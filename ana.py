from collections import deque
from ase.neighborlist import NeighborList
from ase.io.lammpsrun import get_max_index, construct_cell, lammps_data_to_ase_atoms
from pyxtal.descriptor import _qlm
import numpy as np
from ase.atoms import Atoms

def read_xyz(filename, cellname=None):
    fileobj = open(filename)
    lines = fileobj.readlines()
    fileobj.close()
    images = []
    while len(lines) > 0:
        symbols = []
        positions = []
        natoms = int(lines.pop(0))
        lines.pop(0)  # Comment line; ignored
        for _ in range(natoms):
            line = lines.pop(0)
            symbol, x, y, z = line.split()[:4]
            symbol = symbol.lower().capitalize()
            symbols.append(symbol)
            positions.append([float(x), float(y), float(z)])
        images.append(Atoms(symbols=symbols, positions=positions))

    # update the cell vectors
    if cellname is not None:
        fileobj = open(cellname)
        lines = deque(fileobj.readlines())
        fileobj.close()
        cells = []

        while len(lines) > 0:
            line = lines.popleft()
            if 'Lattice vectors (A)' in line:
                print(line)
                celldatarows = [lines.popleft() for _ in range(3)]
                cells.append(np.loadtxt(celldatarows))

        if len(cells) <= len(images):
            images = images[:len(cells)]
            for mat, struc in zip(cells, images):
                struc.set_pbc([1, 1, 1])
                struc.set_cell(mat)
                print(struc)
        else:
            print(len(cells), len(images))
            raise RuntimeError('Number of structure is inconsistent')

    return images

def read_lammps_dump_text(filename, **kwargs):
    # Load all dumped timesteps into memory simultaneously
    fileobj = open(filename)
    lines = deque(fileobj.readlines())
    fileobj.close()
    index_end = get_max_index(-1)

    n_atoms = 0
    images = []

    # avoid references before assignment in case of incorrect file structure
    cell, celldisp, pbc = None, None, False

    while len(lines) > n_atoms:
        line = lines.popleft()

        if "ITEM: TIMESTEP" in line:
            n_atoms = 0
            line = lines.popleft()
            # !TODO: pyflakes complains about this line -> do something
            # ntimestep = int(line.split()[0])  # NOQA

        if "ITEM: NUMBER OF ATOMS" in line:
            line = lines.popleft()
            n_atoms = int(line.split()[0])

        if "ITEM: BOX BOUNDS" in line:
            # save labels behind "ITEM: BOX BOUNDS" in triclinic case
            # (>=lammps-7Jul09)
            tilt_items = line.split()[3:]
            celldatarows = [lines.popleft() for _ in range(3)]
            celldata = np.loadtxt(celldatarows)
            diagdisp = celldata[:, :2].reshape(6, 1).flatten()

            # determine cell tilt (triclinic case!)
            if len(celldata[0]) > 2:
                # for >=lammps-7Jul09 use labels behind "ITEM: BOX BOUNDS"
                # to assign tilt (vector) elements ...
                offdiag = celldata[:, 2]
                # ... otherwise assume default order in 3rd column
                # (if the latter was present)
                if len(tilt_items) >= 3:
                    sort_index = [tilt_items.index(i)
                                  for i in ["xy", "xz", "yz"]]
                    offdiag = offdiag[sort_index]
            else:
                offdiag = (0.0,) * 3

            cell, celldisp = construct_cell(diagdisp, offdiag)

            # Handle pbc conditions
            if len(tilt_items) == 3:
                pbc_items = tilt_items
            elif len(tilt_items) > 3:
                pbc_items = tilt_items[3:6]
            else:
                pbc_items = ["f", "f", "f"]
            pbc = ["p" in d.lower() for d in pbc_items]

        if "ITEM: ATOMS" in line:
            colnames = line.split()[2:]
            datarows = [lines.popleft() for _ in range(n_atoms)]
            data = np.loadtxt(datarows, dtype=str)
            out_atoms = lammps_data_to_ase_atoms(
                data=data,
                colnames=colnames,
                cell=cell,
                celldisp=celldisp,
                atomsobj=Atoms,
                pbc=pbc,
                **kwargs
            )
            images.append(out_atoms)

        if len(images) > index_end >= 0:
            break

    return images

def get_mol_centers(struc, N_mols, N_atoms):
    all_scaled_positions = np.zeros([N_atoms*N_mols, 3])
    centers = np.zeros([N_mols, 3])
    for i in range(N_mols):
        start, end = i*N_atoms, (i+1)*N_atoms
        scaled_pos = struc.get_scaled_positions()[start:end, :]
        dist = scaled_pos - scaled_pos[0]
        shift = np.round(dist)
        scaled_pos -= shift
        centers[i, :] = np.mean(scaled_pos, axis=0)
    return centers

def compute_q4_q6(strucs, N_atoms, ls=[4, 6], N_cut=10):

    N_mols = int(len(strucs[0])/N_atoms)
    neighbors = NeighborList([5.0]*N_mols, self_interaction=False, bothways=True, skin=0.0)
    qs = np.zeros([len(strucs), N_mols, len(ls)])

    for i, struc in enumerate(strucs):
        pos = get_mol_centers(struc, N_mols, N_atoms)
        atom1 = Atoms([6]*N_mols, scaled_positions=pos, cell=struc.cell[:], pbc=[1,1,1])
        neighbors.update(atom1)
        for j in range(N_mols):
            indices, offsets = neighbors.get_neighbors(j)
            Ri = atom1.positions[j]
            tmp = np.zeros([len(indices), 3])
            count = 0
            for k, offset in zip(indices, offsets):
                tmp[count] = atom1.positions[k] + np.dot(offset, atom1.get_cell())
                count += 1
            tmp -= Ri
            ds = np.linalg.norm(tmp, axis=1)
            #print(np.sort(ds))
            tol = np.sort(ds)[N_cut-1]
            tmp = tmp[ds<tol+1e-5]
            for l0, l in enumerate(ls):
                factor = (4 * np.pi) / (2*l + 1)
                qlms = _qlm(tmp, l)
                dot = float(np.sum(qlms*np.conjugate(qlms)))
                val = np.sqrt(factor * dot)
                qs[i, j, l0] = val
                print(i, j, len(tmp), val)
    return qs


def plot_q4_q6(prefixes, labels, title, ls=[4, 6]):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import matplotlib.gridspec as gridspec
    sns.set(font_scale=1.0)
    sns.set_theme(style='white')
    
    fig = plt.figure(figsize=(6.4, 5.0))
    gs = gridspec.GridSpec(nrows=len(ls), ncols=2, 
                           width_ratios=[1, 0.5], height_ratios=[0.5, 1], 
                           wspace=0, hspace=0)
    
    drange = (0, 0.5)
    q4_range = [0.25, 0.52]
    q6_range = [0.25, 0.52]

    for l0, l in enumerate(ls):
        if l0 == 0:
            ax = fig.add_subplot(gs[0, 0])
            ori = 'vertical'
            ax.set_xlim(q4_range)
            ax.set_xticks([])
        else:
            ax = fig.add_subplot(gs[1, 1])
            ori = 'horizontal'
            ax.set_ylim(q6_range)
            ax.set_yticks([])
    
        for f in prefixes:
            tag = f+'-q'+str(l)+'.txt'
            d = np.loadtxt(tag).flatten()
            ax.hist(d, bins=200, range=drange, density=True, label=tag[:-4], alpha=0.5, orientation=ori)
    
    ax = fig.add_subplot(gs[1, 0])
    for f, la in zip(prefixes, labels):
        x = None
        y = None
        for l0, l in enumerate(ls):
            tag = f+'-q'+str(l)+'.txt'
            if l0 == 0:
                x = np.loadtxt(tag)
            else:
                y = np.loadtxt(tag)
        ax.scatter(x, y, label=la, alpha=0.5, s=0.1)
    
    legend = ax.legend(frameon=False, markerscale=16, fontsize=15)
    ax.set_xlabel('q4')
    ax.set_ylabel('q6')
    ax.set_xlim(q4_range)
    ax.set_ylim(q6_range)
    plt.figtext(0.77, 0.68, title.replace('_', '\n'), ha='center', fontsize=16)
    plt.savefig(title+'.png', dpi=300)

###########################################################
N_atoms = 21 #number of atoms per molecule
ls = [4, 6]  #q4 and q6
N_cut = 10   #number of neighbors

#strucs = read_xyz('ACS/geo_final.xyz', 'ACS/md.out')
#print(len(strucs))

labels = ['form I', 'form II']
prefixes = []

# MLP
files = ['MLP/asp_I_300K.lammpstrj', 'MLP/asp_II_300K.lammpstrj']
title = 'MLP_NPT-MD_300K_1atm'

# GAFF
files = ['GAFF/ACSALA/dump.lammpstrj', 'GAFF/ACSALA13/dump.lammpstrj']
title = 'GAFF_NPT-MD_300K_1atm'

for f in files:
    strucs = read_lammps_dump_text(f)
    qs = compute_q4_q6(strucs, N_atoms, N_cut=N_cut, ls=ls)
    pre = f.split('.')[0] + '_Ncut_' + str(N_cut)
    for l0, l in enumerate(ls):
        np.savetxt(pre+'-q'+str(l)+'.txt', qs[:,:,l0])
    prefixes.append(pre)

plot_q4_q6(prefixes, labels, title)
