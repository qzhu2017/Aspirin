from collections import deque
from ase.neighborlist import NeighborList
from ase.io.lammpsrun import get_max_index, construct_cell, lammps_data_to_ase_atoms
#from pyxtal.descriptor import _qlm
import numpy as np
from ase.atoms import Atoms

def _qlm(dists, l=4):
    """
    Calculates the vector associated with an atomic site and
    one of its neighbors

    Args:
        distss: a list of distance vectors
        l:  free integer quantum number
    Returns:
        q: numpy array(complex128), the complex vector qlm normalized
            by the number of nearest neighbors
    """
    # initiate variable as a complex number
    q = np.zeros(2 * l + 1, dtype=np.complex128)
    neighbors_count = len(dists)

    for i, m in enumerate(range(-l, l + 1)):
        for _j, r_vec in enumerate(dists):
            # find the position vector of the site/neighbor pair
            r_mag = np.linalg.norm(r_vec)
            theta = np.arccos(r_vec[2] / r_mag)
            if abs((r_vec[2] / r_mag) - 1.0) < 10.0 ** (-8.0):
                theta = 0.0
            elif abs((r_vec[2] / r_mag) + 1.0) < 10.0 ** (-8.0):
                theta = np.pi

            # phi
            if r_vec[0] < 0.0:
                phi = np.pi + np.arctan(r_vec[1] / r_vec[0])
            elif r_vec[0] > 0.0 and r_vec[1] < 0.0:
                phi = 2 * np.pi + np.arctan(r_vec[1] / r_vec[0])
            elif r_vec[0] > 0.0 and r_vec[1] >= 0.0:
                phi = np.arctan(r_vec[1] / r_vec[0])
            elif r_vec[0] == 0.0 and r_vec[1] > 0.0:
                phi = 0.5 * np.pi
            elif r_vec[0] == 0.0 and r_vec[1] < 0.0:
                phi = 1.5 * np.pi
            else:
                phi = 0.0

            q[i] += sph_harm(m, l, phi, theta)
    # normalize by number of neighbors
    return q / neighbors_count

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
                celldatarows = [lines.popleft() for _ in range(3)]
                cells.append(np.loadtxt(celldatarows))

        if len(cells) <= len(images):
            images = images[:len(cells)]
            for mat, struc in zip(cells, images):
                struc.set_pbc([1, 1, 1])
                struc.set_cell(mat)
                #print(struc)
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
        centers[i, :] = np.mean(scaled_pos[:13], axis=0)
        #centers[i, :] = scaled_pos[0]
    return centers

def get_dimer_centers(struc, N_mols, N_atoms, id=10):
    all_scaled_positions = np.zeros([N_atoms*N_mols, 3])
    centers = []
    for i in range(0, N_mols, 2):
        j = 1
        id1, id2 = i*N_atoms+id, (i+j)*N_atoms+id
        pos1 = struc.get_scaled_positions()[id1]
        pos2 = struc.get_scaled_positions()[id2]
        dist = pos2 - pos1
        shift = np.round(dist)
        pos2 -= shift
        #print('dimer', i, i+j, pos1, pos2, np.linalg.norm((pos2-pos1).dot(struc.cell[:])))
        centers.append((pos1 + pos2)/2)
        #centers[i, :] = scaled_pos[0]
    return np.array(centers)

def compute_shift(strucs, N_atoms, ids=[(0, 2), (1, 3)], ref='mol'):

    dists = []
    N_mols = int(len(strucs[0])/N_atoms)
    for i, struc in enumerate(strucs):
        pos = get_dimer_centers(struc, N_mols, N_atoms)
        _dists = []
        #print(struc.cell)
        for id in ids:
            (a, b) = id
            ds = (pos[a]-pos[b])
            ds -= np.round(ds)
            _dists.append(np.dot(ds, struc.cell.array))
            #print(struc.cell.cellpar()[:3], ds, ds*struc.cell.cellpar()[:3], np.dot(ds, struc.cell.array))#; import sys; sys.exit()
            #dists.append(ds*struc.cell.cellpar()[:3])
        dists.append(np.array(_dists).mean(axis=0))
        #dists.append(np.array(_dists)) #.mean(axis=0))
        #ds = (pos[a]+pos[b]-pos[c]-pos[d])/2
        #if ds[1] > 0.3:
        #    ds[1] -= 0.5
        #elif ds[1] < -0.25:
        #    ds[1] += 0.5
        #dist = np.dot(ds, struc.cell.array)
        #dists.append(dist)
        #print(pos); print(dists); import sys; sys.exit()
        #print(i, ds, dist)
        #import sys; sys.exit()
    return np.array(dists)


def compute_qs(strucs, N_atoms, ls=[4, 6], N_cut=10, r=5.0, ref='mol'):

    N_mols = int(len(strucs[0])/N_atoms)
    if ref == 'mol':
        qs = np.zeros([len(strucs), N_mols, len(ls)])
    else:
        qs = np.zeros([len(strucs), int(N_mols/2), len(ls)])
        #N_cut = 20 #check all

    neighbors = NeighborList([r]*qs.shape[1], self_interaction=False, bothways=True, skin=0.0)

    for i, struc in enumerate(strucs):
        if ref == 'mol':
            pos = get_mol_centers(struc, N_mols, N_atoms)
        else:
            pos = get_dimer_centers(struc, N_mols, N_atoms)
        cell = struc.cell[:]
        #cell[0, 0] /= 1.5
        atom1 = Atoms([6]*len(pos), scaled_positions=pos, cell=cell, pbc=[1,1,1])
        #atom1.write('0.cif', format='cif')
        neighbors.update(atom1)

        for j in range(len(pos)):
            indices, offsets = neighbors.get_neighbors(j)
            Ri = atom1.positions[j]
            tmp = np.zeros([len(indices), 3])
            count = 0
            for k, offset in zip(indices, offsets):
                tmp[count] = atom1.positions[k] + np.dot(offset, atom1.get_cell())
                count += 1
            tmp -= Ri
            ds = np.linalg.norm(tmp, axis=1)

            print(i, j, len(ds), np.sort(ds))

            if N_cut < len(ds):
                tol = np.sort(ds)[N_cut-1]
                tmp = tmp[ds<tol+1e-5]
            for l0, l in enumerate(ls):
                factor = (4 * np.pi) / (2*l + 1)
                qlms = _qlm(tmp, l)
                dot = float(np.sum(qlms*np.conjugate(qlms)))
                val = np.sqrt(factor * dot)
                qs[i, j, l0] = val
                #print(i, j, len(tmp), val)
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
    q4_range = [0.12, 0.52]
    q6_range = [0.12, 0.52]

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
    ax.set_xlabel('q'+str(ls[0]))
    ax.set_ylabel('q'+str(ls[1]))
    ax.set_xlim(q4_range)
    ax.set_ylim(q6_range)
    plt.figtext(0.77, 0.68, title.replace('_', '\n'), ha='center', fontsize=16)
    plt.savefig(title+'.png', dpi=300)

def plot_yz(titles, suffix):
    import matplotlib.pyplot as plt
    import seaborn as sns
    import matplotlib.gridspec as gridspec
    sns.set(font_scale=1.0)
    sns.set_theme(style='white')

    fig = plt.figure(figsize=(5.6, 5.6))
    gs = gridspec.GridSpec(nrows=2, ncols=2,
                           width_ratios=[1, 0.5], height_ratios=[0.5, 1],
                           wspace=0, hspace=0)

    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[0, 0])
    ax3 = fig.add_subplot(gs[1, 1])
    ori2 = 'vertical'
    drange = (-1.5, 1.5)
    nbins = 75
    ori3 = 'horizontal'

    for i, title in enumerate(titles):
        if i == 0:
            col = 'red'
        else:
            col = 'blue'
        d = np.loadtxt(title+'.txt')
        ax1.scatter(d[:,2], d[:,1], s=1.0, c=col, alpha=0.5, label=title)
        ax2.hist(d[:,2], bins=nbins, range=drange, density=True, color=col, label=title+'-c', orientation=ori2, alpha=0.5)
        ax3.hist(d[:,1], bins=nbins, range=drange, density=True, color=col, label=title+'-b', orientation=ori3, alpha=0.5)

    ax1.set_xlim(drange)
    ax1.set_ylim(drange)
    ax2.set_ylim([0.001, 2])
    ax3.set_xlim([0.001, 2])
    ax2.set_xticks([])
    ax3.set_yticks([])
    ax1.legend(loc=1)
    ax2.legend(loc=1)
    ax3.legend(loc=2)
    ax1.grid()
    ax2.grid()
    ax3.grid()
    #ax1.set_aspect('equal')
    ax1.set_aspect(1./ax1.get_data_ratio())
    plt.savefig('fluc'+suffix+'.png')


###########################################################
N_atoms = 21 # number of atoms per molecule
cut = 100 #600
#files = ['MLP/asp_I_300K.lammpstrj', 'MLP/asp_II_300K.lammpstrj']
#files = ['DFTB/ACSALA/', 'DFTB/ACSALA13/']

ids = [(0, 2), (1, 3)]#, (4, 6), (5, 7), (8, 10), (9, 11), (12, 14), (13, 15)]
for T in [300]: #[300, 350, 380]:
    labels = []
    for n in ['I', 'II']:
        label = n+'-'+str(T)+'K' #['I-380K', 'II-380K']
        if n == 'I':
            form = 'ACSALA'
        else:
            form = 'ACSALA13'
            #ids = [(0, 2), (1, 3)]
        #f = 'GAFF/'+form+'-'+str(T)+'K/dump.lammpstrj'
        f = 'GAFF/'+form+'/dump.lammpstrj'
        strucs = read_lammps_dump_text(f)#[:1]
        dists = compute_shift(strucs[cut:], N_atoms, ids)
        #if i == 0: dists[:,1] += 0.1
        print(f, cut, len(dists), np.max(dists, axis=0), np.min(dists, axis=0))
        np.savetxt(label+'.txt', dists)
        labels.append(label)
    plot_yz(labels, '-'+str(T))

    #qs = compute_qs(strucs, N_atoms, N_cut=N_cut, ls=ls, r=r_cut, ref=ref)
    #pre = f.split('.')[0] + '_Ncut_' + str(N_cut)
    #for l0, l in enumerate(ls):
    #    np.savetxt(pre+'-q'+str(l)+'.txt', qs[:,:,l0])
    #prefixes.append(pre)

#title = dir0+'_NPT-MD_300K_1atm'
#plot_q4_q6(prefixes, labels, title, ls)
