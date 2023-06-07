from pyxtal.db import database
from pyxtal import pyxtal
db = database('ht.db')

for code in ['ACSALA', 'ACSALA13']:
    f1 = db.get_pyxtal(code) 
    f1.to_file(code+'-pyxtal.cif')
    c1 = pyxtal()
    sites = [{}]

    # drop the cif file with molecuar centers
    s = f1.mol_sites[0]
    sites[0][str(s.wp.multiplicity)+s.wp.letter] = s.position
    print(sites)
    c1.from_random(3, f1.group.number, ['C'], f1.numMols, lattice=f1.lattice, sites=sites)
    c1.to_file(code+'-center.cif')

    # drop the cif file of the atoms by ID
    id = 0
    s = f1.mol_sites[0]
    coord0s, specie0s = s._get_coords_and_species(first=True)
    sites[0][str(s.wp.multiplicity)+s.wp.letter] = coord0s[id]
    print(sites)
    c1.from_random(3, f1.group.number, [specie0s[id]], f1.numMols, lattice=f1.lattice, sites=sites)
    c1.to_file(code+'-atom'+str(id)+'.cif')
