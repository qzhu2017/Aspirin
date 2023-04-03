# Aspirin
Free energy calculation of aspirin polymorphs

## Useful links
- [Set up environment regarding plumed and lammps](https://github.com/plumed/masterclass-2022)
- [Master class](https://www.plumed.org/doc-master/user-doc/html/masterclass-22-12.html)
- [Plumed examples](https://www.plumed-nest.org/browse.html)
- [Instruction on Collective variables](https://www.plumed.org/doc-v2.8/user-doc/html/colvarintro.html)

## Our own script
```
pip install --upgrade git+https://github.com/qzhu2017/PyXtal.git@master
```

`analysis.py` provides a short utility to 
- parse the MD trajectory files
- compute the q-series parameter for each molecule/dimer center
- plot the histogram for the q values

## Structures
The crystal structures can be found [here](https://github.com/qzhu2017/Aspirin/tree/main/cifs)
