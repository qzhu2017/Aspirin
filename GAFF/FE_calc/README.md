# Aspirin
Free energy calculation of aspirin polymorphs

## Setup
```
$ conda activate plumed-masterclass-2022
$ mpirun -np 24 lmp < npt.in > npt.out &
$ tail -f COLVAR
```

- Analyze the output `COLVAR` to check the if transition happens and how frequent it is
- download the lammps dumpfile to ovito check if the simulation is reasonable
- Play with the parameters `BARRIER` in `plumed.dat`


# Download
```
cp dump.lammpstrj ~/
scp qzhu@zinc.physics.unlv.edu:dump.lammpstrj ./
```
