# Set up

```
conda create --name plumed-masterclass-2022
conda activate plumed-masterclass-2022
conda install --strict-channel-priority -c plumed/label/masterclass-2022 -c conda-forge plumed lammps
```

# Run lammps (to get the trajectory, e.g., dump files)

Goto `ACSACA` or `ACAALA13`, make sure you have correct input file (e.g., `npt.in`) then run
```
mpirun -np 4 lmp_mpi < npt.in > npt.log
```

# Run plumed (to compute the CV from lammps dump)
Goto `ACSACA` or `ACAALA13`, make sure you have correct input files
 - `plumed.dat`: control the plumed output
 - `dump.na`: lammps dump file
 - `out.dcd`: lammps output 

then run
```
plumed driver --plumed plumed.dat --mf_dcd out.dcd > plumed.out
```

You expect to see some files

## Plumed
Here is a typical `plumed.dat`,

```
ENVIRONMENTSIMILARITY ...
 SPECIES=1-1024:21
 SIGMA=0.05
 CRYSTAL_STRUCTURE=CUSTOM
 REFERENCE_1=../env_01.pdb 
 LABEL=s1
 MEAN
... ENVIRONMENTSIMILARITY


hh1: HISTOGRAM ...
   DATA=s1
   GRID_MIN=-0.5
   GRID_MAX=2
   GRID_BIN=1000
   BANDWIDTH=0.01
...

DUMPGRID GRID=hh1 FILE=histo1
PRINT STRIDE=1  ARG=* FILE=COLVAR

```

This will consider the reference atoms
```
1, 22, 43, ..., 
```
and then compute the similarity CV based on the reference environment from `../env_01.pdb`

It will generate two output files

- histo1
- COLVAR


