# Electronic_Coupling

Used to calculat the electronc coupling between two molecules.

Calculations are based on files generated by gaussian program.

Two kinds of fiels are needed: fchk file and log file, both of which can be generated by gaussian.

6 files are needed:
monomer1.fchk
monomer1.log
monomer2.fchk
monomer2.log
dimer.fchk
dimer.log

When doing calculation using gaussian, do add the following keywords for both dimer and monomers:
```
# b3lyp/6-31* iop(3/33=1) punch=mo nosymm pop=full
```


