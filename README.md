# [MD simulation of Pyrazinamidase](./3pl1.pdb)
![Pyrazinamidase](/home/meco/Downloads/Bioinformatics/Pyrazinamidase/images/pyrazinamidase.png  "Pyrazinamidase Structure")
---
## [1. Delete Water molecules-residues HOH in PDB file](./3pl1_clean.pdb)
```
grep -v HOH 3pl1.pdb > 3pl1_clean.pdb
```

---
## [2. Fixing missing or non-bonded atoms, and adding missing hidrogens using PDBFixer](./3pl1_fixed.pdb)
* non-bonded atoms
* filter backone chain

---
## [3. Execute pdb2gmx](./3pl1.gro)
![Water](/home/meco/Downloads/Bioinformatics/Pyrazinamidase/images/water.png  "Water models")
```
gmx pdb2gmx -f 3pl1_fixed.pdb -o 3pl1_processed.gro -water tip3 -ignh
```

---
## [3. Define the box dimentions](./3pl1_newbox.gro)
```
gmx editconf -f 3pl1_processed.gro -o 3pl1_newbox.gro -c -d 1.0 -bt cubic
```

---
## [4. Fill the box with water](./3pl1_solv.gro)
The TIP3P water model as implemented in CHARMM (MacKerell) specifies a 3-site rigid water molecule with charges and Lennard-Jones parameters assigned to each of the 3 atoms. In GROMACS the fix shake command can be used to hold the two O-H bonds and the H-O-H angle rigid. A bond style of harmonic and an angle style of harmonic or charmm should also be used.

![Box](/home/meco/Downloads/Bioinformatics/Pyrazinamidase/images/box.png  "Simulation box with solvated protein")
```
gmx solvate -cp 3pl1_newbox.gro -cs spc216.gro -o 3pl1_solv.gro -p topol.top
```

* Using **AMBER99 SB-IDLN** Force Field: Option 6

These are the additional parameters (in real units) to set for O and H atoms and the water molecule to run a rigid TIP3P model with a cutoff. The K values can be used if a flexible TIP3P model (without fix shake) is desired. If the LJ epsilon and sigma for HH and OH are set to 0.0, it corresponds to the original 1983 TIP3P model (Jorgensen).
```
O mass = 15.9994
H mass = 1.008
O charge = -0.834
H charge = 0.417
LJ e of OO = 0.1521
LJ s of OO = 3.1507
LJ e of HH = 0.0460
LJ s of HH = 0.4000
LJ e of OH = 0.0836
LJ s of OH = 1.7753
K of OH bond = 450
ro of OH bond = 0.9572
K of HOH angle = 55
0 of HOH angle = 104.52
```

These are the parameters to use for TIP3P with a long-range Coulomb solver:
```
O mass = 15.9994
H mass = 1.008
O charge = -0.830
H charge = 0.415
LJ e of OO = 0.102
LJ s of OO = 3.188
LJ e, s of OH, HH = 0.0
K of OH bond = 450
ro of OH bond = 0.9572
K of HOH angle = 55
0 of HOH angle = 104.52
```
---
## [5. Adding Ions to the simulation box](./ions.tpr)
Assemble the solvated electroneutral system.

![ions](/home/meco/Downloads/Bioinformatics/Pyrazinamidase/images/ions.png  "Neutralized System")

### [ions.mdp File](./ions.mdp)
```
; ions.mdp - used as input into grompp to generate ions.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 20         ; Frequency to update the neighbor list and long range forces
cutoff-scheme	= Verlet    ; Buffered neighbor searching 
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
```
```
gmx grompp -f ions.mdp -c 3pl1_solv.gro -p topol.top -o ions.tpr -maxwarn 1
```
```
gmx genion -s ions.tpr -o 3pl1_solv_ions.gro -p topol.top -pname NA -nname CL -neutral
```

---
## [6. Energy Minimization](./em.trr)
Relax the structure through Energy Minimization

### [minim.mdp](./minim.mdp)
```
; minim.mdp - used as input into grompp to generate em.tpr
; Parameters describing what to do, when to stop and what to save
integrator  = steep         ; Algorithm (steep = steepest descent minimization)
emtol       = 1000.0        ; Stop minimization when the maximum force < 1000.0 kJ/mol/nm
emstep      = 0.01          ; Minimization step size
nsteps      = 50000         ; Maximum number of (minimization) steps to perform

; Parameters describing how to find the neighbors of each atom and how to calculate the interactions
nstlist         = 20        ; Frequency to update the neighbor list and long range forces
cutoff-scheme   = Verlet    ; Buffered neighbor searching
coulombtype     = PME       ; Treatment of long range electrostatic interactions
rcoulomb        = 1.0       ; Short-range electrostatic cut-off
rvdw            = 1.0       ; Short-range Van der Waals cut-off
pbc             = xyz       ; Periodic Boundary Conditions in all 3 dimensions
```

#### [Minimazing the Potential](./em.tpr)
```
gmx grompp -f minim.mdp -c 3pl1_solv_ions.gro -p topol.top -o em.tpr
```
```
gmx mdrun -v -deffnm em
```
We'll get the following files:
* **em.log:** ASCII-text log file of the EM process
* **em.edr:** Binary energy file
* **em.trr:** Binary full-precision trajectory
* **em.gro:** Energy minimized structure

### [Visualizing the Potential progression using Grace](./potential.xmg)
```
gmx energy -f em.edr -o potential.xvg
```
```
xmgrace potential.xvg
```

![Potential Progression](/home/meco/Downloads/Bioinformatics/Pyrazinamidase/images/potential.png  "Potential Progression")

## [7. NVT ensemble: Constant Number of particles, Volume, and Temperature](./nvt.trr)
Equilibrate the solvent and ions around the protein. It needs to be bought to the temperature, you wish to simulate and stablish the proper orientation about the solute.

### [nvt.md File](./nvt.mdp)
```
title                   = AMBER Pyrazinamidase NVT Equilibration
define                  = -DPOSRES  ; position restrain the protein
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 50000     ; 2 * 50000 = 100 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 500       ; save coordinates every 1.0 ps
nstvout                 = 500       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = no        ; first dynamics run
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings 
cutoff-scheme           = Verlet    ; Buffered neighbor searching
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is off
pcoupl                  = no        ; no pressure coupling in NVT
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = yes       ; assign velocities from Maxwell distribution
gen_temp                = 300       ; temperature for Maxwell distribution
gen_seed                = -1        ; generate a random seed
```

#### [Minimazing the Temperature](./nvt.tpr)
```
gmx grompp -f nvt.mdp -c em.gro -r em.gro -p topol.top -o nvt.tpr
```
```
gmx mdrun -v -deffnm nvt
```
We'll get the following files:
* **nvt.log:** ASCII-text log file of the NVT process
* **nvt.edr:** Binary nvt file
* **nvt.trr:** Binary full-precision trajectory
* **nvt.gro:** NVT minimized structure

### [Visualizing Temperature progression using Grace](./temperature.xmg)
```
gmx energy -f nvt.edr -o temperature.xvg
```
```
xmgrace temperature.xvg
```

![Temperature](/home/meco/Downloads/Bioinformatics/Pyrazinamidase/images/temperature.png  "Temperature Progression")

---
## [8. NPT ensemble: Constant Number of Particles, Pressure, and Temperature](./npt.trr)
After the model arrives at the correct temperature based on kinetic energy, you'll apply pressure to the system until it reaches the proper density.

### [npt.md File](./npt.mdp)
```
title                   = AMBER Pyrazinamidase NPT Equilibration
define                  = -DPOSRES  ; position restrain the protein
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 250000    ; 2 * 250000 = 500 ps
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 500       ; save coordinates every 1.0 ps
nstvout                 = 500       ; save velocities every 1.0 ps
nstenergy               = 500       ; save energies every 1.0 ps
nstlog                  = 500       ; update log file every 1.0 ps
; Bond parameters
continuation            = yes       ; Restarting after NVT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Nonbonded settings 
cutoff-scheme           = Verlet    ; Buffered neighbor searching
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
DispCorr                = EnerPres  ; account for cut-off vdW scheme
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
refcoord_scaling        = com
; Periodic boundary conditions
pbc                     = xyz       ; 3-D PBC
; Velocity generation
gen_vel                 = no        ; Velocity generation is off 
```
#### [Minimazing the Pressure](./npt.tpr)
```
gmx grompp -f npt.mdp -c nvt.gro -r nvt.gro -t nvt.cpt -p topol.top -o npt.tpr
```
```
gmx mdrun -v -deffnm npt
```
We'll get the following files:
* **npt.log:** ASCII-text log file of the NPT process
* **npt.edr:** Binary npt file
* **npt.trr:** Binary full-precision trajectory
* **npt.gro:** NPT minimized structure

### [Visualizing Pressure progression using Grace](./pressure.xmg)
```
gmx energy -f npt.edr -o pressure.xvg
```
```
xmgrace pressure.xvg
```

![Pressure](/home/meco/Downloads/Bioinformatics/Pyrazinamidase/images/pressure.png  "TPressure Progression")

### [Visualizing Energy progression using Grace](./energy.xmg)
```
gmx energy -f npt.edr -o energy.xvg
```
```
xmgrace energy.xvg
```
![Energy Progression](/home/meco/Downloads/Bioinformatics/Pyrazinamidase/images/energy.png  "Energy Progression")

### [Visualizing Density progression using Grace](./density.xmg)
```
gmx energy -f npt.edr -o density.xvg
```
```
xmgrace density.xvg
```

![Density Progression](/home/meco/Downloads/Bioinformatics/Pyrazinamidase/images/density.png  "Density Progression")

---
## [9. Run MD Simulation](./md_0.trr)
Upon completion of the 2 equilibration phases, the system is now well-equilibrated at the desired temperature and pressure.

### [md.mdp File](./md.mdp)
```
title                   = AMBER Pyrazinamidase MD 
; Run parameters
integrator              = md        ; leap-frog integrator
nsteps                  = 5000000    ; 2 * 5000000 = 10 ns
dt                      = 0.002     ; 2 fs
; Output control
nstxout                 = 0         ; suppress bulky .trr file by specifying 
nstvout                 = 0         ; 0 for output frequency of nstxout,
nstfout                 = 0         ; nstvout, and nstfout
nstenergy               = 5000      ; save energies every 10.0 ps
nstlog                  = 5000      ; update log file every 10.0 ps
nstxout-compressed      = 5000      ; save compressed coordinates every 10.0 ps
compressed-x-grps       = System    ; save the whole system
; Bond parameters
continuation            = yes       ; Restarting after NPT 
constraint_algorithm    = lincs     ; holonomic constraints 
constraints             = h-bonds   ; bonds involving H are constrained
lincs_iter              = 1         ; accuracy of LINCS
lincs_order             = 4         ; also related to accuracy
; Neighborsearching
cutoff-scheme           = Verlet    ; Buffered neighbor searching
nstlist                 = 10        ; 20 fs, largely irrelevant with Verlet scheme
rcoulomb                = 1.0       ; short-range electrostatic cutoff (in nm)
rvdw                    = 1.0       ; short-range van der Waals cutoff (in nm)
; Electrostatics
coulombtype             = PME       ; Particle Mesh Ewald for long-range electrostatics
pme_order               = 4         ; cubic interpolation
fourierspacing          = 0.16      ; grid spacing for FFT
; Temperature coupling is on
tcoupl                  = V-rescale             ; modified Berendsen thermostat
tc-grps                 = Protein Non-Protein   ; two coupling groups - more accurate
tau_t                   = 0.1     0.1           ; time constant, in ps
ref_t                   = 300     300           ; reference temperature, one for each group, in K
; Pressure coupling is on
pcoupl                  = Parrinello-Rahman     ; Pressure coupling on in NPT
pcoupltype              = isotropic             ; uniform scaling of box vectors
tau_p                   = 2.0                   ; time constant, in ps
ref_p                   = 1.0                   ; reference pressure, in bar
compressibility         = 4.5e-5                ; isothermal compressibility of water, bar^-1
```
```
gmx grompp -f md.mdp -c npt.gro -t npt.cpt -p topol.top -o md_0.tpr
```
```
gmx mdrun -v -deffnm md_0
```

---
### [10. Trjconv module to account for any periodicity in the system](./md_0_noPBC.xtc)
Use trjconv module as a post-processing tool string-out coordinates, correct for periodicity, or manually alter the trajectory (time units, frame frequency, etc)
```
gmx trjconv -s md_0.tpr -f md_0.xtc -o md_0_noPBC.xtc -pbc mol -center
```
**NOTICE:** Select 1 for protein and 0 for system 

---
### [11. RMSD: Root Mean Square Diviation](./rmsd.xvg)
```
gmx rms -s md_0.tpr -f md_0_noPBC.xtc -o rmsd.xvg -tu ns
```
```
xmgrace rmsd.xvg
```
**NOTICE:** Select 4 for backbone and 4 for backbone
* **tu**: time units -> ns: nano-seconds

![RMSD](/home/meco/Downloads/Bioinformatics/Pyrazinamidase/images/rmsd.png  "Root Mean Square Deviations")

---
### [12. RMSF: Root Mean Square Fluctuation](./rmsf.xvg)
```
gmx rmsf -s md_0.tpr -f md_0_noPBC.xtc -o rmsf.xvg
```
```
xmgrace rmsf.xvg
```
**NOTICE:** Select 4 for backbone and 4 for backbone
* **tu**: time units -> ns: nano-seconds

```
gmx rmsf -s em.tpr -f md_0_noPBC.xtc -o rmsf_em.xvg
```
```
xmgrace rmsf_em.xvg
```

![RMSF](/home/meco/Downloads/Bioinformatics/Pyrazinamidase/images/rmsf.png  "Root Mean Square Fluctuations")

---
### [13. Radius of Gyration](./gyrate.xvg)
Radius of gyration of a protein is a meassure of its compactness. If the protein is stably folded
```
gmx gyrate -s md_0.tpr -f md_0_noPBC.xtc -o gyrate.xvg
```
```
xmgrace gyrate.xvg
```

![Gyration](/home/meco/Downloads/Bioinformatics/Pyrazinamidase/images/gyrate.png  "Radius of Gyration")

---
### [14. Visualize the simulation using VMD](./md_0.gro)
```
vmd md_0.gro md_0_noPBC.xtc
```

---
## Authors

* **Robinson Montes** - [mecomonteshbtn](https://github.com/mecomonteshbtn)
