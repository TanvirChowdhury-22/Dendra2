Amber/TeraChem QM/MM setups for a Dendra2 construct with a covalently bound chromophore:
	
 	•	prot/ – chromophore protonated (phenol, net 0 in isolation)
	•	deprot/ – chromophore deprotonated (phenolate, net −1 in isolation)

Include solvated topologies, minimize, equilibrate, production inputs using Amber sander with TeraChem as the QM engine, and runner scripts.

For setting the environment var, if sander is not on $PATH, :
SANDER=/full/path/to/sander bash prot/MD/run_qmmm_prot.sh


Top-level layout

	Dendra2/
		prot/                    # protonated build
    		MD/                    # QM/MM mdins + runner script
    		build_files/           # inputs used during system construction
    		complex_prot.*         # pre-solvation complex (prmtop/inpcrd/pdb)
    		solvated_complex_prot.*# protonated solvated system (prmtop/inpcrd/pdb)
    		min_local_solv.*       # protonated solvent-only minimization (before QM/MM)
		deprot/                  # deprotonated build
    		MD/                    # QM/MM mdins + runner script
    		build_files/           # (optional) inputs used during system construction
    		complex_deprot.*       # pre-solvation complex (prmtop/inpcrd/pdb)
    		solvated_complex_deprot.*# deprotonated solvated system (prmtop/inpcrd/pdb)
    		min_local_solv_deprot.*# deprotonated solvent-only minimization (before QM/MM)


QM region and charge
	
 	•	QM region for simulation: :60,57,62,140,207
	•	Spin multiplicity: 1 (closed shell)
	•	QM charge (formal):


	Protonated: −1 (5SQ 0, Arg +1, 2×Glu −2, Thr 0 → −1) 


 parmed -p solvated_complex_prot.prmtop -c min_local_solv.rst

                 ./       |      |        \.
               .:(        |i __ j|        ):`.
             .'   `._     |`::::'|     _.'    `.
           .'        "---.j `::' f.---"         `.
     _____/     ___    ____      __    __  ____   ___    
    |      \   |   |  |     `__'|  \  /  ||    | |   \    
    |  .-.  | .'   `| | .-.  |-/|   \/   || ___| |  . \   
    |  |_|  | |  i  | | |_| /"":|        || |    | | \ \  
    |       / | .^. | |    /::::|        || |__. | |  \ \ 
    |  ----'  | | | | |    \ :: |        ||  __| | |  |  )
    |  |     .' ''' `.|  |\ \   |  i  i  j| |    . | /  /  
    |  |     |   _   ||  | \ \  |  |\/|  || |__. | |.  / .
   [|  |     |  | |  ||  |  \ \ |  |  |  ||    | |    /   ].
  ] `--'     :--' `--::--'   \_|`--' ::--"|____|-"-- /    :[
  |      __  ::-'''`.:' "--.    .----::.----:: ,.---._    :|
  [  .-""  "`'              \  /      "      `'       `-. :].
 ]:.'                        \/                          `.:[
 |/                                                        \|

ParmEd: a Parameter file Editor


Loaded Amber topology file solvated_complex.prmtop with coordinates from min_local_solv.rst

Reading input from STDIN...
> strip !(:60,57,62,140,207)
Removing mask '!(:60,57,62,140,207)' (33385 atoms) from the topology file.
> summary
Amino Acid Residues:   4
Nucleic Acid Residues: 0
Number of cations:     0
Number of anions:      0
Num. of solvent mols:  0
Num. of unknown res:   1
Total charge (e-):     -0.9140
Total mass (amu):      851.8660
Number of atoms:       108
Number of residues:    5
Residue set:           5SQ, ARG, GLU, THR
Residue count:         5SQ: 1, ARG: 1, GLU: 2, THR: 1
System volume (ang^3): 381203.01
System density (g/mL): 0.003711

> quit
Done!



	•	Deprotonated: −2 (phenolate −1 + Arg +1 + 2×Glu −2 + Thr 0 → −2)

parmed -p solvated_complex_deprot.prmtop -c min_local_solv.rst       

                                  _____
                         __...---'-----`---...__
                   _===============================
   ______________,/'      `---..._______...---'
  (______________). .    ,--'                            
   /    /.---'       `. /                  
  '--------_  - - - - _/         P A R M E D
            `~~~~~~~~'

ParmEd: a Parameter file Editor


Loaded Amber topology file solvated_complex_deprot.prmtop with coordinates from min_local_solv.rst

Reading input from STDIN...
> strip !(:60,57,62,140,207)
Removing mask '!(:60,57,62,140,207)' (33389 atoms) from the topology file.
> summary
Amino Acid Residues:   4
Nucleic Acid Residues: 0
Number of cations:     0
Number of anions:      0
Num. of solvent mols:  0
Num. of unknown res:   1
Total charge (e-):     -1.7880
Total mass (amu):      850.8580
Number of atoms:       107
Number of residues:    5
Residue set:           ***, ARG, GLU, THR
Residue count:         ***: 1, ARG: 1, GLU: 2, THR: 1
System volume (ang^3): 381233.69
System density (g/mL): 0.003706

> quit
Done!

QM method
	
 	•	DFT: B3LYP
	•	Basis: 6-311G**
	•	Electrostatics: electrostatic embedding (qmmm_int=1), periodic QM/MM (qm_ewald=1)

Reason for the residues selected in the QM region:

	:57 — Thr
		•	Formal charge 0 in both states.
		•	Proton wires and short-range H-bond networks around the chromophore enable ESPT. Immediate H-bond neighbors should be allowed to polarize/respond self-consistently. [https://pubs.acs.org/doi/10.1021/jp309219m]
	:62 — Arg
		•	Formal charge +1 in both states.
		•	Arg residue adjacent to the chromophore should work as a proton donor/stabilizer. It interacts through a salt bridge with the glutamate. [https://pubs.acs.org/doi/10.1021/jp309219m]
	:140 — Glu
		•	Formal charge −1 in both states.
		•	conserved Glu near the chromophore that can undergo photo-Kolbe decarboxylation and participate in electron/proton-transfer chemistry. [https://pubs.acs.org/doi/10.1021/jp309219m]
	:207 — Glu
		•	Same as 140 — Glu

Protonated version simulation
From the project root:

	chmod +x prot/MD/run_qmmm_prot.sh
	bash prot/MD/run_qmmm_prot.sh

it does in order:

	1.	Min1 – restrained QM/MM minimization (starts from prot/min_local_solv.rst)
	2.	Min2 – unrestrained QM/MM minimization
	3.	Eq1 – NVT heat 0 → 300 K (dt=1 fs, light restraints on MM heavy atoms)
	4.	Eq2 – NPT density at 300 K / 1 atm (dt=2 fs, restraints weakened)
	5.	Eq3 – short NPT hold (dt=2 fs, restraints feathered down)
	6.	Prod NVT – 5 ns, dt=2 fs, no restraints
	7.	Prod NPT – 5 ns, dt=2 fs, no restraints
