The chromophore_variants folder contains the PDB files of the protonated and deprotonated versions of the chains (A-H)

python chromophore_extractor_v2.py 2vzx.pdb --resname 5SQ --resid 64 --oxygen-name OH

python validate_chrom_pdbs_any.py --oxygen_name OH
