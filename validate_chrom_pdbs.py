from pathlib import Path
from Bio.PDB import PDBParser
from math import sqrt

def dist(a, b):
    ax, ay, az = a.get_coord()
    bx, by, bz = b.get_coord()
    return sqrt((ax-bx)**2 + (ay-by)**2 + (az-bz)**2)

OXYGEN_NAME = "O3"
PDB_DIR = Path("chromophore_variants")

parser = PDBParser(QUIET=True)

for pdb_path in sorted(PDB_DIR.glob("*.pdb")):
    name = pdb_path.name
    tag = name.lower()

    structure = parser.get_structure("s", pdb_path)
    residue = next(structure.get_residues())

    oxy = None
    hydrogens = []
    for atom in residue:
        if atom.get_name().strip() == OXYGEN_NAME:
            oxy = atom
        elif getattr(atom, "element", "").upper() == "H":
            hydrogens.append(atom)

    if oxy is None:
        print(f"[FAIL] {name}: oxygen {OXYGEN_NAME} not found")
        continue

    oh_dists = [dist(oxy, h) for h in hydrogens if dist(oxy, h) < 1.25]

    if "deprot" in tag:
        if len(oh_dists) == 0:
            print(f"[PASS] {name} — no H near O")
        else:
            print(f"[FAIL] {name} — has {len(oh_dists)} H near O, distances = {[round(d,2) for d in oh_dists]}")
    elif "prot" in tag:
        if len(oh_dists) == 1 and abs(oh_dists[0] - 0.96) < 0.1:
            print(f"[PASS] {name} — 1 H at {oh_dists[0]:.2f} A")
        else:
            print(f"[FAIL] {name} — {len(oh_dists)} H near O, distances = {[round(d,2) for d in oh_dists]}")
    else:
        print(f"[WARN] {name}: could not determine prot/deprot from filename")
