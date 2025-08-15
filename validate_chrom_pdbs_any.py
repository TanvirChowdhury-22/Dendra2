from pathlib import Path
import argparse
from math import sqrt
from Bio.PDB import PDBParser

ELEMENT2 = {
    "CL":"Cl","BR":"Br","ZN":"Zn","MG":"Mg","NA":"Na","CA":"Ca","FE":"Fe","CU":"Cu","NI":"Ni","MN":"Mn",
    "SI":"Si","SE":"Se","CO":"Co","AL":"Al","PT":"Pt","PD":"Pd","HG":"Hg","AG":"Ag","AU":"Au","PB":"Pb",
    "SN":"Sn","TI":"Ti","CR":"Cr","LI":"Li","CS":"Cs","RB":"Rb","ZR":"Zr","MO":"Mo","RU":"Ru","RH":"Rh",
    "IR":"Ir","OS":"Os","RE":"Re","TA":"Ta","GA":"Ga","GE":"Ge","AS":"As","SR":"Sr","BA":"Ba","K":"K",
    "F":"F","I":"I"
}

def element_of(atom):
    if hasattr(atom, "element") and atom.element:
        e = atom.element.strip()
        if len(e) == 2:
            return ELEMENT2.get(e.upper(), e.capitalize())
        return e[0].upper()
    full = atom.get_fullname()
    letters = "".join(ch for ch in full if ch.isalpha())
    if len(letters) >= 2 and letters[:2].upper() in ELEMENT2:
        return ELEMENT2[letters[:2].upper()]
    return letters[0].upper() if letters else "X"

def is_h(atom):
    return element_of(atom).upper() == "H"

def dist(a, b):
    ax, ay, az = a.get_coord()
    bx, by, bz = b.get_coord()
    return sqrt((ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2)

def analyze_one(pdb_path: Path, oxygen_name: str, ideal_oh: float, tol: float):
    parser = PDBParser(QUIET=True)
    s = parser.get_structure("s", pdb_path)
    # these files are single-residue structures; grab the first residue
    try:
        r = next(s.get_residues())
    except StopIteration:
        return ("FAIL", "no residue found")

    O = None
    Hs = []
    for a in r.get_atoms():
        nm = a.get_name().strip()
        if nm == oxygen_name:
            O = a
        elif is_h(a):
            Hs.append(a)

    if O is None:
        return ("FAIL", f"oxygen {oxygen_name} not found")

    oh_dists = [dist(O, h) for h in Hs if dist(O, h) < 1.25]
    lname = pdb_path.name.lower().rstrip(".pdb")

    if lname.endswith("_prot"):
        if len(oh_dists) == 1 and abs(oh_dists[0] - ideal_oh) < tol:
            return ("PASS", f"1 H at {oh_dists[0]:.2f} Å")
        return ("FAIL", f"{len(oh_dists)} H near {oxygen_name}, distances = {[round(d,2) for d in oh_dists]}")
    if lname.endswith("_deprot"):
        if len(oh_dists) == 0:
            return ("PASS", "no H near O (as expected)")
        return ("FAIL", f"has {len(oh_dists)} H near {oxygen_name}, distances = {[round(d,2) for d in oh_dists]}")

    if "deprot" in lname and "prot" not in lname:
        if len(oh_dists) == 0:
            return ("PASS", "no H near O (as expected)")
        return ("FAIL", f"has {len(oh_dists)} H near {oxygen_name}, distances = {[round(d,2) for d in oh_dists]}")
    if "prot" in lname and "deprot" not in lname:
        if len(oh_dists) == 1 and abs(oh_dists[0] - ideal_oh) < tol:
            return ("PASS", f"1 H at {oh_dists[0]:.2f} Å")
        return ("FAIL", f"{len(oh_dists)} H near {oxygen_name}, distances = {[round(d,2) for d in oh_dists]}")

    return ("INFO", f"{len(oh_dists)} H within 1.25 Å of {oxygen_name}: {[round(d,2) for d in oh_dists]}")

def main():
    ap = argparse.ArgumentParser(description="Validate protonated/deprotonated chromophore PDBs around a chosen oxygen.")
    ap.add_argument("--pdb_dir", type=Path, default=Path("chromophore_variants"),
                    help="Directory with generated PDBs (default: chromophore_variants)")
    ap.add_argument("--oxygen_name", default="O3",
                    help="Oxygen atom name to check (e.g., O3 or OH). Default: O3")
    ap.add_argument("--ideal_oh", type=float, default=0.96,
                    help="Ideal O–H bond length for protonated variant (Å). Default: 0.96")
    ap.add_argument("--tol", type=float, default=0.12,
                    help="Tolerance around ideal O–H (Å) for PASS. Default: 0.12")
    args = ap.parse_args()

    if not args.pdb_dir.is_dir():
        print(f"Directory not found: {args.pdb_dir}")
        return

    files = sorted(args.pdb_dir.glob("*.pdb"))
    if not files:
        print(f"No PDBs found in {args.pdb_dir}")
        return

    for pdb_path in files:
        status, details = analyze_one(pdb_path, args.oxygen_name, args.ideal_oh, args.tol)
        print(f"[{status}] {pdb_path.name} — {details}")

if __name__ == "__main__":
    main()
