import argparse
from pathlib import Path
from math import sqrt
from Bio.PDB import PDBParser, PDBIO, StructureBuilder, Atom

ELEMENT2 = {
    "CL":"Cl","BR":"Br","ZN":"Zn","MG":"Mg","NA":"Na","CA":"Ca","FE":"Fe",
    "CU":"Cu","NI":"Ni","MN":"Mn","SI":"Si","SE":"Se","CO":"Co","AL":"Al",
    "PT":"Pt","PD":"Pd","HG":"Hg","AG":"Ag","AU":"Au","PB":"Pb","SN":"Sn",
    "TI":"Ti","CR":"Cr","LI":"Li","CS":"Cs","RB":"Rb","ZR":"Zr","MO":"Mo",
    "RU":"Ru","RH":"Rh","IR":"Ir","OS":"Os","RE":"Re","TA":"Ta","GA":"Ga",
    "GE":"Ge","AS":"As","SR":"Sr","BA":"Ba","KR":"Kr","XE":"Xe","AR":"Ar",
    "HE":"He","NE":"Ne","RN":"Rn","CD":"Cd","CO":"Co","VA":"Va"
}

def dist(a, b):
    ax, ay, az = a.get_coord()
    bx, by, bz = b.get_coord()
    return sqrt((ax-bx)**2 + (ay-by)**2 + (az-bz)**2)

def unit_vec(from_xyz, to_xyz):
    vx = to_xyz[0] - from_xyz[0]
    vy = to_xyz[1] - from_xyz[1]
    vz = to_xyz[2] - from_xyz[2]
    n = sqrt(vx*vx + vy*vy + vz*vz)
    if n == 0:
        return (0.0, 0.0, 1.0)
    return (vx/n, vy/n, vz/n)

def element_of(atom):
    if hasattr(atom, "element") and atom.element:
        e = atom.element.strip()
        if len(e) == 2:
            return ELEMENT2.get(e.upper(), e.capitalize())
        return e[0].upper()

    fullname = atom.get_fullname()
    letters_only = "".join(ch for ch in fullname if ch.isalpha())
    if len(letters_only) >= 2:
        cand2 = letters_only[:2].upper()
        if cand2 in ELEMENT2:
            return ELEMENT2[cand2]
    if letters_only:
        return letters_only[0].upper()
    return "X"

def is_hydrogen(atom):
    return element_of(atom) == "H"

def heavy_neighbors(all_atoms, center, cutoff=1.8):
    return [a for a in all_atoms
            if a is not center and element_of(a) != "H" and dist(a, center) < cutoff]

def h_neighbors(all_atoms, center, cutoff=1.2):
    return [a for a in all_atoms
            if a is not center and is_hydrogen(a) and dist(a, center) < cutoff]

def build_single_residue_structure(chain_id, residue):
    sb = StructureBuilder.StructureBuilder()
    sb.init_structure("chrom")
    sb.init_model(0)
    sb.init_chain(chain_id)
    
    res_id = residue.get_id()
    sb.init_seg("    ")
    sb.init_residue(residue.get_resname(), res_id[0], res_id[1], res_id[2])
    for atom in residue.get_atoms():
        a = Atom.Atom(
            name=atom.get_name(),
            coord=atom.get_coord(),
            bfactor=atom.get_bfactor(),
            occupancy=atom.get_occupancy() if atom.get_occupancy() is not None else 1.0,
            altloc=atom.get_altloc() if hasattr(atom, "get_altloc") else ' ',
            fullname=atom.get_fullname(),
            serial_number=atom.get_serial_number() if hasattr(atom, "get_serial_number") else 0,
            element=element_of(atom)
        )
        sb.structure[0][chain_id].child_list[0].add(a)
    return sb.get_structure()

def find_phenolic_oxygen(residue, forced_name=None, debug=True):
    atoms = list(residue.get_atoms())
    if forced_name:
        for a in atoms:
            if a.get_name().strip() == forced_name.strip():
                if debug:
                    print(f"[DEBUG] Forced oxygen '{forced_name}' selected.")
                return a
        raise ValueError(f"Oxygen atom named '{forced_name}' not found in residue "
                         f"{residue.get_resname()} {residue.get_id()[1]}")

    candidates = []
    for a in atoms:
        if element_of(a) != "O":
            continue
        hn = len(h_neighbors(atoms, a))
        hnbrs = heavy_neighbors(atoms, a)
        if debug:
            print(f"[DEBUG] Oxygen {a.get_name()} > heavy neighbors: {len(hnbrs)}, H neighbors: {hn}")
        candidates.append((a, len(hnbrs), hn))

    candidates = sorted(candidates, key=lambda t: (abs(t[1]-1), t[2]))
    if not candidates:
        raise ValueError("No oxygen atoms found in residue.")
    chosen = candidates[0][0]
    if debug:
        print(f"[DEBUG] Chosen oxygen: {chosen.get_name()} in residue {residue.get_resname()} {residue.get_id()[1]}")
    return chosen

def add_h_to_oxygen(structure, chain_id, resid_tuple, o_atom, bond_ref_atom=None, oh_length=0.96, debug=False):
    chain = structure[0][chain_id]
    residue = list(chain.get_residues())[0]         
    atoms = list(residue.get_atoms())               

    if h_neighbors(atoms, o_atom):
        if debug:
            print(f"[DEBUG] Oxygen {o_atom.get_name()} already has hydrogen — skipping.")
        return

    if bond_ref_atom is None:
        hnbrs = heavy_neighbors(atoms, o_atom)
        if not hnbrs:
            candidates = [a for a in atoms if element_of(a) != "H" and a is not o_atom]
            if not candidates:
                raise ValueError("No heavy atom found to orient O–H.")
            bond_ref_atom = min(candidates, key=lambda a: dist(a, o_atom))
            if debug:
                print(f"[DEBUG] No close heavy neighbor; fallback nearest is "
                      f"{bond_ref_atom.get_name()} at {dist(bond_ref_atom, o_atom):.2f} A")
        else:
            bond_ref_atom = min(hnbrs, key=lambda a: dist(a, o_atom))
            if debug:
                print(f"[DEBUG] Using nearest heavy neighbor {bond_ref_atom.get_name()} "
                      f"at {dist(bond_ref_atom, o_atom):.2f} A")

    ox, oy, oz = o_atom.get_coord()
    rx, ry, rz = bond_ref_atom.get_coord()
    vx, vy, vz = unit_vec((rx, ry, rz), (ox, oy, oz))
    hx, hy, hz = ox + oh_length * vx, oy + oh_length * vy, oz + oh_length * vz

    base = "HO"
    existing_names = {a.get_name().strip() for a in atoms}
    name = base
    idx = 1
    while name in existing_names:
        name = f"HO{idx}"
        idx += 1

    if debug:
        print(f"[DEBUG] Adding H '{name}' at ({hx:.3f}, {hy:.3f}, {hz:.3f}); O–H = {oh_length:.2f} A")

    h_atom = Atom.Atom(
        name=name,
        coord=(hx, hy, hz),       
        bfactor=1.0,
        occupancy=1.0,
        altloc=' ',
        fullname=f"{name:>4}",
        serial_number=0,
        element='H'
    )
    residue.add(h_atom)

def remove_o_hydrogens(structure, chain_id, resid_tuple, o_atom, debug=True):
    chain = structure[0][chain_id]
    residue = list(chain.get_residues())[0]
    to_remove = []

    if debug:
        print(f"\n[DEBUG] Checking hydrogens near oxygen {o_atom.get_name()} in chain {chain_id}...")

    for a in list(residue.get_atoms()):
        if is_hydrogen(a):
            d = dist(a, o_atom)
            if debug:
                print(f"  Found H atom: {a.get_name()} at distance {d:.3f} A from {o_atom.get_name()}")
            if d < 1.25:
                if debug:
                    print(f"    > Marked for removal (bonded to O)")
                to_remove.append(a.get_name())

    for name in to_remove:
        residue.detach_child(name)
        if debug:
            print(f"  Removed hydrogen: {name}")

    if debug and not to_remove:
        print("  No hydrogens removed — oxygen was already deprotonated.")

def write_pdb(structure, outpath):
    io = PDBIO()
    io.set_structure(structure)
    io.save(str(outpath))

def main():
    ap = argparse.ArgumentParser(description="Extract chromophore and write protonated/deprotonated variants.")
    ap.add_argument("pdb", type=Path, help="Input PDB (e.g., 2vzx.pdb)")
    ap.add_argument("--resname", default="5SQ", help="Chromophore residue name (default: 5SQ)")
    ap.add_argument("--resid", type=int, default=64, help="Residue number (default: 64)")
    ap.add_argument("--oxygen-name", default=None, help="Force oxygen atom name (e.g., OZ or OH).")
    ap.add_argument("--chains", default=None,
                    help="Comma-separated chain IDs (e.g., A,B,C). Default: all chains containing that residue.")
    ap.add_argument("--outdir", type=Path, default=Path("chromophore_variants"), help="Output directory")
    args = ap.parse_args()

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure("input", str(args.pdb))

    targets = []
    chain_filter = set([c.strip() for c in args.chains.split(",")]) if args.chains else None
    for model in structure:
        for chain in model:
            cid = chain.get_id()
            if chain_filter and cid not in chain_filter:
                continue
            for res in chain:
                het, rnum, icode = res.get_id()
                if res.get_resname().strip() == args.resname and rnum == args.resid:
                    targets.append((cid, res))

    if not targets:
        raise SystemExit(f"No matches for {args.resname} {args.resid} in the provided PDB.")

    args.outdir.mkdir(parents=True, exist_ok=True)
    written = []

    for cid, residue in targets:
        base_struct = build_single_residue_structure(cid, residue)

        oxy = find_phenolic_oxygen(list(base_struct[0][cid].get_residues())[0],
                                   forced_name=args.oxygen_name, debug=True)

        deprot = build_single_residue_structure(cid, residue)
        oxy_d = find_phenolic_oxygen(list(deprot[0][cid].get_residues())[0],
                                     forced_name=args.oxygen_name, debug=False)
        remove_o_hydrogens(deprot, cid, None, oxy_d, debug=True)
        out_d = args.outdir / f"chrom_{cid}_{args.resname}{args.resid}_deprot.pdb"
        write_pdb(deprot, out_d)
        written.append(out_d)

        prot = build_single_residue_structure(cid, residue)
        oxy_p = find_phenolic_oxygen(list(prot[0][cid].get_residues())[0],
                                     forced_name=args.oxygen_name, debug=False)
        add_h_to_oxygen(prot, cid, None, oxy_p, debug=True)
        out_p = args.outdir / f"chrom_{cid}_{args.resname}{args.resid}_prot.pdb"
        write_pdb(prot, out_p)
        written.append(out_p)

    print("Wrote:")
    for p in written:
        print("  ", p)


if __name__ == "__main__":
    main()
