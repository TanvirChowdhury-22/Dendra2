import argparse
from pathlib import Path
from math import sqrt, cos, sin, radians
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

def add_h_to_oxygen(structure, chain_id, resid_tuple, o_atom, bond_ref_atom=None,
                    oh_length=0.96, target_angle_deg=109.0, debug=False):
    chain = structure[0][chain_id]
    residue = list(chain.get_residues())[0]
    atoms = list(residue.get_atoms())

    if any(a.element == 'H' and dist(a, o_atom) < 1.25 for a in atoms):
        if debug:
            print(f"[DEBUG] {o_atom.get_name()} already has an H nearby; skipping.")
        return

    if bond_ref_atom is None:
        heavies = [a for a in atoms if a is not o_atom and element_of(a) != "H"]
        if not heavies:
            raise ValueError("No heavy atom found to orient O–H.")
        bond_ref_atom = min(heavies, key=lambda a: dist(a, o_atom))
        if debug:
            print(f"[DEBUG] Using nearest heavy neighbor {bond_ref_atom.get_name()} "
                  f"at {dist(bond_ref_atom, o_atom):.2f} Å")

<<<<<<< HEAD
    ox, oy, oz = o_atom.get_coord()
    rx, ry, rz = bond_ref_atom.get_coord()
    vx, vy, vz = unit_vec((rx, ry, rz), (ox, oy, oz))
    hx, hy, hz = ox + oh_length * vx, oy + oh_length * vy, oz + oh_length * vz
=======
    O = o_atom.get_coord()
    C = bond_ref_atom.get_coord()

    u_OC = unit_vec(O, C)                
    u_CO = (-u_OC[0], -u_OC[1], -u_OC[2])

    ring_neighbors = [a for a in atoms
                      if a is not o_atom and a is not bond_ref_atom and element_of(a) != "H"]
    by_name = {a.get_name().strip(): a for a in ring_neighbors}
    picks = []
    for nm in ("CE1", "CE2", "CD1", "CD2"):
        if nm in by_name:
            picks.append(by_name[nm])
    if len(picks) < 2:
        for a in sorted(ring_neighbors, key=lambda x: dist(x, bond_ref_atom)):
            if a not in picks:
                picks.append(a)
            if len(picks) == 2:
                break

    def norm(v):
        return (v[0]*v[0] + v[1]*v[1] + v[2]*v[2])**0.5

    if len(picks) >= 2:
        v1 = picks[0].get_coord() - C
        v2 = picks[1].get_coord() - C
        # ring normal n = v1 x v2
        nx = v1[1]*v2[2] - v1[2]*v2[1]
        ny = v1[2]*v2[0] - v1[0]*v2[2]
        nz = v1[0]*v2[1] - v1[1]*v2[0]
        nlen = norm((nx, ny, nz))
        if nlen > 1e-6:
            n = (nx/nlen, ny/nlen, nz/nlen)
            # t = n x u_CO (lies in ring plane, ⟂ to O–C)
            tx = n[1]*u_CO[2] - n[2]*u_CO[1]
            ty = n[2]*u_CO[0] - n[0]*u_CO[2]
            tz = n[0]*u_CO[1] - n[1]*u_CO[0]
            tlen = norm((tx, ty, tz))
            t = (tx/tlen, ty/tlen, tz/tlen) if tlen > 1e-6 else (0.0, 0.0, 1.0)
        else:
            t = (0.0, 0.0, 1.0)
    else:
        zx, zy, zz = 0.0, 0.0, 1.0
        dot = u_CO[2]
        px, py, pz = zx - dot*u_CO[0], zy - dot*u_CO[1], zz - dot*u_CO[2]
        plen = norm((px, py, pz))
        t = (1.0, 0.0, 0.0) if plen < 1e-6 else (px/plen, py/plen, pz/plen)

    theta = radians(180.0 - target_angle_deg)
    d = oh_length
    hx = O[0] + d*(cos(theta)*u_CO[0] + sin(theta)*t[0])
    hy = O[1] + d*(cos(theta)*u_CO[1] + sin(theta)*t[1])
    hz = O[2] + d*(cos(theta)*u_CO[2] + sin(theta)*t[2])
>>>>>>> 96af0ea (prot-deprot site correction)

    base = "HO"
    existing = {a.get_name().strip() for a in atoms}
    name = base; k = 1
    while name in existing:
        name = f"HO{k}"; k += 1

    if debug:
        print(f"[DEBUG] Adding H '{name}' at ({hx:.3f}, {hy:.3f}, {hz:.3f}); "
              f"O–H = {oh_length:.2f} Å, ∠(C–O–H) ≈ {target_angle_deg:g}°")

    h_atom = Atom.Atom(
        name=name, coord=(hx, hy, hz),
        bfactor=1.0, occupancy=1.0, altloc=' ',
        fullname=f"{name:>4}", serial_number=0, element='H'
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
