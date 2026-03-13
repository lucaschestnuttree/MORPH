import math
import random
from pymol import cmd
from string import ascii_uppercase

def custom_cba(selection):
    """
    Custom atomic coloring logic scoped to a specific selection.
    Colors atoms by type within the selection without affecting other atoms.
    """
    atom_colors = {"N": "blue", "O": "red", "S": "orange", "H": "white"}
    for atom_type, color in atom_colors.items():
        cmd.color(color, f"{selection} and elem {atom_type}")

def apply_unique_translation(selection: str, r: float, used_directions: set):
    """
    Apply a unique translation of magnitude r in a random 3D direction to all atoms in a selection.
    Ensures the translation does not match previous ones.
    """
    max_attempts = 10
    for _ in range(max_attempts):
        theta = random.uniform(0, 2 * math.pi)
        phi = random.uniform(0, math.pi)

        dx = r * math.sin(phi) * math.cos(theta)
        dy = r * math.sin(phi) * math.sin(theta)
        dz = r * math.cos(phi)

        direction = (round(dx, 2), round(dy, 2), round(dz, 2))
        if direction not in used_directions:
            used_directions.add(direction)
            cmd.alter_state(1, selection, f"(x, y, z) = (x + {dx}, y + {dy}, z + {dz})")
            print(f"Translated '{selection}' by: dx={dx:.2f}, dy={dy:.2f}, dz={dz:.2f}")
            return

    print(f"Warning: Could not find a unique direction for '{selection}' after {max_attempts} attempts.")

def build_mismatch_dict(obj: str,
                        variant_chains=("V", "W", "X", "Y"),
                        use_ca_only=True,
                        b_min=None,
                        b_max=None):
    """
    dict[resi][chain] = True/False indicating mismatch vs chain A at same resi.
    If b_min/b_max given, restrict comparison to CA atoms with b in (b_min, b_max).
    """
    atom_filter = "name CA" if use_ca_only else "polymer.protein"

    b_filter = ""
    if b_min is not None and b_max is not None:
        b_filter = f" and (b > {b_min} and b < {b_max})"

    a_model = cmd.get_model(f"({obj}) and chain A and {atom_filter}{b_filter}")
    a_resn = {at.resi: at.resn for at in a_model.atom}

    mismatch = {}
    for ch in variant_chains:
        if cmd.count_atoms(f"({obj}) and chain {ch}") == 0:
            continue

        v_model = cmd.get_model(f"({obj}) and chain {ch} and {atom_filter}{b_filter}")
        for at in v_model.atom:
            resi = at.resi
            mismatch.setdefault(resi, {})
            a_res = a_resn.get(resi)
            mismatch[resi][ch] = (a_res is None) or (at.resn != a_res)

    return mismatch

def resi_clause(mismatch: dict, chain_letter: str):
    """Build a PyMOL `resi` clause (e.g., '12+45+103') for mismatches in the given chain."""
    resis = [r for r, d in mismatch.items() if d.get(chain_letter)]
    resis = sorted(set(resis), key=lambda x: (len(x), x))
    return "+".join(resis) if resis else None

def color_gnomaf(*args, **kwargs):
    """
    Colors structures based on gnomAD allele frequency bins.
    threshold=1 extends coloring to a wider range.
    """
    selection = kwargs.get("sel")
    if not selection:
        print("Error: No selection (PDB object) provided. Use sel='your_object'.")
        return

    threshold = int(args[0]) if len(args) >= 1 else 0

    # Colors
    cmd.set_color(f"{selection}_bin_0_color", [224, 224, 224]) # Frequency: [1e-9, 1e-8) # Grey Shades
    cmd.set_color(f"{selection}_bin_1_color", [192, 192, 192]) # Frequency: [1e-8, 1e-7) # Grey Shades
    cmd.set_color(f"{selection}_bin_2_color", [160, 160, 160]) # Frequency: [1e-7, 1e-6) # Grey Shades
    cmd.set_color(f"{selection}_bin_3_color", [128, 128, 128]) # Frequency: [1e-6, 1e-5) # Grey Shades
    cmd.set_color(f"{selection}_bin_5_color", [255, 51, 51]) # Frequency: [1e-4, 1e-3) # Red
    cmd.set_color(f"{selection}_bin_6_color", [153, 0, 0]) # Frequency: [1e-3, 1e-2) # Dark Red
    cmd.set_color(f"{selection}_bin_7_color", [51, 51, 255]) # Frequency: [1e-2, 0.51) # Blue
    cmd.set_color(f"{selection}_bin_8_color", [153, 51, 255]) # Frequency: [0.51, 1.0]  # Purple
    cmd.set_color(f"{selection}_bin_paleyellow_color", [255, 255, 204]) # Yellow for absent/0 allele frequency

    cmd.set_color(f"{selection}_light_orange", [255, 204, 153])     # V (min_residue)
    cmd.set_color(f"{selection}_lighter_orange", [255, 229, 204])   # W (min_residue2)
    cmd.set_color(f"{selection}_lightest_orange", [255, 240, 224])  # X (min_residue3)

    # Bin window
    if threshold == 1:
        bin_min, bin_max = (3, 9) # bins 4, 5, 6, 7, 8 [1e-5, 1.0] ~20 seq
        cmd.set_color(f"{selection}_bin_4_color", [255, 204, 204]) # bin 4 [1e-5, 1e-4) ~
    else:
        bin_min, bin_max = (4, 9)  # bins 5, 6, 7, 8  [1e-4, 1.0] ~ 200 seq
        cmd.set_color(f"{selection}_bin_4_color", [96, 96, 96]) # bin 4  [1e-5, 1e-4) 

    # KEEP chain A selections as-is
    cmd.select(f"{selection}_variants",
               f"({selection}) and chain A and (name CA) and (b > {bin_min} and b < {bin_max})")
    cmd.select(f"{selection}_maj_residue",
               f"({selection}) and chain A and (b > {bin_min} and b < {bin_max}) and (sidechain or name CA)")

    # Build mismatch dict restricted by SAME b window (CA-only)
    mismatch = build_mismatch_dict(selection,
                                   variant_chains=("V", "W", "X"),
                                   use_ca_only=True,
                                   b_min=bin_min,
                                   b_max=bin_max)

    resi_V = resi_clause(mismatch, "V")
    resi_W = resi_clause(mismatch, "W")
    resi_X = resi_clause(mismatch, "X")

    # Variant selections: mismatches only (still within your b window)
    cmd.select(f"{selection}_min_residue",
               f"({selection}) and chain V and resi {resi_V} and (b > {bin_min} and b < {bin_max}) and (sidechain or name CA)"
               if resi_V else "none")

    cmd.select(f"{selection}_min_residue2",
               f"({selection}) and chain W and resi {resi_W} and (b > {bin_min} and b < {bin_max}) and (sidechain or name CA)"
               if resi_W else "none")

    cmd.select(f"{selection}_min_residue3",
               f"({selection}) and chain X and resi {resi_X} and (b > {bin_min} and b < {bin_max}) and (sidechain or name CA)"
               if resi_X else "none")

    # Final display settings
    cmd.bg_color("white")
    cmd.set("fog", 0)
    cmd.hide("everything", f"({selection}) and (chain V or chain W or chain X)")
    cmd.set("orthoscopic", 1)
    cmd.set("specular", 0)

    # Main structure colored by allele-frequency bins 
    cmd.show("cartoon", f"({selection}) and chain A")
    for i in range(9):
        cmd.color(f"{selection}_bin_{i}_color", f"({selection}) and chain A and b={i}")
    cmd.color(f"{selection}_bin_paleyellow_color", f"({selection}) and chain A and b=42")

    # Variant spheres
    cmd.create(f"{selection}_variant_spheres", f"{selection}_variants")
    cmd.alter(f"{selection}_variant_spheres", "chain = 'Z'")
    cmd.alter(f"{selection}_variant_spheres", "name = 'O'")
    cmd.alter(f"{selection}_variant_spheres", "resn = 'HOH'")
    cmd.delete(f"{selection}_variants")
    cmd.show("spheres", f"{selection}_variant_spheres")
    cmd.set("sphere_scale", 0.7)

    # Color variant chains
    cmd.color(f"{selection}_light_orange", f"({selection}) and chain V")
    cmd.color(f"{selection}_lighter_orange", f"({selection}) and chain W")
    cmd.color(f"{selection}_lightest_orange", f"({selection}) and chain X")

    # Sidechains coloring
    custom_cba(f"({selection}_maj_residue) and sidechain")
    custom_cba(f"{selection}_min_residue")
    custom_cba(f"{selection}_min_residue2")
    custom_cba(f"{selection}_min_residue3")
    cmd.color("brightorange", f"({selection}_maj_residue) and sidechain and elem C")

def color_gnomcat(*args, **kwargs):
    selection = kwargs.get("sel")
    if not selection:
        print("Error: No selection (PDB object) provided. Use sel='your_object'.")
        return

    # colors
    cmd.set_color(f"{selection}_reder", [255, 102, 102])
    cmd.set_color(f"{selection}_blueer", [102, 102, 255])
    cmd.set_color(f"{selection}_violeter", [178, 102, 255])
    cmd.set_color(f"{selection}_oranger", [255, 178, 102])
    cmd.set_color(f"{selection}_yellower", [255, 255, 102])  

    cmd.set_color(f"{selection}_redish", [255, 204, 204])
    cmd.set_color(f"{selection}_blueish", [204, 204, 255])
    cmd.set_color(f"{selection}_violetish", [229, 204, 255])
    cmd.set_color(f"{selection}_orangish", [255, 229, 204])
    cmd.set_color(f"{selection}_yellowish", [255, 255, 204])  
    cmd.set_color(f"{selection}_light_grey", [224, 224, 224])

    threshold = int(args[0]) if len(args) >= 1 else 0
    if threshold == 1:
        bin_min, bin_max = (4, 9)  # bins 5, 6, 7 ,8: [1e-4, 1.0] ~10 seq
        bin_min2, bin_max2 = (4, 6) # bin 5 [1e-4, 1e-3) 
    else:
        bin_min, bin_max = (5, 9) # bins 6, 7, 8: [1e-3, 1.0] ~100 seq
        bin_min2, bin_max2 = (0, 1) # no bin

    variant_chains = ["V", "W", "X", "Y", "Z"]
    variant_chains = [c for c in variant_chains if cmd.count_atoms(f"({selection}) and chain {c}") > 0]
    if not variant_chains:
        print("No variant chains found with atoms. Nothing to process.")
        return

    atom = "CA"
    colors  = [f"{selection}_yellower", f"{selection}_oranger", f"{selection}_reder", f"{selection}_violeter", f"{selection}_blueer"]
    colors2 = [f"{selection}_yellowish", f"{selection}_orangish", f"{selection}_redish", f"{selection}_violetish", f"{selection}_blueish"]

    used_directions = set()

    cmd.select(f"{selection}_chain_A", f"({selection}) and chain A")
    cmd.show("cartoon", f"{selection}_chain_A")
    cmd.color(f"{selection}_light_grey", f"{selection}_chain_A")
    cmd.delete(f"{selection}_chain_A")

    cmd.select(f"{selection}_maj_residue",
               f"({selection}) and chain A and (b > {bin_min} and b < {bin_max})")
    custom_cba(f"{selection}_maj_residue")

    for i, chain in enumerate(variant_chains):
        color = colors[i % len(colors)]
        color2 = colors2[i % len(colors2)]

        tag = chain  # ← chain-only naming

        tmp_sel = f"{selection}_temp_variants_{tag}"
        cmd.select(tmp_sel,
                   f"({selection}) and chain {chain} and name {atom} "
                   f"and (b > {bin_min} and b < {bin_max})")
        cmd.color(color, tmp_sel)
        cmd.color(color2, f"{tmp_sel} and (b > {bin_min2} and b < {bin_max2})")

        var_obj = f"{selection}_{tag}_variations"
        cmd.create(var_obj, tmp_sel)
        apply_unique_translation(var_obj, r=0.5, used_directions=used_directions)

        min_sel = f"{selection}_{tag}_min_residue"
        cmd.select(min_sel,
                   f"({selection}) and chain {chain} "
                   f"and (b > {bin_min} and b < {bin_max}) "
                   f"and (sidechain or name CA)")
        cmd.color(color, min_sel)
        cmd.color(color2, f"{min_sel} and (b > {bin_min2} and b < {bin_max2})")

        cmd.delete(tmp_sel)

    # final viz
    cmd.bg_color("white")
    cmd.hide("everything", f"({selection}) and (chain V or chain W or chain X or chain Y or chain Z)")
    cmd.set("specular", 0)
    cmd.set("use_shaders", 0)

    for chain in variant_chains:
        tag = chain
        cmd.hide("everything", f"{selection}_{tag}_variations")
        cmd.show("spheres", f"{selection}_{tag}_variations")
        custom_cba(f"{selection}_{tag}_min_residue")

    cmd.set("sphere_scale", 0.7)
    cmd.hide("sticks", f"({selection}) and (chain A or chain V or chain W or chain X or chain Y)")

def pos_var(*args, **kwargs):
    """
    Per-chain selections at a residue position, showing sidechain + CA.

    Usage:
      pos_var sel=myObj, 123        # m=1 default: only mismatching non-A chains
      pos_var sel=myObj, 123, 0     # m=0: show all chains

    Creates selections like:
      myObj_pos_123_A_sc
      myObj_pos_123_V_sc
      ...
    """
    sel = kwargs.get("sel")
    if sel is None:
        print("Error: sel=<object> is required (e.g., pos_var sel=myObj, 123)")
        return

    if len(args) not in (1, 2):
        print("Error: usage is pos_var sel=myObj, <resi> [, m]")
        return

    resi = str(args[0])
    m = int(args[1]) if len(args) == 2 else 1  # default mismatch-only for non-A

    # discover chains
    model = cmd.get_model(f"({sel}) and name CA")
    chains = sorted({a.chain for a in model.atom if a.chain})
    if not chains:
        model = cmd.get_model(f"({sel})")
        chains = sorted({a.chain for a in model.atom if a.chain})

    if not chains:
        print(f"No chains found in '{sel}'.")
        return

    # mismatch dict (only needed if m==1)
    mismatch = {}
    if m == 1:
        variant_chains = [c for c in chains if c != "A"]
        mismatch = build_mismatch_dict(sel, variant_chains=tuple(variant_chains), use_ca_only=True)

    made = []
    skipped_match = []
    skipped_empty = []
    skipped_missing = []

    for ch in chains:
        # residue present?
        if cmd.count_atoms(f"({sel}) and chain {ch} and resi {resi} and name CA") == 0:
            skipped_missing.append(ch)
            continue

        # mismatch filtering for non-A
        if m == 1 and ch != "A":
            if not mismatch.get(resi, {}).get(ch, False):
                skipped_match.append(ch)
                continue

        # selection includes sidechain + CA
        pos_sel = f"{sel}_pos_{resi}_{ch}"
        cmd.select(pos_sel, f"({sel}) and chain {ch} and resi {resi} and (sidechain or name CA)")

        if cmd.count_atoms(pos_sel) == 0:
            cmd.delete(pos_sel)
            skipped_empty.append(ch)
            continue

        # color chain A carbon atoms brightorange
        if ch == "A":
            cmd.color("brightorange", f"{pos_sel} and (sidechain)")

        # custom element coloring (N/O/S/H etc) for everything displayed
        custom_cba(pos_sel)

        cmd.show("sticks", f"{sel}_pos_{resi}_A and (sidechain or name CA)")

        made.append(ch)

    if not made:
        msg = f"pos_var: nothing displayed at resi {resi} in '{sel}' (m={m})."
        print(msg)
        if skipped_match:
            print("Skipped (matched chain A): " + ", ".join(skipped_match))
        if skipped_missing:
            print("Skipped (missing residue): " + ", ".join(skipped_missing))
        if skipped_empty:
            print("Skipped (empty selection): " + ", ".join(skipped_empty))
        return

    print(f"pos_var: displayed resi {resi} for chains: {', '.join(made)} (m={m})")
    if skipped_match:
        print("Not shown (matched chain A): " + ", ".join(skipped_match))

# Register the functions with PyMOL
cmd.extend("color_gnomaf", color_gnomaf)
cmd.extend("color_gnomcat", color_gnomcat)
cmd.extend("pos_var", pos_var)