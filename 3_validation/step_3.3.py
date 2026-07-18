"""
phaC_structural_validation.py

Consolidated structural validation of the candidate PHA synthase
(phaC; locus tag MPDKNC_03777) from a ColabFold/AlphaFold2 prediction.

This script performs every analysis step carried out after the ColabFold
result zip was produced:

    1. Unzip the ColabFold result archive and locate the top-ranked
       (rank_001) model.
    2. Sanity-check the structure: confirm sequence identity and residue
       numbering against the original query (no His-tag / offset issues).
    3. Per-residue pLDDT confidence analysis (global + lipase-box region).
    4. Catalytic triad geometric verification:
           a. Search for His/Asp side-chain atoms near the Cys133 (S-gamma)
              nucleophile (lipase-box catalytic residue).
           b. Secondary search for Asp hydrogen-bonding to any candidate
              His (the true Cys-His-Asp relay geometry), independent of
              that Asp's raw distance to Cys.
    5. Render two annotated figures with PyMOL:
           a. Overall fold colored by per-residue pLDDT.
           b. Active-site close-up of the Cys133-His318-Asp289 triad.

Requirements:
    biopython, pymol-open-source
    (pip install biopython pymol-open-source --break-system-packages)

Usage:
    python3 phaC_structural_validation.py <path_to_result_zip> [output_dir]

Outputs (written to <output_dir>, default "./phaC_structural_results"):
    extracted/                              (unzipped ColabFold archive)
    phaC_overall_plddt.png                  (Figure 1: overall confidence)
    phaC_active_site_triad.png              (Figure 2: active-site close-up)
    phaC_structural_validation_summary.txt  (all numerical results)
"""

import os
import sys
import glob
import zipfile

import numpy as np
from Bio.PDB import PDBParser, MMCIFParser, is_aa

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------

QUERY_SEQUENCE = (
    "MERKTEKDRWMGFFETLTRREEWKPHHPRNAIQEIGRATLWHYPTVESNGALPILMVYSHINKPSILDL"
    "TEQHSMIGEFLRNGYDVFLLDFGIPDERDKDTGLESYLFDYIDGAVKTVLQHQQTGRLTLAGFCLGGTL"
    "AALYAALYPDRIHNLLLFVTPIDFNQLPDFREWIKAIQTGAIDPAVIAPATGIIPAEQIRYGMRLITAP"
    "VYYSPYLSLLHRSHDPAYTEHWYRFNQWTNDHIPMTGEFLRDLLHYFIKQNALMNGGMTLRGREINPAR"
    "IQSNVYMVCSKFDQMVPAAISYPLMELVSSEEKQFIEVPGGHASLIKDGVSSRLMDWLREHS"
)

CYS_NUCLEOPHILE_RESSEQ = 133   # from sequence-based lipase-box motif GxCxG at 131-135
LIPASE_BOX_RANGE = (126, 140)  # padding around residues 129-135
TRIAD_SEARCH_RADIUS = 7.0      # Angstroms: Cys(SG) -> His/Asp search
HIS_ASP_SEARCH_RADIUS = 6.0    # Angstroms: candidate His -> Asp search

THREE_TO_ONE = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D", "CYS": "C", "GLN": "Q",
    "GLU": "E", "GLY": "G", "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S", "THR": "T", "TRP": "W",
    "TYR": "Y", "VAL": "V",
}


# --------------------------------------------------------------------------
# STEP 1: Unzip and locate the top-ranked model
# --------------------------------------------------------------------------

def extract_and_find_top_model(zip_path, out_dir):
    extract_dir = os.path.join(out_dir, "extracted")
    os.makedirs(extract_dir, exist_ok=True)

    print(f"Extracting {zip_path} -> {extract_dir}")
    with zipfile.ZipFile(zip_path, "r") as zf:
        zf.extractall(extract_dir)

    candidates = glob.glob(
        os.path.join(extract_dir, "**", "*rank_001*.pdb"), recursive=True
    )
    if not candidates:
        candidates = glob.glob(
            os.path.join(extract_dir, "**", "*rank_001*.cif"), recursive=True
        )
    if not candidates:
        raise FileNotFoundError(
            f"Could not find a rank_001 model (.pdb or .cif) under {extract_dir}. "
            f"Check the archive contents."
        )

    # Prefer a relaxed model if one exists among the rank_001 matches
    relaxed = [c for c in candidates if "relaxed" in os.path.basename(c) and "unrelaxed" not in os.path.basename(c)]
    top_model_path = relaxed[0] if relaxed else candidates[0]

    print(f"Top-ranked model: {top_model_path}")
    if not relaxed:
        print("  NOTE: no amber-relaxed rank_001 model found -- using the "
              "unrelaxed model. Minor local geometry (bond angles/clashes) "
              "may differ slightly from a relaxed structure, though this "
              "rarely affects overall fold or inter-residue distances "
              "meaningfully.")
    return top_model_path


# --------------------------------------------------------------------------
# STEP 2: Load structure and sanity-check against the query sequence
# --------------------------------------------------------------------------

def load_structure(path):
    parser = MMCIFParser(QUIET=True) if path.lower().endswith(".cif") else PDBParser(QUIET=True)
    return parser.get_structure("phaC", path)


def sanity_check_sequence(structure, query_sequence):
    model = structure[0]
    chain = next(iter(model))

    residues = [res for res in chain if is_aa(res, standard=True)]
    resseqs = [res.id[1] for res in residues]
    pdb_seq = "".join(THREE_TO_ONE.get(res.get_resname(), "X") for res in residues)

    print("\n--- Sequence / numbering sanity check ---")
    print(f"  Structure length: {len(pdb_seq)} residues "
          f"(range {min(resseqs)}-{max(resseqs)})")
    print(f"  Query length:     {len(query_sequence)} residues")

    match = pdb_seq == query_sequence
    print(f"  Sequence identical to query: {match}")
    if not match:
        # find first mismatch for diagnostic purposes
        n = min(len(pdb_seq), len(query_sequence))
        first_diff = next((i for i in range(n) if pdb_seq[i] != query_sequence[i]), None)
        print(f"  WARNING: sequence mismatch detected"
              + (f" at position {first_diff + 1}" if first_diff is not None else " (length mismatch)"))

    cys_positions = [res.id[1] for res in residues if res.get_resname() == "CYS"]
    print(f"  Cysteine residues in structure: {cys_positions}")
    if CYS_NUCLEOPHILE_RESSEQ not in cys_positions:
        print(f"  WARNING: expected catalytic Cys{CYS_NUCLEOPHILE_RESSEQ} not found "
              f"at that residue number -- numbering offset likely (e.g. tag/linker).")
    else:
        print(f"  Confirmed: Cys{CYS_NUCLEOPHILE_RESSEQ} numbering matches query sequence.")

    return match


# --------------------------------------------------------------------------
# STEP 3: Per-residue pLDDT confidence analysis
# --------------------------------------------------------------------------

def get_residue_plddt(structure):
    """ColabFold/AlphaFold2 store per-residue pLDDT in the B-factor column."""
    plddt = {}
    model = structure[0]
    for chain in model:
        for res in chain:
            if not is_aa(res, standard=True):
                continue
            bfactors = [atom.get_bfactor() for atom in res]
            plddt[res.id[1]] = float(np.mean(bfactors))
        break  # only first chain
    return plddt


def summarize_confidence(plddt, log):
    log("=" * 70)
    log("PER-REGION CONFIDENCE (mean pLDDT)")
    log("=" * 70)

    all_vals = np.array(list(plddt.values()))
    log(f"Overall mean pLDDT: {all_vals.mean():.1f}  (n={len(all_vals)} residues)")

    lo, hi = LIPASE_BOX_RANGE
    box_vals = [v for k, v in plddt.items() if lo <= k <= hi]
    if box_vals:
        log(f"Lipase-box region ({lo}-{hi}): mean pLDDT = {np.mean(box_vals):.1f}  "
            f"(min={min(box_vals):.1f}, max={max(box_vals):.1f})")

    cys_plddt = plddt.get(CYS_NUCLEOPHILE_RESSEQ)
    if cys_plddt is not None:
        log(f"Cys{CYS_NUCLEOPHILE_RESSEQ} (catalytic nucleophile): pLDDT = {cys_plddt:.1f}")

    very_high = int((all_vals >= 90).sum())
    confident = int(((all_vals >= 70) & (all_vals < 90)).sum())
    low = int(((all_vals >= 50) & (all_vals < 70)).sum())
    very_low = int((all_vals < 50).sum())
    n = len(all_vals)
    log("\nConfidence distribution:")
    log(f"  Very high (>=90): {very_high:4d} residues ({100*very_high/n:.1f}%)")
    log(f"  Confident (70-90): {confident:4d} residues ({100*confident/n:.1f}%)")
    log(f"  Low (50-70):       {low:4d} residues ({100*low/n:.1f}%)")
    log(f"  Very low (<50):    {very_low:4d} residues ({100*very_low/n:.1f}%)")


# --------------------------------------------------------------------------
# STEP 4: Catalytic triad geometric verification
# --------------------------------------------------------------------------

def find_triad_candidates(structure, cys_resseq=CYS_NUCLEOPHILE_RESSEQ, radius=TRIAD_SEARCH_RADIUS):
    """His/Asp side-chain atoms within `radius` of Cys(SG)."""
    model = structure[0]
    chain = next(iter(model))

    cys_res = next((res for res in chain if res.id[1] == cys_resseq and is_aa(res, standard=True)), None)
    if cys_res is None or "SG" not in cys_res:
        return None  # signal: not found

    sg_coord = cys_res["SG"].get_coord()

    candidates = []
    for res in chain:
        if not is_aa(res, standard=True):
            continue
        resname = res.get_resname()
        if resname not in ("HIS", "ASP"):
            continue
        atoms_to_check = [a for a in (("NE2", "ND1") if resname == "HIS" else ("OD1", "OD2")) if a in res]

        min_dist = None
        for atom_name in atoms_to_check:
            d = float(np.linalg.norm(res[atom_name].get_coord() - sg_coord))
            if min_dist is None or d < min_dist:
                min_dist = d

        if min_dist is not None and min_dist <= radius:
            candidates.append({"resname": resname, "resseq": res.id[1], "distance_to_Cys_SG": min_dist})

    candidates.sort(key=lambda c: c["distance_to_Cys_SG"])
    return candidates


def find_asp_near_his(structure, his_resseq, radius=HIS_ASP_SEARCH_RADIUS):
    """Asp side-chain atoms (OD1/OD2) within `radius` of a candidate His's
    side-chain N atoms (ND1/NE2) -- the true triad relay geometry."""
    model = structure[0]
    chain = next(iter(model))

    his_res = next((res for res in chain if res.id[1] == his_resseq and res.get_resname() == "HIS"), None)
    if his_res is None:
        return []

    his_atoms = [his_res[a].get_coord() for a in ("ND1", "NE2") if a in his_res]

    hits = []
    for res in chain:
        if not is_aa(res, standard=True) or res.get_resname() != "ASP":
            continue
        asp_atoms = [res[a].get_coord() for a in ("OD1", "OD2") if a in res]
        min_dist = None
        for ha in his_atoms:
            for aa in asp_atoms:
                d = float(np.linalg.norm(ha - aa))
                if min_dist is None or d < min_dist:
                    min_dist = d
        if min_dist is not None and min_dist <= radius:
            hits.append((res.id[1], min_dist))

    hits.sort(key=lambda x: x[1])
    return hits


def run_triad_analysis(structure, plddt, log):
    log("\n" + "=" * 70)
    log(f"CATALYTIC TRIAD SEARCH (partners within {TRIAD_SEARCH_RADIUS} A of "
        f"Cys{CYS_NUCLEOPHILE_RESSEQ} SG)")
    log("=" * 70)

    candidates = find_triad_candidates(structure)
    if candidates is None:
        log(f"WARNING: Could not find Cys{CYS_NUCLEOPHILE_RESSEQ} SG atom in structure. "
            f"Check residue numbering matches the query sequence.")
        return None

    if not candidates:
        log("No His/Asp side chains found within range. Consider increasing "
            "TRIAD_SEARCH_RADIUS or visually inspecting the structure.")
        return None

    for c in candidates:
        log(f"  {c['resname']}{c['resseq']:>4d}   distance to Cys{CYS_NUCLEOPHILE_RESSEQ} SG = "
            f"{c['distance_to_Cys_SG']:.2f} A")

    his_hits = [c for c in candidates if c["resname"] == "HIS"]
    asp_hits_direct = [c for c in candidates if c["resname"] == "ASP"]

    result = {
        "cys_resseq": CYS_NUCLEOPHILE_RESSEQ,
        "his_candidates": his_hits,
        "asp_direct_candidates": asp_hits_direct,
        "triad": None,
    }

    if not his_hits:
        log("\nNo candidate His found near the nucleophile -- triad incomplete.")
        return result

    log("\n" + "=" * 70)
    log("SECONDARY SEARCH: Asp hydrogen-bonding to candidate His "
        "(classic triad geometry: Asp---His---Cys)")
    log("=" * 70)

    best_triad = None
    for h in his_hits:
        asp_near = find_asp_near_his(structure, h["resseq"])
        if asp_near:
            log(f"  His{h['resseq']} (Cys{CYS_NUCLEOPHILE_RESSEQ} distance "
                f"{h['distance_to_Cys_SG']:.2f} A) has {len(asp_near)} Asp candidate(s) "
                f"within {HIS_ASP_SEARCH_RADIUS} A:")
            for asp_id, d in asp_near:
                log(f"    Asp{asp_id}: {d:.2f} A from His{h['resseq']} side chain")
            if best_triad is None:
                best_triad = {
                    "cys_resseq": CYS_NUCLEOPHILE_RESSEQ,
                    "his_resseq": h["resseq"],
                    "asp_resseq": asp_near[0][0],
                    "cys_his_distance": h["distance_to_Cys_SG"],
                    "his_asp_distance": asp_near[0][1],
                }
        else:
            log(f"  His{h['resseq']}: no Asp found within {HIS_ASP_SEARCH_RADIUS} A "
                f"(may need visual inspection / wider radius)")

    result["triad"] = best_triad

    if best_triad:
        c, h, a = best_triad["cys_resseq"], best_triad["his_resseq"], best_triad["asp_resseq"]
        log(f"\n-> Candidate catalytic triad confirmed: "
            f"Cys{c} - His{h} - Asp{a}")
        log(f"   Cys{c}(SG)-His{h} distance: {best_triad['cys_his_distance']:.2f} A")
        log(f"   His{h}-Asp{a} distance:      {best_triad['his_asp_distance']:.2f} A")
        log(f"   pLDDT: Cys{c}={plddt.get(c, float('nan')):.1f}  "
            f"His{h}={plddt.get(h, float('nan')):.1f}  "
            f"Asp{a}={plddt.get(a, float('nan')):.1f}")
    else:
        log("\nNo complete triad found (His identified but no nearby Asp). "
            "Recommend visual inspection in PyMOL/ChimeraX.")

    return result


# --------------------------------------------------------------------------
# STEP 5: Render annotated figures with PyMOL
# --------------------------------------------------------------------------

def render_figures(pdb_path, triad, out_dir):
    try:
        import pymol
        from pymol import cmd
    except ImportError:
        print("\nNOTE: pymol-open-source not installed -- skipping figure "
              "rendering. Install with:\n"
              "  pip install pymol-open-source --break-system-packages")
        return None, None

    pymol.finish_launching(["pymol", "-qc"])

    overall_png = os.path.join(out_dir, "phaC_overall_plddt.png")
    triad_png = os.path.join(out_dir, "phaC_active_site_triad.png")

    cmd.load(pdb_path, "phaC")
    cmd.hide("everything")
    cmd.bg_color("white")
    cmd.set("ray_opaque_background", 0)
    cmd.set("antialias", 2)
    cmd.set("cartoon_fancy_helices", 1)
    cmd.set("cartoon_highlight_color", "grey70")

    # --- Figure 1: overall fold colored by pLDDT ---
    cmd.show("cartoon", "phaC")
    cmd.spectrum("b", "red_yellow_green", "phaC", minimum=50, maximum=100)
    cmd.orient("phaC")
    cmd.zoom("phaC", buffer=3)
    cmd.ray(1400, 1100)
    cmd.png(overall_png, dpi=300)
    print(f"Figure written: {overall_png}")

    # --- Figure 2: active-site close-up ---
    if triad is not None:
        c, h, a = triad["cys_resseq"], triad["his_resseq"], triad["asp_resseq"]
        cmd.color("grey80", "phaC")
        cmd.show("cartoon", "phaC")
        cmd.set("cartoon_transparency", 0.6)

        triad_sel = f"resi {c}+{h}+{a}"
        cmd.show("sticks", triad_sel + " and not name C+N+O")
        cmd.show("sticks", triad_sel + " and name CA")
        cmd.util.cbay(triad_sel)
        cmd.color("yellow", f"resi {c} and elem C")
        cmd.color("cyan", f"resi {h} and elem C")
        cmd.color("magenta", f"resi {a} and elem C")

        cmd.label(f"resi {c} and name CA", f'"Cys{c}"')
        cmd.label(f"resi {h} and name CA", f'"His{h}"')
        cmd.label(f"resi {a} and name CA", f'"Asp{a}"')
        cmd.set("label_size", 18)
        cmd.set("label_color", "black")
        cmd.set("label_font_id", 7)
        cmd.set("label_outline_color", "white")

        cmd.distance("d1", f"resi {c} and name SG", f"resi {h} and name NE2")
        cmd.distance("d2", f"resi {h} and name ND1", f"resi {a} and name OD1")
        cmd.hide("labels", "d1")
        cmd.hide("labels", "d2")
        cmd.color("black", "d1")
        cmd.color("black", "d2")

        cmd.orient(triad_sel)
        cmd.zoom(triad_sel, buffer=6)
        cmd.ray(1400, 1100)
        cmd.png(triad_png, dpi=300)
        print(f"Figure written: {triad_png}")
    else:
        triad_png = None
        print("Skipping active-site figure: no complete triad found to render.")

    return overall_png, triad_png


# --------------------------------------------------------------------------
# MAIN
# --------------------------------------------------------------------------

def main():
    if len(sys.argv) < 2:
        print("Usage: python3 phaC_structural_validation.py <path_to_result_zip> [output_dir]")
        sys.exit(1)

    zip_path = sys.argv[1]
    out_dir = sys.argv[2] if len(sys.argv) > 2 else "./phaC_structural_results"
    os.makedirs(out_dir, exist_ok=True)

    log_lines = []

    def log(msg=""):
        print(msg)
        log_lines.append(str(msg))

    log("STEP 1: Extract ColabFold result archive")
    log("-" * 70)
    pdb_path = extract_and_find_top_model(zip_path, out_dir)

    log("\nSTEP 2: Load structure and verify against query sequence")
    log("-" * 70)
    structure = load_structure(pdb_path)
    sanity_check_sequence(structure, QUERY_SEQUENCE)

    log("\nSTEP 3: Per-residue confidence (pLDDT) analysis")
    log("-" * 70)
    plddt = get_residue_plddt(structure)
    summarize_confidence(plddt, log)

    log("\nSTEP 4: Catalytic triad geometric verification")
    log("-" * 70)
    triad_result = run_triad_analysis(structure, plddt, log)
    triad = triad_result["triad"] if triad_result else None

    log("\nSTEP 5: Render annotated structural figures")
    log("-" * 70)
    overall_png, triad_png = render_figures(pdb_path, triad, out_dir)

    summary_path = os.path.join(out_dir, "phaC_structural_validation_summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("\n".join(log_lines))

    log(f"\nAll results written to: {out_dir}")
    log(f"Summary log: {summary_path}")


if __name__ == "__main__":
    main()
