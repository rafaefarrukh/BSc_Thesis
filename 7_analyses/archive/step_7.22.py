"""
step_7.22_minimal_cut_sets.py

Minimal Cut Sets (MCS): a minimal cut set is a minimal-size set of
reaction knockouts that eliminates an undesired phenotype — here, the
undesired phenotype is "growth without PHB production" (i.e. any flux
state where the cell grows above a viability threshold but makes little
or no PHB). Removing a minimal cut set forces every remaining growth
phenotype to also make PHB.

If the `straindesign` package is installed, this script uses its exact
constraint-based MCS enumeration (the standard approach, based on
Ballerstein/von Kamp/Klamt's dual-network formulation). Otherwise it
falls back to a direct enumeration heuristic implemented in cobrapy:

    A candidate knockout set S is a cut set if, after removing S, the LP
        maximize growth
        subject to: PHB flux <= PHB_LEAK_THRESHOLD
    becomes infeasible or has growth below VIABILITY_THRESHOLD — i.e. you
    can no longer grow without also producing PHB.

    Singles are tested exhaustively; if none work alone, pairs are tested
    from a restricted candidate pool (reactions within the PHB pathway
    neighborhood from step_7.12, if available, else all non-essential
    reactions) to keep runtime bounded. This is NOT exhaustive over all
    pairs in a genome-scale model and is not guaranteed to find the true
    minimal cut sets — treat it as a practical candidate-generation pass.

Requires step_7.12 to have been run first for the (recommended,
smaller/faster) pathway-scoped candidate pool; falls back to a size-
limited pool of non-essential reactions otherwise.

Outputs:
    step_7_outputs/step_7.22_minimal_cut_sets.csv
    step_7_outputs/step_7.22_mcs_summary.txt
"""

import os

import pandas as pd

from step_7_utils import (
    MODEL_PATH,
    OUTPUT_DIR,
    describe_reaction,
    ensure_output_dir,
    find_biomass_reaction,
    find_phb_reaction,
    load_model,
)

PHB_REACTION_ID = "PHBS_syn_1"
PHB_LEAK_THRESHOLD = 1e-3         # PHB flux below this counts as "no production"
VIABILITY_THRESHOLD_FRACTION = 0.10  # growth must fall below this fraction of wild-type to count as "cut"
MAX_CANDIDATE_POOL = 60           # cap on candidate reactions considered for pairs (runtime control)
TRY_STRAINDESIGN_FIRST = True

DISTANCE_CSV = os.path.join(OUTPUT_DIR, "step_7.12_phb_pathway_distances.csv")


def get_candidate_pool(model, biomass, phb):
    if os.path.exists(DISTANCE_CSV):
        dist_df = pd.read_csv(DISTANCE_CSV)
        ids = dist_df.sort_values("distance_from_PHB")["reaction_id"].tolist()
    else:
        print(f"NOTE: {DISTANCE_CSV} not found (run step_7.12 first for a smarter, "
              f"smaller candidate pool). Using all non-boundary reactions instead.")
        ids = [r.id for r in model.reactions if not r.boundary]

    ids = [rid for rid in ids if rid not in (biomass.id, phb.id)]
    return ids[:MAX_CANDIDATE_POOL]


def growth_without_phb(model, biomass, phb, knockout_ids):
    """Max growth achievable while PHB flux stays <= PHB_LEAK_THRESHOLD,
    after knocking out the given reactions."""
    with model as m:
        for rid in knockout_ids:
            m.reactions.get_by_id(rid).knock_out()
        m.reactions.get_by_id(phb.id).upper_bound = min(
            m.reactions.get_by_id(phb.id).upper_bound, PHB_LEAK_THRESHOLD
        )
        m.objective = biomass.id
        g = m.slim_optimize()
    return g if g else 0.0


def try_straindesign_mcs(model, biomass, phb):
    try:
        import straindesign as sd
    except ImportError:
        return None

    print("straindesign detected — running exact MCS enumeration...")
    module = sd.SDModule(
        model,
        sd.names.SUPPRESS,
        constraints=[f"{phb.id} <= {PHB_LEAK_THRESHOLD}", f"{biomass.id} >= "
                     f"{VIABILITY_THRESHOLD_FRACTION}"],
    )
    solutions = sd.compute_strain_designs(model, sd_modules=[module], max_solutions=20)
    return solutions


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)
    biomass = find_biomass_reaction(model)
    phb = model.reactions.get_by_id(PHB_REACTION_ID) if PHB_REACTION_ID else find_phb_reaction(model)
    if phb is None:
        raise ValueError("Could not auto-detect a PHB reaction. Set PHB_REACTION_ID at the top of this script.")

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")

    model.objective = biomass.id
    wt_growth = model.slim_optimize()
    viability_cutoff = VIABILITY_THRESHOLD_FRACTION * wt_growth
    baseline_growth_no_phb = growth_without_phb(model, biomass, phb, [])
    print(f"\nWild-type growth: {wt_growth:.6f} /h")
    print(f"Growth achievable WITHOUT PHB production (baseline, no knockouts): "
          f"{baseline_growth_no_phb:.6f} /h")

    if TRY_STRAINDESIGN_FIRST:
        try:
            sd_result = try_straindesign_mcs(model, biomass, phb)
        except Exception as exc:  # noqa: BLE001
            print(f"straindesign MCS computation failed ({exc}); falling back to heuristic search.")
            sd_result = None
        if sd_result is not None:
            csv_path = os.path.join(OUTPUT_DIR, "step_7.22_straindesign_mcs_results.csv")
            pd.DataFrame(sd_result).to_csv(csv_path)
            print(f"\nstraindesign MCS results written to: {csv_path}")
            return
        print("straindesign not available — falling back to the cobrapy-only heuristic search.\n")

    candidates = get_candidate_pool(model, biomass, phb)
    print(f"Candidate pool for MCS search: {len(candidates)} reactions "
          f"(pathway-neighborhood-prioritized if step_7.12 was run)")

    # --- Test singles ------------------------------------------------------------
    print("\nTesting single-reaction cut sets...")
    found_sets = []
    for rid in candidates:
        g_no_phb = growth_without_phb(model, biomass, phb, [rid])
        if g_no_phb <= viability_cutoff:
            # confirm the knockout doesn't just kill the cell outright
            with model as m:
                m.reactions.get_by_id(rid).knock_out()
                m.objective = biomass.id
                g_normal = m.slim_optimize()
            if g_normal and g_normal > viability_cutoff:
                found_sets.append({
                    "cut_set": rid, "size": 1,
                    "growth_without_PHB_after_cut": g_no_phb,
                    "max_growth_after_cut": g_normal,
                })

    # --- If no singles work, try pairs from the same candidate pool -----------
    if not found_sets:
        print("No single-reaction cut sets found; testing pairs (this is slower)...")
        for i, rid1 in enumerate(candidates):
            for rid2 in candidates[i + 1:]:
                g_no_phb = growth_without_phb(model, biomass, phb, [rid1, rid2])
                if g_no_phb <= viability_cutoff:
                    with model as m:
                        m.reactions.get_by_id(rid1).knock_out()
                        m.reactions.get_by_id(rid2).knock_out()
                        m.objective = biomass.id
                        g_normal = m.slim_optimize()
                    if g_normal and g_normal > viability_cutoff:
                        found_sets.append({
                            "cut_set": f"{rid1} + {rid2}", "size": 2,
                            "growth_without_PHB_after_cut": g_no_phb,
                            "max_growth_after_cut": g_normal,
                        })

    df = pd.DataFrame(found_sets).sort_values(["size", "max_growth_after_cut"], ascending=[True, False])
    csv_path = os.path.join(OUTPUT_DIR, "step_7.22_minimal_cut_sets.csv")
    df.to_csv(csv_path, index=False)

    lines = ["STEP 7.22 — Minimal Cut Sets (heuristic search)", "=" * 70]
    lines.append(f"Candidate pool size: {len(candidates)}")
    lines.append(f"Cut sets found: {len(df)}")
    if len(df):
        lines.append("\nTop candidate cut sets (smallest size, highest retained growth):")
        lines.append(df.head(15).to_string(index=False))
    else:
        lines.append(
            "\nNo cut sets found within the searched pool/size limit. Consider "
            "increasing MAX_CANDIDATE_POOL, testing triples, or installing "
            "'straindesign' for exact enumeration."
        )
    summary_path = os.path.join(OUTPUT_DIR, "step_7.22_mcs_summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("\n".join(lines))

    print("\n" + "\n".join(lines))
    print(f"\nCut sets table written to: {csv_path}")
    print(f"Summary written to: {summary_path}")


if __name__ == "__main__":
    main()
