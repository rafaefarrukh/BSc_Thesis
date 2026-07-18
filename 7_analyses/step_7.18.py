"""
step_7.18.py

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
    neighborhood from step_7.11, if available, else all non-essential
    reactions) to keep runtime bounded. This is NOT exhaustive over all
    pairs in a genome-scale model and is not guaranteed to find the true
    minimal cut sets — treat it as a practical candidate-generation pass.

Requires step_7.11 to have been run first for the (recommended,
smaller/faster) pathway-scoped candidate pool; falls back to a size-
limited pool of non-essential reactions otherwise.

--- Notes on the exact (straindesign) path ---

Exact MCS enumeration is a MILP. A few things matter a lot for runtime:
  - Candidate knockouts are now restricted via `ko_cost` to the same
    pathway-scoped candidate pool used by the heuristic fallback,
    instead of leaving every reaction in the (genome-scale) model as a
    free knockout candidate. Without this the search space is far too
    large for exact enumeration to finish in a reasonable time,
    especially with GLPK's big-M MILP formulation.
  - A STRAINDESIGN_TIME_LIMIT_S wall-clock limit is passed to
    compute_strain_designs so a hard MILP instance fails fast and falls
    back to the heuristic search, instead of hanging indefinitely.
  - The growth-viability constraint passed to straindesign is now
    correctly expressed as VIABILITY_THRESHOLD_FRACTION * wt_growth
    (an absolute growth rate), matching the fraction-of-wild-type
    semantics used everywhere else in this script. Previously the raw
    fraction (e.g. 0.10) was passed directly as if it were an absolute
    growth rate, which happened to be harmless here only because
    0.10 < wt_growth, but would silently be wrong in general.
  - compute_strain_designs() returns an SDSolutions object, not a plain
    list. Solutions live in its `.sd` attribute: a list of dicts mapping
    reaction IDs to -1 (knocked out), 0 (untouched), or 1 (added). We
    parse that correctly now instead of iterating over the object
    directly or assuming a `.knockouts` attribute exists.

Outputs:
    step_7_outputs/step_7.18_straindesign_mcs_results.csv  (if straindesign succeeds)
    step_7_outputs/step_7.18_minimal_cut_sets.csv          (heuristic fallback)
    step_7_outputs/step_7.18_mcs_summary.txt
"""

import os
import pandas as pd
import cobra

from step_7_config import (
    MODEL_PATH,
    OUTPUT_DIR,
    BIOMASS_ID,
    PHB_REACTION_ID,
    describe_reaction,
)

PHB_LEAK_THRESHOLD = 1e-3         # PHB flux below this counts as "no production"
VIABILITY_THRESHOLD_FRACTION = 0.10  # growth must fall below this fraction of wild-type to count as "cut"
MAX_CANDIDATE_POOL = 60           # cap on candidate reactions considered for pairs (runtime control)
TRY_STRAINDESIGN_FIRST = True
STRAINDESIGN_MAX_SOLUTIONS = 5    # exact enumeration: how many MCS to look for (fewer = faster)
STRAINDESIGN_TIME_LIMIT_S = 600   # wall-clock budget before giving up and falling back to heuristic

DISTANCE_CSV = os.path.join(OUTPUT_DIR, "step_7.11_phb_pathway_distances.csv")


def get_candidate_pool(model, biomass, phb):
    if os.path.exists(DISTANCE_CSV):
        dist_df = pd.read_csv(DISTANCE_CSV)
        ids = dist_df.sort_values("distance_from_PHB")["reaction_id"].tolist()
    else:
        print(f"NOTE: {DISTANCE_CSV} not found (run step_7.11 first for a smarter, "
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


def try_straindesign_mcs(model, biomass, phb, wt_growth, candidate_pool):
    try:
        import straindesign as sd
    except ImportError:
        return None

    print("straindesign detected — running exact MCS enumeration...")
    print(f"  candidate knockout pool restricted to {len(candidate_pool)} reactions "
          f"(pathway-neighborhood-prioritized if step_7.11 was run)")

    # Restrict knockout candidates to the same pathway-scoped pool used by the
    # heuristic fallback. Without this, every reaction in the genome-scale
    # model is a free knockout candidate, which makes exact enumeration
    # infeasible in practice (this is what was causing multi-hour hangs).
    # NOTE: ko_cost is an argument of compute_strain_designs(), not of
    # SDModule() -- passing it to SDModule raises "Key ko_cost is not
    # supported."
    ko_cost = {rid: 1 for rid in candidate_pool}

    viability_growth = VIABILITY_THRESHOLD_FRACTION * wt_growth

    module = sd.SDModule(
        model,
        sd.names.SUPPRESS,
        constraints=[
            f"{phb.id} <= {PHB_LEAK_THRESHOLD}",
            f"{biomass.id} >= {viability_growth}",
        ],
    )
    solutions = sd.compute_strain_designs(
        model,
        sd_modules=[module],
        ko_cost=ko_cost,
        max_solutions=STRAINDESIGN_MAX_SOLUTIONS,
        time_limit=STRAINDESIGN_TIME_LIMIT_S,
    )
    return solutions


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    model = cobra.io.read_sbml_model(MODEL_PATH)

    biomass = model.reactions.get_by_id(BIOMASS_ID)
    phb = model.reactions.get_by_id(PHB_REACTION_ID)

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")

    model.objective = biomass.id
    wt_growth = model.slim_optimize()
    viability_cutoff = VIABILITY_THRESHOLD_FRACTION * wt_growth
    baseline_growth_no_phb = growth_without_phb(model, biomass, phb, [])
    print(f"\nWild-type growth: {wt_growth:.6f} /h")
    print(f"Growth achievable WITHOUT PHB production (baseline, no knockouts): "
          f"{baseline_growth_no_phb:.6f} /h")

    # Candidate pool is computed up front so both the exact and heuristic
    # paths use the same pathway-scoped restriction.
    candidates = get_candidate_pool(model, biomass, phb)

    if TRY_STRAINDESIGN_FIRST:
        try:
            sd_result = try_straindesign_mcs(model, biomass, phb, wt_growth, candidates)
        except Exception as exc:  # noqa: BLE001
            print(f"straindesign MCS computation failed ({exc}); falling back to heuristic search.")
            sd_result = None

        if sd_result is not None:
            # compute_strain_designs() returns an SDSolutions object, not a
            # plain list. Solutions live in `.sd`: a list of dicts mapping
            # reaction IDs to -1 (knocked out), 0 (untouched), 1 (added).
            solution_dicts = getattr(sd_result, "sd", None)
            if not solution_dicts:
                print("straindesign returned no solutions (search space may be "
                      "too restricted, or the time limit was hit); falling back "
                      "to heuristic search.\n")
            else:
                formatted_results = []
                for design in solution_dicts:
                    knockouts = [rid for rid, val in design.items() if val == -1]
                    formatted_results.append({
                        "knockouts": ", ".join(knockouts),
                        "n_knockouts": len(knockouts),
                    })
                csv_path = os.path.join(OUTPUT_DIR, "step_7.18_straindesign_mcs_results.csv")
                pd.DataFrame(formatted_results).to_csv(csv_path, index=False)
                print(f"\n{len(formatted_results)} exact minimal cut set(s) found via straindesign.")
                print(f"straindesign MCS results written to: {csv_path}")
                return
        else:
            print("straindesign not available — falling back to the cobrapy-only heuristic search.\n")

    print(f"Candidate pool for MCS search: {len(candidates)} reactions "
          f"(pathway-neighborhood-prioritized if step_7.11 was run)")

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

    if found_sets:
        df = pd.DataFrame(found_sets).sort_values(
            ["size", "max_growth_after_cut"], ascending=[True, False]
        )
    else:
        # pd.DataFrame([]) has no columns at all, so sort_values on "size"
        # would raise KeyError -- this branch is what previously crashed
        # when neither single nor pair cut sets were found.
        df = pd.DataFrame(columns=["cut_set", "size", "growth_without_PHB_after_cut", "max_growth_after_cut"])
    csv_path = os.path.join(OUTPUT_DIR, "step_7.18_minimal_cut_sets.csv")
    df.to_csv(csv_path, index=False)

    lines = ["STEP 7.18 — Minimal Cut Sets (heuristic search)", "=" * 70]
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
    summary_path = os.path.join(OUTPUT_DIR, "step_7.18_mcs_summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("\n".join(lines))

    print("\n" + "\n".join(lines))
    print(f"\nCut sets table written to: {csv_path}")
    print(f"Summary written to: {summary_path}")


if __name__ == "__main__":
    main()
