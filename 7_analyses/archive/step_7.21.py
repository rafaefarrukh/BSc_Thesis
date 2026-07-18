"""
step_7.21_strain_design_optknock.py

Strain design: searches for reaction knockouts that force coupling
between growth and PHB production — i.e. knockouts after which a cell
CANNOT grow without also making PHB (the classic OptKnock/RobustKnock
goal), rather than PHB being just one of many optima the cell could
choose.

If the `cameo` package is installed, this script uses its exact MILP
OptKnock implementation (cameo.strain_design.deterministic — true
bilevel optimization). Otherwise it falls back to a greedy heuristic
search implemented directly in cobrapy:

    Coupling score for a candidate knockout set = the MINIMUM PHB flux
    across the feasible space at maximum growth after the knockout
    (found via FVA: minimize PHB flux subject to growth = max growth).
    A high minimum means production is unavoidable at optimal growth —
    i.e. growth and PHB are coupled. The heuristic does:
        1. Score every single-reaction knockout this way.
        2. Greedily extend the best single knockouts with a second
           knockout (from a restricted candidate pool) to see whether
           coupling can be strengthened further.

This heuristic is NOT the exact bilevel OptKnock MILP — it can miss
combinations the true MILP would find, and doesn't guarantee optimality
— but it runs anywhere cobrapy runs and surfaces the same qualitative
answer (which knockouts push the network toward growth-coupled
production) without extra solver dependencies.

Outputs:
    step_7_outputs/step_7.21_single_knockout_coupling_scores.csv
    step_7_outputs/step_7.21_best_double_knockout_coupling.csv
    step_7_outputs/step_7.21_strain_design_summary.txt
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
MIN_VIABLE_GROWTH_FRACTION = 0.10   # candidate knockouts must retain >= this much wild-type growth
CANDIDATE_POOL_SIZE = 40            # how many top single knockouts to try extending to doubles
TRY_CAMEO_FIRST = True


def coupling_score(model, biomass, phb, knockout_ids):
    """Returns (max_growth, min_PHB_at_max_growth) after knocking out the
    given reaction IDs, or (0, None) if the knockout set is lethal."""
    with model as m:
        for rid in knockout_ids:
            m.reactions.get_by_id(rid).knock_out()
        m.objective = biomass.id
        max_growth = m.slim_optimize()
        if not max_growth or max_growth <= 1e-9:
            return 0.0, None
        biomass_m = m.reactions.get_by_id(biomass.id)
        biomass_m.lower_bound = max_growth
        biomass_m.upper_bound = max_growth
        m.objective = phb.id
        m.objective_direction = "min"
        min_phb = m.slim_optimize()
        m.objective_direction = "max"
    return max_growth, (min_phb if min_phb else 0.0)


def try_cameo_optknock(model, biomass, phb):
    try:
        from cameo.strain_design.deterministic.linear_programming import OptKnock
    except ImportError:
        return None

    print("cameo detected — running exact MILP OptKnock (this may take a while)...")
    optknock = OptKnock(model, fraction_of_optimum=0.1)
    result = optknock.run(
        max_knockouts=2, biomass=biomass.id, target=phb.id, essential_reactions=None
    )
    return result.data_frame


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
    print(f"\nWild-type growth: {wt_growth:.6f} /h")

    if TRY_CAMEO_FIRST:
        cameo_result = try_cameo_optknock(model, biomass, phb)
        if cameo_result is not None:
            csv_path = os.path.join(OUTPUT_DIR, "step_7.21_cameo_optknock_results.csv")
            cameo_result.to_csv(csv_path)
            print(f"\ncameo OptKnock results written to: {csv_path}")
            return
        else:
            print("cameo not available — falling back to the cobrapy-only heuristic search.\n")

    # --- Heuristic single-knockout coupling screen -----------------------------
    candidates = [
        r.id for r in model.reactions
        if not r.boundary and r.id not in (biomass.id, phb.id)
    ]
    print(f"Screening {len(candidates)} single-reaction knockouts for growth/PHB coupling "
          f"(this is the slow step)...")

    rows = []
    for i, rid in enumerate(candidates):
        max_growth, min_phb = coupling_score(model, biomass, phb, [rid])
        if max_growth >= MIN_VIABLE_GROWTH_FRACTION * wt_growth:
            rows.append({
                "knockouts": rid, "n_knockouts": 1,
                "max_growth": max_growth, "min_PHB_at_max_growth": min_phb,
            })
        if (i + 1) % 200 == 0:
            print(f"  ...{i + 1}/{len(candidates)} screened")

    single_df = pd.DataFrame(rows).sort_values("min_PHB_at_max_growth", ascending=False)
    single_csv = os.path.join(OUTPUT_DIR, "step_7.21_single_knockout_coupling_scores.csv")
    single_df.to_csv(single_csv, index=False)
    print(f"\nSingle-knockout coupling scores written to: {single_csv}")
    print("\nTop 10 single knockouts by forced minimum PHB flux at max growth:")
    print(single_df.head(10).to_string(index=False))

    # --- Greedy double-knockout extension of the top single candidates --------
    top_singles = single_df.head(CANDIDATE_POOL_SIZE)["knockouts"].tolist()
    print(f"\nExtending top {len(top_singles)} single knockouts with a second knockout "
          f"(greedy search, this is also slow)...")

    double_rows = []
    for base_rid in top_singles:
        for other_rid in top_singles:
            if other_rid == base_rid:
                continue
            pair = tuple(sorted([base_rid, other_rid]))
            max_growth, min_phb = coupling_score(model, biomass, phb, list(pair))
            if max_growth >= MIN_VIABLE_GROWTH_FRACTION * wt_growth:
                double_rows.append({
                    "knockouts": " + ".join(pair), "n_knockouts": 2,
                    "max_growth": max_growth, "min_PHB_at_max_growth": min_phb,
                })

    double_df = pd.DataFrame(double_rows).drop_duplicates(subset="knockouts")
    double_df = double_df.sort_values("min_PHB_at_max_growth", ascending=False)
    double_csv = os.path.join(OUTPUT_DIR, "step_7.21_best_double_knockout_coupling.csv")
    double_df.to_csv(double_csv, index=False)
    print(f"\nDouble-knockout coupling scores written to: {double_csv}")
    if len(double_df):
        print("\nTop 10 double knockouts by forced minimum PHB flux at max growth:")
        print(double_df.head(10).to_string(index=False))

    # --- Summary -----------------------------------------------------------------
    lines = ["STEP 7.21 — Strain Design (heuristic growth/PHB coupling search)", "=" * 70]
    lines.append(f"Wild-type growth: {wt_growth:.6f} /h; wild-type min PHB at max growth: "
                  f"{coupling_score(model, biomass, phb, [])[1]:.6f}")
    if len(single_df):
        best_single = single_df.iloc[0]
        lines.append(
            f"\nBest single knockout: {best_single['knockouts']} -> "
            f"max growth {best_single['max_growth']:.4f} /h, "
            f"forced min PHB flux {best_single['min_PHB_at_max_growth']:.4f}"
        )
    if len(double_df):
        best_double = double_df.iloc[0]
        lines.append(
            f"Best double knockout: {best_double['knockouts']} -> "
            f"max growth {best_double['max_growth']:.4f} /h, "
            f"forced min PHB flux {best_double['min_PHB_at_max_growth']:.4f}"
        )
    summary_path = os.path.join(OUTPUT_DIR, "step_7.21_strain_design_summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("\n".join(lines))
    print("\n" + "\n".join(lines))
    print(f"\nSummary written to: {summary_path}")


if __name__ == "__main__":
    main()
