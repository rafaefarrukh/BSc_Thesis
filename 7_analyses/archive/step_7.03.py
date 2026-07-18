"""
step_7.3_pfba_at_phb_optimum.py

Parsimonious FBA (pFBA) at the PHB production optimum.

Plain FBA solutions are typically highly degenerate — many different flux
distributions achieve the same optimal objective value, including ones
that route flux through unnecessary loops. pFBA resolves this by first
finding the optimal objective value (max PHB flux here), then finding the
flux distribution that achieves that same value with the minimum possible
total enzyme usage (sum of absolute fluxes). This is generally taken as a
more physiologically realistic solution, since cells are expected to
minimize unnecessary protein cost.

Outputs:
    step_7_outputs/step_7.3_pfba_flux_distribution.csv   (all reaction fluxes)
    step_7_outputs/step_7.3_pfba_summary.txt
"""

import csv
import os

from cobra.flux_analysis import pfba

from step_7_utils import (
    MODEL_PATH,
    OUTPUT_DIR,
    describe_reaction,
    ensure_output_dir,
    find_biomass_reaction,
    find_phb_reaction,
    load_model,
)

# Set manually if auto-detection picks the wrong reaction
PHB_REACTION_ID = "PHBS_syn_1"

# How many top (by |flux|) reactions to print/report in the summary
TOP_N_REACTIONS = 25


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)

    biomass = find_biomass_reaction(model)
    phb = model.reactions.get_by_id(PHB_REACTION_ID) if PHB_REACTION_ID else find_phb_reaction(model)
    if phb is None:
        raise ValueError(
            "Could not auto-detect a PHB reaction. Set PHB_REACTION_ID at the "
            "top of this script."
        )

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")

    # Objective = PHB production; pfba() will optimize this first, then
    # minimize total flux subject to achieving that same optimal value.
    model.objective = phb.id
    pfba_solution = pfba(model)

    phb_flux = pfba_solution.fluxes[phb.id]
    growth_flux = pfba_solution.fluxes[biomass.id]
    total_flux = pfba_solution.fluxes.abs().sum()
    active_reactions = (pfba_solution.fluxes.abs() > 1e-9).sum()

    # --- Write full flux distribution ---------------------------------------
    csv_path = os.path.join(OUTPUT_DIR, "step_7.3_pfba_flux_distribution.csv")
    fluxes_sorted = pfba_solution.fluxes.reindex(
        pfba_solution.fluxes.abs().sort_values(ascending=False).index
    )
    with open(csv_path, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["reaction_id", "reaction_name", "flux"])
        for rxn_id, flux in fluxes_sorted.items():
            rxn = model.reactions.get_by_id(rxn_id)
            writer.writerow([rxn_id, rxn.name, flux])

    # --- Summary ---------------------------------------------------------------
    lines = ["=" * 70, "STEP 7.3 — pFBA AT PHB OPTIMUM", "=" * 70]
    lines.append(f"Objective (maximized): {phb.id}")
    lines.append(f"Status:                {pfba_solution.status}")
    lines.append(f"PHB flux:              {phb_flux:.6f}")
    lines.append(f"Growth rate at optimum: {growth_flux:.6f} /h")
    lines.append(f"Sum of |fluxes| (total parsimonious flux): {total_flux:.4f}")
    lines.append(f"Reactions carrying nonzero flux: {active_reactions} / {len(model.reactions)}")
    lines.append("")
    lines.append(f"Top {TOP_N_REACTIONS} reactions by |flux|:")
    for rxn_id, flux in fluxes_sorted.head(TOP_N_REACTIONS).items():
        rxn = model.reactions.get_by_id(rxn_id)
        lines.append(f"  {rxn_id:20s} {flux:12.4f}  {rxn.name}")

    txt_path = os.path.join(OUTPUT_DIR, "step_7.3_pfba_summary.txt")
    with open(txt_path, "w") as fh:
        fh.write("\n".join(lines))

    print("\n".join(lines))
    print(f"\nFull flux distribution written to: {csv_path}")
    print(f"Summary written to: {txt_path}")


if __name__ == "__main__":
    main()
