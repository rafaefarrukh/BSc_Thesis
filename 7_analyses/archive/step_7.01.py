"""
step_7.1_validation_and_stats.py

Basic validation and summary statistics for model_7.xml.

This is a lightweight, independent sanity check to run alongside MEMOTE
(see MEMOTE_Report_Summary___model_7.md for the full scoring report).
It reports:
    - SBML/COBRA validation (structural + FBC consistency errors/warnings)
    - Model size (reactions, metabolites, genes, compartments)
    - Mass/charge balance summary
    - Reactions without GPR
    - Blocked reaction count
    - Basic solvability (can the model produce nonzero flux at all?)

Outputs:
    step_7_outputs/step_7.1_model_summary.txt
    step_7_outputs/step_7.1_mass_charge_imbalances.csv
"""

import csv
import os

import cobra
from cobra.flux_analysis import find_blocked_reactions

from step_7_utils import MODEL_PATH, OUTPUT_DIR, ensure_output_dir, find_biomass_reaction

# Set to False to skip the blocked-reaction scan on very large models
# (it requires one FVA-style optimization per reaction and can be slow)
RUN_BLOCKED_REACTION_SCAN = True


def validate_and_load(model_path):
    """Use cobra's SBML validator (catches structural/FBC issues that a
    plain read_sbml_model call would silently ignore or warn about)."""
    print(f"Validating SBML file: {model_path}")
    model, errors = cobra.io.sbml.validate_sbml_model(model_path)
    if model is None:
        raise RuntimeError(f"Model failed to load. Errors: {errors}")
    return model, errors


def summarize_compartments(model):
    return {cid: cname for cid, cname in model.compartments.items()}


def mass_charge_balance_report(model, out_csv):
    """Check every reaction's mass/charge balance (skips exchange/
    boundary reactions, which are expected to be unbalanced by design)."""
    imbalanced = []
    for rxn in model.reactions:
        if rxn.boundary:
            continue
        balance = rxn.check_mass_balance()
        if balance:
            imbalanced.append((rxn.id, rxn.name, str(balance)))

    with open(out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["reaction_id", "reaction_name", "imbalance"])
        writer.writerows(imbalanced)

    return imbalanced


def main():
    ensure_output_dir(OUTPUT_DIR)
    model, sbml_errors = validate_and_load(MODEL_PATH)

    lines = []
    lines.append("=" * 70)
    lines.append("STEP 7.1 — MODEL VALIDATION AND SUMMARY STATISTICS")
    lines.append("=" * 70)
    lines.append(f"Model file: {MODEL_PATH}")
    lines.append(f"Model ID:   {model.id}")
    lines.append("")

    # --- SBML validation -----------------------------------------------
    lines.append("--- SBML / COBRA validation ---")
    if sbml_errors and any(sbml_errors.values()):
        for category, issues in sbml_errors.items():
            if issues:
                lines.append(f"  [{category}] {len(issues)} issue(s):")
                for issue in issues[:20]:
                    lines.append(f"      - {issue}")
                if len(issues) > 20:
                    lines.append(f"      ... and {len(issues) - 20} more")
    else:
        lines.append("  No structural validation errors reported.")
    lines.append("")

    # --- Basic size ------------------------------------------------------
    lines.append("--- Model size ---")
    lines.append(f"  Reactions:    {len(model.reactions)}")
    lines.append(f"  Metabolites:  {len(model.metabolites)}")
    lines.append(f"  Genes:        {len(model.genes)}")
    compartments = summarize_compartments(model)
    lines.append(f"  Compartments: {len(compartments)} -> {compartments}")
    exchanges = model.exchanges
    lines.append(f"  Exchange reactions: {len(exchanges)}")
    lines.append("")

    # --- GPR coverage ------------------------------------------------------
    no_gpr = [r for r in model.reactions if not r.boundary and not r.gene_reaction_rule]
    lines.append("--- Gene-Protein-Reaction (GPR) coverage ---")
    lines.append(
        f"  Reactions without a GPR rule (excluding boundary rxns): "
        f"{len(no_gpr)} / {len(model.reactions) - len(exchanges)}"
    )
    lines.append("")

    # --- Mass/charge balance ------------------------------------------------
    lines.append("--- Mass/charge balance (non-boundary reactions) ---")
    imbalance_csv = os.path.join(OUTPUT_DIR, "step_7.1_mass_charge_imbalances.csv")
    imbalanced = mass_charge_balance_report(model, imbalance_csv)
    non_boundary = len(model.reactions) - len(exchanges)
    lines.append(
        f"  Mass/charge-imbalanced reactions: {len(imbalanced)} / {non_boundary} "
        f"({100 * len(imbalanced) / max(non_boundary, 1):.1f}%)"
    )
    lines.append(f"  Full list written to: {imbalance_csv}")
    lines.append("")

    # --- Solvability -------------------------------------------------------
    lines.append("--- Solvability ---")
    try:
        biomass = find_biomass_reaction(model)
        model.objective = biomass.id
        solution = model.optimize()
        lines.append(f"  Biomass reaction: {biomass.id}")
        lines.append(f"  Solver status:    {solution.status}")
        lines.append(f"  Max growth rate:  {solution.objective_value:.6f} /h")
    except Exception as exc:  # noqa: BLE001 - report any solver issue plainly
        lines.append(f"  FAILED to optimize model: {exc}")
    lines.append("")

    # --- Blocked reactions ---------------------------------------------------
    if RUN_BLOCKED_REACTION_SCAN:
        lines.append("--- Blocked reactions ---")
        print("Scanning for blocked reactions (this can take a while on large models)...")
        blocked = find_blocked_reactions(model)
        lines.append(
            f"  Blocked reactions: {len(blocked)} / {len(model.reactions)} "
            f"({100 * len(blocked) / len(model.reactions):.1f}%)"
        )
        lines.append("")

    summary_path = os.path.join(OUTPUT_DIR, "step_7.1_model_summary.txt")
    with open(summary_path, "w") as fh:
        fh.write("\n".join(lines))

    print("\n".join(lines))
    print(f"\nSummary written to: {summary_path}")


if __name__ == "__main__":
    main()
