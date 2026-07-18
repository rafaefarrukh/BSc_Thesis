"""
step_7.03.py

PHB-vs-growth production envelope: for a range of fixed growth rates
(0 -> mu_max), compute the minimum and maximum achievable PHB flux
(via FVA at each fixed growth rate). This visualizes the classic
growth/production trade-off space, showing:
    - the maximum theoretical PHB flux at zero growth
    - how the achievable PHB range shrinks as growth increases
    - the growth rate beyond which no PHB production is possible

Outputs:
    step_7.03_production_envelope.csv
    step_7.03_production_envelope.png
"""

import os
import cobra
import matplotlib
matplotlib.use("Agg")  # headless/non-interactive backend for saving to file
import matplotlib.pyplot as plt

from cobra.flux_analysis import production_envelope

from step_7_config import (
    MODEL_PATH,
    OUTPUT_DIR,
    BIOMASS_ID,
    PHB_REACTION_ID,
    describe_reaction,
)

POINTS = 25  # number of growth-rate steps across the envelope


def main():
    # Ensure output directory exists
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    # Load model using configuration parameters
    model = cobra.io.read_sbml_model(MODEL_PATH)

    biomass = model.reactions.get_by_id(BIOMASS_ID)
    phb = model.reactions.get_by_id(PHB_REACTION_ID)

    describe_reaction(biomass, "Biomass reaction")
    describe_reaction(phb, "PHB reaction")

    print("Computing production envelope (this runs an FVA at each growth step)...")
    envelope = production_envelope(model, reactions=[biomass.id], objective=phb.id, points=POINTS)

    csv_path = os.path.join(OUTPUT_DIR, "step_7.03_production_envelope.csv")
    envelope.to_csv(csv_path, index=False)

    # --- Plot -----------------------------------------------------------------
    growth_col = biomass.id
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.fill_between(
        envelope[growth_col],
        envelope["flux_minimum"],
        envelope["flux_maximum"],
        alpha=0.3,
        color="tab:green",
        label="Feasible PHB flux range",
    )
    ax.plot(envelope[growth_col], envelope["flux_maximum"], color="tab:green", lw=2, label="Max PHB flux")
    ax.plot(envelope[growth_col], envelope["flux_minimum"], color="tab:green", lw=1, linestyle="--", label="Min PHB flux")
    ax.set_xlabel("Growth rate (1/h)")
    ax.set_ylabel(f"{phb.id} flux (mmol/gDW/h)")
    ax.set_title("PHB Production vs. Growth — Production Envelope")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, "step_7.03_production_envelope.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    max_phb_at_zero_growth = envelope["flux_maximum"].iloc[0]
    max_growth_row = envelope.loc[envelope[growth_col].idxmax()]

    print(f"Max PHB flux at zero growth: {max_phb_at_zero_growth:.4f}")
    print(
        f"At max growth ({max_growth_row[growth_col]:.4f} /h), "
        f"PHB flux range = [{max_growth_row['flux_minimum']:.4f}, "
        f"{max_growth_row['flux_maximum']:.4f}]"
    )
    print(f"\nData written to: {csv_path}")
    print(f"Plot written to: {png_path}")


if __name__ == "__main__":
    main()
