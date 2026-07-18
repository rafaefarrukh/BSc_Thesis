"""
step_7.12_phb_pathway_distance_and_map.py

Ranks every reaction in the model by its graph distance (number of
reaction "hops" through shared, non-currency metabolites) from PHB
synthesis. The PHB synthase reaction itself is distance 0; a reaction
that shares a non-currency metabolite directly with it (e.g. an
acetoacetyl-CoA reductase, "AACOAR_syn"-style reaction) is distance 1;
and so on outward through the network. This is the same idea as the
example in the request: PHB_syn_1 = 0, AACOAR_syn = 1, etc.

Also renders a "PHB pathway reaction map": a bipartite reaction/metabolite
network diagram restricted to the near neighborhood of PHB synthesis
(distance <= MAP_MAX_DISTANCE), colored by distance, so the core PHB
biosynthetic route and its immediate branch points are visible at a glance.

Requires: networkx (pip install networkx)

Outputs:
    step_7_outputs/step_7.12_phb_pathway_distances.csv   (ALL reactions, ranked)
    step_7_outputs/step_7.12_phb_pathway_map.png
"""

import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import networkx as nx
except ImportError as exc:  # pragma: no cover
    raise SystemExit(
        "This script requires networkx. Install it with: pip install networkx"
    ) from exc

from step_7_utils import (
    MODEL_PATH,
    OUTPUT_DIR,
    base_metabolite_id,
    bfs_reaction_distance,
    build_reaction_graph,
    describe_reaction,
    ensure_output_dir,
    find_phb_reaction,
    is_currency_metabolite,
    load_model,
)

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------
# Set explicitly if the auto-detected PHB reaction isn't the actual PHA
# synthase step (e.g. it picked an exchange/demand reaction instead) —
# for pathway-distance purposes you want the internal synthesis reaction,
# e.g. "PHB_syn_1". Can be a single ID or a list of IDs (multiple starting
# points, e.g. if PHB is produced by more than one isozyme/route).
PHB_SYN_REACTION_ID = "PHBS_syn_1"

MAP_MAX_DISTANCE = 2   # how far out from PHB synthesis to draw the map


def main():
    ensure_output_dir(OUTPUT_DIR)
    model = load_model(MODEL_PATH)

    if PHB_SYN_REACTION_ID:
        start_ids = (
            [PHB_SYN_REACTION_ID] if isinstance(PHB_SYN_REACTION_ID, str) else list(PHB_SYN_REACTION_ID)
        )
        for rid in start_ids:
            describe_reaction(model.reactions.get_by_id(rid), f"PHB synthesis start ({rid})")
    else:
        phb_rxn = find_phb_reaction(model, prefer_exchange=False)
        if phb_rxn is None:
            raise ValueError(
                "Could not auto-detect the PHB synthase reaction. Set "
                "PHB_SYN_REACTION_ID at the top of this script (e.g. 'PHB_syn_1')."
            )
        describe_reaction(phb_rxn, "PHB synthesis start (auto-detected)")
        start_ids = [phb_rxn.id]

    print("\nBuilding reaction-adjacency graph (excluding currency metabolites)...")
    graph = build_reaction_graph(model, exclude_currency=True)

    print("Running BFS from PHB synthesis...")
    distances = bfs_reaction_distance(graph, start_ids)

    unreached = len(graph) - len(distances)
    print(f"  {len(distances)} / {len(graph)} non-boundary reactions reached "
          f"({unreached} unreachable at any distance, excluding currency metabolites)")

    # --- Full ranking CSV -----------------------------------------------------
    rows = []
    for rxn_id, dist in distances.items():
        rxn = model.reactions.get_by_id(rxn_id)
        rows.append({
            "reaction_id": rxn_id,
            "reaction_name": rxn.name,
            "subsystem": getattr(rxn, "subsystem", ""),
            "distance_from_PHB": dist,
        })
    rows.sort(key=lambda r: (r["distance_from_PHB"], r["reaction_id"]))

    import csv
    csv_path = os.path.join(OUTPUT_DIR, "step_7.12_phb_pathway_distances.csv")
    with open(csv_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["reaction_id", "reaction_name", "subsystem", "distance_from_PHB"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nFull ranking written to: {csv_path}")

    print("\nReactions at distance 0-3:")
    for r in rows:
        if r["distance_from_PHB"] <= 3:
            print(f"  [{r['distance_from_PHB']}] {r['reaction_id']:20s} {r['reaction_name']}")

    # --- PHB pathway reaction map (bipartite reaction/metabolite graph) -------
    near_rxn_ids = [rid for rid, d in distances.items() if d <= MAP_MAX_DISTANCE]
    print(f"\nBuilding pathway map for {len(near_rxn_ids)} reactions within distance {MAP_MAX_DISTANCE}...")

    G = nx.Graph()
    for rid in near_rxn_ids:
        rxn = model.reactions.get_by_id(rid)
        G.add_node(rid, kind="reaction", distance=distances[rid], label=rid)
        for met in rxn.metabolites:
            if is_currency_metabolite(met):
                continue
            met_key = f"met::{met.id}"
            if met_key not in G:
                G.add_node(met_key, kind="metabolite", distance=None, label=base_metabolite_id(met))
            G.add_edge(rid, met_key)

    fig, ax = plt.subplots(figsize=(12, 10))
    pos = nx.spring_layout(G, seed=42, k=0.6)

    rxn_nodes = [n for n, d in G.nodes(data=True) if d["kind"] == "reaction"]
    met_nodes = [n for n, d in G.nodes(data=True) if d["kind"] == "metabolite"]

    rxn_colors = [G.nodes[n]["distance"] for n in rxn_nodes]
    nx.draw_networkx_edges(G, pos, ax=ax, alpha=0.3, width=1)
    nx.draw_networkx_nodes(G, pos, nodelist=met_nodes, node_shape="o",
                            node_color="lightgray", node_size=250, ax=ax, label="Metabolite")
    rxn_scatter = nx.draw_networkx_nodes(
        G, pos, nodelist=rxn_nodes, node_shape="s", node_color=rxn_colors,
        cmap="viridis_r", node_size=550, ax=ax,
    )
    nx.draw_networkx_labels(
        G, pos, labels={n: G.nodes[n]["label"] for n in met_nodes}, font_size=6, ax=ax
    )
    nx.draw_networkx_labels(
        G, pos, labels={n: G.nodes[n]["label"] for n in rxn_nodes}, font_size=7,
        font_weight="bold", ax=ax,
    )
    fig.colorbar(rxn_scatter, ax=ax, label="Distance from PHB synthesis (reaction hops)")
    ax.set_title(
        f"PHB Pathway Reaction Map (reactions within {MAP_MAX_DISTANCE} hops; "
        f"squares = reactions, circles = metabolites)"
    )
    ax.axis("off")
    fig.tight_layout()

    png_path = os.path.join(OUTPUT_DIR, "step_7.12_phb_pathway_map.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"Pathway map written to: {png_path}")


if __name__ == "__main__":
    main()
