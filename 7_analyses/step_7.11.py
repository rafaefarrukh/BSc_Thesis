"""
step_7.11.py

Ranks every reaction in the model by its graph distance (number of
reaction "hops" through shared, non-currency metabolites) from PHB
synthesis. The PHB synthase reaction itself is distance 0; a reaction
that shares a non-currency metabolite directly with it (e.g. an
acetoacetyl-CoA reductase, "AACOAR_syn"-style reaction) is distance 1;
and so on outward through the network.

Also renders a "PHB pathway reaction map": a bipartite reaction/metabolite
network diagram restricted to the near neighborhood of PHB synthesis
(distance <= MAP_MAX_DISTANCE), colored by distance, so the core PHB
biosynthetic route and its immediate branch points are visible at a glance.

Requires: networkx (pip install networkx)

Outputs:
    step_7_outputs/step_7.11_phb_pathway_distances.csv   (ALL reactions, ranked)
    step_7_outputs/step_7.11_phb_pathway_map.png
"""

import os
import csv
import collections
import cobra
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    import networkx as nx
except ImportError as exc:
    raise SystemExit(
        "This script requires networkx. Install it with: pip install networkx"
    ) from exc

from step_7_config import (
    MODEL_PATH,
    OUTPUT_DIR,
    PHB_REACTION_ID,
    describe_reaction,
)

# --------------------------------------------------------------------------
# CONFIG
# --------------------------------------------------------------------------
MAP_MAX_DISTANCE = 2   # how far out from PHB synthesis to draw the map

# Common currency metabolites to ignore to prevent artificial shortcuts in the network
CURRENCY_METABOLITES = {
    "atp_c", "atp_e", "adp_c", "amp_c", "nad_c", "nadh_c", "nadp_c", "nadph_c",
    "h2o_c", "h2o_e", "h_c", "h_e", "co2_c", "co2_e", "pi_c", "ppi_c", "coa_c",
    "nh4_c", "nh4_e", "o2_c", "o2_e", "na1_c", "na1_e", "cl_c", "cl_e"
}


def base_metabolite_id(met):
    """Strip compartment suffixes (e.g., '_c', '_e') for cleaner visualization labels."""
    if met.id.endswith("_c") or met.id.endswith("_e"):
        return met.id[:-2]
    return met.id


def is_currency_metabolite(met):
    """Check if a metabolite is a high-degree connection currency molecule."""
    return met.id in CURRENCY_METABOLITES


def build_reaction_graph(model):
    """Builds an adjacency graph where reactions are connected if they share a non-currency metabolite."""
    graph = collections.defaultdict(set)
    met_to_rxns = collections.defaultdict(list)

    # Map metabolites to their participating reactions
    for rxn in model.reactions:
        if rxn.boundary:
            continue
        for met in rxn.metabolites:
            if is_currency_metabolite(met):
                continue
            met_to_rxns[met.id].append(rxn.id)

    # Link reactions sharing a non-currency metabolite
    for rxns in met_to_rxns.values():
        for r1 in rxns:
            for r2 in rxns:
                if r1 != r2:
                    graph[r1].add(r2)
                    graph[r2].add(r1)
    return graph


def bfs_reaction_distance(graph, start_ids):
    """Computes shortest reaction path distances using Breadth-First Search (BFS)."""
    distances = {}
    queue = collections.deque()

    for start_id in start_ids:
        if start_id in graph or start_id:
            distances[start_id] = 0
            queue.append(start_id)

    while queue:
        current = queue.popleft()
        current_dist = distances[current]

        for neighbor in graph.get(current, []):
            if neighbor not in distances:
                distances[neighbor] = current_dist + 1
                queue.append(neighbor)
    return distances


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    model = cobra.io.read_sbml_model(MODEL_PATH)

    phb_rxn = model.reactions.get_by_id(PHB_REACTION_ID)
    describe_reaction(phb_rxn, f"PHB synthesis start ({PHB_REACTION_ID})")
    start_ids = [phb_rxn.id]

    print("\nBuilding reaction-adjacency graph (excluding currency metabolites)...")
    graph = build_reaction_graph(model)

    print("Running BFS from PHB synthesis...")
    distances = bfs_reaction_distance(graph, start_ids)

    unreached = len(model.reactions) - len(distances)
    print(f"  {len(distances)} / {len(model.reactions)} reactions reached "
          f"({unreached} unreachable/boundary reactions excluded)")

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

    csv_path = os.path.join(OUTPUT_DIR, "step_7.11_phb_pathway_distances.csv")
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

    png_path = os.path.join(OUTPUT_DIR, "step_7.11_phb_pathway_map.png")
    fig.savefig(png_path, dpi=200)
    plt.close(fig)

    print(f"Pathway map written to: {png_path}")


if __name__ == "__main__":
    main()
