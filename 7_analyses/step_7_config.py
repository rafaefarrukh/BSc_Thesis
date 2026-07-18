"""
Step 7 Configuration
Contains model path, reaction IDs, and helper functions
"""

# libraries
import os
import sys
import cobra

# paths
MODEL_PATH = "../models/model_5.xml"
OUTPUT_DIR = "step_7_results/"

# reaction ids
BIOMASS_ID = "Growth"
PHB_REACTION_ID = "PHBS_syn_1"
GLUCOSE_EXCHANGE_ID = "EX_glc__D_e"
MALTOSE_EXCHANGE_ID = "EX_malt_e"
OXYGEN_EXCHANGE_ID = "EX_o2_e"
NITROGEN_EXCHANGE_ID = "EX_nh4_e"

# carbon sources
ALTERNATE_CARBON_SOURCES = {
    "Glucose": "EX_glc__D_e",
    "Fructose": "EX_fru_e",
    "Glycerol": "EX_glyc_e",
    "Sucrose": "EX_sucr_e",
    "Maltose": "EX_malt_e",
    "Mannose": "EX_man_e",
    "Mannitol": "EX_mnl_e",
    "Arabinose": "EX_arab__L_e",
    "Melibiose": "EX_melib_e",
    "Acetate": "EX_ac_e",
    "Succinate": "EX_succ_e",
    "Lactate (D)": "EX_lac__D_e",
    "Pyruvate": "EX_pyr_e",
}

# load_model
# model = cobra.io.read_sbml_model(model_path)

# Michaelis Menten Equation
def michaelis_menten(conc, vmax, km):
    """Standard Michaelis-Menten uptake kinetics.
    conc: current extracellular substrate concentration (mM, must be >= 0)
    vmax: maximum uptake rate (mmol/gDW/h)
    km:   half-saturation constant (mM)
    Returns the kinetic uptake rate (mmol/gDW/h), always >= 0.
    """
    if conc <= 0:
        return 0.0
    return vmax * conc / (km + conc)


def growth_coupled_phb(model, biomass, phb, fraction=0.90):
    """
    Standard two-step growth-coupled evaluation used throughout the Step 7
    series:
    (1) find max growth,
    (2) fix growth at `fraction` of that max,  then maximize PHB. Returns (max_growth, phb_flux_at_fixed_growth).
    Must be called with `model` already inside (or about to be wrapped in)
    a `with model as m:` context by the caller if the bound changes should
    not persist — this function itself uses its own `with` block so the
    model is left unmodified afterwards.
    """
    with model as m:
        m.objective = biomass.id
        max_growth = m.slim_optimize()
        max_growth = max_growth if max_growth and max_growth > 1e-9 else 0.0

        phb_flux = 0.0
        if max_growth > 0:
            fixed_growth = fraction * max_growth
            biomass_m = m.reactions.get_by_id(biomass.id)
            biomass_m.lower_bound = fixed_growth
            biomass_m.upper_bound = fixed_growth
            m.objective = phb.id
            phb_flux = m.slim_optimize()
            phb_flux = phb_flux if phb_flux and phb_flux > 1e-9 else 0.0

    return max_growth, phb_flux


def describe_reaction(rxn, label=""):
    """Print a one/two-line human-readable description of a reaction."""
    if rxn is None:
        print(f"{label}: NOT FOUND")
        return
    print(f"{label}: {rxn.id} | {rxn.name} | bounds=({rxn.lower_bound}, {rxn.upper_bound})")
    print(f"    {rxn.reaction}")


def capped_figsize(n_items, per_item=0.15, min_size=8.0, max_size=18.0):
    """
    Compute a figure dimension that grows with the number of items being
    plotted (e.g. reactions on an axis) but is HARD-CAPPED at max_size.

    This exists because a pattern like `figsize=(0.3 * len(reaction_list), ...)`
    is dangerous on genome-scale models: if a reaction list unexpectedly
    contains hundreds of reactions (e.g. a pathway-distance neighborhood
    that turned out larger than expected), matplotlib will try to
    allocate a canvas hundreds of inches wide at ~200 dpi, which can
    require tens of gigabytes of RAM and crash the process. Always route
    reaction/metabolite-count-dependent figure sizes through this
    function instead of computing them inline.
    """
    size = max(min_size, per_item * n_items)
    return min(size, max_size)


def cap_item_list(items, max_n, label="items"):
    """
    Truncate a list to at most max_n entries, printing a warning if
    truncation occurred. Use this to bound the size of anything that will
    later drive a figure size, a correlation matrix, or a combinatorial
    search — silently processing an unexpectedly large list is how a
    'should be quick' script turns into an OOM crash.
    """
    if len(items) > max_n:
        print(f"NOTE: {len(items)} {label} found, capping to the first {max_n} "
              f"to keep runtime/memory bo3unded. Adjust the relevant MAX_* "
              f"constant at the top of the script if you need more.")
        return items[:max_n]
    return items

# Highly-connected "currency" metabolites (protons, water, energy/redox carriers, common cofactors) are excluded when building the reaction adjacency graph — including them would make almost every reaction in the model "adjacent" to almost every other reaction via e.g. ATP/ADP, which destroys the usefulness of a pathway-distance metric.
CURRENCY_METABOLITE_BASE_IDS = {
    "h", "h2o", "o2", "co2", "pi", "ppi",
    "atp", "adp", "amp", "gtp", "gdp", "gmp", "utp", "udp", "ctp", "cdp",
    "nad", "nadh", "nadp", "nadph",
    "fad", "fadh2", "fadh",
    "coa", "accoa",
    "nh4", "so4", "so3",
    "q8", "q8h2", "mqn8", "mql8",
    "na1", "k", "fe2", "fe3", "mg2", "ca2", "cl", "zn2", "cu2", "mn2", "cobalt2",
}


def base_metabolite_id(met):
    """Strip the compartment suffix from a metabolite's ID using its
    actual assigned compartment, e.g. 'atp_c' with compartment 'c' -> 'atp'.
    Falls back to the raw ID if the expected suffix isn't present."""
    suffix = f"_{met.compartment}"
    if met.id.endswith(suffix):
        return met.id[: -len(suffix)]
    return met.id


def is_currency_metabolite(met):
    return base_metabolite_id(met).lower() in CURRENCY_METABOLITE_BASE_IDS


def build_reaction_graph(model, exclude_currency=True):
    """
    Build an undirected reaction-adjacency graph: two reactions are
    connected if they share at least one metabolite (optionally excluding
    common currency metabolites, see CURRENCY_METABOLITE_BASE_IDS).

    Returns a dict: {reaction_id: set(neighbor_reaction_ids)}.
    """
    met_to_rxns = {}
    for rxn in model.reactions:
        if rxn.boundary:
            continue  # exchange/sink/demand reactions aren't pathway steps
        for met in rxn.metabolites:
            if exclude_currency and is_currency_metabolite(met):
                continue
            met_to_rxns.setdefault(met.id, set()).add(rxn.id)

    graph = {rxn.id: set() for rxn in model.reactions if not rxn.boundary}
    for met_id, rxn_ids in met_to_rxns.items():
        for r1 in rxn_ids:
            graph[r1].update(rxn_ids - {r1})

    return graph


def bfs_reaction_distance(graph, start_ids):
    """
    Breadth-first search over a reaction-adjacency graph (as built by
    build_reaction_graph) starting from one or more reaction IDs (distance 0).
    Returns a dict {reaction_id: distance}; reactions unreachable from any
    start node are omitted.
    """
    from collections import deque

    if isinstance(start_ids, str):
        start_ids = [start_ids]

    distances = {rid: 0 for rid in start_ids if rid in graph}
    queue = deque(distances.keys())
    while queue:
        current = queue.popleft()
        for neighbor in graph.get(current, ()):
            if neighbor not in distances:
                distances[neighbor] = distances[current] + 1
                queue.append(neighbor)
    return distances
