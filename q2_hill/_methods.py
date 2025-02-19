# ----------------------------------------------------------------------------
# Copyright (c) 2024, Stephanie Hereira-Pacheco.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import numpy as np
from skbio import TreeNode
import pandas as pd

def alpha_taxa(table: pd.DataFrame, q: float) -> pd.Series:
    # check if is a df
    if isinstance(table, np.ndarray):
        comm = pd.DataFrame(table)

    # abundances to props
    prop = table.div(table.sum(axis=1), axis=0)

    # Calculation of Hill numbers
    if q == 0:
        hill_numbers = (table > 0).sum(axis=1)  # Richness
    elif q == 1:
        hill_numbers = np.exp(-np.nansum(prop * np.log(prop.replace(0, np.nan)), axis=1))  # Shannon
    else:
        hill_numbers = (prop**q).sum(axis=1) ** (1 / (1 - q))  # General formula

    # Return as serie
    return pd.Series(hill_numbers, index=table.index, name=f"q={q}")


# Functions for data prep phylo
def get_descendants(node, tips=False):
    if tips:
        return [leaf.name for leaf in node.traverse() if leaf.is_tip()]
    return [leaf.name for leaf in node.traverse()]

def get_node_tips(phylogeny):
    return [get_descendants(node, tips=True) for node in phylogeny.traverse()][1:]

def traverse_internal(node):
    internal_nodes = []
    for node in node.traverse():
        if not node.is_tip():
            internal_nodes.append(node)
    return internal_nodes

def collect_edges(node):
    edges = []
    for node in node.traverse():
        for child in node.children:
            edges.append((node.number, child.number))
    return edges

def dat_prep_phylo(comm, phylogeny):
    #node tips
    node_tips = get_node_tips(phylogeny)
    tips = list(phylogeny.tips())
    num_tips = len(tips)

    # Assign tip numbers
    for i, tip in enumerate(tips, start=1):
        tip.number = i

    # internal nodes
    internal_nodes = traverse_internal(phylogeny)

    # ids for internal nodes
    for i, node in enumerate(internal_nodes, start=num_tips + 1):
        node.number = i

    # edges
    edges = collect_edges(phylogeny)

    #  ids for tree
    for node in phylogeny.traverse():
        node.id = node.number

    # ids for edges
    node_identifiers = [node.id if node.name is None else node.name for node in phylogeny.traverse()]

    # Renanme edges
    for new_name, node in zip(node_identifiers, phylogeny.traverse()):
        node.name = str(new_name)

    # Creat matrix
    n_cols = phylogeny.count() - 1
    n_rows = comm.shape[1]
    node_names = [node.name for node in phylogeny.traverse()][1:]
    M2 = pd.DataFrame(0, index=comm.columns, columns=node_names)  

    tips_names = [tip.name for tip in phylogeny.tips()]
    for i in range(len(M2.columns) - 1, -1, -1):
        tips_indices = node_tips[i]
        for tip in tips_indices:
            if tip in M2.index:
                M2.loc[tip, M2.columns[i]] = 1

    result = np.dot(comm.values, M2.values).T
    result_df = pd.DataFrame(result, index=M2.columns, columns=comm.T.columns)
    return result_df




def alpha_phylo(table: pd.DataFrame, phylogeny: TreeNode, q: float) -> pd.Series:

    if (table < 0).any().any():
        raise ValueError("Table has negative values")

    if isinstance(table, np.ndarray):
          table = pd.DataFrame(table)

    # Remove empty samples and species
    comm = table.copy()
    comm = comm.loc[(comm.sum(axis=1) > 0), (comm.sum(axis=0) > 0)]

    # Check tree to right structre
    species_in_tree = {leaf.name for leaf in phylogeny.tips()}
    species_in_comm = set(comm.columns)
    common_species = species_in_comm.intersection(species_in_tree)

    # Remove species that are not in tree
    missing_species = species_in_comm - species_in_tree
    if missing_species:
        if show_warning:
            print(f"Warning: Removing {len(missing_species)} species that are not in tree.")
        comm = comm.drop(columns=missing_species, errors="ignore")

    # Order columns by tree 
    tree_leaf_order = [leaf.name for leaf in phylogeny.tips()]

    # Ordenar table
    comm = comm[tree_leaf_order]

    # Check at least 2 species
    if comm.shape[1] < 2:
        raise ValueError("Community has less than 2 species.")

    # Transform to incidence
    if q == 0:
        comm[comm > 0] = 1

    # Convert to props
    comm = comm.div(comm.sum(axis=1), axis=0).fillna(0)

    # Transform with `dat_prep_phylo`
    pabun = dat_prep_phylo(comm, phylogeny)

    # get branch lengths
    plength = np.array([node.length for node in phylogeny.traverse()][1:])

    # Calculate phylo diversity
    PD = {}
    D_t = {}
    for sample in pabun.columns:
        pabun_sample = pabun.loc[:, sample].values
        TT = np.sum(pabun_sample * plength)
        I = np.where(pabun_sample > 0)[0]

        if len(I) == 0:
            PD[sample] = np.nan
            D_t[sample] = np.nan
            continue

        D_t[sample] = TT

        if q == 1:
            PD[sample] = np.exp(-np.sum(plength[I] * (pabun_sample[I] / TT) * np.log(pabun_sample[I] / TT)))
        else:
            PD[sample] = np.sum(plength[I] * (pabun_sample[I] / TT) ** q) ** (1 / (1 - q))

    return pd.Series(PD, name="PD")

