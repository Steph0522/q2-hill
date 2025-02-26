# ---------------------------------------------------------------------------- Copyright (c) 2024, Stephanie 
# Hereira-Pacheco.
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
    return pd.Series(hill_numbers, index=table.index, name=f"TD q={q}")


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

def tss(abund):
     # Si 'abund' es un DataFrame
    if isinstance(abund, pd.DataFrame):
        return abund.div(abund.sum(axis=0), axis=1).replace(np.nan, 0)
    # Si 'abund' es una Serie (vector)
    elif isinstance(abund, pd.Series):
        return abund / abund.sum() if abund.sum() != 0 else 0
    else:
        raise ValueError("El argumento debe ser un DataFrame o una Serie de pandas.")


def tss2(pool_dict):
    total = sum(pool_dict.values())  # Sumar todas las abundancias
    if total != 0:
        return {key: value / total for key, value in pool_dict.items()}  # Normalizar por la suma total


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




def alpha_phylo_hillr(table: pd.DataFrame, phylogeny: TreeNode, q: float) -> pd.Series:

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

    return pd.Series(PD, name=f"PD q={q}")

def alpha_phylo_hilldiv2(matrix: pd.DataFrame, phylogeny: TreeNode, q: float) -> pd.Series:
    matrix = matrix.T
    Li = np.array([node.length for node in phylogeny.traverse()][1:])
    ltips = get_node_tips(phylogeny)
    ltips_idx = [[matrix.index.get_loc(name) for name in group] for group in ltips]

    # ai.multi calculation
    # Aplicar la función tss y luego calcular 'ai.multi'
    ai_multi = np.array([
     [np.sum(tss(matrix.iloc[:, j]).iloc[ltips_idx_i]) for ltips_idx_i in ltips_idx]
        for j in range(matrix.shape[1])
    ]).T

    # Crear un DataFrame con los resultados
    ai_multi_df = pd.DataFrame(ai_multi, columns=matrix.columns, index=[f"tip_group_{i+1}" for i in range(len(ltips))])
    ai = ai_multi_df

    # Calculate reference T for an even distribution of all present MAGs
    pool = pd.DataFrame(matrix).apply(lambda row: 1 if any(row != 0) else 0, axis=1).values
    pool_dict = dict(zip(matrix.index, pool))

    #ai.all
    # Normalizamos pool_dict
    normalized_pool_dict = tss2(pool_dict)
    # Asegúrate de que las especies en ltips estén en el diccionario normalized_pool_dict
    for TipVector in ltips:
        for tip in TipVector:
            if tip not in normalized_pool_dict:
                print(f"Error: '{tip}' no está en el diccionario normalized_pool_dict")


    # Calcula ai.all para cada conjunto de especies en ltips
    ai_all = [sum(normalized_pool_dict[tip] for tip in TipVector) for TipVector in ltips] 
    T = np.sum(Li * ai_all)

    # Calculate present OTUs/ASVs/MAGs
    present = np.sum(matrix != 0, axis=0)


 # Calcular valores de Hill con un solo if-else
    results = np.zeros(ai_multi.shape[1])

    for i in range(ai_multi.shape[1]):
        ai = ai_multi[:, i]

        if q == 1:
            results[i] = np.exp(-np.sum(Li[ai != 0] * (ai[ai != 0] / T) * np.log(ai[ai != 0] / T))) / T
        else:
            results[i] = (np.sum(Li[ai != 0] * (ai[ai != 0] / T) ** q)) ** (1 / (1 - q)) / T

#    return results

    return pd.Series(results, index=matrix.columns, name=f"PD q={q}")



def alpha_phylo(table: pd.DataFrame, phylogeny: TreeNode, q: float, metric: str = "PD") -> pd.Series:
    """
    Calculate phylogenetic diversity using hillR (PD) or hillhiv2 (qDT).
    """
    if metric == "PD":
        return alpha_phylo_hillr(table, phylogeny, q)
    elif metric == "qDT":
        return alpha_phylo_hilldiv2(table, phylogeny, q)
    else:
        raise ValueError("Invalid metric. Use 'PD' or 'qDT'.")


