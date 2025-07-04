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
import skbio
from skbio import DistanceMatrix
from skbio.diversity import beta_diversity
from qiime2 import Metadata
import qiime2
import itertools



################# ALPHA DIVERSITY

########### functions for taxonomic alpha diversity


def alpha_taxa(table: pd.DataFrame, q: float) -> pd.Series:
    if isinstance(table, np.ndarray):
        table = pd.DataFrame(table)
    prop = table.div(table.sum(axis=1), axis=0)
    if q == 0:
        hill_numbers = (table > 0).sum(axis=1)  # Richness
    elif q == 1:
        hill_numbers = np.exp(
            -np.nansum(prop * np.log(prop.replace(0, np.nan)), axis=1)
        )  # Shannon
    else:
        hill_numbers = (prop ** q).sum(axis=1) ** (
            1 / (1 - q)
        )  # General formula
    # Return as Series
    return pd.Series(hill_numbers, index=table.index, name=f"TD q={q}")

########### functions for phylogenetic alpha diversity


def get_descendants(node, tips=False):
    if tips:
        return [leaf.name for leaf in node.traverse() if leaf.is_tip()]
    return [leaf.name for leaf in node.traverse()]


def get_node_tips(phylogeny):
    return [get_descendants(node, tips=True) for node in phylogeny.traverse()][
        1:
    ]


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
    if isinstance(abund, pd.DataFrame):
        return abund.div(abund.sum(axis=0), axis=1).replace(np.nan, 0)
    elif isinstance(abund, pd.Series):
        return abund / abund.sum() if abund.sum() != 0 else 0
    else:
        raise ValueError(
            "It must be a Feature Table."
        )


def tss2(pool_dict):
    total = sum(pool_dict.values())
    if total != 0:
        return {
            key: value / total for key, value in pool_dict.items()
        }


def dat_prep_phylo(comm, phylogeny):
    node_tips = get_node_tips(phylogeny)
    tips = list(phylogeny.tips())
    num_tips = len(tips)
    for i, tip in enumerate(tips, start=1):
        tip.number = i
    internal_nodes = traverse_internal(phylogeny)
    for i, node in enumerate(internal_nodes, start=num_tips + 1):
        node.number = i
    for node in phylogeny.traverse():
        node.id = node.number
    node_identifiers = [
        node.id if node.name is None else node.name
        for node in phylogeny.traverse()
    ]
    for new_name, node in zip(node_identifiers, phylogeny.traverse()):
        node.name = str(new_name)
    node_names = [node.name for node in phylogeny.traverse()][1:]
    M2 = pd.DataFrame(0, index=comm.columns, columns=node_names)
    for tips in range(len(M2.columns) - 1, -1, -1):
        tips_indices = node_tips[tips]
        for tip in tips_indices:
            if tip in M2.index:
                M2.loc[tip, M2.columns[tips]] = 1
    result = np.dot(comm.values, M2.values).T
    result_df = pd.DataFrame(result, index=M2.columns, columns=comm.T.columns)
    return result_df


def alpha_phylo_hillr(
    table: pd.DataFrame, phylogeny: TreeNode, q: float
) -> pd.Series:
    if (table < 0).any().any():
        raise ValueError("Table has negative values")
    if isinstance(table, np.ndarray):
        table = pd.DataFrame(table)
    comm = table.copy()
    comm = comm.loc[(comm.sum(axis=1) > 0), (comm.sum(axis=0) > 0)]
    species_in_tree = {leaf.name for leaf in phylogeny.tips()}
    species_in_comm = set(comm.columns)
    missing_species = species_in_comm - species_in_tree
    comm = comm.drop(columns=missing_species, errors="ignore")
    tree_leaf_order = [leaf.name for leaf in phylogeny.tips()]
    comm = comm[tree_leaf_order]
    if comm.shape[1] < 2:
        raise ValueError("Community has less than 2 species.")
    if q == 0:
        comm[comm > 0] = 1
    comm = comm.div(comm.sum(axis=1), axis=0).fillna(0)
    pabun = dat_prep_phylo(comm, phylogeny)
    plength = np.array([node.length for node in phylogeny.traverse()][1:])
    PD = {}
    D_t = {}
    for sample in pabun.columns:
        pabun_sample = pabun.loc[:, sample].values
        TT = np.sum(pabun_sample * plength)
        Item = np.where(pabun_sample > 0)[0]
        if len(Item) == 0:
            PD[sample] = np.nan
            D_t[sample] = np.nan
            continue
        D_t[sample] = TT
        if q == 1:
            PD[sample] = np.exp(
                -np.sum(
                    plength[Item]
                    * (pabun_sample[Item] / TT)
                    * np.log(pabun_sample[Item] / TT)
                )
            )
        else:
            PD[sample] = np.sum(
                plength[Item] * (pabun_sample[Item] / TT) ** q
            ) ** (1 / (1 - q))
    return pd.Series(PD, name=f"PD q={q}")


def alpha_phylo_hilldiv2(
    matrix: pd.DataFrame, phylogeny: TreeNode, q: float
) -> pd.Series:
    matrix = matrix.T
    Li = np.array([node.length for node in phylogeny.traverse()][1:])
    ltips = get_node_tips(phylogeny)
    ltips_idx = [
        [matrix.index.get_loc(name) for name in group] for group in ltips
    ]
    ai_multi = np.array(
        [
            [
                np.sum(tss(matrix.iloc[:, j]).iloc[ltips_idx_i])
                for ltips_idx_i in ltips_idx
            ]
            for j in range(matrix.shape[1])
        ]
    ).T
    ai_multi_df = pd.DataFrame(
        ai_multi,
        columns=matrix.columns,
        index=[f"tip_group_{i+1}" for i in range(len(ltips))],
    )
    ai = ai_multi_df
    pool = (
        pd.DataFrame(matrix)
        .apply(lambda row: 1 if any(row != 0) else 0, axis=1)
        .values
    )
    pool_dict = dict(zip(matrix.index, pool))
    normalized_pool_dict = tss2(pool_dict)
    for TipVector in ltips:
        for tip in TipVector:
            if tip not in normalized_pool_dict:
                print(f"Error: '{tip}' is not in normalized_pool_dict")
    ai_all = [
        sum(normalized_pool_dict[tip] for tip in TipVector)
        for TipVector in ltips
    ]
    T = np.sum(Li * ai_all)
    results = np.zeros(ai_multi.shape[1])
    for i in range(ai_multi.shape[1]):
        ai = ai_multi[:, i]
        if q == 1:
            results[i] = (
                np.exp(
                    -np.sum(
                        Li[ai != 0]
                        * (ai[ai != 0] / T)
                        * np.log(ai[ai != 0] / T)
                    )
                )
                / T
            )
        else:
            results[i] = (np.sum(Li[ai != 0] * (ai[ai != 0] / T) ** q)) ** (
                1 / (1 - q)
            ) / T
    return pd.Series(results, index=matrix.columns, name=f"PD q={q}")


def alpha_phylo(
    table: pd.DataFrame, phylogeny: TreeNode, q: float, metric: str = "PD"
) -> pd.Series:
    if metric == "PD":
        return alpha_phylo_hillr(table, phylogeny, q)
    elif metric == "qDT":
        return alpha_phylo_hilldiv2(table, phylogeny, q)
    else:
        raise ValueError("Invalid metric. Use 'PD' or 'qDT'.")

########## functions for functional alpha diversity


def compute_distance(
    traits: Metadata, dist: str = "euclidean"
) -> DistanceMatrix:
    traits_df = (
        traits if isinstance(traits, pd.DataFrame) else traits.to_dataframe()
    )
    if dist in ["euclidean", "manhattan"]:
        dist_matrix = beta_diversity(
            dist, traits_df.values, ids=traits_df.index
        )
    elif dist == "gower":
        # Manual Gower distance
        def gower_distance(data):
            def normalize_column(col):
                return (col - col.min()) / (col.max() - col.min())

            data_normalized = data.copy()
            for col in data.columns:
                if data[col].dtype != "object":
                    data_normalized[col] = normalize_column(data[col])
            n_species = len(data)
            n_features = len(data.columns)
            dist_matrix = np.zeros((n_species, n_species))
            for i in range(n_species):
                for j in range(i + 1, n_species):
                    total_distance = 0
                    for k in range(n_features):
                        if data.iloc[:, k].dtype != "object":
                            distance = abs(
                                data_normalized.iloc[i, k]
                                - data_normalized.iloc[j, k]
                            )
                        else:
                            distance = (
                                0 if data.iloc[i, k] == data.iloc[j, k] else 1
                            )
                        total_distance += distance
                    dist_matrix[i, j] = dist_matrix[j, i] = (
                        total_distance / n_features
                    )
            return dist_matrix

        gower_matrix = gower_distance(traits)
        dist_matrix = DistanceMatrix(gower_matrix, ids=traits.index)
    else:
        raise ValueError(
            "Parameter --p-dist must be 'euclidean', 'manhattan' or 'gower'."
        )
    return dist_matrix


def compute_tau(
    dist_matrix: skbio.DistanceMatrix, tau: str or float = "mean"
) -> float:
    if isinstance(tau, (int, float)) and tau > 0:
        return tau  # If it is a number, return it
    n = dist_matrix.shape[0]
    lower_values = dist_matrix.data[
        np.tril_indices(n, k=-1)
    ]  # extract lower values
    nonzero_lower_values = lower_values[lower_values != 0]  # exclude zeros
    if len(nonzero_lower_values) == 0:
        raise ValueError("Values are zero")
    if tau == "min":
        return np.min(nonzero_lower_values)
    elif tau == "max":
        return np.max(nonzero_lower_values)
    elif tau == "mean":
        return np.mean(nonzero_lower_values)
    else:
        raise ValueError(
            "Parameter tau must be a number > 0 or 'min', 'max' o 'mean'."
        )


def alpha_functional_hilldiv2(
    table: pd.DataFrame,
    traits: skbio.DistanceMatrix,
    q_value,
    tau: str or float = "mean",
) -> pd.Series:
    # Check common species between table and distance matrix
    table = table.copy().T
    common_species = table.index.intersection(traits.ids)
    if common_species.empty:
        raise ValueError("Species of table and traits do not match")
    # Determine tau
    tau = (
        compute_tau(traits, tau)
        if not (isinstance(tau, (int, float)) and tau > 0)
        else tau
    )
    # Extract distance matrix and limit to tau
    dij = traits.filter(common_species, strict=False).data
    dij = np.where(dij > tau, tau, dij)
    # Normalization
    comm_norm = table.div(table.sum(axis=0), axis=1).fillna(0).values
    species_names = table.index
    sample_names = table.columns.tolist()
    comm_norm = pd.DataFrame(
        comm_norm, index=species_names, columns=sample_names
    )
    a_multi = np.dot(1 - dij / tau, comm_norm)
    a_multi[a_multi < 0] = 0
    valid_rows = np.sum(a_multi, axis=1) != 0
    a_multi = a_multi[valid_rows, :]
    comm_norm = comm_norm.iloc[valid_rows, :]
    if a_multi.shape[0] == 0 or a_multi.shape[1] == 0:
        raise ValueError("All rows were removed after a-multi filter")
    v_multi = comm_norm / (a_multi + 1e-10)
    results = {}
    for j, sample in enumerate(v_multi.columns):
        v_col = v_multi[sample].values
        a_col = a_multi[:, j]
        mask = ~np.isnan(a_col)
        v_filtered = v_col[mask]
        a_filtered = a_col[mask]
        if q_value == 1:
            # For q=1:
            entropy = np.sum(
                -v_filtered * a_filtered * np.log(a_filtered + 1e-10)
            )
            results[sample] = np.exp(entropy)
        else:
            results[sample] = (np.sum(v_filtered * a_filtered ** q_value)) ** (
                1 / (1 - q_value)
            )
    results_series = pd.Series(results, name=f"qDT q={q_value}")
    return results_series


def fdisp(dij, a, tol=1e-07):
    n = dij.shape[0]
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(i):
            A[i, j] = -0.5 * (dij[i, j] ** 2)
            A[j, i] = A[i, j]
    row_means = np.mean(A, axis=1, keepdims=True)
    col_means = np.mean(A, axis=0, keepdims=True)
    grand_mean = np.mean(A)
    G = A - row_means - col_means + grand_mean
    eigvals, eigvecs = np.linalg.eigh(G)
    eigvals = eigvals[::-1]
    eigvecs = eigvecs[:, ::-1]
    w0 = eigvals[-1] / eigvals[0]
    if w0 > -tol:
        r = np.sum(eigvals > (eigvals[0] * tol))
    else:
        r = len(eigvals)
    eigvals_r = eigvals[:r]
    eigvecs_r = eigvecs[:, :r]
    vectors = eigvecs_r * np.sqrt(np.abs(eigvals_r))[None, :]
    pos = eigvals_r > 0
    n_sites = a.shape[0]
    FDis = np.zeros(n_sites)
    for i in range(n_sites):
        pres = np.where(a[i, :] > 0)[0]
        if pres.size >= 2:
            vec = vectors[pres, :]
            vec_unique = np.unique(vec, axis=0)
            if vec_unique.shape[0] >= 2:
                w = a[i, pres]
                centroid = np.average(vec, axis=0, weights=w)
                if np.any(pos):
                    diff_pos = vec[:, pos] - centroid[pos]
                    dist_pos = np.sum(diff_pos ** 2, axis=1)
                else:
                    dist_pos = np.zeros(vec.shape[0])
                if np.any(~pos):
                    diff_neg = vec[:, ~pos] - centroid[~pos]
                    dist_neg = np.sum(diff_neg ** 2, axis=1)
                else:
                    dist_neg = np.zeros(vec.shape[0])
                zij = np.sqrt(np.abs(dist_pos - dist_neg))
                FDis[i] = np.average(zij, weights=w)
            else:
                FDis[i] = 0
        else:
            FDis[i] = 0
    return {"FDis": FDis, "eig": eigvals_r, "vectors": vectors}


def alpha_functional_hillr(
    comm: pd.DataFrame,
    traits: skbio.DistanceMatrix,
    q: float,
    metric: str = "FD_q",
    base: float = np.e,
    check_data: bool = True,
    div_by_sp: bool = False,
    fdis_calc: bool = True,
    stand_dij: bool = False,
) -> pd.Series:
    common_species = comm.columns.intersection(traits.ids)
    if common_species.empty:
        raise ValueError("Species of community and traits do not match")
    comm = comm.loc[:, common_species]
    traits = traits.filter(common_species, strict=False)
    dij = traits.data.copy()
    if stand_dij:
        dij = dij / np.nanmax(dij)
    comm_np = comm.values.astype(float)
    row_sums = np.sum(comm_np, axis=1, keepdims=True)
    comm_rel = comm_np / row_sums
    SR = np.sum(comm_np > 0, axis=1)
    if q == 0:
        comm_np[comm_np > 0] = 1
        row_sums = np.sum(comm_np, axis=1, keepdims=True)
        comm_rel = comm_np / row_sums
    N = comm_np.shape[0]
    Q = np.zeros(N)
    for k in range(N):
        Q[k] = np.dot(comm_rel[k, :], np.dot(dij, comm_rel[k, :]))
    D_q = np.zeros(N)
    MD_q = np.zeros(N)
    FD_q = np.zeros(N)
    for k in range(N):
        idx = np.where(comm_rel[k, :] > 0)[0]
        df2 = comm_rel[k, idx]
        dis2 = dij[np.ix_(idx, idx)]
        if Q[k] == 0:
            D_q[k] = 0
        else:
            if q == 0:
                D_q[k] = (np.sum(dis2 / Q[k])) ** 0.5
            elif q == 1:
                eps = 1e-12
                outer_prod = np.outer(df2, df2)
                D_q[k] = np.exp(
                    -0.5
                    * np.sum(
                        (dis2 / Q[k]) * outer_prod * np.log(outer_prod + eps)
                    )
                    / np.log(base)
                )
            else:
                din = np.sum((dis2 / Q[k]) * (np.outer(df2, df2) ** q))
                D_q[k] = 0 if din == 0 else din ** (1 / (2 * (1 - q)))
        MD_q[k] = D_q[k] * Q[k]
        FD_q[k] = (D_q[k] ** 2) * Q[k]
    if div_by_sp:
        choose_val = np.where(SR >= 2, (SR * (SR - 1)) / 2, np.nan)
        D_q = D_q / SR
        MD_q = MD_q / SR
        FD_q = FD_q / choose_val
    if metric.upper() == "Q":
        results = Q
    elif metric.upper() == "FD_Q":
        results = FD_q
    elif metric.upper() == "D_Q":
        results = D_q
    elif metric.upper() == "MD_Q":
        results = MD_q
    elif metric.upper() == "FDIS":
        res_dict = fdisp(traits, comm.values, tol=1e-07)
        results = res_dict["FDis"]
    else:
        raise ValueError(
            "Invalid metric. Choose from: Q, FDis, D_q, MD_q, FD_q"
        )
    results_series = pd.Series(results, index=comm.index, name=f"FDis q={q}")
    return results_series


def alpha_functional(
    table: pd.DataFrame,
    traits: pd.DataFrame,
    q: float,
    dist: str = "euclidean",
    metric: str = "FD",
    tau: str or float = "mean",
) -> pd.Series:
    if isinstance(traits, qiime2.Metadata):
        traits = traits.to_dataframe()
    dm = compute_distance(traits, dist=dist)
    if metric == "FD":
        result = alpha_functional_hilldiv2(table, dm, q, tau=tau)
        result.name = f"FD q={q}"
    elif metric == "FD_q":
        result = alpha_functional_hillr(table, dm, q, metric="FD_q")
        result.name = f"FD_q q={q}"
    elif metric == "D_q":
        result = alpha_functional_hillr(table, dm, q, metric="D_q")
        result.name = f"D_q q={q}"
    elif metric == "Q":
        result = alpha_functional_hillr(table, dm, q, metric="Q")
        result.name = f"Q q={q}"
    elif metric == "MD_q":
        result = alpha_functional_hillr(table, dm, q, metric="MD_q")
        result.name = f"MD_q q={q}"
    elif metric == "FDis":
        result = alpha_functional_hillr(table, dm, q, metric="FDis")
        result.name = f"FDis q={q}"
    else:
        raise ValueError("Metric should be: FD, FD_q, D_q, Q, MD_q or FDis")
    return result


################### BETA DIVERSITY

########### functions for taxonomic beta diversity


def hilldiss_S(beta, N):
    return 1 - (
        ((1 / beta) - (1 / N)) / (1 - (1 / N))
    )


def hilldiss_C(beta, N, q):
    if q == 1:
        q_adj = 0.999999
    else:
        q_adj = q
    return 1 - (
        ((1 / beta) ** (q_adj - 1) - (1 / N) ** (q_adj - 1))
        / (1 - (1 / N) ** (q_adj - 1))
    )


def hilldiss_V(beta, N):
    return 1 - ((N - beta) / (N - 1))


def hilldiss_U(beta, N, q):
    if q == 1:
        q_adj = 0.999999
    else:
        q_adj = q
    return 1 - (
        ((1 / beta) ** (1 - q_adj) - (1 / N) ** (1 - q_adj))
        / (1 - (1 / N) ** (1 - q_adj))
    )


def hillpart_taxa(comm: pd.DataFrame, q: float):
    if isinstance(q, (int, float)):
        q = [q]

    if not isinstance(q, list):
        raise ValueError("The parameter q must be numeric")

    data = comm.loc[(comm != 0).any(axis=1)]
    pi = data.div(data.sum(axis=0), axis=1)

    results = []

    for q_value in q:
        if q_value == 0:
            species_per_sample = (pi > 0).sum(axis=0)
            alpha = species_per_sample.mean()
            gamma = (pi.sum(axis=1) > 0).sum()

        elif q_value == 1:
            pi_nonzero = pi[pi != 0]
            alpha = np.exp(
                -(pi_nonzero * np.log(pi_nonzero)).sum(axis=0).mean()
            )

            gamma_vec = pi.mean(axis=1)
            gamma = np.exp(
                -np.sum(
                    gamma_vec[gamma_vec != 0]
                    * np.log(gamma_vec[gamma_vec != 0])
                )
            )

        else:
            alpha_vals = (pi ** q_value).sum(axis=0)
            alpha = (alpha_vals.mean()) ** (1 / (1 - q_value))

            gamma_vec = pi.mean(axis=1)
            gamma = (gamma_vec ** q_value).sum() ** (1 / (1 - q_value))

        beta = gamma / alpha
        results.append([alpha, gamma, beta])

    result_df = pd.DataFrame(
        results,
        columns=["alpha", "gamma", "beta"],
        index=[f"q{q_value}" for q_value in q],
    )

    return result_df

def hilldiss_taxa(comm: pd.DataFrame, q: float, metric: str = "C", similarity: bool = False):
    if not isinstance(q, (int, float)):
        raise ValueError("q parameter must be numeric (int or float).")

    data = comm
    N = data.shape[1]

    if isinstance(q, list):
        q = q[0]

    df_hill = hillpart_taxa(data, q)
    beta = df_hill["beta"].values[0]

    if metric == "S":
        diss = hilldiss_S(beta, N)
    elif metric == "C":
        diss = hilldiss_C(beta, N, q)
    elif metric == "V":
        diss = hilldiss_V(beta, N)
    elif metric == "U":
        diss = hilldiss_U(beta, N, q)
    else:
        raise ValueError(f"Invalid metric: {metric}. Valid options: 'C', 'U', 'V', 'S'.")

    if similarity:
        diss = 1 - diss

    return pd.DataFrame({f"{metric}qN": [diss]})

def beta_taxa(
    table: pd.DataFrame, q: float, metric: str = "C", similarity: bool = False
) -> DistanceMatrix:
    if not isinstance(table, pd.DataFrame):
        raise ValueError("Input table must be a FeatureTable")

    table = table.T
    sample_names = table.columns.tolist()
    pairs = list(itertools.combinations(sample_names, 2))

    results = []
    for s1, s2 in pairs:
        sub_table = table[[s1, s2]]
        diss = hilldiss_taxa(sub_table, q=q, metric=metric, similarity=similarity)
        diss["site1"] = s1
        diss["site2"] = s2
        diss["q"] = q
        results.append(diss)

    df = pd.concat(results, ignore_index=True)
    metric_col = f"{metric}qN"

    if metric_col not in df.columns:
        raise ValueError(f"Invalid metric: '{metric_col}'")

    df_results = df[["q", "site1", "site2", metric_col]]

    dist_matrix = pd.DataFrame(index=sample_names, columns=sample_names, dtype=float)

    for _, row in df_results.iterrows():
        value = row[metric_col]
        dist_matrix.loc[row["site1"], row["site2"]] = value
        dist_matrix.loc[row["site2"], row["site1"]] = value
        dist_matrix.loc[row["site1"], row["site1"]] = 1 if similarity else 0
        dist_matrix.loc[row["site2"], row["site2"]] = 1 if similarity else 0

    dm = DistanceMatrix(dist_matrix.values, ids=dist_matrix.index.tolist())
    return dm




