# Bundled from Bridges et al., "Bridgesetal-CCC"
# (https://github.com/miller-jensen-lab/Bridgesetal-CCC), released under the
# Unlicense (public domain) -- see LICENSE-THIRD-PARTY.md. Unmodified except
# where noted. If you use this pipeline, please cite that repository/paper
# (see README.md).

import numpy as np
import scanpy as sc
import itertools
# milopy / distinctipy / networkx / community_louvain are imported lazily by the
# functions that need them so that lighter consumers (e.g. permutation tests)
# can use this module without the full milopy+R stack installed.


def encode_replicates(adata, rep_dict):
    rep = [None] * adata.shape[0]
    for z in np.arange(adata.shape[0]):
        if np.sum(adata.obsm['hash'][z]) > 0:
            z_arg = np.argmax(adata.obsm['hash'][z])
            if adata.obs['sample'][z] == 'BD1':
                rep[z] = rep_dict[z_arg]
            else:
                rep[z] = 'R{}'.format(z_arg+1)
        else:
            rep[z] = 'R0'
    return rep


def build_samplerep(adata, sample_slot, replicate_slot):
    sample_rep = [None] * adata.shape[0]
    for r in np.arange(adata.shape[0]):
        sample_rep[r] = str(adata.obs[sample_slot][r]) + ' ' + str(adata.obs[replicate_slot][r])
    return sample_rep


def run_milo(adata):
    import milopy
    import milopy.core as milo
    sc.pp.neighbors(adata)
    milo.make_nhoods(adata)
    # Count cells from each sample (just cond, no rep) in each nhood
    adata.obs['rep_code'] = adata.obs['Condition'].cat.codes
    milo.count_nhoods(adata, sample_col="rep_code")
    # Test for differential abundance between conditions
    # need to convert to "continuous" encoding to describe time
    milo.DA_nhoods(adata, design="~ cond_continuous")
    # build graph for viz
    milopy.utils.build_nhood_graph(adata)

    # annotation of nhoods by joint celltypes of interacting pairs
    milopy.utils.annotate_nhoods(adata, anno_col='celltype.Joint')
    adata.uns['nhood_adata'].obs.loc[adata.uns['nhood_adata'].obs["nhood_annotation_frac"] < 0.5, "nhood_annotation"] = "Mixed"

    return adata


def group_nhoods(adata, min_connect, max_difflfc):
    adj_nhood = np.zeros((adata.obsm['nhoods'].shape[1], adata.obsm['nhoods'].shape[1]))

    # only considering single cells belonging to more than one
    overlap_ind = np.where(np.sum(adata.obsm['nhoods'], axis=1) > 1)[0]
    for g in overlap_ind:
        nhood_ind = np.where(adata.obsm['nhoods'][g, :].todense() == 1)[1]
        ij = list(itertools.permutations(nhood_ind, 2))
        for q in ij:
            adj_nhood[q] = adj_nhood[q] + 1

    # still need to filter adj matrix entries to zero by connections (< 3) and LFC match (diff > 0.25?)
    nonzero_ind = np.where(adj_nhood > 0)
    logFC = adata.uns['nhood_adata'].obs['logFC']
    for f in np.arange(len(nonzero_ind[0])):
        if adj_nhood[nonzero_ind[0][f], nonzero_ind[1][f]] < min_connect or abs(logFC[nonzero_ind[0][f]] - logFC[nonzero_ind[1][f]]) > max_difflfc:
            adj_nhood[nonzero_ind[0][f], nonzero_ind[1][f]] = 0

    return adj_nhood


def cluster_nhoods(adata, min_connect, max_difflfc):
    import networkx as nx
    from community import community_louvain
    test_adj = group_nhoods(adata, min_connect, max_difflfc)
    G_test = nx.from_numpy_matrix(test_adj)
    partition2 = community_louvain.best_partition(G_test)
    print(np.max(list(partition2.values())))
    return partition2


def plot_nhood_clusters(adata, cluster_labels, title, alpha=0.1, min_size=10, plot_edges=False):
    from distinctipy import distinctipy
    nhood_adata = adata.uns["nhood_adata"].copy()

    nhood_adata.obs["graph_color"] = cluster_labels
    nhood_adata.obs["graph_color"] = nhood_adata.obs["graph_color"].astype('category')

    clust_col = distinctipy.get_colors(len(np.unique(nhood_adata.obs["graph_color"])))
    clust_pal = {np.unique(nhood_adata.obs["graph_color"])[i]: clust_col[i] for i in range(len(clust_col))}

    nhood_adata.obs.loc[nhood_adata.obs["SpatialFDR"] > alpha, "graph_color"] = np.nan

    # plotting order
    ordered = nhood_adata.obs.sort_values('SpatialFDR', na_position='last').index[::-1]
    nhood_adata = nhood_adata[ordered]

    sc.pl.embedding(nhood_adata, "X_milo_graph",
                    color="graph_color", palette=clust_pal,
                    size=adata.uns["nhood_adata"].obs["Nhood_size"] * min_size,
                    edges=plot_edges, neighbors_key="nhood",
                    frameon=False,
                    title=title
                    )

    return nhood_adata.obs["graph_color"], clust_pal


def plot_durable_clusters(adata, cluster_labels, title, alpha=0.1, beta=0.5, min_size=10, plot_edges=False):
    from distinctipy import distinctipy
    nhood_adata = adata.uns["nhood_adata"].copy()

    nhood_adata.obs["graph_color"] = cluster_labels
    nhood_adata.obs["graph_color"] = nhood_adata.obs["graph_color"].astype('category')

    clust_col = distinctipy.get_colors(len(np.unique(nhood_adata.obs["graph_color"])))
    clust_pal = {np.unique(nhood_adata.obs["graph_color"])[i]: clust_col[i] for i in range(len(clust_col))}

    nhood_adata.obs.loc[nhood_adata.obs["SpatialFDR"] < alpha, "graph_color"] = np.nan
    nhood_adata.obs.loc[nhood_adata.obs["logFC"] > beta, "graph_color"] = np.nan
    nhood_adata.obs.loc[nhood_adata.obs["logFC"] < -beta, "graph_color"] = np.nan

    sc.pl.embedding(nhood_adata, "X_milo_graph",
                    color="graph_color", palette=clust_pal,
                    size=adata.uns["nhood_adata"].obs["Nhood_size"] * min_size,
                    edges=plot_edges, neighbors_key="nhood",
                    frameon=False,
                    title=title
                    )

    return nhood_adata.obs["graph_color"], clust_pal


def get_sc_louvain(adata, cluster_slot='louvain'):
    louvain_onehot = np.zeros((adata.uns['nhood_adata'].obs[cluster_slot].shape[0], (np.unique(adata.uns['nhood_adata'].obs[cluster_slot])[-2] + 1).astype('int')))
    for c in adata.uns['nhood_adata'].obs[cluster_slot].index:
        if adata.uns['nhood_adata'].obs[cluster_slot][c] < float('inf'):
            louvain_onehot[int(c), adata.uns['nhood_adata'].obs[cluster_slot][c].astype('int')] = 1

    # get single-cell louvain neighborhood cluster labels
    sc_onehot = adata.obsm['nhoods']*louvain_onehot
    sc_louvain = np.zeros(sc_onehot.shape[0])
    for t in np.arange(sc_onehot.shape[0]):
        if np.sum(sc_onehot[t, :]) == 0:
            sc_louvain[t] = -1
        else:
            sc_louvain[t] = np.argmax(sc_onehot[t, :])

    return sc_louvain.astype('int')


def highlight_ind(clust, adata):
    highlight_ind_ = []
    for g in clust:
        i = np.where(adata.obs['louvain_str'] == g)[0]
        highlight_ind_.append(i)
    adata_highlight = adata[np.array(list(itertools.chain(*highlight_ind_)))]
    return adata_highlight


def highlight_NICHEScluster(niches_adata, adata, cluster_no):
    codes_1 = niches_adata[niches_adata.obs['sc_louvain'] == cluster_no].obs.index

    v = np.zeros(adata.shape[0])
    w = np.zeros(adata.shape[0])

    for n in codes_1:
        ind = np.where(adata.obs.index.str.contains("-".join(n.split("-", 2)[:1])))
        if len(ind[0]) > 0:
            v[ind[0]] = 1

    for n in codes_1:
        ind = np.where(adata.obs.index.str.contains("-".join(n.split('—')[1].split("-", 2)[:1])))
        if len(ind[0]) > 0:
            w[ind[0]] = 1

    adata.obs['cluster{}_sending'.format(cluster_no)] = v
    adata.obs['cluster{}_receiving'.format(cluster_no)] = w
    highlight_map = {0: 'Other', 1: 'Highlight'}
    adata.obs['cluster{}_sending'.format(cluster_no)] = adata.obs['cluster{}_sending'.format(cluster_no)].map(highlight_map)
    adata.obs['cluster{}_receiving'.format(cluster_no)] = adata.obs['cluster{}_receiving'.format(cluster_no)].map(highlight_map)

    return adata


def permutation_loop_test(niches_adata, sc_adata,
                          cluster_A, role_A,
                          cluster_B, role_B,
                          celltype_background,
                          stratify_sample=None,
                          n_perm=10000,
                          random_state=0):
    """VectorType-matched permutation test for NICHES reciprocal-loop overlap.

    Tests whether the unique cells participating in cluster A (as sender or
    receiver) overlap with the unique cells in cluster B more than expected
    under a null that resamples NICHES pairs from a VectorType-matched pool.

    Parameters
    ----------
    niches_adata : AnnData
        NICHES + Milo network with 'sc_louvain' and 'VectorType' in .obs.
    sc_adata : AnnData
        Single-cell AnnData with 'grouping' (or equivalent) in .obs to define
        the celltype-restricted background. Obs index is the cell barcode used
        by NICHES (16-mer barcode + sample suffix).
    cluster_A, cluster_B : int
        Internal sc_louvain cluster ids being compared.
    role_A, role_B : {'sender', 'receiver'}
        Which side of each cluster's pairs to extract cells from.
    celltype_background : pandas.Series
        Boolean mask indexed like sc_adata.obs (e.g. grouping == 'T cell').
        The overlap is computed within this subset; defines the relevant pool
        of cells for the null draw.
    stratify_sample : str or None
        If given (e.g. 'BD6'), restricts both the NICHES VT pool and the
        celltype background to that sample. Use this when the clusters were
        selected for enrichment in one condition; use None for clusters that
        span samples by construction (e.g. durable d8/d10 clusters).
    n_perm : int
        Number of permutations.
    random_state : int
        Seed for the random number generator.

    Returns
    -------
    dict with: observed, n_A_cells, n_B_cells, n_background,
               null_mean, null_ci95, fold, empirical_p, fisher_p, null_dist,
               vt_A, vt_B, pool_A_size, pool_B_size.
    """
    from scipy.stats import fisher_exact

    rng = np.random.default_rng(random_state)

    if stratify_sample is not None:
        niches_view = niches_adata[niches_adata.obs['sample'] == stratify_sample]
        sc_view = sc_adata[sc_adata.obs['sample'] == stratify_sample]
        bg_mask = celltype_background.loc[sc_view.obs.index]
    else:
        niches_view = niches_adata
        sc_view = sc_adata
        bg_mask = celltype_background

    # Map 16-mer sequence -> set of sc obs_names
    seq_to_obs = {}
    for obs_name in sc_view.obs.index:
        seq_to_obs.setdefault(obs_name.split('-')[0], set()).add(obs_name)

    def cells_in(pair_names, role):
        idx = 0 if role == 'sender' else 1
        out = set()
        for n in pair_names:
            seq = n.split('—')[idx].split('-')[0]
            if seq in seq_to_obs:
                out.update(seq_to_obs[seq])
        return out

    pairs_A = niches_view.obs.index[niches_view.obs['sc_louvain'] == cluster_A].tolist()
    pairs_B = niches_view.obs.index[niches_view.obs['sc_louvain'] == cluster_B].tolist()
    vt_A = niches_view.obs.loc[pairs_A, 'VectorType'].mode().iloc[0]
    vt_B = niches_view.obs.loc[pairs_B, 'VectorType'].mode().iloc[0]
    pool_A = niches_view.obs.index[niches_view.obs['VectorType'] == vt_A].tolist()
    pool_B = niches_view.obs.index[niches_view.obs['VectorType'] == vt_B].tolist()
    nA, nB = len(pairs_A), len(pairs_B)

    bg = set(sc_view.obs.index[bg_mask])
    obs_A_cells = cells_in(pairs_A, role_A) & bg
    obs_B_cells = cells_in(pairs_B, role_B) & bg
    observed = len(obs_A_cells & obs_B_cells)

    null = np.empty(n_perm, dtype=int)
    for i in range(n_perm):
        samp_A = rng.choice(pool_A, size=nA, replace=False)
        samp_B = rng.choice(pool_B, size=nB, replace=False)
        a_cells = cells_in(samp_A, role_A) & bg
        b_cells = cells_in(samp_B, role_B) & bg
        null[i] = len(a_cells & b_cells)

    null_mean = float(null.mean())
    null_ci = tuple(np.percentile(null, [2.5, 97.5]).tolist())
    fold = observed / null_mean if null_mean > 0 else float('inf')
    p_emp = (np.sum(null >= observed) + 1) / (n_perm + 1)

    # Fisher exact for direct comparison
    n_bg = len(bg)
    a = observed
    b = len(obs_A_cells) - a
    c = len(obs_B_cells) - a
    d = n_bg - a - b - c
    _, fisher_p = fisher_exact([[a, b], [c, d]], alternative='greater')

    return {
        'observed': int(observed),
        'n_A_cells': len(obs_A_cells),
        'n_B_cells': len(obs_B_cells),
        'n_background': n_bg,
        'null_mean': null_mean,
        'null_ci95': null_ci,
        'fold': float(fold),
        'empirical_p': float(p_emp),
        'fisher_p': float(fisher_p),
        'null_dist': null,
        'vt_A': vt_A,
        'vt_B': vt_B,
        'pool_A_size': len(pool_A),
        'pool_B_size': len(pool_B),
    }


def nhood_expression_mapping(adata, gene_oi):
    avg_expr = [np.mean(adata[adata.obs['sc_louvain'] == b][:, gene_oi].X).tolist() for b in np.unique(adata.obs['sc_louvain'])]
    expr_nhood_map = {np.unique(adata.obs['sc_louvain'])[c]: avg_expr[c] for c in np.arange(len(np.unique(adata.obs['sc_louvain'])))}

    return expr_nhood_map

