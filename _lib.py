"""
_lib.py -- small helpers local to this pipeline, kept separate from
milo_helpers.py (bundled from Bridges et al. -- see that file's header) so
the original code there stays a clean, unmodified copy.

Six things live here, all workarounds for old/unmaintained dependencies
colliding with newer library versions -- not bugs in this pipeline's own
code, nor in milo_helpers.py's own code either (it was written against
older library versions that are no longer what `uv sync` resolves today):

1. write_wilcoxon_results(): a modern-pandas-safe Excel writer for
   scanpy's rank_genes_groups() results, one sheet per cluster. Needed
   because `pd.ExcelWriter(...).save()` was removed in pandas 2.0 (replaced
   by the writer's context-manager `__exit__` / `.close()`).

2. make_nhoods_fixed(): an exact copy of milopy.core.make_nhoods() (pinned
   commit 30646f5, see pyproject.toml) with one line fixed -- see its
   docstring below.

3. DA_nhoods_fixed(): an exact copy of milopy.core.DA_nhoods() (same pinned
   commit) with the rpy2 global-conversion setup fixed -- see its docstring
   below.

4. count_nhoods_fixed(): an exact copy of milopy.core.count_nhoods() (same
   pinned commit) with one line fixed -- see its docstring below.

5. cluster_nhoods_fixed(): an exact copy of milo_helpers.cluster_nhoods()
   with one line fixed -- see its docstring below. (This is the one
   function here patching milo_helpers.py's own code, not milopy's --
   still done as a copy in this file so milo_helpers.py itself stays
   unmodified.)

6. get_sc_louvain_fixed(): a corrected copy of milo_helpers.get_sc_louvain(),
   handling a real edge case (zero significant neighborhoods) the original
   crashes on -- see its docstring below.

Everything else in run_milo_da.py reuses milo_helpers.py directly, since
plot_nhood_clusters / get_sc_louvain / group_nhoods are already fully
generic and don't need changes.

The rest of this file (select_clusters, filter_clusters_by_celltype,
rank_clusters, dominant_celltype_labels, print_cluster_significance,
auto_select_lr_pairs, benjamini_hochberg, compute_comparison_stats,
compute_group_means, compute_log2fc) isn't a compatibility patch -- it's
shared plotting/analysis logic used by plot_lr_heatmap.py,
plot_lr_comparison_heatmap.py, plot_lr_comparison_sidebyside.py,
plot_milo_volcano.py, and plot_lr_log2fc_heatmap.py, factored out here so
the scripts don't duplicate it.
"""
import logging
import random
import re
import sys

import anndata
import numpy as np
import pandas as pd
import scipy.sparse


def write_wilcoxon_results(excel_path, adata, groups, de_key):
    """Write per-cluster scanpy rank_genes_groups results to one Excel sheet
    per cluster. Same output shape as utils.write_deres(), fixed for
    pandas>=2.0.
    """
    res_cat = ["names", "scores", "logfoldchanges", "pvals", "pvals_adj"]

    with pd.ExcelWriter(excel_path, engine="xlsxwriter") as writer:
        for g in groups:
            d = {cat: adata.uns[de_key][cat][str(g)].tolist() for cat in res_cat}
            df = pd.DataFrame(data=d)
            df.to_excel(writer, sheet_name=str(g))


def make_nhoods_fixed(adata, neighbors_key=None, prop=0.1, seed=42):
    """Exact copy of milopy.core.make_nhoods() (pinned commit 30646f5:
    https://github.com/emdann/milopy/blob/30646f538481151b6101b5e2f133858d2171000a/milopy/core.py),
    with ONE line changed.

    The original does:
        nhood_ixs = adata.obs["nhood_ixs_refined"] == 1
        dist_mat = knn_dists[nhood_ixs, :]
    i.e. indexes a scipy sparse matrix with a raw pandas boolean Series.
    scipy's sparse indexing calls `.nonzero()` on that index internally,
    which used to work because pandas Series had a passthrough `.nonzero()`
    method -- removed in a later pandas release. milopy hasn't been updated
    since. Fixed here by converting to a numpy array (`.values`) first,
    which scipy handles natively via numpy's own `.nonzero()`. Nothing else
    is changed from the original.
    """
    # Get reduced dim used for KNN graph
    if neighbors_key is None:
        try:
            use_rep = adata.uns["neighbors"]["params"]["use_rep"]
        except KeyError:
            logging.warning('Using X_pca as default embedding')
            use_rep = "X_pca"
        try:
            knn_graph = adata.obsp["connectivities"].copy()
        except KeyError:
            raise KeyError(
                'No "connectivities" slot in adata.obsp -- please run scanpy.pp.neighbors(adata) first'
            )
    else:
        try:
            use_rep = adata.uns[neighbors_key]["params"]["use_rep"]
        except KeyError:
            logging.warning('Using X_pca as default embedding')
            use_rep = "X_pca"
        knn_graph = adata.obsp[neighbors_key + "_connectivities"].copy()

    # Get reduced dim
    if use_rep == 'X':
        X_dimred = adata.X
        if scipy.sparse.issparse(X_dimred):
            X_dimred = X_dimred.A
    else:
        X_dimred = adata.obsm[use_rep]

    # Sample size
    n_ixs = int(np.round(adata.n_obs * prop))

    # Binarize
    knn_graph[knn_graph != 0] = 1

    #  Sample random vertices
    random.seed(seed)
    random_vertices = random.sample(range(adata.n_obs), k=n_ixs)
    random_vertices.sort()

    ixs_nn = knn_graph[random_vertices, :]

    # Refine sampling
    non_zero_rows = ixs_nn.nonzero()[0]
    non_zero_cols = ixs_nn.nonzero()[1]

    refined_vertices = np.empty(shape=[len(random_vertices), ])

    from sklearn.metrics.pairwise import euclidean_distances
    for i in range(len(random_vertices)):
        nh_pos = np.median(
            X_dimred[non_zero_cols[non_zero_rows == i], :], 0).reshape(-1, 1)
        nn_ixs = non_zero_cols[non_zero_rows == i]
        # Find closest real point (amongst nearest neighbors)
        dists = euclidean_distances(
            X_dimred[non_zero_cols[non_zero_rows == i], :], nh_pos.T)
        # Update vertex index
        refined_vertices[i] = nn_ixs[dists.argmin()]

    refined_vertices = np.unique(refined_vertices.astype("int"))
    refined_vertices.sort()

    nhoods = knn_graph[:, refined_vertices]
    adata.obsm['nhoods'] = nhoods

    # Add ixs to adata
    adata.obs["nhood_ixs_random"] = adata.obs_names.isin(
        adata.obs_names[random_vertices])
    adata.obs["nhood_ixs_refined"] = adata.obs_names.isin(
        adata.obs_names[refined_vertices])
    adata.obs["nhood_ixs_refined"] = adata.obs["nhood_ixs_refined"].astype(
        "int")
    adata.obs["nhood_ixs_random"] = adata.obs["nhood_ixs_random"].astype("int")
    # Store info on neighbor_key used
    adata.uns["nhood_neighbors_key"] = neighbors_key
    # Store distance to K-th nearest neighbor (used for spatial FDR correction)
    if neighbors_key is None:
        knn_dists = adata.obsp["distances"]
    else:
        knn_dists = adata.obsp[neighbors_key + "_distances"]
    # THE FIX: .values converts the pandas boolean Series to a numpy array
    # before using it to index the sparse matrix (see docstring above).
    nhood_ixs = (adata.obs["nhood_ixs_refined"] == 1).values
    dist_mat = knn_dists[nhood_ixs, :]
    k_distances = dist_mat.max(1).toarray().ravel()
    adata.obs["nhood_kth_distance"] = 0
    adata.obs.loc[adata.obs["nhood_ixs_refined"]
                  == 1, "nhood_kth_distance"] = k_distances


def DA_nhoods_fixed(adata, design, model_contrasts=None, subset_samples=None, add_intercept=True):
    """Exact copy of milopy.core.DA_nhoods() (pinned commit 30646f5:
    https://github.com/emdann/milopy/blob/30646f538481151b6101b5e2f133858d2171000a/milopy/core.py),
    with the rpy2 setup fixed.

    The original starts with:
        rpy2.robjects.numpy2ri.activate()
        rpy2.robjects.pandas2ri.activate()
    relying on rpy2's old *global* implicit-conversion mode for every R call
    that follows (passing numpy arrays / pandas DataFrames straight into R
    function calls). Newer rpy2 (the global-activation mechanism was
    deprecated, then hardened into a hard failure) makes `.activate()`
    itself raise `DeprecationWarning` as an exception instead of activating
    anything -- so the original function fails on its very first line, before
    doing any real work. Fixed here by using the modern, non-deprecated
    `rpy2.robjects.conversion.localconverter()` context manager instead,
    wrapped around the same body that needs numpy/pandas<->R conversion.
    Nothing else is changed from the original logic (including its use of
    milopy.core's own private `_try_import_bioc_library` and
    `_graph_spatialFDR` helpers, called here via the `milopy.core` module
    directly rather than copied).
    """
    import rpy2.robjects.numpy2ri
    import rpy2.robjects.pandas2ri
    from rpy2.robjects import conversion, default_converter
    from rpy2.robjects.packages import importr, STAP
    import milopy.core as milo

    # THE FIX: a scoped converter via localconverter(), instead of the
    # deprecated-and-now-broken global .activate() calls.
    conv = (default_converter
            + rpy2.robjects.numpy2ri.converter
            + rpy2.robjects.pandas2ri.converter)
    with conversion.localconverter(conv):
        # Loaded (importr attaches the package into R's own search path, like
        # library(edgeR)/library(limma)) so milo_da_fit()'s R code below can
        # call DGEList/calcNormFactors/estimateDisp/glmQLFit/glmQLFTest/
        # topTags by bare name -- these bindings aren't called directly from
        # Python anymore, see the "THE FIX, continued" comment further down.
        milo._try_import_bioc_library("edgeR")
        limma = milo._try_import_bioc_library("limma")
        stats = importr("stats")

        nhood_adata = adata.uns["nhood_adata"]
        covariates = [x.strip(" ") for x in set(
            re.split('\\+|\\*', design.lstrip("~ ")))]

        # Add covariates used for testing to nhood_adata.var
        sample_col = nhood_adata.uns["sample_col"]
        try:
            nhoods_var = adata.obs[covariates + [sample_col]].drop_duplicates()
        except KeyError:
            missing_cov = [
                x for x in covariates if x not in nhood_adata.var.columns]
            raise KeyError(
                'Covariates {c} are not columns in adata.obs'.format(
                    c=" ".join(missing_cov))
            )
        nhoods_var = nhoods_var[covariates + [sample_col]]
        nhoods_var.index = nhoods_var[sample_col].astype("str")

        try:
            assert nhoods_var.loc[nhood_adata.var_names].shape[0] == len(
                nhood_adata.var_names)
        except Exception:
            raise ValueError(
                "Covariates cannot be unambiguously assigned to each sample -- each sample value should match a single covariate value")
        nhood_adata.var = nhoods_var.loc[nhood_adata.var_names]
        # Get design dataframe
        try:
            design_df = nhood_adata.var[covariates]
        except KeyError:
            missing_cov = [
                x for x in covariates if x not in nhood_adata.var.columns]
            raise KeyError(
                'Covariates {c} are not columns in adata.uns["nhood_adata"].var'.format(
                    c=" ".join(missing_cov))
            )

        # Get count matrix
        count_mat = nhood_adata.X.toarray()
        lib_size = count_mat.sum(0)

        # Filter out samples with zero counts
        keep_smp = lib_size > 0

        # Subset samples
        if subset_samples is not None:
            keep_smp = keep_smp & nhood_adata.var_names.isin(subset_samples)

        design_df = design_df[keep_smp]
        for i, e in enumerate(design_df.columns):
            if design_df.dtypes[i].name == 'category':
                design_df[e] = design_df[e].cat.remove_unused_categories()

        # Filter out nhoods with zero counts (they can appear after sample filtering)
        keep_nhoods = count_mat[:, keep_smp].sum(1) > 0

        # Define model matrix
        if not add_intercept or model_contrasts is not None:
            design = design + ' + 0'
        model = stats.model_matrix(object=stats.formula(
            design), data=design_df)

        # Fit NB-GLM and test -- ALL done inside one R function (THE FIX,
        # continued): edgeR::DGEList() and glmQLFit() return R S4 objects
        # (DGEList / DGEGLM). Round-tripping those through Python as
        # intermediate variables (dge <- ...; dge <- calcNormFactors(dge);
        # fit <- glmQLFit(dge, ...), the way the original milopy code and
        # my first attempt at this fix both did) hits a real rpy2 S4
        # conversion bug: the returned object comes back into Python as a
        # generic rpy2.rlike.container.OrdDict instead of something that
        # converts cleanly back to R, so the *next* R call in the chain
        # fails (either with edgeR's own "colSums(x) : 'x' must be numeric
        # or complex", under some rpy2 versions, or a hard
        # NotImplementedError converting the OrdDict back to R, under
        # others -- same underlying cause, different symptom). Only plain
        # numeric matrices (counts, lib_size, model) and a final
        # data.frame result -- both well-supported, ordinary conversions --
        # ever cross the Python/R boundary here; the DGEList/fit objects
        # are created and consumed entirely inside R.
        _counts_for_r = count_mat[keep_nhoods, :][:, keep_smp]
        _libsize_for_r = lib_size[keep_smp]

        _da_r_code = '''
        milo_da_fit <- function(counts, lib_size, model, contrasts=NULL, coef=NULL) {
            dge <- DGEList(counts=counts, lib.size=lib_size)
            dge <- calcNormFactors(dge, method="TMM")
            dge <- estimateDisp(dge, model)
            fit <- glmQLFit(dge, model, robust=TRUE)
            if (!is.null(contrasts)) {
                qlf <- glmQLFTest(fit, contrast=contrasts)
            } else {
                qlf <- glmQLFTest(fit, coef=coef)
            }
            as.data.frame(topTags(qlf, sort.by='none', n=Inf))
        }
        '''
        milo_da_fit = STAP(_da_r_code, "milo_da_fit").milo_da_fit

        n_coef = model.shape[1]
        if model_contrasts is not None:
            r_str = '''
            get_model_cols <- function(design_df, design){
                m = model.matrix(object=formula(design), data=design_df)
                return(colnames(m))
            }
            '''
            get_model_cols = STAP(r_str, "get_model_cols")
            model_mat_cols = get_model_cols.get_model_cols(design_df, design)
            model_df = pd.DataFrame(model)
            model_df.columns = model_mat_cols
            try:
                mod_contrast = limma.makeContrasts(
                    contrasts=model_contrasts, levels=model_df)
            except Exception:
                raise ValueError(
                    "Model contrasts must be in the form 'A-B' or 'A+B'")
            res = milo_da_fit(_counts_for_r, _libsize_for_r, model, contrasts=mod_contrast)
        else:
            res = milo_da_fit(_counts_for_r, _libsize_for_r, model, coef=n_coef)
        res = conversion.rpy2py(res)
        if not isinstance(res, pd.DataFrame):
            res = pd.DataFrame(res)

    # Save outputs (plain pandas from here on -- no R conversion needed)
    res.index = nhood_adata.obs_names[keep_nhoods]
    if any([x in nhood_adata.obs.columns for x in res.columns]):
        nhood_adata.obs = nhood_adata.obs.drop(res.columns, axis=1)
    nhood_adata.obs = pd.concat([nhood_adata.obs, res], axis=1)

    # Run Graph spatial FDR correction (pure Python/pandas, reused as-is)
    milo._graph_spatialFDR(adata)


def count_nhoods_fixed(adata, sample_col):
    """Exact copy of milopy.core.count_nhoods() (pinned commit 30646f5:
    https://github.com/emdann/milopy/blob/30646f538481151b6101b5e2f133858d2171000a/milopy/core.py),
    with ONE line fixed.

    The original does:
        sample_dummies = pd.get_dummies(adata.obs[sample_col])
        ...
        sample_dummies = scipy.sparse.csr_matrix(sample_dummies.values)
        nhood_count_mat = adata.obsm["nhoods"].T.dot(sample_dummies)
    `pd.get_dummies()` used to default to `uint8` columns; current pandas
    defaults to `bool`. That bool matrix ends up as `nhood_adata.X`
    unchanged, and is later passed straight into `edgeR::DGEList(counts=...)`
    via rpy2 in DA_nhoods_fixed() -- where it stops surviving as something
    edgeR's internal `colSums()` accepts ("'x' must be numeric or complex").
    Fixed here by casting the dummy matrix to a real numeric dtype
    (`float64`) right after `get_dummies()`, before anything downstream
    ever sees it. Nothing else is changed from the original.
    """
    try:
        nhoods = adata.obsm["nhoods"]
    except KeyError:
        raise KeyError(
            'Cannot find "nhoods" slot in adata.obsm -- please run milopy.make_nhoods(adata)'
        )
    #  Make nhood abundance matrix
    # THE FIX: force a numeric dtype -- pd.get_dummies() now defaults to
    # bool, which doesn't survive the trip through scipy sparse + rpy2 as
    # something edgeR's R-side colSums() will accept.
    sample_dummies = pd.get_dummies(adata.obs[sample_col]).astype(np.float64)
    all_samples = sample_dummies.columns
    sample_dummies = scipy.sparse.csr_matrix(sample_dummies.values)
    nhood_count_mat = adata.obsm["nhoods"].T.dot(sample_dummies)
    nhood_var = pd.DataFrame(index=all_samples)
    nhood_adata = anndata.AnnData(X=nhood_count_mat, var=nhood_var)
    nhood_adata.uns["sample_col"] = sample_col
    # Save nhood index info
    nhood_adata.obs["index_cell"] = adata.obs_names[adata.obs["nhood_ixs_refined"] == 1]
    nhood_adata.obs["kth_distance"] = adata.obs.loc[adata.obs["nhood_ixs_refined"]
                                                    == 1, "nhood_kth_distance"].values
    adata.uns["nhood_adata"] = nhood_adata


def cluster_nhoods_fixed(adata, min_connect, max_difflfc):
    """Exact copy of milo_helpers.cluster_nhoods(), with ONE line fixed.

    The original does:
        G_test = nx.from_numpy_matrix(test_adj)
    `networkx.from_numpy_matrix` was renamed to `networkx.from_numpy_array`
    in networkx 3.0 (2023) and the old name was removed entirely in a later
    release -- current networkx raises `AttributeError: module 'networkx'
    has no attribute 'from_numpy_matrix'`. Fixed here by calling the new
    name instead; the adjacency-matrix semantics are identical, so nothing
    else about the clustering changes. Reuses milo_helpers.py's own
    `group_nhoods()` unchanged, since that helper doesn't call anything
    broken.
    """
    import networkx as nx
    from community import community_louvain
    from milo_helpers import group_nhoods

    test_adj = group_nhoods(adata, min_connect, max_difflfc)
    # THE FIX: from_numpy_matrix (removed) -> from_numpy_array (networkx >= 3.0).
    G_test = nx.from_numpy_array(test_adj)
    partition2 = community_louvain.best_partition(G_test)
    print(np.max(list(partition2.values())))
    return partition2


def get_sc_louvain_fixed(adata, cluster_slot='louvain'):
    """Corrected copy of milo_helpers.get_sc_louvain().

    The original computes the number of louvain clusters via:
        (np.unique(adata.uns['nhood_adata'].obs[cluster_slot])[-2] + 1)
    relying on `inf` (used elsewhere to mark non-significant neighborhoods)
    always sorting to the very last position of np.unique()'s output, so
    `[-2]` grabs the largest *finite* cluster label. This raises
    `IndexError: index -2 is out of bounds for axis 0 with size 1` whenever
    cluster_slot has only one unique value -- which happens whenever EVERY
    neighborhood is non-significant (SpatialFDR > alpha for all of them),
    since then cluster_slot is entirely `inf`.

    That's a real, valid Milo result -- "no differentially abundant
    neighborhoods at this --alpha" -- not a bug in the data, and not
    something that should crash the whole run. This version detects it
    directly (rather than via the fragile sort-position trick), returns
    -1 for every cell (consistent with the existing "-1 = not in any
    significant neighborhood" convention used for individual cells below),
    and prints a note so the situation is visible instead of passing
    silently. Everything else -- the one-hot / argmax assignment of each
    cell to its most-overlapping significant neighborhood -- is unchanged.
    """
    cluster_vals = adata.uns['nhood_adata'].obs[cluster_slot]
    finite_vals = cluster_vals[cluster_vals < float('inf')]
    if finite_vals.empty:
        print(
            f"[get_sc_louvain_fixed] No neighborhoods were significant "
            f"('{cluster_slot}' is entirely non-finite for all "
            f"{cluster_vals.shape[0]} neighborhoods) -- every cell is being "
            "marked -1 (unassigned). This means Milo found no "
            "differentially abundant neighborhoods at the --alpha "
            "threshold used; it is not a bug in this pipeline. Consider "
            "rerunning with a larger --alpha if this is unexpected."
        )
        return np.full(adata.n_obs, -1, dtype='int')

    # THE FIX: derive the cluster count directly from the finite labels,
    # instead of the fragile "second-to-last after sorting" trick.
    n_clusters = int(finite_vals.max()) + 1
    louvain_onehot = np.zeros((cluster_vals.shape[0], n_clusters))
    for c in cluster_vals.index:
        if cluster_vals[c] < float('inf'):
            louvain_onehot[int(c), int(cluster_vals[c])] = 1

    # get single-cell louvain neighborhood cluster labels
    sc_onehot = adata.obsm['nhoods'] * louvain_onehot
    sc_louvain = np.zeros(sc_onehot.shape[0])
    for t in np.arange(sc_onehot.shape[0]):
        if np.sum(sc_onehot[t, :]) == 0:
            sc_louvain[t] = -1
        else:
            sc_louvain[t] = np.argmax(sc_onehot[t, :])

    return sc_louvain.astype('int')


def rank_clusters(adata, clusters, groupby, rank_by):
    """Return `clusters` reordered best-first by `rank_by` (caller truncates).

    'size' = most cell pairs assigned, always available. 'spatialfdr' =
    most significant on average (lowest mean SpatialFDR among that
    cluster's own significant neighborhoods) -- needs
    adata.uns['nhood_adata'] with a 'louvain' column, i.e. the direct h5ad
    output of run_milo_da.py.
    """
    if rank_by == "size":
        counts = (
            adata.obs[adata.obs[groupby].astype(str).isin(clusters)][groupby]
            .astype(str).value_counts()
        )
        return sorted(clusters, key=lambda c: counts.get(c, 0), reverse=True)

    if "nhood_adata" not in adata.uns or "louvain" not in adata.uns["nhood_adata"].obs.columns:
        sys.exit(
            "--rank-clusters-by spatialfdr needs adata.uns['nhood_adata'] "
            "with a 'louvain' column -- this must be the direct h5ad output "
            "of run_milo_da.py (not a re-saved/subsetted copy). Use "
            "--rank-clusters-by size instead if that's not available."
        )
    nhood_obs = adata.uns["nhood_adata"].obs
    finite = nhood_obs[nhood_obs["louvain"] < float("inf")].copy()
    finite["louvain_str"] = finite["louvain"].astype(int).astype(str)
    finite = finite[finite["louvain_str"].isin(clusters)]
    mean_fdr = finite.groupby("louvain_str")["SpatialFDR"].mean()
    return sorted(clusters, key=lambda c: mean_fdr.get(c, float("inf")))


def select_clusters(adata, groupby, clusters_arg=None, top_n_clusters=None,
                     rank_clusters_by="size", include_unassigned=False,
                     celltype_filter=None, celltype_col="VectorType",
                     min_cluster_size=None):
    """Shared cluster-selection logic for all five cluster-based scripts
    (plot_lr_heatmap.py, plot_lr_comparison_heatmap.py,
    plot_lr_comparison_sidebyside.py, plot_milo_volcano.py,
    plot_lr_log2fc_heatmap.py): resolve --clusters / --top-n-clusters /
    --rank-clusters-by / --include-unassigned / --celltype-filter /
    --min-cluster-size into a concrete, sorted list of cluster labels
    (strings) to display. Exits with an explanatory message on bad/empty
    input.

    --celltype-filter and --min-cluster-size (if given) are both applied
    FIRST, before --top-n-clusters truncation and even when --clusters is
    set explicitly -- so "--celltype-filter Macrophage --top-n-clusters 10"
    means "the top 10 macrophage clusters", not "whichever of the
    (unfiltered) top 10 clusters happen to be macrophage ones" (those give
    different, and sometimes empty, results). Same reasoning for
    --min-cluster-size: it's a quality floor ("don't show me a cluster
    Milo called significant on 4 cell pairs"), not a scope choice, so it's
    enforced everywhere -- including against clusters named explicitly via
    --clusters, on the theory that if you type a tiny cluster's number by
    hand you probably didn't realize how small it was.
    """
    all_vals = adata.obs[groupby].astype(str).unique().tolist()

    if clusters_arg:
        clusters = [str(c) for c in clusters_arg]
        missing_cl = [c for c in clusters if c not in all_vals]
        if missing_cl:
            sys.exit(
                f"--clusters value(s) not found in '{groupby}': "
                f"{missing_cl}. Available: {sorted(all_vals)}"
            )
    else:
        clusters = [c for c in all_vals if c != "-1"]
        try:
            clusters = sorted(clusters, key=lambda x: int(x))
        except ValueError:
            clusters = sorted(clusters)

    if celltype_filter:
        before = list(clusters)
        clusters = filter_clusters_by_celltype(adata, groupby, clusters, celltype_col, celltype_filter)
        print(
            f"--celltype-filter {celltype_filter!r} (dominant '{celltype_col}' "
            f"contains it, case-insensitive): kept {len(clusters)}/{len(before)} "
            f"cluster(s): {clusters}"
        )

    if min_cluster_size:
        before = list(clusters)
        counts = adata.obs[groupby].astype(str).value_counts()
        dropped = [c for c in before if counts.get(c, 0) < min_cluster_size]
        clusters = [c for c in before if c not in dropped]
        if dropped:
            print(
                f"--min-cluster-size {min_cluster_size}: dropped "
                f"{len(dropped)}/{len(before)} cluster(s) with fewer cell "
                "pairs: " + ", ".join(f"{c} (n={counts.get(c, 0)})" for c in dropped)
            )

    if not clusters_arg:
        if top_n_clusters is not None:
            ranked = rank_clusters(adata, clusters, groupby, rank_clusters_by)
            clusters = ranked[:top_n_clusters]
            try:
                clusters = sorted(clusters, key=lambda x: int(x))
            except ValueError:
                clusters = sorted(clusters)
            print(f"Narrowed to top {len(clusters)} cluster(s) by {rank_clusters_by}: {clusters}")

        if include_unassigned and "-1" in all_vals:
            clusters = clusters + ["-1"]

    if not clusters:
        sys.exit(
            "No clusters to display after selection/filtering. If you "
            "passed --celltype-filter, double check it against the "
            f"dominant '{celltype_col}' values actually present (printed "
            "above, or inspect adata.obs[celltype_col].unique()); if you "
            "passed --min-cluster-size, try lowering it (printed above are "
            "the sizes of whatever got dropped); otherwise this usually "
            "means no neighborhoods were significant at the --alpha used "
            "in Step 3 -- rerun with a larger --alpha."
        )
    return clusters


def filter_clusters_by_celltype(adata, groupby, clusters, celltype_col, substring):
    """Keep only clusters whose DOMINANT `celltype_col` value (same "most
    common value" computed by dominant_celltype_labels()) contains
    `substring`, case-insensitive. Backs --celltype-filter, so you can ask
    for e.g. all macrophage-involving clusters by name instead of having to
    look up and type out cluster numbers by hand.
    """
    if celltype_col not in adata.obs.columns:
        sys.exit(
            f"--celltype-filter needs --celltype-col '{celltype_col}' in "
            f"adata.obs. Available: {list(adata.obs.columns)}"
        )
    needle = substring.lower()
    keep = []
    for cl in clusters:
        sub = adata.obs.loc[adata.obs[groupby].astype(str) == cl, celltype_col]
        if len(sub) == 0:
            continue
        dominant = sub.astype(str).value_counts().idxmax()
        if needle in dominant.lower():
            keep.append(cl)
    return keep


def dominant_celltype_labels(adata, groupby, clusters, celltype_col):
    """Map each cluster label to '<cluster>: <most common celltype_col value>'."""
    if celltype_col not in adata.obs.columns:
        print(
            f"Note: --celltype-col '{celltype_col}' not found in adata.obs -- "
            "showing raw cluster labels instead. Available columns: "
            f"{list(adata.obs.columns)}"
        )
        return {c: c for c in clusters}
    labels = {}
    for cl in clusters:
        sub = adata.obs.loc[adata.obs[groupby].astype(str) == cl, celltype_col]
        if len(sub) == 0:
            labels[cl] = cl
            continue
        dominant = sub.astype(str).value_counts().idxmax()
        labels[cl] = f"{cl}: {dominant}"
    return labels


def print_cluster_significance(adata, clusters, groupby):
    """Print per-displayed-cluster size + Milo DA stats (mean logFC,
    mean/min SpatialFDR among the cluster's own significant neighborhoods)
    -- this is the actual reference-vs-treatment significance of the
    *neighborhoods*, distinct from any per-L-R-pair test plotted on top of
    it (a cluster-vs-rest Wilcoxon in plot_lr_heatmap.py, or this cluster's
    own reference-vs-treatment Wilcoxon in plot_lr_comparison_heatmap.py).
    """
    if "nhood_adata" not in adata.uns or "louvain" not in adata.uns["nhood_adata"].obs.columns:
        print(
            "Note: can't print per-cluster Milo significance -- "
            "adata.uns['nhood_adata'] with a 'louvain' column is needed "
            "(must be the direct h5ad output of run_milo_da.py)."
        )
        return
    nhood_obs = adata.uns["nhood_adata"].obs
    finite = nhood_obs[nhood_obs["louvain"] < float("inf")].copy()
    finite["louvain_str"] = finite["louvain"].astype(int).astype(str)

    print("\nPer-cluster Milo differential-abundance stats (reference- vs "
          "treatment-condition significance of the underlying neighborhoods):")
    header = f"  {'cluster':>10} | {'n_cellpairs':>11} | {'n_nhoods':>8} | {'mean_logFC':>10} | {'mean_SpatialFDR':>15} | {'max_SpatialFDR':>14}"
    print(header)
    for cl in clusters:
        n_cellpairs = int((adata.obs[groupby].astype(str) == cl).sum())
        sub = finite[finite["louvain_str"] == cl]
        if len(sub) == 0:
            print(f"  {cl:>10} | {n_cellpairs:>11} | {'0':>8} | {'--':>10} | {'--':>15} | {'--':>14}")
            continue
        print(
            f"  {cl:>10} | {n_cellpairs:>11} | {len(sub):>8} | "
            f"{sub['logFC'].mean():>10.3f} | {sub['SpatialFDR'].mean():>15.4g} | "
            f"{sub['SpatialFDR'].max():>14.4g}"
        )
    print(
        "(mean_logFC > 0 means that cluster's neighborhoods are enriched in "
        "condition 1 relative to condition 0, as encoded by "
        "prep_niches_input.py's --treatment-value/--reference-value; "
        "max_SpatialFDR at or below your Step 3 --alpha confirms every "
        "neighborhood in that cluster individually cleared the "
        "significance bar, not just the cluster's average.)\n"
    )


def auto_select_lr_pairs(adata, clusters, wilcox_key, top_n):
    """Auto-select L-R pairs from the cluster-identity (cluster-vs-rest)
    Wilcoxon results run_milo_da.py stores under adata.uns[wilcox_key] when
    run with --wilcox-out -- the top `top_n` pairs per cluster, by score.
    """
    if wilcox_key not in adata.uns:
        sys.exit(
            f"--lr-pairs not given and adata.uns[{wilcox_key!r}] not found "
            "for auto-selection. Rerun Step 3 (run_milo_da.py) with "
            "--wilcox-out to store Wilcoxon results in the h5ad, or pass "
            "--lr-pairs explicitly."
        )
    import scanpy as sc
    pairs = []
    for cl in clusters:
        df = sc.get.rank_genes_groups_df(adata, group=str(cl), key=wilcox_key)
        top = df.sort_values("scores", ascending=False).head(top_n)["names"].tolist()
        for lr in top:
            if lr not in pairs:
                pairs.append(lr)
    return pairs


def benjamini_hochberg(pvals):
    """Benjamini-Hochberg FDR correction, plain numpy (no statsmodels
    dependency, which isn't otherwise needed by this pipeline -- see
    pyproject.toml). NaNs pass through untouched and are excluded from the
    correction itself (not counted toward `n`).
    """
    pvals = np.asarray(pvals, dtype=float)
    out = np.full(pvals.shape, np.nan)
    finite = np.isfinite(pvals)
    p = pvals[finite]
    n = len(p)
    if n == 0:
        return out
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # Enforce monotonicity, the standard BH step-up correction.
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    ranked = np.clip(ranked, 0, 1)
    adj = np.empty(n)
    adj[order] = ranked
    out[finite] = adj
    return out


def compute_comparison_stats(adata, clusters, groupby, condition_col, lr_pairs,
                              min_cells_per_group=3):
    """For every (L-R pair, cluster) combination, compare NICHES scores
    between the two `condition_col` groups (0 = reference, 1 = treatment, as
    encoded by prep_niches_input.py) among only the cell pairs in that
    cluster. Used by plot_lr_comparison_heatmap.py.

    Returns two pandas DataFrames, both indexed by `lr_pairs` with
    `clusters` as columns:

    - effect: Cohen's d = (mean_treatment - mean_reference) / pooled_std.
      Deliberately NOT scanpy's built-in "logfoldchanges" -- that assumes
      the input is already log1p-transformed (it un-logs with expm1 before
      taking a ratio), which NICHES L-R scores aren't guaranteed to be
      (they're continuous, non-negative products of ligand/receptor
      expression). Using that assumption blindly would give a number that
      looks like a fold-change but isn't one. Cohen's d needs no such
      assumption and is scale-free, so L-R pairs with very different score
      magnitudes are still comparable on one shared color scale.

    - padj: Benjamini-Hochberg-adjusted Mann-Whitney U (= Wilcoxon
      rank-sum, the same test family run_milo_da.py's --wilcox-out already
      uses elsewhere in this pipeline) p-value, corrected across every
      (L-R pair, cluster) combination tested together. Rank-based, so it
      needs no log-space assumption either.

    Cells are NaN in both outputs wherever a cluster has fewer than
    `min_cells_per_group` cell pairs in one of the two condition groups
    (too few to test) -- a printed note lists which combinations were
    skipped and why.
    """
    from scipy.stats import mannwhitneyu

    effect = pd.DataFrame(index=lr_pairs, columns=clusters, dtype=float)
    pval = pd.DataFrame(index=lr_pairs, columns=clusters, dtype=float)
    skipped = []

    for cl in clusters:
        mask_cl = (adata.obs[groupby].astype(str) == cl).values
        sub = adata[mask_cl]
        cond_sub = sub.obs[condition_col]
        mask1 = (cond_sub == 1).values
        mask0 = (cond_sub == 0).values
        if mask1.sum() < min_cells_per_group or mask0.sum() < min_cells_per_group:
            skipped.append((cl, int(mask0.sum()), int(mask1.sum())))
            continue

        X = sub[:, lr_pairs].X
        if scipy.sparse.issparse(X):
            X = X.toarray()
        X = np.asarray(X)

        for j, lr in enumerate(lr_pairs):
            v1 = X[mask1, j]
            v0 = X[mask0, j]
            n1, n0 = len(v1), len(v0)
            mean1, mean0 = v1.mean(), v0.mean()
            var1 = v1.var(ddof=1) if n1 > 1 else 0.0
            var0 = v0.var(ddof=1) if n0 > 1 else 0.0
            pooled_std = (
                np.sqrt(((n1 - 1) * var1 + (n0 - 1) * var0) / (n1 + n0 - 2))
                if (n1 + n0) > 2 else np.nan
            )
            d = (mean1 - mean0) / pooled_std if pooled_std and pooled_std > 0 else np.nan
            try:
                _, p = mannwhitneyu(v1, v0, alternative="two-sided")
            except ValueError:
                p = np.nan
            effect.loc[lr, cl] = d
            pval.loc[lr, cl] = p

    if skipped:
        print(
            "Note: skipped the reference-vs-treatment comparison for these clusters "
            f"(< {min_cells_per_group} cell pairs in one condition group -- "
            "left blank in the heatmap): " +
            ", ".join(f"{cl} (n0={n0}, n1={n1})" for cl, n0, n1 in skipped)
        )

    padj = pd.DataFrame(
        benjamini_hochberg(pval.values.ravel()).reshape(pval.shape),
        index=pval.index, columns=pval.columns,
    )
    return effect, padj


def compute_group_means(adata, clusters, groupby, condition_col, lr_pairs):
    """Mean NICHES score per (L-R pair, cluster), split by the two
    `condition_col` groups (0 = reference, 1 = treatment). Used by
    plot_lr_comparison_sidebyside.py to show raw reference/treatment magnitudes next to
    each other -- unlike compute_comparison_stats(), which only reports the
    *difference* between them (Cohen's d), this keeps each condition's own
    mean score so the two can be plotted as two directly comparable panels.

    Returns two pandas DataFrames (mean0, mean1), both indexed by `lr_pairs`
    with `clusters` as columns. A cell is NaN if that cluster has zero cell
    pairs in that condition group (nothing to average).
    """
    mean0 = pd.DataFrame(index=lr_pairs, columns=clusters, dtype=float)
    mean1 = pd.DataFrame(index=lr_pairs, columns=clusters, dtype=float)

    for cl in clusters:
        mask_cl = (adata.obs[groupby].astype(str) == cl).values
        sub = adata[mask_cl]
        cond_sub = sub.obs[condition_col]
        mask1 = (cond_sub == 1).values
        mask0 = (cond_sub == 0).values

        X = sub[:, lr_pairs].X
        if scipy.sparse.issparse(X):
            X = X.toarray()
        X = np.asarray(X)

        if mask0.sum() > 0:
            mean0[cl] = X[mask0, :].mean(axis=0)
        if mask1.sum() > 0:
            mean1[cl] = X[mask1, :].mean(axis=0)

    return mean0, mean1


def compute_log2fc(mean0, mean1, cap=4.0, zero_thresh=1e-3):
    """log2(treatment/reference) fold change from two mean-score DataFrames
    (same shape/index/columns as compute_group_means()'s output), for
    plot_lr_log2fc_heatmap.py -- ported from Bridges et al.'s
    figures/human_cd40ag_validation.py `_log2fc()` (its Fig. 7G L-R
    heatmap, a pre- vs post-sotigalimab comparison in the original),
    generalized from a per-cell function to a vectorized one over an
    arbitrary matrix, with one correctness fix (see below).

    A near-zero baseline makes a literal log2 ratio meaningless (a tiny
    denominator blows the ratio up arbitrarily), so both near-zero cases are
    capped rather than computed directly:
      - reference ~0, treatment clearly positive  -> +cap (looks like
        signaling "turned on" in the treatment condition, capped rather
        than +inf)
      - treatment ~0, reference clearly positive  -> -cap (mirror case:
        "turned off")
      - both ~0                        -> 0.0 (no signal in either
        condition -- not "unchanged", just nothing to compare)
    Anything else is a normal log2(treatment/reference).

    THE FIX vs. the original `_log2fc()`: it only ever handled the first
    case explicitly (reference~0, treatment>threshold -> +4.0); the mirror
    case (treatment~0, reference>threshold) fell through to a literal
    `np.log2(0/reference)`, which is -inf, not a capped value -- a real bug
    in the original, just one that likely never triggered in its
    15-pair x 3-interaction table
    (didn't happen to hit a clean zero). Handled symmetrically here.

    NaN cells in `mean0`/`mean1` (from compute_group_means()'s own "zero
    cell pairs" case) stay NaN -- there's nothing to compute a ratio from
    at all, which is different from "both means are ~0" (a real, if small,
    signal was measured in both conditions).
    """
    reference = mean0.values.astype(float)
    treatment = mean1.values.astype(float)

    ref_zero = reference < zero_thresh
    trt_zero = treatment < zero_thresh
    both_zero = ref_zero & trt_zero

    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = np.log2(treatment / reference)

    out = np.where(both_zero, 0.0, ratio)
    out = np.where(ref_zero & ~trt_zero, cap, out)
    out = np.where(trt_zero & ~ref_zero, -cap, out)
    out = np.where(np.isnan(reference) | np.isnan(treatment), np.nan, out)

    return pd.DataFrame(out, index=mean0.index, columns=mean0.columns)
