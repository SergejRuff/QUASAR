# Utilities

Functions that support the deconvolution workflows but are useful on
their own: pseudo-bulk simulation, evaluation metrics, normalisation,
and `.h5ad` interoperability.

``` r

library(quasar)
library(QuasarDeconData)

sc_ref <- load_cov_imm("train")
pbulk  <- load_cov_pbulk(1)
truth  <- pbulk$ground_truth_proportions
```

## Pseudo-bulk simulation

Each ported method ships its own simulator matching its Python original.
[`quasar_sim_bulk()`](https://sergejruff.github.io/QUASAR/reference/quasar_sim_bulk.md)
is a general-purpose alternative, independent of any particular method.

``` r

sim <- quasar_sim_bulk(
  sc_ref,
  n_bulk_samples          = 1000,
  cells_per_bulk          = 500,
  cell_type_column        = "cell_type",
  patient_id_column       = "donor_id",
  return_signature_matrix = TRUE
)

dim(sim$bulk_expression_profiles)
head(sim$ground_truth_proportions)
```

Two details are worth knowing.

Supplying `patient_id_column` draws each pseudo-bulk from a single donor
where that donor has cells of the required type, which preserves
donor-specific expression heterogeneity instead of averaging it away.

The returned `ground_truth_proportions` reflect the **realised** cell
allocations after rounding, not the raw Dirichlet draws. Rows sum to one
and match the cells actually pooled — which is what you want when
scoring, since a method cannot recover proportions that were never in
the mixture.

`return_signature_matrix = TRUE` additionally computes per-cell-type
mean profiles and a global genes × cell-types signature matrix.

### Spatial mode

``` r

spots <- quasar_sim_bulk(
  sc_ref,
  n_bulk_samples   = 500,
  mode             = "spatial",
  cell_type_column = "cell_type"
)

range(spots$cells_per_sample)
table(rowSums(spots$ground_truth_proportions > 0))
```

Each spot draws how many cell types it contains, then pools a small
number of cells drawn from `cells_per_spot_range`. Allocation uses
rounding rather than flooring, because flooring a handful of cells would
empty most spots. The `sparse` and `rare` arguments are ignored here —
sparsity is already imposed by per-spot cell-type selection.

## Evaluation metrics

``` r

metrics <- quasar_prop_metrics(
  estimated_proportions    = predictions,
  ground_truth_proportions = truth
)
```

Estimates and ground truth are aligned automatically on shared sample
and cell-type names, so column order and extra cell types in either
matrix are handled without manual subsetting.

Per cell type, across samples:

``` r

metrics$cell_type_rmse            # root mean squared error
metrics$cell_type_mad             # median absolute deviation of the error
metrics$cell_type_nmae            # MAE normalised by ground-truth range
metrics$pearson_celltype_cor
metrics$spearman_celltype_cor
metrics$per_celltype_jsd
metrics$mean_ground_truth_for_nmae
```

Per sample, across cell types:

``` r

metrics$per_sample_jsd
```

`mean_ground_truth_for_nmae` is included because per-cell-type errors
are only interpretable against the abundance of that cell type. An RMSE
of 0.01 means something very different for a population at 40% than for
one at 1%.

A note on the correlation functions: if either vector has near-zero
standard deviation, noise of magnitude `epsilon` is added before
computing, so constant vectors return a finite value rather than `NA`.

## Normalisation

[`normalize_bulks()`](https://sergejruff.github.io/QUASAR/reference/normalize_bulks.md)
selects the most variable genes shared between a reference and a target
matrix, then applies sample-wise normalisation. It is independent of the
method-specific `*_process()` functions and useful when preparing data
for something else.

``` r

norm <- normalize_bulks(
  ref_matrix    = sim$bulk_expression_profiles,
  target_matrix = pbulk$bulk_expression_profiles,
  top_perc      = 0.98,
  normalization = "samplewise_ss"
)

dim(norm$pseudobulk_norm)
head(norm$genes)
```

Filtering happens independently on each matrix, then the retained gene
sets are intersected so both end up with identical features. Available
methods are `"samplewise_ss"` (log1p then standardisation),
`"samplewise_minmax"`, `"samplewise_cpm"`, and `"none"`. Both raw
filtered matrices are returned alongside the normalised ones.

## Interoperability with Python

References and results move between R and Python via `.h5ad`.

``` r

sce <- quasar_h5adimport("reference.h5ad", as = "sce")
seu <- quasar_h5adimport("reference.h5ad", as = "seurat")
```

The importer handles modern per-column categorical and nullable
encodings as well as two legacy factor layouts, sparse CSR/CSC matrices,
and `obsm` matrices. A missing `adata.raw` downgrades to `adata.X` with
a warning rather than failing. `load.X = FALSE` reads only the metadata,
which is much faster when you just want to inspect `obs` and `var`.

``` r

quasar_h5adexporter(
  object           = pbulk$bulk_expression_profiles,
  filename         = "bulk.h5ad",
  fractions        = truth,
  fractions_to     = "both",
  bulk_orientation = "genes_x_samples"
)
```

The exporter accepts a `Seurat` object, a `SingleCellExperiment`, or a
plain bulk matrix. An optional `fractions` matrix can be written to
`adata.obsm['fractions']`, spread across `obs` columns, or both — useful
for shipping deconvolution ground truth alongside the expression data.

Two lower-level readers are exported for cases where you want one
component rather than a whole object:

``` r

obs <- quasar_h5_read_dataframe("reference.h5ad", "obs")
X   <- quasar_h5_read_matrix("reference.h5ad", "X")
```
