# Simulate DISSECT training mixtures from single-cell data

Generates simulated bulk or spatial transcriptomics mixtures from a
single-cell reference, following the simulation logic used in the
original DISSECT workflow. Returned proportions are stored as samples x
cell types.

## Usage

``` r
dissect_simulate(
  sc_data,
  celltype_col = "celltype",
  batch_col = NULL,
  type = "bulk",
  n_samples = NULL,
  cells_per_sample = 500,
  prop_sparse = 0.5,
  concentration = NULL,
  save_expr = TRUE,
  min_genes = 200,
  min_cells = 3,
  mt_cutoff = 5,
  min_expr = 0,
  downsample = NULL,
  seed = 42
)
```

## Arguments

- sc_data:

  A \`Seurat\` or \`SingleCellExperiment\` object containing single-cell
  counts.

- celltype_col:

  Character scalar. Column name in cell metadata containing cell-type
  labels.

- batch_col:

  Character scalar or \`NULL\`. Optional column name in cell metadata
  containing batch labels. If multiple batches are present, simulation
  is performed per batch and concatenated.

- type:

  Character scalar. Either \`"bulk"\` or \`"st"\`.

- n_samples:

  Integer scalar or \`NULL\`. Number of simulated samples. If \`NULL\`,
  defaults to \`1000 \* n_celltypes\`.

- cells_per_sample:

  Integer scalar. Number of cells per simulated sample for bulk
  simulation.

- prop_sparse:

  Numeric scalar in \`\[0, 1\]\`. Fraction of sparse simulated bulk
  samples.

- concentration:

  Numeric vector or \`NULL\`. Dirichlet concentration parameter for bulk
  simulation. If \`NULL\`, a vector of ones is used.

- save_expr:

  Logical scalar. If \`TRUE\`, stores per-cell-type expression
  contributions for each simulated sample.

- min_genes:

  Integer scalar. Minimum number of detected genes required for a cell
  to be kept during preprocessing.

- min_cells:

  Integer scalar. Minimum number of cells required for a gene to be kept
  during preprocessing.

- mt_cutoff:

  Numeric scalar. Maximum allowed mitochondrial percentage per cell.

- min_expr:

  Numeric scalar. Minimum mean \`log1p\` expression threshold

- downsample:

  Numeric scalar or \`NULL\`. Optional downsampling factor for spatial
  transcriptomics simulation.

- seed:

  Integer scalar. Random seed used for simulation.

## Value

A named list with components:

- X:

  A simulated expression matrix with samples in rows and genes in
  columns.

- props:

  A proportions matrix with samples in rows and cell types in columns.

- cells:

  An integer matrix of sampled cell counts with samples in rows and cell
  types in columns.

- gene_names:

  Character vector of gene names.

- celltype_names:

  Character vector of cell-type names.

- layers:

  Named list of per-cell-type expression matrices, or \`NULL\` if
  \`save_expr = FALSE\`.

- batch:

  Character vector giving the simulated batch assignment of each sample.

- type:

  Simulation type, either \`"bulk"\` or \`"st"\`.

## Details

The preprocessing and simulation logic are designed to follow the
original DISSECT implementation as closely as possible while using
in-memory R data structures instead of file-based Python objects.

## References

Khatri, R., Machart, P., & Bonn, S. (2024). DISSECT: deep
semi-supervised consistency regularization for accurate cell type
fraction and gene expression estimation. *Genome Biology*, 25(1), 112.

Original DISSECT software repository:
<https://github.com/imsb-uke/DISSECT>

## Examples

``` r
if (FALSE) { # \dontrun{
sim <- dissect_simulate(
  sc_data = sce,
  celltype_col = "celltype",
  batch_col = "batch",
  type = "bulk"
)
} # }
```
