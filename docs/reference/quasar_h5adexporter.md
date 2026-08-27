# Export a Seurat, SingleCellExperiment or bulk matrix to .h5ad

Single-cell objects contribute a features x cells assay; bulk data is a
plain samples x genes matrix (`n_obs x n_vars`, the anndata convention).
An optional cell-type `fractions` matrix (samples x cell types) can be
written alongside — useful for deconvolution ground truth.

## Usage

``` r
quasar_h5adexporter(
  object,
  filename,
  assay = NULL,
  layer = "counts",
  obs = NULL,
  var = NULL,
  fractions = NULL,
  fractions_to = c("obsm", "obs", "both"),
  obsm = NULL,
  bulk_orientation = c("samples_x_genes", "genes_x_samples"),
  export_reductions = TRUE,
  overwrite = TRUE,
  verbose = TRUE
)
```

## Arguments

- object:

  a `Seurat`, `SingleCellExperiment`, or a matrix / `Matrix` /
  `data.frame` (treated as bulk).

- filename:

  output path (overwritten unless `overwrite=FALSE`).

- assay:

  assay to export. Seurat: assay name (default = active assay). SCE:
  assay name (default `"counts"` if present, else the first). Ignored
  for bulk.

- layer:

  Seurat layer/slot to pull the matrix from (default `"counts"`).

- obs, var:

  optional per-observation / per-variable metadata (`data.frame` or
  matrix). Mainly for bulk, where the object carries no metadata; row
  counts must match the matrix.

- fractions:

  optional cell-type fraction matrix, samples x cell types (rows =
  observations). Column names become cell-type labels.

- fractions_to:

  where to place the fractions: `"obsm"` (as `adata.obsm['fractions']`),
  `"obs"` (spread across obs columns), or `"both"`.

- obsm:

  optional named list of extra per-observation matrices (bulk only; for
  single-cell use `export_reductions`).

- bulk_orientation:

  orientation of a bulk matrix: `"samples_x_genes"` (default) or
  `"genes_x_samples"`.

- export_reductions:

  logical, export Seurat reductions / SCE reducedDims into `adata.obsm`.

- overwrite:

  logical, overwrite an existing file.

- verbose:

  logical, print a start line and a summary block.

## Value

(invisibly) the output `filename`.

## Examples

``` r
if (FALSE) { # \dontrun{
## --- tiny bulk export (6 samples x 20 genes) + ground-truth fractions ----
set.seed(1)
bulk <- matrix(rpois(6 * 20, 5), nrow = 6, ncol = 20,
               dimnames = list(paste0("sample_", 1:6),
                               paste0("ENSG", 1:20)))
frac <- matrix(runif(6 * 3), nrow = 6,
               dimnames = list(rownames(bulk),
                               c("Tcell", "Bcell", "Mono")))
frac <- frac / rowSums(frac)                         # rows sum to 1
quasar_h5adexporter(bulk, tempfile(fileext = ".h5ad"),
                    fractions = frac, fractions_to = "both")

## --- tiny SingleCellExperiment export (15 genes x 8 cells) ---------------
sce <- SingleCellExperiment::SingleCellExperiment(
  assays = list(counts = matrix(rpois(15 * 8, 2), nrow = 15,
                dimnames = list(paste0("gene", 1:15),
                                paste0("cell", 1:8)))))
quasar_h5adexporter(sce, tempfile(fileext = ".h5ad"))
} # }
```
