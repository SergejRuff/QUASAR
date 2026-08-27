# Import an .h5ad file as a Seurat or SingleCellExperiment object

A single entry point that reads adata.X (or adata.raw.X), adata.obs,
adata.var and adata.obsm and assembles either a `Seurat` or a
`SingleCellExperiment`. A missing `adata.raw` downgrades gracefully to
`adata.X` with a warning rather than failing.

## Usage

``` r
quasar_h5adimport(
  filename,
  as = c("seurat", "sce"),
  use.raw = TRUE,
  load.X = TRUE,
  load.obsm = TRUE,
  assay = "RNA",
  verbose = TRUE
)
```

## Arguments

- filename:

  path to the .h5ad file.

- as:

  one of `"seurat"` or `"sce"`; the object type to return.

- use.raw:

  logical, default `TRUE`. Use `adata.raw.X` / `adata.raw.var` when
  present, otherwise `adata.X` / `adata.var`.

- load.X:

  logical, whether to load the expression matrix. If `FALSE` an all-zero
  sparse matrix of the right shape is used (much faster).

- load.obsm:

  logical, whether to load `adata.obsm` as reduced dimensions /
  DimReducs.

- assay:

  Seurat assay name to populate (ignored for `as = "sce"`).

- verbose:

  logical, print a start line and a summary block.

## Value

a `Seurat` or `SingleCellExperiment` object.

## Examples

``` r
if (FALSE) { # \dontrun{
## Round-trip on tiny data: write a small bulk matrix, then read it back.
f   <- tempfile(fileext = ".h5ad")
mat <- matrix(rpois(15 * 8, 2), nrow = 15,
              dimnames = list(paste0("gene", 1:15), paste0("cell", 1:8)))
quasar_h5adexporter(mat, f, bulk_orientation = "genes_x_samples")

sce <- quasar_h5adimport(f, as = "sce",    use.raw = FALSE)
seu <- quasar_h5adimport(f, as = "seurat", use.raw = FALSE)
} # }
```
