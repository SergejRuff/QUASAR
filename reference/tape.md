# Run the full TAPE workflow

Convenience wrapper that simulates pseudobulks from single-cell data,
processes simulated and real bulk data, trains the TAPE autoencoder, and
predicts cell-type fractions with optional adaptive refinement.

## Usage

``` r
tape(
  sc_data,
  real_bulk,
  variance_threshold = 0.98,
  scaler = "mms",
  d_prior = NULL,
  mode = "overall",
  adaptive = TRUE,
  sparse = TRUE,
  batch_size = 128L,
  epochs = 128L,
  seed = 0L,
  samplenum = 5000L,
  n = 500L,
  celltype_col = "CellType",
  assay = "RNA",
  slot = "counts"
)
```

## Arguments

- sc_data:

  A \`Seurat\` object, a \`SingleCellExperiment\` object, or a
  matrix/data.frame with cells in rows and genes in columns.

- real_bulk:

  Numeric matrix with genes in rows and samples in columns.

- variance_threshold:

  Numeric scalar in \`\[0, 1\]\` used for variance-based gene filtering.

- scaler:

  Character scalar. Either \`"mms"\` or \`"ss"\`.

- d_prior:

  Numeric vector or \`NULL\`. Dirichlet prior used in simulation.

- mode:

  Character scalar. Either \`"overall"\` or \`"high-resolution"\`.

- adaptive:

  Logical scalar. Whether to use adaptive refinement.

- sparse:

  Logical scalar. Whether to simulate sparse mixtures.

- batch_size:

  Integer scalar. Batch size for training.

- epochs:

  Integer scalar. Number of training epochs.

- seed:

  Integer scalar. Random seed.

- samplenum:

  Integer scalar. Number of simulated pseudobulk samples.

- n:

  Integer scalar. Number of cells per simulated pseudobulk.

- celltype_col:

  Character scalar. Metadata column containing cell-type labels.

- assay:

  Character scalar giving the assay name for \`Seurat\` input.

- slot:

  Character scalar giving the assay slot or assay name to extract.

## Value

A named list with components:

- sigm:

  Predicted signature matrix output, depending on \`mode\` and
  \`adaptive\`.

- pred:

  Predicted cell-type proportions with samples in rows and cell types in
  columns.

## Details

This function provides an R implementation of TAPE (Tissue-AdaPtive
autoEncoder) as described by Chen et al. (2022), using torch for model
training, adaptive refinement, and prediction. The implementation
follows the TAPE methodology for pseudobulk simulation, preprocessing,
autoencoder training, and optional tissue-adaptive prediction within an
R-based workflow.

The model architecture and workflow are implemented in torch in R
following the PyTorch implementation provided in the original TAPE
repository.

## References

Chen, Y., Wang, Y., Chen, Y., Cheng, Y., Wei, Y., Li, Y., Wang, J., Wei,
Y., Chan, T.-F., & Li, Y. (2022). Deep autoencoder for interpretable
tissue-adaptive deconvolution and cell-type-specific gene analysis.
*Nature Communications*, 13(1), 6735.

Original TAPE software repository:
<https://github.com/poseidonchan/TAPE>

## Examples

``` r
if (FALSE) { # \dontrun{
res <- tape(
  sc_data = sce,
  real_bulk = bulk_mat,
  celltype_col = "CellType"
)
} # }
```
