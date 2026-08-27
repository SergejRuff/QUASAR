# Predict cell fractions and signature matrices with TAPE

Applies a trained TAPE model to processed bulk data, optionally with
adaptive refinement in either overall or high-resolution mode.

## Usage

``` r
tape_predict(
  model,
  test_x,
  genename,
  celltypes,
  samplename,
  adaptive = TRUE,
  mode = "overall",
  chunk_size = 16L
)
```

## Arguments

- model:

  A trained model returned by \[tape_train()\].

- test_x:

  Numeric matrix with samples in rows and genes in columns.

- genename:

  Character vector of gene names.

- celltypes:

  Character vector of cell-type names.

- samplename:

  Character vector of sample names.

- adaptive:

  Logical scalar. Whether to run adaptive refinement.

- mode:

  Character scalar. Either \`"overall"\` or \`"high-resolution"\`.

- chunk_size:

  Integer scalar. Number of samples adapted simultaneously in
  high-resolution mode.

## Value

A named list with components:

- sigm:

  A named list of per-cell-type signature matrices in high-resolution
  mode, a signature matrix data frame in overall adaptive mode, and
  \`NULL\` when \`adaptive\` is \`FALSE\`.

- pred:

  A data frame of predicted proportions with samples in rows and cell
  types in columns.

## Details

High-resolution mode adapts every sample with its own copy of the
trained model, and processes \`chunk_size\` samples concurrently as an
ensemble of independent models. The supplied model is read only in this
mode and is returned unmodified.

## Examples

``` r
if (FALSE) { # \dontrun{
res <- tape_predict(
  model = model,
  test_x = processed$test_x,
  genename = processed$genename,
  celltypes = processed$celltypes,
  samplename = processed$samplename,
  mode = "high-resolution",
  chunk_size = 16L
)
} # }
```
