# Process simulated and bulk data for OmicsTweezer

Aligns simulated pseudobulk training data with real bulk expression
data, filters genes by variance, applies log transformation and
per-sample scaling, and returns matrices ready for OmicsTweezer training
and prediction.

## Usage

``` r
omics_process(
  simudata,
  real_bulk,
  variance_threshold = 0.98,
  scaler = c("ss", "mms", "none")
)
```

## Arguments

- simudata:

  A list returned by \[omics_simulate()\]. Must contain \`X\` and
  \`obs\`.

- real_bulk:

  Numeric matrix or \`data.frame\` with genes in rows and bulk samples
  in columns.

- variance_threshold:

  Numeric scalar. Fraction of genes used to define the variance cutoff
  separately in simulated and real bulk data.

- scaler:

  Character scalar specifying the per-sample scaling method. One of
  \`"ss"\` for standard scaling, \`"mms"\` for min-max scaling, or
  \`"none"\`.

## Value

A named list with components:

- train_x:

  Processed training matrix with samples in rows and genes in columns.

- train_y:

  Training target matrix with samples in rows and cell types in columns.

- test_x:

  Processed real bulk matrix with samples in rows and genes in columns.

- genename:

  Character vector of retained shared genes.

- celltypes:

  Character vector of cell-type names.

- samplename:

  Character vector of real bulk sample names.

## Details

Both simulated and real bulk matrices are filtered by variance before
being restricted to their shared gene set. The log1p transformation is
applied before optional per-sample scaling.
