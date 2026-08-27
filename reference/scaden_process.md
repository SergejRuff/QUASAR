# Process simulated and bulk data for Scaden

Aligns simulated pseudobulk training data with real bulk data, filters
genes, log-transforms, and applies per-sample scaling. Two
gene-filtering strategies are available via \`mode\`.

## Usage

``` r
scaden_process(
  sim_data,
  bulk_data,
  mode = c("tape", "scaden"),
  variance_threshold = 0.98,
  scaler = c("mms", "ss"),
  var_cutoff = 0.1,
  top_var_genes = NULL,
  verbose = TRUE
)
```

## Arguments

- sim_data:

  A list returned by \[scaden_sim_pb()\], with \`train_x\` and
  \`train_y\`.

- bulk_data:

  Bulk expression matrix or \`data.frame\`. Gene names must appear in
  row or column names; orientation is inferred.

- mode:

  Character scalar. \`"tape"\` (default) reproduces TAPE's
  \`ProcessInputData\`; \`"scaden"\` reproduces original Scaden's
  preprocessing.

- variance_threshold:

  Numeric scalar in \`\[0, 1)\`. \`mode = "tape"\` only. Fraction of
  genes defining the variance cutoff, applied separately to each matrix.

- scaler:

  Character scalar. \`mode = "tape"\` only. \`"mms"\` for per-sample
  min-max, \`"ss"\` for per-sample standardization.

- var_cutoff:

  Numeric scalar. \`mode = "scaden"\` only. Genes in \`bulk_data\` with
  variance at or below this value are dropped. Set to \`NULL\` to skip.

- top_var_genes:

  Optional integer scalar. \`mode = "scaden"\` only. Keep only the most
  variable genes after \`var_cutoff\`.

- verbose:

  Logical scalar. Print progress.

## Value

A named list with \`train_x\`, \`train_y\`, \`test_x\`, \`genename\`,
\`celltypes\`, \`samplename\`, \`mode\`, and \`elapsed\`.

## Details

The two modes differ only in gene filtering and scaler choice:

- \`"tape"\`:

  Top-quantile variance cutoff applied to \*\*both\*\* matrices, then
  intersect. \`log(x + 1)\`, then \`mms\` or \`ss\`.

- \`"scaden"\`:

  Absolute variance cutoff applied to the \*\*bulk matrix only\*\*, then
  intersect with the training genes. \`log2(x + 1)\`, then per-sample
  min-max (\`log_min_max\`); \`scaler\` is ignored.

The log base is in fact immaterial: both scalers are invariant to
multiplication by a positive constant, so \`log2\` and \`log\` give
identical output up to floating-point rounding. The bases are kept
literal to their respective sources.

Both scalers operate per sample (\`fit_transform(x.T).T\` in sklearn),
use population statistics (\`ddof = 0\`), and map constant samples to
exactly zero, matching sklearn's \`\_handle_zeros_in_scale\`. The
\`"tape"\` variance filter uses \`ddof = 1\` to match
\`pandas.DataFrame.var\`.
