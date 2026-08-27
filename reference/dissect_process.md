# Prepare DISSECT input matrices

Prepares real and simulated matrices for DISSECT fraction estimation by
applying the same preprocessing order used in the original workflow:
transformation, variance filtering, deduplication, normalisation, gene
intersection, and size balancing.

## Usage

``` r
dissect_process(
  bulk,
  reference = NULL,
  sim_data = NULL,
  test_dataset_type = "bulk",
  duplicated = "first",
  normalize_simulated = "cpm",
  normalize_test = "cpm",
  var_cutoff = 0.1,
  test_in_mix = 1
)
```

## Arguments

- bulk:

  Numeric matrix with genes in rows and samples in columns.

- reference:

  Numeric matrix with genes in rows and cell types in columns, or
  \`NULL\`. This is a convenience path and is only used when
  \`sim_data\` is not supplied.

- sim_data:

  A list returned by \[dissect_simulate()\], or \`NULL\`.

- test_dataset_type:

  Character scalar. Either \`"bulk"\` or \`"microarray"\`.

- duplicated:

  Character scalar. How duplicated gene names should be resolved. One of
  \`"first"\`, \`"sum"\`, or \`"mean"\`.

- normalize_simulated:

  Character scalar or \`NULL\`. Currently \`"cpm"\` or \`NULL\`.

- normalize_test:

  Character scalar or \`NULL\`. Currently \`"cpm"\` or \`NULL\`.

- var_cutoff:

  Numeric scalar or \`NULL\`. Variance threshold applied to the bulk
  input before transposition.

- test_in_mix:

  Integer scalar. Number of real samples used in the online mixing step.

## Value

A named list with components:

- X_real_train:

  Real training matrix with samples in rows and genes in columns.

- X_sim:

  Simulated matrix with samples in rows and genes in columns.

- y_sim:

  Simulated proportions with samples in rows and cell types in columns.

- X_real_test:

  Real test matrix with samples in rows and genes in columns.

- sample_names:

  Character vector of sample names.

- celltypes:

  Character vector of cell-type names.

- genes:

  Character vector of common genes used for modelling.

- reference:

  Reference matrix if supplied, otherwise \`NULL\`.

- sim_data:

  Simulation object if supplied, otherwise \`NULL\`.

## Details

The intended DISSECT workflow uses \`sim_data\` generated from
\[dissect_simulate()\]. A direct \`reference\` matrix is supported as a
convenience interface for proportion estimation.

## References

Khatri, R., Machart, P., & Bonn, S. (2024). DISSECT: deep
semi-supervised consistency regularization for accurate cell type
fraction and gene expression estimation. *Genome Biology*, 25(1), 112.

Original DISSECT software repository:
<https://github.com/imsb-uke/DISSECT>

## Examples

``` r
if (FALSE) { # \dontrun{
proc <- dissect_process(
  bulk = bulk_mat,
  sim_data = sim
)
} # }
```
