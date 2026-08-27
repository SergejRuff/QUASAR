# Normalize reference and target bulk matrices

Selects the most variable genes shared between a reference (pseudobulk)
matrix and a target matrix, then applies sample-wise normalization.

## Usage

``` r
normalize_bulks(
  ref_matrix,
  target_matrix,
  top_n = 0,
  top_perc = 0.98,
  normalization = c("samplewise_ss", "samplewise_minmax", "samplewise_cpm", "none")
)
```

## Arguments

- ref_matrix:

  Numeric matrix containing the reference/pseudobulk data. Rows
  correspond to genes and columns correspond to samples.

- target_matrix:

  Numeric matrix containing the target bulk data. Rows correspond to
  genes and columns correspond to samples.

- top_n:

  Integer. Number of most variable genes to retain. If greater than
  zero, this takes precedence over \`top_perc\`.

- top_perc:

  Numeric between 0 and 1. Fraction of the most variable genes to
  retain. Default is 0.98.

- normalization:

  Character specifying the normalization method:

  "samplewise_ss"

  : Log1p transformation followed by sample-wise standardization (mean
    and population standard deviation).

  "samplewise_minmax"

  : Sample-wise min-max scaling.

  "samplewise_cpm"

  : Counts per million normalization per sample.

  "none"

  : No normalization.

## Value

A list containing:

- pseudobulk_norm:

  Normalized filtered reference matrix.

- target_norm:

  Normalized filtered target matrix.

- pseudobulk_raw:

  Filtered reference matrix before normalization.

- target_raw:

  Filtered target matrix before normalization.

- genes:

  Genes retained after filtering and intersection.

- normalization_method:

  Selected normalization method.

## Details

The function first filters genes based on their variance within each
matrix. Only genes present in both filtered matrices are retained. The
filtered matrices can then be normalized using one of several
sample-wise approaches.

Genes are ranked by variance and only the most variable genes are
retained. Filtering is performed independently for the reference and
target matrices, followed by intersection of retained genes to ensure
both matrices contain identical features.

## Examples

``` r
# Create example reference and target matrices
set.seed(1)

ref <- matrix(
  rpois(300, lambda = 10),
  nrow = 100,
  ncol = 3,
  dimnames = list(
    paste0("gene", 1:100),
    paste0("ref_sample", 1:3)
  )
)

target <- matrix(
  rpois(300, lambda = 15),
  nrow = 100,
  ncol = 3,
  dimnames = list(
    paste0("gene", 1:100),
    paste0("target_sample", 1:3)
  )
)

result <- normalize_bulks(
  ref_matrix = ref,
  target_matrix = target,
  top_perc = 0.5,
  normalization = "samplewise_ss"
)

# Normalized matrices
dim(result$pseudobulk_norm)
#> [1] 25  3
dim(result$target_norm)
#> [1] 25  3

# Retained genes
head(result$genes)
#> [1] "gene36" "gene96" "gene57" "gene87" "gene32" "gene99"
```
