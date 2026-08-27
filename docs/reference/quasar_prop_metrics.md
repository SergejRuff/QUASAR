# Calculate proportion-estimation performance metrics

Computes multiple evaluation metrics for estimated cell-type proportions
relative to ground-truth proportions after aligning both matrices by
common sample names and common cell-type names.

## Usage

``` r
quasar_prop_metrics(
  estimated_proportions,
  ground_truth_proportions,
  epsilon = 1e-12
)
```

## Arguments

- estimated_proportions:

  Numeric matrix or data.frame of estimated cell-type proportions, with
  samples in rows and cell types in columns.

- ground_truth_proportions:

  Numeric matrix or data.frame of ground-truth cell-type proportions,
  with samples in rows and cell types in columns.

- epsilon:

  Numeric scalar. Small constant used to stabilize correlation and
  divergence calculations for nearly constant vectors and zero entries.

## Value

A named list with the following components:

- cell_type_rmse:

  Named numeric vector of root mean squared error (RMSE) for each cell
  type across samples.

- cell_type_mad:

  Named numeric vector of median absolute deviation (MAD) of the
  cell-type-specific error vector across samples.

- cell_type_nmae:

  Named numeric vector of normalized mean absolute error (NMAE) for each
  cell type across samples.

- mean_ground_truth_for_nmae:

  Named numeric vector of mean ground-truth proportions for each cell
  type.

- pearson_celltype_cor:

  Named numeric vector of Pearson correlations between estimated and
  true proportions for each cell type across samples.

- spearman_celltype_cor:

  Named numeric vector of Spearman correlations between estimated and
  true proportions for each cell type across samples.

- per_celltype_jsd:

  Named numeric vector of Jensen-Shannon distances for each cell type,
  where the distribution is formed across samples.

- per_sample_jsd:

  Named numeric vector of Jensen-Shannon distances for each sample,
  where the distribution is formed across cell types.

## Details

Both inputs must be numeric matrices or coercible to numeric matrices,
with samples in rows and cell types in columns. Only the intersecting
samples and cell types are retained for metric calculation.

Let \\\hat{p}\_{ic}\\ denote the estimated proportion for sample \\i\\
and cell type \\c\\, and let \\p\_{ic}\\ denote the corresponding
ground-truth proportion. After restricting both inputs to common samples
and cell types, the following metrics are calculated.

**1. Root mean squared error (RMSE) per cell type**

For each cell type \\c\\, RMSE is calculated across all shared samples:
\$\$ \mathrm{RMSE}\_c = \sqrt{\frac{1}{n}\sum\_{i=1}^{n}(\hat{p}\_{ic} -
p\_{ic})^2} \$\$ where \\n\\ is the number of shared samples.

**2. Median absolute deviation (MAD) of the error per cell type**

For each cell type \\c\\, the error vector is first defined as \\e\_{ic}
= \hat{p}\_{ic} - p\_{ic}\\. The reported MAD is then: \$\$
\mathrm{MAD}\_c = \mathrm{median}\_i\left(\left\|e\_{ic} -
\mathrm{median}\_j(e\_{jc})\right\|\right) \$\$ This quantifies the
spread of the cell-type-specific errors across samples.

**3. Normalized mean absolute error (NMAE) per cell type**

For each cell type \\c\\, the mean absolute error is normalized by the
range of the ground-truth proportions across samples: \$\$
\mathrm{NMAE}\_c = \frac{1}{n}\sum\_{i=1}^{n} \frac{\|\hat{p}\_{ic} -
p\_{ic}\|}{\max_i(p\_{ic}) - \min_i(p\_{ic})} \$\$ when the ground-truth
range is non-zero.

If the ground-truth range for a cell type is zero, the function uses the
following fallback:

- returns \\0\\ if estimated and true values are identical for all
  samples,

- otherwise returns the unnormalized mean absolute error.

**4. Mean ground-truth proportion per cell type**

For each cell type \\c\\, the mean ground-truth proportion is: \$\$
\bar{p}\_c = \frac{1}{n}\sum\_{i=1}^{n} p\_{ic} \$\$

**5. Pearson correlation per cell type**

For each cell type \\c\\, the Pearson correlation is calculated between
the estimated and true proportions across samples: \$\$ r_c =
\mathrm{cor}(\hat{p}\_{\cdot c}, p\_{\cdot c}) \$\$

If the estimated or true vector has near-zero standard deviation, random
noise with standard deviation \`epsilon\` is added before correlation
calculation to avoid undefined results for constant vectors.

**6. Spearman correlation per cell type**

For each cell type \\c\\, the Spearman rank correlation is calculated
across samples: \$\$ \rho_c = \mathrm{cor}(\hat{p}\_{\cdot c}, p\_{\cdot
c}, \mathrm{method} = "spearman") \$\$

As for Pearson correlation, a small perturbation is added when one of
the vectors has near-zero standard deviation.

**7. Jensen-Shannon distance (JSD) per cell type**

For each cell type \\c\\, the function compares the distribution of that
cell type across samples between ground truth and estimation. Let \\P_c
= (p\_{1c}, \dots, p\_{nc})\\ and \\Q_c = (\hat{p}\_{1c}, \dots,
\hat{p}\_{nc})\\. These vectors are first normalized to sum to 1: \$\$
P_c^\* = \frac{P_c}{\sum_i p\_{ic}}, \qquad Q_c^\* = \frac{Q_c}{\sum_i
\hat{p}\_{ic}} \$\$

The midpoint distribution is: \$\$ M_c = \frac{1}{2}(P_c^\* + Q_c^\*)
\$\$

The Kullback-Leibler divergence is computed as: \$\$ \mathrm{KL}(P
\parallel Q) = \sum_k P_k \log\left(\frac{P_k}{Q_k}\right) \$\$ with
\`epsilon\` added internally to both \\P\\ and \\Q\\ for numerical
stability.

The reported Jensen-Shannon distance is: \$\$ \mathrm{JSD}\_c = \sqrt{
\frac{1}{2}\mathrm{KL}(P_c^\* \parallel M_c) +
\frac{1}{2}\mathrm{KL}(Q_c^\* \parallel M_c) } \$\$

**8. Jensen-Shannon distance (JSD) per sample**

For each sample \\i\\, the function compares the distribution of
cell-type proportions across cell types between ground truth and
estimation. Let \\P_i = (p\_{i1}, \dots, p\_{im})\\ and \\Q_i =
(\hat{p}\_{i1}, \dots, \hat{p}\_{im})\\, where \\m\\ is the number of
shared cell types. After normalization to sum to 1, the same
Jensen-Shannon distance formula is applied: \$\$ \mathrm{JSD}\_i =
\sqrt{ \frac{1}{2}\mathrm{KL}(P_i^\* \parallel M_i) +
\frac{1}{2}\mathrm{KL}(Q_i^\* \parallel M_i) } \$\$ with \\M_i =
\frac{1}{2}(P_i^\* + Q_i^\*)\\.

## Author

Sergej Ruff

## Examples

``` r
if (FALSE) { # \dontrun{
metrics <- quasar_prop_metrics(
  estimated_proportions = pred_mat,
  ground_truth_proportions = truth_mat
)

metrics$cell_type_rmse
metrics$pearson_celltype_cor
metrics$per_sample_jsd
} # }
```
