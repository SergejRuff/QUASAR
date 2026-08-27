# Estimate cell-type-specific expression with DISSECT

Trains the DISSECT expression model and returns estimated expression
profiles for each sample-cell-type combination.

## Usage

``` r
dissect_expr(
  bulk = NULL,
  fractions,
  sim_data,
  normalize_simulated = "cpm",
  normalize_test = "cpm",
  n_steps_expr = 5000,
  expr_scaling = "p99",
  latent_dim = 128,
  batch_size = 128,
  lr = 0.001,
  beta_vae = 0.01,
  lambda_cons = 0.1,
  seed = 42,
  device = c("auto", "cpu", "cuda"),
  cuda_index = NULL
)
```

## Arguments

- bulk:

  Numeric matrix with genes in rows and samples in columns.

- fractions:

  A data frame or matrix of estimated cell-type fractions with samples
  in rows and cell types in columns, typically obtained from
  \[dissect_prop()\].

- sim_data:

  A list returned by \[dissect_simulate()\] with \`save_expr = TRUE\`.

- normalize_simulated:

  Character scalar or \`NULL\`. Currently \`"cpm"\` or \`NULL\`.

- normalize_test:

  Character scalar or \`NULL\`. Currently \`"cpm"\` or \`NULL\`.

- n_steps_expr:

  Integer scalar or \`NULL\`. Number of training steps for the
  expression model.

- expr_scaling:

  Character scalar. Scaling method used before model fitting, typically
  \`"p99"\` or \`"max"\`.

- latent_dim:

  Integer scalar. Latent dimension for the conditional VAE.

- batch_size:

  Integer scalar. Batch size.

- lr:

  Numeric scalar. Learning rate.

- beta_vae:

  Numeric scalar. Weight of the KL divergence term.

- lambda_cons:

  Numeric scalar. Weight of the consistency loss term.

- seed:

  Integer scalar or \`NULL\`. Optional random seed.

- device:

  Character scalar. One of \`"auto"\`, \`"cpu"\`, or \`"cuda"\`.

- cuda_index:

  Integer scalar or \`NULL\`. Optional CUDA device index.

## Value

A named list with components:

- expression_layered:

  Named list of matrices, one per cell type, with samples in rows and
  genes in columns.

- expression_combined:

  Data frame containing estimated expression for all sample-cell-type
  combinations, along with \`cell_type\` and \`sample\` columns.

- scaled_counts:

  Scaled expression matrix corresponding to the final combined
  estimates.

- celltypes:

  Character vector of cell-type names.

- genes:

  Character vector of genes used in the model.

## Details

This implementation keeps the DISSECT expression methodology the same,
but avoids materializing the full expanded \`(sample x celltype)\`
training design matrices in memory. Instead, minibatches are assembled
on the fly.

dissect_expr implements the DISSECT cell-type-specific expression
estimation workflow described by Khatri, Machart, and Bonn (2024) in
torch in R. It follows the DISSECT expression-estimation methodology
based on simulated per-cell-type layers, estimated sample fractions, and
a conditional variational autoencoder. The training logic is unchanged,
but expanded input matrices are built batch-wise rather than all at
once.

Progress bars and elapsed time are reported during training.

## References

Khatri, R., Machart, P., & Bonn, S. (2024). DISSECT: deep
semi-supervised consistency regularization for accurate cell type
fraction and gene expression estimation. *Genome Biology*, 25(1), 112.

Original DISSECT software repository:
<https://github.com/imsb-uke/DISSECT>

## Examples

``` r
if (FALSE) { # \dontrun{
expr <- dissect_expr(
  bulk = bulk_mat,
  fractions = prop$fractions,
  sim_data = sim,
  device = "auto"
)
} # }
```
