# Run the full DISSECT workflow

Convenience wrapper that runs simulation, data processing, fraction
estimation, and optionally expression estimation in sequence.

## Usage

``` r
dissect(
  sc_data = NULL,
  bulk,
  reference = NULL,
  celltype_col = "celltype",
  batch_col = NULL,
  type = "bulk",
  n_samples = NULL,
  cells_per_sample = 500,
  prop_sparse = 0.5,
  concentration = NULL,
  save_expr = TRUE,
  min_genes = 200,
  min_cells = 3,
  mt_cutoff = 5,
  min_expr = 0,
  downsample = NULL,
  test_dataset_type = "bulk",
  duplicated = "first",
  normalize_simulated = "cpm",
  normalize_test = "cpm",
  var_cutoff = 0.1,
  test_in_mix = 1,
  n_hidden_layers = 4,
  hidden_units = c(512, 256, 128, 64),
  hidden_activation = "relu6",
  output_activation = "softmax",
  loss = "kldivergence",
  n_steps = 5000,
  lr = 1e-05,
  batch_size = 64,
  dropout = NULL,
  alpha_range = c(0.1, 0.9),
  normalization_per_batch = "log1p-MinMax",
  models = c(1, 2, 3, 4, 5),
  mix = "srm",
  n_steps_expr = 5000,
  expr_scaling = "p99",
  latent_dim = 128,
  batch_size_expr = 128,
  lr_expr = 0.001,
  beta_vae = 0.01,
  lambda_cons = 0.1,
  seed = 42,
  device = c("auto", "cpu", "cuda"),
  cuda_index = NULL
)
```

## Arguments

- sc_data:

  A \`Seurat\` or \`SingleCellExperiment\` object, or \`NULL\`. If
  supplied, simulated mixtures are generated with
  \[dissect_simulate()\].

- bulk:

  Numeric matrix with genes in rows and samples in columns.

- reference:

  Numeric matrix with genes in rows and cell types in columns, or
  \`NULL\`.

- celltype_col:

  Character scalar. Column name in the single-cell metadata containing
  cell-type labels.

- batch_col:

  Character scalar or \`NULL\`. Optional column name in the single-cell
  metadata containing batch labels.

- type:

  Character scalar. Simulation type, \`"bulk"\` or \`"st"\`.

- n_samples:

  Integer scalar or \`NULL\`. Number of simulated samples.

- cells_per_sample:

  Integer scalar. Number of cells per bulk sample.

- prop_sparse:

  Numeric scalar. Fraction of sparse bulk samples.

- concentration:

  Numeric vector or \`NULL\`. Dirichlet concentration parameter for bulk
  simulation.

- save_expr:

  Logical scalar. If \`TRUE\`, stores simulated per-cell-type expression
  and enables downstream expression estimation.

- min_genes:

  Integer scalar. Minimum genes per cell for preprocessing.

- min_cells:

  Integer scalar. Minimum cells per gene for preprocessing.

- mt_cutoff:

  Numeric scalar. Mitochondrial percentage cutoff.

- min_expr:

  Numeric scalar. Minimum mean \`log1p\` expression threshold.

- downsample:

  Numeric scalar or \`NULL\`. Downsampling factor for ST simulation.

- test_dataset_type:

  Character scalar. Either \`"bulk"\` or \`"microarray"\`.

- duplicated:

  Character scalar. How duplicated genes should be handled.

- normalize_simulated:

  Character scalar or \`NULL\`. Currently \`"cpm"\` or \`NULL\`.

- normalize_test:

  Character scalar or \`NULL\`. Currently \`"cpm"\` or \`NULL\`.

- var_cutoff:

  Numeric scalar or \`NULL\`. Variance threshold for bulk preprocessing.

- test_in_mix:

  Integer scalar. Number of real samples used in online mixing.

- n_hidden_layers:

  Integer scalar. Number of hidden layers in the fraction model.

- hidden_units:

  Integer vector. Units per hidden layer.

- hidden_activation:

  Character scalar. Hidden-layer activation.

- output_activation:

  Character scalar. Output-layer activation.

- loss:

  Character scalar. Fraction-model loss.

- n_steps:

  Integer scalar. Fraction-model training steps.

- lr:

  Numeric scalar. Fraction-model learning rate.

- batch_size:

  Integer scalar. Fraction-model batch size.

- dropout:

  Numeric vector or \`NULL\`. Fraction-model dropout rates.

- alpha_range:

  Numeric vector of length two. Mixing coefficient range.

- normalization_per_batch:

  Character scalar or \`NULL\`. Batch normalisation mode for the
  fraction model.

- models:

  Integer vector. Ensemble model identifiers.

- mix:

  Character scalar. Mixing strategy.

- n_steps_expr:

  Integer scalar or \`NULL\`. Expression-model training steps.

- expr_scaling:

  Character scalar. Expression scaling mode.

- latent_dim:

  Integer scalar. Expression-model latent dimension.

- batch_size_expr:

  Integer scalar. Expression-model batch size.

- lr_expr:

  Numeric scalar. Expression-model learning rate.

- beta_vae:

  Numeric scalar. Weight of the KL term in the expression model.

- lambda_cons:

  Numeric scalar. Weight of the consistency term in the expression
  model.

- seed:

  Integer scalar. Random seed.

- device:

  Character scalar. One of \`"auto"\`, \`"cpu"\`, or \`"cuda"\`.

- cuda_index:

  Integer scalar or \`NULL\`. Optional CUDA device index.

## Value

A named list with components:

- fractions:

  Estimated cell-type fractions.

- scores:

  Estimated pre-activation fraction scores.

- expression:

  Expression estimation results, or \`NULL\` if not run.

- sim_data:

  Simulation results, or \`NULL\` if no single-cell object was provided.

- processed:

  Processed matrices used for training and prediction.

## Details

This is the main end-to-end wrapper for the R implementation of DISSECT.
If \`sc_data\` is provided, the workflow uses simulated training data
generated from the single-cell reference. If \`save_expr = TRUE\`, the
expression model is trained after fraction estimation.

Elapsed time for each major step and the total runtime are reported.

## References

Khatri, R., Machart, P., & Bonn, S. (2024). DISSECT: deep
semi-supervised consistency regularization for accurate cell type
fraction and gene expression estimation. *Genome Biology*, 25(1), 112.

Original DISSECT software repository:
<https://github.com/imsb-uke/DISSECT>

## Examples

``` r
if (FALSE) { # \dontrun{
res <- dissect(
  sc_data = sce,
  bulk = bulk_mat,
  celltype_col = "celltype",
  batch_col = "batch",
  device = "auto"
)
} # }
```
