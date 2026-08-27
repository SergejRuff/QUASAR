# Estimate cell-type proportions with DISSECT

Trains the DISSECT fraction model and returns estimated cell-type
proportions and pre-activation scores for each sample.

## Usage

``` r
dissect_prop(
  bulk = NULL,
  reference = NULL,
  sim_data = NULL,
  processed = NULL,
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
  device = c("auto", "cpu", "cuda"),
  cuda_index = NULL
)
```

## Arguments

- bulk:

  Numeric matrix with genes in rows and samples in columns, or \`NULL\`
  if \`processed\` is supplied.

- reference:

  Numeric matrix with genes in rows and cell types in columns, or
  \`NULL\`.

- sim_data:

  A list returned by \[dissect_simulate()\], or \`NULL\`.

- processed:

  A list returned by \[dissect_process()\], or \`NULL\`.

- test_dataset_type:

  Character scalar. Either \`"bulk"\` or a supported alternative
  matching the original workflow.

- duplicated:

  Character scalar controlling duplicated gene handling.

- normalize_simulated:

  Character scalar or \`NULL\`. Currently \`"cpm"\` or \`NULL\`.

- normalize_test:

  Character scalar or \`NULL\`. Currently \`"cpm"\` or \`NULL\`.

- var_cutoff:

  Numeric scalar or \`NULL\`. Variance threshold for bulk preprocessing.

- test_in_mix:

  Integer scalar. Number of real samples to use in online mixing.

- n_hidden_layers:

  Integer scalar. Number of hidden layers in the fraction model.

- hidden_units:

  Integer vector. Number of units in each hidden layer.

- hidden_activation:

  Character scalar. Hidden-layer activation, applied to every ensemble
  model. Note that the reference Python mutates its config in place
  inside the model loop, so only the first model receives \`relu6\` and
  the remaining four fall back to \`relu\`; this implementation follows
  the documented behaviour instead.

- output_activation:

  Character scalar. Output activation function.

- loss:

  Character scalar. One of \`"kldivergence"\`, \`"l2"\`, or \`"l1"\`.

- n_steps:

  Integer scalar. Number of training steps.

- lr:

  Numeric scalar. Learning rate.

- batch_size:

  Integer scalar. Batch size.

- dropout:

  Numeric vector or \`NULL\`. Dropout rates for hidden layers.

- alpha_range:

  Numeric vector of length two. Range for the mixing coefficient used in
  online mixtures.

- normalization_per_batch:

  Character scalar or \`NULL\`. Currently \`"log1p-MinMax"\` or
  \`NULL\`.

- models:

  Integer vector. Ensemble model identifiers. As in the original Python
  code, only the number of models is used.

- mix:

  Character scalar. Mixing strategy, either \`"srm"\` or \`"rrm"\`.

- device:

  Character scalar. One of \`"auto"\`, \`"cpu"\`, or \`"cuda"\`.

- cuda_index:

  Integer scalar or \`NULL\`. Optional CUDA device index.

## Value

A named list with components:

- fractions:

  A data frame of estimated cell-type proportions with samples in rows
  and cell types in columns.

- scores:

  A data frame of pre-activation output scores with samples in rows and
  cell types in columns.

- processed:

  The processed input object used for training.

## Details

dissect_prop implements the DISSECT cell-type fraction estimation
strategy described by Khatri, Machart, and Bonn (2024) in torch in R. It
reproduces the semi-supervised consistency-regularized training workflow
for fraction estimation within an R interface.

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
prop <- dissect_prop(
  processed = proc,
  n_steps = 5000,
  device = "auto"
)
} # }
```
