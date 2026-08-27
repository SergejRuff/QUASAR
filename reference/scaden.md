# Train a Scaden ensemble in torch

Trains a three-model multilayer perceptron ensemble corresponding to the
Scaden architectures \`m256\`, \`m512\`, and \`m1024\`. The function
expects a processed training matrix with samples in rows and genes in
columns, and a matching target matrix with cell-type proportions.

## Usage

``` r
scaden(
  train_x,
  train_y,
  lr = 1e-04,
  batch_size = 128,
  epochs = 20,
  seed = 123,
  device = NULL
)
```

## Arguments

- train_x:

  Numeric matrix with samples in rows and genes in columns.

- train_y:

  Numeric matrix with samples in rows and cell types in columns.

- lr:

  Numeric scalar. Learning rate passed to Adam.

- batch_size:

  Integer scalar. Minibatch size.

- epochs:

  Integer scalar. Number of training epochs for each ensemble model.

- seed:

  Integer scalar random seed.

- device:

  Optional torch device. If \`NULL\`, CUDA is used when available,
  otherwise CPU is used.

## Value

A named list with components:

- model256:

  Trained \`m256\` torch model.

- model512:

  Trained \`m512\` torch model.

- model1024:

  Trained \`m1024\` torch model.

- loss256:

  Numeric vector of batch-level losses for the \`m256\` model.

- loss512:

  Numeric vector of batch-level losses for the \`m512\` model.

- loss1024:

  Numeric vector of batch-level losses for the \`m1024\` model.

- epoch_loss256:

  Numeric vector of mean epoch losses for the \`m256\` model.

- epoch_loss512:

  Numeric vector of mean epoch losses for the \`m512\` model.

- epoch_loss1024:

  Numeric vector of mean epoch losses for the \`m1024\` model.

- total_training_time_sec:

  Numeric scalar giving total ensemble training time in seconds.

- architectures:

  List of model architectures used.

- lr:

  Learning rate used for training.

- batch_size:

  Batch size used for training.

- epochs:

  Number of training epochs.

- inputdim:

  Number of input genes.

- outputdim:

  Number of output cell types.

- celltypes:

  Character vector of output cell-type names.

- train_samples:

  Character vector of training sample names.

- device:

  Torch device used for training.

## Details

The ensemble consists of three feed-forward neural networks with softmax
output. Predictions are averaged across the three trained models in
\[scaden_predict()\].

## Examples

``` r
if (FALSE) { # \dontrun{
fit <- scaden(
  train_x = proc$train_x,
  train_y = proc$train_y,
  lr = 1e-4,
  batch_size = 128,
  epochs = 20
)
} # }
```
