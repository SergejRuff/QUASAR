# Train the TAPE autoencoder

Trains the TAPE autoencoder on simulated pseudobulk expression and
corresponding cell-type proportions.

## Usage

``` r
tape_train(train_x, train_y, batch_size = 128L, epochs = 128L, seed = 0L)
```

## Arguments

- train_x:

  Numeric matrix with samples in rows and genes in columns.

- train_y:

  Numeric matrix with samples in rows and cell types in columns.

- batch_size:

  Integer scalar. Batch size used for training.

- epochs:

  Integer scalar. Number of training epochs.

- seed:

  Integer scalar. Random seed used for reproducibility.

## Value

A trained model list containing \`encoder\` and \`decoder\`.

## Details

This function implements the TAPE autoencoder training procedure
described by Chen et al. (2022) in torch in R.

## Examples

``` r
if (FALSE) { # \dontrun{
model <- tape_train(
  train_x = processed$train_x,
  train_y = processed$train_y,
  epochs = 128
)
} # }
```
