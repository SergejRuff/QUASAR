# Predict cell-type proportions with a trained Scaden ensemble

Uses a fitted Scaden ensemble returned by \[scaden()\] to predict
cell-type proportions for new bulk samples.

## Usage

``` r
scaden_predict(fit, test_x, device = NULL)
```

## Arguments

- fit:

  A fitted Scaden object returned by \[scaden()\].

- test_x:

  Numeric matrix with samples in rows and genes in columns.

- device:

  Optional torch device. If \`NULL\`, the device stored in \`fit\` is
  used.

## Value

A named list with components:

- average_output:

  Average prediction across the three ensemble models.

- prediction_model256:

  Predictions from the \`m256\` model.

- prediction_model512:

  Predictions from the \`m512\` model.

- prediction_model1024:

  Predictions from the \`m1024\` model.

- model256:

  The fitted \`m256\` model.

- model512:

  The fitted \`m512\` model.

- model1024:

  The fitted \`m1024\` model.

## Details

Predictions are generated independently for all three ensemble members
and then averaged sample-wise.

## Examples

``` r
if (FALSE) { # \dontrun{
pred <- scaden_predict(
  fit = fit,
  test_x = proc$test_x
)
} # }
```
