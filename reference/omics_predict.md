# Predict cell-type proportions with one OmicsTweezer model

Uses a trained OmicsTweezer encoder and predictor to estimate cell-type
proportions for real bulk samples.

## Usage

``` r
omics_predict(
  model,
  test_x,
  celltypes,
  samplename = NULL,
  batch_size = 128L,
  device = c("auto", "cpu", "cuda"),
  cuda_index = NULL,
  verbose = TRUE
)
```

## Arguments

- model:

  A trained model returned by \[omics_train()\]. Must contain
  \`encoder\` and \`predictor\`.

- test_x:

  Numeric matrix with samples in rows and genes in columns.

- celltypes:

  Character vector of output cell-type names.

- samplename:

  Optional character vector of sample names for the returned prediction
  table.

- batch_size:

  Integer scalar. Prediction batch size.

- device:

  Character scalar. One of \`"auto"\`, \`"cpu"\`, or \`"cuda"\`. If
  \`"auto"\`, CUDA is used when available, otherwise CPU is used.

- cuda_index:

  Optional integer scalar giving the CUDA device index to use when
  \`device = "cuda"\`. CUDA indices are zero-based.

- verbose:

  Logical scalar. If \`TRUE\`, prints the selected prediction device.

## Value

A data frame with samples in rows and predicted cell-type proportions in
columns.

## Details

Predictions are generated in evaluation mode and returned on the CPU as
a regular R \`data.frame\`.
