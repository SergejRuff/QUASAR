# Train one OmicsTweezer model

Trains one OmicsTweezer neural network with supervised cell-type
proportion prediction on simulated pseudobulks and a Wasserstein-style
domain adaptation term using real bulk samples.

## Usage

``` r
omics_train(
  train_x,
  train_y,
  test_x,
  dims,
  drops,
  epochs = 128L,
  batch_size = 128L,
  device = c("auto", "cpu", "cuda"),
  cuda_index = NULL,
  learning_rate = 1e-04,
  loss_weight = 1,
  seed = 2021L,
  verbose = TRUE
)
```

## Arguments

- train_x:

  Numeric matrix with samples in rows and genes in columns.

- train_y:

  Numeric matrix with samples in rows and cell types in columns.

- test_x:

  Numeric matrix with target-domain bulk samples in rows and genes in
  columns.

- dims:

  Integer vector of length 4 specifying hidden-layer dimensions.

- drops:

  Numeric vector of length 4 specifying dropout probabilities.

- epochs:

  Integer scalar. Number of training epochs.

- batch_size:

  Integer scalar. Minibatch size.

- device:

  Character scalar. One of \`"auto"\`, \`"cpu"\`, or \`"cuda"\`. If
  \`"auto"\`, CUDA is used when available, otherwise CPU is used.

- cuda_index:

  Optional integer scalar giving the CUDA device index to use when
  \`device = "cuda"\`. CUDA indices are zero-based.

- learning_rate:

  Numeric scalar. Learning rate passed to Adam.

- loss_weight:

  Numeric scalar. Weight applied to the Wasserstein-style domain
  adaptation loss.

- seed:

  Integer scalar random seed.

- verbose:

  Logical scalar. If \`TRUE\`, prints an epoch-level progress bar.

## Value

A named list with components:

- encoder:

  Trained torch encoder module.

- predictor:

  Trained torch predictor module.

- history:

  Data frame with epoch-level total, MSE, OT, and weighted OT losses.

## Details

Each batch contains source-domain simulated pseudobulks with known
proportions and target-domain real bulk profiles without labels. The
training objective combines supervised MSE loss with a domain adaptation
term computed from the encoded source and target embeddings.

## References

TAPE project PyTorch implementation:
<https://github.com/poseidonchan/TAPE/blob/main/TAPE/model.py>
