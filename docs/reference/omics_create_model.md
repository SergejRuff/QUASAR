# Create an OmicsTweezer model

Builds the encoder and predictor modules for one OmicsTweezer
architecture. The encoder contains two fully connected blocks, and the
predictor contains two fully connected blocks followed by a linear
output layer and softmax.

## Usage

``` r
omics_create_model(feature_num, celltype_num, dims, drops)
```

## Arguments

- feature_num:

  Integer scalar. Number of input genes.

- celltype_num:

  Integer scalar. Number of output cell types.

- dims:

  Integer vector of length 4 specifying hidden-layer dimensions.

- drops:

  Numeric vector of length 4 specifying dropout probabilities.

## Value

A named list with components:

- encoder:

  Torch encoder module.

- predictor:

  Torch predictor module returning cell-type proportions.

## Details

In R torch, softmax over the feature axis of an \`N x F\` tensor is
specified with \`dim = 2\`.
