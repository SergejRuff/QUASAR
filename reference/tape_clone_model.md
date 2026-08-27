# Clone a trained TAPE model

Internal helper that creates an independent copy of a trained TAPE model
by copying encoder and decoder state dictionaries in memory.

## Usage

``` r
tape_clone_model(model)
```

## Arguments

- model:

  A trained TAPE model.

## Value

A new model list with the same weights as \`model\`.
