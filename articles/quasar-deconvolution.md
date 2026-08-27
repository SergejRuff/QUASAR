# Get started

## Installation

``` r

# install.packages("pak")
pak::pak("SergejRuff/QUASAR")
```

The example datasets used throughout the documentation live in a
companion data package:

``` r

pak::pak("SergejRuff/QuasarDeconData")
```

All deep-learning methods in this package are built on `torch`. If this
is a fresh installation, run the following once:

``` r

torch::install_torch()
```

Every training and prediction function accepts `device = "auto"`,
`"cpu"`, or `"cuda"`. `"auto"` selects CUDA when a compatible `torch`
build is available and falls back to CPU otherwise, so the same code
runs on a laptop and on a GPU node without modification.

``` r

library(quasar)
library(QuasarDeconData)
```

## What this package does

Bulk transcriptomics measures gene expression from mixtures of cell
types, so a difference between two samples can reflect a change in
cellular composition, a change within individual cell types, or both.
Cell-type deconvolution estimates those proportions from bulk data and
makes the two explanations separable.

The package provides three things:

- **QUASAR**, an uncertainty-aware deconvolution method that returns
  confidence and prediction intervals alongside point estimates.
- **Native R implementations** of Scaden, TAPE, DISSECT, and
  OmicsTweezer, written in `torch` for R with GPU support and direct
  `Seurat` / `SingleCellExperiment` compatibility.
- **Supporting utilities** for pseudo-bulk simulation, normalisation,
  evaluation metrics, and `.h5ad` interoperability.

## Where to go next

| If you want to | Go to |
|----|----|
| Estimate proportions with uncertainty intervals | **QUASAR** |
| Compare deconvolution methods, or use one in R instead of Python | **Deep learning methods** |
| Load the COVID-19, PBMC, or benchmark datasets | **Data** |
| Simulate pseudo-bulks, score predictions, read or write `.h5ad` | **Utilities** |
| Look up a specific function | **Reference** |

If you are new to the package, the **Data** page is the shortest route
to a working example: it shows how to load a single-cell reference and a
matching bulk matrix, which is the input every method expects.

## The shape of a deconvolution workflow

All four ported methods follow the same four stages, which is worth
internalising before reading any individual method page:

1.  **Simulate** pseudo-bulk training data from an annotated single-cell
    reference, with known cell-type proportions.
2.  **Process** the simulated and real bulk matrices onto a shared gene
    set, with variance filtering, log transformation, and per-sample
    scaling.
3.  **Train** a model to predict proportions from expression.
4.  **Predict** proportions for the real bulk samples.

Each method exposes those four stages as separate functions, plus a
wrapper that runs all of them in one call. Once you have seen one
method, the others read the same way — the differences are in the model
and in what happens during training, not in the surrounding workflow.

## Citation

If you use QUASAR, cite the QUASAR paper. If you use one of the R
implementations of Scaden, TAPE, DISSECT, or OmicsTweezer, please cite
**both** the QUASAR paper — for the R implementation — **and** the
original publication of the method. Full references are given on each
method’s page.

``` r

citation("quasar")
```
