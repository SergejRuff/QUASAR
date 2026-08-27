# Run the complete OmicsTweezer workflow

Runs OmicsTweezer cell-type deconvolution from a single-cell reference
and real bulk expression profiles. The workflow simulates pseudobulks,
processes simulated and real bulk data, trains three neural network
architectures, and averages their predictions.

## Usage

``` r
omics_tweezer(
  sc_data,
  real_bulk,
  ot_weight = 1,
  scale_minmax = FALSE,
  samplenum = 5000L,
  n_cells_per_bulk = 500L,
  sparse = FALSE,
  variance_threshold = 0.98,
  scaler = c("ss", "mms", "none"),
  batch_size = 128L,
  epochs = 30L,
  learning_rate = 1e-04,
  seed = 2021L,
  celltype_col = "CellType",
  assay = "RNA",
  slot = "counts",
  verbose = TRUE,
  device = c("auto", "cpu", "cuda"),
  cuda_index = NULL
)
```

## Arguments

- sc_data:

  A single-cell reference object. Supported inputs are a \`Seurat\`
  object, a \`SingleCellExperiment\` object, a matrix, or a
  \`data.frame\`.

- real_bulk:

  Numeric matrix or \`data.frame\` with genes in rows and real bulk
  samples in columns.

- ot_weight:

  Numeric scalar. Weight applied to the Wasserstein-style domain
  adaptation loss.

- scale_minmax:

  Logical scalar. If \`TRUE\`, applies an outer global min-max scaling
  step to simulated and real bulk count matrices before processing.

- samplenum:

  Integer scalar. Number of simulated pseudobulk samples.

- n_cells_per_bulk:

  Integer scalar. Target number of single cells per simulated
  pseudobulk.

- sparse:

  Logical scalar. If \`TRUE\`, applies sparse perturbation during
  pseudobulk simulation.

- variance_threshold:

  Numeric scalar. Fraction of genes used to define the variance cutoff
  during processing.

- scaler:

  Character scalar specifying the per-sample scaling method. One of
  \`"ss"\`, \`"mms"\`, or \`"none"\`.

- batch_size:

  Integer scalar. Minibatch size for training and prediction.

- epochs:

  Integer scalar. Number of training epochs for each architecture.

- learning_rate:

  Numeric scalar. Learning rate passed to Adam.

- seed:

  Integer scalar random seed.

- celltype_col:

  Character scalar giving the metadata column containing cell-type
  labels for \`Seurat\` and \`SingleCellExperiment\` input.

- assay:

  Character scalar giving the assay name to extract from a \`Seurat\`
  object.

- slot:

  Character scalar giving the assay slot or assay name to extract.

- verbose:

  Logical scalar. If \`TRUE\`, prints workflow progress messages.

- device:

  Character scalar. One of \`"auto"\`, \`"cpu"\`, or \`"cuda"\`. If
  \`"auto"\`, CUDA is used when available, otherwise CPU is used.

- cuda_index:

  Optional integer scalar giving the CUDA device index to use when
  \`device = "cuda"\`. CUDA indices are zero-based.

## Value

A named list with components:

- pred:

  Data frame of averaged predicted cell-type proportions.

- per_model:

  Named list of predictions from the \`m256\`, \`m512\`, and \`m1024\`
  architectures.

- simudata:

  Simulated pseudobulk data returned by \[omics_simulate()\].

- processed:

  Processed matrices returned by \[omics_process()\].

- models:

  Named list of trained OmicsTweezer models.

- loss_history:

  Data frame containing epoch-level training losses for all
  architectures.

## Details

The ensemble uses three architectures: \`m256\`, \`m512\`, and
\`m1024\`. Predictions are generated independently for each architecture
and then averaged sample-wise.
