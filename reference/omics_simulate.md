# Simulate pseudobulk training data for OmicsTweezer

Generates artificial pseudobulk samples from single-cell reference data
and returns them in the genes x samples format used by the OmicsTweezer
workflow.

## Usage

``` r
omics_simulate(
  sc_data,
  d_prior = NULL,
  n = 500L,
  samplenum = 5000L,
  random_state = NULL,
  sparse = FALSE,
  sparse_prob = 0.5,
  rare = FALSE,
  rare_percentage = 0.4,
  celltype_col = "CellType",
  assay = "RNA",
  slot = "counts"
)
```

## Arguments

- sc_data:

  A single-cell reference object. Supported inputs are a \`Seurat\`
  object, a \`SingleCellExperiment\` object, a matrix, or a
  \`data.frame\`.

- d_prior:

  Optional numeric vector of Dirichlet concentration parameters. If
  \`NULL\`, a symmetric Dirichlet prior of ones is used.

- n:

  Integer scalar. Target number of single cells per simulated pseudobulk
  before rounding to integer cell counts.

- samplenum:

  Integer scalar. Number of pseudobulk samples to generate.

- random_state:

  Optional integer scalar used to control reproducible simulation.

- sparse:

  Logical scalar. If \`TRUE\`, a subset of simulated samples is forced
  to contain zero fractions for a subset of cell types.

- sparse_prob:

  Numeric scalar. Controls both the proportion of sparse samples and the
  proportion of cell types zeroed within those samples.

- rare:

  Logical scalar. If \`TRUE\`, a subset of cell types is perturbed to
  have very small fractions in a subset of samples.

- rare_percentage:

  Numeric scalar in \`\[0, 1\]\`. Fraction of cell types to treat as
  rare when \`rare = TRUE\`.

- celltype_col:

  Character scalar giving the metadata column that contains cell-type
  labels for \`Seurat\` and \`SingleCellExperiment\` input.

- assay:

  Character scalar giving the assay name to extract from a \`Seurat\`
  object.

- slot:

  Character scalar giving the assay slot or assay name to extract. For
  \`Seurat\`, this is interpreted as a slot or layer within \`assay\`.
  For \`SingleCellExperiment\`, this is interpreted as an assay name.

## Value

A named list with components:

- X:

  Numeric matrix with genes in rows and simulated pseudobulk samples in
  columns.

- obs:

  Data frame with samples in rows and realized cell-type proportions in
  columns.

- var:

  Data frame indexed by gene names.

## Details

Counts are sampled with replacement within each cell type according to
Dirichlet-generated mixture proportions. Optional sparse and rare-cell
perturbations can be applied to mimic heterogeneous cellular
compositions.

The simulation follows the OmicsTweezer-style pseudobulk generation
strategy. Single-cell counts are first grouped by cell type, then
sampled according to Dirichlet-generated proportions and aggregated into
pseudobulk profiles.
