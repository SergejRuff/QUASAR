# Simulate pseudobulk training data for Scaden

Faithful port of \`generate_simulated_data()\`. Counts are sampled with
replacement within each cell type according to Dirichlet-generated
mixture proportions. Optionally merges a set of cell types into a single
unknown class, following Scaden's post-simulation merge semantics.

## Usage

``` r
scaden_sim_pb(
  sc_data,
  metadata = NULL,
  celltypes = NULL,
  celltype_col = NULL,
  assay = "RNA",
  slot = "counts",
  d_prior = NULL,
  n = 500,
  samplenum = 5000,
  seed = NULL,
  sparse = TRUE,
  sparse_prob = 0.5,
  rare = FALSE,
  rare_percentage = 0.4,
  unknown_celltypes = character(0),
  unknown_label = "Unknown",
  verbose = TRUE
)
```

## Arguments

- sc_data:

  A \`Seurat\`, \`SingleCellExperiment\`, matrix, \`data.frame\`, or
  sparse \`Matrix\` reference.

- metadata:

  Optional metadata for matrix-like input.

- celltypes:

  Optional character vector of cell-type labels for matrix-like input.
  Used only when \`metadata\` is not supplied.

- celltype_col:

  Character scalar giving the metadata column holding cell-type labels.

- assay:

  Character scalar giving the Seurat assay name.

- slot:

  Character scalar giving the assay slot/layer (Seurat) or assay name
  (SingleCellExperiment).

- d_prior:

  Optional numeric vector of Dirichlet concentration parameters. If
  \`NULL\`, a symmetric prior of ones is used. Length must equal the
  number of cell types \*before\* unknown merging.

- n:

  Integer scalar. Target cells per pseudobulk before flooring.

- samplenum:

  Integer scalar. Number of pseudobulks to generate.

- seed:

  Optional integer scalar for reproducible simulation.

- sparse:

  Logical scalar. Zero out a subset of cell types in a subset of
  samples.

- sparse_prob:

  Numeric scalar in \`\[0, 1)\`. Controls both the proportion of sparse
  samples and the proportion of cell types zeroed within them.

- rare:

  Logical scalar. Perturb a subset of cell types to very small fractions
  in a subset of samples.

- rare_percentage:

  Numeric scalar in \`\[0, 1\]\`. Fraction of cell types treated as rare
  when \`rare = TRUE\`.

- unknown_celltypes:

  Character vector of cell-type labels to merge into a single unknown
  class. Defaults to \`character(0)\` (no merging).

- unknown_label:

  Character scalar used as the merged label.

- verbose:

  Logical scalar. Print progress.

## Value

A named list with \`train_x\` (samples x genes), \`train_y\` (samples x
cell types), \`cell_num\`, \`celltypes\`, \`genes\`, \`elapsed\`, and
\`settings\`.

## Details

Unknown merging happens \*\*after\*\* simulation: each original cell
type gets its own Dirichlet component and is sampled independently, and
the proportion columns of the unknown types are then summed into one
column. This mirrors Scaden's \`y.groupby(y.columns, axis=1).sum()\` and
differs from relabelling cells before simulation, which would give the
merged class a single Dirichlet component and pool its cells uniformly.
