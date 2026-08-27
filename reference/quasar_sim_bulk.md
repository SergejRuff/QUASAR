# Simulate pseudo-bulk or spatial transcriptomics profiles from single-cell data

Generates pseudo-bulk expression profiles by sampling cells from a
single-cell reference according to Dirichlet-drawn cell-type fractions.
With `mode = "spatial"` the same machinery produces spot-level profiles
instead, pooling only a handful of cells per spot and restricting each
spot to a small number of cell types. Optionally returns per-cell-type
signature profiles and a global signature matrix.

## Usage

``` r
quasar_sim_bulk(
  ...,
  n_bulk_samples = NULL,
  cells_per_bulk = 500,
  mode = c("bulk", "spatial"),
  cells_per_spot_range = c(5L, 11L),
  celltypes_per_spot_range = c(1L, 5L),
  spatial_background = 1e-06,
  cell_type_column = "cell_type",
  patient_id_column = NULL,
  select_ct = NULL,
  dirich_alpha = 1,
  seet = 1,
  sparse = FALSE,
  sparse_prob = 0.5,
  rare = FALSE,
  rare_percentage = 0.4,
  return_used_samples = FALSE,
  return_patient_metadata = FALSE,
  return_signature_matrix = FALSE,
  verbose = TRUE
)
```

## Arguments

- ...:

  The single-cell reference, supplied in one of two forms:

  - A single object: a `Seurat` object (counts taken from the `"RNA"`
    assay) or a `SingleCellExperiment` (counts taken from `counts()`).

  - Two objects: a genes \\\times\\ cells count matrix followed by a
    cell-metadata `data.frame` whose row names match the column names of
    the count matrix.

- n_bulk_samples:

  Integer number of pseudo-bulk samples to generate. If `NULL`
  (default), `1000 * (number of cell types)` samples are generated.

- cells_per_bulk:

  Integer number of cells pooled into each pseudo-bulk sample. Default
  `500`. Ignored when `mode = "spatial"`, where the pool size is drawn
  per spot from `cells_per_spot_range`.

- mode:

  Character scalar selecting the simulation type. `"bulk"` (default)
  pools `cells_per_bulk` cells per sample from a symmetric Dirichlet
  over all cell types. `"spatial"` emulates spot-level transcriptomics:
  each spot keeps only a few cell types and pools a small, randomly
  drawn number of cells.

- cells_per_spot_range:

  Integer vector of length two giving the inclusive range from which the
  number of cells per spot is drawn uniformly when `mode = "spatial"`.
  Default `c(5, 11)`.

- celltypes_per_spot_range:

  Integer vector of length two giving the inclusive range from which the
  number of cell types present in a spot is drawn uniformly when
  `mode = "spatial"`. Default `c(1, 5)`, capped at the number of
  available cell types.

- spatial_background:

  Numeric scalar giving the Dirichlet concentration assigned to cell
  types that are *not* selected for a spot when `mode = "spatial"`. A
  small positive value leaves trace amounts rather than exact zeros,
  which is what the reference implementation does. Default `1e-6`.

- cell_type_column:

  Name of the metadata column holding cell-type labels. Default
  `"cell_type"`.

- patient_id_column:

  Optional name of a metadata column holding patient/donor IDs. When
  supplied, each pseudo-bulk is drawn from a single donor where that
  donor has cells of the required type (falling back to the full
  cell-type pool otherwise). Default `NULL` (no patient structure).

- select_ct:

  Optional character vector restricting which cell types are used (and
  fixing their order in the output). Default `NULL` (all cell types).

- dirich_alpha:

  Concentration parameter of the symmetric Dirichlet used to draw
  cell-type fractions. Smaller values give more skewed mixtures. Default
  `1` (uniform over the simplex).

- seet:

  Integer random seed (passed to `set.seed`) for reproducibility.
  Default `1`.

- sparse:

  Logical; if `TRUE`, randomly zero out a fraction of cell-type entries
  before renormalisation, producing samples in which some cell types are
  absent. Default `FALSE`.

- sparse_prob:

  Probability that a given cell-type fraction is dropped when
  `sparse = TRUE`. Default `0.5`.

- rare:

  Logical; if `TRUE`, force a random subset of cell-type fractions to
  small values in `[0, 0.03]` before renormalisation, to emulate rare
  populations. Default `FALSE`.

- rare_percentage:

  Probability that a given cell-type fraction is made rare when
  `rare = TRUE`. Default `0.4`.

- return_used_samples:

  Logical; if `TRUE`, include the per-sample, per-cell-type cell indices
  that were drawn. Default `FALSE`.

- return_patient_metadata:

  Logical; if `TRUE` and `patient_id_column` is set, include a data
  frame mapping each pseudo-bulk to its source donor. Default `FALSE`.

- return_signature_matrix:

  Logical; if `TRUE`, also compute per-cell-type signature profiles (one
  genes \\\times\\ bulk matrix per cell type) and a global genes
  \\\times\\ cell-types signature matrix. Default `FALSE`.

- verbose:

  Logical; if `TRUE` (default), print the header, progress bars, and
  timing summary. If `FALSE`, nothing is printed.

## Value

A named list containing:

- `bulk_expression_profiles`:

  Genes \\\times\\ `n_bulk_samples` matrix of summed pseudo-bulk counts.

- `ground_truth_proportions`:

  `n_bulk_samples` \\\times\\ cell-types matrix of realized proportions
  (rows sum to 1).

- `cells_per_sample`:

  Integer vector giving the realized number of cells pooled into each
  sample. Constant in bulk mode, variable in spatial mode.

- `used_samples_by_ct`:

  (if `return_used_samples`) list of length `n_bulk_samples`, each a
  per-cell-type list of drawn cell indices.

- `bulk_patient_metadata`:

  (if `return_patient_metadata` and patient mode) data frame mapping
  `sample_id` to `patient_id`.

- `bulk_signature_profiles`:

  (if `return_signature_matrix`) list of per-cell-type genes \\\times\\
  bulk mean-expression matrices.

- `global_signature_matrix`:

  (if `return_signature_matrix`) genes \\\times\\ cell-types matrix of
  mean signatures.

- `timing`:

  list with `pseudobulk_seconds`, `signature_seconds`, and
  `total_seconds`.

## Details

In `mode = "bulk"`, cell-type fractions are drawn from a symmetric
Dirichlet, optionally modified by the `sparse`/`rare` masks,
renormalised per sample, and turned into integer cell counts by
`floor(fraction * cells_per_bulk)` (each bulk is guaranteed at least one
cell).

In `mode = "spatial"`, each spot first draws how many cell types it
contains, then a Dirichlet whose concentration is `dirich_alpha` for the
selected types and `spatial_background` for the rest. The number of
cells in the spot is drawn uniformly from `cells_per_spot_range` and
allocated by `round(fraction * n_cells)` rather than `floor`, because
flooring a handful of cells would empty most spots. The `sparse` and
`rare` arguments are ignored in this mode, since sparsity is already
imposed by the per-spot cell-type selection.

In both modes cells are then sampled with replacement from the reference
and summed, so the returned `ground_truth_proportions` reflect the
\*realized\* allocations rather than the raw Dirichlet draws.

## Examples

``` r

## Tiny synthetic reference: 1000 genes, 3 cell types, 100 cells each
set.seed(1)
n_genes        <- 1000
cell_types     <- c("Tcell", "Bcell", "Mono")
cells_per_type <- 100
n_cells        <- length(cell_types) * cells_per_type

counts <- matrix(
  rpois(n_genes * n_cells, lambda = 5),
  nrow = n_genes, ncol = n_cells
)
rownames(counts) <- paste0("gene_", seq_len(n_genes))
colnames(counts) <- paste0("cell_", seq_len(n_cells))

meta <- data.frame(
  cell_type = rep(cell_types, each = cells_per_type),
  row.names = colnames(counts)
)

## Generate 50 pseudo-bulks of 100 cells each, with signatures
res <- quasar_sim_bulk(
  counts, meta,
  n_bulk_samples          = 50,
  cells_per_bulk          = 100,
  return_signature_matrix = TRUE,
  verbose                 = TRUE
)
#> -------------------------------
#> Generating pseudo-bulks
#> -------------------------------
#> Creating pseudo-bulks
#> [                              ]   0.00%
#> Samples       : 0 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█                             ]   2.00%
#> Samples       : 1 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█                             ]   4.00%
#> Samples       : 2 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [██                            ]   6.00%
#> Samples       : 3 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [██                            ]   8.00%
#> Samples       : 4 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [███                           ]  10.00%
#> Samples       : 5 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [████                          ]  12.00%
#> Samples       : 6 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [████                          ]  14.00%
#> Samples       : 7 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█████                         ]  16.00%
#> Samples       : 8 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█████                         ]  18.00%
#> Samples       : 9 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [██████                        ]  20.00%
#> Samples       : 10 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [███████                       ]  22.00%
#> Samples       : 11 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [███████                       ]  24.00%
#> Samples       : 12 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [████████                      ]  26.00%
#> Samples       : 13 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [████████                      ]  28.00%
#> Samples       : 14 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█████████                     ]  30.00%
#> Samples       : 15 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [██████████                    ]  32.00%
#> Samples       : 16 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [██████████                    ]  34.00%
#> Samples       : 17 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [███████████                   ]  36.00%
#> Samples       : 18 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [███████████                   ]  38.00%
#> Samples       : 19 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [████████████                  ]  40.00%
#> Samples       : 20 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█████████████                 ]  42.00%
#> Samples       : 21 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█████████████                 ]  44.00%
#> Samples       : 22 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [██████████████                ]  46.00%
#> Samples       : 23 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [██████████████                ]  48.00%
#> Samples       : 24 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [███████████████               ]  50.00%
#> Samples       : 25 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [████████████████              ]  52.00%
#> Samples       : 26 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [████████████████              ]  54.00%
#> Samples       : 27 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█████████████████             ]  56.00%
#> Samples       : 28 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█████████████████             ]  58.00%
#> Samples       : 29 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [██████████████████            ]  60.00%
#> Samples       : 30 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [███████████████████           ]  62.00%
#> Samples       : 31 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [███████████████████           ]  64.00%
#> Samples       : 32 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [████████████████████          ]  66.00%
#> Samples       : 33 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [████████████████████          ]  68.00%
#> Samples       : 34 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█████████████████████         ]  70.00%
#> Samples       : 35 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [██████████████████████        ]  72.00%
#> Samples       : 36 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [██████████████████████        ]  74.00%
#> Samples       : 37 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [███████████████████████       ]  76.00%
#> Samples       : 38 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [███████████████████████       ]  78.00%
#> Samples       : 39 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [████████████████████████      ]  80.00%
#> Samples       : 40 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█████████████████████████     ]  82.00%
#> Samples       : 41 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█████████████████████████     ]  84.00%
#> Samples       : 42 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [██████████████████████████    ]  86.00%
#> Samples       : 43 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [██████████████████████████    ]  88.00%
#> Samples       : 44 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [███████████████████████████   ]  90.00%
#> Samples       : 45 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [████████████████████████████  ]  92.00%
#> Samples       : 46 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [████████████████████████████  ]  94.00%
#> Samples       : 47 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█████████████████████████████ ]  96.00%
#> Samples       : 48 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [█████████████████████████████ ]  98.00%
#> Samples       : 49 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating pseudo-bulks
#> [██████████████████████████████] 100.00%
#> Samples       : 50 / 50
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Building signature matrices
#> [                              ]   0.00%
#> Steps         : 0 / 6
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Building signature matrices
#> [█████                         ]  16.67%
#> Steps         : 1 / 6
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Building signature matrices
#> [██████████                    ]  33.33%
#> Steps         : 2 / 6
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Building signature matrices
#> [███████████████               ]  50.00%
#> Steps         : 3 / 6
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Building signature matrices
#> [████████████████████          ]  66.67%
#> Steps         : 4 / 6
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Building signature matrices
#> [█████████████████████████     ]  83.33%
#> Steps         : 5 / 6
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Building signature matrices
#> [██████████████████████████████] 100.00%
#> Steps         : 6 / 6
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Generated 50 pseudobulks (100 cells each) from 1000 genes
#> Pseudobulk creation time: 0.12 secs
#> Computed bulk_signature_profiles (list of genes×bulk matrices per cell type)
#> Signature matrix time: 0.01 secs
#> Returned global_signature_matrix (genes × celltypes)
#> Total quasar_sim_bulk time: 0.13 secs
#> 

dim(res$bulk_expression_profiles)   # 1000 x 50
#> [1] 1000   50
head(res$ground_truth_proportions)  # rows sum to 1
#>               Tcell     Bcell       Mono
#> sample_1 0.04081633 0.4897959 0.46938776
#> sample_2 0.25252525 0.3838384 0.36363636
#> sample_3 0.71717172 0.2222222 0.06060606
#> sample_4 0.16161616 0.3333333 0.50505051
#> sample_5 0.02040816 0.5408163 0.43877551
#> sample_6 0.43877551 0.3877551 0.17346939
dim(res$global_signature_matrix)    # 1000 x 3
#> [1] 1000    3

## Spatial spots: few cells and few cell types per spot
spots <- quasar_sim_bulk(
  counts, meta,
  n_bulk_samples = 200,
  mode           = "spatial",
  verbose        = TRUE
)
#> -------------------------------
#> Generating spatial spots
#> -------------------------------
#> Creating spatial spots
#> [                              ]   0.00%
#> Spots         : 0 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [                              ]   1.00%
#> Spots         : 2 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█                             ]   2.00%
#> Spots         : 4 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█                             ]   3.00%
#> Spots         : 6 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█                             ]   4.00%
#> Spots         : 8 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██                            ]   5.00%
#> Spots         : 10 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██                            ]   6.00%
#> Spots         : 12 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██                            ]   7.00%
#> Spots         : 14 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██                            ]   8.00%
#> Spots         : 16 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███                           ]   9.00%
#> Spots         : 18 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███                           ]  10.00%
#> Spots         : 20 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███                           ]  11.00%
#> Spots         : 22 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████                          ]  12.00%
#> Spots         : 24 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████                          ]  13.00%
#> Spots         : 26 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████                          ]  14.00%
#> Spots         : 28 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████                          ]  15.00%
#> Spots         : 30 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████                         ]  16.00%
#> Spots         : 32 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████                         ]  17.00%
#> Spots         : 34 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████                         ]  18.00%
#> Spots         : 36 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████                        ]  19.00%
#> Spots         : 38 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████                        ]  20.00%
#> Spots         : 40 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████                        ]  21.00%
#> Spots         : 42 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████                       ]  22.00%
#> Spots         : 44 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████                       ]  23.00%
#> Spots         : 46 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████                       ]  24.00%
#> Spots         : 48 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████                      ]  25.00%
#> Spots         : 50 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████                      ]  26.00%
#> Spots         : 52 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████                      ]  27.00%
#> Spots         : 54 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████                      ]  28.00%
#> Spots         : 56 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████                     ]  29.00%
#> Spots         : 58 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████                     ]  30.00%
#> Spots         : 60 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████                     ]  31.00%
#> Spots         : 62 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████                    ]  32.00%
#> Spots         : 64 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████                    ]  33.00%
#> Spots         : 66 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████                    ]  34.00%
#> Spots         : 68 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████                    ]  35.00%
#> Spots         : 70 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████                   ]  36.00%
#> Spots         : 72 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████                   ]  37.00%
#> Spots         : 74 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████                   ]  38.00%
#> Spots         : 76 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████                  ]  39.00%
#> Spots         : 78 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████                  ]  40.00%
#> Spots         : 80 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████                  ]  41.00%
#> Spots         : 82 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████                 ]  42.00%
#> Spots         : 84 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████                 ]  43.00%
#> Spots         : 86 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████                 ]  44.00%
#> Spots         : 88 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████                ]  45.00%
#> Spots         : 90 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████                ]  46.00%
#> Spots         : 92 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████                ]  47.00%
#> Spots         : 94 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████                ]  48.00%
#> Spots         : 96 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████████               ]  49.00%
#> Spots         : 98 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████████               ]  50.00%
#> Spots         : 100 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████████               ]  51.00%
#> Spots         : 102 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████              ]  52.00%
#> Spots         : 104 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████              ]  53.00%
#> Spots         : 106 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████              ]  54.00%
#> Spots         : 108 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████              ]  55.00%
#> Spots         : 110 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████████             ]  56.00%
#> Spots         : 112 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████████             ]  57.00%
#> Spots         : 114 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████████             ]  58.00%
#> Spots         : 116 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████████            ]  59.00%
#> Spots         : 118 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████████            ]  60.00%
#> Spots         : 120 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████████            ]  61.00%
#> Spots         : 122 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████████████           ]  62.00%
#> Spots         : 124 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████████████           ]  63.00%
#> Spots         : 126 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████████████           ]  64.00%
#> Spots         : 128 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████████          ]  65.00%
#> Spots         : 130 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████████          ]  66.00%
#> Spots         : 132 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████████          ]  67.00%
#> Spots         : 134 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████████          ]  68.00%
#> Spots         : 136 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████████████         ]  69.00%
#> Spots         : 138 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████████████         ]  70.00%
#> Spots         : 140 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████████████         ]  71.00%
#> Spots         : 142 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████████████        ]  72.00%
#> Spots         : 144 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████████████        ]  73.00%
#> Spots         : 146 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████████████        ]  74.00%
#> Spots         : 148 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████████████        ]  75.00%
#> Spots         : 150 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████████████████       ]  76.00%
#> Spots         : 152 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████████████████       ]  77.00%
#> Spots         : 154 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████████████████       ]  78.00%
#> Spots         : 156 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████████████      ]  79.00%
#> Spots         : 158 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████████████      ]  80.00%
#> Spots         : 160 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████████████      ]  81.00%
#> Spots         : 162 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████████████████     ]  82.00%
#> Spots         : 164 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████████████████     ]  83.00%
#> Spots         : 166 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████████████████     ]  84.00%
#> Spots         : 168 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████████████████    ]  85.00%
#> Spots         : 170 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████████████████    ]  86.00%
#> Spots         : 172 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████████████████    ]  87.00%
#> Spots         : 174 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████████████████    ]  88.00%
#> Spots         : 176 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████████████████████   ]  89.00%
#> Spots         : 178 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████████████████████   ]  90.00%
#> Spots         : 180 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [███████████████████████████   ]  91.00%
#> Spots         : 182 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████████████████  ]  92.00%
#> Spots         : 184 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████████████████  ]  93.00%
#> Spots         : 186 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████████████████  ]  94.00%
#> Spots         : 188 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [████████████████████████████  ]  95.00%
#> Spots         : 190 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████████████████████ ]  96.00%
#> Spots         : 192 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████████████████████ ]  97.00%
#> Spots         : 194 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [█████████████████████████████ ]  98.00%
#> Spots         : 196 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████████████████████]  99.00%
#> Spots         : 198 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Creating spatial spots
#> [██████████████████████████████] 100.00%
#> Spots         : 200 / 200
#> Elapsed       : 00:00
#> ETA           : 00:00
#> Generated 200 spatial spots (4-12 cells each, median 8) from 1000 genes
#> Median cell types per spot: 2
#> Spot creation time: 0.22 secs
#> Total quasar_sim_bulk time: 0.22 secs
#> 

dim(spots$bulk_expression_profiles)  # 1000 x 200
#> [1] 1000  200
range(spots$cells_per_sample)        # within cells_per_spot_range
#> [1]  4 12
table(rowSums(spots$ground_truth_proportions > 0))
#> 
#>  1  2  3 
#> 74 85 41 

```
