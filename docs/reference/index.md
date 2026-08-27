# Package index

## QUASAR

Uncertainty-aware cell-type deconvolution. Returns Dirichlet-based point
estimates together with confidence and prediction intervals.

- [`quasar_sim_bulk()`](https://sergejruff.github.io/QUASAR/reference/quasar_sim_bulk.md)
  : Simulate pseudo-bulk or spatial transcriptomics profiles from
  single-cell data

## Scaden

Native R implementation of Scaden (Menden et al., 2020). A three-model
multilayer perceptron ensemble trained on simulated pseudo-bulks.

- [`scaden_sim_pb()`](https://sergejruff.github.io/QUASAR/reference/scaden_sim_pb.md)
  : Simulate pseudobulk training data for Scaden
- [`scaden_process()`](https://sergejruff.github.io/QUASAR/reference/scaden_process.md)
  : Process simulated and bulk data for Scaden
- [`scaden()`](https://sergejruff.github.io/QUASAR/reference/scaden.md)
  : Train a Scaden ensemble in torch
- [`scaden_predict()`](https://sergejruff.github.io/QUASAR/reference/scaden_predict.md)
  : Predict cell-type proportions with a trained Scaden ensemble

## TAPE

Native R implementation of TAPE (Chen et al., 2022). An autoencoder with
a tissue-adaptive refinement stage, in overall or high-resolution mode.

- [`tape_simulate()`](https://sergejruff.github.io/QUASAR/reference/tape_simulate.md)
  : Simulate TAPE training pseudobulks from single-cell data
- [`tape_process()`](https://sergejruff.github.io/QUASAR/reference/tape_process.md)
  : Process simulated and bulk data for TAPE
- [`tape_train()`](https://sergejruff.github.io/QUASAR/reference/tape_train.md)
  : Train the TAPE autoencoder
- [`tape_predict()`](https://sergejruff.github.io/QUASAR/reference/tape_predict.md)
  : Predict cell fractions and signature matrices with TAPE
- [`tape()`](https://sergejruff.github.io/QUASAR/reference/tape.md) :
  Run the full TAPE workflow
- [`tape_clone_model()`](https://sergejruff.github.io/QUASAR/reference/tape_clone_model.md)
  : Clone a trained TAPE model

## DISSECT

Native R implementation of DISSECT (Khatri et al., 2024).
Semi-supervised consistency regularization for cell-type fractions and
cell-type-specific expression.

- [`dissect_simulate()`](https://sergejruff.github.io/QUASAR/reference/dissect_simulate.md)
  : Simulate DISSECT training mixtures from single-cell data
- [`dissect_process()`](https://sergejruff.github.io/QUASAR/reference/dissect_process.md)
  : Prepare DISSECT input matrices
- [`dissect_prop()`](https://sergejruff.github.io/QUASAR/reference/dissect_prop.md)
  : Estimate cell-type proportions with DISSECT
- [`dissect_expr()`](https://sergejruff.github.io/QUASAR/reference/dissect_expr.md)
  : Estimate cell-type-specific expression with DISSECT
- [`dissect()`](https://sergejruff.github.io/QUASAR/reference/dissect.md)
  : Run the full DISSECT workflow

## OmicsTweezer

Native R implementation of OmicsTweezer (Yang et al., 2025). Domain
adaptation between simulated pseudo-bulks and real bulk samples.

- [`omics_simulate()`](https://sergejruff.github.io/QUASAR/reference/omics_simulate.md)
  : Simulate pseudobulk training data for OmicsTweezer
- [`omics_process()`](https://sergejruff.github.io/QUASAR/reference/omics_process.md)
  : Process simulated and bulk data for OmicsTweezer
- [`omics_create_model()`](https://sergejruff.github.io/QUASAR/reference/omics_create_model.md)
  : Create an OmicsTweezer model
- [`omics_train()`](https://sergejruff.github.io/QUASAR/reference/omics_train.md)
  : Train one OmicsTweezer model
- [`omics_predict()`](https://sergejruff.github.io/QUASAR/reference/omics_predict.md)
  : Predict cell-type proportions with one OmicsTweezer model
- [`omics_tweezer()`](https://sergejruff.github.io/QUASAR/reference/omics_tweezer.md)
  : Run the complete OmicsTweezer workflow
- [`omics_set_seed()`](https://sergejruff.github.io/QUASAR/reference/omics_set_seed.md)
  : Set reproducible random seeds for OmicsTweezer

## Evaluation

Metrics for comparing estimated proportions against known ground truth.

- [`quasar_prop_metrics()`](https://sergejruff.github.io/QUASAR/reference/quasar_prop_metrics.md)
  : Calculate proportion-estimation performance metrics

## Normalisation

Gene filtering and sample-wise normalisation shared across methods.

- [`normalize_bulks()`](https://sergejruff.github.io/QUASAR/reference/normalize_bulks.md)
  : Normalize reference and target bulk matrices

## Interoperability

Reading and writing AnnData `.h5ad` files.

- [`quasar_h5adimport()`](https://sergejruff.github.io/QUASAR/reference/quasar_h5adimport.md)
  : Import an .h5ad file as a Seurat or SingleCellExperiment object
- [`quasar_h5adexporter()`](https://sergejruff.github.io/QUASAR/reference/quasar_h5adexporter.md)
  : Export a Seurat, SingleCellExperiment or bulk matrix to .h5ad
- [`quasar_h5_read_dataframe()`](https://sergejruff.github.io/QUASAR/reference/quasar_h5_read_dataframe.md)
  : Read an anndata dataframe group (obs / var) into a data.frame
- [`quasar_h5_read_matrix()`](https://sergejruff.github.io/QUASAR/reference/quasar_h5_read_matrix.md)
  : Read an anndata matrix group / dataset (X, raw/X, layers, obsm) as
  an R matrix
