#' Simulate pseudobulk training data for Scaden
#'
#' Generates artificial pseudobulk samples from single-cell reference data and
#' returns them in a format directly compatible with [scaden()].
#'
#' The function accepts `Seurat`, `SingleCellExperiment`, base matrices,
#' `data.frame`s, and sparse `Matrix` objects. Counts are sampled with
#' replacement within each cell type according to Dirichlet-generated mixture
#' proportions. Optional sparse and rare-cell perturbations can be applied to
#' mimic more heterogeneous compositions.
#'
#' @param sc_data A single-cell reference object. Supported inputs are a
#'   `Seurat` object, a `SingleCellExperiment` object, a matrix, a
#'   `data.frame`, or a sparse `Matrix`.
#' @param metadata Optional metadata for matrix-like input. Must contain the
#'   column specified by `celltype_col` when supplied.
#' @param celltypes Optional character vector of cell-type labels for matrix-like
#'   input. Used only when `metadata` is not supplied.
#' @param celltype_col Character scalar giving the metadata column that contains
#'   cell-type labels. Must be provided for `Seurat`, `SingleCellExperiment`,
#'   or matrix input with `metadata`.
#' @param assay Character scalar giving the assay name to extract from a
#'   `Seurat` object.
#' @param slot Character scalar giving the assay slot or assay name to extract.
#'   For `Seurat`, this is interpreted as a slot within `assay`. For
#'   `SingleCellExperiment`, this is interpreted as an assay name.
#' @param d_prior Optional numeric vector of Dirichlet concentration parameters.
#'   If `NULL`, a symmetric Dirichlet prior of ones is used.
#' @param n Integer scalar. Target number of single cells per simulated
#'   pseudobulk before rounding to integer cell counts.
#' @param samplenum Integer scalar. Number of pseudobulk samples to generate.
#' @param seed Optional integer scalar used to control reproducible simulation.
#' @param sparse Logical scalar. If `TRUE`, a subset of simulated samples is
#'   forced to contain zero fractions for a subset of cell types.
#' @param sparse_prob Numeric scalar in `[0, 1)`. Controls both the proportion
#'   of sparse samples and the proportion of cell types zeroed within those
#'   samples.
#' @param rare Logical scalar. If `TRUE`, a subset of cell types is perturbed to
#'   have very small fractions in a subset of samples.
#' @param rare_percentage Numeric scalar in `[0, 1]`. Fraction of cell types to
#'   treat as rare when `rare = TRUE`.
#' @param unknown_celltypes Character vector of cell-type labels that should be
#'   merged into a single unknown class before simulation.
#' @param unknown_label Character scalar used as the merged label for unknown
#'   cell types.
#' @param select_ct Optional character vector restricting simulation to a subset
#'   of cell types after unknown-cell merging.
#' @param verbose Logical scalar. If `TRUE`, prints progress information and
#'   elapsed time.
#'
#' @return A named list with components:
#' \describe{
#'   \item{train_x}{Numeric matrix with samples in rows and genes in columns.
#'   This is the training input expected by [scaden()].}
#'   \item{train_y}{Numeric matrix with samples in rows and cell types in
#'   columns. This is the training target expected by [scaden()].}
#'   \item{cell_num}{Integer matrix giving realized sampled cell counts per
#'   sample and cell type.}
#'   \item{celltypes}{Character vector of cell-type names used in the
#'   simulation.}
#'   \item{genes}{Character vector of gene names.}
#'   \item{elapsed}{Character scalar giving total runtime in `HH:MM:SS`
#'   format.}
#'   \item{settings}{List of simulation settings used to generate the output.}
#' }
#'
#' @details
#' The output is structured to plug directly into [scaden_process()] and
#' [scaden()]. In particular, `train_x` is returned as a samples x genes matrix
#' and `train_y` as a samples x cell-types matrix.
#' 
#' @references
#' Menden, K., Marouf, M., Oller, S., Dalmia, A., Magruder, D. S., Kloiber, K.,
#' Heutink, P., & Bonn, S. (2020). Deep learning-based cell composition analysis
#' from tissue expression profiles. \emph{Science Advances}, 6(30), eaba2619.
#'
#' TAPE project PyTorch implementation:
#' \url{https://github.com/poseidonchan/TAPE/blob/main/TAPE/model.py}
#'
#' @importFrom methods slot slotNames
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom Matrix Matrix sparseMatrix
#' @importFrom MCMCpack rdirichlet
#' @importFrom stats runif
#'
#' @examples
#' \dontrun{
#' sim <- scaden_sim_pb(
#'   sc_data = sce,
#'   celltype_col = "celltype",
#'   samplenum = 5000,
#'   n = 500,
#'   seed = 123
#' )
#' }
#'
#' @export
scaden_sim_pb <- function(
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
  unknown_celltypes = "unknown",
  unknown_label = "Unknown",
  select_ct = NULL,
  verbose = TRUE
) {

  # -- Local helpers ----------------------------------------------------------

  fmt_time <- function() {
    secs <- max(0L, as.integer(round(difftime(Sys.time(), start_time, units = "secs"))))
    sprintf("%02d:%02d:%02d", secs %/% 3600L, (secs %% 3600L) %/% 60L, secs %% 60L)
  }

  step <- 0L
  pb   <- NULL

  tick <- function(msg, ...) {
    step <<- step + 1L
    if (verbose && !is.null(pb)) {
      setTxtProgressBar(pb, step)
      cat(sprintf("  [%s] %s\n", fmt_time(), sprintf(msg, ...)))
    }
  }

  require_celltype_col <- function(meta) {
    if (is.null(celltype_col) || !nzchar(celltype_col))
      stop("`celltype_col` must be provided explicitly.")
    if (!celltype_col %in% colnames(meta))
      stop(sprintf("`celltype_col` '%s' not found in metadata.", celltype_col))
  }

  # -- Setup ------------------------------------------------------------------

  start_time <- Sys.time()

  if (!is.null(seed)) {
    stopifnot("`seed` must be a single finite number" = length(seed) == 1L && is.finite(seed))
    set.seed(as.integer(seed))
  }

  # -- Extract counts matrix (genes x cells) and cell-type vector -------------

  if (inherits(sc_data, "Seurat")) {
    stopifnot("Seurat object missing requested assay" = assay %in% names(sc_data@assays))
    assay_obj <- sc_data[[assay]]
    stopifnot("Slot not found in assay" = slot %in% methods::slotNames(assay_obj))
    counts <- as.matrix(methods::slot(assay_obj, slot))
    meta   <- sc_data@meta.data
    require_celltype_col(meta)
    cell_types <- as.character(meta[[celltype_col]])

  } else if (inherits(sc_data, "SingleCellExperiment")) {
    avail <- SummarizedExperiment::assayNames(sc_data)
    stopifnot("Assay not found in SCE object" = slot %in% avail)
    counts <- as.matrix(SummarizedExperiment::assay(sc_data, slot))
    meta   <- as.data.frame(SummarizedExperiment::colData(sc_data))
    require_celltype_col(meta)
    cell_types <- as.character(meta[[celltype_col]])

  } else if (is.matrix(sc_data) || is.data.frame(sc_data) || inherits(sc_data, "Matrix")) {
    x <- as.matrix(sc_data)

    if (!is.null(metadata)) {
      metadata <- as.data.frame(metadata)
      require_celltype_col(metadata)
      ct_vec <- as.character(metadata[[celltype_col]])

      if (ncol(x) == nrow(metadata)) {
        counts <- x
      } else if (nrow(x) == nrow(metadata)) {
        counts <- t(x)
      } else {
        stop("Dimensions of `sc_data` and `metadata` do not match.")
      }
      cell_types <- ct_vec

    } else if (!is.null(celltypes)) {
      if (length(celltypes) == ncol(x)) {
        counts <- x
      } else if (length(celltypes) == nrow(x)) {
        counts <- t(x)
      } else {
        stop("Length of `celltypes` does not match rows or columns of `sc_data`.")
      }
      cell_types <- as.character(celltypes)

    } else {
      stop("For matrix input, provide either `metadata` + `celltype_col` or `celltypes`.")
    }

  } else {
    stop("Unsupported `sc_data` type. Expected Seurat, SingleCellExperiment, matrix, data.frame, or Matrix.")
  }

  # Ensure sparse storage, default row/col names, sanitise cell types
  counts <- Matrix::Matrix(counts, sparse = TRUE)
  if (is.null(rownames(counts))) rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
  if (is.null(colnames(counts))) colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))
  cell_types <- gsub("/", "_", cell_types)
  stopifnot("Cell-type vector length != number of cells" = length(cell_types) == ncol(counts))

  # -- Apply unknown / select filters ----------------------------------------

  if (length(unknown_celltypes) > 0L)
    cell_types[cell_types %in% unknown_celltypes] <- unknown_label

  if (!is.null(select_ct)) {
    keep <- cell_types %in% select_ct
    if (!any(keep)) stop("No cells found for `select_ct`.")
    counts     <- counts[, keep, drop = FALSE]
    cell_types <- cell_types[keep]
  }

  ct_levels  <- sort(unique(cell_types))
  n_ct       <- length(ct_levels)
  if (n_ct < 1L) stop("No cell types available after filtering.")

  # -- Validate numeric parameters --------------------------------------------

  stopifnot("`n` must be a positive number"         = is.numeric(n)         && length(n) == 1L         && n > 0,
            "`samplenum` must be a positive number"  = is.numeric(samplenum) && length(samplenum) == 1L && samplenum > 0)
  n         <- as.integer(n)
  samplenum <- as.integer(samplenum)

  if (is.null(d_prior)) {
    d_prior <- rep(1, n_ct)
  } else {
    stopifnot("`d_prior` length must equal number of cell types" = length(d_prior) == n_ct,
              "All `d_prior` entries must be > 0"                = all(d_prior > 0))
  }
  if (sparse) stopifnot("`sparse_prob` must be in [0, 1)" = sparse_prob >= 0 && sparse_prob < 1)
  if (rare)   stopifnot("`rare_percentage` must be in [0, 1]" = rare_percentage >= 0 && rare_percentage <= 1)

  # -- Initialise progress bar ------------------------------------------------
  # Fixed steps: prepare | dirichlet | cell counts | sample | build + done
  # Conditional steps: +1 if sparse, +1 if rare
  n_steps <- 5L + as.integer(sparse) + as.integer(rare)

  if (verbose) {
    cat(sprintf("Simulating pseudobulks for Scaden ... (%d samples, %d cell types)\n",
                samplenum, n_ct))
    pb <- txtProgressBar(min = 0, max = n_steps, style = 3)
  }

  # -- Step: Prepare reference ------------------------------------------------

  tick("Preparing reference — %d cells, %d genes, %d cell types",
       ncol(counts), nrow(counts), n_ct)

  # -- Step: Generate Dirichlet proportions -----------------------------------

  prop <- MCMCpack::rdirichlet(samplenum, d_prior)
  prop <- prop / rowSums(prop)
  tick("Generated Dirichlet fractions")

  # -- Step (conditional): Sparse perturbation --------------------------------

  if (sparse) {
    n_sparse <- as.integer(samplenum * sparse_prob)
    n_zero   <- min(as.integer(n_ct * sparse_prob), n_ct - 1L)
    if (n_sparse > 0L && n_zero > 0L) {
      for (i in seq_len(n_sparse)) {
        prop[i, sample.int(n_ct, n_zero)] <- 0
      }
      rs <- rowSums(prop); rs[rs == 0] <- 1
      prop <- prop / rs
    }
    tick("Applied sparse perturbation (prob = %.3f)", sparse_prob)
  }

  # -- Step (conditional): Rare-cell perturbation -----------------------------

  if (rare) {
    n_rare_ct <- as.integer(n_ct * rare_percentage)
    if (n_rare_ct > 0L) {
      set.seed(0L)
      rare_idx <- sample.int(n_ct, n_rare_ct)
      prop     <- prop / rowSums(prop)
      n_rows   <- min(as.integer(0.5 * samplenum) + as.integer(rare_percentage * 0.5 * samplenum),
                      samplenum)
      for (i in seq_len(n_rows)) {
        buf <- stats::runif(n_rare_ct, 0, 0.03)
        prop[i, rare_idx] <- 0
        remain <- sum(prop[i, ])
        if (remain > 0) prop[i, ] <- (1 - sum(buf)) * prop[i, ] / remain
        prop[i, rare_idx] <- buf
      }
    }
    tick("Applied rare-cell perturbation (pct = %.3f)", rare_percentage)
  }

  # -- Step: Convert proportions to integer cell counts -----------------------

  cell_num <- matrix(as.integer(floor(n * prop)),
                     nrow = samplenum, ncol = n_ct,
                     dimnames = list(NULL, ct_levels))

  # Guarantee at least one cell per sample
  empty <- rowSums(cell_num) == 0L
  if (any(empty)) {
    best_col <- max.col(prop[empty, , drop = FALSE], ties.method = "first")
    for (ii in seq_along(best_col)) cell_num[which(empty)[ii], best_col[ii]] <- 1L
  }

  prop_realized <- cell_num / rowSums(cell_num)
  tick("Computed cell counts per sample")

  # -- Step: Sample cells -----------------------------------------------------

  ct_indices <- split(seq_along(cell_types), cell_types)[ct_levels]

  ids_all <- vector("list", n_ct)
  grp_all <- vector("list", n_ct)

  for (k in seq_along(ct_levels)) {
    total_k <- sum(cell_num[, k])
    if (total_k > 0L && length(ct_indices[[k]]) > 0L) {
      ids_all[[k]] <- sample(ct_indices[[k]], total_k, replace = TRUE)
      grp_all[[k]] <- rep.int(seq_len(samplenum), times = cell_num[, k])
    }
  }

  tick("Sampled cells for %d cell types", n_ct)

  # -- Step: Build pseudobulks ------------------------------------------------

  ids <- unlist(ids_all, use.names = FALSE)
  grp <- unlist(grp_all, use.names = FALSE)

  S    <- Matrix::sparseMatrix(i = ids, j = grp, x = 1,
                               dims = c(ncol(counts), samplenum))
  bulk <- as.matrix(counts %*% S)

  sample_names <- paste0("sample_", seq_len(samplenum))
  colnames(bulk) <- sample_names

  train_x <- t(bulk)                               # samples x genes
  train_y <- prop_realized                          # samples x cell types
  rownames(train_x) <- rownames(train_y) <- sample_names
  colnames(train_y) <- ct_levels

  tick("Simulation finished")

  if (verbose) close(pb)

  # -- Return -----------------------------------------------------------------

  list(
    train_x                  = train_x,
    train_y                  = train_y,
    cell_num                 = cell_num,
    celltypes                = ct_levels,
    genes                    = rownames(counts),
    elapsed                  = fmt_time(),
    settings = list(
      d_prior           = d_prior,
      n                 = n,
      samplenum         = samplenum,
      seed              = seed,
      sparse            = sparse,
      sparse_prob       = sparse_prob,
      rare              = rare,
      rare_percentage   = rare_percentage,
      unknown_celltypes = unknown_celltypes,
      unknown_label     = unknown_label,
      select_ct         = select_ct,
      assay             = assay,
      slot              = slot,
      celltype_col      = celltype_col
    )
  )
}































#' Process simulated and bulk data for Scaden
#'
#' Aligns simulated pseudobulk training data with optional bulk test data,
#' filters genes, applies normalization, and returns matrices ready for
#' [scaden()] and [scaden_predict()].
#'
#' @param sim_data A list returned by [scaden_sim_pb()]. Must contain at least
#'   `train_x` and `train_y`.
#' @param bulk_data Optional bulk expression matrix or `data.frame`. Gene names
#'   must be present in either row names or column names. The function will
#'   infer orientation and convert the matrix to samples x genes.
#' @param var_cutoff Numeric scalar. Genes in `bulk_data` with variance less
#'   than or equal to this threshold are removed before intersecting with the
#'   simulated training genes.
#' @param top_var_genes Optional integer scalar. If provided, only the top most
#'   variable genes in `bulk_data` are kept after applying `var_cutoff`.
#' @param normalize_train Character scalar specifying the transformation applied
#'   to the training matrix. One of `"log1p_minmax"`, `"log1p"`, or `"none"`.
#' @param normalize_test Character scalar specifying the transformation applied
#'   to the test matrix. One of `"log1p_minmax"`, `"log1p"`, or `"none"`.
#' @param verbose Logical scalar. If `TRUE`, prints progress messages and
#'   elapsed time.
#'
#' @return A named list with components:
#' \describe{
#'   \item{train_x}{Processed training matrix with samples in rows and genes in
#'   columns.}
#'   \item{train_y}{Training target matrix with samples in rows and cell types
#'   in columns.}
#'   \item{test_x}{Processed bulk test matrix with samples in rows and genes in
#'   columns, or `NULL` if `bulk_data` was not provided.}
#'   \item{genes}{Character vector of genes retained after filtering and
#'   intersection.}
#'   \item{train_samples}{Character vector of training sample names.}
#'   \item{test_samples}{Character vector of test sample names, or `NULL`.}
#'   \item{celltypes}{Character vector of cell-type names.}
#'   \item{normalize_train}{Training normalization method used.}
#'   \item{normalize_test}{Test normalization method used.}
#'   \item{elapsed}{Character scalar giving total runtime in `HH:MM:SS`
#'   format.}
#' }
#'
#' @details
#' The processing logic follows the Scaden-style setup where the simulated
#' training matrix and real bulk matrix are restricted to a shared gene set.
#' Optionally, genes can first be filtered by variance on the bulk data.
#'
#' 
#' @references
#' Menden, K., Marouf, M., Oller, S., Dalmia, A., Magruder, D. S., Kloiber, K.,
#' Heutink, P., & Bonn, S. (2020). Deep learning-based cell composition analysis
#' from tissue expression profiles. \emph{Science Advances}, 6(30), eaba2619.
#'
#' TAPE project PyTorch implementation:
#' \url{https://github.com/poseidonchan/TAPE/blob/main/TAPE/model.py}
#' 
#' @importFrom matrixStats rowMins rowMaxs colVars
#'
#' @examples
#' \dontrun{
#' proc <- scaden_process(
#'   sim_data = sim,
#'   bulk_data = bulk_mat,
#'   var_cutoff = 0.1,
#'   normalize_train = "log1p_minmax",
#'   normalize_test = "log1p_minmax"
#' )
#' }
#'
#' @export
scaden_process <- function(sim_data,
  bulk_data = NULL,
  var_cutoff = 0.1,
  top_var_genes = NULL,
  normalize_train = c("log1p_minmax", "log1p", "none"),
  normalize_test  = c("log1p_minmax", "log1p", "none"),
  verbose = TRUE) {

  normalize_train <- match.arg(normalize_train)
  normalize_test  <- match.arg(normalize_test)

  start_time <- Sys.time()

  # -- Local helpers ----------------------------------------------------------

  fmt_time <- function() {
    secs <- max(0L, as.integer(round(difftime(Sys.time(), start_time, units = "secs"))))
    sprintf("%02d:%02d:%02d", secs %/% 3600L, (secs %% 3600L) %/% 60L, secs %% 60L)
  }

  log_msg <- function(...) {
    if (verbose) cat(sprintf("[%s] %s\n", fmt_time(), sprintf(...)))
  }

  normalize <- function(x, method) {
    if (method == "none") return(x)
    x <- log1p(x)
    if (method == "log1p") return(x)
    # log1p_minmax: row-wise min-max scaling
    mins  <- matrixStats::rowMins(x)
    range <- pmax(matrixStats::rowMaxs(x) - mins, 1e-8)
    sweep(sweep(x, 1, mins, "-"), 1, range, "/")
  }

  orient_to_samples_x_genes <- function(mat, target_genes, label) {
    rn <- rownames(mat)
    cn <- colnames(mat)
    hit_rows <- if (!is.null(rn)) length(intersect(rn, target_genes)) else 0L
    hit_cols <- if (!is.null(cn)) length(intersect(cn, target_genes)) else 0L

    if (hit_cols > hit_rows) return(mat)            # already samples x genes
    if (hit_rows > hit_cols) return(t(mat))         # genes x samples → transpose

    # Tie-break on dimension match
    if (!is.null(rn) && nrow(mat) == length(target_genes)) return(t(mat))
    if (!is.null(cn) && ncol(mat) == length(target_genes)) return(mat)

    stop(sprintf(
      "Cannot infer orientation of `%s`. Ensure gene names appear in row or column names.",
      label
    ))
  }

  # -- Validate sim_data ------------------------------------------------------

  if (!is.list(sim_data) || is.null(sim_data$train_x) || is.null(sim_data$train_y))
    stop("`sim_data` must be a list with `train_x` and `train_y` (output of `scaden_sim_pb()`).")

  train_x <- as.matrix(sim_data$train_x)
  train_y <- as.matrix(sim_data$train_y)

  if (is.null(colnames(train_x)))
    stop("Training matrix must have gene names as column names.")
  if (is.null(rownames(train_x)))
    rownames(train_x) <- paste0("sample_", seq_len(nrow(train_x)))
  if (is.null(rownames(train_y)))
    rownames(train_y) <- rownames(train_x)
  if (is.null(colnames(train_y)))
    colnames(train_y) <- paste0("celltype_", seq_len(ncol(train_y)))

  log_msg("Processing Scaden training data ...")
  log_msg("Input: %d samples, %d genes, %d cell types",
          nrow(train_x), ncol(train_x), ncol(train_y))

  # -- Process bulk (test) data -----------------------------------------------

  genes_keep <- colnames(train_x)
  test_x     <- NULL

  if (!is.null(bulk_data)) {
    bulk_x <- orient_to_samples_x_genes(as.matrix(bulk_data), genes_keep, "bulk_data")

    if (is.null(colnames(bulk_x)))
      stop("`bulk_data` must have gene names in row or column names.")
    if (is.null(rownames(bulk_x)))
      rownames(bulk_x) <- paste0("sample_", seq_len(nrow(bulk_x)))

    log_msg("Bulk data: %d samples, %d genes", nrow(bulk_x), ncol(bulk_x))

    # -- Gene filtering (variance on bulk) ------------------------------------

    if (!is.null(var_cutoff) || !is.null(top_var_genes)) {
      gene_var <- setNames(matrixStats::colVars(bulk_x), colnames(bulk_x))

      if (!is.null(var_cutoff)) {
        before <- length(gene_var)
        gene_var <- gene_var[gene_var > var_cutoff]
        log_msg("Variance cutoff > %g: %d -> %d genes", var_cutoff, before, length(gene_var))
      }

      if (!is.null(top_var_genes) && length(gene_var) > 0L) {
        top_n <- min(as.integer(top_var_genes), length(gene_var))
        gene_var <- sort(gene_var, decreasing = TRUE)[seq_len(top_n)]
        log_msg("Top variable genes: kept %d", top_n)
      }

      bulk_genes <- names(gene_var)
    } else {
      bulk_genes <- colnames(bulk_x)
    }

    genes_keep <- intersect(genes_keep, bulk_genes)
    log_msg("Shared genes after filtering: %d", length(genes_keep))

    if (length(genes_keep) == 0L)
      stop("No genes remain after intersection — check that training and bulk gene names match.")

    test_x <- normalize(bulk_x[, genes_keep, drop = FALSE], normalize_test)
    log_msg("Normalised test data (%s): %d samples x %d genes",
            normalize_test, nrow(test_x), ncol(test_x))
  }

  # -- Subset and normalise training data -------------------------------------

  train_x <- normalize(train_x[, genes_keep, drop = FALSE], normalize_train)
  log_msg("Normalised training data (%s): %d samples x %d genes",
          normalize_train, nrow(train_x), ncol(train_x))

  log_msg("Done")

  # -- Return -----------------------------------------------------------------

  list(
    train_x         = train_x,
    train_y         = train_y,
    test_x          = test_x,
    genes           = genes_keep,
    train_samples   = rownames(train_x),
    test_samples    = if (!is.null(test_x)) rownames(test_x),
    celltypes       = colnames(train_y),
    normalize_train = normalize_train,
    normalize_test  = normalize_test,
    elapsed         = fmt_time()
  )
}






























#' Set reproducible random seeds for Scaden training
#'
#' Internal helper that sets both the R random seed and the torch random seed.
#'
#' @param seed Integer scalar random seed.
#'
#' @return Invisibly returns `seed`.
#'
#' @importFrom torch torch_manual_seed
#' @keywords internal
#' @noRd
reproducibility <- function(seed = 9) {
  set.seed(seed)
  torch_manual_seed(seed)
  invisible(seed)
}

#' Train a Scaden ensemble in torch
#'
#' Trains a three-model multilayer perceptron ensemble corresponding to the
#' Scaden architectures `m256`, `m512`, and `m1024`. The function expects a
#' processed training matrix with samples in rows and genes in columns, and a
#' matching target matrix with cell-type proportions.
#'
#' @param train_x Numeric matrix with samples in rows and genes in columns.
#' @param train_y Numeric matrix with samples in rows and cell types in columns.
#' @param lr Numeric scalar. Learning rate passed to Adam.
#' @param batch_size Integer scalar. Minibatch size.
#' @param epochs Integer scalar. Number of training epochs for each ensemble
#'   model.
#' @param seed Integer scalar random seed.
#' @param device Optional torch device. If `NULL`, CUDA is used when available,
#'   otherwise CPU is used.
#'
#' @return A named list with components:
#' \describe{
#'   \item{model256}{Trained `m256` torch model.}
#'   \item{model512}{Trained `m512` torch model.}
#'   \item{model1024}{Trained `m1024` torch model.}
#'   \item{loss256}{Numeric vector of batch-level losses for the `m256` model.}
#'   \item{loss512}{Numeric vector of batch-level losses for the `m512` model.}
#'   \item{loss1024}{Numeric vector of batch-level losses for the `m1024`
#'   model.}
#'   \item{epoch_loss256}{Numeric vector of mean epoch losses for the `m256`
#'   model.}
#'   \item{epoch_loss512}{Numeric vector of mean epoch losses for the `m512`
#'   model.}
#'   \item{epoch_loss1024}{Numeric vector of mean epoch losses for the `m1024`
#'   model.}
#'   \item{total_training_time_sec}{Numeric scalar giving total ensemble
#'   training time in seconds.}
#'   \item{architectures}{List of model architectures used.}
#'   \item{lr}{Learning rate used for training.}
#'   \item{batch_size}{Batch size used for training.}
#'   \item{epochs}{Number of training epochs.}
#'   \item{inputdim}{Number of input genes.}
#'   \item{outputdim}{Number of output cell types.}
#'   \item{celltypes}{Character vector of output cell-type names.}
#'   \item{train_samples}{Character vector of training sample names.}
#'   \item{device}{Torch device used for training.}
#' }
#'
#' @details
#' The ensemble consists of three feed-forward neural networks with softmax
#' output. Predictions are averaged across the three trained models in
#' [scaden_predict()]. This function implements the Scaden framework described by Menden et al.
#' (2020). The present implementation follows the PyTorch realization used in
#' the TAPE project, specifically the Scaden-related model definition provided
#' in `TAPE/model.py`, while preserving the three-model ensemble design and
#' prediction strategy of Scaden.
#' 
#' 
#' @references
#' Menden, K., Marouf, M., Oller, S., Dalmia, A., Magruder, D. S., Kloiber, K.,
#' Heutink, P., & Bonn, S. (2020). Deep learning-based cell composition analysis
#' from tissue expression profiles. \emph{Science Advances}, 6(30), eaba2619.
#'
#' TAPE project PyTorch implementation:
#' \url{https://github.com/poseidonchan/TAPE/blob/main/TAPE/model.py}
#'
#' @importFrom torch torch_manual_seed cuda_is_available torch_device
#' @importFrom torch nn_sequential nn_linear nn_dropout nn_relu nn_softmax
#' @importFrom torch nn_init_xavier_uniform_ torch_tensor torch_float
#' @importFrom torch dataloader tensor_dataset optim_adam nnf_l1_loss as_array
#' @importFrom coro loop
#'
#' @examples
#' \dontrun{
#' fit <- scaden(
#'   train_x = proc$train_x,
#'   train_y = proc$train_y,
#'   lr = 1e-4,
#'   batch_size = 128,
#'   epochs = 20
#' )
#' }
#'
#' @export
scaden <- function(train_x,
                   train_y,
                   lr = 1e-4,
                   batch_size = 128,
                   epochs = 20,
                   seed = 123,
                   device = NULL) {
  
  scaden_ascii <- r"(Running
     ____                _            
    / ___|  ___ __ _  __| | ___ _ __  
    \___ \ / __/ _` |/ _` |/ _ \ '_ \ 
     ___) | (_| (_| | (_| |  __/ | | |
    |____/ \___\__,_|\__,_|\___|_| |_|
  )"

  cat("\033[35m", "\033[1m", scaden_ascii, "\033[0m", sep = "")

  reproducibility(seed)

  architectures <- list(
    m256 = list(c(256, 128, 64, 32), c(0, 0, 0, 0)),
    m512 = list(c(512, 256, 128, 64), c(0, 0.3, 0.2, 0.1)),
    m1024 = list(c(1024, 512, 256, 128), c(0, 0.6, 0.3, 0.1))
  )

  if (is.null(device)) {
    device <- if (cuda_is_available()) torch_device("cuda") else torch_device("cpu")
  }

  train_x <- as.matrix(train_x)
  train_y <- as.matrix(train_y)

  inputdim <- ncol(train_x)
  outputdim <- ncol(train_y)

  if (is.null(rownames(train_x))) {
    rownames(train_x) <- paste0("sample_", seq_len(nrow(train_x)))
  }
  if (is.null(colnames(train_y))) {
    colnames(train_y) <- paste0("celltype_", seq_len(ncol(train_y)))
  }
  if (is.null(rownames(train_y))) {
    rownames(train_y) <- rownames(train_x)
  }

  make_mlp <- function(input_dim, output_dim, hidden_units, dropout_rates) {
    nn_sequential(
      nn_linear(input_dim, hidden_units[1]),
      nn_dropout(dropout_rates[1]),
      nn_relu(),
      nn_linear(hidden_units[1], hidden_units[2]),
      nn_dropout(dropout_rates[2]),
      nn_relu(),
      nn_linear(hidden_units[2], hidden_units[3]),
      nn_dropout(dropout_rates[3]),
      nn_relu(),
      nn_linear(hidden_units[3], hidden_units[4]),
      nn_dropout(dropout_rates[4]),
      nn_relu(),
      nn_linear(hidden_units[4], output_dim),
      nn_softmax(dim = 2)
    )
  }

  initialize_weight <- function(m) {
    if (inherits(m, "nn_linear")) {
      with_no_grad({
        nn_init_xavier_uniform_(m$weight)
        if (!is.null(m$bias)) {
          m$bias$fill_(0)
        }
      })
    }
  }

  subtrain <- function(model, optimizer, train_loader, epochs, device, model_name = "model") {
    model$train()
    loss <- c()
    epoch_loss <- numeric(epochs)

    format_secs <- function(x) {
      x <- max(0, as.integer(round(x)))
      h <- x %/% 3600
      m <- (x %% 3600) %/% 60
      s <- x %% 60
      sprintf("%02d:%02d:%02d", h, m, s)
    }

    make_bar <- function(i, n, width = 30L) {
      filled <- if (n > 0) floor(width * i / n) else 0L
      paste0(
        paste(rep("=", filled), collapse = ""),
        paste(rep("-", width - filled), collapse = "")
      )
    }

    start_time <- Sys.time()

    for (i in seq_len(epochs)) {
      batch_losses <- c()

      coro::loop(for (b in train_loader) {
        data <- b[[1]]$to(device = device)
        label <- b[[2]]$to(device = device)

        optimizer$zero_grad()
        batch_loss <- nnf_l1_loss(model(data), label)
        batch_loss$backward()
        optimizer$step()

        loss_value <- as.numeric(batch_loss$to(device = torch_device("cpu"))$item())
        loss <- c(loss, loss_value)
        batch_losses <- c(batch_losses, loss_value)
      })

      epoch_loss[i] <- mean(batch_losses)

      elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
      avg_per_epoch <- elapsed / i
      eta <- avg_per_epoch * (epochs - i)

      cat(sprintf(
        "\r%s | [%s] %3d/%3d | mean loss %.6f | elapsed %s | ETA %s",
        model_name,
        make_bar(i, epochs, width = 30L),
        i, epochs,
        epoch_loss[i],
        format_secs(elapsed),
        format_secs(eta)
      ))
      flush.console()
    }

    cat("\n")

    list(model = model, loss = loss, epoch_loss = epoch_loss)
  }

  x_tensor <- torch_tensor(train_x, dtype = torch_float())
  y_tensor <- torch_tensor(train_y, dtype = torch_float())

  train_loader <- dataloader(
    dataset = tensor_dataset(x_tensor, y_tensor),
    batch_size = batch_size,
    shuffle = TRUE
  )

  model256 <- make_mlp(
    input_dim = inputdim,
    output_dim = outputdim,
    hidden_units = architectures[["m256"]][[1]],
    dropout_rates = architectures[["m256"]][[2]]
  )$to(device = device)

  model512 <- make_mlp(
    input_dim = inputdim,
    output_dim = outputdim,
    hidden_units = architectures[["m512"]][[1]],
    dropout_rates = architectures[["m512"]][[2]]
  )$to(device = device)

  model1024 <- make_mlp(
    input_dim = inputdim,
    output_dim = outputdim,
    hidden_units = architectures[["m1024"]][[1]],
    dropout_rates = architectures[["m1024"]][[2]]
  )$to(device = device)

  model256$apply(initialize_weight)
  model512$apply(initialize_weight)
  model1024$apply(initialize_weight)

  total_start <- Sys.time()

  optimizer256 <- optim_adam(model256$parameters, lr = lr, eps = 1e-07)
  cat("train model256 now\n")
  fit256 <- subtrain(model256, optimizer256, train_loader, epochs, device, "model256")

  optimizer512 <- optim_adam(model512$parameters, lr = lr, eps = 1e-07)
  cat("train model512 now\n")
  fit512 <- subtrain(model512, optimizer512, train_loader, epochs, device, "model512")

  optimizer1024 <- optim_adam(model1024$parameters, lr = lr, eps = 1e-07)
  cat("train model1024 now\n")
  fit1024 <- subtrain(model1024, optimizer1024, train_loader, epochs, device, "model1024")

  total_time <- as.numeric(difftime(Sys.time(), total_start, units = "secs"))
  cat(sprintf("Training of Scaden is done in %.1f seconds\n", total_time))

  list(
    model256 = fit256$model,
    model512 = fit512$model,
    model1024 = fit1024$model,
    loss256 = fit256$loss,
    loss512 = fit512$loss,
    loss1024 = fit1024$loss,
    epoch_loss256 = fit256$epoch_loss,
    epoch_loss512 = fit512$epoch_loss,
    epoch_loss1024 = fit1024$epoch_loss,
    total_training_time_sec = total_time,
    architectures = architectures,
    lr = lr,
    batch_size = batch_size,
    epochs = epochs,
    inputdim = inputdim,
    outputdim = outputdim,
    celltypes = colnames(train_y),
    train_samples = rownames(train_x),
    device = device
  )
}
































#' Predict cell-type proportions with a trained Scaden ensemble
#'
#' Uses a fitted Scaden ensemble returned by [scaden()] to predict cell-type
#' proportions for new bulk samples.
#'
#' @param fit A fitted Scaden object returned by [scaden()].
#' @param test_x Numeric matrix with samples in rows and genes in columns.
#' @param device Optional torch device. If `NULL`, the device stored in `fit` is
#'   used.
#'
#' @return A named list with components:
#' \describe{
#'   \item{average_output}{Average prediction across the three ensemble models.}
#'   \item{prediction_model256}{Predictions from the `m256` model.}
#'   \item{prediction_model512}{Predictions from the `m512` model.}
#'   \item{prediction_model1024}{Predictions from the `m1024` model.}
#'   \item{model256}{The fitted `m256` model.}
#'   \item{model512}{The fitted `m512` model.}
#'   \item{model1024}{The fitted `m1024` model.}
#' }
#'
#' @details
#' Predictions are generated independently for all three ensemble members and
#' then averaged sample-wise.
#' 
#' 
#' @references
#' Menden, K., Marouf, M., Oller, S., Dalmia, A., Magruder, D. S., Kloiber, K.,
#' Heutink, P., & Bonn, S. (2020). Deep learning-based cell composition analysis
#' from tissue expression profiles. \emph{Science Advances}, 6(30), eaba2619.
#'
#' TAPE project PyTorch implementation:
#' \url{https://github.com/poseidonchan/TAPE/blob/main/TAPE/model.py}
#'
#' @importFrom torch torch_tensor torch_float with_no_grad as_array torch_device
#'
#' @examples
#' \dontrun{
#' pred <- scaden_predict(
#'   fit = fit,
#'   test_x = proc$test_x
#' )
#' }
#'
#' @export
scaden_predict <- function(fit, test_x, device = NULL) {
  if (is.null(device)) {
    device <- fit$device
  }

  test_x <- as.matrix(test_x)

  if (is.null(rownames(test_x))) {
    rownames(test_x) <- paste0("sample_", seq_len(nrow(test_x)))
  }

  x_tensor <- torch_tensor(test_x, dtype = torch_float())$to(device = device)

  fit$model256$eval()
  fit$model512$eval()
  fit$model1024$eval()

  pred256 <- with_no_grad({
    fit$model256(x_tensor)
  })

  pred512 <- with_no_grad({
    fit$model512(x_tensor)
  })

  pred1024 <- with_no_grad({
    fit$model1024(x_tensor)
  })

  average_output <- (pred256 + pred512 + pred1024) / 3

  pred256 <- as.matrix(as_array(pred256$to(device = torch_device("cpu"))))
  pred512 <- as.matrix(as_array(pred512$to(device = torch_device("cpu"))))
  pred1024 <- as.matrix(as_array(pred1024$to(device = torch_device("cpu"))))
  average_output <- as.matrix(as_array(average_output$to(device = torch_device("cpu"))))

  rownames(pred256) <- rownames(test_x)
  rownames(pred512) <- rownames(test_x)
  rownames(pred1024) <- rownames(test_x)
  rownames(average_output) <- rownames(test_x)

  colnames(pred256) <- fit$celltypes
  colnames(pred512) <- fit$celltypes
  colnames(pred1024) <- fit$celltypes
  colnames(average_output) <- fit$celltypes

  list(
    average_output = average_output,
    prediction_model256 = pred256,
    prediction_model512 = pred512,
    prediction_model1024 = pred1024,
    model256 = fit$model256,
    model512 = fit$model512,
    model1024 = fit$model1024
  )
}