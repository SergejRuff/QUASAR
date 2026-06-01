
#' Format elapsed seconds as HH:MM:SS
#'
#' Internal helper that converts a duration in seconds to a fixed-width
#' `HH:MM:SS` character string.
#'
#' @param seconds Numeric scalar. Duration in seconds.
#'
#' @return Character scalar giving the formatted duration.
#'
#' @keywords internal
#' @noRd
ot_format_hms <- function(seconds) {
  seconds <- max(0L, as.integer(round(seconds)))
  sprintf("%02d:%02d:%02d", seconds %/% 3600L,
          (seconds %% 3600L) %/% 60L, seconds %% 60L)
}

#' Print an OmicsTweezer progress message
#'
#' Internal helper that prints formatted progress messages. If `start_time` is
#' supplied, the message is prefixed with elapsed runtime in `HH:MM:SS` format.
#'
#' @param msg Character scalar. A format string passed to [sprintf()].
#' @param ... Additional values passed to [sprintf()].
#' @param start_time Optional start time, usually created with [Sys.time()].
#'
#' @return Invisibly returns `NULL`.
#'
#' @keywords internal
#' @noRd
ot_log <- function(msg, ..., start_time = NULL) {
  prefix <- ""
  if (!is.null(start_time)) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    prefix <- sprintf("[%s] ", ot_format_hms(elapsed))
  }
  message(prefix, sprintf(msg, ...))
}

#' Create a text progress bar
#'
#' Internal helper that creates a fixed-width textual progress bar.
#'
#' @param i Integer scalar. Current progress value.
#' @param n Integer scalar. Total number of steps.
#' @param width Integer scalar. Width of the progress bar.
#'
#' @return Character scalar containing the progress bar.
#'
#' @keywords internal
#' @noRd
ot_make_bar <- function(i, n, width = 30L) {
  if (n <= 0L) return(paste(rep("-", width), collapse = ""))
  filled <- max(0L, min(width, as.integer(floor(width * i / n))))
  paste0(paste(rep("=", filled), collapse = ""),
         paste(rep("-", width - filled), collapse = ""))
}

#' Resolve torch device for OmicsTweezer
#'
#' Internal helper that resolves CPU or CUDA execution for OmicsTweezer
#' training and prediction.
#'
#' @param device Character scalar. One of `"auto"`, `"cpu"`, or `"cuda"`.
#'   If `"auto"`, CUDA is used when available, otherwise CPU is used.
#' @param cuda_index Optional integer scalar giving the CUDA device index to use.
#'   CUDA indices are zero-based.
#' @param verbose Logical scalar. If `TRUE`, prints the selected device.
#' @param stage Optional character scalar describing the workflow stage for the
#'   progress message.
#'
#' @return A named list with components:
#' \describe{
#'   \item{device}{Resolved torch device object.}
#'   \item{device_name}{Character scalar describing the selected device.}
#'   \item{using_cuda}{Logical scalar indicating whether CUDA is used.}
#' }
#'
#' @importFrom torch cuda_is_available cuda_device_count torch_device
#'
#' @keywords internal
#' @noRd
omics_resolve_device <- function(device = c("auto", "cpu", "cuda"),
                                 cuda_index = NULL,
                                 verbose = TRUE,
                                 stage = NULL) {
  device <- match.arg(device)

  cuda_ok <- torch::cuda_is_available()

  if (identical(device, "auto")) {
    device <- if (cuda_ok) "cuda" else "cpu"
  }

  if (identical(device, "cuda")) {
    if (!cuda_ok) {
      stop("device = 'cuda' was requested, but torch::cuda_is_available() is FALSE")
    }

    n_cuda <- torch::cuda_device_count()

    if (!is.null(cuda_index)) {
      cuda_index <- as.integer(cuda_index)

      if (length(cuda_index) != 1L ||
          is.na(cuda_index) ||
          cuda_index < 0L ||
          cuda_index >= n_cuda) {
        stop(sprintf(
          "cuda_index must be a single integer in [0, %d]",
          max(0L, n_cuda - 1L)
        ))
      }

      dev <- torch::torch_device("cuda", cuda_index)
      dev_name <- sprintf("cuda:%d", cuda_index)
    } else {
      dev <- torch::torch_device("cuda")
      dev_name <- "cuda"
    }

    out <- list(
      device = dev,
      device_name = dev_name,
      using_cuda = TRUE
    )
  } else {
    out <- list(
      device = torch::torch_device("cpu"),
      device_name = "cpu",
      using_cuda = FALSE
    )
  }

  if (isTRUE(verbose)) {
    if (is.null(stage)) {
      cat(sprintf("Using device: %s\n", out$device_name))
    } else {
      cat(sprintf("Using device for %s: %s\n", stage, out$device_name))
    }
  }

  out
}

#' Set reproducible random seeds for OmicsTweezer
#'
#' Sets both the R random seed and the torch random seed.
#'
#' @param seed Integer scalar random seed.
#'
#' @return Invisibly returns `seed`.
#'
#' @importFrom torch torch_manual_seed
#'
#' @export
omics_set_seed <- function(seed = 42L) {
  set.seed(seed)
  torch::torch_manual_seed(seed)
  message(sprintf("seed is fixed, seed is %d", seed))
  invisible(seed)
}


#' Extract single-cell counts and cell-type labels
#'
#' Extracts a genes x cells expression matrix and matching cell-type labels from
#' a `Seurat` object, a `SingleCellExperiment` object, a matrix, or a
#' `data.frame`. The returned count matrix is transposed to cells x genes for
#' internal OmicsTweezer processing.
#'
#' @param sc_data A single-cell reference object. Supported inputs are a
#'   `Seurat` object, a `SingleCellExperiment` object, a matrix, or a
#'   `data.frame`.
#' @param celltype_col Character scalar giving the metadata column containing
#'   cell-type labels for `Seurat` and `SingleCellExperiment` input.
#' @param assay Character scalar giving the assay name to extract from a
#'   `Seurat` object.
#' @param slot Character scalar giving the assay slot or layer to extract. For
#'   `Seurat`, this is interpreted as a slot or layer within `assay`. For
#'   `SingleCellExperiment`, this is interpreted as an assay name.
#'
#' @return A named list with components:
#' \describe{
#'   \item{counts}{Numeric matrix with cells in rows and genes in columns.}
#'   \item{celltypes}{Character vector of cell-type labels, one per cell.}
#'   \item{genes}{Character vector of gene names.}
#' }
#'
#' @details
#' Matrix and `data.frame` input are expected to follow the genes x cells
#' convention. For these inputs, column names are interpreted as cell-type
#' labels.
#'
#' @importFrom methods slot slotNames
#' @importFrom SeuratObject LayerData
#' @importFrom SummarizedExperiment assay assayNames colData
#'
#' @keywords internal
#' @noRd
omics_extract_sc_data <- function(sc_data,
                                  celltype_col = "CellType",
                                  assay = "RNA",
                                  slot = "counts") {

  if (inherits(sc_data, "Seurat")) {
    if (!assay %in% names(sc_data@assays))
      stop("Assay '", assay, "' not found in Seurat object.")
    meta <- sc_data@meta.data
    if (!celltype_col %in% colnames(meta))
      stop("Column '", celltype_col, "' not found in Seurat meta.data.")

    mat <- NULL
    if ("LayerData" %in% getNamespaceExports("SeuratObject")) {
      mat <- tryCatch(
        SeuratObject::LayerData(sc_data, assay = assay, layer = slot),
        error = function(e) NULL)
    }
    if (is.null(mat)) {
      assay_obj <- sc_data[[assay]]
      mat <- switch(slot,
                    counts     = assay_obj$counts,
                    data       = assay_obj$data,
                    scale.data = assay_obj$scale.data,
                    if (slot %in% methods::slotNames(assay_obj))
                      methods::slot(assay_obj, slot)
                    else
                      stop("Could not extract slot/layer '", slot,
                           "' from assay '", assay, "'."))
    }
    mat <- as.matrix(mat)            # genes x cells
    ct  <- as.character(meta[[celltype_col]])
    genes <- rownames(mat)

  } else if (inherits(sc_data, "SingleCellExperiment")) {
    avail <- SummarizedExperiment::assayNames(sc_data)
    if (!slot %in% avail)
      stop("Assay '", slot, "' not found in SingleCellExperiment object.")
    meta <- as.data.frame(SummarizedExperiment::colData(sc_data))
    if (!celltype_col %in% colnames(meta))
      stop("Column '", celltype_col, "' not found in colData.")
    mat   <- as.matrix(SummarizedExperiment::assay(sc_data, slot))   # genes x cells
    ct    <- as.character(meta[[celltype_col]])
    genes <- rownames(mat)

  } else if (is.matrix(sc_data) || is.data.frame(sc_data)) {
    mat <- as.matrix(sc_data)       
    ct  <- colnames(mat)             
    if (is.null(ct))
      stop("Matrix input must have column names = cell-type labels ",
           "(genes x cells convention).")
    genes <- rownames(mat)

  } else {
    stop("sc_data must be a Seurat object, a SingleCellExperiment object, ",
         "or a matrix/data.frame (genes x cells).")
  }

  if (is.null(genes)) {
    genes <- paste0("gene_", seq_len(nrow(mat)))
    rownames(mat) <- genes
  }


  counts <- t(mat)
  list(counts = counts, celltypes = ct, genes = genes)
}


#' Simulate pseudobulk training data for OmicsTweezer
#'
#' Generates artificial pseudobulk samples from single-cell reference data and
#' returns them in the genes x samples format used by the OmicsTweezer workflow.
#'
#' Counts are sampled with replacement within each cell type according to
#' Dirichlet-generated mixture proportions. Optional sparse and rare-cell
#' perturbations can be applied to mimic heterogeneous cellular compositions.
#'
#' @param sc_data A single-cell reference object. Supported inputs are a
#'   `Seurat` object, a `SingleCellExperiment` object, a matrix, or a
#'   `data.frame`.
#' @param d_prior Optional numeric vector of Dirichlet concentration parameters.
#'   If `NULL`, a symmetric Dirichlet prior of ones is used.
#' @param n Integer scalar. Target number of single cells per simulated
#'   pseudobulk before rounding to integer cell counts.
#' @param samplenum Integer scalar. Number of pseudobulk samples to generate.
#' @param random_state Optional integer scalar used to control reproducible
#'   simulation.
#' @param sparse Logical scalar. If `TRUE`, a subset of simulated samples is
#'   forced to contain zero fractions for a subset of cell types.
#' @param sparse_prob Numeric scalar. Controls both the proportion of sparse
#'   samples and the proportion of cell types zeroed within those samples.
#' @param rare Logical scalar. If `TRUE`, a subset of cell types is perturbed to
#'   have very small fractions in a subset of samples.
#' @param rare_percentage Numeric scalar in `[0, 1]`. Fraction of cell types to
#'   treat as rare when `rare = TRUE`.
#' @param celltype_col Character scalar giving the metadata column that contains
#'   cell-type labels for `Seurat` and `SingleCellExperiment` input.
#' @param assay Character scalar giving the assay name to extract from a
#'   `Seurat` object.
#' @param slot Character scalar giving the assay slot or assay name to extract.
#'   For `Seurat`, this is interpreted as a slot or layer within `assay`. For
#'   `SingleCellExperiment`, this is interpreted as an assay name.
#'
#' @return A named list with components:
#' \describe{
#'   \item{X}{Numeric matrix with genes in rows and simulated pseudobulk samples
#'   in columns.}
#'   \item{obs}{Data frame with samples in rows and realized cell-type
#'   proportions in columns.}
#'   \item{var}{Data frame indexed by gene names.}
#' }
#'
#' @details
#' The simulation follows the OmicsTweezer-style pseudobulk generation strategy.
#' Single-cell counts are first grouped by cell type, then sampled according to
#' Dirichlet-generated proportions and aggregated into pseudobulk profiles.
#'
#'
#' @importFrom MCMCpack rdirichlet
#' @importFrom Matrix sparseMatrix
#' @importFrom stats runif
#' @importFrom utils txtProgressBar setTxtProgressBar
#'
#' @export
omics_simulate <- function(sc_data,
                           d_prior = NULL,
                           n = 500L,
                           samplenum = 5000L,
                           random_state = NULL,
                           sparse = TRUE,
                           sparse_prob = 0.5,
                           rare = FALSE,
                           rare_percentage = 0.4,
                           celltype_col = "CellType",
                           assay = "RNA",
                           slot = "counts") {

  sim_start <- Sys.time()
  ot_log("Simulation started", start_time = sim_start)

  if (!is.null(random_state)) set.seed(random_state)

  ref       <- omics_extract_sc_data(sc_data, celltype_col, assay, slot)
  sc_counts <- ref$counts                       # cells x genes
  celltypes <- ref$celltypes
  genename  <- ref$genes

  # match Python: groupby('celltype') preserves first-seen order
  ct_levels    <- unique(celltypes)
  num_celltype <- length(ct_levels)
  if (is.null(d_prior)) d_prior <- rep(1, num_celltype)
  if (length(d_prior) != num_celltype)
    stop("Length of `d_prior` must equal the number of cell types ",
         "(found ", num_celltype, ").")

  ct_groups <- split(seq_len(nrow(sc_counts)),
                     factor(celltypes, levels = ct_levels))

  ot_log("Drawing Dirichlet proportions", start_time = sim_start)

  prop <- MCMCpack::rdirichlet(samplenum, d_prior)
  prop <- prop / rowSums(prop)

  # ---- sparse perturbation (matches Python exactly) -----------------------
  if (sparse) {
    ot_log("Applying sparse perturbation (prob = %.2f)", sparse_prob,
           start_time = sim_start)
    n_sparse_samples <- as.integer(nrow(prop) * sparse_prob)
    n_zero_ct        <- as.integer(ncol(prop) * sparse_prob)
    if (n_sparse_samples > 0L && n_zero_ct > 0L) {
      n_zero_ct <- min(n_zero_ct, ncol(prop) - 1L)
      for (i in seq_len(n_sparse_samples)) {
        idx_zero <- sample.int(ncol(prop), n_zero_ct, replace = FALSE)
        prop[i, idx_zero] <- 0
      }
      prop <- prop / rowSums(prop)
    }
  }

  # ---- rare perturbation (matches Python's quirky row count exactly) ------
  if (rare) {
    ot_log("Applying rare-cell perturbation", start_time = sim_start)
    set.seed(0)  # original hard-codes np.random.seed(0) here
    n_rare_ct <- as.integer(ncol(prop) * rare_percentage)
    if (n_rare_ct > 0L) {
      rare_idx  <- sample.int(ncol(prop), n_rare_ct, replace = FALSE)
      prop      <- prop / rowSums(prop)
      n_rows_rare <- as.integer(0.5 * nrow(prop)) +
        as.integer(as.integer(rare_percentage * 0.5 * nrow(prop)))
      n_rows_rare <- min(n_rows_rare, nrow(prop))
      for (i in seq_len(n_rows_rare)) {
        buf <- stats::runif(length(rare_idx), 0, 0.03)
        prop[i, rare_idx] <- 0
        remain <- sum(prop[i, ])
        if (remain > 0) {
          prop[i, ] <- (1 - sum(buf)) * prop[i, ] / remain
        }
        prop[i, rare_idx] <- buf
      }
    }
  }


  cell_num <- floor(n * prop)
  zero_rows <- which(rowSums(cell_num) == 0)
  if (length(zero_rows) > 0L) {
    best_ct <- max.col(prop[zero_rows, , drop = FALSE], ties.method = "first")
    for (i in seq_along(zero_rows)) cell_num[zero_rows[i], best_ct[i]] <- 1L
  }
  prop <- cell_num / rowSums(cell_num)


  ot_log("Sampling cells and building pseudobulks", start_time = sim_start)
  sample_mat <- matrix(0, nrow = samplenum, ncol = ncol(sc_counts))
  pb <- utils::txtProgressBar(min = 0, max = length(ct_levels), style = 3)
  for (j in seq_along(ct_levels)) {
    idx_pool <- ct_groups[[j]]
    counts_j <- as.integer(cell_num[, j])
    if (length(idx_pool) == 0L || sum(counts_j) == 0L) {
      utils::setTxtProgressBar(pb, j); next
    }
    sample_ids <- rep.int(seq_len(samplenum), times = counts_j)
    local_ids  <- sample.int(length(idx_pool), size = length(sample_ids),
                             replace = TRUE)
    if (!requireNamespace("Matrix", quietly = TRUE))
      stop("Package 'Matrix' is required.")
    W <- Matrix::sparseMatrix(i = local_ids, j = sample_ids, x = 1,
                              dims = c(length(idx_pool), samplenum))
    contrib    <- Matrix::t(Matrix::t(sc_counts[idx_pool, , drop = FALSE]) %*% W)
    sample_mat <- sample_mat + as.matrix(contrib)
    utils::setTxtProgressBar(pb, j)
  }
  close(pb)

  sample_names <- paste0("sample_", seq_len(samplenum))
  rownames(sample_mat) <- sample_names
  colnames(sample_mat) <- genename
  prop_df <- as.data.frame(prop)
  colnames(prop_df) <- ct_levels
  rownames(prop_df) <- sample_names

  ot_log("Simulation finished", start_time = sim_start)
  list(X   = t(sample_mat),                     # genes x samples (per user spec)
       obs = prop_df,                           # samples x celltypes
       var = data.frame(row.names = genename))
}


#' Process simulated and bulk data for OmicsTweezer
#'
#' Aligns simulated pseudobulk training data with real bulk expression data,
#' filters genes by variance, applies log transformation and per-sample scaling,
#' and returns matrices ready for OmicsTweezer training and prediction.
#'
#' @param simudata A list returned by [omics_simulate()]. Must contain `X` and
#'   `obs`.
#' @param real_bulk Numeric matrix or `data.frame` with genes in rows and bulk
#'   samples in columns.
#' @param variance_threshold Numeric scalar. Fraction of genes used to define
#'   the variance cutoff separately in simulated and real bulk data.
#' @param scaler Character scalar specifying the per-sample scaling method. One
#'   of `"ss"` for standard scaling, `"mms"` for min-max scaling, or `"none"`.
#'
#' @return A named list with components:
#' \describe{
#'   \item{train_x}{Processed training matrix with samples in rows and genes in
#'   columns.}
#'   \item{train_y}{Training target matrix with samples in rows and cell types
#'   in columns.}
#'   \item{test_x}{Processed real bulk matrix with samples in rows and genes in
#'   columns.}
#'   \item{genename}{Character vector of retained shared genes.}
#'   \item{celltypes}{Character vector of cell-type names.}
#'   \item{samplename}{Character vector of real bulk sample names.}
#' }
#'
#' @details
#' Both simulated and real bulk matrices are filtered by variance before being
#' restricted to their shared gene set. The log1p transformation is applied
#' before optional per-sample scaling.
#'
#' @importFrom stats var sd
#'
#' @export
omics_process <- function(simudata, real_bulk,
                          variance_threshold = 0.98,
                          scaler = c("ss", "mms", "none")) {
  scaler <- match.arg(scaler)
  proc_start <- Sys.time()
  ot_log("Processing started (scaler = %s)", scaler, start_time = proc_start)


  train_x <- as.data.frame(t(simudata$X))                        
  train_y <- as.matrix(simudata$obs)                             
  test_x  <- as.data.frame(t(as.matrix(real_bulk)))             

  ot_log("Variance filtering", start_time = proc_start)
  train_var    <- apply(train_x, 2, var)
  cutoff_train <- sort(train_var, decreasing = TRUE)[
    as.integer(ncol(train_x) * variance_threshold)]
  train_x <- train_x[, train_var > cutoff_train, drop = FALSE]

  test_var    <- apply(test_x, 2, var)
  cutoff_test <- sort(test_var, decreasing = TRUE)[
    as.integer(ncol(test_x) * variance_threshold)]
  test_x <- test_x[, test_var > cutoff_test, drop = FALSE]

  inter   <- intersect(colnames(train_x), colnames(test_x))
  train_x <- train_x[, inter, drop = FALSE]
  test_x  <- test_x[,  inter, drop = FALSE]
  ot_log("Intersected gene number: %d", length(inter), start_time = proc_start)

  genename   <- inter
  celltypes  <- colnames(train_y)
  samplename <- rownames(test_x)

  ot_log("log1p transform", start_time = proc_start)
  train_x <- log1p(as.matrix(train_x))
  test_x  <- log1p(as.matrix(test_x))


  ot_log("Scaling (%s, per sample)", scaler, start_time = proc_start)
  if (scaler == "ss") {
    train_mean <- rowMeans(train_x)
    train_sd   <- apply(train_x, 1, sd)
    sd_safe    <- ifelse(train_sd == 0, 1, train_sd)
    train_x    <- sweep(sweep(train_x, 1, train_mean, "-"), 1, sd_safe, "/")
    train_x[train_sd == 0, ] <- 0

    test_mean <- rowMeans(test_x)
    test_sd   <- apply(test_x, 1, sd)
    sd_safe   <- ifelse(test_sd == 0, 1, test_sd)
    test_x    <- sweep(sweep(test_x, 1, test_mean, "-"), 1, sd_safe, "/")
    test_x[test_sd == 0, ] <- 0

  } else if (scaler == "mms") {
    train_min <- apply(train_x, 1, min)
    train_max <- apply(train_x, 1, max)
    rng       <- train_max - train_min
    rng_safe  <- ifelse(rng == 0, 1, rng)
    train_x   <- sweep(sweep(train_x, 1, train_min, "-"), 1, rng_safe, "/")
    train_x[rng == 0, ] <- 0

    test_min <- apply(test_x, 1, min)
    test_max <- apply(test_x, 1, max)
    rng      <- test_max - test_min
    rng_safe <- ifelse(rng == 0, 1, rng)
    test_x   <- sweep(sweep(test_x, 1, test_min, "-"), 1, rng_safe, "/")
    test_x[rng == 0, ] <- 0
  }


  ot_log("Processing finished", start_time = proc_start)
  list(train_x = train_x, train_y = train_y, test_x = test_x,
       genename = genename, celltypes = celltypes, samplename = samplename)
}


#' Create an OmicsTweezer encoder block
#'
#' Internal helper that creates one fully connected encoder block consisting of
#' a linear layer, leaky ReLU activation, and dropout.
#'
#' @param in_dim Integer scalar. Number of input features.
#' @param out_dim Integer scalar. Number of output features.
#' @param do_rate Numeric scalar. Dropout probability.
#'
#' @return A torch `nn_sequential` module.
#'
#' @importFrom torch nn_sequential nn_linear nn_leaky_relu nn_dropout
#'
#' @keywords internal
#' @noRd
omics_encoder_block <- function(in_dim, out_dim, do_rate) {
  torch::nn_sequential(
    torch::nn_linear(in_dim, out_dim),
    torch::nn_leaky_relu(negative_slope = 0.2),
    torch::nn_dropout(p = do_rate)
  )
}

#' Create an OmicsTweezer model
#'
#' Builds the encoder and predictor modules for one OmicsTweezer architecture.
#' The encoder contains two fully connected blocks, and the predictor contains
#' two fully connected blocks followed by a linear output layer and softmax.
#'
#' @param feature_num Integer scalar. Number of input genes.
#' @param celltype_num Integer scalar. Number of output cell types.
#' @param dims Integer vector of length 4 specifying hidden-layer dimensions.
#' @param drops Numeric vector of length 4 specifying dropout probabilities.
#'
#' @return A named list with components:
#' \describe{
#'   \item{encoder}{Torch encoder module.}
#'   \item{predictor}{Torch predictor module returning cell-type proportions.}
#' }
#'
#' @details
#' In R torch, softmax over the feature axis of an `N x F` tensor is specified
#' with `dim = 2`.
#'
#' @importFrom torch nn_sequential nn_linear nn_softmax
#'
#' @export
omics_create_model <- function(feature_num, celltype_num, dims, drops) {
  if (length(dims) != 4L || length(drops) != 4L)
    stop("`dims` and `drops` must each have length 4 ",
         "(matches m256/m512/m1024 architectures).")

  encoder <- torch::nn_sequential(
    omics_encoder_block(feature_num, dims[1], drops[1]),
    omics_encoder_block(dims[1],     dims[2], drops[2])
  )
  predictor <- torch::nn_sequential(
    omics_encoder_block(dims[2], dims[3], drops[3]),
    omics_encoder_block(dims[3], dims[4], drops[4]),
    torch::nn_linear(dims[4], celltype_num),
    torch::nn_softmax(dim = 2)
  )
  list(encoder = encoder, predictor = predictor)
}


#' Compute mean squared prediction loss
#'
#' Internal helper that computes the mean squared error between predicted and
#' target cell-type proportions.
#'
#' @param preds Torch tensor containing predicted cell-type proportions.
#' @param gt Torch tensor containing target cell-type proportions.
#'
#' @return A scalar torch tensor containing the mean squared error.
#'
#' @keywords internal
#' @noRd
omics_mse_loss <- function(preds, gt) {
  ((preds - gt)^2)$mean()
}

#' Compute naive Wasserstein-style distance
#'
#' Internal helper that computes the absolute difference between the scalar
#' means of source and target embeddings.
#'
#' @param emb_source Torch tensor containing source-domain embeddings.
#' @param emb_target Torch tensor containing target-domain embeddings.
#'
#' @return A scalar torch tensor containing the absolute mean difference.
#'
#' @details
#' This helper mirrors the simplified Wasserstein-style distance used in the
#' translated OmicsTweezer training loop.
#'
#' @keywords internal
#' @noRd
omics_naive_w_distance <- function(emb_source, emb_target) {
  (emb_source$mean() - emb_target$mean())$abs()
}


#' Train one OmicsTweezer model
#'
#' Trains one OmicsTweezer neural network with supervised cell-type proportion
#' prediction on simulated pseudobulks and a Wasserstein-style domain adaptation
#' term using real bulk samples.
#'
#' @param train_x Numeric matrix with samples in rows and genes in columns.
#' @param train_y Numeric matrix with samples in rows and cell types in columns.
#' @param test_x Numeric matrix with target-domain bulk samples in rows and
#'   genes in columns.
#' @param dims Integer vector of length 4 specifying hidden-layer dimensions.
#' @param drops Numeric vector of length 4 specifying dropout probabilities.
#' @param epochs Integer scalar. Number of training epochs.
#' @param batch_size Integer scalar. Minibatch size.
#' @param device Character scalar. One of `"auto"`, `"cpu"`, or `"cuda"`.
#'   If `"auto"`, CUDA is used when available, otherwise CPU is used.
#' @param cuda_index Optional integer scalar giving the CUDA device index to use
#'   when `device = "cuda"`. CUDA indices are zero-based.
#' @param learning_rate Numeric scalar. Learning rate passed to Adam.
#' @param loss_weight Numeric scalar. Weight applied to the Wasserstein-style
#'   domain adaptation loss.
#' @param seed Integer scalar random seed.
#' @param verbose Logical scalar. If `TRUE`, prints an epoch-level progress bar.
#'
#' @return A named list with components:
#' \describe{
#'   \item{encoder}{Trained torch encoder module.}
#'   \item{predictor}{Trained torch predictor module.}
#'   \item{history}{Data frame with epoch-level total, MSE, OT, and weighted OT
#'   losses.}
#' }
#'
#' @details
#' Each batch contains source-domain simulated pseudobulks with known
#' proportions and target-domain real bulk profiles without labels. The training
#' objective combines supervised MSE loss with a domain adaptation term computed
#' from the encoded source and target embeddings.
#'
#' @references
#' TAPE project PyTorch implementation:
#' \url{https://github.com/poseidonchan/TAPE/blob/main/TAPE/model.py}
#'
#' @importFrom torch cuda_is_available optim_adam torch_tensor torch_float
#' @importFrom utils flush.console
#'
#' @export
omics_train <- function(train_x, train_y, test_x,
                        dims, drops,
                        epochs = 128L, batch_size = 128L,
                        device = c("auto", "cpu", "cuda"),
                        cuda_index = NULL,
                        learning_rate = 1e-4,
                        loss_weight = 1.0,
                        seed = 2021L,
                        verbose = TRUE) {

  train_start <- Sys.time()
  ot_log("Training started (epochs = %d, batch = %d, lr = %.1e)",
         epochs, batch_size, learning_rate, start_time = train_start)

  omics_set_seed(seed)

  feature_num  <- ncol(train_x)
  celltype_num <- ncol(train_y)
  model <- omics_create_model(feature_num, celltype_num, dims, drops)

  dev_info <- omics_resolve_device(
    device = device,
    cuda_index = cuda_index,
    verbose = verbose,
    stage = "training"
  )

  dev <- dev_info$device
  model$encoder   <- model$encoder$to(device = dev)
  model$predictor <- model$predictor$to(device = dev)

  optimizer <- torch::optim_adam(
    c(model$encoder$parameters, model$predictor$parameters),
    lr = learning_rate
  )

  n_src <- nrow(train_x); n_tgt <- nrow(test_x)
  src_x_t <- torch::torch_tensor(train_x, dtype = torch::torch_float(), device = dev)
  src_y_t <- torch::torch_tensor(train_y, dtype = torch::torch_float(), device = dev)
  tgt_x_t <- torch::torch_tensor(test_x,  dtype = torch::torch_float(), device = dev)

  model$encoder$train(); model$predictor$train()
  history <- data.frame(
    epoch = seq_len(epochs),
    total_loss = NA_real_,
    mse_loss = NA_real_,
    ot_loss = NA_real_,
    weighted_ot_loss = NA_real_
  )
  model_name <- sprintf("OmicsTweezer[d=%d]", dims[1])

  for (ep in seq_len(epochs)) {
    src_perm <- sample.int(n_src)
    tgt_perm <- sample.int(n_tgt)
    tgt_cursor <- 1L
    total_loss_epoch <- 0
    mse_loss_epoch <- 0
    ot_loss_epoch <- 0
    weighted_ot_loss_epoch <- 0
    n_batches <- 0L

    for (start in seq(1L, n_src, by = batch_size)) {
      end <- min(start + batch_size - 1L, n_src)
      src_idx <- src_perm[start:end]

      # cycle target iterator (matches the try/StopIteration pattern)
      bs <- length(src_idx)
      if (tgt_cursor + bs - 1L > n_tgt) {
        tgt_perm   <- sample.int(n_tgt)
        tgt_cursor <- 1L
      }
      tgt_idx    <- tgt_perm[tgt_cursor:(tgt_cursor + bs - 1L)]
      tgt_cursor <- tgt_cursor + bs

      src_x <- src_x_t[src_idx, , drop = FALSE]
      src_y <- src_y_t[src_idx, , drop = FALSE]
      tgt_x <- tgt_x_t[tgt_idx, , drop = FALSE]

      emb_src   <- model$encoder(src_x)
      emb_tgt   <- model$encoder(tgt_x)
      frac_pred <- model$predictor(emb_src)

      pred_loss <- omics_mse_loss(frac_pred, src_y)
      w_dist    <- omics_naive_w_distance(emb_src, emb_tgt)
      loss      <- pred_loss + loss_weight * w_dist

      optimizer$zero_grad()
      loss$backward(retain_graph = TRUE)        # kept 1:1 with original
      optimizer$step()

      mse_val <- as.numeric(pred_loss$item())
      ot_val  <- as.numeric(w_dist$item())
      total_val <- as.numeric(loss$item())

      mse_loss_epoch <- mse_loss_epoch + mse_val
      ot_loss_epoch <- ot_loss_epoch + ot_val
      weighted_ot_loss_epoch <- weighted_ot_loss_epoch + loss_weight * ot_val
      total_loss_epoch <- total_loss_epoch + total_val

      n_batches <- n_batches + 1L
    }

    n_batches_safe <- max(1L, n_batches)

    history$total_loss[ep] <- total_loss_epoch / n_batches_safe
    history$mse_loss[ep] <- mse_loss_epoch / n_batches_safe
    history$ot_loss[ep] <- ot_loss_epoch / n_batches_safe
    history$weighted_ot_loss[ep] <- weighted_ot_loss_epoch / n_batches_safe

    if (verbose) {
      elapsed <- as.numeric(difftime(Sys.time(), train_start, units = "secs"))
      eta     <- (elapsed / ep) * (epochs - ep)
      cat(sprintf(
        "\r%s | [%s] %3d/%3d | total %.6f | mse %.6f | ot %.6f | elapsed %s | ETA %s",
        model_name, ot_make_bar(ep, epochs), ep, epochs,
        history$total_loss[ep],
        history$mse_loss[ep],
        history$ot_loss[ep],
        ot_format_hms(elapsed),
        ot_format_hms(eta)
      ))
      utils::flush.console()
    }
  }
  if (verbose) cat("\n")
  ot_log("Training finished", start_time = train_start)

  list(
    encoder = model$encoder,
    predictor = model$predictor,
    history = history
  )
}




#' Predict cell-type proportions with one OmicsTweezer model
#'
#' Uses a trained OmicsTweezer encoder and predictor to estimate cell-type
#' proportions for real bulk samples.
#'
#' @param model A trained model returned by [omics_train()]. Must contain
#'   `encoder` and `predictor`.
#' @param test_x Numeric matrix with samples in rows and genes in columns.
#' @param celltypes Character vector of output cell-type names.
#' @param samplename Optional character vector of sample names for the returned
#'   prediction table.
#' @param batch_size Integer scalar. Prediction batch size.
#' @param device Character scalar. One of `"auto"`, `"cpu"`, or `"cuda"`.
#'   If `"auto"`, CUDA is used when available, otherwise CPU is used.
#' @param cuda_index Optional integer scalar giving the CUDA device index to use
#'   when `device = "cuda"`. CUDA indices are zero-based.
#' @param verbose Logical scalar. If `TRUE`, prints the selected prediction
#'   device.
#'
#' @return A data frame with samples in rows and predicted cell-type proportions
#' in columns.
#'
#' @details
#' Predictions are generated in evaluation mode and returned on the CPU as a
#' regular R `data.frame`.
#'
#' @importFrom torch cuda_is_available torch_tensor torch_float with_no_grad
#' @importFrom torch torch_tensor torch_float with_no_grad
#'
#' @export
omics_predict <- function(model, test_x, celltypes, samplename = NULL,
                          batch_size = 128L,device = c("auto", "cpu", "cuda"),
                          cuda_index = NULL,verbose=TRUE) {

  # dev <- if (torch::cuda_is_available()) "cuda" else "cpu"
  # model$encoder$eval(); model$predictor$eval()
  dev_info <- omics_resolve_device(
    device = device,
    cuda_index = cuda_index,
    verbose = verbose,
    stage = "prediction"
  )

  dev <- dev_info$device

  model$encoder   <- model$encoder$to(device = dev)
  model$predictor <- model$predictor$to(device = dev)

  model$encoder$eval()
  model$predictor$eval()

  n <- nrow(test_x)
  out <- matrix(0, nrow = n, ncol = length(celltypes))
  test_x_t <- torch::torch_tensor(test_x, dtype = torch::torch_float(), device = dev)

  torch::with_no_grad({
    for (start in seq(1L, n, by = batch_size)) {
      end <- min(start + batch_size - 1L, n)
      x   <- test_x_t[start:end, , drop = FALSE]
      logits <- model$predictor(model$encoder(x))
      out[start:end, ] <- as.matrix(logits$cpu())
    }
  })

  pred <- as.data.frame(out)
  colnames(pred) <- celltypes
  if (!is.null(samplename) && length(samplename) == nrow(pred))
    rownames(pred) <- samplename
  pred
}


#' Run the complete OmicsTweezer workflow
#'
#' Runs OmicsTweezer cell-type deconvolution from a single-cell reference and
#' real bulk expression profiles. The workflow simulates pseudobulks, processes
#' simulated and real bulk data, trains three neural network architectures, and
#' averages their predictions.
#'
#' @param sc_data A single-cell reference object. Supported inputs are a
#'   `Seurat` object, a `SingleCellExperiment` object, a matrix, or a
#'   `data.frame`.
#' @param real_bulk Numeric matrix or `data.frame` with genes in rows and real
#'   bulk samples in columns.
#' @param real_bulk_obs Optional matrix or `data.frame` with ground-truth
#'   cell-type proportions for real bulk samples. Used only for reporting in the
#'   returned object.
#' @param ot_weight Numeric scalar. Weight applied to the Wasserstein-style
#'   domain adaptation loss.
#' @param scale_minmax Logical scalar. If `TRUE`, applies an outer global
#'   min-max scaling step to simulated and real bulk count matrices before
#'   processing.
#' @param samplenum Integer scalar. Number of simulated pseudobulk samples.
#' @param n_cells_per_bulk Integer scalar. Target number of single cells per
#'   simulated pseudobulk.
#' @param sparse Logical scalar. If `TRUE`, applies sparse perturbation during
#'   pseudobulk simulation.
#' @param variance_threshold Numeric scalar. Fraction of genes used to define
#'   the variance cutoff during processing.
#' @param scaler Character scalar specifying the per-sample scaling method. One
#'   of `"ss"`, `"mms"`, or `"none"`.
#' @param batch_size Integer scalar. Minibatch size for training and prediction.
#' @param epochs Integer scalar. Number of training epochs for each architecture.
#' @param learning_rate Numeric scalar. Learning rate passed to Adam.
#' @param seed Integer scalar random seed.
#' @param celltype_col Character scalar giving the metadata column containing
#'   cell-type labels for `Seurat` and `SingleCellExperiment` input.
#' @param assay Character scalar giving the assay name to extract from a
#'   `Seurat` object.
#' @param slot Character scalar giving the assay slot or assay name to extract.
#' @param verbose Logical scalar. If `TRUE`, prints workflow progress messages.
#' @param device Character scalar. One of `"auto"`, `"cpu"`, or `"cuda"`.
#'   If `"auto"`, CUDA is used when available, otherwise CPU is used.
#' @param cuda_index Optional integer scalar giving the CUDA device index to use
#'   when `device = "cuda"`. CUDA indices are zero-based.
#'
#' @return A named list with components:
#' \describe{
#'   \item{pred}{Data frame of averaged predicted cell-type proportions. If
#'   `real_bulk_obs` is supplied, this component is itself a list containing
#'   `pred` and `ground_truth`.}
#'   \item{per_model}{Named list of predictions from the `m256`, `m512`, and
#'   `m1024` architectures.}
#'   \item{simudata}{Simulated pseudobulk data returned by [omics_simulate()].}
#'   \item{processed}{Processed matrices returned by [omics_process()].}
#'   \item{models}{Named list of trained OmicsTweezer models.}
#'   \item{loss_history}{Data frame containing epoch-level training losses for
#'   all architectures.}
#' }
#'
#' @details
#' The ensemble uses three architectures: `m256`, `m512`, and `m1024`.
#' Predictions are generated independently for each architecture and then
#' averaged sample-wise.
#'
#'
#' @export
omics_tweezer <- function(sc_data,
                          real_bulk,
                          real_bulk_obs = NULL,
                          ot_weight      = 1.0,
                          scale_minmax   = FALSE,        
                          samplenum      = 5000L,
                          n_cells_per_bulk = 500L,
                          sparse         = TRUE,
                          variance_threshold = 0.98,
                          scaler         = c("ss", "mms", "none"),
                          batch_size     = 128L,
                          epochs         = 30L,
                          learning_rate  = 1e-4,
                          seed           = 2021L,
                          celltype_col   = "CellType",
                          assay          = "RNA",
                          slot           = "counts",
                          verbose        = TRUE,device         = c("auto", "cpu", "cuda"),
cuda_index     = NULL) {
  scaler <- match.arg(scaler)
  device <- match.arg(device)
  total_start <- Sys.time()

  ascii <- r"(
    /\     |   ____            _         
   /  \    |  / __ \          (_)        
  / /\ \   | | |  | |_ __ ___  _  ___ ___
 |_|  |_|  | | |  | | '_ ` _ \| |/ __/ __|
   0===0   | | |__| | | | | | | | (__\__ \
    O=o    |  \____/|_|_|_| |_|_|\___|___/
     O     |     |__   __|
    o=O    |        | |__      _____  ___ _______ _ __ 
   0===0   |        | |\ \ /\ / / _ \/ _ \_  / _ \ '__|
    O=o    |        | | \ V  V /  __/  __// /  __/ |   
     O     |        |_|  \_/\_/ \___|\___/___\___|_|   
    o=O    |
   0===0   | Optimal transport-based cell-type deconvolution
  )"
  if (verbose) cat("\033[36m\033[1m", ascii, "\033[0m\n", sep = "")
  ot_log("OmicsTweezer workflow started", start_time = total_start)

  # 1. simulate pseudobulks
  stage <- Sys.time()
  ot_log("Stage: simulation", start_time = total_start)
  simudata <- omics_simulate(
    sc_data        = sc_data,
    samplenum      = samplenum,
    n              = n_cells_per_bulk,
    sparse         = sparse,
    random_state   = seed,
    celltype_col   = celltype_col,
    assay          = assay,
    slot           = slot
  )
  ot_log("Stage finished: simulation | %s",
         ot_format_hms(as.numeric(difftime(Sys.time(), stage, units = "secs"))),
         start_time = total_start)

  # 2. optional outer MinMax on raw counts (mirrors train_predict(scale=TRUE))
  if (scale_minmax) {
    ot_log("Outer MinMax scaling (mirrors train_predict scale=TRUE)",
           start_time = total_start)
    rescale <- function(M) {
      v <- as.vector(M); lo <- min(v); hi <- max(v); rng <- hi - lo
      if (rng == 0) return(M * 0)
      (M - lo) / rng
    }
    simudata$X <- rescale(simudata$X)
    real_bulk  <- rescale(as.matrix(real_bulk))
  }

  # 3. process
  stage <- Sys.time()
  ot_log("Stage: processing", start_time = total_start)
  processed <- omics_process(simudata, real_bulk,
                             variance_threshold = variance_threshold,
                             scaler = scaler)
  ot_log("Stage finished: processing | %s",
         ot_format_hms(as.numeric(difftime(Sys.time(), stage, units = "secs"))),
         start_time = total_start)
  ot_log("Training data shape: %d x %d", nrow(processed$train_x),
         ncol(processed$train_x), start_time = total_start)
  ot_log("Test data shape: %d x %d", nrow(processed$test_x),
         ncol(processed$test_x), start_time = total_start)

  # 4. ensemble: three architectures, averaged
  archs <- list(
    m256  = list(dims = c(256L, 128L, 64L,  32L),  drops = c(0,   0,   0,   0)),
    m512  = list(dims = c(512L, 256L, 128L, 64L),  drops = c(0, 0.3, 0.2, 0.1)),
    m1024 = list(dims = c(1024L, 512L, 256L, 128L), drops = c(0, 0.6, 0.3, 0.1))
  )

  per_model <- list()
  models <- list()
  loss_history <- list()
  for (nm in names(archs)) {
    stage <- Sys.time()
    ot_log("Stage: training %s", nm, start_time = total_start)
    m <- omics_train(
      train_x       = processed$train_x,
      train_y       = processed$train_y,
      test_x        = processed$test_x,
      dims          = archs[[nm]]$dims,
      drops         = archs[[nm]]$drops,
      epochs        = epochs,
      batch_size    = batch_size,
      learning_rate = learning_rate,
      loss_weight   = ot_weight,
      seed          = seed,
      verbose       = verbose,
      cuda_index    = cuda_index,device        = device

    )
    ot_log("Stage finished: %s | %s", nm,
           ot_format_hms(as.numeric(difftime(Sys.time(), stage, units = "secs"))),
           start_time = total_start)

    pred <- omics_predict(m, processed$test_x,
                          celltypes  = processed$celltypes,
                          samplename = processed$samplename,
                          batch_size = batch_size,device     = device,verbose = verbose,
  cuda_index = cuda_index)
    per_model[[nm]] <- pred
    models[[nm]]    <- m
    loss_history[[nm]] <- data.frame(
      model = nm,
      m$history,
      row.names = NULL
    )
  }

  loss_history <- do.call(rbind, loss_history)
  rownames(loss_history) <- NULL

  # 5. average predictions
  avg <- Reduce(`+`, per_model) / length(per_model)
  rownames(avg) <- processed$samplename
  colnames(avg) <- processed$celltypes

  if (!is.null(real_bulk_obs))
    avg <- list(pred = avg,
                ground_truth = as.matrix(real_bulk_obs)[rownames(avg), , drop = FALSE])

  ot_log("Workflow finished | total %s",
         ot_format_hms(as.numeric(difftime(Sys.time(), total_start, units = "secs"))),
         start_time = total_start)

  return(list(
    pred = avg,
    per_model = per_model,
    simudata = simudata,
    processed = processed,
    models = models,
    loss_history = loss_history
  ))
}





