

#' Extract single-cell counts and cell-type labels for TAPE
#'
#' Internal helper that extracts a cells x genes expression matrix and the
#' corresponding cell-type labels from a `Seurat`,
#' `SingleCellExperiment`, or matrix-like object.
#'
#' @param sc_data A `Seurat` object, a `SingleCellExperiment` object, or a
#'   matrix/data.frame. For matrix-like input, rows are assumed to be cells and
#'   columns genes.
#' @param celltype_col Character scalar giving the metadata column that contains
#'   cell-type labels.
#' @param assay Character scalar giving the assay name for `Seurat` input.
#' @param slot Character scalar giving the assay slot or assay name to extract.
#'
#' @return A named list with components:
#' \describe{
#'   \item{counts}{A numeric matrix with cells in rows and genes in columns.}
#'   \item{celltypes}{A character vector of cell-type labels, one per cell.}
#'   \item{genes}{A character vector of gene names.}
#' }
#'
#' @importFrom SummarizedExperiment assay colData
#' @importFrom SeuratObject LayerData
#' @keywords internal
#' @noRd
tape_extract_sc_data <- function(sc_data,
                                 celltype_col = "CellType",
                                 assay = "RNA",
                                 slot = "counts") {
  if (inherits(sc_data, "Seurat")) {
    if (!assay %in% names(sc_data@assays)) {
      stop("Assay '", assay, "' not found in Seurat object.")
    }

    meta <- sc_data@meta.data
    if (!celltype_col %in% colnames(meta)) {
      stop("Column '", celltype_col, "' not found in Seurat meta.data.")
    }

    mat <- NULL

 
    if ("LayerData" %in% getNamespaceExports("SeuratObject")) {
      mat <- tryCatch(
        SeuratObject::LayerData(sc_data, assay = assay, layer = slot),
        error = function(e) NULL
      )
    }


    if (is.null(mat)) {
      assay_obj <- sc_data[[assay]]

      if (slot %in% c("counts", "data", "scale.data")) {
        mat <- switch(
          slot,
          counts = assay_obj$counts,
          data = assay_obj$data,
          scale.data = assay_obj$scale.data
        )
      } else if (slot %in% methods::slotNames(assay_obj)) {
        mat <- methods::slot(assay_obj, slot)
      } else {
        stop(
          "Could not extract slot/layer '", slot, "' from assay '", assay, "'."
        )
      }
    }

    mat <- as.matrix(mat)
    mat <- t(mat)  # genes x cells -> cells x genes
    ct <- as.character(meta[[celltype_col]])

    genes <- colnames(mat)

  } else if (inherits(sc_data, "SingleCellExperiment")) {
    avail <- SummarizedExperiment::assayNames(sc_data)
    if (!slot %in% avail) {
      stop("Assay '", slot, "' not found in SingleCellExperiment object.")
    }

    meta <- as.data.frame(SummarizedExperiment::colData(sc_data))
    if (!celltype_col %in% colnames(meta)) {
      stop("Column '", celltype_col, "' not found in colData.")
    }

    mat <- t(as.matrix(SummarizedExperiment::assay(sc_data, slot)))
    ct <- as.character(meta[[celltype_col]])
    genes <- colnames(mat)

  } else if (is.matrix(sc_data) || is.data.frame(sc_data)) {
    mat <- as.matrix(sc_data)
    ct <- rownames(mat)

    if (is.null(ct)) {
      stop("Matrix input must have rownames as cell type labels.")
    }

    genes <- colnames(mat)

  } else {
    stop("sc_data must be Seurat, SingleCellExperiment, or matrix/data.frame.")
  }

  if (is.null(genes)) {
    genes <- paste0("gene_", seq_len(ncol(mat)))
    colnames(mat) <- genes
  }

  list(
    counts = mat,
    celltypes = ct,
    genes = genes
  )
}


















#' Simulate TAPE training pseudobulks from single-cell data
#'
#' Generates simulated bulk RNA-seq mixtures from single-cell reference data
#' using Dirichlet-sampled cell-type proportions, optional sparsity, and
#' optional rare-cell perturbations.
#'
#' @param sc_data A `Seurat` object, a `SingleCellExperiment` object, or a
#'   matrix/data.frame with cells in rows and genes in columns.
#' @param d_prior Numeric vector or `NULL`. Dirichlet prior for cell-type
#'   fractions. If `NULL`, a vector of ones is used.
#' @param n Integer scalar. Number of cells per simulated bulk sample.
#' @param samplenum Integer scalar. Number of pseudobulk samples to generate.
#' @param random_state Integer scalar or `NULL`. Random seed for reproducible
#'   simulation.
#' @param sparse Logical scalar. Whether to generate sparse cell-type mixtures.
#' @param sparse_prob Numeric scalar in `[0, 1]`. Controls both the proportion
#'   of sparse samples and the proportion of cell types set to zero within those
#'   samples, matching the provided implementation.
#' @param rare Logical scalar. Whether to perturb selected cell types to very
#'   small fractions.
#' @param rare_percentage Numeric scalar in `[0, 1]`. Fraction of cell types
#'   selected for rare-cell perturbation.
#' @param celltype_col Character scalar giving the metadata column with
#'   cell-type labels.
#' @param assay Character scalar giving the assay name for `Seurat` input.
#' @param slot Character scalar giving the assay slot or assay name to extract.
#'
#' @return A named list with components:
#' \describe{
#'   \item{X}{A numeric matrix of simulated pseudobulks with samples in rows and
#'   genes in columns.}
#'   \item{obs}{A data frame of simulated cell-type proportions with samples in
#'   rows and cell types in columns.}
#'   \item{var}{A data frame whose row names are the simulated gene names.}
#' }
#'
#' @details
#' This function generates simulated pseudobulk mixtures for the R/torch
#' implementation of TAPE, following the simulation strategy used in the
#' original TAPE workflow.
#'
#' @importFrom MCMCpack rdirichlet
#' @importFrom stats runif
#' @examples
#' 
#' set.seed(1)
#' n_genes        <- 1000
#' cell_types     <- c("Tcell", "Bcell", "Mono")
#' cells_per_type <- 100
#' n_cells        <- length(cell_types) * cells_per_type
#'
#' counts <- matrix(
#'   rpois(n_genes * n_cells, lambda = 5),
#'   nrow = n_genes, ncol = n_cells
#' )
#' rownames(counts) <- paste0("gene", seq_len(n_genes))
#' colnames(counts) <- paste0("cell", seq_len(n_cells))
#'
#' meta <- data.frame(
#'   cell_type = rep(cell_types, each = cells_per_type),
#'   row.names = colnames(counts)
#' )
#' sc <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
#'
#' 
#' simudata <- tape_simulate(
#'   sc_data = sc,
#'   samplenum = 50,
#'   n = 100,
#'   celltype_col = "cell_type"
#' )
#' 
#' print(names(simudata))
#' 
#' @export
tape_simulate <- function(sc_data,
                          d_prior = NULL,
                          n = 500,
                          samplenum = 5000,
                          random_state = NULL,
                          sparse = TRUE,
                          sparse_prob = 0.5,
                          rare = FALSE,
                          rare_percentage = 0.4,
                          celltype_col = "CellType",
                          assay = "RNA",
                          slot = "counts") {

  sim_start <- Sys.time()
  tape_log("Simulation started", start_time = sim_start)

  if (!is.null(random_state)) set.seed(random_state)

  ref <- tape_extract_sc_data(
    sc_data = sc_data,
    celltype_col = celltype_col,
    assay = assay,
    slot = slot
  )
  sc_counts <- ref$counts   # cells x genes
  celltypes <- ref$celltypes
  genename  <- ref$genes

  celltype_levels <- sort(unique(celltypes), method = "radix")
  num_celltype <- length(celltype_levels)

  if (is.null(d_prior)) d_prior <- rep(1, num_celltype)
  if (length(d_prior) != num_celltype) {
    stop("Length of `d_prior` must equal the number of cell types.")
  }

  celltype_groups <- split(
    seq_len(nrow(sc_counts)),
    factor(celltypes, levels = celltype_levels)
  )

  tape_log("Generating Dirichlet proportions", start_time = sim_start)

  # Dirichlet sampling
  prop <- MCMCpack::rdirichlet(samplenum, d_prior)
  prop <- prop / rowSums(prop)

  # Sparse fractions
  if (sparse) {
    tape_log("Applying sparse perturbation", start_time = sim_start)
    n_sparse_samples <- as.integer(nrow(prop) * sparse_prob)
    n_zero_ct <- as.integer(ncol(prop) * sparse_prob)
    if (n_sparse_samples > 0 && n_zero_ct > 0) {
      n_zero_ct <- min(n_zero_ct, ncol(prop) - 1L)
      for (i in seq_len(n_sparse_samples)) {
        idx_zero <- sample.int(ncol(prop), n_zero_ct, replace = FALSE)
        prop[i, idx_zero] <- 0
      }
      prop <- prop / rowSums(prop)
    }
  }

  # Rare cell types
  if (rare) {
    tape_log("Applying rare-cell perturbation", start_time = sim_start)
    set.seed(0)
    n_rare_ct <- as.integer(ncol(prop) * rare_percentage)
    if (n_rare_ct > 0) {
      rare_idx <- sample.int(ncol(prop), n_rare_ct, replace = FALSE)
      prop <- prop / rowSums(prop)
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

  # Discretize to cell counts
  cell_num <- floor(n * prop)

  # Guard against all-zero rows
  row_sums <- rowSums(cell_num)
  zero_rows <- which(row_sums == 0)
  if (length(zero_rows) > 0) {
    best_ct <- max.col(prop[zero_rows, , drop = FALSE], ties.method = "first")
    for (i in seq_along(zero_rows)) {
      cell_num[zero_rows[i], best_ct[i]] <- 1
    }
  }

  prop <- cell_num / rowSums(cell_num)

  # Sample pseudo-bulk via sparse matrix aggregation
  tape_log("Sampling cells and building pseudobulks", start_time = sim_start)

  sample_mat <- matrix(0, nrow = samplenum, ncol = ncol(sc_counts))
  pb <- utils::txtProgressBar(min = 0, max = length(celltype_levels), style = 3)

  for (j in seq_along(celltype_levels)) {
    idx_pool <- celltype_groups[[j]]
    counts_j <- as.integer(cell_num[, j])
    if (length(idx_pool) == 0 || sum(counts_j) == 0) {
      utils::setTxtProgressBar(pb, j)
      next
    }

    sample_ids <- rep.int(seq_len(samplenum), times = counts_j)
    local_ids  <- sample.int(length(idx_pool), size = length(sample_ids), replace = TRUE)

    W <- Matrix::sparseMatrix(
      i = local_ids,
      j = sample_ids,
      x = 1,
      dims = c(length(idx_pool), samplenum)
    )
    contrib <- Matrix::t(Matrix::t(sc_counts[idx_pool, , drop = FALSE]) %*% W)
    sample_mat <- sample_mat + as.matrix(contrib)
    utils::setTxtProgressBar(pb, j)
  }
  close(pb)

  sample_names <- paste0("sample_", seq_len(samplenum))
  rownames(sample_mat) <- sample_names
  colnames(sample_mat) <- genename

  prop_df <- as.data.frame(prop)
  colnames(prop_df) <- celltype_levels
  rownames(prop_df) <- sample_names

  tape_log("Simulation finished", start_time = sim_start)

  return(list(
    X   = sample_mat,
    obs = prop_df,
    var = data.frame(row.names = genename)
  ))
}
















#' Process simulated and bulk data for TAPE
#'
#' Prepares simulated training data and real bulk data for TAPE model training
#' and prediction. The function performs variance filtering, gene intersection,
#' log transformation, and per-sample scaling.
#'
#' @param simudata A list returned by [tape_simulate()] containing simulated
#'   pseudobulk expression and proportions.
#' @param real_bulk Numeric matrix with genes in rows and samples in columns.
#' @param variance_threshold Numeric scalar in `[0, 1]`. Fractional rank used
#'   to define the variance cutoff independently for simulated and real bulk
#'   data.
#' @param scaler Character scalar. Either `"mms"` for min-max scaling or `"ss"`
#'   for standard scaling.
#'
#' @return A named list with components:
#' \describe{
#'   \item{train_x}{Numeric matrix of simulated training expression with samples
#'   in rows and genes in columns.}
#'   \item{train_y}{Numeric matrix of simulated proportions with samples in rows
#'   and cell types in columns.}
#'   \item{test_x}{Numeric matrix of processed real bulk expression with samples
#'   in rows and genes in columns.}
#'   \item{genename}{Character vector of intersected gene names.}
#'   \item{celltypes}{Character vector of cell-type names.}
#'   \item{samplename}{Character vector of real bulk sample names.}
#' }
#'
#' @details
#' The current implementation applies scaling row-wise after `log1p`
#' transformation, matching the intended sample-wise scaling direction.
#'
#' @importFrom stats var sd
#' @examples
#' 
#' ## Build a tiny single-cell reference (genes x cells) ---------------------
#' set.seed(1)
#' n_genes        <- 1000
#' cell_types     <- c("Tcell", "Bcell", "Mono")
#' cells_per_type <- 100
#' n_cells        <- length(cell_types) * cells_per_type
#'
#' counts <- matrix(
#'   rpois(n_genes * n_cells, lambda = 5),
#'   nrow = n_genes, ncol = n_cells
#' )
#' rownames(counts) <- paste0("gene", seq_len(n_genes))
#' colnames(counts) <- paste0("cell", seq_len(n_cells))
#'
#' meta <- data.frame(
#'   cell_type = rep(cell_types, each = cells_per_type),
#'   row.names = colnames(counts)
#' )
#' sc <- SeuratObject::CreateSeuratObject(counts = counts, meta.data = meta)
#'
#' ## Simulate training pseudobulks (samples x genes) ------------------------
#' train_sim <- tape_simulate(
#'   sc, celltype_col = "cell_type",
#'   samplenum = 50, n = 100, sparse = FALSE, random_state = 1
#' )
#'
#' ## Simulate a separate set to stand in for the real bulk ------------------
#' test_sim <- tape_simulate(
#'   sc, celltype_col = "cell_type",
#'   samplenum = 10, n = 100, sparse = FALSE, random_state = 2
#' )
#' real_bulk <- t(test_sim$X)   # tape_process expects genes x samples
#'
#' ## Variance filter, gene intersection, log1p, per-sample scaling ----------
#' processed <- tape_process(train_sim, real_bulk, scaler = "mms")
#' head(processed$train_x)
#' head(processed$test_x)
#' 
#' @export
tape_process <- function(simudata, real_bulk,
                         variance_threshold = 0.98,
                         scaler = "mms") {

 
  proc_start <- Sys.time()
  tape_log("Processing started", start_time = proc_start)
 
  train_x <- simudata$X      
  train_y <- as.matrix(simudata$obs)
  test_x  <- t(as.matrix(real_bulk))    # samples x genes
 
  tape_log("Applying variance filtering", start_time = proc_start)
 
  pick_cutoff <- function(v, vt) {
    k <- as.integer(length(v) * vt) + 1L          # 0-based int(n*vt) -> 1-based
    k <- min(k, length(v))                         # guard (vt < 1 keeps k <= n)
    sort(v, decreasing = TRUE)[k]
  }
 
  train_var <- apply(train_x, 2, var)
  cutoff_train <- pick_cutoff(train_var, variance_threshold)
  train_x <- train_x[, train_var > cutoff_train, drop = FALSE]
 
  test_var <- apply(test_x, 2, var)
  cutoff_test <- pick_cutoff(test_var, variance_threshold)
  test_x <- test_x[, test_var > cutoff_test, drop = FALSE]
 

  inter <- intersect(colnames(train_x), colnames(test_x))
  train_x <- train_x[, inter, drop = FALSE]
  test_x  <- test_x[, inter, drop = FALSE]
  tape_log("Intersected gene number: %d", length(inter), start_time = proc_start)
 
  genename   <- inter
  celltypes  <- colnames(train_y)
  samplename <- rownames(test_x)
 

  tape_log("Applying log1p transform", start_time = proc_start)
  train_x <- log1p(as.matrix(train_x))
  test_x  <- log1p(as.matrix(test_x))
 
  tape_log("Applying %s scaling", scaler, start_time = proc_start)
 

  if (scaler == "mms") {
    # sklearn MinMaxScaler().fit_transform(x.T).T  ==  (x - rowMin)/(rowMax-rowMin)
    scale_minmax <- function(x) {
      mn  <- apply(x, 1, min)
      mx  <- apply(x, 1, max)
      rng <- mx - mn
      rng_safe <- rng; rng_safe[rng_safe == 0] <- 1
      x <- sweep(x, 1, mn, "-")
      x <- sweep(x, 1, rng_safe, "/")
      if (any(rng == 0)) x[rng == 0, ] <- 0   # constant rows -> 0 (matches sklearn)
      x
    }
    train_x <- scale_minmax(train_x)
    test_x  <- scale_minmax(test_x)
 
  } else if (scaler == "ss") {
    # sklearn StandardScaler().fit_transform(x.T).T  ==  (x - rowMean)/popSD
    # popSD divides by n (ddof = 0); base-R sd() (ddof = 1) is NOT used.
    scale_standard <- function(x) {
      mu  <- rowMeans(x)
      sdv <- sqrt(rowMeans((x - mu)^2))        # population std, ddof = 0
      sd_safe <- sdv; sd_safe[sd_safe == 0] <- 1
      x <- sweep(x, 1, mu, "-")
      x <- sweep(x, 1, sd_safe, "/")
      if (any(sdv == 0)) x[sdv == 0, ] <- 0    # constant rows -> 0 (matches sklearn)
      x
    }
    train_x <- scale_standard(train_x)
    test_x  <- scale_standard(test_x)
 
  } else {
    stop("`scaler` must be 'mms' or 'ss'.")
  }
 
  tape_log("Processing finished", start_time = proc_start)
 
  list(
    train_x    = train_x,
    train_y    = as.matrix(train_y),
    test_x     = test_x,
    genename   = genename,
    celltypes  = celltypes,
    samplename = samplename
  )
}








#' Create the TAPE autoencoder modules
#'
#' Internal helper that constructs the encoder and decoder used by TAPE.
#'
#' @param input_dim Integer scalar. Number of genes.
#' @param output_dim Integer scalar. Number of cell types.
#'
#' @return A named list with `encoder` and `decoder` torch modules.
#'
#' @importFrom torch nn_sequential nn_dropout nn_linear nn_celu
#' @keywords internal
#' @noRd
tape_create_model <- function(input_dim, output_dim) {
  encoder <- nn_sequential(
    nn_dropout(p = 0.5),
    nn_linear(input_dim, 512),
    nn_celu(),
    nn_dropout(p = 0.5),
    nn_linear(512, 256),
    nn_celu(),
    nn_dropout(p = 0.5),
    nn_linear(256, 128),
    nn_celu(),
    nn_dropout(p = 0.5),
    nn_linear(128, 64),
    nn_celu(),
    nn_linear(64, output_dim)
  )

  decoder <- nn_sequential(
    nn_linear(output_dim, 64, bias = FALSE),
    nn_linear(64, 128, bias = FALSE),
    nn_linear(128, 256, bias = FALSE),
    nn_linear(256, 512, bias = FALSE),
    nn_linear(512, input_dim, bias = FALSE)
  )

  list(encoder = encoder, decoder = decoder)
}

#' Compute the TAPE signature matrix
#'
#' Internal helper that multiplies decoder weights to obtain the non-negative
#' signature matrix.
#'
#' @param model A model list returned by [tape_create_model()] or trained by
#'   [tape_train()].
#'
#' @return A `torch_tensor` containing the signature matrix.
#'
#' @importFrom torch torch_mm nnf_relu
#' @keywords internal
#' @noRd
tape_sigmatrix <- function(model) {
  w0 <- model$decoder[[1]]$weight$t()
  w1 <- model$decoder[[2]]$weight$t()
  w2 <- model$decoder[[3]]$weight$t()
  w3 <- model$decoder[[4]]$weight$t()
  w4 <- model$decoder[[5]]$weight$t()

  w <- torch_mm(torch_mm(torch_mm(torch_mm(w0, w1), w2), w3), w4)
  nnf_relu(w)
}

#' Forward pass for the TAPE autoencoder
#'
#' Internal helper that computes reconstructed input, latent proportions, and
#' signature matrix.
#'
#' @param model A model list returned by [tape_create_model()] or [tape_train()].
#' @param x A `torch_tensor` with samples in rows and genes in columns.
#' @param state Character scalar. Either `"train"` or `"test"`.
#'
#' @return A named list with `x_recon`, `z`, and `sigm`.
#'
#' @importFrom torch torch_mm
#' @keywords internal
#' @noRd
tape_forward <- function(model, x, state = "train") {
  sigm <- tape_sigmatrix(model)
  z <- model$encoder(x)

  if (state == "test") {
    z <- nnf_relu(z)
    z_sum <- z$sum(dim = 2, keepdim = TRUE)
    z <- z / (z_sum)
  }

  x_recon <- torch_mm(z, sigm)
  list(x_recon = x_recon, z = z, sigm = sigm)
}

#' Get all trainable parameters of the TAPE model
#'
#' @param model A model list returned by [tape_create_model()] or [tape_train()].
#'
#' @return A list of torch parameters.
#'
#' @keywords internal
#' @noRd
tape_get_all_parameters <- function(model) {
  c(model$encoder$parameters, model$decoder$parameters)
}

#' Get encoder parameters of the TAPE model
#'
#' @param model A model list returned by [tape_create_model()] or [tape_train()].
#'
#' @return A list of torch parameters from the encoder.
#'
#' @keywords internal
#' @noRd
tape_get_encoder_parameters <- function(model) {
  model$encoder$parameters
}

#' Get decoder parameters of the TAPE model
#'
#' @param model A model list returned by [tape_create_model()] or [tape_train()].
#'
#' @return A list of torch parameters from the decoder.
#'
#' @keywords internal
#' @noRd
tape_get_decoder_parameters <- function(model) {
  model$decoder$parameters
}

#' Set random seeds for TAPE
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
tape_set_seed <- function(seed = 1) {
  torch_manual_seed(seed)
  set.seed(seed)
}

#' @keywords internal
#' @noRd
tape_format_hms <- function(seconds) {
  seconds <- max(0L, as.integer(round(seconds)))
  hh <- seconds %/% 3600L
  mm <- (seconds %% 3600L) %/% 60L
  ss <- seconds %% 60L
  sprintf("%02d:%02d:%02d", hh, mm, ss)
}

#' @keywords internal
#' @noRd
tape_log <- function(msg, ..., start_time = NULL) {
  prefix <- ""
  if (!is.null(start_time)) {
    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    prefix <- paste0("[", tape_format_hms(elapsed), "] ")
  }
  message(prefix, sprintf(msg, ...))
}

#' @keywords internal
#' @noRd
tape_make_bar <- function(i, n, width = 30L) {
  if (n <= 0L) {
    return(paste(rep("-", width), collapse = ""))
  }
  filled <- as.integer(floor(width * i / n))
  paste0(
    paste(rep("=", filled), collapse = ""),
    paste(rep("-", width - filled), collapse = "")
  )
}

#' Train the TAPE autoencoder
#'
#' Trains the TAPE autoencoder on simulated pseudobulk expression and
#' corresponding cell-type proportions.
#'
#' @param train_x Numeric matrix with samples in rows and genes in columns.
#' @param train_y Numeric matrix with samples in rows and cell types in columns.
#' @param batch_size Integer scalar. Batch size used for training.
#' @param epochs Integer scalar. Number of training epochs.
#' @param seed Integer scalar. Random seed used for reproducibility.
#'
#' @return A trained model list containing `encoder` and `decoder`.
#' 
#' @details
#' This function implements the TAPE autoencoder training procedure described
#' by Chen et al. (2022) in torch in R.
#'
#' @importFrom torch cuda_is_available torch_tensor optim_adam nnf_l1_loss
#' @examples
#' \dontrun{
#' model <- tape_train(
#'   train_x = processed$train_x,
#'   train_y = processed$train_y,
#'   epochs = 128
#' )
#' }
#' @export
tape_train <- function(train_x, train_y,
                       batch_size = 128L,
                       epochs = 128L,
                       seed = 0L) {
  train_start <- Sys.time()
  tape_log("Training started", start_time = train_start)

  tape_set_seed(seed)

  input_dim <- ncol(train_x)
  output_dim <- ncol(train_y)
  model <- tape_create_model(input_dim, output_dim)

  dev <- if (cuda_is_available()) "cuda" else "cpu"
  model$encoder <- model$encoder$to(device = dev)
  model$decoder <- model$decoder$to(device = dev)

  optimizer <- optim_adam(tape_get_all_parameters(model), lr = 1e-4)

  n <- nrow(train_x)
  train_x_t <- torch_tensor(train_x, device = dev)
  train_y_t <- torch_tensor(train_y, device = dev)

  model$encoder$train()
  model$decoder$train()

  epoch_loss <- numeric(epochs)
  model_name <- "TAPE"

  for (i in seq_len(epochs)) {
    perm <- sample.int(n)
    batch_losses <- numeric(0)

    for (start in seq(1, n, by = batch_size)) {
      end <- min(start + batch_size - 1L, n)
      idx <- perm[start:end]

      data <- train_x_t[idx, , drop = FALSE]
      label <- train_y_t[idx, , drop = FALSE]

      optimizer$zero_grad()
      out <- tape_forward(model, data, state = "train")
      loss <- nnf_l1_loss(out$z, label) + nnf_l1_loss(out$x_recon, data)
      loss$backward()
      optimizer$step()

      batch_losses <- c(batch_losses, as.numeric(loss$item()))
    }

    epoch_loss[i] <- mean(batch_losses)

    elapsed <- as.numeric(difftime(Sys.time(), train_start, units = "secs"))
    avg_per_epoch <- elapsed / i
    eta <- avg_per_epoch * (epochs - i)

    cat(sprintf(
      "\r%s | [%s] %3d/%3d | mean loss %.6f | elapsed %s | ETA %s",
      model_name,
      tape_make_bar(i, epochs, width = 30L),
      i, epochs,
      epoch_loss[i],
      tape_format_hms(elapsed),
      tape_format_hms(eta)
    ))
    utils::flush.console()
  }

  cat("\n")
  tape_log("Training finished", start_time = train_start)

  model
}


#' Clone a trained TAPE model
#'
#' Internal helper that creates an independent copy of a trained TAPE model by
#' copying encoder and decoder state dictionaries in memory.
#'
#' @param model A trained TAPE model.
#'
#' @return A new model list with the same weights as `model`.
#'
#' @importFrom torch cuda_is_available
#' @export
tape_clone_model <- function(model) {
  input_dim <- model$encoder[[2]]$in_features
  output_dim <- model$encoder[[13]]$out_features

  new_model <- tape_create_model(input_dim, output_dim)

  dev <- if (cuda_is_available()) "cuda" else "cpu"
  new_model$encoder <- new_model$encoder$to(device = dev)
  new_model$decoder <- new_model$decoder$to(device = dev)

  new_model$encoder$load_state_dict(model$encoder$state_dict())
  new_model$decoder$load_state_dict(model$decoder$state_dict())

  new_model
}



#' Run the adaptive refinement stage for TAPE
#'
#' Internal helper that performs alternating decoder and encoder refinement on
#' real bulk expression.
#'
#' @param model A trained TAPE model.
#' @param data_np Numeric matrix with samples in rows and genes in columns.
#' @param step Integer scalar. Number of optimization steps per sub-stage.
#' @param max_iter Integer scalar. Number of alternating refinement rounds.
#' @param verbose Logical scalar. Whether to print per-step progress.
#'
#' @return A named list with `sigm` and `pred`.
#'
#' @importFrom torch cuda_is_available torch_tensor optim_adam nnf_l1_loss
#' @keywords internal
#' @noRd
tape_adaptive_stage <- function(model, data_np, step = 300L, max_iter = 5L,
                                verbose = TRUE) {
  adapt_start <- Sys.time()
 
  dev <- if (cuda_is_available()) "cuda" else "cpu"
  data <- torch_tensor(data_np, device = dev)
 
  model$encoder$eval()
  model$decoder$eval()
 
  out <- with_no_grad(tape_forward(model, data, state = "test"))
  ori_sigm <- out$sigm$detach()$clone()
  ori_pred <- out$z$detach()$clone()
  rm(out)
 
  optD <- optim_adam(tape_get_decoder_parameters(model), lr = 1e-4)
  optE <- optim_adam(tape_get_encoder_parameters(model), lr = 1e-4)
 
  masks <- tape_make_masks(
    dims = c(ncol(data_np), 512L, 256L, 128L),
    n_rows = nrow(data_np),
    dev = dev
  )
 
  total_steps <- max_iter * step * 2L
  step_counter <- 0L
  model_name <- "TAPE-adapt"
 
  report <- function(value) {
    elapsed <- as.numeric(difftime(Sys.time(), adapt_start, units = "secs"))
    eta <- (elapsed / step_counter) * (total_steps - step_counter)
    cat(sprintf(
      "\r%s | [%s] %4d/%4d | mean loss %.6f | elapsed %s | ETA %s",
      model_name, tape_make_bar(step_counter, total_steps, width = 30L),
      step_counter, total_steps, value,
      tape_format_hms(elapsed), tape_format_hms(eta)
    ))
    utils::flush.console()
  }
 
  for (k in seq_len(max_iter)) {
    model$encoder$train()
    model$decoder$train()
 
    z_fixed <- with_no_grad(tape_masked_encode(model$encoder, data, masks))
 
    for (i in seq_len(step)) {
      optD$zero_grad()
      sigm <- tape_sigmatrix(model)
      loss <- nnf_l1_loss(torch_mm(z_fixed, sigm), data) +
        nnf_l1_loss(sigm, ori_sigm)
      loss$backward()
      optD$step()
 
      if (verbose) {
        step_counter <- step_counter + 1L
        report(as.numeric(loss$item()))
      }
    }
 
    rm(z_fixed)
    sigm_fixed <- with_no_grad(tape_sigmatrix(model))
 
    for (i in seq_len(step)) {
      optE$zero_grad()
      z <- tape_masked_encode(model$encoder, data, masks)
      loss <- nnf_l1_loss(ori_pred, z) +
        nnf_l1_loss(torch_mm(z, sigm_fixed), data)
      loss$backward()
      optE$step()
 
      if (verbose) {
        step_counter <- step_counter + 1L
        report(as.numeric(loss$item()))
      }
    }
 
    rm(sigm_fixed)
  }
 
  if (verbose) cat("\n")
 
  model$encoder$eval()
  model$decoder$eval()
  out <- with_no_grad(tape_forward(model, data, state = "test"))
 
  list(
    sigm = as.matrix(out$sigm$cpu()$detach()),
    pred = as.matrix(out$z$cpu()$detach())
  )
}













#' Predict cell fractions and signature matrices with TAPE
#'
#' Applies a trained TAPE model to processed bulk data, optionally with adaptive
#' refinement in either overall or high-resolution mode.
#'
#' @param model A trained model returned by [tape_train()].
#' @param test_x Numeric matrix with samples in rows and genes in columns.
#' @param genename Character vector of gene names.
#' @param celltypes Character vector of cell-type names.
#' @param samplename Character vector of sample names.
#' @param adaptive Logical scalar. Whether to run adaptive refinement.
#' @param mode Character scalar. Either `"overall"` or `"high-resolution"`.
#' @param chunk_size Integer scalar. Number of samples adapted simultaneously in
#'   high-resolution mode.
#'
#' @return A named list with components:
#' \describe{
#'   \item{sigm}{A named list of per-cell-type signature matrices in
#'   high-resolution mode, a signature matrix data frame in overall adaptive
#'   mode, and `NULL` when `adaptive` is `FALSE`.}
#'   \item{pred}{A data frame of predicted proportions with samples in rows and
#'   cell types in columns.}
#' }
#'
#' @details
#' High-resolution mode adapts every sample with its own copy of the trained
#' model, and processes `chunk_size` samples concurrently as an ensemble of
#' independent models. The supplied model is read only in this mode and is
#' returned unmodified.
#'
#' @importFrom torch cuda_is_available  torch_tensor with_no_grad
#' @examples
#' \dontrun{
#' res <- tape_predict(
#'   model = model,
#'   test_x = processed$test_x,
#'   genename = processed$genename,
#'   celltypes = processed$celltypes,
#'   samplename = processed$samplename,
#'   mode = "high-resolution",
#'   chunk_size = 16L
#' )
#' }
#' @export
tape_predict <- function(model, test_x, genename, celltypes, samplename,
                         adaptive = TRUE,
                         mode = "overall",
                         chunk_size = 16L) {
 
  pred_start <- Sys.time()
  tape_log("Prediction started | mode = %s | adaptive = %s", mode, adaptive,
           start_time = pred_start)
 
  if (adaptive) {
    if (mode == "high-resolution") {
      n_samples <- nrow(test_x)
      n_ct <- length(celltypes)
      n_genes <- length(genename)
 
      TestSigmList <- array(0, dim = c(n_samples, n_ct, n_genes))
      TestPred <- matrix(0, nrow = n_samples, ncol = n_ct)
 
      chunk_size <- max(1L, as.integer(chunk_size))
      starts <- seq(1L, n_samples, by = chunk_size)
 
      tape_log("Start adaptive training at high-resolution | chunk size %d",
               chunk_size, start_time = pred_start)
      pb <- utils::txtProgressBar(min = 0, max = n_samples, style = 3)
 
      for (s in starts) {
        e <- min(s + chunk_size - 1L, n_samples)
 
        res <- tape_adaptive_chunk(
          model = model,
          chunk_x = test_x[s:e, , drop = FALSE],
          step = 300L,
          max_iter = 3L
        )
 
        TestPred[s:e, ] <- res$pred
        TestSigmList[s:e, , ] <- res$sigm
 
        utils::setTxtProgressBar(pb, e)
        tape_log("High-resolution sample %d/%d finished", e, n_samples,
                 start_time = pred_start)
 
        rm(res)
        gc(verbose = FALSE)

      }
 
      close(pb)
 
      TestPred <- as.data.frame(TestPred)
      colnames(TestPred) <- celltypes
      rownames(TestPred) <- samplename
 
      CellTypeSigm <- list()
      for (j in seq_along(celltypes)) {
        sigm_df <- as.data.frame(TestSigmList[, j, ])
        colnames(sigm_df) <- genename
        rownames(sigm_df) <- samplename
        CellTypeSigm[[celltypes[j]]] <- sigm_df
      }
 
      tape_log("Prediction finished", start_time = pred_start)
      return(list(sigm = CellTypeSigm, pred = TestPred))
 
    } else if (mode == "overall") {
      tape_log("Start adaptive training for all samples", start_time = pred_start)
 
      res <- tape_adaptive_stage(model, test_x, step = 300L, max_iter = 3L,
                                 verbose = TRUE)
 
      sigm <- as.data.frame(res$sigm)
      colnames(sigm) <- genename
      rownames(sigm) <- celltypes
 
      pred <- as.data.frame(res$pred)
      colnames(pred) <- celltypes
      rownames(pred) <- samplename
 
      tape_log("Prediction finished", start_time = pred_start)
      return(list(sigm = sigm, pred = pred))
    }
  } else {
    tape_log("Predicting cell fractions without adaptive training",
             start_time = pred_start)
 
    dev <- if (cuda_is_available()) "cuda" else "cpu"
    model$encoder$eval()
    model$decoder$eval()
    data <- torch_tensor(test_x, device = dev)
    out <- with_no_grad(tape_forward(model, data, state = "test"))
 
    pred <- as.data.frame(as.matrix(out$z$cpu()$detach()))
    colnames(pred) <- celltypes
    rownames(pred) <- samplename
 
    tape_log("Prediction finished", start_time = pred_start)
    return(list(sigm = NULL, pred = pred))
  }
}




#' Run the full TAPE workflow
#'
#' Convenience wrapper that simulates pseudobulks from single-cell data,
#' processes simulated and real bulk data, trains the TAPE autoencoder, and
#' predicts cell-type fractions with optional adaptive refinement.
#'
#' @param sc_data A `Seurat` object, a `SingleCellExperiment` object, or a
#'   matrix/data.frame with cells in rows and genes in columns.
#' @param real_bulk Numeric matrix with genes in rows and samples in columns.
#' @param variance_threshold Numeric scalar in `[0, 1]` used for variance-based
#'   gene filtering.
#' @param scaler Character scalar. Either `"mms"` or `"ss"`.
#' @param d_prior Numeric vector or `NULL`. Dirichlet prior used in simulation.
#' @param mode Character scalar. Either `"overall"` or `"high-resolution"`.
#' @param adaptive Logical scalar. Whether to use adaptive refinement.
#' @param sparse Logical scalar. Whether to simulate sparse mixtures.
#' @param batch_size Integer scalar. Batch size for training.
#' @param epochs Integer scalar. Number of training epochs.
#' @param seed Integer scalar. Random seed.
#' @param samplenum Integer scalar. Number of simulated pseudobulk samples.
#' @param n Integer scalar. Number of cells per simulated pseudobulk.
#' @param celltype_col Character scalar. Metadata column containing cell-type
#'   labels.
#' @param assay Character scalar giving the assay name for `Seurat` input.
#' @param slot Character scalar giving the assay slot or assay name to extract.
#'
#' @details
#' This function provides an R implementation of TAPE (Tissue-AdaPtive
#' autoEncoder) as described by Chen et al. (2022), using torch for model
#' training, adaptive refinement, and prediction. The implementation follows
#' the TAPE methodology for pseudobulk simulation, preprocessing, autoencoder
#' training, and optional tissue-adaptive prediction within an R-based
#' workflow.
#'
#' The model architecture and workflow are implemented in torch in R following
#' the PyTorch implementation provided in the original TAPE repository.
#'
#' @references
#' Chen, Y., Wang, Y., Chen, Y., Cheng, Y., Wei, Y., Li, Y., Wang, J., Wei, Y.,
#' Chan, T.-F., & Li, Y. (2022). Deep autoencoder for interpretable
#' tissue-adaptive deconvolution and cell-type-specific gene analysis.
#' \emph{Nature Communications}, 13(1), 6735.
#'
#' Original TAPE software repository:
#' \url{https://github.com/poseidonchan/TAPE}
#' 
#' @return A named list with components:
#' \describe{
#'   \item{sigm}{Predicted signature matrix output, depending on `mode` and
#'   `adaptive`.}
#'   \item{pred}{Predicted cell-type proportions with samples in rows and cell
#'   types in columns.}
#' }
#'
#' @examples
#' \dontrun{
#' res <- tape(
#'   sc_data = sce,
#'   real_bulk = bulk_mat,
#'   celltype_col = "CellType"
#' )
#' }
#' @export
tape <- function(sc_data, real_bulk,
                 variance_threshold = 0.98,
                 scaler = "mms",
                 d_prior = NULL,
                 mode = "overall",
                 adaptive = TRUE,
                 sparse = TRUE,
                 batch_size = 128L,
                 epochs = 128L,
                 seed = 0L,
                 samplenum = 5000L,
                 n = 500L,
                 celltype_col = "CellType",
                 assay = "RNA",
                 slot = "counts") {
  # sc_data: Seurat, SingleCellExperiment, or matrix (cells x genes with celltypes as rownames)
  # real_bulk: matrix genes x samples
  #
  # Returns: list(sigm, pred)

  tape_ascii <- r"(
  Running...
   _______  _______  _______  _______
  |       ||   _   ||       ||       |
  |_     _||  |_|  ||    _  ||    ___|
    |   |  |       ||   |_| ||   |___
    |   |  |       ||    ___||    ___|
    |   |  |   _   ||   |    |   |___
    |___|  |__| |__||___|    |_______|
  )"
  cat("\033[35m", "\033[1m", tape_ascii, "\033[0m", sep = "")

  total_start <- Sys.time()
  tape_log("TAPE workflow started", start_time = total_start)

  stage_start <- Sys.time()
  tape_log("Stage: simulation", start_time = total_start)
  simudata <- tape_simulate(sc_data,
                            samplenum = samplenum,
                            n = n,
                            d_prior = d_prior,
                            sparse = sparse,
                            random_state = seed,
                            celltype_col = celltype_col,
                            assay = assay,
                            slot = slot)
  tape_log("Stage finished: simulation | step_time = %s",
           tape_format_hms(as.numeric(difftime(Sys.time(), stage_start, units = "secs"))),
           start_time = total_start)

  stage_start <- Sys.time()
  tape_log("Stage: processing", start_time = total_start)
  processed <- tape_process(simudata, real_bulk,
                            variance_threshold = variance_threshold,
                            scaler = scaler)
  tape_log("Stage finished: processing | step_time = %s",
           tape_format_hms(as.numeric(difftime(Sys.time(), stage_start, units = "secs"))),
           start_time = total_start)

  tape_log("Training data shape: %d x %d", nrow(processed$train_x), ncol(processed$train_x), start_time = total_start)
  tape_log("Test data shape: %d x %d", nrow(processed$test_x), ncol(processed$test_x), start_time = total_start)

  stage_start <- Sys.time()
  tape_log("Stage: training", start_time = total_start)
  tape_set_seed(seed)
  model <- tape_train(processed$train_x, processed$train_y,
                      batch_size = batch_size,
                      epochs = epochs,
                      seed = seed)
  tape_log("Stage finished: training | step_time = %s",
           tape_format_hms(as.numeric(difftime(Sys.time(), stage_start, units = "secs"))),
           start_time = total_start)

  stage_start <- Sys.time()
  tape_log("Stage: prediction", start_time = total_start)
  result <- tape_predict(model, processed$test_x,
                         genename = processed$genename,
                         celltypes = processed$celltypes,
                         samplename = processed$samplename,
                         adaptive = adaptive,
                         mode = mode)
  tape_log("Stage finished: prediction | step_time = %s",
           tape_format_hms(as.numeric(difftime(Sys.time(), stage_start, units = "secs"))),
           start_time = total_start)

  tape_log("TAPE workflow finished | total_time = %s",
           tape_format_hms(as.numeric(difftime(Sys.time(), total_start, units = "secs"))),
           start_time = total_start)

  result
}





#' Refract latent scores for an ensemble of TAPE models
#'
#' Internal helper that rectifies and row-normalizes batched latent scores.
#'
#' @param z A `torch_tensor` of shape `(B, 1, celltypes)`.
#'
#' @return A `torch_tensor` of the same shape.
#'
#' @importFrom torch nnf_relu
#' @keywords internal
#' @noRd
tape_batched_refract <- function(z) {
  z <- nnf_relu(z)
  z / z$sum(dim = 3, keepdim = TRUE)
}
 
 
#' Adapt a chunk of samples as an ensemble of independent TAPE models
#'
#' Internal helper that performs per-sample adaptive refinement for several
#' samples at once, each with its own copy of the trained model.
#'
#' @param model A trained TAPE model. The model is read only and is left
#'   unmodified.
#' @param chunk_x Numeric matrix with samples in rows and genes in columns.
#' @param step Integer scalar. Number of optimization steps per sub-stage.
#' @param max_iter Integer scalar. Number of alternating refinement rounds.
#'
#' @return A named list with components:
#' \describe{
#'   \item{pred}{Numeric matrix of predicted proportions, samples by cell types.}
#'   \item{sigm}{Numeric array of adapted signature matrices with dimensions
#'   samples by cell types by genes.}
#' }
#'
#' @details
#' Ensemble member `b` sees only sample `b`, so its gradients depend on no other
#' sample provided the loss is reduced per member before being summed. Each loss
#' term is averaged over the elements belonging to one member, reproducing the
#' gradient scale a lone model would receive, and the per-member losses are then
#' summed so that the batch is a concatenation of independent problems.
#'
#' Parameters and optimizers are constructed fresh for every chunk, reproducing
#' the per-sample model reload and the per-sample optimizer state of the
#' reference implementation.
#'
#' @importFrom torch cuda_is_available torch_tensor torch_bmm with_no_grad
#' @keywords internal
#' @noRd
tape_adaptive_chunk <- function(model, chunk_x, step = 300L, max_iter = 3L) {
  dev <- if (cuda_is_available()) "cuda" else "cpu"
 
  B <- nrow(chunk_x)
  data <- torch_tensor(chunk_x, device = dev)$unsqueeze(2)
 
  params <- tape_batchify(model, B)
  EW <- params$EW
  EB <- params$EB
  DW <- params$DW
 
  masks <- lapply(
    tape_make_masks(
      dims = c(ncol(chunk_x), 512L, 256L, 128L),
      n_rows = 1L,
      dev = dev
    ),
    function(m) m$unsqueeze(1)
  )
 
  anchors <- with_no_grad({
    list(
      pred = tape_batched_refract(
        tape_batched_encode(data, EW, EB, masks = NULL)
      )$clone(),
      sigm = tape_batched_sigmatrix(DW)$clone()
    )
  })
 
  optD <- tape_fused_adam(DW, lr = 1e-4)
  optE <- tape_fused_adam(c(EW, EB), lr = 1e-4)
 
  per_member_mean <- function(x) {
    x$abs()$mean(dim = 3L)$mean(dim = 2L)
  }
 
  for (k in seq_len(max_iter)) {
    z_fixed <- with_no_grad(tape_batched_encode(data, EW, EB, masks))
 
    for (i in seq_len(step)) {
      optD$zero_grad()
      sigm <- tape_batched_sigmatrix(DW)
      loss <- (
        per_member_mean(torch_bmm(z_fixed, sigm) - data) +
          per_member_mean(sigm - anchors$sigm)
      )$sum()
      loss$backward()
      optD$step()
    }
 
    rm(z_fixed)
    sigm_fixed <- with_no_grad(tape_batched_sigmatrix(DW))
 
    for (i in seq_len(step)) {
      optE$zero_grad()
      z <- tape_batched_encode(data, EW, EB, masks)
      loss <- (
        per_member_mean(anchors$pred - z) +
          per_member_mean(torch_bmm(z, sigm_fixed) - data)
      )$sum()
      loss$backward()
      optE$step()
    }
 
    rm(sigm_fixed)
  }
 
  final <- with_no_grad({
    list(
      pred = tape_batched_refract(
        tape_batched_encode(data, EW, EB, masks = NULL)
      ),
      sigm = tape_batched_sigmatrix(DW)
    )
  })
 
  list(
    pred = as.matrix(final$pred$squeeze(2)$cpu()$detach()),
    sigm = as.array(final$sigm$cpu()$detach())
  )
}


#' Compute signature matrices for an ensemble of TAPE models
#'
#' Internal helper that collapses batched decoder weights into one non-negative
#' signature matrix per ensemble member.
#'
#' @param DW A list of five batched decoder weight tensors.
#'
#' @return A `torch_tensor` of shape `(B, celltypes, genes)`.
#'
#' @importFrom torch torch_bmm nnf_relu
#' @keywords internal
#' @noRd
tape_batched_sigmatrix <- function(DW) {
  w <- torch_bmm(DW[[1]], DW[[2]])
  w <- torch_bmm(w, DW[[3]])
  w <- torch_bmm(w, DW[[4]])
  w <- torch_bmm(w, DW[[5]])
  nnf_relu(w)
}
 
 
#' Encode with an ensemble of TAPE encoders
#'
#' Internal helper that applies `B` independent encoders to `B` samples in one
#' batched forward pass.
#'
#' @param x A `torch_tensor` of shape `(B, 1, genes)`.
#' @param EW A list of five batched encoder weight tensors.
#' @param EB A list of five batched encoder bias tensors.
#' @param masks A list of four masks broadcastable over the batch dimension, or
#'   `NULL` to encode without dropout.
#'
#' @return A `torch_tensor` of shape `(B, 1, celltypes)`.
#'
#' @importFrom torch torch_bmm nnf_celu
#' @keywords internal
#' @noRd
tape_batched_encode <- function(x, EW, EB, masks = NULL) {
  h <- x
 
  for (l in 1:4) {
    if (!is.null(masks)) {
      h <- h * masks[[l]]
    }
    h <- nnf_celu(torch_bmm(h, EW[[l]]) + EB[[l]])
  }
 
  torch_bmm(h, EW[[5]]) + EB[[5]]
}


#' Expand TAPE parameters into an independent ensemble
#'
#' Internal helper that replicates the encoder and decoder weights of a trained
#' TAPE model into leading-dimension tensors, giving `B` independent copies that
#' can be adapted simultaneously.
#'
#' @param model A trained TAPE model. The model is read only.
#' @param B Integer scalar. Number of ensemble members.
#'
#' @return A named list with `EW`, `EB`, and `DW`, each a list of leaf tensors
#'   with `requires_grad` enabled.
#'
#' @details
#' Encoder weights are stored transposed with shape `(B, in, out)` and biases
#' with shape `(B, 1, out)`, so that a forward pass becomes a batched matrix
#' product. Transposition permutes elements, and Adam updates each element
#' independently, so an optimizer running over these tensors performs exactly the
#' updates that `B` separate optimizers would perform on `B` separate models.
#'
#' @keywords internal
#' @noRd
tape_batchify <- function(model, B) {
  B <- as.integer(B)
  linears <- tape_encoder_linear_index()
 
  EW <- lapply(linears, function(i) {
    model$encoder[[i]]$weight$detach()$t()$unsqueeze(1)$
      expand(c(B, -1L, -1L))$contiguous()$requires_grad_(TRUE)
  })
 
  EB <- lapply(linears, function(i) {
    model$encoder[[i]]$bias$detach()$view(c(1L, 1L, -1L))$
      expand(c(B, 1L, -1L))$contiguous()$requires_grad_(TRUE)
  })
 
  DW <- lapply(1:5, function(i) {
    model$decoder[[i]]$weight$detach()$t()$unsqueeze(1)$
      expand(c(B, -1L, -1L))$contiguous()$requires_grad_(TRUE)
  })
 
  list(EW = EW, EB = EB, DW = DW)
}


#' Positions of the linear layers inside the TAPE encoder
#'
#' Internal helper returning the indices of the five `nn_linear` modules within
#' the encoder built by [tape_create_model()].
#'
#' @return An integer vector of length five.
#'
#' @keywords internal
#' @noRd
tape_encoder_linear_index <- function() {
  c(2L, 5L, 8L, 11L, 13L)
}
 
 
#' Draw fixed dropout masks for the TAPE adaptive stage
#'
#' Internal helper that draws the four encoder dropout masks a single time.
#'
#' @param dims Integer vector of length four giving the feature dimension seen by
#'   each encoder dropout layer, in forward order.
#' @param n_rows Integer scalar. Number of rows the masks are drawn for. Use the
#'   number of samples in the batch for overall adaptation, or `1` when masks are
#'   broadcast across an ensemble of single-sample models.
#' @param dev Character scalar. Torch device the masks are allocated on.
#' @param seed Integer scalar. Seed used for the draw.
#'
#' @return A list of four `torch_tensor` multiplicative masks.
#'
#' @details
#' The reference implementation reseeds both the CPU and CUDA generators before
#' every optimization step of its adaptive stage, so the same four masks are
#' drawn at each of the `2 * step * max_iter` steps. Drawing them once is an
#' exact restatement of that behaviour rather than an approximation, and it
#' removes any dependence on whether a given torch build propagates a manual seed
#' to the CUDA generator.
#'
#' @importFrom torch torch_manual_seed torch_ones nnf_dropout
#' @keywords internal
#' @noRd
tape_make_masks <- function(dims, n_rows, dev, seed = 0L) {
  torch_manual_seed(seed)
 
  lapply(dims, function(d) {
    nnf_dropout(
      torch_ones(c(as.integer(n_rows), as.integer(d)), device = dev),
      p = 0.5,
      training = TRUE
    )
  })
}
 
 
#' Encode with externally supplied dropout masks
#'
#' Internal helper that runs the TAPE encoder while substituting precomputed
#' multiplicative masks for its `nn_dropout` modules.
#'
#' @param encoder The encoder module of a TAPE model.
#' @param x A `torch_tensor` with samples in rows and genes in columns.
#' @param masks A list of four masks as returned by [tape_make_masks()], or
#'   `NULL` to encode without dropout.
#'
#' @return A `torch_tensor` of latent cell-type scores.
#'
#' @importFrom torch nnf_celu
#' @keywords internal
#' @noRd
tape_masked_encode <- function(encoder, x, masks = NULL) {
  linears <- tape_encoder_linear_index()
  h <- x
 
  for (l in 1:4) {
    if (!is.null(masks)) {
      h <- h * masks[[l]]
    }
    h <- nnf_celu(encoder[[linears[l]]](h))
  }
 
  encoder[[linears[5]]](h)
}

#' Fused Adam optimizer for TAPE
#'
#' Minimal Adam optimizer that performs the entire parameter update in a single
#' fused CUDA kernel via `torch__fused_adam_()`.
#'
#' @param params A `torch_tensor` or a list of `torch_tensor` leaves with
#'   `requires_grad` enabled.
#' @param lr Numeric scalar. Learning rate.
#' @param beta1 Numeric scalar. Exponential decay rate for the first moment.
#' @param beta2 Numeric scalar. Exponential decay rate for the second moment.
#' @param eps Numeric scalar. Term added to the denominator for stability.
#' @param weight_decay Numeric scalar. L2 penalty.
#'
#' @return A list with components `zero_grad` and `step`, both functions of no
#'   arguments, plus `params` for inspection.
#'
#' @details
#' Implements the same update as [torch::optim_adam()] with default arguments,
#' but replaces the R-level loop over parameters and its roughly ten separate
#' elementwise kernels with one fused kernel per call. On the large batched
#' tensors used by [tape_adaptive_chunk()] the optimizer step is bandwidth
#' bound, so the reduction in passes over memory dominates.
#'
#' Step counters are held as device-resident float tensors, one per parameter,
#' and are incremented before the fused call so that bias correction sees a step
#' count of one on the first update.
#'
#' @importFrom torch torch_zeros_like torch_zeros torch_float
#'   torch__foreach_add_ torch__fused_adam_
#' @keywords internal
#' @noRd
tape_fused_adam <- function(params,
                            lr = 1e-4,
                            beta1 = 0.9,
                            beta2 = 0.999,
                            eps = 1e-8,
                            weight_decay = 0) {

  if (inherits(params, "torch_tensor")) {
    params <- list(params)
  }

  params <- unname(params)

  if (length(params) == 0L) {
    stop("`params` is empty.")
  }

  exp_avgs    <- lapply(params, torch_zeros_like)
  exp_avg_sqs <- lapply(params, torch_zeros_like)

  state_steps <- lapply(params, function(p) {
    torch_zeros(1, dtype = torch_float(), device = p$device)
  })

  max_exp_avg_sqs <- list()

  zero_grad <- function() {
    for (p in params) {
      if (!is.null(p$grad)) {
        p$grad$zero_()
      }
    }
    invisible(NULL)
  }

  step <- function() {
    grads <- lapply(params, function(p) p$grad)

    if (any(vapply(grads, is.null, logical(1)))) {
      stop("tape_fused_adam: at least one parameter has no gradient.")
    }

    torch__foreach_add_(state_steps, 1)

    torch__fused_adam_(
      params,
      grads,
      exp_avgs,
      exp_avg_sqs,
      max_exp_avg_sqs,
      state_steps,
      lr,
      beta1,
      beta2,
      weight_decay,
      eps,
      FALSE,
      FALSE
    )

    invisible(NULL)
  }

  list(
    params = params,
    zero_grad = zero_grad,
    step = step
  )
}