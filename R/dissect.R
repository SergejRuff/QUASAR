

#' Extract single-cell counts and metadata
#'
#' Internal helper to extract a count matrix and selected metadata from either
#' a Seurat object or a SingleCellExperiment object.
#'
#' @param sc_data A `Seurat` or `SingleCellExperiment` object.
#' @param celltype_col Character scalar. Column name in cell metadata containing
#'   cell-type labels.
#' @param batch_col Character scalar or `NULL`. Optional column name in cell
#'   metadata containing batch labels.
#'
#' @return A named list with components:
#' \describe{
#'   \item{expr}{A genes x cells count matrix.}
#'   \item{celltypes}{A character vector of cell-type labels, one per cell.}
#'   \item{batches}{A character vector of batch labels, one per cell, or
#'   `NULL` if no batch column was supplied.}
#' }
#'
#' @importFrom SummarizedExperiment assay colData
#' @keywords internal
#' @noRd
dissect_extract_sc_data <- function(sc_data,
                                    celltype_col = "celltype",
                                    batch_col = NULL) {
  if (inherits(sc_data, "Seurat")) {
    expr <- as.matrix(sc_data[["RNA"]]$counts)
    meta <- sc_data@meta.data
  } else if (inherits(sc_data, "SingleCellExperiment")) {
    expr <- as.matrix(SummarizedExperiment::assay(sc_data, "counts"))
    meta <- as.data.frame(SummarizedExperiment::colData(sc_data))
  } else {
    stop("sc_data must be a Seurat or SingleCellExperiment object")
  }

  if (!(celltype_col %in% colnames(meta))) {
    stop("celltype_col not found in single-cell metadata")
  }

  celltypes <- as.character(meta[[celltype_col]])
  celltypes <- gsub("/", "_", celltypes)

  batches <- NULL
  if (!is.null(batch_col) && batch_col %in% colnames(meta)) {
    batches <- as.character(meta[[batch_col]])
  }

  list(
    expr = expr,
    celltypes = celltypes,
    batches = batches
  )
}















#' Per-batch log1p-MinMax normalisation
#'
#' Internal helper implementing the same per-batch normalisation used in the
#' fraction model: log1p transformation followed by feature-wise MinMax scaling
#' within each sample.
#'
#' @param x A `torch_tensor` with shape batch x features.
#' @param epsilon Numeric scalar. Small constant added to the denominator to
#'   avoid division by zero.
#'
#' @return A `torch_tensor` with the same shape as `x`.
#'
#' @importFrom torch torch_log1p
#' @keywords internal
#' @noRd
dissect_normalize_per_batch_torch <- function(x, epsilon = 1e-8) {
  x1 <- torch_log1p(x) / log(2)
  min_vals <- x1$min(dim = 2, keepdim = TRUE)[[1]]
  max_vals <- x1$max(dim = 2, keepdim = TRUE)[[1]]
  (x1 - min_vals) / (max_vals - min_vals + epsilon)
}

















#' Aggregate duplicated row names
#'
#' Internal helper for resolving duplicated gene names in matrices with genes in
#' rows.
#'
#' @param x A matrix with genes in rows.
#' @param duplicated Character scalar specifying how duplicated row names should
#'   be handled. One of `"first"`, `"sum"`, or `"mean"`.
#'
#' @return A matrix with unique row names.
#'
#' @keywords internal
#' @noRd
dissect_aggregate_rows_by_name <- function(x, duplicated = "first") {
  rn <- rownames(x)
  if (is.null(rn) || !anyDuplicated(rn)) {
    return(x)
  }

  if (duplicated == "first") {
    return(x[!duplicated(rn), , drop = FALSE])
  }

  summed <- rowsum(x, group = rn, reorder = FALSE)
  if (duplicated == "sum") {
    return(summed)
  }
  if (duplicated == "mean") {
    counts <- as.numeric(table(rn)[rownames(summed)])
    return(summed / counts)
  }

  stop("duplicated must be one of: 'first', 'sum', 'mean'")
}




















#' Simulate DISSECT training mixtures from single-cell data
#'
#' Generates simulated bulk or spatial transcriptomics mixtures from a
#' single-cell reference, following the simulation logic used in the original
#' DISSECT workflow. Returned proportions are stored as samples x cell types.
#'
#' @param sc_data A `Seurat` or `SingleCellExperiment` object containing
#'   single-cell counts.
#' @param celltype_col Character scalar. Column name in cell metadata containing
#'   cell-type labels.
#' @param batch_col Character scalar or `NULL`. Optional column name in cell
#'   metadata containing batch labels. If multiple batches are present,
#'   simulation is performed per batch and concatenated.
#' @param type Character scalar. Either `"bulk"` or `"st"`.
#' @param n_samples Integer scalar or `NULL`. Number of simulated samples. If
#'   `NULL`, defaults to `1000 * n_celltypes`.
#' @param cells_per_sample Integer scalar. Number of cells per simulated sample
#'   for bulk simulation.
#' @param prop_sparse Numeric scalar in `[0, 1]`. Fraction of sparse simulated
#'   bulk samples.
#' @param concentration Numeric vector or `NULL`. Dirichlet concentration
#'   parameter for bulk simulation. If `NULL`, a vector of ones is used.
#' @param save_expr Logical scalar. If `TRUE`, stores per-cell-type expression
#'   contributions for each simulated sample.
#' @param min_genes Integer scalar. Minimum number of detected genes required
#'   for a cell to be kept during preprocessing.
#' @param min_cells Integer scalar. Minimum number of cells required for a gene
#'   to be kept during preprocessing.
#' @param mt_cutoff Numeric scalar. Maximum allowed mitochondrial percentage per
#'   cell.
#' @param min_expr Numeric scalar. Minimum mean `log1p` expression threshold
#' @param downsample Numeric scalar or `NULL`. Optional downsampling factor for
#'   spatial transcriptomics simulation.
#' @param seed Integer scalar. Random seed used for simulation.
#'
#' @return A named list with components:
#' \describe{
#'   \item{X}{A simulated expression matrix with samples in rows and genes in
#'   columns.}
#'   \item{props}{A proportions matrix with samples in rows and cell types in
#'   columns.}
#'   \item{cells}{An integer matrix of sampled cell counts with samples in rows
#'   and cell types in columns.}
#'   \item{gene_names}{Character vector of gene names.}
#'   \item{celltype_names}{Character vector of cell-type names.}
#'   \item{layers}{Named list of per-cell-type expression matrices, or `NULL`
#'   if `save_expr = FALSE`.}
#'   \item{batch}{Character vector giving the simulated batch assignment of each
#'   sample.}
#'   \item{type}{Simulation type, either `"bulk"` or `"st"`.}
#' }
#'
#' @details
#' The preprocessing and simulation logic are designed to follow the original
#' DISSECT implementation as closely as possible while using in-memory R data
#' structures instead of file-based Python objects.
#' 
#' 
#' @references
#' Khatri, R., Machart, P., & Bonn, S. (2024). DISSECT: deep semi-supervised
#' consistency regularization for accurate cell type fraction and gene
#' expression estimation. \emph{Genome Biology}, 25(1), 112.
#'
#' Original DISSECT software repository:
#' \url{https://github.com/imsb-uke/DISSECT}
#'
#' @importFrom stats median rgamma rhyper setNames
#' @examples
#' \dontrun{
#' sim <- dissect_simulate(
#'   sc_data = sce,
#'   celltype_col = "celltype",
#'   batch_col = "batch",
#'   type = "bulk"
#' )
#' }
#'
#' @export
dissect_simulate <- function(sc_data,
                             celltype_col = "celltype",
                             batch_col = NULL,
                             type = "bulk",
                             n_samples = NULL,
                             cells_per_sample = 500,
                             prop_sparse = 0.5,
                             concentration = NULL,
                             save_expr = TRUE,
                             min_genes = 200,
                             min_cells = 3,
                             mt_cutoff = 5,
                             min_expr = 0,
                             downsample = NULL,
                             seed = 42) {
  set.seed(seed)

  if (!(type %in% c("bulk", "st"))) {
    stop("type must be 'bulk' or 'st'")
  }

  sc <- dissect_extract_sc_data(
    sc_data = sc_data,
    celltype_col = celltype_col,
    batch_col = batch_col
  )

  expr <- sc$expr
  celltypes <- sc$celltypes
  batches <- sc$batches

  if (is.null(rownames(expr))) {
    rownames(expr) <- paste0("gene_", seq_len(nrow(expr)))
  }
  rownames(expr) <- make.unique(rownames(expr))

  if (!inherits(expr, "sparseMatrix")) {
    expr <- Matrix::Matrix(expr, sparse = TRUE)
  }

  keep_cells <- Matrix::colSums(expr > 0) >= min_genes
  expr <- expr[, keep_cells, drop = FALSE]
  celltypes <- celltypes[keep_cells]
  if (!is.null(batches)) {
    batches <- batches[keep_cells]
  }

  keep_genes <- Matrix::rowSums(expr > 0) >= min_cells
  expr <- expr[keep_genes, , drop = FALSE]

  mt_mask <- grepl("^MT-", rownames(expr))
  if (any(mt_mask)) {
    lib <- Matrix::colSums(expr)
    lib[lib == 0] <- 1
    pct_mt <- Matrix::colSums(expr[mt_mask, , drop = FALSE]) / lib * 100
  } else {
    pct_mt <- rep(0, ncol(expr))
  }

  keep_mt <- pct_mt < mt_cutoff
  expr <- expr[, keep_mt, drop = FALSE]
  celltypes <- celltypes[keep_mt]
  if (!is.null(batches)) {
    batches <- batches[keep_mt]
  }

  lib <- Matrix::colSums(expr)
  lib_nonzero <- lib[lib > 0]
  median_lib <- if (length(lib_nonzero) > 0) stats::median(lib_nonzero) else 1
  lib_safe <- lib
  lib_safe[lib_safe == 0] <- 1
  expr <- expr %*% Matrix::Diagonal(x = median_lib / lib_safe)

  gene_means <- log1p(Matrix::rowMeans(expr))
  expr <- expr[gene_means > min_expr, , drop = FALSE]

  expr_cg <- Matrix::t(expr)
  gene_names <- colnames(expr_cg)
  if (is.null(gene_names)) {
    gene_names <- rownames(expr)
  }

  celltype_names <- sort(unique(celltypes), method = "radix")
  n_celltypes <- length(celltype_names)

  if (is.null(n_samples)) {
    n_samples <- 1000L * n_celltypes
  }

  if (type == "bulk") {
    if (is.null(concentration)) {
      concentration <- rep(1, n_celltypes)
    }
    if (length(concentration) != n_celltypes) {
      stop("concentration must have length equal to the number of cell types")
    }

    n_sparse <- as.integer(n_samples * prop_sparse)
    n_complete <- n_samples - n_sparse
    min_prc <- 1 / cells_per_sample

    props_complete <- matrix(0, n_complete, n_celltypes)
    if (n_complete > 0) {
      draws_complete <- matrix(
        rgamma(n_complete * n_celltypes,
               shape = rep(concentration, each = n_complete),
               rate = 1),
        nrow = n_complete,
        ncol = n_celltypes,
        byrow = FALSE
      )
      draws_complete <- draws_complete / rowSums(draws_complete)
      draws_complete[draws_complete < min_prc] <- 0

      zero_rows <- which(rowSums(draws_complete) == 0)
      if (length(zero_rows) > 0) {
        max_ct <- which.max(concentration)
        draws_complete[cbind(zero_rows, max_ct)] <- 1
      }

      props_complete <- draws_complete / rowSums(draws_complete)
    }

    cells_complete <- round(props_complete * cells_per_sample)
    if (n_complete > 0) {
      zero_rows <- which(rowSums(cells_complete) == 0)
      if (length(zero_rows) > 0) {
        max_idx <- max.col(props_complete[zero_rows, , drop = FALSE], ties.method = "first")
        cells_complete[cbind(zero_rows, max_idx)] <- 1L
      }
    }
    props_complete <- cells_complete / rowSums(cells_complete)

    props_sparse <- matrix(0, n_sparse, n_celltypes)
    if (n_sparse > 0) {
      for (i in seq_len(n_sparse)) {
        alpha <- rep(1, n_celltypes)
        no_keep <- sample.int(n_celltypes - 1L, 1L)
        no_keeps <- sample.int(n_celltypes, no_keep, replace = FALSE)
        alpha[no_keeps] <- 1e-6

        draw <- rgamma(n_celltypes, shape = alpha, rate = 1)
        draw <- draw / sum(draw)
        draw[draw < min_prc] <- 0

        if (sum(draw) == 0) {
          draw[which.max(alpha)] <- 1
        }

        props_sparse[i, ] <- draw / sum(draw)
      }
    }

    cells_sparse <- round(props_sparse * cells_per_sample)
    if (n_sparse > 0) {
      zero_rows <- which(rowSums(cells_sparse) == 0)
      if (length(zero_rows) > 0) {
        max_idx <- max.col(props_sparse[zero_rows, , drop = FALSE], ties.method = "first")
        cells_sparse[cbind(zero_rows, max_idx)] <- 1L
      }
    }
    props_sparse <- cells_sparse / rowSums(cells_sparse)

    props <- rbind(props_complete, props_sparse)
    cells <- rbind(cells_complete, cells_sparse)
  } else {
    props <- matrix(0, n_samples, n_celltypes)
    cells <- matrix(0L, n_samples, n_celltypes)

    for (i in seq_len(n_samples)) {
      keep_sparse <- rep(1e-6, n_celltypes)
      n_keep <- sample.int(min(5L, n_celltypes), 1L)
      keeps <- sample.int(n_celltypes, n_keep, replace = FALSE)
      keep_sparse[keeps] <- 1

      draw <- rgamma(n_celltypes, shape = keep_sparse, rate = 1)
      draw <- draw / sum(draw)
      props[i, ] <- draw

      n_cells_spot <- sample(5:11, 1)
      cells[i, ] <- as.integer(round(draw * n_cells_spot))

      if (sum(cells[i, ]) == 0) {
        cells[i, which.max(draw)] <- 1L
      }
    }

    props <- cells / rowSums(cells)
  }

  colnames(props) <- celltype_names
  colnames(cells) <- celltype_names

  if (is.null(batches) || length(unique(batches)) == 1) {
    batch_levels <- NULL
  } else {
    batch_levels <- unique(batches)
  }

  if (is.null(batch_levels)) {
    batch_expr_list <- list(expr_cg)
    batch_celltype_list <- list(celltypes)
    batch_names <- "batch1"
  } else {
    batch_expr_list <- lapply(batch_levels, function(b) expr_cg[batches == b, , drop = FALSE])
    batch_celltype_list <- lapply(batch_levels, function(b) celltypes[batches == b])
    batch_names <- batch_levels
  }

  X_batches <- vector("list", length(batch_expr_list))
  layer_batches <- if (save_expr) vector("list", length(batch_expr_list)) else NULL

  for (b in seq_along(batch_expr_list)) {
    expr_now <- batch_expr_list[[b]]
    ct_now <- batch_celltype_list[[b]]

    if (!inherits(expr_now, "sparseMatrix")) {
      expr_now <- Matrix::Matrix(expr_now, sparse = TRUE)
    }

    n_sim <- nrow(cells)
    n_genes <- ncol(expr_now)

    ct_indices <- lapply(celltype_names, function(ct) which(ct_now == ct))
    names(ct_indices) <- celltype_names
    ct_sizes <- vapply(ct_indices, length, integer(1))

    expr_by_ct <- vector("list", n_celltypes)
    for (j in seq_len(n_celltypes)) {
      if (ct_sizes[j] > 0) {
        expr_by_ct[[j]] <- expr_now[ct_indices[[j]], , drop = FALSE]
      } else {
        expr_by_ct[[j]] <- NULL
      }
    }

    X_out <- matrix(0, nrow = n_sim, ncol = n_genes)

    if (save_expr) {
      layer_out <- setNames(vector("list", length(celltype_names)), celltype_names)
      for (ct in celltype_names) {
        layer_out[[ct]] <- matrix(0, nrow = n_sim, ncol = n_genes)
      }
    } else {
      layer_out <- NULL
    }

    for (j in seq_len(n_celltypes)) {
      pool_size <- ct_sizes[j]
      if (pool_size == 0L) next

      counts_j <- cells[, j]

      if (type == "bulk") {
        counts_eff <- counts_j
      } else {
        counts_eff <- pmin(counts_j, pool_size)
      }

      total_draws <- sum(counts_eff)
      if (total_draws <= 0L) next

      sample_ids <- rep.int(seq_len(n_sim), counts_eff)
      local_ids <- sample.int(pool_size, size = total_draws, replace = TRUE)

      W_j <- Matrix::sparseMatrix(
        i = local_ids,
        j = sample_ids,
        x = 1,
        dims = c(pool_size, n_sim)
      )

      contrib_j <- Matrix::t(Matrix::t(expr_by_ct[[j]]) %*% W_j)
      contrib_j <- as.matrix(contrib_j)

      X_out <- X_out + contrib_j

      if (save_expr) {
        if (type == "bulk") {
          layer_out[[celltype_names[j]]] <- contrib_j
        } else {
          denom <- counts_eff
          denom_safe <- denom
          denom_safe[denom_safe == 0] <- 1
          layer_j <- sweep(contrib_j, 1, denom_safe, "/")
          layer_j[denom == 0, ] <- 0
          layer_out[[celltype_names[j]]] <- layer_j
        }
      }
    }

    X_batches[[b]] <- X_out
    if (save_expr) {
      layer_batches[[b]] <- layer_out
    }
  }

  X <- do.call(rbind, X_batches)
  props_out <- do.call(rbind, replicate(length(X_batches), props, simplify = FALSE))
  cells_out <- do.call(rbind, replicate(length(X_batches), cells, simplify = FALSE))
  batch_out <- rep(batch_names, each = nrow(props))

  rownames(X) <- paste0(batch_out, "_sample_", sequence(tabulate(match(batch_out, batch_out))))
  colnames(X) <- gene_names
  rownames(props_out) <- rownames(X)
  rownames(cells_out) <- rownames(X)

  layers <- NULL
  if (save_expr) {
    layers <- setNames(vector("list", length(celltype_names)), celltype_names)
    for (ct in celltype_names) {
      layers[[ct]] <- do.call(rbind, lapply(layer_batches, function(x) x[[ct]]))
      rownames(layers[[ct]]) <- rownames(X)
      colnames(layers[[ct]]) <- gene_names
    }
  }

  if (!is.null(downsample) && type == "st") {
    for (i in seq_len(nrow(X))) {
      lib_i <- sum(X[i, ])
      if (lib_i <= 0) next
      target <- as.integer(round(downsample * lib_i))
      if (target <= 0) { X[i, ] <- 0; next }
      if (target >= lib_i) next
      # sc.pp.downsample_counts defaults to replace = FALSE: reads are drawn
      # without replacement (multivariate hypergeometric), so a gene can never
      # gain reads. rmultinom() samples with replacement.
      x <- as.integer(X[i, ])
      remaining <- sum(x)
      drawn <- target
      out <- integer(length(x))
      for (g in seq_along(x)) {
        if (drawn <= 0L) break
        if (x[g] > 0L) {
          out[g] <- stats::rhyper(1, m = x[g], n = remaining - x[g], k = drawn)
          drawn <- drawn - out[g]
        }
        remaining <- remaining - x[g]
      }
      X[i, ] <- out
    }
  }

  list(
    X = X,
    props = props_out,
    cells = cells_out,
    gene_names = gene_names,
    celltype_names = celltype_names,
    layers = layers,
    batch = batch_out,
    type = type
  )
}














#' Prepare DISSECT input matrices
#'
#' Prepares real and simulated matrices for DISSECT fraction estimation by
#' applying the same preprocessing order used in the original workflow:
#' transformation, variance filtering, deduplication, normalisation, gene
#' intersection, and size balancing.
#'
#' @param bulk Numeric matrix with genes in rows and samples in columns.
#' @param reference Numeric matrix with genes in rows and cell types in columns,
#'   or `NULL`. This is a convenience path and is only used when `sim_data` is
#'   not supplied.
#' @param sim_data A list returned by [dissect_simulate()], or `NULL`.
#' @param test_dataset_type Character scalar. Either `"bulk"` or
#'   `"microarray"`.
#' @param duplicated Character scalar. How duplicated gene names should be
#'   resolved. One of `"first"`, `"sum"`, or `"mean"`.
#' @param normalize_simulated Character scalar or `NULL`. Currently `"cpm"` or
#'   `NULL`.
#' @param normalize_test Character scalar or `NULL`. Currently `"cpm"` or
#'   `NULL`.
#' @param var_cutoff Numeric scalar or `NULL`. Variance threshold applied to the
#'   bulk input before transposition.
#' @param test_in_mix Integer scalar. Number of real samples used in the online
#'   mixing step.
#'
#' @return A named list with components:
#' \describe{
#'   \item{X_real_train}{Real training matrix with samples in rows and genes in
#'   columns.}
#'   \item{X_sim}{Simulated matrix with samples in rows and genes in columns.}
#'   \item{y_sim}{Simulated proportions with samples in rows and cell types in
#'   columns.}
#'   \item{X_real_test}{Real test matrix with samples in rows and genes in
#'   columns.}
#'   \item{sample_names}{Character vector of sample names.}
#'   \item{celltypes}{Character vector of cell-type names.}
#'   \item{genes}{Character vector of common genes used for modelling.}
#'   \item{reference}{Reference matrix if supplied, otherwise `NULL`.}
#'   \item{sim_data}{Simulation object if supplied, otherwise `NULL`.}
#' }
#'
#' @details
#' The intended DISSECT workflow uses `sim_data` generated from
#' [dissect_simulate()]. A direct `reference` matrix is supported as a
#' convenience interface for proportion estimation.
#' 
#' 
#' @references
#' Khatri, R., Machart, P., & Bonn, S. (2024). DISSECT: deep semi-supervised
#' consistency regularization for accurate cell type fraction and gene
#' expression estimation. \emph{Genome Biology}, 25(1), 112.
#'
#' Original DISSECT software repository:
#' \url{https://github.com/imsb-uke/DISSECT}
#' @importFrom stats var
#' @examples
#' \dontrun{
#' proc <- dissect_process(
#'   bulk = bulk_mat,
#'   sim_data = sim
#' )
#' }
#'
#' @export
dissect_process <- function(bulk,
                            reference = NULL,
                            sim_data = NULL,
                            test_dataset_type = "bulk",
                            duplicated = "first",
                            normalize_simulated = "cpm",
                            normalize_test = "cpm",
                            var_cutoff = 0.1,
                            test_in_mix = 1) {
  if (is.null(rownames(bulk))) {
    rownames(bulk) <- paste0("gene_", seq_len(nrow(bulk)))
  }
  if (is.null(colnames(bulk))) {
    colnames(bulk) <- paste0("sample_", seq_len(ncol(bulk)))
  }

  bulk <- as.matrix(bulk)

  if (test_dataset_type == "microarray") {
    bulk <- 2^bulk - 1
    bulk[bulk < 0] <- 0
  }

  if (!is.null(var_cutoff)) {
    keep_var <- apply(bulk, 1, var) > var_cutoff
    bulk <- bulk[keep_var, , drop = FALSE]
  }

  X_real <- t(bulk)
  sample_names <- rownames(X_real)

  if (anyDuplicated(colnames(X_real))) {
    x_real_gs <- t(X_real)
    x_real_gs <- dissect_aggregate_rows_by_name(x_real_gs, duplicated = duplicated)
    X_real <- t(x_real_gs)
  }

  if (!is.null(normalize_test) && normalize_test == "cpm") {
    lib <- rowSums(X_real)
    lib[lib == 0] <- 1
    X_real <- sweep(X_real, 1, lib, "/")
    X_real <- sweep(X_real, 1, rep(1e6, nrow(X_real)), "*")
  }

  if (!is.null(sim_data)) {
    X_sim <- sim_data$X
    y_sim <- sim_data$props
  } else {
    if (is.null(reference)) {
      stop("Provide either sim_data or reference")
    }
    if (is.null(rownames(reference))) {
      rownames(reference) <- paste0("gene_", seq_len(nrow(reference)))
    }
    if (is.null(colnames(reference))) {
      colnames(reference) <- paste0("celltype_", seq_len(ncol(reference)))
    }

    reference <- dissect_aggregate_rows_by_name(as.matrix(reference), duplicated = duplicated)

    X_sim <- t(reference)
    y_sim <- diag(ncol(reference))
    colnames(y_sim) <- colnames(reference)
    rownames(y_sim) <- colnames(reference)
  }

  if (!is.null(normalize_simulated) && normalize_simulated == "cpm") {
    lib <- rowSums(X_sim)
    lib[lib == 0] <- 1
    X_sim <- sweep(X_sim, 1, lib, "/")
    X_sim <- sweep(X_sim, 1, rep(1e6, nrow(X_sim)), "*")
  }

  if (anyDuplicated(colnames(X_sim))) {
    x_sim_gs <- t(X_sim)
    x_sim_gs <- dissect_aggregate_rows_by_name(x_sim_gs, duplicated = duplicated)
    X_sim <- t(x_sim_gs)
  }

  y_sim <- as.matrix(y_sim)
  y_sim[y_sim < 0.005] <- 0
  y_sum <- rowSums(y_sim)
  y_sum[y_sum == 0] <- 1
  y_sim <- y_sim / y_sum

  genes_common <- intersect(colnames(X_real), colnames(X_sim))
  X_real <- X_real[, genes_common, drop = FALSE]
  X_sim <- X_sim[, genes_common, drop = FALSE]
  X_real_test <- X_real

  if (!is.null(test_in_mix) && test_in_mix > 0) {
    n_mix <- min(test_in_mix, nrow(X_real))
    X_real_train <- X_real[seq_len(n_mix), , drop = FALSE]
  } else {
    X_real_train <- X_real
  }

  real_size <- nrow(X_real_train)
  sim_size <- nrow(X_sim)

  if (real_size < sim_size) {
    idx <- rep(seq_len(real_size), length.out = sim_size)
    X_real_train <- X_real_train[idx, , drop = FALSE]
  }

  real_size <- nrow(X_real_train)
  if (sim_size < real_size) {
    idx <- rep(seq_len(sim_size), length.out = real_size)
    X_sim <- X_sim[idx, , drop = FALSE]
    y_sim <- y_sim[idx, , drop = FALSE]
  }

  list(
    X_real_train = X_real_train,
    X_sim = X_sim,
    y_sim = y_sim,
    X_real_test = X_real_test,
    sample_names = sample_names,
    celltypes = colnames(y_sim),
    genes = genes_common,
    reference = reference,
    sim_data = sim_data
  )
}


















#' Resolve DISSECT compute device
#'
#' Internal helper to select CPU or CUDA for torch-based training and inference.
#'
#' @param device Character scalar. One of `"auto"`, `"cpu"`, or `"cuda"`.
#' @param cuda_index Integer scalar or `NULL`. Optional CUDA device index.
#' @param verbose Logical scalar. If `TRUE`, prints the selected device.
#' @param stage Character scalar or `NULL`. Optional label for status output.
#'
#' @return A list with:
#' \describe{
#'   \item{device}{A `torch_device` object.}
#'   \item{device_name}{A character label such as `"cpu"` or `"cuda:0"`.}
#'   \item{using_cuda}{Logical scalar indicating whether CUDA is used.}
#' }
#'
#' @importFrom torch cuda_is_available cuda_device_count torch_device
#' @keywords internal
#' @noRd
dissect_resolve_device <- function(device = c("auto", "cpu", "cuda"),
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
      if (length(cuda_index) != 1L || is.na(cuda_index) || cuda_index < 0L || cuda_index >= n_cuda) {
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


















#' Estimate cell-type proportions with DISSECT
#'
#' Trains the DISSECT fraction model and returns estimated cell-type
#' proportions and pre-activation scores for each sample.
#'
#' @param bulk Numeric matrix with genes in rows and samples in columns, or
#'   `NULL` if `processed` is supplied.
#' @param reference Numeric matrix with genes in rows and cell types in
#'   columns, or `NULL`.
#' @param sim_data A list returned by [dissect_simulate()], or `NULL`.
#' @param processed A list returned by [dissect_process()], or `NULL`.
#' @param test_dataset_type Character scalar. Either `"bulk"` or a supported
#'   alternative matching the original workflow.
#' @param duplicated Character scalar controlling duplicated gene handling.
#' @param normalize_simulated Character scalar or `NULL`. Currently `"cpm"` or
#'   `NULL`.
#' @param normalize_test Character scalar or `NULL`. Currently `"cpm"` or
#'   `NULL`.
#' @param var_cutoff Numeric scalar or `NULL`. Variance threshold for bulk
#'   preprocessing.
#' @param test_in_mix Integer scalar. Number of real samples to use in online
#'   mixing.
#' @param n_hidden_layers Integer scalar. Number of hidden layers in the
#'   fraction model.
#' @param hidden_units Integer vector. Number of units in each hidden layer.
#' @param hidden_activation Character scalar. Hidden-layer activation, applied
#'   to every ensemble model. Note that the reference Python mutates its config
#'   in place inside the model loop, so only the first model receives `relu6`
#'   and the remaining four fall back to `relu`; this implementation follows the
#'   documented behaviour instead.
#' @param output_activation Character scalar. Output activation function.
#' @param loss Character scalar. One of `"kldivergence"`, `"l2"`, or `"l1"`.
#' @param n_steps Integer scalar. Number of training steps.
#' @param lr Numeric scalar. Learning rate.
#' @param batch_size Integer scalar. Batch size.
#' @param dropout Numeric vector or `NULL`. Dropout rates for hidden layers.
#' @param alpha_range Numeric vector of length two. Range for the mixing
#'   coefficient used in online mixtures.
#' @param normalization_per_batch Character scalar or `NULL`. Currently
#'   `"log1p-MinMax"` or `NULL`.
#' @param models Integer vector. Ensemble model identifiers. As in the original
#'   Python code, only the number of models is used.
#' @param mix Character scalar. Mixing strategy, either `"srm"` or `"rrm"`.
#' @param device Character scalar. One of `"auto"`, `"cpu"`, or `"cuda"`.
#' @param cuda_index Integer scalar or `NULL`. Optional CUDA device index.
#'
#' @return A named list with components:
#' \describe{
#'   \item{fractions}{A data frame of estimated cell-type proportions with
#'   samples in rows and cell types in columns.}
#'   \item{scores}{A data frame of pre-activation output scores with samples in
#'   rows and cell types in columns.}
#'   \item{processed}{The processed input object used for training.}
#' }
#'
#' @details
#' dissect_prop implements the DISSECT cell-type fraction estimation
#' strategy described by Khatri, Machart, and Bonn (2024) in torch in R.
#' It reproduces the semi-supervised consistency-regularized training
#' workflow for fraction estimation within an R interface.
#'
#' Progress bars and elapsed time are reported during training.
#' 
#' @references
#' Khatri, R., Machart, P., & Bonn, S. (2024). DISSECT: deep semi-supervised
#' consistency regularization for accurate cell type fraction and gene
#' expression estimation. \emph{Genome Biology}, 25(1), 112.
#'
#' Original DISSECT software repository:
#' \url{https://github.com/imsb-uke/DISSECT}
#' @importFrom stats runif
#' @importFrom torch torch_tensor torch_float32 nn_module nn_module_list nn_linear nn_dropout nnf_relu nnf_softmax torch_sigmoid torch_tanh torch_clamp optim_adam with_no_grad torch_mean torch_sum torch_log torch_abs torch_manual_seed  cuda_synchronize torch_device
#' @examples
#' \dontrun{
#' prop <- dissect_prop(
#'   processed = proc,
#'   n_steps = 5000,
#'   device = "auto"
#' )
#' }
#'
#' @export
dissect_prop <- function(bulk = NULL,
                         reference = NULL,
                         sim_data = NULL,
                         processed = NULL,
                         test_dataset_type = "bulk",
                         duplicated = "first",
                         normalize_simulated = "cpm",
                         normalize_test = "cpm",
                         var_cutoff = 0.1,
                         test_in_mix = 1,
                         n_hidden_layers = 4,
                         hidden_units = c(512, 256, 128, 64),
                         hidden_activation = "relu6",
                         output_activation = "softmax",
                         loss = "kldivergence",
                         n_steps = 5000,
                         lr = 1e-5,
                         batch_size = 64,
                         dropout = NULL,
                         alpha_range = c(0.1, 0.9),
                         normalization_per_batch = "log1p-MinMax",
                         models = c(1, 2, 3, 4, 5),
                         mix = "srm",
                         device = c("auto", "cpu", "cuda"),
                         cuda_index = NULL) {

  self <- NULL

  format_hms <- function(seconds) {
    seconds <- max(0, as.integer(round(seconds)))
    hh <- seconds %/% 3600
    mm <- (seconds %% 3600) %/% 60
    ss <- seconds %% 60
    sprintf("%02d:%02d:%02d", hh, mm, ss)
  }

  make_bar <- function(i, n, width = 30L) {
    if (n <= 0) {
      return(paste(rep("-", width), collapse = ""))
    }
    filled <- as.integer(floor(width * i / n))
    paste0(
      paste(rep("=", filled), collapse = ""),
      paste(rep("-", width - filled), collapse = "")
    )
  }

  if (is.null(processed)) {
    if (is.null(bulk)) {
      stop("Provide bulk or processed")
    }
    processed <- dissect_process(
      bulk = bulk,
      reference = reference,
      sim_data = sim_data,
      test_dataset_type = test_dataset_type,
      duplicated = duplicated,
      normalize_simulated = normalize_simulated,
      normalize_test = normalize_test,
      var_cutoff = var_cutoff,
      test_in_mix = test_in_mix
    )
  }

  backend <- dissect_resolve_device(
    device = device,
    cuda_index = cuda_index,
    verbose = TRUE,
    stage = "proportion estimation"
  )
  device_obj <- backend$device

  X_real_train <- processed$X_real_train
  X_sim <- processed$X_sim
  y_sim <- processed$y_sim
  X_real_test <- processed$X_real_test
  sample_names <- processed$sample_names
  celltypes <- processed$celltypes

  use_sig_matrix <- !is.null(processed$reference) && is.null(processed$sim_data)

  n_features <- ncol(X_sim)
  n_celltypes <- ncol(y_sim)
  n_rows <- nrow(X_sim)

  if (n_rows < 1L) {
    stop("No rows available for training")
  }
  if (length(hidden_units) != n_hidden_layers) {
    stop("hidden_units must have length n_hidden_layers")
  }
  if (!is.null(dropout) && length(dropout) != n_hidden_layers) {
    stop("dropout must have length n_hidden_layers when provided")
  }
  if (nrow(X_real_train) != n_rows) {
    stop("X_real_train and X_sim must have the same number of rows after processing")
  }
  if (!(mix %in% c("srm", "rrm"))) {
    stop("mix must be one of 'srm' or 'rrm'")
  }
  if (!(loss %in% c("kldivergence", "l2", "l1"))) {
    stop("loss must be one of 'kldivergence', 'l2', 'l1'")
  }
  if (!is.null(normalization_per_batch) &&
      normalization_per_batch != "log1p-MinMax") {
    stop("Only 'log1p-MinMax' or NULL are supported for normalization_per_batch")
  }

  X_real_train_t <- torch_tensor(X_real_train, dtype = torch_float32(), device = device_obj)
  X_sim_t <- torch_tensor(X_sim, dtype = torch_float32(), device = device_obj)
  y_sim_t <- torch_tensor(y_sim, dtype = torch_float32(), device = device_obj)
  X_real_test_t <- torch_tensor(X_real_test, dtype = torch_float32(), device = device_obj)

  apply_hidden_activation <- function(x) {
    if (identical(hidden_activation, "relu6")) {
      torch_clamp(nnf_relu(x), max = 6)
    } else {
      nnf_relu(x)
    }
  }

  apply_output_activation <- function(x) {
    if (identical(output_activation, "softmax")) {
      nnf_softmax(x, dim = 2)
    } else if (identical(output_activation, "sigmoid")) {
      torch_sigmoid(x)
    } else if (identical(output_activation, "tanh")) {
      torch_tanh(x)
    } else if (identical(output_activation, "relu")) {
      nnf_relu(x)
    } else if (identical(output_activation, "linear") || is.null(output_activation)) {
      x
    } else {
      stop("Unsupported output_activation: ", output_activation)
    }
  }

  do_norm <- function(x) {
    if (is.null(normalization_per_batch)) {
      x
    } else {
      dissect_normalize_per_batch_torch(x)
    }
  }

  DissectNet <- nn_module(
    initialize = function() {
      self$hidden <- nn_module_list()
      self$drop <- nn_module_list()

      in_features <- n_features
      for (i in seq_len(n_hidden_layers)) {
        self$hidden$append(nn_linear(in_features, hidden_units[i]))
        if (is.null(dropout)) {
          self$drop$append(nn_dropout(p = 0))
        } else {
          self$drop$append(nn_dropout(p = dropout[i]))
        }
        in_features <- hidden_units[i]
      }

      self$out <- nn_linear(in_features, n_celltypes)
    },
    forward_hidden = function(x) {
      for (i in seq_len(n_hidden_layers)) {
        x <- self$hidden[[i]](x)
        x <- apply_hidden_activation(x)
        if (!is.null(dropout)) {
          x <- self$drop[[i]](x)
        }
      }
      x
    },
    forward_scores = function(x) {
      x <- self$forward_hidden(x)
      self$out(x)
    },
    forward = function(x) {
      scores <- self$forward_scores(x)
      apply_output_activation(scores)
    }
  )

  ensemble_preds <- vector("list", length(models))
  ensemble_scores <- vector("list", length(models))

  prop_start_time <- Sys.time()
  cat(sprintf("[00:00:00] Running proportion estimation (%d model%s) on %s\n",
              length(models), if (length(models) == 1L) "" else "s", backend$device_name))

  for (m in seq_along(models)) {
    model_seed <- m - 1L
    set.seed(model_seed)
    torch_manual_seed(model_seed)

    model <- DissectNet()
    model <- model$to(device = device_obj)

    optimizer <- optim_adam(model$parameters, lr = lr)
    model$train()

    shuffle_buffer <- min(1000L, n_rows)
    buffer <- seq_len(shuffle_buffer)
    next_source <- if (shuffle_buffer < n_rows) shuffle_buffer + 1L else 1L

    next_stream_index <- function() {
      out <- next_source
      next_source <<- next_source + 1L
      if (next_source > n_rows) {
        next_source <<- 1L
      }
      out
    }

    next_batch_idx <- function(bs) {
      out <- integer(bs)
      for (k in seq_len(bs)) {
        pick <- sample.int(length(buffer), 1L)
        out[k] <- buffer[pick]
        buffer[pick] <<- next_stream_index()
      }
      out
    }

    model_start_time <- Sys.time()
    last_loss <- NA_real_

    cat(sprintf("[00:00:00] Model %d/%d started\n", m, length(models)))

    for (step in seq_len(n_steps)) {
      step_idx <- step - 1L

      batch_idx <- next_batch_idx(batch_size)

      x_sim <- X_sim_t[batch_idx, ]
      y_sim_batch <- y_sim_t[batch_idx, ]
      x_real <- X_real_train_t[batch_idx, ]

      alpha_val <- runif(1, alpha_range[1], alpha_range[2])
      alpha <- torch_tensor(alpha_val, dtype = torch_float32(), device = device_obj)

      if (mix == "rrm") {
        shuf_idx <- sample.int(x_real$size(1), x_real$size(1), replace = FALSE)
        x_real_s <- x_real[shuf_idx, ]
        x_mix <- alpha * x_real + (1 - alpha) * x_real_s
      } else {
        x_mix <- alpha * x_real + (1 - alpha) * x_sim
      }

      x_real_n <- do_norm(x_real)
      x_sim_n <- do_norm(x_sim)
      x_mix_n <- do_norm(x_mix)
      if (mix == "rrm") {
        x_real_s_n <- do_norm(x_real_s)
      }

      if (step_idx == 0L) {
        invisible(model(x_sim_n))
        invisible(model(x_real_n))
        invisible(model(x_mix_n))
      }

      optimizer$zero_grad()

      y_hat_sim <- model(x_sim_n)
      y_hat_real <- model(x_real_n)
      y_hat_mix <- model(x_mix_n)

      if (loss == "kldivergence") {
        reg_loss <- torch_mean(torch_sum(
          y_sim_batch * (torch_log(y_sim_batch + 1e-10) - torch_log(y_hat_sim + 1e-10)),
          dim = 2
        ))
      } else if (loss == "l2") {
        reg_loss <- torch_mean((y_sim_batch - y_hat_sim)^2)
      } else {
        reg_loss <- torch_mean(torch_abs(y_sim_batch - y_hat_sim))
      }

      if (mix == "rrm") {
        y_hat_real_s <- model(x_real_s_n)
        y_mix_target <- alpha * y_hat_real + (1 - alpha) * y_hat_real_s
      } else {
        y_mix_target <- alpha * y_hat_real + (1 - alpha) * y_hat_sim
      }

      if (loss == "l1") {
        cons_loss <- torch_mean(torch_abs(y_mix_target - y_hat_mix))
      } else {
        cons_loss <- torch_mean((y_mix_target - y_hat_mix)^2)
      }

      if (step_idx < 2000L) {
        total_loss <- reg_loss
      } else if (step_idx < 4000L) {
        lambda <- 15
        if (use_sig_matrix || identical(test_dataset_type, "spatial_sparse")) {
          lambda <- 0.15
        }
        total_loss <- reg_loss + lambda * cons_loss
      } else {
        lambda <- 10
        if (use_sig_matrix || identical(test_dataset_type, "spatial_sparse")) {
          lambda <- 0.10
        }
        total_loss <- reg_loss + lambda * cons_loss
      }

      total_loss$backward()
      optimizer$step()

      if (backend$using_cuda) {
        torch::cuda_synchronize()
      }

      last_loss <- as.numeric(total_loss$item())
      elapsed_model <- as.numeric(difftime(Sys.time(), model_start_time, units = "secs"))

      cat(sprintf(
        "\r[%s] Model %d/%d step %d/%d | loss=%.6f | elapsed=%s",
        make_bar(step, n_steps),
        m, length(models),
        step, n_steps,
        last_loss,
        format_hms(elapsed_model)
      ))
      utils::flush.console()
    }

    if (backend$using_cuda) {
      torch::cuda_synchronize()
    }

    model_elapsed <- as.numeric(difftime(Sys.time(), model_start_time, units = "secs"))
    cat(sprintf(
      "\n[%s] Model %d/%d finished | final_loss=%.6f | time=%s\n",
      format_hms(as.numeric(difftime(Sys.time(), prop_start_time, units = "secs"))),
      m, length(models),
      last_loss,
      format_hms(model_elapsed)
    ))

    model$eval()
    with_no_grad({
      x_test <- do_norm(X_real_test_t)
      pred <- model(x_test)
      score <- model$forward_scores(x_test)

      if (backend$using_cuda) {
        torch::cuda_synchronize()
      }

      pred_cpu <- pred$to(device = torch_device("cpu"))
      score_cpu <- score$to(device = torch_device("cpu"))

      ensemble_preds[[m]] <- as.matrix(pred_cpu)
      ensemble_scores[[m]] <- as.matrix(score_cpu)
    })
  }

  pred_ens <- Reduce("+", ensemble_preds) / length(ensemble_preds)
  score_ens <- Reduce("+", ensemble_scores) / length(ensemble_scores)

  score_df <- as.data.frame(score_ens)
  colnames(score_df) <- celltypes
  rownames(score_df) <- sample_names

  pred_df <- as.data.frame(pred_ens)
  colnames(pred_df) <- celltypes
  rownames(pred_df) <- sample_names

  prop_total_elapsed <- as.numeric(difftime(Sys.time(), prop_start_time, units = "secs"))
  cat(sprintf("[%-8s] Proportion estimation done\n", format_hms(prop_total_elapsed)))

  list(
    fractions = pred_df,
    scores = score_df,
    processed = processed
  )
}








#' Estimate cell-type-specific expression with DISSECT
#'
#' Trains the DISSECT expression model and returns estimated expression profiles
#' for each sample-cell-type combination.
#'
#' This implementation keeps the DISSECT expression methodology the same, but
#' avoids materializing the full expanded `(sample x celltype)` training design
#' matrices in memory. Instead, minibatches are assembled on the fly.
#'
#' @param bulk Numeric matrix with genes in rows and samples in columns.
#' @param fractions A data frame or matrix of estimated cell-type fractions with
#'   samples in rows and cell types in columns, typically obtained from
#'   [dissect_prop()].
#' @param sim_data A list returned by [dissect_simulate()] with `save_expr = TRUE`.
#' @param normalize_simulated Character scalar or `NULL`. Currently `"cpm"` or
#'   `NULL`.
#' @param normalize_test Character scalar or `NULL`. Currently `"cpm"` or
#'   `NULL`.
#' @param n_steps_expr Integer scalar or `NULL`. Number of training steps for
#'   the expression model.
#' @param expr_scaling Character scalar. Scaling method used before model
#'   fitting, typically `"p99"` or `"max"`.
#' @param latent_dim Integer scalar. Latent dimension for the conditional VAE.
#' @param batch_size Integer scalar. Batch size.
#' @param lr Numeric scalar. Learning rate.
#' @param beta_vae Numeric scalar. Weight of the KL divergence term.
#' @param lambda_cons Numeric scalar. Weight of the consistency loss term.
#' @param seed Integer scalar or `NULL`. Optional random seed.
#' @param device Character scalar. One of `"auto"`, `"cpu"`, or `"cuda"`.
#' @param cuda_index Integer scalar or `NULL`. Optional CUDA device index.
#'
#' @return A named list with components:
#' \describe{
#'   \item{expression_layered}{Named list of matrices, one per cell type, with
#'   samples in rows and genes in columns.}
#'   \item{expression_combined}{Data frame containing estimated expression for
#'   all sample-cell-type combinations, along with `cell_type` and `sample`
#'   columns.}
#'   \item{scaled_counts}{Scaled expression matrix corresponding to the final
#'   combined estimates.}
#'   \item{celltypes}{Character vector of cell-type names.}
#'   \item{genes}{Character vector of genes used in the model.}
#' }
#'
#' @details
#' dissect_expr implements the DISSECT cell-type-specific expression
#' estimation workflow described by Khatri, Machart, and Bonn (2024) in
#' torch in R. It follows the DISSECT expression-estimation methodology
#' based on simulated per-cell-type layers, estimated sample fractions,
#' and a conditional variational autoencoder. The training logic is unchanged, but
#' expanded input matrices are built batch-wise rather than all at once.
#'
#' Progress bars and elapsed time are reported during training.
#'
#' 
#' @references
#' Khatri, R., Machart, P., & Bonn, S. (2024). DISSECT: deep semi-supervised
#' consistency regularization for accurate cell type fraction and gene
#' expression estimation. \emph{Genome Biology}, 25(1), 112.
#'
#' Original DISSECT software repository:
#' \url{https://github.com/imsb-uke/DISSECT}
#' 
#' @importFrom stats setNames quantile runif
#' @importFrom torch torch_tensor torch_float32 nn_module nn_linear nnf_relu torch_cat torch_randn_like torch_exp torch_mean torch_sum optim_adam with_no_grad torch_manual_seed cuda_synchronize torch_device
#' @examples
#' \dontrun{
#' expr <- dissect_expr(
#'   bulk = bulk_mat,
#'   fractions = prop$fractions,
#'   sim_data = sim,
#'   device = "auto"
#' )
#' }
#'
#' @export
dissect_expr <- function(bulk = NULL,
                         fractions,
                         sim_data,
                         normalize_simulated = "cpm",
                         normalize_test = "cpm",
                         n_steps_expr = 5000,
                         expr_scaling = "p99",
                         latent_dim = 128,
                         batch_size = 128,
                         lr = 1e-3,
                         beta_vae = 0.01,
                         lambda_cons = 0.1,
                         seed = 42,
                         device = c("auto", "cpu", "cuda"),
                         cuda_index = NULL) {

  format_hms <- function(seconds) {
    seconds <- max(0, as.integer(round(seconds)))
    hh <- seconds %/% 3600
    mm <- (seconds %% 3600) %/% 60
    ss <- seconds %% 60
    sprintf("%02d:%02d:%02d", hh, mm, ss)
  }

  self <- NULL

  make_bar <- function(i, n, width = 30L) {
    if (n <= 0) {
      return(paste(rep("-", width), collapse = ""))
    }
    filled <- as.integer(floor(width * i / n))
    paste0(
      paste(rep("=", filled), collapse = ""),
      paste(rep("-", width - filled), collapse = "")
    )
  }

  if (is.null(sim_data$layers)) {
    stop("Expression estimation requires sim_data generated with save_expr = TRUE")
  }
  if (is.null(bulk)) {
    stop("bulk must be provided for expression estimation")
  }

  if (!is.null(seed)) {
    set.seed(seed)
    torch_manual_seed(seed)
  }

  backend <- dissect_resolve_device(
    device = device,
    cuda_index = cuda_index,
    verbose = TRUE,
    stage = "expression estimation"
  )
  device_obj <- backend$device

  expr_total_start <- Sys.time()
  cat("[00:00:00] Preparing expression-estimation inputs\n")

  if (is.null(rownames(bulk))) {
    rownames(bulk) <- paste0("gene_", seq_len(nrow(bulk)))
  }
  if (is.null(colnames(bulk))) {
    colnames(bulk) <- paste0("sample_", seq_len(ncol(bulk)))
  }

  ct_names <- sim_data$celltype_names
  layers <- sim_data$layers

  X_real <- t(as.matrix(bulk))
  X_sim <- as.matrix(sim_data$X)
  props_sim <- as.matrix(sim_data$props)
  fractions <- as.matrix(fractions)

  if (!all(ct_names %in% colnames(fractions))) {
    stop("All simulated cell types must be present as columns in `fractions`.")
  }
  fractions <- fractions[, ct_names, drop = FALSE]

  if (!is.null(rownames(X_real)) && !is.null(rownames(fractions))) {
    if (all(rownames(X_real) %in% rownames(fractions))) {
      fractions <- fractions[rownames(X_real), , drop = FALSE]
    }
  }

  if (!is.null(normalize_test) && normalize_test == "cpm") {
    lib <- rowSums(X_real)
    lib[lib == 0] <- 1
    X_real <- sweep(X_real, 1, lib, "/")
    X_real <- sweep(X_real, 1, rep(1e6, nrow(X_real)), "*")
  }
  X_real <- log1p(X_real)
  X_real <- X_real[, !duplicated(colnames(X_real)), drop = FALSE]

  if (!is.null(normalize_simulated) && normalize_simulated == "cpm") {
    lib <- rowSums(X_sim)
    lib[lib == 0] <- 1
    X_sim <- sweep(X_sim, 1, lib, "/")
    X_sim <- sweep(X_sim, 1, rep(1e6, nrow(X_sim)), "*")
  }
  X_sim <- log1p(X_sim)

  common_genes <- intersect(colnames(X_real), colnames(X_sim))
  if (length(common_genes) == 0L) {
    stop("No common genes between `bulk` and simulated data.")
  }

  X_real <- X_real[, common_genes, drop = FALSE]
  X_sim <- X_sim[, common_genes, drop = FALSE]

  layer_gene_idx <- setNames(vector("list", length(ct_names)), ct_names)
  for (ct in ct_names) {
    idx <- match(common_genes, colnames(layers[[ct]]))
    if (anyNA(idx)) {
      stop("Gene matching failed for layer: ", ct)
    }
    layer_gene_idx[[ct]] <- idx
  }

  layer_libs <- setNames(vector("list", length(ct_names)), ct_names)
  for (ct in ct_names) {
    lib <- rowSums(layers[[ct]])   # full gene set -- layers is not subset here
    lib[lib == 0] <- 1
    layer_libs[[ct]] <- lib
  }

  mean_sum_sim <- mean(rowSums(X_sim))
  real_lib <- rowSums(X_real)
  real_lib[real_lib == 0] <- 1
  X_real_normalized <- X_real / real_lib * mean_sum_sim
  X_sim_normalized <- X_sim

  if (expr_scaling == "p99") {
    max_val <- as.numeric(stats::quantile(
      c(as.vector(X_sim_normalized), as.vector(X_real_normalized)),
      0.99
    ))
  } else {
    max_val <- max(c(X_sim_normalized, X_real_normalized))
  }
  if (max_val == 0) {
    max_val <- 1
  }

  X_sim_scaled <- X_sim_normalized / max_val
  X_real_scaled <- X_real_normalized / max_val

  n_genes <- length(common_genes)
  n_ct <- length(ct_names)
  n_sim <- nrow(X_sim_scaled)
  n_real <- nrow(X_real_scaled)

  hot_encoding <- diag(n_ct)
  colnames(hot_encoding) <- ct_names

  X_sim_scaled_t <- torch_tensor(X_sim_scaled, dtype = torch_float32(), device = device_obj)
  X_real_scaled_t <- torch_tensor(X_real_scaled, dtype = torch_float32(), device = device_obj)

  rm(X_sim_scaled, X_real_scaled, X_sim_normalized, X_real_normalized, X_real, X_sim)
  gc()

  Encoder <- nn_module(
    initialize = function() {
      self$dense_proj2 <- nn_linear(n_genes, 200)
      self$dense_proj3 <- nn_linear(200, 200)
      self$dense_mean <- nn_linear(200, latent_dim)
      self$dense_log_var <- nn_linear(200, latent_dim)
    },
    forward = function(inputs) {
      x <- nnf_relu(self$dense_proj2(inputs))
      x <- nnf_relu(self$dense_proj3(x))
      z_mean <- self$dense_mean(x)
      z_log_var <- self$dense_log_var(x)
      epsilon <- torch_randn_like(z_mean)
      z <- z_mean + torch_exp(0.5 * z_log_var) * epsilon
      list(z_mean = z_mean, z_log_var = z_log_var, z = z)
    }
  )

  Decoder <- nn_module(
    initialize = function() {
      self$dense_proj5 <- nn_linear(latent_dim + n_ct, 200)
      self$dense_proj6 <- nn_linear(200, 200)
      self$dense_output <- nn_linear(200, n_genes)
    },
    forward = function(inputs, labels) {
      x <- torch_cat(list(inputs, labels), dim = 2)
      x <- nnf_relu(self$dense_proj5(x))
      x <- nnf_relu(self$dense_proj6(x))
      nnf_relu(self$dense_output(x))
    }
  )

  VAE <- nn_module(
    initialize = function() {
      self$encoder <- Encoder()
      self$decoder <- Decoder()
    },
    forward = function(inputs, labels) {
      enc <- self$encoder(inputs)
      reconstructed <- self$decoder(enc$z, labels)
      kl_loss <- -0.5 * torch_mean(
        enc$z_log_var - enc$z_mean^2 - torch_exp(enc$z_log_var) + 1
      )
      list(reconstructed = reconstructed, kl = kl_loss)
    }
  )

  vae <- VAE()
  vae <- vae$to(device = device_obj)

  optimizer <- optim_adam(vae$parameters, lr = lr)
  mse_loss_fn <- function(pred, target) torch_sum(torch_mean((pred - target)^2, dim = 2))

  expand_flat_indices <- function(flat_idx, n_ct) {
    sample_idx <- ((flat_idx - 1L) %/% n_ct) + 1L
    ct_idx <- ((flat_idx - 1L) %% n_ct) + 1L
    list(sample_idx = sample_idx, ct_idx = ct_idx)
  }

  build_sim_batch <- function(flat_idx) {
    pos <- expand_flat_indices(flat_idx, n_ct)
    sim_idx <- pos$sample_idx
    ct_idx <- pos$ct_idx

    x_sim_batch <- X_sim_scaled_t[sim_idx, , drop = FALSE]

    labels_mat <- hot_encoding[ct_idx, , drop = FALSE]
    gt_mat <- matrix(0, nrow = length(flat_idx), ncol = n_genes)
    frac_vec <- numeric(length(flat_idx))

    for (j in seq_len(n_ct)) {
      take <- which(ct_idx == j)
      if (!length(take)) next

      ct_j <- ct_names[j]
      layer_block <- layers[[ct_j]][sim_idx[take], layer_gene_idx[[ct_j]], drop = FALSE]

      if (!is.null(normalize_simulated) && normalize_simulated == "cpm") {
        layer_block <- sweep(layer_block, 1, layer_libs[[ct_j]][sim_idx[take]], "/") * 1e6
      }

      layer_block <- log1p(layer_block)
      gt_mat[take, ] <- layer_block / max_val
      frac_vec[take] <- props_sim[sim_idx[take], j] / max_val
    }

    list(
      x_genes = x_sim_batch,
      labels = torch_tensor(labels_mat, dtype = torch_float32(), device = device_obj),
      gt = torch_tensor(gt_mat, dtype = torch_float32(), device = device_obj),
      frac = torch_tensor(frac_vec, dtype = torch_float32(), device = device_obj)
    )
  }

  build_real_batch <- function(flat_idx) {
    pos <- expand_flat_indices(flat_idx, n_ct)
    real_idx <- pos$sample_idx
    ct_idx <- pos$ct_idx

    labels_mat <- hot_encoding[ct_idx, , drop = FALSE]
    frac_vec <- fractions[cbind(real_idx, ct_idx)]

    list(
      x_genes = X_real_scaled_t[real_idx, , drop = FALSE],
      labels = torch_tensor(labels_mat, dtype = torch_float32(), device = device_obj),
      frac = torch_tensor(frac_vec, dtype = torch_float32(), device = device_obj)
    )
  }

  total_steps <- if (is.null(n_steps_expr)) 5000L * n_ct else as.integer(n_steps_expr)
  n_epochs <- as.integer((total_steps * 128L) / (n_genes + n_ct))
  steps_per_epoch <- as.integer((n_sim * n_ct) / batch_size)
  total_updates <- max(1L, n_epochs * steps_per_epoch)

  sim_pos <- 1L
  real_pos <- 1L

  next_seq_batch <- function(n_total, batch_size_now, pos_ref) {
    idx <- ((pos_ref - 1L) + seq_len(batch_size_now) - 1L) %% n_total + 1L
    next_pos <- ((pos_ref - 1L) + batch_size_now) %% n_total + 1L
    list(idx = idx, next_pos = next_pos)
  }

  prep_elapsed <- as.numeric(difftime(Sys.time(), expr_total_start, units = "secs"))
  cat(sprintf(
    "[%s] Running expression estimation (%d epoch%s, %d step%s/epoch) on %s\n",
    format_hms(prep_elapsed),
    n_epochs, if (n_epochs == 1L) "" else "s",
    steps_per_epoch, if (steps_per_epoch == 1L) "" else "s",
    backend$device_name
  ))

  train_start <- Sys.time()
  update_idx <- 0L
  last_loss <- NA_real_

  vae$train()
  if (n_epochs > 0 && steps_per_epoch > 0) {
    for (epoch in seq_len(n_epochs)) {
      for (step in seq_len(steps_per_epoch)) {
        update_idx <- update_idx + 1L

        sim_batch_idx <- next_seq_batch(n_sim * n_ct, batch_size, sim_pos)
        sim_pos <- sim_batch_idx$next_pos

        real_batch_idx <- next_seq_batch(n_real * n_ct, batch_size, real_pos)
        real_pos <- real_batch_idx$next_pos

        sim_batch <- build_sim_batch(sim_batch_idx$idx)
        real_batch <- build_real_batch(real_batch_idx$idx)

        beta <- torch_tensor(runif(1, 0.1, 0.9), dtype = torch_float32(), device = device_obj)

        batch_size_actual <- min(sim_batch$x_genes$size(1), real_batch$x_genes$size(1))
        x_train_sim_genes <- sim_batch$x_genes[1:batch_size_actual, , drop = FALSE]
        x_train_sim_gt <- sim_batch$gt[1:batch_size_actual, , drop = FALSE]
        labels_sim <- sim_batch$labels[1:batch_size_actual, , drop = FALSE]
        true_fractions_sim <- sim_batch$frac[1:batch_size_actual]

        x_train_real_genes <- real_batch$x_genes[1:batch_size_actual, , drop = FALSE]
        labels_real <- real_batch$labels[1:batch_size_actual, , drop = FALSE]
        true_fractions_real <- real_batch$frac[1:batch_size_actual]

        x_train_mix <- beta * x_train_real_genes + (1 - beta) * x_train_sim_genes
        true_fractions_mix <- beta * true_fractions_real + (1 - beta) * true_fractions_sim

        optimizer$zero_grad()

        reconstructed_sim <- vae(x_train_sim_genes, labels_sim)
        recon_loss <- mse_loss_fn(reconstructed_sim$reconstructed, x_train_sim_gt)

        reconstructed_mix <- vae(x_train_mix, labels_sim)
        reconstructed_real <- vae(x_train_real_genes, labels_real)

        recon_mix_weighted <- reconstructed_mix$reconstructed * true_fractions_mix$unsqueeze(2)
        recon_real_weighted <- reconstructed_real$reconstructed * true_fractions_real$unsqueeze(2)
        recon_sim_weighted <- reconstructed_sim$reconstructed * true_fractions_sim$unsqueeze(2)
        expected_mix <- beta * recon_real_weighted + (1 - beta) * recon_sim_weighted

        consistency_loss <- mse_loss_fn(recon_mix_weighted, expected_mix)
        total_loss <- recon_loss + beta_vae * reconstructed_sim$kl + lambda_cons * consistency_loss

        total_loss$backward()
        optimizer$step()

        if (backend$using_cuda) {
          torch::cuda_synchronize()
        }

        last_loss <- as.numeric(total_loss$item())
        elapsed_expr <- as.numeric(difftime(Sys.time(), train_start, units = "secs"))

        cat(sprintf(
          "\r[%s] Epoch %d/%d step %d/%d | loss=%.6f | elapsed=%s",
          make_bar(update_idx, total_updates),
          epoch, n_epochs,
          step, steps_per_epoch,
          last_loss,
          format_hms(elapsed_expr)
        ))
        utils::flush.console()
      }
    }
  }

  if (backend$using_cuda) {
    torch::cuda_synchronize()
  }

  train_elapsed <- as.numeric(difftime(Sys.time(), train_start, units = "secs"))
  cat(sprintf(
    "\n[%s] Expression estimation finished | final_loss=%.6f\n",
    format_hms(train_elapsed),
    last_loss
  ))

  pred_start <- Sys.time()
  cat(sprintf("[%s] Reconstructing cell-type-specific expression\n",
              format_hms(as.numeric(difftime(Sys.time(), expr_total_start, units = "secs")))))

  vae$eval()

  sample_names <- rownames(fractions)
  if (is.null(sample_names)) {
    sample_names <- paste0("sample_", seq_len(n_real))
  }

  expression_layered <- setNames(vector("list", length(ct_names)), ct_names)
  for (ct in ct_names) {
    expression_layered[[ct]] <- matrix(0, nrow = n_real, ncol = n_genes)
    rownames(expression_layered[[ct]]) <- sample_names
    colnames(expression_layered[[ct]]) <- common_genes
  }

  pred_batch_size <- max(batch_size, 256L)

  with_no_grad({
    for (j in seq_len(n_ct)) {
      for (start in seq(1L, n_real, by = pred_batch_size)) {
        end <- min(start + pred_batch_size - 1L, n_real)
        idx <- start:end

        labels_mat <- hot_encoding[rep.int(j, length(idx)), , drop = FALSE]
        labels_t <- torch_tensor(labels_mat, dtype = torch_float32(), device = device_obj)
        real_genes_t <- X_real_scaled_t[idx, , drop = FALSE]

        est_block <- vae(real_genes_t, labels_t)$reconstructed
        est_block <- est_block * max_val

        if (backend$using_cuda) {
          torch::cuda_synchronize()
        }

        est_block <- as.matrix(est_block$to(device = torch_device("cpu")))
        expression_layered[[ct_names[j]]][idx, ] <- est_block
      }
    }
  })

  pred_elapsed <- as.numeric(difftime(Sys.time(), pred_start, units = "secs"))
  cat(sprintf("[%s] Reconstruction finished | step_time=%s\n",
              format_hms(as.numeric(difftime(Sys.time(), expr_total_start, units = "secs"))),
              format_hms(pred_elapsed)))

  est <- matrix(0, nrow = n_real * n_ct, ncol = n_genes)
  for (j in seq_len(n_ct)) {
    rows_j <- seq.int(from = j, to = n_real * n_ct, by = n_ct)
    est[rows_j, ] <- expression_layered[[ct_names[j]]]
  }

  est[est < 0] <- 0
  colnames(est) <- common_genes

  ct_labels <- rep(ct_names, times = n_real)
  sample_labels <- rep(sample_names, each = n_ct)

  # scanpy's sc.pp.scale(max_value = 10) clips only the upper tail and uses a
  # population sd (ddof = 0); base R scale() uses ddof = 1 and no clipping.
  ctr <- colMeans(est)
  sdv <- sqrt(colMeans(sweep(est, 2, ctr, "-")^2))
  sdv[sdv == 0] <- 1
  scaled_counts <- sweep(sweep(est, 2, ctr, "-"), 2, sdv, "/")
  scaled_counts[scaled_counts > 10] <- 10
  colnames(scaled_counts) <- common_genes

  expression_combined <- as.data.frame(est)
  expression_combined$cell_type <- ct_labels
  expression_combined$sample <- sample_labels

  total_elapsed <- as.numeric(difftime(Sys.time(), expr_total_start, units = "secs"))
  cat(sprintf("[%s] Expression workflow done | total_time=%s\n",
              format_hms(total_elapsed), format_hms(total_elapsed)))

  list(
    expression_layered = expression_layered,
    expression_combined = expression_combined,
    scaled_counts = scaled_counts,
    celltypes = ct_names,
    genes = common_genes
  )
}











#' Run the full DISSECT workflow
#'
#' Convenience wrapper that runs simulation, data processing, fraction
#' estimation, and optionally expression estimation in sequence.
#'
#' @param sc_data A `Seurat` or `SingleCellExperiment` object, or `NULL`. If
#'   supplied, simulated mixtures are generated with [dissect_simulate()].
#' @param bulk Numeric matrix with genes in rows and samples in columns.
#' @param reference Numeric matrix with genes in rows and cell types in
#'   columns, or `NULL`.
#' @param celltype_col Character scalar. Column name in the single-cell metadata
#'   containing cell-type labels.
#' @param batch_col Character scalar or `NULL`. Optional column name in the
#'   single-cell metadata containing batch labels.
#' @param type Character scalar. Simulation type, `"bulk"` or `"st"`.
#' @param n_samples Integer scalar or `NULL`. Number of simulated samples.
#' @param cells_per_sample Integer scalar. Number of cells per bulk sample.
#' @param prop_sparse Numeric scalar. Fraction of sparse bulk samples.
#' @param concentration Numeric vector or `NULL`. Dirichlet concentration
#'   parameter for bulk simulation.
#' @param save_expr Logical scalar. If `TRUE`, stores simulated per-cell-type
#'   expression and enables downstream expression estimation.
#' @param min_genes Integer scalar. Minimum genes per cell for preprocessing.
#' @param min_cells Integer scalar. Minimum cells per gene for preprocessing.
#' @param mt_cutoff Numeric scalar. Mitochondrial percentage cutoff.
#' @param min_expr Numeric scalar. Minimum mean `log1p` expression threshold.
#' @param downsample Numeric scalar or `NULL`. Downsampling factor for ST
#'   simulation.
#' @param test_dataset_type Character scalar. Either `"bulk"` or
#'   `"microarray"`.
#' @param duplicated Character scalar. How duplicated genes should be handled.
#' @param normalize_simulated Character scalar or `NULL`. Currently `"cpm"` or
#'   `NULL`.
#' @param normalize_test Character scalar or `NULL`. Currently `"cpm"` or
#'   `NULL`.
#' @param var_cutoff Numeric scalar or `NULL`. Variance threshold for bulk
#'   preprocessing.
#' @param test_in_mix Integer scalar. Number of real samples used in online
#'   mixing.
#' @param n_hidden_layers Integer scalar. Number of hidden layers in the
#'   fraction model.
#' @param hidden_units Integer vector. Units per hidden layer.
#' @param hidden_activation Character scalar. Hidden-layer activation.
#' @param output_activation Character scalar. Output-layer activation.
#' @param loss Character scalar. Fraction-model loss.
#' @param n_steps Integer scalar. Fraction-model training steps.
#' @param lr Numeric scalar. Fraction-model learning rate.
#' @param batch_size Integer scalar. Fraction-model batch size.
#' @param dropout Numeric vector or `NULL`. Fraction-model dropout rates.
#' @param alpha_range Numeric vector of length two. Mixing coefficient range.
#' @param normalization_per_batch Character scalar or `NULL`. Batch
#'   normalisation mode for the fraction model.
#' @param models Integer vector. Ensemble model identifiers.
#' @param mix Character scalar. Mixing strategy.
#' @param n_steps_expr Integer scalar or `NULL`. Expression-model training
#'   steps.
#' @param expr_scaling Character scalar. Expression scaling mode.
#' @param latent_dim Integer scalar. Expression-model latent dimension.
#' @param batch_size_expr Integer scalar. Expression-model batch size.
#' @param lr_expr Numeric scalar. Expression-model learning rate.
#' @param beta_vae Numeric scalar. Weight of the KL term in the expression
#'   model.
#' @param lambda_cons Numeric scalar. Weight of the consistency term in the
#'   expression model.
#' @param seed Integer scalar. Random seed.
#' @param device Character scalar. One of `"auto"`, `"cpu"`, or `"cuda"`.
#' @param cuda_index Integer scalar or `NULL`. Optional CUDA device index.
#'
#' @return A named list with components:
#' \describe{
#'   \item{fractions}{Estimated cell-type fractions.}
#'   \item{scores}{Estimated pre-activation fraction scores.}
#'   \item{expression}{Expression estimation results, or `NULL` if not run.}
#'   \item{sim_data}{Simulation results, or `NULL` if no single-cell object was
#'   provided.}
#'   \item{processed}{Processed matrices used for training and prediction.}
#' }
#'
#' @details
#' This is the main end-to-end wrapper for the R implementation of DISSECT. If
#' `sc_data` is provided, the workflow uses simulated training data generated
#' from the single-cell reference. If `save_expr = TRUE`, the expression model
#' is trained after fraction estimation.
#'
#' Elapsed time for each major step and the total runtime are reported.
#' 
#' 
#' @references
#' Khatri, R., Machart, P., & Bonn, S. (2024). DISSECT: deep semi-supervised
#' consistency regularization for accurate cell type fraction and gene
#' expression estimation. \emph{Genome Biology}, 25(1), 112.
#'
#' Original DISSECT software repository:
#' \url{https://github.com/imsb-uke/DISSECT}
#'
#' @examples
#' \dontrun{
#' res <- dissect(
#'   sc_data = sce,
#'   bulk = bulk_mat,
#'   celltype_col = "celltype",
#'   batch_col = "batch",
#'   device = "auto"
#' )
#' }
#'
#' @export
dissect <- function(sc_data = NULL,
                    bulk,
                    reference = NULL,
                    celltype_col = "celltype",
                    batch_col = NULL,
                    type = "bulk",
                    n_samples = NULL,
                    cells_per_sample = 500,
                    prop_sparse = 0.5,
                    concentration = NULL,
                    save_expr = TRUE,
                    min_genes = 200,
                    min_cells = 3,
                    mt_cutoff = 5,
                    min_expr = 0,
                    downsample = NULL,
                    test_dataset_type = "bulk",
                    duplicated = "first",
                    normalize_simulated = "cpm",
                    normalize_test = "cpm",
                    var_cutoff = 0.1,
                    test_in_mix = 1,
                    n_hidden_layers = 4,
                    hidden_units = c(512, 256, 128, 64),
                    hidden_activation = "relu6",
                    output_activation = "softmax",
                    loss = "kldivergence",
                    n_steps = 5000,
                    lr = 1e-5,
                    batch_size = 64,
                    dropout = NULL,
                    alpha_range = c(0.1, 0.9),
                    normalization_per_batch = "log1p-MinMax",
                    models = c(1, 2, 3, 4, 5),
                    mix = "srm",
                    n_steps_expr = 5000,
                    expr_scaling = "p99",
                    latent_dim = 128,
                    batch_size_expr = 128,
                    lr_expr = 1e-3,
                    beta_vae = 0.01,
                    lambda_cons = 0.1,
                    seed = 42,
                    device = c("auto", "cpu", "cuda"),
                    cuda_index = NULL) {

  format_hms <- function(seconds) {
    seconds <- max(0, as.integer(round(seconds)))
    hh <- seconds %/% 3600
    mm <- (seconds %% 3600) %/% 60
    ss <- seconds %% 60
    sprintf("%02d:%02d:%02d", hh, mm, ss)
  }

  dissect_ascii <- r"(
    Running
       ____  _                      __ 
      / __ \(_)_____________  _____/ /_
     / / / / / ___/ ___/ _ \/ ___/ __/
    / /_/ / (__  |__  )  __/ /__/ /_  
   /_____/_/____/____/\___/\___/\__/  
  )"
  cat("\033[35m", "\033[1m", dissect_ascii, "\033[0m", sep = "")
  cat("\n")
  total_start_time <- Sys.time()
  cat("[00:00:00] DISSECT workflow started\n")

  backend <- dissect_resolve_device(
    device = device,
    cuda_index = cuda_index,
    verbose = TRUE,
    stage = "workflow"
  )

  sim_data <- NULL

  if (!is.null(sc_data)) {
    step_start <- Sys.time()
    cat("[00:00:00] Running simulation\n")
    sim_data <- dissect_simulate(
      sc_data = sc_data,
      celltype_col = celltype_col,
      batch_col = batch_col,
      type = type,
      n_samples = n_samples,
      cells_per_sample = cells_per_sample,
      prop_sparse = prop_sparse,
      concentration = concentration,
      save_expr = save_expr,
      min_genes = min_genes,
      min_cells = min_cells,
      mt_cutoff = mt_cutoff,
      min_expr = min_expr,
      downsample = downsample,
      seed = seed
    )
    step_elapsed <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
    total_elapsed <- as.numeric(difftime(Sys.time(), total_start_time, units = "secs"))
    cat(sprintf("[%s] Simulation finished | step_time=%s\n",
                format_hms(total_elapsed), format_hms(step_elapsed)))
  }

  step_start <- Sys.time()
  cat(sprintf("[%s] Running preprocessing\n",
              format_hms(as.numeric(difftime(Sys.time(), total_start_time, units = "secs")))))
  processed <- dissect_process(
    bulk = bulk,
    reference = reference,
    sim_data = sim_data,
    test_dataset_type = test_dataset_type,
    duplicated = duplicated,
    normalize_simulated = normalize_simulated,
    normalize_test = normalize_test,
    var_cutoff = var_cutoff,
    test_in_mix = test_in_mix
  )
  step_elapsed <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
  total_elapsed <- as.numeric(difftime(Sys.time(), total_start_time, units = "secs"))
  cat(sprintf("[%s] Preprocessing finished | step_time=%s\n",
              format_hms(total_elapsed), format_hms(step_elapsed)))

  step_start <- Sys.time()
  cat(sprintf("[%s] Running proportion estimation\n",
              format_hms(as.numeric(difftime(Sys.time(), total_start_time, units = "secs")))))
  prop_result <- dissect_prop(
    processed = processed,
    test_dataset_type = test_dataset_type,
    duplicated = duplicated,
    normalize_simulated = normalize_simulated,
    normalize_test = normalize_test,
    var_cutoff = var_cutoff,
    test_in_mix = test_in_mix,
    n_hidden_layers = n_hidden_layers,
    hidden_units = hidden_units,
    hidden_activation = hidden_activation,
    output_activation = output_activation,
    loss = loss,
    n_steps = n_steps,
    lr = lr,
    batch_size = batch_size,
    dropout = dropout,
    alpha_range = alpha_range,
    normalization_per_batch = normalization_per_batch,
    models = models,
    mix = mix,
    device = device,
    cuda_index = cuda_index
  )
  step_elapsed <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
  total_elapsed <- as.numeric(difftime(Sys.time(), total_start_time, units = "secs"))
  cat(sprintf("[%s] Proportion estimation wrapper finished | step_time=%s\n",
              format_hms(total_elapsed), format_hms(step_elapsed)))

  expr_result <- NULL
  if (!is.null(sim_data) && save_expr) {
    step_start <- Sys.time()
    cat(sprintf("[%s] Running expression estimation\n",
                format_hms(as.numeric(difftime(Sys.time(), total_start_time, units = "secs")))))
    expr_result <- tryCatch({
      dissect_expr(
        bulk = bulk,
        fractions = prop_result$fractions,
        sim_data = sim_data,
        normalize_simulated = normalize_simulated,
        normalize_test = normalize_test,
        n_steps_expr = n_steps_expr,
        expr_scaling = expr_scaling,
        latent_dim = latent_dim,
        batch_size = batch_size_expr,
        lr = lr_expr,
        beta_vae = beta_vae,
        lambda_cons = lambda_cons,
        seed = seed,
        device = device,
        cuda_index = cuda_index
      )
    }, error = function(e) {
      warning("Expression estimation failed: ", conditionMessage(e),
              "\nReturning proportion results only.")
      NULL
    })
    step_elapsed <- as.numeric(difftime(Sys.time(), step_start, units = "secs"))
    total_elapsed <- as.numeric(difftime(Sys.time(), total_start_time, units = "secs"))
    cat(sprintf("[%s] Expression estimation %s | step_time=%s\n",
                format_hms(total_elapsed),
                if (is.null(expr_result)) "FAILED" else "finished",
                format_hms(step_elapsed)))
  }

  total_elapsed <- as.numeric(difftime(Sys.time(), total_start_time, units = "secs"))
  cat(sprintf("[%s] DISSECT workflow finished | total_time=%s\n",
              format_hms(total_elapsed), format_hms(total_elapsed)))

  list(
    fractions = prop_result$fractions,
    scores = prop_result$scores,
    expression = expr_result,
    sim_data = sim_data,
    processed = processed
  )
}