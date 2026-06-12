#' Simulate pseudo-bulk RNA-seq profiles from single-cell data
#'
#' Generates pseudo-bulk expression profiles by sampling cells from a
#' single-cell reference according to Dirichlet-drawn cell-type fractions.
#' Optionally returns per-cell-type signature profiles and a global signature matrix.
#'
#' @param ... The single-cell reference, supplied in one of two forms:
#'   \itemize{
#'     \item A single object: a \code{Seurat} object (counts taken from the
#'       \code{"RNA"} assay) or a \code{SingleCellExperiment} (counts taken from
#'       \code{counts()}).
#'     \item Two objects: a genes \eqn{\times} cells count matrix followed by a
#'       cell-metadata \code{data.frame} whose row names match the column names
#'       of the count matrix.
#'   }
#' @param n_bulk_samples Integer number of pseudo-bulk samples to generate. If
#'   \code{NULL} (default), \code{1000 * (number of cell types)} samples are
#'   generated.
#' @param cells_per_bulk Integer number of cells pooled into each pseudo-bulk
#'   sample. Default \code{500}.
#' @param cell_type_column Name of the metadata column holding cell-type labels.
#'   Default \code{"cell_type"}.
#' @param patient_id_column Optional name of a metadata column holding
#'   patient/donor IDs. When supplied, each pseudo-bulk is drawn from a single
#'   donor where that donor has cells of the required type (falling back to the
#'   full cell-type pool otherwise). Default \code{NULL} (no patient structure).
#' @param select_ct Optional character vector restricting which cell types are
#'   used (and fixing their order in the output). Default \code{NULL} (all cell
#'   types).
#' @param dirich_alpha Concentration parameter of the symmetric Dirichlet used
#'   to draw cell-type fractions. Smaller values give more skewed mixtures.
#'   Default \code{1} (uniform over the simplex).
#' @param seet Integer random seed (passed to \code{set.seed}) for
#'   reproducibility. Default \code{1}.
#' @param sparse Logical; if \code{TRUE}, randomly zero out a fraction of
#'   cell-type entries before renormalisation, producing samples in which some
#'   cell types are absent. Default \code{FALSE}.
#' @param sparse_prob Probability that a given cell-type fraction is dropped
#'   when \code{sparse = TRUE}. Default \code{0.5}.
#' @param rare Logical; if \code{TRUE}, force a random subset of cell-type
#'   fractions to small values in \code{[0, 0.03]} before renormalisation, to
#'   emulate rare populations. Default \code{FALSE}.
#' @param rare_percentage Probability that a given cell-type fraction is made
#'   rare when \code{rare = TRUE}. Default \code{0.4}.
#' @param return_used_samples Logical; if \code{TRUE}, include the per-sample,
#'   per-cell-type cell indices that were drawn. Default \code{FALSE}.
#' @param return_patient_metadata Logical; if \code{TRUE} and
#'   \code{patient_id_column} is set, include a data frame mapping each
#'   pseudo-bulk to its source donor. Default \code{FALSE}.
#' @param return_signature_matrix Logical; if \code{TRUE}, also compute
#'   per-cell-type signature profiles (one genes \eqn{\times} bulk matrix per
#'   cell type) and a global genes \eqn{\times} cell-types signature matrix.
#'   Default \code{FALSE}.
#' @param verbose Logical; if \code{TRUE} (default), print the header, progress
#'   bars, and timing summary. If \code{FALSE}, nothing is printed.
#'
#' @details
#' Cell-type fractions are drawn from a symmetric Dirichlet, optionally modified
#' by the \code{sparse}/\code{rare} masks, renormalised per sample, and turned
#' into integer cell counts (each bulk is guaranteed at least one cell).
#' Cells are then sampled with replacement from the reference and summed to form each pseudo-bulk,
#' so the returned \code{ground_truth_proportions} reflect the *realized* allocations rather
#' than the raw Dirichlet draws.
#'
#'
#' @return A named list containing:
#'   \describe{
#'     \item{\code{bulk_expression_profiles}}{Genes \eqn{\times}
#'       \code{n_bulk_samples} matrix of summed pseudo-bulk counts.}
#'     \item{\code{ground_truth_proportions}}{\code{n_bulk_samples} \eqn{\times}
#'       cell-types matrix of realized proportions (rows sum to 1).}
#'     \item{\code{used_samples_by_ct}}{(if \code{return_used_samples}) list of
#'       length \code{n_bulk_samples}, each a per-cell-type list of drawn cell
#'       indices.}
#'     \item{\code{bulk_patient_metadata}}{(if \code{return_patient_metadata}
#'       and patient mode) data frame mapping \code{sample_id} to
#'       \code{patient_id}.}
#'     \item{\code{bulk_signature_profiles}}{(if \code{return_signature_matrix})
#'       list of per-cell-type genes \eqn{\times} bulk mean-expression
#'       matrices.}
#'     \item{\code{global_signature_matrix}}{(if \code{return_signature_matrix})
#'       genes \eqn{\times} cell-types matrix of mean signatures.}
#'     \item{\code{timing}}{list with \code{pseudobulk_seconds},
#'       \code{signature_seconds}, and \code{total_seconds}.}
#'   }
#'
#' @examples
#' 
#' ## Tiny synthetic reference: 1000 genes, 3 cell types, 100 cells each
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
#' rownames(counts) <- paste0("gene_", seq_len(n_genes))
#' colnames(counts) <- paste0("cell_", seq_len(n_cells))
#'
#' meta <- data.frame(
#'   cell_type = rep(cell_types, each = cells_per_type),
#'   row.names = colnames(counts)
#' )
#'
#' ## Generate 50 pseudo-bulks of 100 cells each, with signatures
#' res <- quasar_sim_bulk(
#'   counts, meta,
#'   n_bulk_samples          = 50,
#'   cells_per_bulk          = 100,
#'   return_signature_matrix = TRUE,
#'   verbose                 = TRUE
#' )
#'
#' dim(res$bulk_expression_profiles)   # 1000 x 50
#' head(res$ground_truth_proportions)  # rows sum to 1
#' dim(res$global_signature_matrix)    # 1000 x 3
#' 
#'
#' @importFrom MCMCpack rdirichlet
#' @importFrom Matrix sparseMatrix Diagonal Matrix rowMeans
#' @importFrom stats runif
#' @importFrom utils head flush.console
#' @export

quasar_sim_bulk <- function(...,
                            n_bulk_samples = NULL,
                            cells_per_bulk = 500,
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
                            verbose = TRUE) {

  start_time_all <- Sys.time()
  args <- list(...)

  # --------------------------------------------------
  # Yellow multi-line progress bar 
  # --------------------------------------------------
  .quasar_progress <- function(done, total, start_time, header, count_label,
                               first = FALSE, block_lines = 5L) {
    percent   <- if (total > 0) (done / total) * 100 else 100
    bar_width <- 30
    filled    <- round(percent / (100 / bar_width))
    filled    <- max(0L, min(bar_width, as.integer(filled)))
    bar       <- paste0("[", strrep("\u2588", filled),
                        strrep(" ", bar_width - filled), "]")

    elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
    rate    <- if (done > 0) elapsed / done else 0
    eta     <- rate * (total - done)

    fmt_t <- function(s) {
      if (!is.finite(s) || s < 0) return("--:--")
      s <- as.integer(round(s))
      sprintf("%02d:%02d", s %/% 60L, s %% 60L)
    }

    lines <- c(
      sprintf("\033[1m%s\033[0m\033[33m", header),
      sprintf("%s %6.2f%%", bar, percent),
      sprintf("%-14s: %d / %d", count_label, done, total),
      sprintf("%-14s: %s", "Elapsed", fmt_t(elapsed)),
      sprintf("%-14s: %s", "ETA", fmt_t(eta))
    )

    colored  <- paste0("\033[33m", lines, "\033[0m")
    move_up  <- if (first) "" else sprintf("\033[%dA\r", block_lines)
    body     <- paste0("\033[2K", colored, "\n", collapse = "")
    cat(move_up, body, sep = "")
    flush.console()
  }

  if (verbose) {
    cat("\033[33m", "-------------------------------", "\033[0m", "\n", sep = "")
    cat("\033[33m", "\033[1m", "Generating pseudo-bulks", "\033[0m", "\n", sep = "")
    cat("\033[33m", "-------------------------------", "\033[0m", "\n", sep = "")
  }

  # --------------------------------------------------
  # Load count matrix & metadata
  # --------------------------------------------------
  if (length(args) == 1) {
    obj <- args[[1]]
    if (inherits(obj, "Seurat")) {
      count_matrix <- obj[["RNA"]]$counts
      cell_metadata <- obj@meta.data
    } else if (inherits(obj, "SingleCellExperiment")) {
      count_matrix <- SingleCellExperiment::counts(obj)
      cell_metadata <- as.data.frame(SummarizedExperiment::colData(obj))
    } else {
      stop("Single input must be Seurat or SingleCellExperiment object")
    }
  } else if (length(args) == 2) {
    count_matrix <- args[[1]]
    cell_metadata <- as.data.frame(args[[2]])
    common_cells <- intersect(colnames(count_matrix), rownames(cell_metadata))
    if (!length(common_cells))
      stop("No matching cells between counts and metadata")
    count_matrix <- count_matrix[, common_cells, drop = FALSE]
    cell_metadata <- cell_metadata[common_cells, , drop = FALSE]
  } else {
    stop("Invalid input format")
  }

  if (!cell_type_column %in% colnames(cell_metadata))
    stop("Missing cell type column in metadata")

  cell_types <- as.character(cell_metadata[[cell_type_column]])
  names(cell_types) <- colnames(count_matrix)

  # --------------------------------------------------
  # Restrict cell types
  # --------------------------------------------------
  if (!is.null(select_ct)) {
    keep_idx <- which(cell_types %in% select_ct)
    if (!length(keep_idx))
      stop("No cells found for given select_ct")
    count_matrix <- count_matrix[, keep_idx, drop = FALSE]
    cell_metadata <- cell_metadata[keep_idx, , drop = FALSE]
    cell_types <- cell_types[keep_idx]
    celltype_levels <- select_ct
  } else {
    celltype_levels <- unique(cell_types)
  }

  # --------------------------------------------------
  # Patient mode
  # --------------------------------------------------
  patient_mode <- !is.null(patient_id_column)
  if (patient_mode) {
    if (!patient_id_column %in% colnames(cell_metadata))
      stop("Missing patient ID column")
    patient_ids <- as.character(cell_metadata[[patient_id_column]])
    patient_levels <- unique(patient_ids)
  }

  set.seed(seet)

  # --------------------------------------------------
  # Precompute indices per cell type
  # --------------------------------------------------
  celltype_indices <- lapply(celltype_levels, function(ct) which(cell_types == ct))
  names(celltype_indices) <- celltype_levels

  if (is.null(n_bulk_samples)) {
    if (verbose) {
      cat("\033[33m",
          sprintf("n_bulk_samples = NULL - ct * 1000 samples will be generated: %d samples.\n",
                  length(celltype_levels) * 1000),
          "\033[0m", sep = "")
    }
    n_bulk_samples <- length(celltype_levels) * 1000
  }

  # Patient pools (optional)
  if (patient_mode) {
    pools_by_pat <- lapply(celltype_levels, function(ct) {
      idx <- celltype_indices[[ct]]
      split(idx, patient_ids[idx])
    })
    names(pools_by_pat) <- celltype_levels
    chosen_patients <- sample(patient_levels, n_bulk_samples, replace = TRUE)
  }

  # --------------------------------------------------
  # Dirichlet fractions (CT x bulk)
  # --------------------------------------------------
  frac_matrix <- t(
    MCMCpack::rdirichlet(
      n_bulk_samples,
      rep(dirich_alpha, length(celltype_levels))
    )
  )

  # --------------------------------------------------
  # Sparsity / rare masking
  # --------------------------------------------------
  if (sparse) {
    drop_mask <- matrix(
      runif(length(frac_matrix)) < sparse_prob,
      nrow = nrow(frac_matrix)
    )
    frac_matrix[drop_mask] <- 0
  }

  if (rare) {
    rare_mask <- matrix(
      runif(length(frac_matrix)) < rare_percentage,
      nrow = nrow(frac_matrix)
    )
    frac_matrix[rare_mask] <- runif(sum(rare_mask), 0, 0.03)
  }


  col_sums <- colSums(frac_matrix)
  valid <- col_sums > 0
  frac_matrix[, valid] <-
    sweep(frac_matrix[, valid, drop = FALSE], 2, col_sums[valid], "/")

  # --------------------------------------------------
  # Allocate cells per CT per bulk
  # --------------------------------------------------
  counts_per_ct_all <- floor(frac_matrix * cells_per_bulk)

  # Enforce >=1 cell per bulk
  empty_bulks <- which(colSums(counts_per_ct_all) == 0)
  if (length(empty_bulks)) {
    for (i in empty_bulks) {
      k <- sample(seq_len(nrow(counts_per_ct_all)), 1)
      counts_per_ct_all[k, i] <- 1
    }
  }


  start_time_pseudobulk <- Sys.time()

  # --------------------------------------------------
  # FAST sampling with preallocation
  # --------------------------------------------------
  K <- length(celltype_levels)
  C <- ncol(count_matrix)
  N <- n_bulk_samples
  G <- nrow(count_matrix)

  used_samples_by_ct <- vector("list", N)

  total_picks <- sum(counts_per_ct_all)
  ids_all <- integer(total_picks)
  grp_all <- integer(total_picks)
  pos_all <- 1L

  total_picks_ct <- rowSums(counts_per_ct_all)
  ids_ct <- lapply(total_picks_ct, function(n) integer(n))
  grp_ct <- lapply(total_picks_ct, function(n) integer(n))
  pos_ct <- rep.int(1L, K)

  # Contiguous ranges in ids_all/grp_all per bulk sample
  # (ids are filled in increasing i-order, so each sample owns a contiguous block)
  samp_sizes <- colSums(counts_per_ct_all)
  ends   <- cumsum(samp_sizes)
  starts <- c(1, utils::head(ends, -1L) + 1)

  # Dense output, filled chunk-by-chunk *inside* the loop so the heavy
  # matmul/densification is tracked by the progress bar (not done in one
  # opaque blocking call at the end).
  bulk_expression <- matrix(0, nrow = G, ncol = N)

  # ~100 refreshes; each chunk does sampling + its sparse matmul block
  pb_chunk <- max(1L, floor(N / 100))
  chunk_start <- 1L

  if (verbose) {
    .quasar_progress(0L, N, start_time_pseudobulk,
                     header = "Creating pseudo-bulks",
                     count_label = "Samples", first = TRUE)
  }

  for (i in seq_len(N)) {

    picks_ct_list <- vector("list", K)
    names(picks_ct_list) <- celltype_levels

    if (patient_mode) {
      pat_i <- chosen_patients[i]
    }

    for (k in seq_along(celltype_levels)) {
      ct <- celltype_levels[k]
      nct <- counts_per_ct_all[k, i]

      if (nct <= 0) {
        picks_ct_list[[k]] <- integer(0)
        next
      }

      # Select pool (patient-specific if available, else fallback to full CT pool)
      pool <- NULL
      if (patient_mode) {
        tmp <- pools_by_pat[[ct]][[pat_i]]
        if (!is.null(tmp) && length(tmp) > 0) pool <- tmp
      }
      if (is.null(pool) || length(pool) == 0) {
        pool <- celltype_indices[[ct]]
      }

      if (is.null(pool) || length(pool) == 0) {
        picks_ct_list[[k]] <- integer(0)
        next
      }

      picks <- sample(pool, nct, replace = TRUE)
      picks_ct_list[[k]] <- picks

      # Fill global selection
      rng_all <- pos_all:(pos_all + nct - 1L)
      ids_all[rng_all] <- picks
      grp_all[rng_all] <- i
      pos_all <- pos_all + nct

      # Fill per-CT selection
      p0 <- pos_ct[k]
      rng_ct <- p0:(p0 + nct - 1L)
      ids_ct[[k]][rng_ct] <- picks
      grp_ct[[k]][rng_ct] <- i
      pos_ct[k] <- p0 + nct
    }

    used_samples_by_ct[[i]] <- picks_ct_list

   
    if (i %% pb_chunk == 0L || i == N) {
      cs  <- chunk_start
      ce  <- i
      rng <- starts[cs]:ends[ce]

      S_block <- Matrix::sparseMatrix(
        i    = ids_all[rng],
        j    = grp_all[rng] - cs + 1L,
        x    = 1,
        dims = c(C, ce - cs + 1L)
      )

      bulk_expression[, cs:ce] <- as.matrix(count_matrix %*% S_block)
      chunk_start <- i + 1L

      if (verbose) {
        .quasar_progress(i, N, start_time_pseudobulk,
                         header = "Creating pseudo-bulks",
                         count_label = "Samples")
      }
    }
  }

  rownames(bulk_expression) <- rownames(count_matrix)
  colnames(bulk_expression) <- paste0("sample_", seq_len(N))

  # --------------------------------------------------
  # Ground truth proportions from realized allocations
  # --------------------------------------------------
  gt_props <- sweep(counts_per_ct_all, 2, colSums(counts_per_ct_all), "/")
  gt_props <- t(gt_props)
  colnames(gt_props) <- celltype_levels
  rownames(gt_props) <- paste0("sample_", seq_len(N))

  # ==================================================
  # TIMING END: pseudobulk creation
  # ==================================================
  end_time_pseudobulk <- Sys.time()
  time_pseudobulk_sec <- as.numeric(difftime(end_time_pseudobulk, start_time_pseudobulk, units = "secs"))

  result <- list(
    bulk_expression_profiles = bulk_expression,
    ground_truth_proportions = gt_props
  )

  if (return_used_samples) {
    result$used_samples_by_ct <- used_samples_by_ct
  }

  if (patient_mode && return_patient_metadata) {
    result$bulk_patient_metadata <- data.frame(
      sample_id = paste0("sample_", seq_len(N)),
      patient_id = chosen_patients
    )
  }

  # ==================================================
  # Signature matrix computation (separately timed)
  # ==================================================
  time_signature_sec <- NA_real_

  if (return_signature_matrix) {

    start_time_sig <- Sys.time()

    inv_n_list <- lapply(seq_len(K), function(k) {
      nct <- counts_per_ct_all[k, ]
      inv <- ifelse(nct > 0, 1 / nct, 0)
      Matrix::Diagonal(n = N, x = inv)
    })

    sig_list <- vector("list", K)
    names(sig_list) <- celltype_levels


    sig_total_steps <- 2L * K
    sig_done <- 0L
    sig_update_every <- max(1L, floor(sig_total_steps / 100))

    if (verbose) {
      .quasar_progress(0L, sig_total_steps, start_time_sig,
                       header = "Building signature matrices",
                       count_label = "Steps", first = TRUE)
    }

    for (k in seq_len(K)) {
      if (total_picks_ct[k] == 0) {
        sig_list[[k]] <- Matrix::Matrix(0, nrow = nrow(count_matrix), ncol = N, sparse = TRUE)
      } else {
        Wk <- Matrix::sparseMatrix(
          i = ids_ct[[k]],
          j = grp_ct[[k]],
          x = 1,
          dims = c(C, N)
        )

        Xk_sum  <- count_matrix %*% Wk         # genes x bulk (sum with multiplicity)
        Xk_mean <- Xk_sum %*% inv_n_list[[k]]  # genes x bulk (mean)
        sig_list[[k]] <- Xk_mean
      }

      sig_done <- sig_done + 1L
      if (verbose && (sig_done %% sig_update_every == 0L || sig_done == sig_total_steps)) {
        .quasar_progress(sig_done, sig_total_steps, start_time_sig,
                         header = "Building signature matrices",
                         count_label = "Steps")
      }
    }

    # ---------------------------------------------
    # global signature (genes x celltypes)
    # ---------------------------------------------
    global_sig <- matrix(0, nrow = nrow(count_matrix), ncol = K,
                         dimnames = list(rownames(count_matrix), celltype_levels))

    for (k in seq_len(K)) {
      m <- sig_list[[k]]
      if (inherits(m, "dgCMatrix")) {
        global_sig[, k] <- Matrix::rowMeans(m)
      } else {
        global_sig[, k] <- rowMeans(as.matrix(m))
      }

      sig_done <- sig_done + 1L
      if (verbose && (sig_done %% sig_update_every == 0L || sig_done == sig_total_steps)) {
        .quasar_progress(sig_done, sig_total_steps, start_time_sig,
                         header = "Building signature matrices",
                         count_label = "Steps")
      }
    }

    result$bulk_signature_profiles <- sig_list
    result$global_signature_matrix <- global_sig

    end_time_sig <- Sys.time()
    time_signature_sec <- as.numeric(difftime(end_time_sig, start_time_sig, units = "secs"))
  }

  end_time_all <- Sys.time()
  time_total_sec <- as.numeric(difftime(end_time_all, start_time_all, units = "secs"))

  result$timing <- list(
    pseudobulk_seconds = time_pseudobulk_sec,
    signature_seconds  = time_signature_sec,
    total_seconds      = time_total_sec
  )

  if (verbose) {
    cat("\033[33m",
        sprintf("Generated %d pseudobulks (%d cells each) from %d genes\n",
                N, cells_per_bulk, nrow(bulk_expression)),
        "\033[0m", sep = "")
    cat("\033[33m",
        sprintf("Pseudobulk creation time: %.2f secs\n", time_pseudobulk_sec),
        "\033[0m", sep = "")

    if (return_signature_matrix) {
      cat("\033[33m",
          sprintf("Computed bulk_signature_profiles (list of genes\u00d7bulk matrices per cell type)\nSignature matrix time: %.2f secs\n",
                  time_signature_sec),
          "\033[0m", sep = "")
      cat("\033[33m", "Returned global_signature_matrix (genes \u00d7 celltypes)\n", "\033[0m", sep = "")
    }

    cat("\033[33m",
        sprintf("Total quasar_sim_bulk time: %.2f secs\n", time_total_sec),
        "\033[0m", sep = "")
  }

  result
}