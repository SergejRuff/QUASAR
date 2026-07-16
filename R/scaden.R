#' Simulate pseudobulk training data for Scaden
#'
#' Faithful port of `generate_simulated_data()`. Counts are sampled with
#' replacement within each cell type according to Dirichlet-generated mixture
#' proportions. Optionally merges a set of cell types into a single unknown
#' class, following Scaden's post-simulation merge semantics.
#'
#' @param sc_data A `Seurat`, `SingleCellExperiment`, matrix, `data.frame`, or
#'   sparse `Matrix` reference.
#' @param metadata Optional metadata for matrix-like input.
#' @param celltypes Optional character vector of cell-type labels for
#'   matrix-like input. Used only when `metadata` is not supplied.
#' @param celltype_col Character scalar giving the metadata column holding
#'   cell-type labels.
#' @param assay Character scalar giving the Seurat assay name.
#' @param slot Character scalar giving the assay slot/layer (Seurat) or assay
#'   name (SingleCellExperiment).
#' @param d_prior Optional numeric vector of Dirichlet concentration
#'   parameters. If `NULL`, a symmetric prior of ones is used. Length must equal
#'   the number of cell types *before* unknown merging.
#' @param n Integer scalar. Target cells per pseudobulk before flooring.
#' @param samplenum Integer scalar. Number of pseudobulks to generate.
#' @param seed Optional integer scalar for reproducible simulation.
#' @param sparse Logical scalar. Zero out a subset of cell types in a subset of
#'   samples.
#' @param sparse_prob Numeric scalar in `[0, 1)`. Controls both the proportion
#'   of sparse samples and the proportion of cell types zeroed within them.
#' @param rare Logical scalar. Perturb a subset of cell types to very small
#'   fractions in a subset of samples.
#' @param rare_percentage Numeric scalar in `[0, 1]`. Fraction of cell types
#'   treated as rare when `rare = TRUE`.
#' @param unknown_celltypes Character vector of cell-type labels to merge into a
#'   single unknown class. Defaults to `character(0)` (no merging).
#' @param unknown_label Character scalar used as the merged label.
#' @param verbose Logical scalar. Print progress.
#'
#' @return A named list with `train_x` (samples x genes), `train_y`
#'   (samples x cell types), `cell_num`, `celltypes`, `genes`, `elapsed`, and
#'   `settings`.
#'
#' @details
#' Unknown merging happens **after** simulation: each original cell type gets
#' its own Dirichlet component and is sampled independently, and the proportion
#' columns of the unknown types are then summed into one column. This mirrors
#' Scaden's `y.groupby(y.columns, axis=1).sum()` and differs from relabelling
#' cells before simulation, which would give the merged class a single Dirichlet
#' component and pool its cells uniformly.
#'
#' @importFrom methods slot slotNames as
#' @importFrom SummarizedExperiment assay assayNames colData
#' @importFrom Matrix sparseMatrix
#' @importFrom MCMCpack rdirichlet
#' @importFrom stats runif
#' @importFrom utils flush.console
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
  unknown_celltypes = character(0),
  unknown_label = "Unknown",
  verbose = TRUE
) {

  start_time <- Sys.time()

  fmt_time <- function() {
    secs <- max(0L, as.integer(round(difftime(Sys.time(), start_time, units = "secs"))))
    sprintf("%02d:%02d:%02d", secs %/% 3600L, (secs %% 3600L) %/% 60L, secs %% 60L)
  }

  draw_progress <- function(done, total, width = 32L) {
    frac  <- if (total > 0L) min(1, done / total) else 1
    nfill <- as.integer(round(frac * width))
    cat(sprintf("\r  |%s%s| %3d%%  %d/%d  %s",
                strrep("\u2588", nfill),
                strrep("\u2591", width - nfill),
                as.integer(round(100 * frac)),
                done, total, fmt_time()))
    utils::flush.console()
  }

  require_celltype_col <- function(meta) {
    if (is.null(celltype_col) || !nzchar(celltype_col))
      stop("`celltype_col` must be provided explicitly.")
    if (!celltype_col %in% colnames(meta))
      stop(sprintf("`celltype_col` '%s' not found in metadata.", celltype_col))
  }

  if (!is.null(seed)) {
    stopifnot("`seed` must be a single finite number" =
                length(seed) == 1L && is.finite(seed))
    set.seed(as.integer(seed))
  }

  if (verbose) cat("Simulating pseudobulks for Scaden ...\n")


  if (inherits(sc_data, "Seurat")) {
    stopifnot("Seurat object missing requested assay" = assay %in% names(sc_data@assays))
    assay_obj <- sc_data[[assay]]
    avail <- if (inherits(assay_obj, "Assay5"))
               SeuratObject::Layers(assay_obj)
             else
               methods::slotNames(assay_obj)
    stopifnot("Slot/layer not found in assay" = slot %in% avail)
    counts <- SeuratObject::GetAssayData(sc_data, assay = assay, layer = slot)
    meta   <- sc_data@meta.data
    require_celltype_col(meta)
    cell_types <- as.character(meta[[celltype_col]])

  } else if (inherits(sc_data, "SingleCellExperiment")) {
    avail <- SummarizedExperiment::assayNames(sc_data)
    stopifnot("Assay not found in SCE object" = slot %in% avail)
    counts <- SummarizedExperiment::assay(sc_data, slot)
    meta   <- as.data.frame(SummarizedExperiment::colData(sc_data))
    require_celltype_col(meta)
    cell_types <- as.character(meta[[celltype_col]])

  } else if (is.matrix(sc_data) || is.data.frame(sc_data) || inherits(sc_data, "Matrix")) {
    x <- if (is.data.frame(sc_data)) as.matrix(sc_data) else sc_data

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
    stop("Unsupported `sc_data` type. Expected Seurat, SingleCellExperiment, ",
         "matrix, data.frame, or Matrix.")
  }

  counts <- methods::as(counts, "CsparseMatrix")
  if (is.null(rownames(counts))) rownames(counts) <- paste0("gene_", seq_len(nrow(counts)))
  if (is.null(colnames(counts))) colnames(counts) <- paste0("cell_", seq_len(ncol(counts)))
  stopifnot("Cell-type vector length != number of cells" =
              length(cell_types) == ncol(counts))


  na_cells <- which(is.na(cell_types))
  if (length(na_cells) > 0L) {
    counts     <- counts[, -na_cells, drop = FALSE]
    cell_types <- cell_types[-na_cells]
  }


  ct_levels <- sort(unique(cell_types), method = "radix")
  n_ct      <- length(ct_levels)
  if (n_ct < 1L) stop("No cell types available.")

  stopifnot(
    "`n` must be a positive number"         = is.numeric(n) && length(n) == 1L && n > 0,
    "`samplenum` must be a positive number" = is.numeric(samplenum) && length(samplenum) == 1L && samplenum > 0
  )
  n         <- as.integer(n)
  samplenum <- as.integer(samplenum)

  if (is.null(d_prior)) {
    d_prior <- rep(1, n_ct)
  } else {
    stopifnot(
      "`d_prior` length must equal number of cell types" = length(d_prior) == n_ct,
      "All `d_prior` entries must be > 0"                = all(d_prior > 0)
    )
  }
  if (sparse) stopifnot("`sparse_prob` must be in [0, 1)" = sparse_prob >= 0 && sparse_prob < 1)
  if (rare)   stopifnot("`rare_percentage` must be in [0, 1]" = rare_percentage >= 0 && rare_percentage <= 1)

  unknown_celltypes <- as.character(unknown_celltypes)
  stopifnot("`unknown_label` must be a single string" =
              is.character(unknown_label) && length(unknown_label) == 1L)

  if (verbose)
    cat(sprintf("  Reference: %d cells | %d genes | %d cell types | %d pseudobulks\n",
                ncol(counts), nrow(counts), n_ct, samplenum))


  prop <- MCMCpack::rdirichlet(samplenum, d_prior)
  prop <- prop / rowSums(prop)


  if (sparse) {
    n_sparse <- as.integer(samplenum * sparse_prob)
    # The min() clamp has no Python counterpart but can never bind:
    # int(K * 0.5) <= K - 1 for all K >= 2.
    n_zero   <- min(as.integer(n_ct * sparse_prob), n_ct - 1L)
    if (n_sparse > 0L && n_zero > 0L) {
      for (i in seq_len(n_sparse)) {
        prop[i, sample.int(n_ct, n_zero)] <- 0
      }
      prop <- prop / rowSums(prop)
    }
  }



  if (rare) {
    n_rare_ct <- as.integer(n_ct * rare_percentage)
    if (n_rare_ct > 0L) {
      set.seed(0L)   # Python hard-codes np.random.seed(0) here
      rare_idx <- sample.int(n_ct, n_rare_ct)
      prop     <- prop / rowSums(prop)
      n_rows   <- min(as.integer(0.5 * samplenum) +
                        as.integer(as.integer(rare_percentage * 0.5 * samplenum)),
                      samplenum)
      for (i in seq_len(n_rows)) {
        buf <- stats::runif(n_rare_ct, 0, 0.03)
        prop[i, rare_idx] <- 0
        remain <- sum(prop[i, ])
        if (remain > 0) prop[i, ] <- (1 - sum(buf)) * prop[i, ] / remain
        prop[i, rare_idx] <- buf
      }
    }
  }

  # -- Cell counts and realized proportions -----------------------------------

  cell_num <- matrix(as.integer(floor(n * prop)),
                     nrow = samplenum, ncol = n_ct,
                     dimnames = list(NULL, ct_levels))

  # No Python counterpart, and provably inert: rows of `prop` sum to 1, so
  # max(prop[i, ]) >= 1/n_ct, and floor(n * prop) can only be all-zero if every
  # proportion is < 1/n -- impossible for n_ct <= n.
  empty <- rowSums(cell_num) == 0L
  if (any(empty)) {
    best_col <- max.col(prop[empty, , drop = FALSE], ties.method = "first")
    for (ii in seq_along(best_col)) cell_num[which(empty)[ii], best_col[ii]] <- 1L
  }

  prop_realized <- cell_num / rowSums(cell_num)
  colnames(prop_realized) <- ct_levels

  # -- Sample cells -----------------------------------------------------------

  ct_indices <- split(seq_along(cell_types), factor(cell_types, levels = ct_levels))

  ids_all <- vector("list", n_ct)
  grp_all <- vector("list", n_ct)

  for (k in seq_along(ct_levels)) {
    total_k <- sum(cell_num[, k])
    pool    <- ct_indices[[k]]
    if (total_k > 0L && length(pool) > 0L) {
      # Index into `pool` rather than calling sample(pool, ...) directly:
      # sample() on a length-1 numeric silently means sample.int(pool).
      ids_all[[k]] <- pool[sample.int(length(pool), total_k, replace = TRUE)]
      grp_all[[k]] <- rep.int(seq_len(samplenum), times = cell_num[, k])
    }
  }

  ids <- unlist(ids_all, use.names = FALSE)
  grp <- unlist(grp_all, use.names = FALSE)

  S <- Matrix::sparseMatrix(i = ids, j = grp, x = 1,
                            dims = c(ncol(counts), samplenum))

  gene_names   <- rownames(counts)
  sample_names <- paste0("sample_", seq_len(samplenum))

  train_x <- matrix(0, nrow = samplenum, ncol = nrow(counts),
                    dimnames = list(sample_names, gene_names))

  block_size <- max(1L, ceiling(samplenum / 100L))
  starts     <- seq.int(1L, samplenum, by = block_size)

  if (verbose) {
    cat("  Generating pseudobulks:\n")
    draw_progress(0L, samplenum)
  }

  for (s0 in starts) {
    cols  <- s0:min(s0 + block_size - 1L, samplenum)
    block <- as.matrix(counts %*% S[, cols, drop = FALSE])
    train_x[cols, ] <- t(block)
    if (verbose) draw_progress(cols[length(cols)], samplenum)
  }

  if (verbose) cat("\n")

  train_y <- prop_realized

  # -- Merge unknown cell types (post-simulation, as in Scaden) ---------------

  if (length(unknown_celltypes) > 0L) {
    hit <- ct_levels %in% unknown_celltypes
    if (any(hit)) {
      new_names <- ct_levels
      new_names[hit] <- unknown_label

      if (unknown_label %in% ct_levels[!hit])
        warning(sprintf(
          "`unknown_label` '%s' is also an existing cell type; it will be merged too.",
          unknown_label))

      # pandas groupby(axis=1).sum() sorts the resulting columns.
      grp_ct  <- factor(new_names, levels = sort(unique(new_names), method = "radix"))
      train_y  <- t(rowsum(t(train_y),  grp_ct, reorder = TRUE))
      cell_num <- t(rowsum(t(cell_num), grp_ct, reorder = TRUE))
      ct_levels <- levels(grp_ct)

      if (verbose)
        cat(sprintf("  Merged %d cell type(s) into '%s' (%d cell types remain)\n",
                    sum(hit), unknown_label, length(ct_levels)))
    }
  }

  rownames(train_y)  <- sample_names
  colnames(train_y)  <- ct_levels
  rownames(cell_num) <- sample_names
  colnames(cell_num) <- ct_levels

  elapsed <- fmt_time()
  if (verbose)
    cat(sprintf("  Done in %s  (%d pseudobulks, %d cell types)\n",
                elapsed, samplenum, length(ct_levels)))

  list(
    train_x   = train_x,
    train_y   = train_y,
    cell_num  = cell_num,
    celltypes = ct_levels,
    genes     = gene_names,
    elapsed   = elapsed,
    settings  = list(
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
      assay             = assay,
      slot              = slot,
      celltype_col      = celltype_col
    )
  )
}


#' Process simulated and bulk data for Scaden
#'
#' Aligns simulated pseudobulk training data with real bulk data, filters genes,
#' log-transforms, and applies per-sample scaling. Two gene-filtering strategies
#' are available via `mode`.
#'
#' @param sim_data A list returned by [scaden_sim_pb()], with `train_x` and
#'   `train_y`.
#' @param bulk_data Bulk expression matrix or `data.frame`. Gene names must
#'   appear in row or column names; orientation is inferred.
#' @param mode Character scalar. `"tape"` (default) reproduces TAPE's
#'   `ProcessInputData`; `"scaden"` reproduces original Scaden's preprocessing.
#' @param variance_threshold Numeric scalar in `[0, 1)`. `mode = "tape"` only.
#'   Fraction of genes defining the variance cutoff, applied separately to each
#'   matrix.
#' @param scaler Character scalar. `mode = "tape"` only. `"mms"` for per-sample
#'   min-max, `"ss"` for per-sample standardization.
#' @param var_cutoff Numeric scalar. `mode = "scaden"` only. Genes in
#'   `bulk_data` with variance at or below this value are dropped. Set to `NULL`
#'   to skip.
#' @param top_var_genes Optional integer scalar. `mode = "scaden"` only. Keep
#'   only the most variable genes after `var_cutoff`.
#' @param verbose Logical scalar. Print progress.
#'
#' @return A named list with `train_x`, `train_y`, `test_x`, `genename`,
#'   `celltypes`, `samplename`, `mode`, and `elapsed`.
#'
#' @details
#' The two modes differ only in gene filtering and scaler choice:
#'
#' \describe{
#'   \item{`"tape"`}{Top-quantile variance cutoff applied to **both** matrices,
#'   then intersect. `log(x + 1)`, then `mms` or `ss`.}
#'   \item{`"scaden"`}{Absolute variance cutoff applied to the **bulk matrix
#'   only**, then intersect with the training genes. `log2(x + 1)`, then
#'   per-sample min-max (`log_min_max`); `scaler` is ignored.}
#' }
#'
#' The log base is in fact immaterial: both scalers are invariant to
#' multiplication by a positive constant, so `log2` and `log` give identical
#' output up to floating-point rounding. The bases are kept literal to their
#' respective sources.
#'
#' Both scalers operate per sample (`fit_transform(x.T).T` in sklearn), use
#' population statistics (`ddof = 0`), and map constant samples to exactly zero,
#' matching sklearn's `_handle_zeros_in_scale`. The `"tape"` variance filter uses
#' `ddof = 1` to match `pandas.DataFrame.var`.
#'
#' @importFrom matrixStats colVars rowMins rowMaxs
#' @importFrom stats setNames
#'
#' @export
scaden_process <- function(sim_data,
                           bulk_data,
                           mode = c("tape", "scaden"),
                           variance_threshold = 0.98,
                           scaler = c("mms", "ss"),
                           var_cutoff = 0.1,
                           top_var_genes = NULL,
                           verbose = TRUE) {

  mode        <- match.arg(mode)
  scaler_given <- !missing(scaler)
  scaler      <- match.arg(scaler)
  start_time  <- Sys.time()

  if (identical(mode, "scaden")) {
    if (scaler_given && !identical(scaler, "mms"))
      warning("Original Scaden offers only `log_min_max`; `scaler` is ignored ",
              "when `mode = 'scaden'`.")
    scaler <- "mms"
  }

  fmt_time <- function() {
    secs <- max(0L, as.integer(round(difftime(Sys.time(), start_time, units = "secs"))))
    sprintf("%02d:%02d:%02d", secs %/% 3600L, (secs %% 3600L) %/% 60L, secs %% 60L)
  }

  orient_to_samples_x_genes <- function(mat, target_genes, label) {
    rn <- rownames(mat); cn <- colnames(mat)
    hit_rows <- if (!is.null(rn)) length(intersect(rn, target_genes)) else 0L
    hit_cols <- if (!is.null(cn)) length(intersect(cn, target_genes)) else 0L
    if (hit_cols > hit_rows) return(mat)
    if (hit_rows > hit_cols) return(t(mat))
    if (!is.null(rn) && nrow(mat) == length(target_genes)) return(t(mat))
    if (!is.null(cn) && ncol(mat) == length(target_genes)) return(mat)
    stop(sprintf(paste0("Cannot infer orientation of `%s`. Ensure gene names ",
                        "appear in row or column names."), label))
  }

  # TAPE: var_cutoff = x.var(axis=0).sort_values(ascending=False)[int(p * ncol)]
  # The Python index is 0-based, so R needs + 1L.
  quantile_var_filter <- function(M, label) {
    idx <- as.integer(ncol(M) * variance_threshold) + 1L
    if (idx > ncol(M))
      stop(sprintf("`variance_threshold` too high for %s: index %d exceeds %d genes.",
                   label, idx, ncol(M)))
    v <- stats::setNames(matrixStats::colVars(M), colnames(M))  # ddof = 1
    M[, v > sort(v, decreasing = TRUE)[idx], drop = FALSE]
  }


  row_ss <- function(M) {
    mu      <- rowMeans(M)
    sdv     <- sqrt(rowMeans((M - mu)^2))       # ddof = 0
    sd_safe <- ifelse(sdv == 0, 1, sdv)         # _handle_zeros_in_scale
    out     <- sweep(sweep(M, 1, mu, "-"), 1, sd_safe, "/")
    out[sdv == 0, ] <- 0
    out
  }


  row_mms <- function(M) {
    mn       <- matrixStats::rowMins(M)
    rng      <- matrixStats::rowMaxs(M) - mn
    rng_safe <- ifelse(rng == 0, 1, rng)        # _handle_zeros_in_scale
    out      <- sweep(sweep(M, 1, mn, "-"), 1, rng_safe, "/")
    out[rng == 0, ] <- 0
    out
  }

  scale_rows <- function(M) if (scaler == "ss") row_ss(M) else row_mms(M)



  if (!is.list(sim_data) || is.null(sim_data$train_x) || is.null(sim_data$train_y))
    stop("`sim_data` must be a list with `train_x` and `train_y`.")

  train_x <- as.matrix(sim_data$train_x)
  train_y <- as.matrix(sim_data$train_y)

  if (is.null(colnames(train_x)))
    stop("Training matrix must have gene names as column names.")
  if (is.null(rownames(train_x)))
    rownames(train_x) <- paste0("sample_", seq_len(nrow(train_x)))
  if (is.null(rownames(train_y))) rownames(train_y) <- rownames(train_x)
  if (is.null(colnames(train_y)))
    colnames(train_y) <- paste0("celltype_", seq_len(ncol(train_y)))

  test_x <- orient_to_samples_x_genes(as.matrix(bulk_data), colnames(train_x), "bulk_data")
  if (is.null(colnames(test_x)))
    stop("`bulk_data` must have gene names in row or column names.")
  if (is.null(rownames(test_x)))
    rownames(test_x) <- paste0("sample_", seq_len(nrow(test_x)))

  if (verbose)
    cat(sprintf("Processing Scaden data (mode = %s, scaler = %s)\n  train %d x %d | bulk %d x %d\n",
                mode, scaler, nrow(train_x), ncol(train_x), nrow(test_x), ncol(test_x)))


  if (identical(mode, "tape")) {
    # Quantile cutoff applied independently to each matrix.
    train_x <- quantile_var_filter(train_x, "training data")
    test_x  <- quantile_var_filter(test_x,  "bulk data")
    bulk_genes <- colnames(test_x)

  } else {
    # Scaden: absolute cutoff on the bulk matrix only.
    gene_var <- stats::setNames(matrixStats::colVars(test_x), colnames(test_x))

    if (!is.null(var_cutoff)) {
      before   <- length(gene_var)
      gene_var <- gene_var[gene_var > var_cutoff]
      if (verbose)
        cat(sprintf("  variance > %g:      %d \u2192 %d\n", var_cutoff, before, length(gene_var)))
    }
    if (!is.null(top_var_genes) && length(gene_var) > 0L) {
      before   <- length(gene_var)
      top_n    <- min(as.integer(top_var_genes), length(gene_var))
      gene_var <- sort(gene_var, decreasing = TRUE)[seq_len(top_n)]
      if (verbose)
        cat(sprintf("  top variable genes: %d \u2192 %d\n", before, top_n))
    }
    bulk_genes <- names(gene_var)
    test_x     <- test_x[, bulk_genes, drop = FALSE]
  }


  inter <- intersect(colnames(train_x), bulk_genes)
  if (length(inter) == 0L)
    stop("No genes remain after intersection - check that gene names match.")

  train_x <- train_x[, inter, drop = FALSE]
  test_x  <- test_x[,  inter, drop = FALSE]

  if (verbose) cat(sprintf("  intersected genes:  %d\n", length(inter)))

  genename   <- inter
  celltypes  <- colnames(train_y)
  samplename <- rownames(test_x)


  if (identical(mode, "tape")) {
    train_x <- log(train_x + 1)
    test_x  <- log(test_x + 1)
  } else {
    train_x <- log2(train_x + 1)
    test_x  <- log2(test_x + 1)
  }


  train_x <- scale_rows(train_x)
  test_x  <- scale_rows(test_x)

  if (verbose) cat(sprintf("Done in %s\n", fmt_time()))

  list(
    train_x    = train_x,
    train_y    = train_y,
    test_x     = test_x,
    genename   = genename,
    celltypes  = celltypes,
    samplename = samplename,
    mode       = mode,
    elapsed    = fmt_time()
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
#' [scaden_predict()].
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

  # Plain (default terminal color), consistent with the progress bar below.
  cat(scaden_ascii, "\n", sep = "")

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

    # Block-style bar matching scaden_sim_pb's draw_progress (plain color):
    # solid blocks for progress, light blocks for the remainder.
    make_bar <- function(i, n, width = 32L) {
      frac  <- if (n > 0) min(1, i / n) else 1
      nfill <- as.integer(round(frac * width))
      paste0(strrep("\u2588", nfill), strrep("\u2591", width - nfill))
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

      pct <- if (epochs > 0) as.integer(round(100 * i / epochs)) else 100L

      cat(sprintf(
        "\r%-9s |%s| %3d%% | %3d/%3d | loss %.6f | elapsed %s | ETA %s",
        model_name,
        make_bar(i, epochs, width = 32L),
        pct,
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