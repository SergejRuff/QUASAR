select_top_variable_genes <- function(mat, top_n = 0, top_perc = 0) {

  # If empty input, nothing to do
  if (is.null(dim(mat)) || nrow(mat) == 0L) return(mat)

  # Make sure we have rownames (needed because you subset by names(vars))
  if (is.null(rownames(mat))) rownames(mat) <- as.character(seq_len(nrow(mat)))

  # Only compute variances if we actually filter
  if (top_n > 0 || top_perc > 0) {

    # Robust variance computation
    vars <- tryCatch(
      matrixStats::rowVars(as.matrix(mat)),
      error = function(e) apply(as.matrix(mat), 1, var)
    )

    # Ensure names exist and are usable
    if (is.null(names(vars))) names(vars) <- rownames(mat)

    # Replace NA/Inf with 0 so sorting/subsetting is stable
    vars[!is.finite(vars)] <- 0

    # Sort once
    vars_sorted <- sort(vars, decreasing = TRUE, na.last = TRUE)

    if (top_n > 0) {
      k <- min(as.integer(top_n), length(vars_sorted))
      k <- max(k, 1L)  # never select 0
      keep_genes <- names(vars_sorted)[seq_len(k)]
      mat <- mat[keep_genes, , drop = FALSE]

    } else if (top_perc > 0) {
      # Interpret top_perc as fraction of genes to keep; clamp to (0,1]
      top_perc <- min(as.numeric(top_perc), 1)
      top_count <- ceiling(top_perc * length(vars_sorted))
      top_count <- max(top_count, 1L)  # never select 0
      keep_genes <- names(vars_sorted)[seq_len(top_count)]
      mat <- mat[keep_genes, , drop = FALSE]
    }
  }

  return(mat)
}

normalize_bulks <- function(ref_matrix,
                            target_matrix,
                            top_n    = 0,
                            top_perc = 0.98,
                            normalization = c("samplewise_ss",
                                              "samplewise_minmax",
                                              "samplewise_cpm",
                                              "none")) {
  normalization <- match.arg(normalization)
  ref_matrix    <- as.matrix(ref_matrix)
  target_matrix <- as.matrix(target_matrix)

  ref_filt    <- select_top_variable_genes(ref_matrix,    top_n = top_n, top_perc = top_perc)
  target_filt <- select_top_variable_genes(target_matrix, top_n = top_n, top_perc = top_perc)
  keep_genes  <- intersect(rownames(ref_filt), rownames(target_filt))
  ref_filt    <- ref_filt[keep_genes, , drop = FALSE]
  target_filt <- target_filt[keep_genes, , drop = FALSE]

  ref_raw_filtered    <- ref_filt
  target_raw_filtered <- target_filt

  pop_sd_cols <- function(M) {
    mu     <- colMeans(M)
    M_cent <- sweep(M, 2, mu, "-")
    sqrt(colMeans(M_cent^2))
  }

  normalize_matrix <- function(mat, method) {

    if (method == "none") {
      out <- mat
      rownames(out) <- rownames(mat)
      colnames(out) <- colnames(mat)
      return(out)
    }
    if (method == "samplewise_ss") {
      mat_log     <- log1p(mat)
      sample_mean <- colMeans(mat_log)
      sample_sd   <- pop_sd_cols(mat_log)
      sd_safe     <- ifelse(sample_sd == 0, 1, sample_sd)
      out <- sweep(sweep(mat_log, 2, sample_mean, "-"), 2, sd_safe, "/")
      out[, sample_sd == 0] <- 0
      out[!is.finite(out)]  <- 0
      rownames(out) <- rownames(mat); colnames(out) <- colnames(mat)
      return(out)
    }
    if (method == "samplewise_minmax") {
      out <- apply(mat, 2, function(x) {
        rng <- range(x, na.rm = TRUE)
        if (!is.finite(rng[1]) || !is.finite(rng[2]) || diff(rng) == 0)
          return(rep(0, length(x)))
        (x - rng[1]) / (rng[2] - rng[1])
      })
      out <- as.matrix(out)
      rownames(out) <- rownames(mat); colnames(out) <- colnames(mat)
      return(out)
    }
    if (method == "samplewise_cpm") {
      lib_sizes <- colSums(mat, na.rm = TRUE)
      out <- sweep(mat, 2, lib_sizes, FUN = "/") * 1e6
      out[, lib_sizes == 0] <- 0
      out <- as.matrix(out)
      rownames(out) <- rownames(mat); colnames(out) <- colnames(mat)
      return(out)
    }
    stop("Unknown normalization method: ", method)
  }

  list(
    pseudobulk_norm      = normalize_matrix(ref_filt,    normalization),
    target_norm          = normalize_matrix(target_filt, normalization),
    pseudobulk_raw       = ref_raw_filtered,
    target_raw           = target_raw_filtered,
    genes                = keep_genes,
    normalization_method = normalization
  )
}