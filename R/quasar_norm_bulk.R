#' Select top variable genes
#'
#' Internal helper used by [normalize_bulks()] to retain the most variable
#' genes based on row variance.
#'
#' @keywords internal
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


#' Normalize reference and target bulk matrices
#'
#' Selects the most variable genes shared between a reference (pseudobulk)
#' matrix and a target matrix, then applies sample-wise normalization.
#'
#' The function first filters genes based on their variance within each
#' matrix. Only genes present in both filtered matrices are retained. The filtered matrices can
#' then be normalized using one of several sample-wise approaches.
#'
#' @param ref_matrix Numeric matrix containing the reference/pseudobulk data.
#'   Rows correspond to genes and columns correspond to samples.
#' @param target_matrix Numeric matrix containing the target bulk data.
#'   Rows correspond to genes and columns correspond to samples.
#' @param top_n Integer. Number of most variable genes to retain. If greater
#'   than zero, this takes precedence over `top_perc`.
#' @param top_perc Numeric between 0 and 1. Fraction of the most variable
#'   genes to retain. Default is 0.98.
#' @param normalization Character specifying the normalization method:
#'   \describe{
#'     \item{"samplewise_ss"}{Log1p transformation followed by sample-wise
#'     standardization (mean and population standard deviation).}
#'     \item{"samplewise_minmax"}{Sample-wise min-max scaling.}
#'     \item{"samplewise_cpm"}{Counts per million normalization per sample.}
#'     \item{"none"}{No normalization.}
#'   }
#'
#' @return A list containing:
#' \describe{
#'   \item{pseudobulk_norm}{Normalized filtered reference matrix.}
#'   \item{target_norm}{Normalized filtered target matrix.}
#'   \item{pseudobulk_raw}{Filtered reference matrix before normalization.}
#'   \item{target_raw}{Filtered target matrix before normalization.}
#'   \item{genes}{Genes retained after filtering and intersection.}
#'   \item{normalization_method}{Selected normalization method.}
#' }
#'
#' @details
#' Genes are ranked by variance and only the most variable genes are retained.
#' Filtering is performed independently for the reference and target matrices,
#' followed by intersection of retained genes to ensure both matrices contain
#' identical features.
#'
#' @examples
#' # Create example reference and target matrices
#' set.seed(1)
#'
#' ref <- matrix(
#'   rpois(300, lambda = 10),
#'   nrow = 100,
#'   ncol = 3,
#'   dimnames = list(
#'     paste0("gene", 1:100),
#'     paste0("ref_sample", 1:3)
#'   )
#' )
#'
#' target <- matrix(
#'   rpois(300, lambda = 15),
#'   nrow = 100,
#'   ncol = 3,
#'   dimnames = list(
#'     paste0("gene", 1:100),
#'     paste0("target_sample", 1:3)
#'   )
#' )
#'
#' result <- normalize_bulks(
#'   ref_matrix = ref,
#'   target_matrix = target,
#'   top_perc = 0.5,
#'   normalization = "samplewise_ss"
#' )
#'
#' # Normalized matrices
#' dim(result$pseudobulk_norm)
#' dim(result$target_norm)
#'
#' # Retained genes
#' head(result$genes)
#'
#' @export
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