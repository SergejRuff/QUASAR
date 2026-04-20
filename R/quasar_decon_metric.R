#' Calculate proportion-estimation performance metrics
#'
#' Computes multiple evaluation metrics for estimated cell-type proportions
#' relative to ground-truth proportions after aligning both matrices by common
#' sample names and common cell-type names.
#'
#' Both inputs must be numeric matrices or coercible to numeric matrices, with
#' samples in rows and cell types in columns. Only the intersecting samples and
#' cell types are retained for metric calculation.
#'
#' @param estimated_proportions Numeric matrix or data.frame of estimated
#'   cell-type proportions, with samples in rows and cell types in columns.
#' @param ground_truth_proportions Numeric matrix or data.frame of ground-truth
#'   cell-type proportions, with samples in rows and cell types in columns.
#' @param epsilon Numeric scalar. Small constant used to stabilize correlation
#'   and divergence calculations for nearly constant vectors and zero entries.
#'
#' @return A named list with the following components:
#' \describe{
#'   \item{cell_type_rmse}{Named numeric vector of root mean squared error
#'   (RMSE) for each cell type across samples.}
#'   \item{cell_type_mad}{Named numeric vector of median absolute deviation
#'   (MAD) of the cell-type-specific error vector across samples.}
#'   \item{cell_type_nmae}{Named numeric vector of normalized mean absolute
#'   error (NMAE) for each cell type across samples.}
#'   \item{mean_ground_truth_for_nmae}{Named numeric vector of mean ground-truth
#'   proportions for each cell type.}
#'   \item{pearson_celltype_cor}{Named numeric vector of Pearson correlations
#'   between estimated and true proportions for each cell type across samples.}
#'   \item{spearman_celltype_cor}{Named numeric vector of Spearman correlations
#'   between estimated and true proportions for each cell type across samples.}
#'   \item{per_celltype_jsd}{Named numeric vector of Jensen-Shannon distances
#'   for each cell type, where the distribution is formed across samples.}
#'   \item{per_sample_jsd}{Named numeric vector of Jensen-Shannon distances for
#'   each sample, where the distribution is formed across cell types.}
#' }
#'
#' @details
#' Let \eqn{\hat{p}_{ic}} denote the estimated proportion for sample \eqn{i} and
#' cell type \eqn{c}, and let \eqn{p_{ic}} denote the corresponding ground-truth
#' proportion. After restricting both inputs to common samples and cell types,
#' the following metrics are calculated.
#'
#' \strong{1. Root mean squared error (RMSE) per cell type}
#'
#' For each cell type \eqn{c}, RMSE is calculated across all shared samples:
#' \deqn{
#' \mathrm{RMSE}_c = \sqrt{\frac{1}{n}\sum_{i=1}^{n}(\hat{p}_{ic} - p_{ic})^2}
#' }
#' where \eqn{n} is the number of shared samples.
#'
#' \strong{2. Median absolute deviation (MAD) of the error per cell type}
#'
#' For each cell type \eqn{c}, the error vector is first defined as
#' \eqn{e_{ic} = \hat{p}_{ic} - p_{ic}}. The reported MAD is then:
#' \deqn{
#' \mathrm{MAD}_c = \mathrm{median}_i\left(\left|e_{ic} -
#' \mathrm{median}_j(e_{jc})\right|\right)
#' }
#' This quantifies the spread of the cell-type-specific errors across samples.
#'
#' \strong{3. Normalized mean absolute error (NMAE) per cell type}
#'
#' For each cell type \eqn{c}, the mean absolute error is normalized by the
#' range of the ground-truth proportions across samples:
#' \deqn{
#' \mathrm{NMAE}_c =
#' \frac{1}{n}\sum_{i=1}^{n}
#' \frac{|\hat{p}_{ic} - p_{ic}|}{\max_i(p_{ic}) - \min_i(p_{ic})}
#' }
#' when the ground-truth range is non-zero.
#'
#' If the ground-truth range for a cell type is zero, the function uses the
#' following fallback:
#' \itemize{
#'   \item returns \eqn{0} if estimated and true values are identical for all
#'   samples,
#'   \item otherwise returns the unnormalized mean absolute error.
#' }
#'
#' \strong{4. Mean ground-truth proportion per cell type}
#'
#' For each cell type \eqn{c}, the mean ground-truth proportion is:
#' \deqn{
#' \bar{p}_c = \frac{1}{n}\sum_{i=1}^{n} p_{ic}
#' }
#'
#' \strong{5. Pearson correlation per cell type}
#'
#' For each cell type \eqn{c}, the Pearson correlation is calculated between the
#' estimated and true proportions across samples:
#' \deqn{
#' r_c = \mathrm{cor}(\hat{p}_{\cdot c}, p_{\cdot c})
#' }
#'
#' If the estimated or true vector has near-zero standard deviation, random
#' noise with standard deviation `epsilon` is added before correlation
#' calculation to avoid undefined results for constant vectors.
#'
#' \strong{6. Spearman correlation per cell type}
#'
#' For each cell type \eqn{c}, the Spearman rank correlation is calculated
#' across samples:
#' \deqn{
#' \rho_c = \mathrm{cor}(\hat{p}_{\cdot c}, p_{\cdot c},
#' \mathrm{method} = "spearman")
#' }
#'
#' As for Pearson correlation, a small perturbation is added when one of the
#' vectors has near-zero standard deviation.
#'
#' \strong{7. Jensen-Shannon distance (JSD) per cell type}
#'
#' For each cell type \eqn{c}, the function compares the distribution of that
#' cell type across samples between ground truth and estimation. Let
#' \eqn{P_c = (p_{1c}, \dots, p_{nc})} and
#' \eqn{Q_c = (\hat{p}_{1c}, \dots, \hat{p}_{nc})}. These vectors are first
#' normalized to sum to 1:
#' \deqn{
#' P_c^* = \frac{P_c}{\sum_i p_{ic}}, \qquad
#' Q_c^* = \frac{Q_c}{\sum_i \hat{p}_{ic}}
#' }
#'
#' The midpoint distribution is:
#' \deqn{
#' M_c = \frac{1}{2}(P_c^* + Q_c^*)
#' }
#'
#' The Kullback-Leibler divergence is computed as:
#' \deqn{
#' \mathrm{KL}(P \parallel Q) = \sum_k P_k \log\left(\frac{P_k}{Q_k}\right)
#' }
#' with `epsilon` added internally to both \eqn{P} and \eqn{Q} for numerical
#' stability.
#'
#' The reported Jensen-Shannon distance is:
#' \deqn{
#' \mathrm{JSD}_c =
#' \sqrt{
#' \frac{1}{2}\mathrm{KL}(P_c^* \parallel M_c) +
#' \frac{1}{2}\mathrm{KL}(Q_c^* \parallel M_c)
#' }
#' }
#'
#' \strong{8. Jensen-Shannon distance (JSD) per sample}
#'
#' For each sample \eqn{i}, the function compares the distribution of cell-type
#' proportions across cell types between ground truth and estimation. Let
#' \eqn{P_i = (p_{i1}, \dots, p_{im})} and
#' \eqn{Q_i = (\hat{p}_{i1}, \dots, \hat{p}_{im})}, where \eqn{m} is the number
#' of shared cell types. After normalization to sum to 1, the same
#' Jensen-Shannon distance formula is applied:
#' \deqn{
#' \mathrm{JSD}_i =
#' \sqrt{
#' \frac{1}{2}\mathrm{KL}(P_i^* \parallel M_i) +
#' \frac{1}{2}\mathrm{KL}(Q_i^* \parallel M_i)
#' }
#' }
#' with \eqn{M_i = \frac{1}{2}(P_i^* + Q_i^*)}.
#'
#' @examples
#' \dontrun{
#' metrics <- quasar_prop_metrics(
#'   estimated_proportions = pred_mat,
#'   ground_truth_proportions = truth_mat
#' )
#'
#' metrics$cell_type_rmse
#' metrics$pearson_celltype_cor
#' metrics$per_sample_jsd
#' }
#'
#' @author Sergej Ruff
#' @export
quasar_prop_metrics <- function(estimated_proportions, ground_truth_proportions, epsilon = 1e-12) {
  
 
  ground_truth_proportions <- as.matrix(ground_truth_proportions)
  estimated_proportions <- as.matrix(estimated_proportions)
  
  if (!is.numeric(ground_truth_proportions)) stop("Ground truth must be numeric")
  if (!is.numeric(estimated_proportions)) stop("Estimated proportions must be numeric")
  

  common_samples <- intersect(rownames(estimated_proportions), rownames(ground_truth_proportions))
  common_celltypes <- intersect(colnames(estimated_proportions), colnames(ground_truth_proportions))
  
  if (length(common_celltypes) == 0) {
    stop("No common cell types found.\n",
         "Estimated has: ", paste(colnames(estimated_proportions), collapse=", "), "\n",
         "Ground truth has: ", paste(colnames(ground_truth_proportions), collapse=", "))
  }
  if (length(common_samples) == 0) {
    stop("No common samples between estimated and ground truth proportions")
  }
  

  estimated_proportions <- estimated_proportions[common_samples, common_celltypes, drop = FALSE]
  ground_truth_proportions <- ground_truth_proportions[common_samples, common_celltypes, drop = FALSE]
  

  q_pm_check_input_type(estimated_proportions, colnames(ground_truth_proportions), 2)
  q_pm_check_input_type(ground_truth_proportions, colnames(estimated_proportions), 2)
  

  
  # RMSE per cell type
  cell_type_rmse <- sqrt(colMeans((estimated_proportions - ground_truth_proportions)^2))
  

  cell_type_mad <- apply(abs(estimated_proportions - ground_truth_proportions), 2,
                         function(x) median(abs(x - median(x))))
  

  mean_ground_truth <- colMeans(ground_truth_proportions)
  
  # Normalized MAE per cell type
  compute_nmae <- function(Y, X) {
    if (length(Y) != length(X)) stop("Vectors must be of equal length.")
    na_mask <- !(is.na(X) | is.na(Y))
    X <- X[na_mask]
    Y <- Y[na_mask]
    n <- length(X)
    if (n == 0) stop("No valid (non-NA) pairs remain for NMAE calculation.")
    range_X <- max(X) - min(X)
    if (range_X == 0) {
      if (all(X == Y)) {
        return(0)
      } else {
        return(mean(abs(Y - X))) # fallback when GT is constant
      }
    } else {
      return(mean(abs(Y - X) / range_X))
    }
  }
  
  cell_type_nmae <- sapply(common_celltypes, function(ct) {
    compute_nmae(estimated_proportions[, ct], ground_truth_proportions[, ct])
  })
  
  # Pearson correlation per cell type
  pearson_celltype_cor <- sapply(common_celltypes, function(ct) {
    est <- estimated_proportions[, ct]
    gt  <- ground_truth_proportions[, ct]
    if (sd(est) < epsilon) est <- est + rnorm(length(est), mean = 0, sd = epsilon)
    if (sd(gt)  < epsilon) gt  <- gt  + rnorm(length(gt),  mean = 0, sd = epsilon)
    cor(est, gt, method = "pearson")
  })
  
  # Spearman correlation per cell type
  spearman_celltype_cor <- sapply(common_celltypes, function(ct) {
    est <- estimated_proportions[, ct]
    gt  <- ground_truth_proportions[, ct]
    if (sd(est) < epsilon) est <- est + rnorm(length(est), mean = 0, sd = epsilon)
    if (sd(gt)  < epsilon) gt  <- gt  + rnorm(length(gt),  mean = 0, sd = epsilon)
    cor(est, gt, method = "spearman")
  })
  

  
  KL_div <- function(P, Q, eps = epsilon) {

    P <- P + eps
    Q <- Q + eps
    sum(P * log(P / Q))
  }
  
  JSD_func <- function(P, Q, eps = epsilon) {

    P <- P / sum(P)
    Q <- Q / sum(Q)
    M <- 0.5 * (P + Q)
    sqrt(0.5 * KL_div(P, M, eps) + 0.5 * KL_div(Q, M, eps))
  }
  
  # Per cell type JSD: distribution = proportions across samples for that cell type
  per_celltype_jsd <- sapply(common_celltypes, function(ct) {
    true_vec <- ground_truth_proportions[, ct]
    est_vec  <- estimated_proportions[, ct]
    JSD_func(true_vec, est_vec, epsilon)
  })
  
  # Per sample JSD: distribution = proportions across cell types for that sample
  per_sample_jsd <- sapply(common_samples, function(sm) {
    true_vec <- ground_truth_proportions[sm, ]
    est_vec  <- estimated_proportions[sm, ]
    JSD_func(true_vec, est_vec, epsilon)
  })
  
  message("\n[SUCCESS]: All metrics calculated!\n")
  
  return(list(
    cell_type_rmse = cell_type_rmse,
    cell_type_mad = cell_type_mad,
    cell_type_nmae = cell_type_nmae,
    mean_ground_truth_for_nmae = mean_ground_truth,
    pearson_celltype_cor = pearson_celltype_cor,
    spearman_celltype_cor = spearman_celltype_cor,
    per_celltype_jsd = per_celltype_jsd,
    per_sample_jsd = per_sample_jsd
  ))
}

#' internal: check if datatype is correct
#'
#' @param file file
#' @param columns object single or list specifying column
#' @param option 1 for character, 2 for numeric
#'
#' @author Sergej Ruff
#' @noRd
q_pm_check_input_type <- function(file, columns, option) {
  if (!option %in% c(1, 2)) {
    stop("Invalid option. Please choose 1 for character or 2 for numeric.")
  }
  expected_class <- if (option == 1) "character" else c("numeric", "integer")
  all_names <- colnames(file)
  for (col in columns) {
    if (!(col %in% all_names)) {
      stop("Column '", col, "' not found in file. Available column names: ", paste(all_names, collapse = ", "))
    }
    if (!inherits(file[, col], expected_class)) {
      stop("Error: Column '", col, "' must be of type ", expected_class, ".")
    }
  }
}