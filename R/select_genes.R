#' Identify top differentially expressed genes per cell type
#'
#' Internal helper that performs one-versus-rest differential expression
#' analysis for each cell type in a `Seurat` or `SingleCellExperiment` object
#' and returns the top marker genes per group.
#'
#' For `SingleCellExperiment` input, the object is first converted to a Seurat
#' object using the count matrix and `colData`. Cell identities are then set
#' from `ident.column`, low-abundance cell types can be removed, optional
#' gene-level filtering is applied, and `Seurat::FindMarkers()` is run for each
#' remaining identity class against all other cells.
#'
#' @param object A `Seurat` object or `SingleCellExperiment` object containing
#'   single-cell expression data.
#' @param ident.column Character scalar giving the metadata column that defines
#'   the grouping variable used for differential expression. If `NULL`, defaults
#'   to `"cluster"` for `SingleCellExperiment` input and `"ident"` otherwise.
#' @param logfc.threshold Numeric scalar passed to `Seurat::FindMarkers()` to
#'   require a minimum log-fold-change threshold for testing.
#' @param p.val.threshold Numeric scalar. Adjusted p-value threshold used to
#'   retain marker genes after differential expression testing.
#' @param top.n Integer scalar or `NULL`. If not `NULL`, only the top `top.n`
#'   genes ranked by decreasing `avg_log2FC` are retained per cell type.
#' @param only.pos Logical scalar. If `TRUE`, only positively enriched marker
#'   genes are returned.
#' @param test.use Character scalar. Statistical test passed to
#'   `Seurat::FindMarkers()`, for example `"wilcox"`.
#' @param assay Character scalar or `NULL`. Optional assay name to use for
#'   differential expression. If `NULL`, the current default assay is used.
#' @param slot Character scalar specifying which assay slot to use. Must be one
#'   of `"counts"`, `"data"`, or `"scale.data"`.
#' @param min.cells.per.group Integer scalar. Minimum number of cells required
#'   for a cell type to be retained for testing.
#' @param min.feature.count Numeric scalar. If greater than zero, genes with
#'   total counts below this threshold are removed before differential
#'   expression analysis.
#'
#' @return A named list of data frames, one per tested cell type. Each element
#'   contains the corresponding `Seurat::FindMarkers()` result table filtered by
#'   adjusted p-value and ordered by decreasing `avg_log2FC`.
#'
#' @details
#' The function performs one-versus-rest differential expression for each cell
#' type separately. Cell types with fewer than `min.cells.per.group` cells are
#' removed before testing. If `slot = "counts"` is requested for a
#' `SingleCellExperiment` input, the object is normalized after conversion and
#' the analysis is performed on the `"data"` slot.
#'
#' When `min.feature.count > 0`, genes are filtered using total counts from the
#' counts matrix before testing. If no significant marker genes remain for any
#' cell type after filtering by `p.val.threshold`, the function stops with an
#' error.
#'
#' @author Sergej Ruff
#' @keywords internal
#' @noRd
getTopDEGs <- function(
    object,
    ident.column = NULL,
    logfc.threshold = 0.25,
    p.val.threshold = 0.05,
    top.n = NULL,
    only.pos = TRUE,
    test.use = "wilcox",
    assay = NULL,
    slot = "data",
    min.cells.per.group = 3,
    min.feature.count = 0
) {




  if (inherits(object, "SingleCellExperiment")) {


    counts <- SingleCellExperiment::counts(object)
    metadata <- as.data.frame(SummarizedExperiment::colData(object))


    if (is.null(ident.column)) ident.column <- "cluster"

    object <- Seurat::CreateSeuratObject(
      counts = counts,
      meta.data = metadata
    )


    if (slot == "counts") {
      object <- Seurat::NormalizeData(object)
      slot <- "data"  # Switch to normalized data
    }
  }


  if (!inherits(object, "Seurat")) {
    stop("Input must be a SingleCellExperiment or Seurat object.")
  }


  if (is.null(ident.column)) ident.column <- "ident"


  if (!ident.column %in% colnames(object@meta.data)) {
    stop("Column '", ident.column, "' not found in metadata.")
  }
  object <- Seurat::SetIdent(object, value = ident.column)




  cell.counts <- table(Seurat::Idents(object))
  # object <- object[rowSums(object) > min.feature.count,]
  keep.celltypes <- names(cell.counts)[cell.counts >= min.cells.per.group]


  if (length(keep.celltypes) == 0) {
    stop("No cell types have at least ", min.cells.per.group, " cells.")
  }

  if (length(keep.celltypes) < length(cell.counts)) {
    removed <- length(cell.counts) - length(keep.celltypes)
    warning("Removed ", removed, " cell type(s) with fewer than ",
            min.cells.per.group, " cells.")
    object <- subset(object, idents = keep.celltypes)
  }


  if (!is.null(assay)) {
    if (!assay %in% names(object@assays)) {
      stop("Assay '", assay, "' not found.")
    }
    Seurat::DefaultAssay(object) <- assay
  }


  assay <- Seurat::DefaultAssay(object)
  if (!slot %in% c("counts", "data", "scale.data")) {
    stop("Slot must be one of: 'counts', 'data', or 'scale.data'")
  }


  if (min.feature.count > 0) {
    if (slot == "counts") {
      expr.matrix <- Seurat::GetAssayData(object, assay = assay, slot = slot)
    } else {

      expr.matrix <- Seurat::GetAssayData(object, assay = assay, slot = "counts")
    }

    keep.genes <- Matrix::rowSums(expr.matrix) >= min.feature.count
    if (sum(keep.genes) == 0) {
      stop("No genes pass the min.feature.count threshold of ", min.feature.count)
    }

    if (sum(!keep.genes) > 0) {
      message("Removing ", sum(!keep.genes), " genes with total counts < ", min.feature.count)
      object <- object[keep.genes, ]
    }
  }


  if (nrow(Seurat::GetAssayData(object, assay = assay, slot = slot)) == 0) {
    if (slot == "data") {
      message("Normalizing data...")
      object <- Seurat::NormalizeData(object)
    } else if (slot == "scale.data") {
      message("Normalizing and scaling data...")
      object <- Seurat::NormalizeData(object)
      object <- Seurat::ScaleData(object)
    }
  }


  cell.types <- Seurat::Idents(object) |> unique() |> as.character()
  de_results <- list()


  for (ct in cell.types) {
    tryCatch({
      markers <- Seurat::FindMarkers(
        object = object,
        ident.1 = ct,
        ident.2 = NULL,
        logfc.threshold = logfc.threshold,
        test.use = test.use,
        only.pos = only.pos,
        min.pct = 0.1,
        slot = slot
      )


      markers <- markers[markers$p_val_adj < p.val.threshold, ]
      if (nrow(markers) > 0) {
        markers <- markers[order(markers$avg_log2FC, decreasing = TRUE), ]
        if (!is.null(top.n)) {
          markers <- head(markers, top.n)
        }
        de_results[[ct]] <- markers
      }
    }, error = function(e) {
      warning("Failed to analyze cell type '", ct, "': ", e$message)
    })
  }

  if (length(de_results) == 0) {
    stop("No differentially expressed genes found in any cell type.")
  }

  return(de_results)
}




#' Select informative genes for QUASAR deconvolution
#'
#' Selects a gene set for QUASAR by combining cell-type-specific
#' differentially expressed genes and, optionally, highly variable genes
#' derived from a cell-type reference matrix.
#'
#' The function accepts a `Seurat` or `SingleCellExperiment` object, filters
#' low-abundance cell types and optionally low-count genes, performs optional
#' normalization, identifies marker genes per cell type, optionally ranks genes
#' by variance across cell-type reference profiles, and returns a final gene set
#' of size up to `total.genes`.
#'
#' @param object A `Seurat` object or `SingleCellExperiment` object containing
#'   single-cell expression data.
#' @param ident.column Character scalar giving the metadata column that defines
#'   the cell-type labels. If `NULL`, defaults to `"cluster"` for
#'   `SingleCellExperiment` input and `"ident"` otherwise.
#' @param total.genes Integer scalar. Maximum total number of genes to retain in
#'   the final output.
#' @param genes.per.cluster Integer scalar. Maximum number of differentially
#'   expressed genes to retain per cell type. If set to `0`, no DE genes are
#'   selected and `variable.genes` must be `TRUE`.
#' @param variable.genes Logical scalar. If `TRUE`, additional genes are ranked
#'   by variance across the constructed reference matrix and added after
#'   differential expression genes until `total.genes` is reached.
#' @param logfc.threshold Numeric scalar. Minimum log-fold-change threshold used
#'   in differential expression testing.
#' @param p.val.threshold Numeric scalar. Adjusted p-value threshold used to
#'   retain differentially expressed genes.
#' @param only.pos Logical scalar. If `TRUE`, only positively enriched marker
#'   genes are retained.
#' @param assay Character scalar. Assay name to use.
#' @param slot Character scalar. Assay slot used for downstream calculations,
#'   typically `"data"` or `"counts"`.
#' @param min.cells.per.group Integer scalar. Minimum number of cells required
#'   for a cell type to be kept in the analysis.
#' @param min.feature.count Numeric scalar. If greater than zero, genes with
#'   total counts below this threshold are removed before gene selection.
#' @param chunk.percentage Numeric scalar in `(0, 100]`. Percentage of cells
#'   processed per chunk when constructing the reference matrix for variable-gene
#'   ranking.
#' @param verbose Logical scalar. If `TRUE`, progress messages and timings are
#'   printed.
#'
#' @return A named list with components:
#' \describe{
#'   \item{genes}{Character vector of selected gene names.}
#'   \item{celltypes}{Character vector of retained cell-type labels.}
#'   \item{n_diff_genes}{Integer scalar giving the number of unique
#'   differentially expressed genes retained before adding variable genes.}
#'   \item{n_var_genes}{Integer scalar giving the number of variable genes added
#'   to complete the final gene set.}
#' }
#'
#' @details
#' The function proceeds in four main stages.
#'
#' First, the input object is standardized to a Seurat object if needed, the
#' identity column is validated, and cell types with fewer than
#' `min.cells.per.group` cells are removed.
#'
#' Second, if `min.feature.count > 0`, genes with low total counts are removed.
#' If `slot = "data"` and the corresponding assay data are empty, the object is
#' normalized before downstream analysis.
#'
#' Third, if `variable.genes = TRUE`, a cell-type reference matrix is built in
#' chunks using `quasar_create_ref()`. Gene-wise variance across the resulting
#' cell-type profiles is then computed and used to rank candidate variable
#' genes.
#'
#' Fourth, if `genes.per.cluster != 0`, differential expression analysis is
#' performed separately for each retained cell type using
#' `Seurat::FindMarkers()`. Significant genes are filtered by
#' `p.val.threshold`, ordered by decreasing `avg_log2FC` and increasing
#' p-value, and up to `genes.per.cluster` genes are retained per cell type.
#'
#' The final selected gene set is constructed as follows:
#' \itemize{
#'   \item If `variable.genes = TRUE`, all retained DE genes are included first,
#'   and the remaining slots up to `total.genes` are filled with the top-ranked
#'   variable genes not already selected.
#'   \item If `variable.genes = FALSE`, only DE genes are used.
#'   \item If `genes.per.cluster = 0`, then `variable.genes` must be `TRUE`.
#' }
#'
#' @examples
#' \dontrun{
#' selected <- quasar_gene_selector(
#'   object = sce,
#'   ident.column = "celltype",
#'   total.genes = 1000,
#'   genes.per.cluster = 30,
#'   variable.genes = TRUE
#' )
#'
#' selected$genes
#' selected$celltypes
#' }
#'
#' @author Sergej Ruff
#' @export
quasar_gene_selector <- function(
    object,
    ident.column = NULL,
    total.genes = 1000,
    genes.per.cluster = 30,
    variable.genes = FALSE,
    logfc.threshold = 0.25,
    p.val.threshold = 0.05,
    only.pos = TRUE,
    assay = "RNA",
    slot = "data",
    min.cells.per.group = 3,
    min.feature.count = 0,
    chunk.percentage = 10,
    verbose = TRUE
) {

  start_time <- Sys.time()
  if (verbose) message("\n=== Starting QUASAR gene selection ===")


  if (chunk.percentage <= 0 || chunk.percentage > 100) {
    stop("chunk.percentage must be between 1 and 100")
  }


  if (inherits(object, "SingleCellExperiment")) {
    if (verbose) message("\n[1] Converting SingleCellExperiment to Seurat...")
    conv_start <- Sys.time()

    counts <- SingleCellExperiment::counts(object)
    metadata <- as.data.frame(SummarizedExperiment::colData(object))

    if (is.null(ident.column)) ident.column <- "cluster"

    object <- Seurat::CreateSeuratObject(
      counts = counts,
      meta.data = metadata
    )

    if (slot == "counts") {
      object <- Seurat::NormalizeData(object)
      slot <- "data"
    }

    rm(counts, metadata)
    gc()

    if (verbose) message("Conversion completed in ", format(Sys.time() - conv_start, digits = 2))
  }

  if (!inherits(object, "Seurat")) {
    stop("Input must be a SingleCellExperiment or Seurat object")
  }


  if (verbose) message("\n[2] Setting identity column...")
  id_start <- Sys.time()

  if (is.null(ident.column)) ident.column <- "ident"
  if (!ident.column %in% colnames(object@meta.data)) {
    stop("Column '", ident.column, "' not found in metadata")
  }

  if (verbose) message("Using identity column: ", ident.column)
  if (verbose) message("Step completed in ", format(Sys.time() - id_start, digits = 2))


  if (min.feature.count > 0) {
    if (verbose) message("\n[3] Filtering low-count genes...")
    filter_start <- Sys.time()

    counts <- Seurat::GetAssayData(object, assay = assay, slot = "counts")
    keep.genes <- Matrix::rowSums(counts) > min.feature.count
    n_removed <- sum(!keep.genes)
    object <- object[keep.genes, ]

    rm(counts, keep.genes)
    gc()

    if (verbose) message("Removed ", n_removed, " genes with counts < ", min.feature.count)
    if (verbose) message("Step completed in ", format(Sys.time() - filter_start, digits = 2))
  }


  if (verbose) message("\n[4] Filtering cell types...")
  filter_start <- Sys.time()

  cell_counts <- table(object@meta.data[[ident.column]])
  valid_celltypes <- names(cell_counts)[cell_counts >= min.cells.per.group]

  if (length(valid_celltypes) == 0) {
    stop("No cell types have at least ", min.cells.per.group, " cells")
  }

  cells_to_keep <- which(object@meta.data[[ident.column]] %in% valid_celltypes)
  object <- subset(object, cells = cells_to_keep)
  object <- Seurat::SetIdent(object, value = ident.column)

  n_removed <- length(cell_counts) - length(valid_celltypes)
  if (n_removed > 0 && verbose) {
    message("Removed ", n_removed, " cell types with < ", min.cells.per.group, " cells")
  }

  rm(cell_counts, cells_to_keep)
  gc()

  if (verbose) message("Keeping ", length(valid_celltypes), " cell types")
  if (verbose) message("Step completed in ", format(Sys.time() - filter_start, digits = 2))


  if (verbose) message("\n[5] Setting assay...")
  assay_start <- Sys.time()

  if (!assay %in% names(object@assays)) {
    stop("Assay '", assay, "' not found in the object")
  }
  Seurat::DefaultAssay(object) <- assay

  if (verbose) message("Using assay: ", assay)
  if (verbose) message("Step completed in ", format(Sys.time() - assay_start, digits = 2))


  if (slot == "data" && nrow(Seurat::GetAssayData(object, slot = "data")) == 0) {
    if (verbose) message("\n[6] Normalizing data...")
    norm_start <- Sys.time()

    object <- Seurat::NormalizeData(object, verbose = FALSE)

    if (verbose) message("Step completed in ", format(Sys.time() - norm_start, digits = 2))
  }


  var_genes_all <- character()
  if (variable.genes) {
    if (verbose) message("\n[7] Running quasar_create_ref for variable genes...")
    var_start <- Sys.time()

    expr_matrix <- Seurat::GetAssayData(object, assay = assay, slot = slot)
    meta_data <- object@meta.data
    n_cells <- ncol(expr_matrix)

    # Calculate dynamic chunk size (at least 1 cell)
    chunk_size <- max(1, floor(n_cells * chunk.percentage / 100))
    n_chunks <- ceiling(n_cells / chunk_size)

    if (verbose) {
      message("Processing ", n_cells, " cells in ", n_chunks, " chunks (~",
              chunk_size, " cells per chunk)")
    }

    ref_matrix <- matrix(0, nrow = nrow(expr_matrix), ncol = length(valid_celltypes))
    rownames(ref_matrix) <- rownames(expr_matrix)
    colnames(ref_matrix) <- valid_celltypes

    celltype_counts <- table(meta_data[[ident.column]])

    for (i in seq_len(n_chunks)) {
      if (verbose && n_chunks > 1) {
        message("Processing chunk ", i, " of ", n_chunks)
      }

      start_idx <- (i-1)*chunk_size + 1
      end_idx <- min(i*chunk_size, n_cells)

      # Use sparse matrices to avoid dense coercion warnings
      chunk_expr <- Seurat::GetAssayData(
        object,
        assay = assay,
        slot = slot,
        layer = "data"
      )[, start_idx:end_idx, drop = FALSE]

      chunk_meta <- meta_data[start_idx:end_idx, , drop = FALSE]

      ref_chunk <- quasar_create_ref(
        as.matrix(chunk_expr),
        as.matrix(chunk_meta),
        ct_col = ident.column,
        verbose = FALSE
      )$ref_matrix


      for (ct in colnames(ref_chunk)) {
        ref_matrix[, ct] <- ref_matrix[, ct] + ref_chunk[, ct] * celltype_counts[ct]
      }

      rm(chunk_expr, ref_chunk)
      gc()
    }


    for (ct in valid_celltypes) {
      ref_matrix[, ct] <- ref_matrix[, ct] / celltype_counts[ct]
    }


    gene_vars <- matrixStats::rowVars(ref_matrix)
    names(gene_vars) <- rownames(ref_matrix)
    var_genes_all <- names(sort(gene_vars, decreasing = TRUE))

    rm(expr_matrix, meta_data, ref_matrix, gene_vars)
    gc()

    if (verbose) message("Found ", length(var_genes_all), " variable genes")
    if (verbose) message("Step completed in ", format(Sys.time() - var_start, digits = 2))
  }


  de_genes <- character()
  de_results <- list()

  if (genes.per.cluster != 0) {
    if (verbose) message("\n[8] Performing DE analysis...")
    de_start <- Sys.time()

    for (celltype in valid_celltypes) {
      if (verbose) message(paste("Running DE for cell type:",celltype))
      tryCatch({
        markers <- Seurat::FindMarkers(
          object = object,
          ident.1 = celltype,
          logfc.threshold = logfc.threshold,
          min.pct = min.feature.count,
          test.use = "wilcox",
          only.pos = only.pos,
          min.cells.group = min.cells.per.group,
          assay = assay,
          slot = slot,
          verbose = FALSE
        )

        if (!is.null(markers) && nrow(markers) > 0) {
          markers <- markers[markers$p_val_adj <= p.val.threshold, ]
          markers <- markers[order(-markers$avg_log2FC, markers$p_val), ]
          markers <- markers[1:min(genes.per.cluster, nrow(markers)), ]
          de_results[[celltype]] <- markers
          de_genes <- unique(c(de_genes, rownames(markers)))
        }
      }, error = function(e) {
        if (verbose) message("Skipping ", celltype, ": ", e$message)
      })
    }

    n_diff_genes <- length(de_genes)

    if (verbose) {
      message("Found ", n_diff_genes, " DE genes across ", length(de_results), " cell types")
      message("Step completed in ", format(Sys.time() - de_start, digits = 2))
    }
  } else {
    n_diff_genes <- 0
  }


  if (verbose) message("\n[9] Selecting final gene set...")
  select_start <- Sys.time()

  if (variable.genes) {
    var_genes_all <- setdiff(var_genes_all, de_genes)
    n_needed <- max(0, total.genes - n_diff_genes)
    var_genes <- var_genes_all[1:min(n_needed, length(var_genes_all))]
    n_var_genes <- length(var_genes)
    selected_genes <- unique(c(de_genes, var_genes))[1:min(total.genes, n_diff_genes + n_var_genes)]
  } else if (genes.per.cluster == 0) {
    stop("When genes.per.cluster = 0, variable.genes must be TRUE")
  } else {
    selected_genes <- de_genes[1:min(total.genes, length(de_genes))]
    n_var_genes <- 0
  }

  if (verbose) {
    message("Selected ", length(selected_genes), " genes (", n_diff_genes, " DE + ", n_var_genes, " variable)")
    message("Step completed in ", format(Sys.time() - select_start, digits = 2))
  }


  rm(var_genes_all, de_genes)
  gc()


  if (verbose) {
    message("\n=== QUASAR gene selection completed ===")
    message("Total time: ", format(Sys.time() - start_time, digits = 2))
    message("Final counts:")
    message("- Cell types: ", length(valid_celltypes))
    message("- DE genes: ", n_diff_genes)
    message("- Variable genes: ", n_var_genes)
    message("- Total selected genes: ", length(selected_genes))
  }

  return(list(
    genes = selected_genes,
    celltypes = valid_celltypes,
    n_diff_genes = n_diff_genes,
    n_var_genes = n_var_genes
  ))
}
