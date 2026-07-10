

# ----------------------------- small utilities ------------------------------

quasar_msg <- function(..., verbose = TRUE) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  cat("[Quasar Importer/Exporter] ", paste0(...), "\n", sep = "")
  invisible(NULL)
}

# one-line header + an aligned key/value block, e.g.
#   [quasar] import done in 17.6s
#            source     adata.raw.X
#            shape      71439 cells x 26727 features
quasar_summary <- function(header, items, seconds, verbose = TRUE) {
  if (!isTRUE(verbose)) return(invisible(NULL))
  cat(sprintf("[Quasar Importer/Exporter] %s in %.1fs\n", header, seconds))
  w <- max(nchar(names(items)))
  for (k in names(items))
    cat(sprintf("         %-*s  %s\n", w, k, items[[k]]))
  invisible(NULL)
}

quasar_require <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Package '", pkg, "' is required for this operation but is not installed.",
         call. = FALSE)
  invisible(TRUE)
}

# Seurat rejects underscores in feature names and dislikes empty names.
quasar_seurat_feature_names <- function(x) {
  x <- as.character(x)
  x <- gsub("_", "-", x)
  empty <- is.na(x) | !nzchar(x)
  if (any(empty)) x[empty] <- paste0("feature-", which(empty))
  make.unique(x)
}


# ============================================================================
#  HDF5 READERS  
# ============================================================================

# Decode one column coming out of rhdf5::h5read(): plain array, a categorical
# group {categories, codes}, or a nullable/masked group {mask, values}.
quasar_parse_column <- function(x) {
  if (is.array(x)) return(as.vector(x))
  if (is.list(x)) {
    n <- names(x)
    if (!is.null(n) && all(n %in% c("categories", "codes"))) {
      codes <- as.vector(x$codes) + 1L
      codes[codes == 0] <- NA                       # anndata -1 sentinel -> NA
      lev <- quasar_parse_column(x$categories)
      return(lev[codes])
    }
    if (!is.null(n) && all(n %in% c("mask", "values"))) {
      v   <- as.vector(x$values)
      bad <- x$mask %in% c(TRUE, 1L, "TRUE", "1")   # TRUE == missing
      v[bad] <- NA
      return(v)
    }
  }
  if (is.atomic(x)) as.vector(x) else x
}

# rhdf5 turns column names containing "/" into nested lists; flatten those back
# out, while leaving genuine categorical / masked groups intact.
quasar_flatten_columns <- function(collist) {
  repeat {
    nested <- logical(length(collist))
    for (i in seq_along(collist)) {
      nm <- names(collist)[i]
      if (identical(nm, "__categories")) next
      el <- collist[[i]]
      if (is.list(el)) {
        en <- names(el)
        is_cat  <- !is.null(en) && all(en %in% c("categories", "codes"))
        is_mask <- !is.null(en) && all(en %in% c("mask", "values"))
        if (!is_cat && !is_mask) nested[i] <- TRUE
      }
    }
    idx <- which(nested)
    if (!length(idx)) break
    i   <- idx[1]
    sub <- collist[[i]]
    for (j in seq_along(sub))
      collist[[paste0(names(collist)[i], "/", names(sub)[j])]] <- sub[[j]]
    collist[[i]] <- NULL
  }
  collist
}

# Legacy: factor levels kept in /uns/<col>_categories, column holds integer codes.
quasar_apply_uns_categories <- function(filename, cols) {
  ls <- rhdf5::h5ls(filename)
  unsc <- ls[ls$group == "/uns" & endsWith(ls$name, "_categories"), , drop = FALSE]
  if (!nrow(unsc)) return(cols)
  unsc$var <- sub("_categories$", "", unsc$name)
  unsc <- unsc[unsc$var %in% names(cols), , drop = FALSE]
  for (k in seq_len(nrow(unsc))) {
    fn <- unsc$var[k]
    if (!is.numeric(cols[[fn]])) next
    lev   <- as.vector(rhdf5::h5read(filename, paste0("/uns/", unsc$name[k])))
    codes <- cols[[fn]] + 1L
    codes[codes == 0] <- NA
    cols[[fn]] <- lev[codes]
  }
  cols
}

#' Read an anndata dataframe group (obs / var) into a data.frame
#'
#' Handles the modern dataframe encoding (per-column categorical groups,
#' nullable/masked integer & boolean groups, string and numeric arrays) as well
#' as two legacy factor layouts (a shared \code{__categories} group and
#' \code{/uns/<col>_categories}). The index dataset becomes the rownames.
#'
#' @param filename path to the .h5ad file.
#' @param name group name, e.g. \code{"obs"}, \code{"var"} or \code{"raw/var"}.
#' @param index_as_column keep the index as a column in addition to using it as
#'   rownames. Default \code{FALSE} (cleaner colData / meta.data).
#'
#' @return a data.frame with rownames taken from the anndata index.
#' @export
quasar_h5_read_dataframe <- function(filename, name, index_as_column = FALSE) {
  raw      <- rhdf5::h5read(filename, name, read.attributes = TRUE)
  grp_attr <- attributes(raw)
  raw      <- quasar_flatten_columns(raw)

  # decode every real column
  cols <- list()
  for (nm in names(raw)) {
    if (identical(nm, "__categories")) next
    cols[[nm]] <- quasar_parse_column(raw[[nm]])
  }

  # legacy shared "__categories" group
  if (!is.null(raw[["__categories"]])) {
    for (fn in names(raw[["__categories"]])) {
      lev   <- as.vector(raw[["__categories"]][[fn]])
      codes <- cols[[fn]] + 1L
      codes[codes == 0] <- NA
      cols[[fn]] <- lev[codes]
    }
  }
  cols <- quasar_apply_uns_categories(filename, cols)

  # assemble, dropping anything with a mismatched length
  len  <- vapply(cols, length, integer(1))
  keep <- len == max(len)
  if (!all(keep))
    warning("dropping columns of unexpected length in '", name, "'", call. = FALSE)
  df <- as.data.frame(cols[keep], check.names = FALSE, stringsAsFactors = FALSE)

  # resolve the index dataset
  idx_name <- grp_attr[["_index"]]
  if (is.null(idx_name) || !(idx_name %in% colnames(df))) {
    idx_name <- if ("_index" %in% colnames(df)) "_index"
                else if ("index" %in% colnames(df)) "index"
                else NULL
  }
  if (!is.null(idx_name))
    rownames(df) <- make.unique(as.character(df[[idx_name]]))

  # column ordering per the anndata attribute
  ord <- grp_attr[["column-order"]]
  if (!is.null(ord)) {
    ord <- ord[ord %in% colnames(df)]
    if (index_as_column && !is.null(idx_name) && !(idx_name %in% ord))
      ord <- c(idx_name, ord)
    if (length(ord)) df <- df[, ord, drop = FALSE]
  }

  # drop the index column unless asked to keep it (preserving nrow / rownames)
  if (!index_as_column && !is.null(idx_name)) {
    rn <- rownames(df)
    df[[idx_name]] <- NULL
    if (ncol(df) == 0) df <- data.frame(row.names = rn)
  }
  df
}

#' Read an anndata matrix group / dataset (X, raw/X, layers, obsm) as an R matrix
#'
#' Sparse CSR/CSC groups become a \pkg{Matrix} sparse matrix; dense datasets a
#' base matrix. Orientation matches the R single-cell convention
#' (features x observations for X). A dataframe stored where a matrix is expected
#' is read and coerced to a numeric matrix.
#'
#' @param filename path to the .h5ad file.
#' @param name group / dataset name, e.g. \code{"X"}, \code{"raw/X"},
#'   \code{"obsm/X_pca"}.
#'
#' @return a dense base matrix or a sparse \code{dgCMatrix}.
#' @export
#' @import Matrix
quasar_h5_read_matrix <- function(filename, name) {
  if (!startsWith(name, "/")) name <- paste0("/", name)
  attr <- rhdf5::h5readAttributes(filename, name)
  enc  <- attr[["encoding-type"]]
  if (is.null(enc)) enc <- attr$h5sparse_format

  # dataframe-as-matrix
  if (!is.null(enc) && enc == "dataframe") {
    df <- quasar_h5_read_dataframe(filename, name, index_as_column = FALSE)
    return(as.matrix(df))
  }

  is_sparse <- !is.null(enc) && (startsWith(enc, "csr") || startsWith(enc, "csc"))
  if (is_sparse) {
    ls  <- rhdf5::h5ls(filename)
    nnz <- suppressWarnings(as.numeric(ls$dim[ls$group == name & ls$name == "data"]))
    if (length(nnz) == 1 && !is.na(nnz) && nnz >= 2^31)
      stop("Expression matrix has >= 2^31-1 non-zero values; too large for a ",
           "standard R/Matrix sparse matrix. Use Python, or load in chunks.",
           call. = FALSE)
  }

  obj <- rhdf5::h5read(filename, name)

  # dense dataset: rhdf5 returns it already transposed to features x observations
  if (!is.list(obj)) return(obj)
  if (is.null(obj$data)) return(as.matrix(as.data.frame(obj)))

  shape <- attr$shape
  if (is.null(shape)) shape <- attr$h5sparse_shape         # (n_obs, n_vars)
  data  <- as.numeric(obj$data)
  i1    <- obj$indices + 1L

  if (!is.null(enc) && startsWith(enc, "csr")) {
    # CSR: indptr per-observation, indices are feature indices
    Matrix::sparseMatrix(i = i1, p = obj$indptr, x = data,
                         dims = rev(as.integer(shape)))     # features x obs
  } else {
    m <- Matrix::sparseMatrix(i = i1, p = obj$indptr, x = data,
                              dims = as.integer(shape))      # obs x var
    Matrix::t(m)                                             # features x obs
  }
}


# =============================================================================
#  IMPORT :  h5ad  ->  Seurat / SingleCellExperiment
# =============================================================================

#' Import an .h5ad file as a Seurat or SingleCellExperiment object
#'
#' A single entry point that reads adata.X (or adata.raw.X), adata.obs,
#' adata.var and adata.obsm and assembles either a \code{Seurat} or a
#' \code{SingleCellExperiment}. A missing \code{adata.raw} downgrades gracefully
#' to \code{adata.X} with a warning rather than failing.
#'
#' @param filename path to the .h5ad file.
#' @param as one of \code{"seurat"} or \code{"sce"}; the object type to return.
#' @param use.raw logical, default \code{TRUE}. Use \code{adata.raw.X} /
#'   \code{adata.raw.var} when present, otherwise \code{adata.X} / \code{adata.var}.
#' @param load.X logical, whether to load the expression matrix. If \code{FALSE}
#'   an all-zero sparse matrix of the right shape is used (much faster).
#' @param load.obsm logical, whether to load \code{adata.obsm} as reduced
#'   dimensions / DimReducs.
#' @param assay Seurat assay name to populate (ignored for \code{as = "sce"}).
#' @param verbose logical, print a start line and a summary block.
#'
#' @return a \code{Seurat} or \code{SingleCellExperiment} object.
#' @export
#' @import Matrix
#' @examples
#' \dontrun{
#' ## Round-trip on tiny data: write a small bulk matrix, then read it back.
#' f   <- tempfile(fileext = ".h5ad")
#' mat <- matrix(rpois(15 * 8, 2), nrow = 15,
#'               dimnames = list(paste0("gene", 1:15), paste0("cell", 1:8)))
#' quasar_h5adexporter(mat, f, bulk_orientation = "genes_x_samples")
#'
#' sce <- quasar_h5adimport(f, as = "sce",    use.raw = FALSE)
#' seu <- quasar_h5adimport(f, as = "seurat", use.raw = FALSE)
#' }
quasar_h5adimport <- function(filename,
                              as        = c("seurat", "sce"),
                              use.raw   = TRUE,
                              load.X    = TRUE,
                              load.obsm = TRUE,
                              assay     = "RNA",
                              verbose   = TRUE) {
  as <- match.arg(as)
  if (!file.exists(filename)) stop("File not found: ", filename)
  quasar_require("rhdf5"); quasar_require("Matrix")
  t0 <- Sys.time()
  quasar_msg("import started <- ", basename(filename), verbose = verbose)

  h5 <- rhdf5::h5ls(filename)
  has_raw <- any(h5$group == "/raw" & h5$name == "X")
  if (use.raw && !has_raw) {
    warning("use.raw=TRUE but no adata.raw in '", filename,
            "'; falling back to adata.X", call. = FALSE)
    use.raw <- FALSE
  }
  expr_path <- if (use.raw) "raw/X"   else "X"
  var_path  <- if (use.raw) "raw/var" else "var"

  # --- metadata ---
  obs <- quasar_h5_read_dataframe(filename, "obs")
  var <- quasar_h5_read_dataframe(filename, var_path)
  if (as == "seurat") rownames(var) <- quasar_seurat_feature_names(rownames(var))

  # --- expression (features x observations) ---
  if (load.X) {
    X <- quasar_h5_read_matrix(filename, expr_path)
  } else {
    X <- Matrix::sparseMatrix(i = integer(0), p = 0L, x = numeric(0),
                              dims = c(nrow(var), nrow(obs)))
  }
  rownames(X) <- rownames(var)
  colnames(X) <- rownames(obs)

  # --- obsm ---
  obsm <- list()
  if (load.obsm && any(h5$group == "/obsm")) {
    for (n in h5$name[h5$group == "/obsm"]) {
      m <- quasar_h5_read_matrix(filename, paste0("obsm/", n))
      if (is.array(m) || inherits(m, "Matrix")) {
        if (nrow(m) != nrow(obs)) m <- t(m)          # anndata stores obs x k
        rownames(m) <- rownames(obs)
        obsm[[n]] <- m
      }
    }
  }

  out <- if (as == "sce")
    quasar_build_sce(X, obs, var, obsm)
  else
    quasar_build_seurat(X, obs, var, obsm, assay)

  result <- if (as == "sce")
    sprintf("SingleCellExperiment, %d features x %d cells (+%d reducedDims)",
            nrow(out), ncol(out), length(obsm))
  else
    sprintf("Seurat, %d features x %d cells (+%d reductions)",
            nrow(out), ncol(out), length(obsm))

  quasar_summary("import done", c(
    source     = paste0("adata.", if (use.raw) "raw.X" else "X"),
    shape      = sprintf("%d cells x %d features", nrow(obs), nrow(var)),
    expression = if (load.X) class(X)[1] else "not loaded",
    obsm       = if (load.obsm) sprintf("%d matrices", length(obsm)) else "skipped",
    result     = result
  ), as.numeric(Sys.time() - t0, units = "secs"), verbose = verbose)

  out
}

quasar_build_sce <- function(X, obs, var, obsm) {
  quasar_require("SingleCellExperiment")
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays  = list(X = X),
    colData = S4Vectors::DataFrame(obs, check.names = FALSE),
    rowData = S4Vectors::DataFrame(var, check.names = FALSE)
  )
  for (n in names(obsm))
    SingleCellExperiment::reducedDim(sce, n) <- as.matrix(obsm[[n]])
  sce
}

quasar_build_seurat <- function(X, obs, var, obsm, assay) {
  quasar_require("Seurat"); quasar_require("SeuratObject")

  # A dense, near-full matrix is almost certainly scaled data; feeding it in as
  # 'counts' would make Seurat sparsify it (and possibly overflow). Route it to
  # the data slot instead.
  scaled <- FALSE
  if (is.matrix(X)) {
    nz   <- sum(X != 0, na.rm = TRUE)
    dens <- nz / length(X)
    scaled <- dens > 0.8 || nz > (2^31 - 1)
  }
  aobj <- if (scaled) SeuratObject::CreateAssayObject(data   = X)
          else        SeuratObject::CreateAssayObject(counts = X)

  seu <- SeuratObject::CreateSeuratObject(counts = aobj, assay = assay)
  seu <- SeuratObject::AddMetaData(seu, obs)
  if (ncol(var))
    seu[[assay]] <- SeuratObject::AddMetaData(seu[[assay]], metadata = var)

  # reductions (rename 'spatial' to stay clear of Seurat's own handling)
  nm <- names(obsm); nm[nm == "spatial"] <- "X_spatial"; names(obsm) <- nm
  for (n in names(obsm)) {
    emb <- as.matrix(obsm[[n]])
    if (ncol(emb) == 0) next
    key <- paste0(gsub("_", "", n), "_")               # Seurat key rules
    colnames(emb) <- paste0(key, seq_len(ncol(emb)))
    rownames(emb) <- colnames(seu)
    seu[[key]] <- SeuratObject::CreateDimReducObject(
      embeddings = emb, key = key, assay = assay)
  }
  seu
}


# =============================================================================
#  EXPORT :  Seurat / SCE / bulk matrix  ->  h5ad
# =============================================================================

#' Export a Seurat, SingleCellExperiment or bulk matrix to .h5ad
#'
#' Single-cell objects contribute a features x cells assay; bulk data is a plain
#' samples x genes matrix (\code{n_obs x n_vars}, the anndata convention). An
#' optional cell-type \code{fractions} matrix (samples x cell types) can be
#' written alongside — useful for deconvolution ground truth.
#'
#' @param object a \code{Seurat}, \code{SingleCellExperiment}, or a
#'   matrix / \code{Matrix} / \code{data.frame} (treated as bulk).
#' @param filename output path (overwritten unless \code{overwrite=FALSE}).
#' @param assay assay to export. Seurat: assay name (default = active assay).
#'   SCE: assay name (default \code{"counts"} if present, else the first).
#'   Ignored for bulk.
#' @param layer Seurat layer/slot to pull the matrix from (default \code{"counts"}).
#' @param obs,var optional per-observation / per-variable metadata
#'   (\code{data.frame} or matrix). Mainly for bulk, where the object carries no
#'   metadata; row counts must match the matrix.
#' @param fractions optional cell-type fraction matrix, samples x cell types
#'   (rows = observations). Column names become cell-type labels.
#' @param fractions_to where to place the fractions: \code{"obsm"} (as
#'   \code{adata.obsm['fractions']}), \code{"obs"} (spread across obs columns),
#'   or \code{"both"}.
#' @param obsm optional named list of extra per-observation matrices (bulk only;
#'   for single-cell use \code{export_reductions}).
#' @param bulk_orientation orientation of a bulk matrix:
#'   \code{"samples_x_genes"} (default) or \code{"genes_x_samples"}.
#' @param export_reductions logical, export Seurat reductions / SCE reducedDims
#'   into \code{adata.obsm}.
#' @param overwrite logical, overwrite an existing file.
#' @param verbose logical, print a start line and a summary block.
#'
#' @return (invisibly) the output \code{filename}.
#' @export
#' @import Matrix
#' @examples
#' \dontrun{
#' ## --- tiny bulk export (6 samples x 20 genes) + ground-truth fractions ----
#' set.seed(1)
#' bulk <- matrix(rpois(6 * 20, 5), nrow = 6, ncol = 20,
#'                dimnames = list(paste0("sample_", 1:6),
#'                                paste0("ENSG", 1:20)))
#' frac <- matrix(runif(6 * 3), nrow = 6,
#'                dimnames = list(rownames(bulk),
#'                                c("Tcell", "Bcell", "Mono")))
#' frac <- frac / rowSums(frac)                         # rows sum to 1
#' quasar_h5adexporter(bulk, tempfile(fileext = ".h5ad"),
#'                     fractions = frac, fractions_to = "both")
#'
#' ## --- tiny SingleCellExperiment export (15 genes x 8 cells) ---------------
#' sce <- SingleCellExperiment::SingleCellExperiment(
#'   assays = list(counts = matrix(rpois(15 * 8, 2), nrow = 15,
#'                 dimnames = list(paste0("gene", 1:15),
#'                                 paste0("cell", 1:8)))))
#' quasar_h5adexporter(sce, tempfile(fileext = ".h5ad"))
#' }
quasar_h5adexporter <- function(object, filename,
                                assay             = NULL,
                                layer             = "counts",
                                obs               = NULL,
                                var               = NULL,
                                fractions         = NULL,
                                fractions_to      = c("obsm", "obs", "both"),
                                obsm              = NULL,
                                bulk_orientation  = c("samples_x_genes",
                                                      "genes_x_samples"),
                                export_reductions = TRUE,
                                overwrite         = TRUE,
                                verbose           = TRUE) {
  quasar_require("rhdf5"); quasar_require("Matrix")
  fractions_to     <- match.arg(fractions_to)
  bulk_orientation <- match.arg(bulk_orientation)
  t0 <- Sys.time()
  quasar_msg("export started -> ", basename(filename), verbose = verbose)

  b <- quasar_extract(object, assay = assay, layer = layer,
                      obs = obs, var = var,
                      bulk_orientation = bulk_orientation,
                      export_reductions = export_reductions, obsm = obsm)

  # --- optional fractions ---
  if (!is.null(fractions)) {
    fr <- as.matrix(fractions)
    if (nrow(fr) != b$n_obs) {
      if (ncol(fr) == b$n_obs) fr <- t(fr)
      else stop("`fractions` has ", nrow(fr), " rows but there are ",
                b$n_obs, " observations.")
    }
    if (is.null(colnames(fr)))
      colnames(fr) <- paste0("frac_", seq_len(ncol(fr)))
    if (!is.null(rownames(fr)) && all(rownames(b$obs) %in% rownames(fr)))
      fr <- fr[rownames(b$obs), , drop = FALSE]
    rownames(fr) <- rownames(b$obs)

    if (fractions_to %in% c("obsm", "both")) b$obsm[["fractions"]] <- fr
    if (fractions_to %in% c("obs", "both"))
      b$obs <- cbind(b$obs, as.data.frame(fr, check.names = FALSE))
  }

  quasar_write_h5ad(b, filename, overwrite = overwrite)

  xdesc <- if (inherits(b$X_vxo, "sparseMatrix"))
    sprintf("CSR sparse [%d x %d], nnz=%d",
            b$n_obs, b$n_vars, Matrix::nnzero(b$X_vxo))
  else
    sprintf("dense [%d x %d]", b$n_obs, b$n_vars)

  quasar_summary("export done", c(
    source = b$source,
    X      = xdesc,
    obs    = sprintf("%d cols", ncol(b$obs)),
    var    = sprintf("%d cols", ncol(b$var)),
    obsm   = if (length(b$obsm)) paste(names(b$obsm), collapse = ", ") else "none",
    output = filename
  ), as.numeric(Sys.time() - t0, units = "secs"), verbose = verbose)

  invisible(filename)
}


# --------------------------- extraction dispatch -----------------------------


quasar_extract <- function(object, assay, layer, obs, var,
                           bulk_orientation, export_reductions, obsm) {
  if (inherits(object, "Seurat"))
    return(quasar_extract_seurat(object, assay, layer, export_reductions))
  if (inherits(object, "SingleCellExperiment"))
    return(quasar_extract_sce(object, assay, export_reductions))
  if (is.matrix(object) || inherits(object, "Matrix") || is.data.frame(object))
    return(quasar_extract_bulk(object, obs, var, bulk_orientation, obsm))
  stop("Unsupported object of class {", paste(class(object), collapse = ", "),
       "}; expected Seurat, SingleCellExperiment or a matrix/data.frame (bulk).")
}

quasar_assay_feature_meta <- function(a) {
  sl <- methods::slotNames(a)
  mf <- if ("meta.features" %in% sl) methods::slot(a, "meta.features")
        else if ("meta.data" %in% sl) methods::slot(a, "meta.data")
        else NULL
  if (is.null(mf) || !NROW(mf)) return(NULL)
  as.data.frame(mf, check.names = FALSE)
}

quasar_seurat_layer <- function(object, assay, layer) {
  m <- tryCatch(SeuratObject::LayerData(object, assay = assay, layer = layer),
                error = function(e) NULL)
  if (is.null(m) || !length(m))
    m <- tryCatch(SeuratObject::GetAssayData(object, assay = assay, slot = layer),
                  error = function(e) NULL)
  if (is.null(m) || !length(m))
    stop("Could not fetch layer/slot '", layer, "' from assay '", assay, "'.")
  m
}

quasar_extract_seurat <- function(object, assay, layer, export_reductions) {
  quasar_require("SeuratObject")
  if (is.null(assay)) assay <- SeuratObject::DefaultAssay(object)

  mat <- quasar_seurat_layer(object, assay, layer)         # features x cells
  obs <- as.data.frame(object@meta.data, check.names = FALSE)
  var <- quasar_assay_feature_meta(object[[assay]])
  if (is.null(var) || nrow(var) != nrow(mat))
    var <- data.frame(row.names = rownames(mat))
  rownames(obs) <- colnames(mat)
  rownames(var) <- rownames(mat)

  obsm <- list()
  if (export_reductions) {
    for (rn in SeuratObject::Reductions(object)) {
      emb <- SeuratObject::Embeddings(object[[rn]])
      if (!is.null(emb) && nrow(emb) == ncol(mat)) obsm[[rn]] <- emb
    }
  }
  list(X_vxo = mat, obs = obs, var = var, obsm = obsm,
       n_obs = ncol(mat), n_vars = nrow(mat), source = "seurat")
}

quasar_extract_sce <- function(object, assay, export_reductions) {
  quasar_require("SingleCellExperiment"); quasar_require("SummarizedExperiment")
  an <- SummarizedExperiment::assayNames(object)
  if (is.null(assay)) assay <- if ("counts" %in% an) "counts" else an[1]

  mat <- SummarizedExperiment::assay(object, assay)        # features x cells
  obs <- as.data.frame(SummarizedExperiment::colData(object), check.names = FALSE)
  var <- as.data.frame(SummarizedExperiment::rowData(object), check.names = FALSE)
  if (is.null(rownames(mat))) rownames(mat) <- paste0("gene_", seq_len(nrow(mat)))
  if (is.null(colnames(mat))) colnames(mat) <- paste0("cell_", seq_len(ncol(mat)))
  rownames(obs) <- colnames(mat)
  rownames(var) <- rownames(mat)

  obsm <- list()
  if (export_reductions) {
    for (rn in SingleCellExperiment::reducedDimNames(object)) {
      emb <- as.matrix(SingleCellExperiment::reducedDim(object, rn))
      if (nrow(emb) == ncol(mat)) obsm[[rn]] <- emb
    }
  }
  list(X_vxo = mat, obs = obs, var = var, obsm = obsm,
       n_obs = ncol(mat), n_vars = nrow(mat), source = "sce")
}

quasar_coerce_meta <- function(meta, n, ids, what) {
  if (is.null(meta)) {
    df <- data.frame(matrix(nrow = n, ncol = 0))
    rownames(df) <- ids
    return(df)
  }
  df <- as.data.frame(meta, check.names = FALSE)
  if (nrow(df) != n)
    stop("`", what, "` has ", nrow(df), " rows but expected ", n, ".")
  rownames(df) <- ids
  df
}

quasar_extract_bulk <- function(object, obs, var, bulk_orientation, obsm) {
  X <- if (is.data.frame(object)) as.matrix(object) else object

  # Canonical internal orientation = features (genes) x observations (samples).
  if (bulk_orientation == "samples_x_genes") {
    Xv    <- t(X)                                   # -> genes x samples
    samp  <- rownames(X); genes <- colnames(X)
  } else {
    Xv    <- X                                      # already genes x samples
    genes <- rownames(X); samp <- colnames(X)
  }
  n_vars <- nrow(Xv); n_obs <- ncol(Xv)
  if (is.null(samp))  samp  <- paste0("sample_", seq_len(n_obs))
  if (is.null(genes)) genes <- paste0("gene_",   seq_len(n_vars))
  rownames(Xv) <- genes; colnames(Xv) <- samp

  obs_df <- quasar_coerce_meta(obs, n_obs, samp, "obs")
  var_df <- quasar_coerce_meta(var, n_vars, genes, "var")
  if (ncol(var_df) == 0)                             # a conventional gene column
    var_df <- data.frame(gene = genes, row.names = genes, check.names = FALSE)

  om <- list()
  if (!is.null(obsm)) for (nn in names(obsm)) {
    mm <- as.matrix(obsm[[nn]])
    if (nrow(mm) == n_obs) { rownames(mm) <- samp; om[[nn]] <- mm }
  }
  list(X_vxo = Xv, obs = obs_df, var = var_df, obsm = om,
       n_obs = n_obs, n_vars = n_vars, source = "bulk")
}


# ------------------------------ HDF5 writers ---------------------------------


quasar_h5_attr <- function(fid, path, attrs) {
  oid <- rhdf5::H5Oopen(fid, path)
  on.exit(rhdf5::H5Oclose(oid), add = TRUE)
  array_attrs <- c("column-order", "shape")
  for (nm in names(attrs)) {
    val <- attrs[[nm]]
    if (length(val) == 0) next                       # skip empty (e.g. no cols)
    scalar <- !(nm %in% array_attrs) && length(val) == 1
    ok <- FALSE
    if (scalar)
      ok <- tryCatch({
        rhdf5::h5writeAttribute(val, h5obj = oid, name = nm, asScalar = TRUE)
        TRUE
      }, error = function(e) FALSE)
    if (!ok)
      rhdf5::h5writeAttribute(val, h5obj = oid, name = nm)
  }
}

quasar_h5_write_string_array <- function(fid, path, x) {
  x <- as.character(x)
  x[is.na(x)] <- ""                                  # fixed-length can't hold NA
  size <- max(1L, max(nchar(x, type = "bytes")) + 1L)
  made <- tryCatch({
    rhdf5::h5createDataset(fid, path, dims = length(x),
                           storage.mode = "character", size = size,
                           encoding = "UTF-8")
    TRUE
  }, error = function(e) FALSE)
  if (!made)
    rhdf5::h5createDataset(fid, path, dims = length(x),
                           storage.mode = "character", size = size)
  rhdf5::h5write(x, fid, path)
  quasar_h5_attr(fid, path, list(`encoding-type` = "string-array",
                                 `encoding-version` = "0.2.0"))
}

quasar_h5_write_array <- function(fid, path, x) {
  if (is.integer(x) && anyNA(x)) x <- as.double(x)   # ints have no NA on disk
  if (is.logical(x))             x <- as.integer(x)
  rhdf5::h5write(x, fid, path)
  quasar_h5_attr(fid, path, list(`encoding-type` = "array",
                                 `encoding-version` = "0.2.0"))
}

quasar_h5_write_categorical <- function(fid, path, x) {
  cats <- levels(x)
  if (length(cats) == 0) {                           # degenerate -> plain strings
    quasar_h5_write_string_array(fid, path, as.character(x)); return(invisible())
  }
  rhdf5::h5createGroup(fid, path)
  codes <- as.integer(x) - 1L
  codes[is.na(codes)] <- -1L                         # anndata NA sentinel
  quasar_h5_write_string_array(fid, paste0(path, "/categories"), cats)
  rhdf5::h5write(as.integer(codes), fid, paste0(path, "/codes"))
  quasar_h5_attr(fid, paste0(path, "/codes"),
                 list(`encoding-type` = "array", `encoding-version` = "0.2.0"))
  quasar_h5_attr(fid, path, list(`encoding-type` = "categorical",
                                 `encoding-version` = "0.2.0",
                                 ordered = as.integer(is.ordered(x))))
}

quasar_h5_write_column <- function(fid, path, x) {
  if (is.factor(x))          quasar_h5_write_categorical(fid, path, x)
  else if (is.character(x))  quasar_h5_write_string_array(fid, path, x)
  else if (is.numeric(x) || is.logical(x)) quasar_h5_write_array(fid, path, x)
  else {
    warning("column '", basename(path), "' has unsupported type '",
            class(x)[1], "'; coercing to character.", call. = FALSE)
    quasar_h5_write_string_array(fid, path, as.character(x))
  }
}

quasar_h5_write_dataframe <- function(fid, path, df, index) {
  rhdf5::h5createGroup(fid, path)
  cols <- colnames(df)
  quasar_h5_write_string_array(fid, paste0(path, "/_index"), as.character(index))
  for (cn in cols) quasar_h5_write_column(fid, paste0(path, "/", cn), df[[cn]])
  quasar_h5_attr(fid, path, list(`encoding-type`    = "dataframe",
                                 `encoding-version` = "0.2.0",
                                 `_index`           = "_index",
                                 `column-order`     = as.character(cols)))
}

quasar_h5_write_ints <- function(fid, path, x) {
  if (length(x) && max(x) > .Machine$integer.max)
    stop("Sparse matrix exceeds 32-bit index range (>2^31-1 nnz or dimension); ",
         "not supported by quasar_h5adexporter.", call. = FALSE)
  rhdf5::h5write(as.integer(x), fid, path)
}


quasar_h5_write_X <- function(fid, mat, n_obs, n_vars) {
  if (inherits(mat, "sparseMatrix")) {
    mat <- methods::as(mat, "CsparseMatrix")
    rhdf5::h5createGroup(fid, "X")
    rhdf5::h5write(as.numeric(mat@x), fid, "X/data")
    quasar_h5_write_ints(fid, "X/indices", mat@i)    # feature idx -> CSR col idx
    quasar_h5_write_ints(fid, "X/indptr",  mat@p)    # per-observation pointers
    quasar_h5_attr(fid, "X", list(`encoding-type`    = "csr_matrix",
                                  `encoding-version` = "0.1.0",
                                  shape = as.integer(c(n_obs, n_vars))))
  } else {
    mat <- as.matrix(mat)                            # R (var,obs) -> py (obs,var)
    rhdf5::h5write(mat, fid, "X")
    quasar_h5_attr(fid, "X", list(`encoding-type` = "array",
                                  `encoding-version` = "0.2.0"))
  }
}

quasar_write_h5ad <- function(b, filename, overwrite) {
  if (file.exists(filename)) {
    if (!overwrite) stop("File exists: ", filename, " (set overwrite=TRUE).")
    unlink(filename)
  }
  rhdf5::h5createFile(filename)
  fid <- rhdf5::H5Fopen(filename)
  on.exit({ try(rhdf5::H5Fclose(fid), silent = TRUE); rhdf5::h5closeAll() },
          add = TRUE)

  quasar_h5_write_X(fid, b$X_vxo, b$n_obs, b$n_vars)
  quasar_h5_write_dataframe(fid, "obs", b$obs, rownames(b$obs))
  quasar_h5_write_dataframe(fid, "var", b$var, rownames(b$var))

  if (length(b$obsm)) {
    rhdf5::h5createGroup(fid, "obsm")
    for (nn in names(b$obsm)) {
      m <- as.matrix(b$obsm[[nn]])                   # R (k,obs) -> py (obs,k)
      rhdf5::h5write(t(m), fid, paste0("obsm/", nn))
      quasar_h5_attr(fid, paste0("obsm/", nn),
                     list(`encoding-type` = "array", `encoding-version` = "0.2.0"))
    }
  }

  quasar_h5_attr(fid, "/", list(`encoding-type`    = "anndata",
                                `encoding-version` = "0.1.0"))
}