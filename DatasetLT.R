
.ensure_packages <- function(pkgs) {
  stopifnot(is.character(pkgs), length(pkgs) > 0)
  
  need <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(need)) {
    message("Installing missing package(s): ", paste(need, collapse = ", "))
    install.packages(need)
  }
  
  invisible(lapply(pkgs, library, character.only = TRUE, quietly = TRUE))
}

## Make sure R6 is ready before we use it ―–––––––––––––––––––––––––––––
.ensure_packages("R6")

#' DatasetLT – Lightweight Multi-Layer Dataset Container
#'
#' @title DatasetLT
#' @description
#' `DatasetLT` is an **R6** class for managing multiple *assays*
#' (numeric matrices) together with optional embeddings, graphs and
#' **named per-assay column metadata**.  
#' The class enforces strict consistency across all layers:
#'
#' * **Column metadata** must match (set-equality) an assay’s columns on
#'   creation, and is *auto-reordered* to the assay’s column order.
#' * Replacing an assay attempts to keep matching metadata and **drops**
#'   mismatches with a warning.
#' * All row-wise operations propagate to every layer; operations are
#'   halted if any layer would lose required samples.
#'
#' @section Public “adder” methods:
#' * [`addAssay()`] – add one or many assays (matrices or data.frames).
#' * [`addEmbedding()`] – add embedding matrices per assay.
#' * [`addGraph()`] – add named square graphs per assay.
#' * [`addColumnMetadata()`] – add/replace per-assay named column
#'   metadata (strict match + auto-reorder).
#'
#' @section Public “getter” methods:
#' * [`getAssay()`], [`getEmbedding()`], [`getGraph()`] – retrieve
#'   layers (respecting active filters unless `force = TRUE`).
#' * [`getColumnMetadata()`] – retrieve column metadata with optional
#'   variable subset.
#'
#' @section Filtering helpers:
#' * `setActiveAssays()`, `setActiveSamples()`, `clearFilters()`
#' * `filterSamples()`, `cleanup()`
#'
#' @section I/O helpers:
#' * `writeCsv()` – **one CSV per assay**. Writes either the *active
#'   intersection* (default `intersect = TRUE`) or the *unfiltered*
#'   stored matrices (`intersect = FALSE`). Sample IDs are written as the
#'   first column in each CSV.
#'
#' @section Fields:
#' * **assays**: named list of numeric matrices  
#' * **embeddings**: list( assay → list( name → matrix ) )  
#' * **graphs**: list( assay → list( name → square matrix ) )  
#' * **columnMetadata**: list( assay → list( name → data.frame ) )  
#' * **activeAssays / activeSamples / hasActiveFlag**: filter state  
#'
#' @note
#' The implementation purposefully avoids dependencies other than
#' **R6**, keeping the package lightweight and self-contained.
#'
#' @seealso [`R6::R6Class`]
#' @export
#' @import R6
DatasetLT <- R6::R6Class(
  "DatasetLT",
  
  public = list(
    
    ## ----------------------------------------------------------------
    ## FIELDS ----------------------------------------------------------
    ## ----------------------------------------------------------------
    assays           = NULL,   # list( assay -> numeric matrix )
    embeddings       = NULL,   # list( assay -> list( name -> matrix ) )
    graphs           = NULL,   # list( assay -> list( name -> square matrix ) )
    columnMetadata   = NULL,   # list( assay -> list( name -> data.frame ) )
    activeAssays     = character(),
    activeSamples    = character(),
    hasActiveFlag    = FALSE,
    
    ## ----------------------------------------------------------------
    ## CONSTRUCTOR -----------------------------------------------------
    ## ----------------------------------------------------------------
    
    #' @description Create a new, empty `DatasetLT` (optionally with an
    #'   initial named list of assays).
    #' @param assays `NULL` **or** a *named* list of numeric matrices /
    #'   data.frames (row = samples, col = features).
    initialize = function(assays = NULL) {
      self$assays         <- list()
      self$embeddings     <- list()
      self$graphs         <- list()
      self$columnMetadata <- list()
      
      if (!is.null(assays)) {
        stopifnot(
          is.list(assays),
          !is.null(base::names(assays)),
          all(base::names(assays) != "")
        )
        self$addAssay(names = base::names(assays), input = assays)
      }
    },
    
    ## ----------------------------------------------------------------
    ## ADDERS ----------------------------------------------------------
    ## ----------------------------------------------------------------
    
    #' @description Add **one or many assays**.
    #' @param names Character vector of assay names (optional when
    #'   `input` is a *named list*).
    #' @param input A numeric matrix/data.frame *or* a list of such.
    #' @return Invisibly, the modified object (`self`).
    addAssay = function(names, input) {
      if (private$.is_bulk(input)) {
        keys <- private$.resolve_names_param(names, input, what = "assay list")
        private$.assert_unique(keys, what = "assay names")
        
        for (i in seq_along(input)) {
          key <- keys[i]
          mat <- private$toNumericMatrix(
            input[[i]],
            what = sprintf("assay '%s'", key)
          )
          private$.assert_colnames(mat, what = sprintf("assay '%s'", key))
          private$.assign_assay_and_reconcile_colmeta(key, mat)
        }
        return(invisible(self))
      }
      
      stopifnot(length(names) == 1, is.character(names))
      key <- names[1]
      mat <- private$toNumericMatrix(
        input, what = sprintf("assay '%s'", key)
      )
      private$.assert_colnames(mat, what = sprintf("assay '%s'", key))
      private$.assign_assay_and_reconcile_colmeta(key, mat)
      invisible(self)
    },
    
    #' @description Add **one or many embeddings** under a single assay.
    #' @inheritParams addAssay
    #' @param assay Name of the assay that owns the embeddings.
    addEmbedding = function(assay, input, names=NULL) {
      stopifnot(
        is.character(assay), length(assay) == 1,
        assay %in% base::names(self$assays)
      )
      
      if (is.null(self$embeddings[[assay]]))
        self$embeddings[[assay]] <- list()
      
      if (private$.is_bulk(input)) {
        keys <- private$.resolve_names_param(
          names, input,
          what = sprintf("embedding list for assay '%s'", assay)
        )
        private$.assert_unique(keys, what = "embedding names")
        for (i in seq_along(input)) {
          mat <- private$toNumericMatrix(
            input[[i]],
            what = sprintf("embedding '%s' (assay '%s')", keys[i], assay)
          )
          if (!is.null(self$embeddings[[assay]][[keys[i]]])) {
            warning(
              sprintf("Embedding '%s' for assay '%s' existed and was replaced.",
                      keys[i], assay),
              call. = FALSE
            )
          }
          self$embeddings[[assay]][[keys[i]]] <- mat
        }
        return(invisible(self))
      }
      
      stopifnot(length(names) == 1, is.character(names))
      key <- names[1]
      mat <- private$toNumericMatrix(
        input, what = sprintf("embedding '%s' (assay '%s')", key, assay)
      )
      if (!is.null(self$embeddings[[assay]][[key]])) {
        warning(
          sprintf("Embedding '%s' for assay '%s' existed and was replaced.",
                  key, assay),
          call. = FALSE
        )
      }
      self$embeddings[[assay]][[key]] <- mat
      invisible(self)
    },
    
    #' @description Add **one or many graphs** (square adjacency /
    #'   distance matrices) under a single assay.
    #' @inheritParams addEmbedding
    addGraph = function(assay,  input, names=NULL) {
      stopifnot(
        is.character(assay), length(assay) == 1,
        assay %in% base::names(self$assays)
      )
      
      if (is.null(self$graphs[[assay]]))
        self$graphs[[assay]] <- list()
      
      check_graph <- function(mat, nm) {
        if (nrow(mat) != ncol(mat))
          stop(sprintf("Graph '%s' (assay '%s') must be square: %d x %d.",
                       nm, assay, nrow(mat), ncol(mat)))
        if (is.null(rownames(mat)) || is.null(colnames(mat)))
          stop(sprintf("Graph '%s' (assay '%s') lacks row/col names.", nm, assay))
        if (!identical(rownames(mat), colnames(mat)))
          stop(sprintf("Graph '%s' (assay '%s') row/col names must match and be in the same order.",
                       nm, assay))
        if (anyDuplicated(rownames(mat)))
          stop(sprintf("Graph '%s' (assay '%s') has duplicate row/col names.", nm, assay))
      }
      
      if (private$.is_bulk(input)) {
        keys <- private$.resolve_names_param(
          names, input,
          what = sprintf("graph list for assay '%s'", assay)
        )
        private$.assert_unique(keys, what = "graph names")
        for (i in seq_along(input)) {
          mat <- private$toNumericMatrix(
            input[[i]],
            what = sprintf("graph '%s' (assay '%s')", keys[i], assay)
          )
          check_graph(mat, keys[i])
          if (!is.null(self$graphs[[assay]][[keys[i]]])) {
            warning(
              sprintf("Graph '%s' for assay '%s' existed and was replaced.",
                      keys[i], assay),
              call. = FALSE
            )
          }
          self$graphs[[assay]][[keys[i]]] <- mat
        }
        return(invisible(self))
      }
      
      stopifnot(length(names) == 1, is.character(names))
      key <- names[1]
      mat <- private$toNumericMatrix(
        input, what = sprintf("graph '%s' (assay '%s')", key, assay)
      )
      check_graph(mat, key)
      if (!is.null(self$graphs[[assay]][[key]])) {
        warning(
          sprintf("Graph '%s' for assay '%s' existed and was replaced.",
                  key, assay),
          call. = FALSE
        )
      }
      self$graphs[[assay]][[key]] <- mat
      invisible(self)
      },
    ## ----------------------------- SUMMARY ----------------------------------
    #' @description
    #' Print a concise inventory of **all assays, embeddings and graphs**.
    #'
    #' @details
    #' The output is human-readable, but the function invisibly returns a
    #' list so you can capture the structure programmatically:
    #' ```
    #' out <- dslt$printLayerNames()
    #' str(out)
    #' ```
    #' @return (Invisibly) a named `list` with three elements
    #'   `assays`, `embeddings`, `graphs`.
    #' @examples
    #' dslt <- DatasetLT$new()
    #' dslt$addAssay("rna", matrix(0, 1, 1,
    #'                             dimnames = list("c1", "g1")))
    #' dslt$addEmbedding("rna", "pca",
    #'                   matrix(0, 1, 2, dimnames = list("c1", c("PC1","PC2"))))
    #' dslt$printLayerNames()
    printLayerNames = function() {
      assays  <- base::names(self$assays)
      
      embeds  <- lapply(self$embeddings, base::names)
      graphs  <- lapply(self$graphs,     base::names)
      
      cat("Assays:\n")
      if (length(assays)) {
        cat("  •", paste(assays, collapse = ", "), "\n")
      } else {
        cat("  (none)\n")
      }
      
      cat("Embeddings:\n")
      if (length(embeds)) {
        for (a in names(embeds)) {
          cat(sprintf("  • %s: %s\n", a,
                      if (length(embeds[[a]]))
                        paste(embeds[[a]], collapse = ", ")
                      else "(none)"))
        }
      } else {
        cat("  (none)\n")
      }
      
      cat("Graphs:\n")
      if (length(graphs)) {
        for (a in names(graphs)) {
          cat(sprintf("  • %s: %s\n", a,
                      if (length(graphs[[a]]))
                        paste(graphs[[a]], collapse = ", ")
                      else "(none)"))
        }
      } else {
        cat("  (none)\n")
      }
      
      invisible(list(assays = assays,
                     embeddings = embeds,
                     graphs = graphs))
    },

#' @description Add or replace **named column metadata** for an
#'   assay (strict matching + auto-reorder).
#' @param assay Assay name.
#' @param name  Metadata object name.
#' @param df    A `data.frame`/matrix whose **rownames equal** the
#'   assay’s column names.
addColumnMetadata = function(assay, name, df) {
  stopifnot(
    is.character(assay), length(assay) == 1,
    is.character(name),  length(name)  == 1,
    assay %in% base::names(self$assays)
  )
  df <- private$toDataFrameWithRownames(
    df,
    what = sprintf("columnMetadata '%s' for assay '%s'", name, assay)
  )
  feats <- colnames(self$assays[[assay]])
  miss  <- setdiff(feats, rownames(df))
  extra <- setdiff(rownames(df), feats)
  if (length(miss) || length(extra)) {
    msg <- c(
      if (length(miss))  sprintf("missing features: %s",  paste(miss,  collapse = ", ")),
      if (length(extra)) sprintf("extra features: %s",     paste(extra, collapse = ", "))
    )
    stop("columnMetadata '", name, "' does not match assay '",
         assay, "': ", paste(msg, collapse = "; "))
  }
  df <- df[feats, , drop = FALSE]  # auto-reorder
  if (is.null(self$columnMetadata[[assay]]))
    self$columnMetadata[[assay]] <- list()
  if (!is.null(self$columnMetadata[[assay]][[name]])) {
    warning(
      sprintf("columnMetadata '%s' for assay '%s' existed and was replaced.",
              name, assay),
      call. = FALSE
    )
  }
  self$columnMetadata[[assay]][[name]] <- df
  invisible(self)
},

## ----------------------------------------------------------------
## GETTERS ---------------------------------------------------------
## ----------------------------------------------------------------

#' @description Retrieve column metadata (ordered to assay columns).
#' @param vars Optional character vector of column names to return.
getColumnMetadata = function(assay, name, vars = NULL) {
  stopifnot(
    is.character(assay), length(assay) == 1,
    is.character(name),  length(name)  == 1,
    assay %in% base::names(self$assays)
  )
  if (is.null(self$columnMetadata[[assay]]) ||
      is.null(self$columnMetadata[[assay]][[name]]))
    stop("No columnMetadata named '", name,
         "' found for assay '", assay, "'.")
  df <- self$columnMetadata[[assay]][[name]]
  feats <- colnames(self$assays[[assay]])
  if (!setequal(rownames(df), feats))
    stop("Stored columnMetadata '", name,
         "' no longer matches assay '", assay,
         "' columns. Please re-add metadata.")
  df <- df[feats, , drop = FALSE]
  if (!is.null(vars)) {
    stopifnot(is.character(vars))
    miss <- setdiff(vars, colnames(df))
    if (length(miss))
      stop("Unknown metadata variable(s): ", paste(miss, collapse = ", "))
    df <- df[, vars, drop = FALSE]
  }
  df
},

#' @description List metadata names (optionally for one assay).
listColumnMetadata = function(assay = NULL) {
  if (is.null(assay)) lapply(self$columnMetadata, base::names)
  else base::names(self$columnMetadata[[assay]])
},

#' @description Remove a column metadata entry (convenience).
removeColumnMetadata = function(assay, name) {
  stopifnot(
    is.character(assay), length(assay) == 1,
    is.character(name),  length(name)  == 1
  )
  if (!is.null(self$columnMetadata[[assay]])) {
    self$columnMetadata[[assay]][[name]] <- NULL
    if (length(self$columnMetadata[[assay]]) == 0)
      self$columnMetadata[[assay]] <- NULL
  }
  invisible(self)
},

## ------------------------- FILTER HELPERS ------------------------

setActiveAssays = function(names) {
  stopifnot(is.character(names), length(names) > 0)
  miss <- setdiff(names, base::names(self$assays))
  if (length(miss))
    stop("Missing assay(s): ", paste(miss, collapse = ", "))
  private$.assert_unique(names, what = "active assay names")
  self$activeAssays <- names
  self$hasActiveFlag <- TRUE
  invisible(self)
},

setActiveSamples = function(samples) {
  stopifnot(is.character(samples))
  self$activeSamples <- samples
  invisible(self)
},

listAssays = function(active = FALSE) {
  if (active) {
    if (self$hasActiveFlag) self$activeAssays else base::names(self$assays)
  } else {
    base::names(self$assays)
  }
},

currentSamples = function(force = FALSE) {
  private$.base_samples(ignore_active_samples = force)
},

clearFilters = function() {
  self$activeAssays  <- character()
  self$activeSamples <- character()
  self$hasActiveFlag <- FALSE
  invisible(self)
},

## ------------------------------ GETTERS --------------------------

getAssay = function(name, force = FALSE) {
  private$checkActive(name, force)
  mat <- self$assays[[name]]
  if (is.null(mat)) stop("Assay '", name, "' not found.")
  mat[private$intersectSamples(mat), , drop = FALSE]
},

getEmbedding = function(assay, name, force = FALSE) {
  private$checkActive(assay, force)
  embed <- self$embeddings[[assay]][[name]]
  if (is.null(embed))
    stop("Embedding '", name, "' not found for assay '", assay, "'.")
  embed[private$intersectSamples(embed), , drop = FALSE]
},

getGraph = function(assay, name, force = FALSE) {
  private$checkActive(assay, force)
  g <- self$graphs[[assay]][[name]]
  if (is.null(g))
    stop("Graph '", name, "' not found for assay '", assay, "'.")
  idx <- private$intersectSamples(g)
  g[idx, idx, drop = FALSE]
},

## -------------------- CLEANUP & FILTERING -----------------------

cleanup = function(samples = NULL, persist = TRUE) {
  if (length(self$assays) == 0)
    stop("Cleanup aborted: no assays have been added.")
  
  keep <- Reduce(intersect, lapply(self$assays, rownames))
  if (!is.null(samples))
    keep <- intersect(keep, samples)
  if (length(keep) == 0)
    stop("Cleanup aborted: resulting sample set is empty.")
  
  # Validate embeddings/graphs cover 'keep'
  missing_msgs <- private$.validate_layers_cover(keep)
  if (length(missing_msgs)) {
    stop(paste0(
      "Cleanup aborted: some layers do not contain all surviving samples.\n",
      paste0(" - ", missing_msgs, collapse = "\n")
    ))
  }
  
  # Subset all layers by rows
  self$assays <- lapply(self$assays, \(m) m[keep, , drop = FALSE])
  for (a in base::names(self$embeddings))
    self$embeddings[[a]] <- lapply(self$embeddings[[a]], \(m) m[keep, , drop = FALSE])
  for (a in base::names(self$graphs))
    self$graphs[[a]] <- lapply(self$graphs[[a]], \(g) g[keep, keep, drop = FALSE])
  
  if (persist) self$activeSamples <- keep
  message("Cleanup complete: ", length(keep), " samples retained.")
  invisible(self)
},

filterSamples = function(samples, persist = TRUE) {
  stopifnot(is.character(samples), length(samples) > 0)
  ord_subset <- \(m) {
    keep <- samples[samples %in% rownames(m)]
    m[keep, , drop = FALSE]
  }
  
  self$assays <- lapply(self$assays, ord_subset)
  for (a in base::names(self$embeddings))
    self$embeddings[[a]] <- lapply(self$embeddings[[a]], ord_subset)
  for (a in base::names(self$graphs))
    self$graphs[[a]] <- lapply(self$graphs[[a]], \(g) {
      keep <- samples[samples %in% rownames(g)]
      g[keep, keep, drop = FALSE]
    })
  
  retained <- unique(unlist(lapply(self$assays, rownames), use.names = FALSE))
  if (persist) self$activeSamples <- intersect(samples, retained)
  message("Filter complete: requested ", length(samples),
          " samples; ", length(retained), " present after filtering.")
  invisible(self)
},
## --------------------------- I/O : embeddings ---------------------------
#' @description
#' Write **every embedding matrix** to disk as one CSV per
#' *embedding*, with sample IDs in the first column.
#'
#' @details
#' * If `intersect = TRUE` (default) the active‐sample intersection
#'   is applied (same rows you’d get from `getEmbedding()`).
#' * If `intersect = FALSE`, the **unfiltered** stored matrix is
#'   written.
#' * Filenames follow the pattern  
#'   &nbsp;&nbsp;• **directory target**: `<dir>/<assay>_<embed>.csv`  
#'   &nbsp;&nbsp;• **prefix target**   : `<prefix>_<assay>_<embed>.csv`
#' * You may limit the export to specific assays and/or embeddings.
#'
#' @param file Character scalar: an existing/new directory **or**
#'   filename prefix.  If it ends in `.csv`, the extension is
#'   stripped and used as the prefix.
#' @param assays Optional character vector of assay names.  Defaults
#'   to active assays (when `intersect = TRUE`) or all assays.
#' @param embeddings Optional named list whose names are assay names
#'   and whose values are the embedding names to keep for that assay
#'   (character vector).  If `NULL`, **all** embeddings of the
#'   selected assays are written.
#' @param intersect Logical (`TRUE`, default) – respect active
#'   filters?
#' @param overwrite Logical – overwrite existing files?  Default
#'   `FALSE`.
#' @param na Character to use for missing values.  Default `"NA"`.
#' @return (Invisibly) a character vector of paths written.
writeEmbeddingsCsv = function(file,
                              assays      = NULL,
                              embeddings  = NULL,
                              intersect   = TRUE,
                              overwrite   = T,
                              na          = "NA") {
  stopifnot(is.character(file), length(file) == 1, nzchar(file))
  stopifnot(is.logical(intersect), length(intersect) == 1)
  stopifnot(is.logical(overwrite), length(overwrite) == 1)
  stopifnot(is.character(na), length(na) == 1)
  
  ## ---------- assays to export -----------------------------------------
  all_assays <- base::names(self$assays)
  if (is.null(assays)) {
    assays <- if (intersect) self$listAssays(active = TRUE) else all_assays
  } else {
    stopifnot(is.character(assays), length(assays) > 0)
    miss <- setdiff(assays, all_assays)
    if (length(miss))
      stop("Unknown assay(s): ", paste(miss, collapse = ", "))
  }
  if (length(assays) == 0)
    stop("No assays selected to write.")
  
  ## ---------- filter embedding sets ------------------------------------
  if (!is.null(embeddings)) {
    stopifnot(is.list(embeddings), all(names(embeddings) != ""))
    bad <- setdiff(names(embeddings), assays)
    if (length(bad))
      stop("Embedding filter supplied for non-selected assay(s): ",
           paste(bad, collapse = ", "))
  }
  
  ## ---------- output naming pattern ------------------------------------
  is_dir_target <- dir.exists(file) || grepl("[/\\\\]$", file)
  written_paths <- character()
  
  gen_path <- function(a, e) {
    if (is_dir_target) {
      outdir <- sub("[/\\\\]$", "", file)
      if (!dir.exists(outdir) && nzchar(outdir))
        dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
      file.path(if (nzchar(outdir)) outdir else ".",
                paste0(a, "_", e, ".csv"))
    } else {
      prefix <- sub("\\.csv$", "", file, ignore.case = TRUE)
      paste0(prefix, "_", a, "_", e, ".csv")
    }
  }
  
  ## ---------- write each selected embedding ----------------------------
  for (a in assays) {
    emb_list <- self$embeddings[[a]]
    if (is.null(emb_list))
      next                           # no embeddings to write
    
    embeds_to_write <- if (is.null(embeddings) || is.null(embeddings[[a]])) {
      names(emb_list)
    } else {
      sel <- embeddings[[a]]
      miss <- setdiff(sel, names(emb_list))
      if (length(miss))
        stop("Assay '", a, "' lacks embedding(s): ",
             paste(miss, collapse = ", "))
      sel
    }
    
    for (e in embeds_to_write) {
      mat <- if (intersect) {
        self$getEmbedding(a, e, force = FALSE)
      } else {
        emb_list[[e]]
      }
      
      target <- gen_path(a, e)
      if (file.exists(target) && !overwrite)
        stop("File exists and overwrite = FALSE: ", target)
      
      utils::write.csv(mat, file = target, row.names = TRUE, na = na)
      written_paths <- c(written_paths, target)
    }
  }
  
  invisible(written_paths)
},
## ----------------------------- I/O --------------------------------

#' @description
#' Write **one CSV per assay** to disk.
#'
#' @details
#' * If `intersect = TRUE` (default), each CSV contains the **active
#'   intersection** of samples across active assays (and respects
#'   `activeSamples`), i.e. the same rows returned by `getAssay()`.
#' * If `intersect = FALSE`, **no filtering** is applied; the full,
#'   stored matrix per assay is written.
#' * The **first column** in each CSV contains the **sample IDs**
#'   (row names).
#' * The `file` argument can be:
#'   - a **directory** path (existing or to be created) → files are
#'     written as `<dir>/<assay>.csv`;
#'   - a **file-like prefix** (e.g., `"out.csv"` or `"out"`) →
#'     files are written as `"prefix_<assay>.csv"`.
#'
#' @param file Character scalar. Output directory **or** filename
#'   prefix. If it ends with `.csv`, the extension is stripped and
#'   `_<assay>.csv` is appended.
#' @param assays Optional character vector of assay names to write.
#'   Defaults to active assays (if any) when `intersect = TRUE`,
#'   otherwise all assays.
#' @param intersect Logical (default `TRUE`). If `TRUE`, write the
#'   active intersection; if `FALSE`, write unfiltered stored
#'   matrices (ignoring active flags).
#' @param overwrite Logical. Overwrite existing files? Default `FALSE`.
#' @param na Character string used for missing values in CSV.
#'   Default `"NA"`.
#' @return (Invisibly) a named character vector of file paths that
#'   were written, indexed by assay name.
writeCsv = function(file,
                    assays = NULL,
                    intersect = TRUE,
                    overwrite = T,
                    na = "NA") {
  stopifnot(is.character(file), length(file) == 1, nzchar(file))
  stopifnot(is.logical(intersect), length(intersect) == 1)
  stopifnot(is.logical(overwrite), length(overwrite) == 1)
  stopifnot(is.character(na), length(na) == 1)
  
  # Determine assays to write
  all_assays <- base::names(self$assays)
  if (is.null(assays)) {
    assays <- if (intersect) self$listAssays(active = TRUE) else all_assays
  } else {
    stopifnot(is.character(assays), length(assays) > 0)
    miss <- setdiff(assays, all_assays)
    if (length(miss))
      stop("Unknown assay(s): ", paste(miss, collapse = ", "))
  }
  if (length(assays) == 0)
    stop("No assays selected to write.")
  
  # Resolve output pattern
  is_dir_target <- dir.exists(file) || grepl("[/\\\\]$", file)
  paths <- stats::setNames(character(length(assays)), assays)
  
  if (is_dir_target) {
    # Ensure directory exists (create if trailing slash given)
    outdir <- sub("[/\\\\]$", "", file)
    if (!dir.exists(outdir) && nzchar(outdir)) {
      dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
    }
    for (a in assays) {
      paths[[a]] <- file.path(if (nzchar(outdir)) outdir else ".", paste0(a, ".csv"))
    }
  } else {
    # Treat as filename prefix
    prefix <- sub("\\.csv$", "", file, ignore.case = TRUE)
    for (a in assays) {
      paths[[a]] <- paste0(prefix, "_", a, ".csv")
    }
  }
  
  # Write each assay
  for (a in assays) {
    mat <- if (intersect) {
      # respect active assays/samples intersection
      self$getAssay(a, force = FALSE)
    } else {
      # raw stored matrix
      self$assays[[a]]
    }
    
    if (is.null(mat))
      stop("Assay '", a, "' has no data to write.")
    
    target <- paths[[a]]
    if (file.exists(target) && !overwrite) {
      stop("File exists and overwrite = FALSE: ", target)
    }

    utils::write.csv(mat, file = target, row.names = TRUE, na = na)
  }
  
  invisible(paths)
}
  ),

private = list(
  
  ## --------------------------- HELPERS ---------------------------
  
  toNumericMatrix = function(obj, what = "matrix") {
    if (is.data.frame(obj)) {
      if (!all(vapply(obj, is.numeric, logical(1))))
        stop(what, " contains non-numeric columns.")
      obj <- as.matrix(obj)
    }
    if (!is.matrix(obj))
      stop(what, " must be a matrix or data.frame.")
    if (!is.numeric(obj))
      stop(what, " must be numeric.")
    if (is.null(rownames(obj)))
      stop(what, " lacks row names.")
    if (anyDuplicated(rownames(obj)))
      stop(what, " has duplicate row names (sample IDs).")
    obj
  },
  
  .assert_colnames = function(mat, what = "matrix") {
    if (is.null(colnames(mat)))
      stop(what, " lacks column names (feature IDs).")
    if (anyDuplicated(colnames(mat)))
      stop(what, " has duplicate column names (feature IDs).")
    invisible(TRUE)
  },
  
  toDataFrameWithRownames = function(obj, what = "data.frame") {
    if (is.matrix(obj))
      obj <- as.data.frame(obj, stringsAsFactors = FALSE)
    if (!is.data.frame(obj))
      stop(what, " must be a data.frame or matrix.")
    if (is.null(rownames(obj)))
      stop(what, " lacks row names (feature IDs).")
    if (anyDuplicated(rownames(obj)))
      stop(what, " has duplicate row names (feature IDs).")
    obj
  },
  
  .is_bulk = function(x) {
    is.list(x) && !is.data.frame(x)
  },
  
  .resolve_names_param = function(names, input, what = "list") {
    if (!is.null(names)) {
      stopifnot(is.character(names), length(names) == length(input))
      return(names)
    }
    keys <- base::names(input)
    if (is.null(keys) || any(keys == ""))
      stop(what, " must be a named list (or provide 'names').")
    keys
  },
  
  .assert_unique = function(keys, what = "names") {
    if (anyDuplicated(keys))
      stop("Duplicate ", what, " are not allowed: ",
           paste(unique(keys[duplicated(keys)]), collapse = ", "))
  },
  
  .base_samples = function(ignore_active_samples = FALSE) {
    if (length(self$assays) == 0)
      stop("No assays have been added yet.")
    
    assays_to_use <- if (self$hasActiveFlag && length(self$activeAssays) > 0) {
      self$activeAssays
    } else {
      base::names(self$assays)
    }
    
    if (length(assays_to_use) == 0)
      stop("No active assays are set. Call setActiveAssays() or add assays.")
    
    base <- Reduce(intersect, lapply(assays_to_use, \(a) rownames(self$assays[[a]])))
    if (!ignore_active_samples && length(self$activeSamples))
      base <- intersect(base, self$activeSamples)
    base
  },
  
  intersectSamples = function(layer) {
    base <- private$.base_samples(ignore_active_samples = FALSE)
    absent <- setdiff(base, rownames(layer))
    if (length(absent))
      stop("Layer lacks some intersected samples: ", paste(absent, collapse = ", "))
    base
  },
  
  checkActive = function(assay, force) {
    active <- if (self$hasActiveFlag) self$activeAssays else base::names(self$assays)
    if (!(assay %in% active) && !force)
      stop("Assay '", assay, "' is not active (use force = TRUE).")
  },
  
  .validate_layers_cover = function(keep) {
    msgs <- character()
    
    for (a in base::names(self$embeddings)) {
      for (nm in base::names(self$embeddings[[a]])) {
        m <- self$embeddings[[a]][[nm]]
        miss <- setdiff(keep, rownames(m))
        if (length(miss))
          msgs <- c(msgs, sprintf("Embedding '%s' (assay '%s') missing: %s",
                                  nm, a, paste(miss, collapse = ", ")))
      }
    }
    
    for (a in base::names(self$graphs)) {
      for (nm in base::names(self$graphs[[a]])) {
        g <- self$graphs[[a]][[nm]]
        miss <- setdiff(keep, rownames(g))
        if (length(miss))
          msgs <- c(msgs, sprintf("Graph '%s' (assay '%s') missing: %s",
                                  nm, a, paste(miss, collapse = ", ")))
      }
    }
    
    msgs
  },
  
  .assign_assay_and_reconcile_colmeta = function(key, mat) {
    if (!is.null(self$assays[[key]])) {
      warning(sprintf("Assay '%s' existed and was replaced.", key), call. = FALSE)
      if (!is.null(self$columnMetadata[[key]])) {
        new_cols <- colnames(mat)
        kept <- 0L; dropped <- character()
        for (nm in base::names(self$columnMetadata[[key]])) {
          df <- self$columnMetadata[[key]][[nm]]
          if (!is.null(df) && setequal(rownames(df), new_cols)) {
            self$columnMetadata[[key]][[nm]] <- df[new_cols, , drop = FALSE]
            kept <- kept + 1L
          } else {
            self$columnMetadata[[key]][[nm]] <- NULL
            dropped <- c(dropped, nm)
          }
        }
        if (length(dropped)) {
          warning(sprintf(
            "columnMetadata for assay '%s' removed due to feature set mismatch: %s",
            key, paste(dropped, collapse = ", ")
          ), call. = FALSE)
        }
        if (kept == 0L && length(dropped))
          if (is.null(base::names(self$columnMetadata[[key]])) ||
              length(self$columnMetadata[[key]]) == 0)
            self$columnMetadata[[key]] <- NULL
      }
    }
    self$assays[[key]] <- mat
    invisible(NULL)
  }
)
)
