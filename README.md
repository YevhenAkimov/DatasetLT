DatasetLT demo
==============
Yevhen Akimov  
2025-08-05  

DatasetLT – Lightweight Multi-Layer Dataset Container
=====================================================

`DatasetLT` is a **single-file, R6-based** helper that keeps any number of
numeric assays (RNA-seq, ATAC-seq, proteomics …) and their related
embeddings / graphs / feature-level metadata in one tidy object – no
hefty Bioconductor stack required.

Key features
------------

* Only **R6** is required – installed automatically on first use.  
* Strict row/column consistency (rows = samples, cols = features).  
* “Active” filters: work on one assay at a time without touching others.  
* Single-line CSV export – one file per layer.  

Installation
------------

```r
# install.packages("devtools")   # if needed
devtools::install_github("YevhenAkimov/DatasetLT")
```

The helper auto-installs **R6** if it is missing.

Quick start – minimal toy workflow
----------------------------------

```r
library(DatasetLT)      # loads & installs R6 automatically
set.seed(42)
path <- "/Users/yevhenakimov/export2/"   # change to any writable dir

## 1 ── Create two tiny assays -----------------------------------------
rna  <- matrix(rnorm(30),  nrow = 6,
               dimnames = list(paste0("cell", 1:6),
                               paste0("gene", 1:5)))

atac <- matrix(rpois(30, 5), nrow = 6,
               dimnames = list(paste0("cell", 1:6),
                               paste0("peak", 1:5)))

## 2 ── Instantiate the container with those assays --------------------
dslt <- DatasetLT$new(list(rna = rna, atac = atac))

## 3 ── Add embeddings for each assay ----------------------------------
dslt$addEmbedding(
  "rna",
  names = c("pca", "tsne"),
  input = list(
    pca  = matrix(rnorm(12), 6,
                  dimnames = list(paste0("cell", 1:6),
                                  paste0("PC", 1:2))),
    tsne = matrix(rnorm(12), 6,
                  dimnames = list(paste0("cell", 1:6),
                                  paste0("tSNE", 1:2)))
  )
)

dslt$addEmbedding(
  "atac",
  names = "lsi",
  input = matrix(rnorm(12), 6,
                 dimnames = list(paste0("cell", 1:6),
                                 paste0("LSI", 1:2)))
)

## 4 ── Inspect what we now have ---------------------------------------
dslt$printLayerNames()
#> Assays:
#>   • rna, atac
#> Embeddings:
#>   • rna: pca, tsne
#>   • atac: lsi
#> Graphs:
#>   (none)

## 5 ── Focus on RNA only ----------------------------------------------
dslt$setActiveAssays("rna")

## 6 ── Export the *active* layers -------------------------------------
#     creates 'export2/' if absent
dslt$writeCsv(path)               # → …/export2/rna.csv
dslt$writeEmbeddingsCsv(path)     # → …/export2/rna_pca.csv , rna_tsne.csv

## 7 ── Export *everything*, ignoring active filters -------------------
dslt$writeCsv(path,          intersect = FALSE)   # → rna.csv, atac.csv
dslt$writeEmbeddingsCsv(path, intersect = FALSE)  # → rna_pca.csv, rna_tsne.csv, atac_lsi.csv
```

How the filters work
--------------------

| overlay        | respected when `intersect = TRUE` | ignored when `intersect = FALSE` |
|----------------|-----------------------------------|-----------------------------------|
| **activeAssays** | yes                               | yes (all assays written)          |
| **activeSamples**| yes                               | yes (all stored rows written)     |

Rows physically removed by `filterSamples()` / `cleanup()` are gone
regardless of the `intersect` flag.

Documentation
-------------

* **Vignette** – `vignette("DatasetLT", package = "DatasetLT")`
* **Help** – `?DatasetLT` once the package is loaded.

Issues & contributions
----------------------

Open an issue or pull request at  
<https://github.com/YevhenAkimov/DatasetLT> – feedback is welcome!
