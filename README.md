# DatasetLT

Lightweight multi-layer dataset container (R6) for storing **assays** (numeric matrices) together with per-assay **embeddings**, **graphs**, and **feature-level column metadata** — with strict consistency rules.
Designed for storing multimodal lineage tracing experiments 

- Rows = **samples**, columns = **features**. Duplicate names are rejected.
- Active filters let you work on a subset without mutating stored data (unless you call the destructive helpers).
- CSV writers output **one file per assay/embedding**.

---
  
## Load
```r
source("https://raw.githubusercontent.com/YevhenAkimov/DatasetLT/main/DatasetLT.R")
```
---
  
## Lets generate demo data 
```r
barcodes   <- paste0("barcode", 1:6)
barcodes_drugs <- barcodes[1:4]
drug_responses <- paste0("drug_responses", 1:4)
FACS_subpopulations   <- paste0("FACS_subpopulations", 1:5)
growth_curves <- paste0("growth_curves", 1:6)
## barcodes selection for drugs
drug_responses_analysis  <- matrix(rnorm(length(barcodes_drugs) * length(drug_responses)),
               nrow = length(barcodes_drugs),
               dimnames = list(barcodes_drugs, drug_responses))
## barcodes selection upon FACS sorting
FACS_subpopulations_analysis <- matrix(rnorm(length(barcodes) * length(FACS_subpopulations)),
               nrow = length(barcodes),
               dimnames = list(barcodes, FACS_subpopulations))
growth_curves_analysis <- matrix(rnorm(length(barcodes) * length(growth_curves)),
                                       nrow = length(barcodes),
                                       dimnames = list(barcodes, growth_curves))


# --- generate demo embeddings (rownames must be sample IDs(eg barcodes)) -------------------------------

drug_pca  <- matrix(rnorm(length(barcodes) * 2),
                    nrow = length(barcodes),
                    dimnames = list(barcodes, c("PC1","PC2")))
drug_tsne <- matrix(rnorm(length(barcodes) * 2),
                    nrow = length(barcodes),
                    dimnames = list(barcodes, c("tSNE1","tSNE2")))
facs_pca <- matrix(rnorm(length(barcodes) * 2),
                   nrow = length(barcodes),
                   dimnames = list(barcodes, c("PC1","PC2")))

# --- generate demo graphs (square, rownames == colnames == sample IDs(eg barcodes)) -------------------
demo_graph=cor(t(drug_responses_analysis))

```

## Lets add demo data to DatasetLT container
```r
# add assays to container 
dslt <- DatasetLT$new(list(drugs = drug_responses_analysis, facs = FACS_subpopulations_analysis,growth=growth_curves_analysis))

## add embeddings to assays
# add drugs assay-related embeddings, if input is list - names embeddings will be set according to names of the list elements, or names can be specified
dslt$addEmbedding("drugs", input = list(pca = drug_pca, tsne = drug_tsne))
# add facs assay-related embeddings,  if input is matrix - name should be provided in names argument
dslt$addEmbedding("facs",          input = facs_pca,names = "pca"  )
# add  drugs assay-related graphs
dslt$addGraph("drugs", input = demo_graph, names = "correlation")
```
## Usage example
```r
# --- inspect dataset -------------------------------------------------------
dslt$printLayerNames()
# We can set active assays to work with, that ensure consistent output of data - only intersect and ordered  barcodes are returned when data is queried
dslt$setActiveAssays(c("drugs","facs") )  

## Get  facs assay
facs_mat   <- dslt$getAssay("facs")
## Get  drugs assay
drugs_mat <- dslt$getAssay("drugs")
# check that we have identical barcodes for both assays
identical(rownames(facs_mat),rownames(drugs_mat))
# [1] TRUE
# for original inputs:
identical(rownames(drug_responses_analysis),rownames(FACS_subpopulations_analysis))
#[1] FALSE

# Get pca embedding for facs assay
facs_pca   <- dslt$getEmbedding(assay="facs", name="pca")
# Get pca embedding for drugs assay
drugs_graph <- dslt$getGraph(assay="drugs", name="pca")



# --- active samples vs destructive filtering --------------------------------
```r
# Assume we have additional filters for samples (e.g., barcodes), by setting active samples we implement additional level of filtering
dslt$setActiveSamples(c("barcode2","barcode3","barcode5"))  # non-destructive preference
dslt$currentSamples()                               # returns sample IDs that are going to be returned on queries
#[1] "barcode2" "barcode3"
# "barcode5" was not included due to active assays - only samples present in both active assays and active samples are returned
facs_mat   <- dslt$getAssay("facs")
rownames(facs_mat)
#[1] "barcode2" "barcode3"

#  we can clear all filter flags. Tn this case all assays will be treated as active 
dslt$clearFilters()
# 

```
## write CSVs 
When intersect = TRUE (default), writers respect active assays + active samples.
When intersect = FALSE, writers dump stored matrices (ignore actives).
```r
# assays
dslt$setActiveAssays(c("facs"))
path1="/path_to_dir/"
path2="/path_to_dir/"

dslt$writeCsv(path1, intersect = TRUE,  overwrite = TRUE)  # active intersection
dslt$writeCsv(path2, intersect = FALSE, overwrite = TRUE)  # all stored matrices


# embeddings
dslt$writeEmbeddingsCsv(path1, intersect = TRUE,  overwrite = TRUE)   # active
dslt$writeEmbeddingsCsv(path2, intersect = FALSE, overwrite = TRUE)   # stored
## save pca for "facs" assay only

```

---
  
  ## What to Pass to Each “Adder”
  
  - **`addAssay(names, input)`**  
  - `input`: numeric matrix or numeric data.frame (rows = **samples** with unique rownames, cols = **features** with unique colnames).  
- Can also pass a *named list* of such matrices; `names` may be omitted if the list is named.

- **`addEmbedding(assay, names, input)`**  
  - Embedding matrices with **rownames == sample IDs** (any number of columns).  
- Single matrix or *named list* for multiple embeddings.

- **`addGraph(assay, names, input)`**  
  - Square matrix with **identical row/col names** (sample IDs, same order).  
- Single matrix or *named list* for multiple graphs.

- **`addColumnMetadata(assay, name, df)`**  
  - `df`: data.frame/matrix whose **rownames equal the assay’s column names** (features).  
- Auto-reordered to assay’s column order; replacing an assay drops mismatched metadata with a warning.

---
  

---
  
  ## Filters & Intersections
  
  - **Active assays** (`setActiveAssays()`): define the assays considered for getters (unless `force=TRUE`) and for `intersect=TRUE` exports.  
- **Active samples** (`setActiveSamples()`): preferred sample IDs; used when intersecting rows. **Non-destructive**.
- **Destructive helpers**:
  - `filterSamples(samples, persist=TRUE)`: subset stored rows to requested IDs (drops missing).
- `cleanup(samples=NULL, persist=TRUE)`: enforce strict across-assay intersection; **aborts** if some embedding/graph would lose required samples.

  Rows removed by `filterSamples()` / `cleanup()` are **physically gone** regardless of `intersect`.

---
  
  ## Notes
  
  - Assays must be **numeric** and have **unique** row/column names.
- Embeddings must **cover samples**; graphs must be **square** with matching row/col names.
- Feature metadata’s **rownames must match the assay’s column names** (set-equality); the class will auto-reorder.
- Replacing an assay tries to keep only matching column-metadata entries; **mismatches are dropped (warning)**.
- `cleanup()` will **abort** if any embedding/graph doesn’t have all surviving samples — fix inputs or remove those layers first.

---
  
