## Build mouse brain RCTD references from the Allen Brain Cell Atlas
## (Yao et al. 2023, Nature 624:317; WMB-10Xv3 raw counts, release 20230630,
## cell metadata 20241115).
##
## Three references are written:
##   allen_forebrain_layered_rctd.rds  35 types: ABC classes without pons,
##     medulla and cerebellum; split into the types 55 um spots can resolve
##     (upper vs deep cortical IT, hippocampal CA/SUB, glia, vascular)
##   allen_wholebrain_layered_rctd.rds 41 types: as above with hindbrain and
##     cerebellum, for sections that contain them
##   allen_wholebrain_class_rctd.rds   the 34 ABC classes
##
## Finer, subclass-level references (~50 types: individual cortical layers,
## CA1/CA3/DG) were tested on Visium SD FFPE mouse brain and rejected: the many
## near-collinear glutamatergic types made RCTD (full mode) assign a generic
## excitatory profile to one sink type (CA2-CA3 Glut, 18% mean weight, high in
## cortex and striatum). Hindbrain types attracted weight in the hypothalamus
## of a coronal forebrain section, hence the forebrain variant.
##
## Usage:
##   Rscript buildAllenWholeBrainRctdRef.R <downloadDir> <outDir> [nPerType]
## <downloadDir> must contain cell_metadata_with_cluster_annotation.csv and the
## WMB-10Xv3-<region>-raw.h5ad files, all from
## https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com
library(data.table)
library(hdf5r)
library(Matrix)

args <- commandArgs(trailingOnly = TRUE)
downloadDir <- args[1]
outDir <- args[2]
nPerType <- if (length(args) >= 3) as.integer(args[3]) else 300
seed <- 1
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

## subclass (3-digit id) -> label; unlisted subclasses keep their class.
## Only the splits that 55 um spots can resolve: hippocampal CA/SUB vs
## cortical IT-ET, upper vs deep cortical IT, glia and vascular subtypes.
splitMap <- c(
  "016" = "CA-SUB Glut", "017" = "CA-SUB Glut", "023" = "CA-SUB Glut",
  "025" = "CA-SUB Glut",
  "316" = "Astro", "317" = "Astro", "318" = "Astro", "319" = "Astro",
  "320" = "Astro",
  "321" = "Ependymal", "322" = "Ependymal", "323" = "Ependymal",
  "324" = "Ependymal",
  "325" = "Choroid plexus",
  "326" = "OPC", "327" = "Oligo",
  "334" = "Microglia",
  "335" = "BAM-immune", "336" = "BAM-immune", "337" = "BAM-immune",
  "338" = "BAM-immune"
)
layeredMap <- c(
  splitMap,
  "007" = "IT CTX upper", "008" = "IT CTX upper", "009" = "IT CTX upper",
  "011" = "IT CTX upper", "018" = "IT CTX upper", "019" = "IT CTX upper",
  "020" = "IT CTX upper", "006" = "IT CTX upper", "021" = "IT CTX upper",
  "003" = "IT-ET CTX deep", "004" = "IT-ET CTX deep", "005" = "IT-ET CTX deep",
  "022" = "IT-ET CTX deep", "024" = "IT-ET CTX deep"
)

meta <- fread(
  file.path(downloadDir, "cell_metadata_with_cluster_annotation.csv"),
  select = c(
    "cell_label", "library_method", "feature_matrix_label",
    "region_of_interest_acronym", "class", "subclass"
  )
)
meta <- meta[library_method == "10Xv3"]
h5Files <- file.path(
  downloadDir,
  paste0(unique(meta$feature_matrix_label), "-raw.h5ad")
)
names(h5Files) <- unique(meta$feature_matrix_label)
h5Files <- h5Files[file.exists(h5Files)]
meta <- meta[feature_matrix_label %in% names(h5Files)]
message("cells available: ", nrow(meta), " in ", length(h5Files), " matrices")

subclassId <- substr(meta$subclass, 1, 3)
meta[, classLabel := sub("^[0-9]+ ", "", class)]
meta[, layered := ifelse(subclassId %in% names(layeredMap), layeredMap[subclassId], classLabel)]

sample_per <- function(dt, col, n) {
  set.seed(seed)
  dt[, .SD[sample(.N, min(.N, n))], by = col]$cell_label
}
pickClass <- sample_per(meta, "classLabel", nPerType)
pickLayered <- sample_per(meta, "layered", nPerType)
picked <- unique(c(pickClass, pickLayered))
message("cells to read: ", length(picked))

read_rows <- function(h5File, cells) {
  h <- H5File$new(h5File, mode = "r")
  on.exit(h$close_all())
  labels <- h[["obs/cell_label"]][]
  rows <- sort(match(cells, labels))
  rows <- rows[!is.na(rows)]
  genes <- h[["var/gene_identifier"]][]
  indptr <- h[["X/indptr"]][]
  dataDs <- h[["X/data"]]
  indDs <- h[["X/indices"]]
  i <- list()
  x <- list()
  for (k in seq_along(rows)) {
    r <- rows[k]
    from <- indptr[r] + 1
    to <- indptr[r + 1]
    if (to >= from) {
      i[[k]] <- indDs[from:to] + 1
      x[[k]] <- dataDs[from:to]
    } else {
      i[[k]] <- integer(0)
      x[[k]] <- numeric(0)
    }
  }
  m <- sparseMatrix(
    i = unlist(i),
    j = rep(seq_along(rows), lengths(i)),
    x = round(unlist(x)),
    dims = c(length(genes), length(rows)),
    dimnames = list(genes, labels[rows])
  )
  m
}

mats <- list()
for (nm in names(h5Files)) {
  cells <- meta[feature_matrix_label == nm & cell_label %in% picked, cell_label]
  if (length(cells) == 0) {
    next
  }
  message("reading ", length(cells), " cells from ", nm)
  mats[[nm]] <- read_rows(h5Files[[nm]], cells)
}
genes <- rownames(mats[[1]])
stopifnot(all(vapply(mats, function(m) identical(rownames(m), genes), TRUE)))
counts <- do.call(cbind, mats)
counts <- counts[Matrix::rowSums(counts) > 0, ]
meta <- meta[match(colnames(counts), cell_label)]
message("matrix: ", nrow(counts), " genes x ", ncol(counts), " cells")

build_ref <- function(cells, labelCol, file) {
  idx <- match(cells, colnames(counts))
  idx <- idx[!is.na(idx)]
  ## spacexr forbids "/" in cell-type names
  ct <- factor(gsub("/", "-", meta[[labelCol]][idx]))
  names(ct) <- colnames(counts)[idx]
  m <- counts[, idx]
  ref <- spacexr::Reference(
    m,
    ct,
    nUMI = Matrix::colSums(m),
    min_UMI = 100,
    n_max_cells = nPerType
  )
  saveRDS(ref, file.path(outDir, file))
  tab <- as.data.frame(table(cell_type = ref@cell_types))
  fwrite(tab, file.path(outDir, sub("\\.rds$", "_celltypes.tsv", file)), sep = "\t")
  message("wrote ", file, ": ", nlevels(ref@cell_types), " types, ", length(ref@cell_types), " cells")
}
hindbrain <- c("P Glut", "P GABA", "MY Glut", "MY GABA", "CB Glut", "CB GABA")
build_ref(pickClass, "classLabel", "allen_wholebrain_class_rctd.rds")
build_ref(pickLayered, "layered", "allen_wholebrain_layered_rctd.rds")
build_ref(
  intersect(pickLayered, meta[!layered %in% hindbrain, cell_label]),
  "layered",
  "allen_forebrain_layered_rctd.rds"
)

mapping <- unique(meta[, .(class, subclass, layered)])[order(subclass)]
fwrite(mapping, file.path(outDir, "allen_brain_layered_mapping.tsv"), sep = "\t")
writeLines(c(
  "Mouse brain RCTD references (spacexr Reference objects).",
  "Source: Allen Brain Cell Atlas, Yao et al. 2023 (Nature 624:317-332),",
  "  WMB-10Xv3 raw counts (expression_matrices 20230630),",
  "  cell metadata with cluster annotation (metadata/WMB-10X 20241115).",
  "  https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com",
  "Genes: Ensembl gene IDs (match via features.tsv gene_id).",
  sprintf(
    "Sampling: up to %d 10Xv3 cells per type, seed %d, min_UMI 100.",
    nPerType,
    seed
  ),
  "allen_forebrain_layered_rctd.rds: 35 types, no pons/medulla/cerebellum;",
  "  for coronal forebrain sections (e.g. 10x Visium mouse brain).",
  "allen_wholebrain_layered_rctd.rds: 41 types incl. hindbrain and cerebellum.",
  "allen_wholebrain_class_rctd.rds: the 34 ABC classes.",
  "Labels per subclass: allen_brain_layered_mapping.tsv.",
  "Subclass-level references (~50 types) were rejected: one hippocampal type",
  "  became a sink for generic excitatory signal in Visium SD FFPE spots.",
  paste("Built", format(Sys.Date()), "with ezRun inst/scripts/buildAllenWholeBrainRctdRef.R")
), file.path(outDir, "README_allen_brain.txt"))
