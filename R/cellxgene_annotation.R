## custorm the standardize gene function
## We want to solve the problem that can't standatize gene name of data matrix if we only have counts matrix
## And in this case, the output will have the counts matrix which is as same as the data matrix.
StandardizeGeneSymbols_customer = function(
  obj,
  assay = NULL,
  layers = c("counts", "data"),
  EnsemblGeneTable = NULL,
  EnsemblGeneFile = NULL
) {
  if (is.null(assay)) {
    assay <- DefaultAssay(obj)
  }
  #If file is given
  if (is.null(EnsemblGeneTable)) {
    if (is.null(EnsemblGeneFile)) {
      stop("Please provide EnsemblID table or file")
    }
    EnsemblGeneTable <- fread(EnsemblGeneFile)
  }

  #Translate Ensembl IDs if necessary

  genes.in <- rownames(GetAssayData(obj, assay = assay, layer = layers[1]))
  ngenes <- length(genes.in)

  ens.count <- length(intersect(genes.in, EnsemblGeneTable[["Gene stable ID"]]))
  gname.count <- length(intersect(genes.in, EnsemblGeneTable[["Gene name"]]))

  ncbi.count <- 0
  if ("NCBI gene (formerly Entrezgene) ID" %in% colnames(EnsemblGeneTable)) {
    ncbi.count <- length(intersect(
      genes.in,
      EnsemblGeneTable[["NCBI gene (formerly Entrezgene) ID"]]
    ))
  }

  max <- max(ens.count, gname.count, ncbi.count)
  if (max < length(genes.in) / 2) {
    warning(
      "Over 50% of genes in input object not found in reference gene table"
    )
  }

  gname.format <- FALSE
  if (max == gname.count) {
    gname.format <- TRUE
  }

  if (max == ens.count) {
    #Input object has Ensembl IDs
    to = "Gene name"
    from = "Gene stable ID"
    genes.tr <- EnsemblGeneTable[[to]][match(
      genes.in,
      EnsemblGeneTable[[from]]
    )]
    names(genes.tr) <- genes.in
    genes.tr <- genes.tr[!is.na(genes.tr) & genes.tr != ""]
  } else if (max == ncbi.count) {
    to = "Gene name"
    from = "NCBI gene (formerly Entrezgene) ID"
    genes.tr <- EnsemblGeneTable[[to]][match(
      genes.in,
      EnsemblGeneTable[[from]]
    )]
    names(genes.tr) <- genes.in
    genes.tr <- genes.tr[!is.na(genes.tr) & genes.tr != ""]
  } else {
    genes.tr <- genes.in
    names(genes.tr) <- genes.in
  }

  ###### 1. First match dictionary
  geneRef_dict <- EnsemblGeneTable[["Gene name"]]
  names(geneRef_dict) <- EnsemblGeneTable[["Gene Synonym"]]
  geneRef_dict <- geneRef_dict[!is.null(names(geneRef_dict))]

  message(paste("Number of genes in input object:", ngenes))
  genesAllowList1 <- genes.tr[
    !is.na(genes.tr) &
      genes.tr != "" &
      genes.tr %in% EnsemblGeneTable[["Gene name"]]
  ] #keep genes with standard Gene.name
  l <- length(genesAllowList1)

  message(sprintf(
    "Number of genes with standard symbols: %i (%.2f%%)",
    l,
    l / ngenes * 100
  ))

  if (l < ngenes & gname.format) {
    message(paste("Examples of non-standard Gene.names:"))
    ns <- head(genes.tr[!genes.tr %in% EnsemblGeneTable[["Gene name"]]])
    message(paste(unname(ns), collapse = ","))
  }

  ###### 2. Search among synonyms
  genesAllowList2 <- genes.tr[
    !genes.tr %in% EnsemblGeneTable[["Gene name"]] &
      genes.tr %in% EnsemblGeneTable[["Gene Synonym"]]
  ] # keep genes with accepted Gene.name synonym
  genesAllowList2.gn <- geneRef_dict[genesAllowList2] # translate Gene.Synonym to standard Gene.name

  message(paste(
    "Additional number of genes with accepted Gene name synonym: ",
    length(genesAllowList2.gn)
  ))

  ##### 2b. Search by replacing dash (-) with dot (.)
  genes.dash <- gsub("\\.", "-", genes.tr)
  genesAllowList2b <- genes.tr[
    !genes.tr %in% EnsemblGeneTable[["Gene name"]] &
      genes.dash %in% EnsemblGeneTable[["Gene name"]]
  ]
  genesAllowList2b.gn <- gsub("\\.", "-", genesAllowList2b)
  message(paste(
    "Additional number of genes after replacing dots: ",
    length(genesAllowList2b.gn)
  ))

  #Names of genesAllowList contain IDs in matrix - elements contain the new names
  genesAllowList <- c(genesAllowList1, genesAllowList2.gn, genesAllowList2b.gn)

  ###### 3. Check for duplicates
  is.dup <- duplicated(genesAllowList)
  genesAllowList <- genesAllowList[!is.dup]
  message(sprintf(
    "Number of duplicated Gene.name: %i (%.2f%%)",
    sum(is.dup),
    sum(is.dup) / ngenes * 100
  ))

  l <- length(genesAllowList)
  message(sprintf("Final number of genes: %i (%.2f%%)", l, l / ngenes * 100))

  ###### 4. Subset matrix for allowed genes, and translate names
  matrix <- list()
  for (s in layers) {
    if (length(rownames(GetAssayData(obj, assay = assay, layer = s))) == 0) {
      next
    }
    matrix[[s]] <- GetAssayData(obj, assay = assay, layer = s)
    rows.select <- rownames(matrix[[s]])[
      rownames(matrix[[s]]) %in% names(genesAllowList)
    ]
    matrix[[s]] <- matrix[[s]][rows.select, ]
    rownames(matrix[[s]]) <- unname(genesAllowList[rows.select])
  }
  for (s in layers) {
    if (s == "counts" || s == "data") {
      obj <- suppressWarnings(RenameAssays(
        obj,
        assay.name = assay,
        new.assay.name = "tmp"
      ))
      obj[[assay]] <- CreateAssayObject(counts = matrix[[s]], assay = assay)
      DefaultAssay(obj) <- assay
      obj[["tmp"]] <- NULL
    } else {
      obj <- SetAssayData(obj, assay = assay, new.data = matrix[[s]], layer = s)
    }
  }

  return(obj)
}


cellxgene_annotation <- function(scData, param) {
  if (
    !ezIsSpecified(param$cellxgeneUrl) || !ezIsSpecified(param$cellxgeneLabel)
  ) {
    return(NULL)
  }
  # run cellxgene_annotation

  library(Seurat)
  library(sctransform)
  library(qs)
  library(STACAS)
  data(EnsemblGeneTable.Hs)
  data(EnsemblGeneTable.Mm)
  library(SeuratData)
  library(ggplot2)
  library(stringr)
  library(httr)
  library(harmony)
  library(schard)
  #library(glmGamPoi)

  cache_dir = "/srv/GT/databases/scRefData/CellxGene"
  # cache_dir = "/scratch/yang/tmp"
  cell_label_author = param$cellxgeneLabel
  species <- sub("/.*", "", param$refBuild)

  scRef <- getCuratedCellxGeneRef(
    param$cellxgeneUrl,
    cache_dir = cache_dir,
    cell_label_author = cell_label_author,
    species = species
  )

  ### StandardizeGeneSymbols
  if (species == "Homo_sapiens") {
    scData <- StandardizeGeneSymbols_customer(
      scData,
      layers = c("counts"),
      EnsemblGeneTable = EnsemblGeneTable.Hs
    )
  } else if (species == "Mus_musculus") {
    scData <- StandardizeGeneSymbols_customer(
      scData,
      layers = c("counts"),
      EnsemblGeneTable = EnsemblGeneTable.Mm
    )
  } else {
    stop(
      "We only support mouse and human dataset for using cellxgene annotation"
    )
  }
  scData <- NormalizeData(scData)
  ## mapping
  if ("harmony2" %in% names(scRef@reductions)) {
    scData.anchors <- FindTransferAnchors(
      reference = scRef,
      query = scData,
      dims = 1:30,
      reference.reduction = "harmony2",
      normalization.method = "LogNormalize"
    )
  } else {
    scData.anchors <- FindTransferAnchors(
      reference = scRef,
      query = scData,
      dims = 1:30,
      reference.reduction = "pca",
      normalization.method = "LogNormalize"
    )
  }

  if (!cell_label_author %in% colnames(scRef@meta.data)) {
    stop(
      "The specified column name for cell labels does not exist in the reference object's metadata."
    )
  }

  predictions.scData <- TransferData(
    anchorset = scData.anchors,
    refdata = scRef@meta.data[[cell_label_author]],
    dims = 1:30
  )
  cellxgeneResults <- data.frame(
    predicted.id.cellxgene.authorlabel = predictions.scData$predicted.id,
    prediction.score.max = predictions.scData$prediction.score.max,
    row.names = colnames(scData)
  )
  return(cellxgeneResults)
}


getCuratedCellxGeneRef <- function(
  ref_dataset_id,
  cache_dir,
  cell_label_author,
  species
) {
  print("start cellxgene")
  lockFile <- paste0(
    cache_dir,
    "/",
    gsub("\\.rds$", "", basename(ref_dataset_id)),
    "__",
    cell_label_author,
    ".lock"
  )
  refData_building_timeout_minutes <- 120

  i <- 0
  while (file.exists(lockFile) && i < refData_building_timeout_minutes) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste(
      "reference building still in progress after",
      refData_building_timeout_minutes,
      "min"
    ))
  }
  cached_curated_ref_data <- sub(".lock$", "-curated.qsd", lockFile)

  if (file.exists(cached_curated_ref_data)) {
    scRef <- qs::qread(cached_curated_ref_data)
    return(scRef)
  }

  ezWrite(Sys.info(), con = lockFile)
  on.exit(file.remove(lockFile), add = TRUE)

  ## Download reference dataset
  timeout_sec <- 3600
  tmp_download_ref <- paste0(cache_dir, "/", basename(ref_dataset_id))
  response <- httr::GET(
    ref_dataset_id,
    write_disk(tmp_download_ref, overwrite = TRUE),
    timeout(timeout_sec)
  )
  if (http_status(response)$category == "Success") {
    ezLog("Download completed successfully.")
  } else {
    ezLog(paste("Download failed. Status code:", http_status(response)$status))
    stop("Failed to download reference dataset")
  }

  # Detect file format and load accordingly
  file_extension <- tolower(tools::file_ext(ref_dataset_id))

  curated_seurat_object <- switch(
    file_extension,
    "rds" = {
      message("Loading RDS format Seurat object...")
      readRDS(tmp_download_ref)
    },
    "h5ad" = {
      message("Loading H5AD format and converting to Seurat object...")
      if (!requireNamespace("schard", quietly = TRUE)) {
        stop(
          "Package 'schard' is needed for h5ad format. Please install it first."
        )
      }
      schard::h5ad2seurat(tmp_download_ref)
    },
    stop(paste(
      "Unsupported file format:",
      file_extension,
      ". Only .rds and .h5ad formats are supported."
    ))
  )

  file.remove(tmp_download_ref)

  ### Standardize the ref dataset gene symbols with STACAS
  if (species == "Homo_sapiens") {
    if (
      "counts" %in%
        names(curated_seurat_object@assays$RNA@data) &&
        nrow(curated_seurat_object@assays$RNA@counts) > 0
    ) {
      # If counts exit
      curated_seurat_object <- StandardizeGeneSymbols_customer(
        curated_seurat_object,
        layers = c("counts"),
        EnsemblGeneTable = EnsemblGeneTable.Hs
      )
    } else {
      # If there is no counts matix
      curated_seurat_object <- StandardizeGeneSymbols_customer(
        curated_seurat_object,
        layers = c("data"),
        EnsemblGeneTable = EnsemblGeneTable.Hs
      )
    }
  } else if (species == "Mus_musculus") {
    if (
      "counts" %in%
        names(curated_seurat_object@assays$RNA@data) &&
        nrow(curated_seurat_object@assays$RNA@counts) > 0
    ) {
      # If counts exit
      curated_seurat_object <- StandardizeGeneSymbols_customer(
        curated_seurat_object,
        layers = c("counts"),
        EnsemblGeneTable = EnsemblGeneTable.Mm
      )
    } else {
      # If there is no counts matix
      curated_seurat_object <- StandardizeGeneSymbols_customer(
        curated_seurat_object,
        layers = c("data"),
        EnsemblGeneTable = EnsemblGeneTable.Mm
      )
    }
  } else {
    stop(
      "We only support mouse and human dataset for using cellxgene annotation"
    )
  }

  ## Downsample reference dataset
  ### split the object by donor_id

  curated_seurat_object.list <- SplitObject(
    curated_seurat_object,
    split.by = "donor_id"
  )

  rm(curated_seurat_object)

  ### choose the 10 biggest sample
  # calculate cell number of every sample
  cell_counts <- sapply(curated_seurat_object.list, ncol)
  print("Cell number of every donor:")
  print(cell_counts)
  if (!any(cell_counts > 500)) {
    stop(
      "Please choose a sample with a high enough number of cells, at least one donor_id with more than 500 cells."
    )
  }
  # number of samples need to be selected
  num_samples_to_select <- min(length(cell_counts), 10)

  # names of top10 samples
  selected_samples <- names(sort(cell_counts, decreasing = TRUE))[
    1:num_samples_to_select
  ]
  # get selected samples
  curated_seurat_object.list <- curated_seurat_object.list[selected_samples]

  # Remove the null object in the list
  curated_seurat_object.list <- curated_seurat_object.list[
    !sapply(curated_seurat_object.list, is.null)
  ]

  # delete sample with cell number smaller than 500
  curated_seurat_object.list <- curated_seurat_object.list[sapply(
    curated_seurat_object.list,
    function(x) ncol(x) >= 500
  )]

  #print("The info of the ref seurat object.list :")
  #print(head(curated_seurat_object.list[[1]]@meta.data))
  print("The column name of author's cell labels:")
  print(cell_label_author)
  print("Whether this column is in the ref data set:")
  print(
    cell_label_author %in% colnames(curated_seurat_object.list[[1]]@meta.data)
  )
  print("Selecting samples is done")
  print("Start to down sampling")

  ### cut the cell population to at most 3k
  set.seed(123)
  # check if the column exists
  if (
    !cell_label_author %in% colnames(curated_seurat_object.list[[1]]@meta.data)
  ) {
    stop(paste(
      "Column",
      cell_label_author,
      "does not exist in the Seurat object."
    ))
  }

  curated_seurat_object.list <- lapply(
    curated_seurat_object.list,
    function(sample) {
      # get the cell population of every sample
      cell_counts <- table(sample[[cell_label_author]])
      print("cell_counts:")
      print(cell_counts)

      #  Here, if the cell population less than 3000 we keep it. If the cell population bigger than 3000, we sample 3000 cells of total.
      sample_sizes <- ifelse(cell_counts < 3000, cell_counts, 3000)
      print("sample_sizes:")
      print(sample_sizes)
      sampled_cells <- unlist(lapply(names(sample_sizes), function(ct) {
        #cells_in_type <- WhichCells(sample,  expression = get(cell_label_author) == ct)
        cells_in_type <- Cells(sample)[
          (which(sample@meta.data[[cell_label_author]] %in% ct))
        ]
        print("cells_in_type:")
        print(head(cells_in_type))
        print(paste(
          "Cell type:",
          ct,
          "Number of cells:",
          length(cells_in_type)
        ))
        sample(cells_in_type, size = sample_sizes[ct], replace = FALSE)
      }))

      print("Sampled Cells:")
      print(head(sampled_cells))

      if (length(sampled_cells) == 0) {
        stop(
          "No cells were sampled. Please check your sample_sizes and cell labels."
        )
      }

      # check every assay
      assay_names <- Assays(sample)
      print("Available assays:")
      print(assay_names)

      # and info of every assay
      if (length(assay_names) == 1) {
        # go out from if clause
        message("Only one assay name is present. No other assay name.")
      } else {
        for (assay_name in assay_names) {
          assay_data <- GetAssayData(sample, assay = assay_name)
          print(paste("Assay name:", assay_name))
          print(paste("Number of features:", nrow(assay_data)))
          print(paste("Number of cells:", ncol(assay_data)))
          print(paste("Key for assay:", Key(sample[[assay_name]])))
        }
      }

      # THe default assay
      print(paste("Original default assay:", DefaultAssay(sample)))

      subset(sample, cells = sampled_cells)
    }
  )

  ## due to the otherwise huge cost of memory usage use harmony method
  scRef <- Reduce(function(x, y) merge(x, y), curated_seurat_object.list)
  scRef <- NormalizeData(scRef)
  scRef <- FindVariableFeatures(
    scRef,
    selection.method = "vst",
    nfeatures = 2000
  )
  #scRef <- SCTransform(scRef,assay = "RNA", new.assay.name = "SCT")
  scRef <- ScaleData(scRef)
  scRef <- RunPCA(scRef, features = VariableFeatures(scRef))
  if (length(unique(scRef$donor_id)) > 1) {
    scRef <- RunHarmony(scRef, "donor_id", assay.use = "RNA")
    scRef <- RunUMAP(scRef, reduction = "harmony", dims = 1:30)
    scRef[['harmony2']] <- CreateDimReducObject(
      embeddings = scRef[['harmony']]@cell.embeddings,
      key = "harmony2_",
      loadings = scRef[['pca']]@feature.loadings,
      assay = "RNA"
    )
    qs::qsave(scRef, cached_curated_ref_data)
  } else {
    scRef <- RunUMAP(scRef, reduction = "pca", dims = 1:30)
    qs::qsave(scRef, cached_curated_ref_data)
  }

  return(scRef)
}
