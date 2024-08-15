
cellxgene_annotation <- function(scData, param) {
  
  # param <- list(cellxgene ='71be997d-ff75-41b9-8a9f-1288c865f921', column_name_of_cell_label = 'Manually_curated_celltype')
  # # ##param.test.2 <- list(cellxgene ='37b21763-7f0f-41ae-9001-60bad6e2841d', column_name_of_cell_label = 'Manually_curated_celltype')
  # system.time({scData <- UpdateSeuratObject(LoadData("pbmc3k"))})
  # test.result <- cellxgene_annotation(scData = scData, param = param)
  
  if (!ezIsSpecified(param$cellxgene) || !ezIsSpecified(param$column_name_of_cell_label)){
    return(NULL)
  }
  # run cellxgene_annotation
  
  library(duckplyr)
  library(CuratedAtlasQueryR)
  library(SingleCellExperiment)
  #library(rlang)
  library(Seurat)
  library(sctransform)
  library(qs)
  library(STACAS)
  data(EnsemblGeneTable.Hs)
  library(SeuratData)
  library(ggplot2)
  library(stringr)
  #library(glmGamPoi)
  
  
  cell_label_author <- param$column_name_of_cell_label
  print(cell_label_author)
  ref_dataset_id <- param$cellxgene
  cache_dir = "/srv/GT/databases/scRefData/CellxGene"
  
  lockFile <- paste0(cache_dir, "/", ref_dataset_id, ".lock")
  refData_building_timeout_minutes <- 120

  i <- 0
  while (file.exists(lockFile) && i < refData_building_timeout_minutes) {
    ### somebody else builds and we wait
    Sys.sleep(60)
    i <- i + 1
  }
  if (file.exists(lockFile)) {
    stop(paste("reference building still in progress after", refData_building_timeout_minutes, "min"))
  }
  cached_curated_ref_data <- paste0(ref_dataset_id, "-curated.qsd")
  
  if (file.exists(cached_curated_ref_data)) {
    scRef <- qs::qread(cached_curated_ref_data)
  } else {
    scRef <- buildCuratedCellxGeneRef(ref_dataset_id, cached_dir=cache_dir)
    qs::qsave(cached_curated_ref_data)
  }
  
  ### StandardizeGeneSymbols
  scData <- StandardizeGeneSymbols(scData, slots = c( "counts"), EnsemblGeneTable = EnsemblGeneTable.Hs)
  
  # ## Do I need to do SCTransform here? or the scData has alrady done by SCTransform?  No, don't need to .
  # scData <- SCTransform(scData)
  
  ## mapping
  scData.anchors <- FindTransferAnchors(reference = scRef, query = scData, dims = 1:30,
                                        reference.reduction = "pca")
  
  if (!cell_label_author %in% colnames(scRef@meta.data)) {
    stop("The specified column name for cell labels does not exist in the reference object's metadata.")
  }
  
  predictions.scData <- TransferData(anchorset = scData.anchors, refdata = scRef@meta.data[[cell_label_author]], dims = 1:30)
  cellxgeneResults <- data.frame(
    predicted.id.cellxgene.authorlabel=predictions.scData$predicted.id,
    prediction.score.max=predictions.scData$prediction.score.max,
    row.names=colnames(scData))
  return(cellxgeneResults)
  
}

  
  
buildCuratedCellxGeneRef <- function(ref_dataset_id, cached_dir=cache_dir){

  
  ## get the unharmonised meta data
  metadata <- get_metadata(cache_directory = cached_dir)
  curated_single_cell_experiment_object <- metadata |>
    dplyr::filter(
      dataset_id  == param$cellxgene
    ) 
  
  print("The meta data of the data set:")
  print(curated_single_cell_experiment_object)
  
  
  ### check whether get a proper metadata
  if (is.null(head(dplyr::pull(.data = curated_single_cell_experiment_object, var = 1),1))) {
    stop("Failed to get unharmonised metadata. Please select a correct data set id. Please be noticed that do not select from collections.")
  }
  
  unharmonised_metadata <- get_unharmonised_metadata(curated_single_cell_experiment_object, cache_directory=cached_dir)
  
  ### get the author version cell labels
  dplyr::pull(unharmonised_metadata) |> head(2)
  df <- unharmonised_metadata$unharmonised[[1]]
  ### check whether this dataset have the donor_id column
  donor_id_exists <- if ("donor_id" %in% colnames(df)) TRUE else FALSE
  if (donor_id_exists) {
    df2 <- df |> 
      dplyr::select(cell_, cell_label_author, donor_id) |> 
      collect()
  }else{
    df2 <- df |> 
      dplyr::select(cell_, cell_label_author) |> 
      collect() |> 
      mutate(sample_id = str_extract(cell_, "-\\d+$") %>% str_remove("-"))
  }
  
  
  ## get the processed ref data
  ### Download panceas dataset, ref one
  
  curated_seurat_object <- metadata |>
    dplyr::filter(
      dataset_id == ref_dataset_id
    ) |>
    get_seurat()
  ### Standardize the ref dataset gene symbols with STACAS
  curated_seurat_object <- StandardizeGeneSymbols(curated_seurat_object,slots = c("counts"), EnsemblGeneTable = EnsemblGeneTable.Hs)
  ### merge unharmonised meta data with seurat object
  metadata <- curated_seurat_object@meta.data
  metadata2 <- metadata %>%
    left_join(df2, by = c("original_cell_id" = "cell_"))
  rownames(metadata2) <- rownames(metadata)
  curated_seurat_object@meta.data <- metadata2
  
  ## Downsample reference dataset
  ### split the object by sample
  curated_seurat_object[["RNA"]] <- curated_seurat_object[["originalexp"]]
  if(donor_id_exists){
    curated_seurat_object.list <- SplitObject(curated_seurat_object, split.by = "donor_id")
  }else{
    curated_seurat_object.list <- SplitObject(curated_seurat_object, split.by = "sample_id")
  }
  
  
  print("The info of the ref seurat object:")
  print(head(curated_seurat_object@meta.data))
  
  rm(curated_seurat_object)
  
  ### choose the 10 biggest sample (test 15 samples and crushed)
  # calculate cell number of every sample
  cell_counts <- sapply(curated_seurat_object.list, ncol)
  print("Cell number of every sample:")
  print(cell_counts)
  # number of samples need to be selected
  num_samples_to_select <- min(length(cell_counts), 10)
  
  # names of top10 samples
  selected_samples <- names(sort(cell_counts, decreasing = TRUE))[1:num_samples_to_select]
  # get selected samples
  curated_seurat_object.list <- curated_seurat_object.list[selected_samples]
  # delete sample with cell number smaller than 500
  curated_seurat_object.list <- curated_seurat_object.list[sapply(curated_seurat_object.list, function(x) ncol(x) >= 500)]
  
  print("The info of the ref seurat object.list :")
  print(head(curated_seurat_object.list[[1]]@meta.data))
  print("The column name of author's cell labels:")
  print(cell_label_author)
  print("Whether this column is in the ref data set:")
  print(cell_label_author %in% colnames(curated_seurat_object.list[[1]]@meta.data))
  print("Selecting samples is done")
  print("Start to down sampling")
  
  ### cut the cell population to at most 3k 
  set.seed(123)
  # check if the column exists
  if (!cell_label_author %in% colnames(curated_seurat_object.list[[1]]@meta.data)) {
    stop(paste("Column", cell_label_author, "does not exist in the Seurat object."))
  }
  
  curated_seurat_object.list <- lapply(curated_seurat_object.list, function(sample) {
    
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
      cells_in_type <- Cells(sample)[(which(sample@meta.data[[cell_label_author]] %in% ct))]
      print("cells_in_type:")
      print(head(cells_in_type))
      print(paste("Cell type:", ct, "Number of cells:", length(cells_in_type)))
      sample(cells_in_type, size = sample_sizes[ct], replace = FALSE)
      
    }))
    
    print("Sampled Cells:")
    print(head(sampled_cells))
    
    if (length(sampled_cells) == 0) {
      stop("No cells were sampled. Please check your sample_sizes and cell labels.")
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
    
  })
  
  ### Performing integration on datasets normalized with SCTransform
  curated_seurat_object.list <- lapply(X = curated_seurat_object.list, FUN = SCTransform, method = "glmGamPoi")
  features <- SelectIntegrationFeatures(object.list = curated_seurat_object.list, nfeatures = 3000)
  curated_seurat_object.list <- PrepSCTIntegration(object.list = curated_seurat_object.list, anchor.features = features)
  curated_seurat_object.list <- lapply(X = curated_seurat_object.list, FUN = RunPCA, features = features)
  
  anchors <- FindIntegrationAnchors(object.list = curated_seurat_object.list, normalization.method = "SCT",
                                    anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
  seurat.combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
  
  seurat.combined.sct <- RunPCA(seurat.combined.sct, verbose = FALSE)
  seurat.combined.sct <- RunUMAP(seurat.combined.sct, reduction = "pca", dims = 1:30)
  return(seurat.combined.sct)
}

