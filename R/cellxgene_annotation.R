# 
# ezIsSpecified = function(x){
#   !is.null(x) && length(x) > 0 && x[1] != "" && !is.na(x[1]) && x[1] != "NA"
# }
# 
# ezWrite = function(..., sep="", collapse=" ", con=stdout()){
#   args = list(...)  ## see function message
#   #args = sapply(args, print)# as.character)
#   text = paste(sapply(args, paste, collapse=collapse), collapse=sep)
#   writeLines(text, con=con)
# }
# 
# 
# param <- list(cellxgene ='https://datasets.cellxgene.cziscience.com/99dc51ce-83ae-4f1f-ae3c-89f3509168fc.rds', column_name_of_cell_label = 'BICCN_subclass_label')
# #param.test.2 <- list(cellxgene ='https://datasets.cellxgene.cziscience.com/d39144df-fa59-4b63-b07b-9b34613b5c84.rds', column_name_of_cell_label = 'Manually_curated_celltype',refBuild = 'Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_42-2023-01-30')
# system.time({scData <- UpdateSeuratObject(LoadData("pbmc3k"))})
# test.result <- cellxgene_annotation(scData = scData, param = param)
# 
# cache_dir = "/scratch/yang/tmp"
# data(EnsemblGeneTable.Mm)
# scRef <- getCuratedCellxGeneRef(param$cellxgene, cache_dir=cache_dir, cell_label_author = param$column_name_of_cell_label, species = 'Mus_musculus')
# 

cellxgene_annotation <- function(scData, param) {
  

  if (!ezIsSpecified(param$cellxgene) || !ezIsSpecified(param$column_name_of_cell_label)){
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
  #library(glmGamPoi)
  

  cache_dir = "/srv/GT/databases/scRefData/CellxGene"
  #cache_dir = "/scratch/yang/tmp"
  cell_label_author = param$column_name_of_cell_label
  species <- sub("/.*", "", param$refBuild)
  
  scRef <- getCuratedCellxGeneRef(param$cellxgene, cache_dir=cache_dir, cell_label_author = cell_label_author, species = species)
  
  ### StandardizeGeneSymbols
  if( species == "Homo_sapiens" ){
    scData <- StandardizeGeneSymbols(scData, slots = c( "counts"), EnsemblGeneTable = EnsemblGeneTable.Hs)
  }else if(species == "Mus_musculus"){
    scData <- StandardizeGeneSymbols(scData, slots = c( "counts"), EnsemblGeneTable = EnsemblGeneTable.Mm)
  }else{
    stop("We only support mouse and human dataset for using cellxgene annotation")
  }
  
  ## mapping
  scData.anchors <- FindTransferAnchors(reference = scRef, query = scData, dims = 1:30,
                                        reference.reduction = "pca", normalization.method = "SCT" )
  
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

  
  
getCuratedCellxGeneRef <- function(ref_dataset_id, cache_dir, cell_label_author, species){

  lockFile <- paste0(cache_dir, "/", gsub("\\.rds$", "", basename(ref_dataset_id)), "__", cell_label_author, ".lock")
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
  cached_curated_ref_data <- sub(".lock$", "-curated.qds", lockFile)
  
  if (file.exists(cached_curated_ref_data)) {
    scRef <- qs::qread(cached_curated_ref_data)
    return(scRef)
  }
  
  ezWrite(Sys.info(), con = lockFile)
  on.exit(file.remove(lockFile), add = TRUE)
  
  
  
  ## get the processed ref data
  ### Download panceas dataset, ref one
  
  timeout_sec <- 3600
  tmp_download_ref <- paste0(cache_dir, "/", basename(ref_dataset_id))
  response <- httr::GET(ref_dataset_id, write_disk(tmp_download_ref, overwrite = TRUE), timeout(timeout_sec))
  if (http_status(response)$category == "Success") {
    cat("Download completed successfully.\n")
  } else {
    cat("Download failed. Status code:", http_status(response)$status, "\n")
  }
  curated_seurat_object <- readRDS(tmp_download_ref)
  file.remove(tmp_download_ref)
  


  ### Standardize the ref dataset gene symbols with STACAS
  if( species == "Homo_sapiens" ){
    curated_seurat_object <- StandardizeGeneSymbols(curated_seurat_object,slots = c("counts"), EnsemblGeneTable = EnsemblGeneTable.Hs)
  }else if(species == "Mus_musculus"){
    curated_seurat_object <- StandardizeGeneSymbols(curated_seurat_object,slots = c("counts"), EnsemblGeneTable = EnsemblGeneTable.Mm)
  }else{
    stop("We only support mouse and human dataset for using cellxgene annotation")
  }
  
  
  ## Downsample reference dataset
  ### split the object by donor_id

  curated_seurat_object.list <- SplitObject(curated_seurat_object, split.by = "donor_id")
  
  #print("The info of the ref seurat object:")
  #print(head(curated_seurat_object@meta.data))
  
  rm(curated_seurat_object)
  
  ### choose the 10 biggest sample 
  # calculate cell number of every sample
  cell_counts <- sapply(curated_seurat_object.list, ncol)
  print("Cell number of every donor:")
  print(cell_counts)
  if (!any(cell_counts > 500)) {
    stop("Please choose a sample with a high enough number of cells, at least one donor_id with more than 500 cells.")
  }
  # number of samples need to be selected
  num_samples_to_select <- min(length(cell_counts), 10)
  
  # names of top10 samples
  selected_samples <- names(sort(cell_counts, decreasing = TRUE))[1:num_samples_to_select]
  # get selected samples
  curated_seurat_object.list <- curated_seurat_object.list[selected_samples]
  # delete sample with cell number smaller than 500
  curated_seurat_object.list <- curated_seurat_object.list[sapply(curated_seurat_object.list, function(x) ncol(x) >= 500)]
  
  #print("The info of the ref seurat object.list :")
  #print(head(curated_seurat_object.list[[1]]@meta.data))
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
  features <- SelectIntegrationFeatures(object.list = curated_seurat_object.list, nfeatures = 2000)
  curated_seurat_object.list <- PrepSCTIntegration(object.list = curated_seurat_object.list, anchor.features = features)
  curated_seurat_object.list <- lapply(X = curated_seurat_object.list, FUN = RunPCA, features = features)
  
  
  
  
  if ( length(curated_seurat_object.list) > 1){
    anchors <- FindIntegrationAnchors(object.list = curated_seurat_object.list, normalization.method = "SCT",
                                      anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
    
    seurat.combined.sct <- IntegrateData(anchorset = anchors, normalization.method = "SCT", dims = 1:30)
    
    seurat.combined.sct <- RunPCA(seurat.combined.sct, verbose = FALSE)
    scRef <- RunUMAP(seurat.combined.sct, reduction = "pca", dims = 1:30)
    qs::qsave(scRef,cached_curated_ref_data)
  } else {
    # If there is only one element in the list(only one donor)
    scRef <- curated_seurat_object.list[[1]]
    scRef <- RunPCA(scRef, verbose = FALSE)
    scRef <- RunUMAP(scRef, reduction = "pca", dims = 1:30)
    qs::qsave(scRef,cached_curated_ref_data)
  }


  return(scRef)
}

