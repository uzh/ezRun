###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

EzAppSTARsolo = 
    setRefClass("EzAppSTARsolo",
            contains = "EzApp",
            methods = list(
                initialize = function()
                {
                    "Initializes the application using its specific defaults."
                    runMethod <<- ezMethodSTARsolo
                    name <<- "EzAppSTARsolo"
                    appDefaults <<- rbind(controlSeqs=ezFrame(Type="charVector",
                                                              DefaultValue="",
                                                              Description="control sequences to add"),
                                          ## STARsolo parameters
                                          soloType=ezFrame(Type="character",
                                                           DefaultValue="CB_UMI_Simple",
                                                           Description="CB_UMI_Simple (a.k.a. Droplet), CB_UMI_Complex (e.g. Droplet)."),
                                          soloCBwhitelist=ezFrame(Type="character",
                                                           DefaultValue="SC3Pv3",
                                                           Description="Barcode whitelist. Choose: SC3Pv1, SC3Pv2, SC3Pv3."),
                                          soloCBstart=ezFrame(Type="character",
                                                           DefaultValue="1",
                                                           Description="Cell barcode start base."),
                                          soloCBlen=ezFrame(Type="character",
                                                           DefaultValue="16",
                                                           Description="Cell barcode length."),
                                          soloUMIstart=ezFrame(Type="character",
                                                           DefaultValue="17",
                                                           Description="UMI start base."),
                                          soloUMIlen=ezFrame(Type="character",
                                                             DefaultValue="auto",
                                                             Description="UMI length. 10 if SC3Pv3 not selected, 12 otherwise."),
                                          soloUMIfiltering=ezFrame(Type="character",
                                                                   DefaultValue="MultiGeneUMI",
                                                                   Description="type of UMI filtering.
                                                                            -               ... basic filtering: remove UMIs with N and homopolymers (similar to CellRanger 2.2.0)
                                                                            MultiGeneUMI    ... remove lower-count UMIs that map to more than one gene (introduced in CellRanger 3.x.x)"),
                                          soloCBmatchWLtype=ezFrame(Type="character",
                                                                    DefaultValue="1MM_multi_pseudocounts",
                                                                    Description="matching the Cell Barcodes to the WhiteList.
                                                                                Exact                   only exact matches allowed
                                                                                1MM                     only one match in whitelist with 1 mismatched base allowed. Allowed CBs have to have at least one read with exact match.
                                                                                1MM_multi               multiple matches in whitelist with 1 mismatched base allowed, posterior probability calculation is used choose one of the matches.
                                                                                                        Allowed CBs have to have at least one read with exact match. Similar to CellRanger 2.2.0
                                                                                1MM_multi_pseudocounts  same as 1MM_Multi, but pseudocounts of 1 are added to all whitelist barcodes.
                                                                                                        Similar to CellRanger 3.x.x"),
                                          soloCellFilter=ezFrame(Type="character",
                                                                 DefaultValue="None",
                                                                 Description="cell filtering type and parameters
                                                                        CellRanger2.2   ... simple filtering of CellRanger 2.2, followed by thre numbers: number of expected cells, robust maximum percentile for UMI count, maximum to minimum ratio for UMI count
                                                                        TopCells        ... only report top cells by UMI count, followed by the excat number of cells
                                                                        None            ... do not output filtered cells")
                                          # soloBarcodeReadLength=ezFrame(Type="character",
                                          #                               DefaultValue="1",
                                          #                               Description="1 if equal to sum of soloCBlen+soloUMIlen; 0 if not defined, do not check."),
                                          # soloCBposition=ezFrame(Type="character",
                                          #                        DefaultValue="-",
                                          #                        Description="Position of Cell Barcode(s) on the barcode read.
                                          #                                   Presently only works with --soloType CB_UMI_Complex, and barcodes are assumed to be on Read2.
                                          #                                   Format for each barcode: startAnchor_startDistance_endAnchor_endDistance
                                          #                                   start(end)Anchor defines the anchor base for the CB: 
                                          #                                   0: read start; 1: read end; 2: adapter start; 3: adapter end
                                          #                                   start(end)Distance is the distance from the CB start(end) to the Anchor base
                                          #                                   String for different barcodes are separated by space.
                                          #                                   Example: inDrop (Zilionis et al, Nat. Protocols, 2017):
                                          #                                   --soloCBposition  0_0_2_-1  3_1_3_8"),
                                          # soloUMIposition=ezFrame(Type="character",
                                          #                         DefaultValue="-",
                                          #                         Description="position of the UMI on the barcode read, same as soloCBposition.
                                          #                                   Example: inDrop (Zilionis et al, Nat. Protocols, 2017):
                                          #                                   --soloCBposition  3_9_3_14"),
                                          # soloAdapterSequence=ezFrame(Type="character",
                                          #                               DefaultValue="-",
                                          #                               Description="Adapter sequence to anchor barcodes."),
                                          # soloAdapterMismatchesNmax=ezFrame(Type="character",
                                          #                             DefaultValue="1",
                                          #                             Description="maximum number of mismatches allowed in adapter sequence."),
                                          # soloStrand=ezFrame(Type="character",
                                          #                    DefaultValue="Forward",
                                          #                    Description="strandedness of the solo libraries.
                                          #                               Unstranded  ... no strand information
                                          #                               Forward     ... read strand same as the original RNA molecule
                                          #                               Reverse     ... read strand opposite to the original RNA molecule"),
                                          # soloFeatures=ezFrame(Type="character",
                                          #                    DefaultValue="Gene",
                                          #                    Description="genomic features for which the UMI counts per Cell Barcode are collected.
                                          #                               Gene            ... genes: reads match the gene transcript
                                          #                               SJ              ... splice junctions: reported in SJ.out.tab
                                          #                               GeneFull        ... full genes: count all reads overlapping genes' exons and introns
                                          #                               Transcript3p   ... quantification of transcript for 3' protocols"),
                                          # soloUMIdedup=ezFrame(Type="character",
                                          #                    DefaultValue="1MM_All",
                                          #                    Description="type of UMI deduplication (collapsing) algorithm.
                                          #                               1MM_All             ... all UMIs with 1 mismatch distance to each other are collapsed (i.e. counted once)
                                          #                               1MM_Directional     ... follows the 'directional' method from the UMI-tools by Smith, Heger and Sudbery (Genome Research 2017).
                                          #                               Exact               ... only exactly matching UMIs are collapsed"),
                                                            )
                                                          }
            )
)

ezMethodSTARsolo = function(input=NA, output=NA, param=NA){
    # analysis vars
    ## check if a new index is needed, then create it with CellRanger and the last version of STAR.
    refDir <- getSTARReference(param)
    sampleName = input$getNames()
    sampleDirs = strsplit(input$getColumn("RawDataDir"), ",")[[sampleName]]
    sampleDirs = file.path(input$dataRoot, sampleDirs)

    # create STARsolo command
    cmd = makeSTARsoloCmd(param, refDir, sampleName, sampleDirs)
    
    # run STARsolo shell command
    ezSystem(cmd)
    
    # filter raw counts with EmptyDrop
    require(Matrix)
    require(readr)
    require(DropletUtils)
    
    outputDir = paste0(sampleName,'/Solo.out/Gene/raw')
    ## Remove alignments file if not setted differently
    samFn = paste0(sampleName,'/Aligned.out.sam')
    if(param[['keepAlignment']]=='True'){
        # convert to BAM
        bamFn = paste0(sampleName,'/Aligned.out.bam')
        #samtoolsBin = '/usr/local/ngseq/bin/samtools'
        ezSystem(paste('samtools','view','-Sb',samFn,'>',bamFn))
        # remove .sam
        ezSystem(paste('rm',samFn))
        # sort .bam
        bamFnSorted = paste0(sampleName,'/Aligned.out.sorted.bam') 
        ezSortIndexBam(bamFn,bamFnSorted)
    }else{
        # remove .sam
        ezSystem(paste('rm',samFn))
    }
    ## Prepare STARsolo output to be processed as CellRanger output
    ### add "Gene Expression" column to features.tsv
    featuresFn = paste0(outputDir,'/features.tsv')
    tmp =  paste0(outputDir,'/tmp.tsv')
    cmd = paste("awk \'{print $0, \"\tGene Expression\"}\'",featuresFn,">",tmp)
    ezSystem(cmd)
    ezSystem(paste("mv",tmp,featuresFn)) # rename tmp into features
    ### gzip all files
    ezSystem(paste("pigz -p 4",paste0(outputDir,"/*")))

    ## Filter raw STARsolo counts with EmptyDrops
    sce = read10xCounts(samples = outputDir, col.names = TRUE, version = '3')
    rawCounts = sce@assays@data@listData$counts
    ### cell calling
    called = defaultDrops(rawCounts)
    ### Data of called cells.
    countsFiltered = rawCounts[,called]
    geneID = sce@rowRanges@elementMetadata@listData[["ID"]]
    geneSymbol = sce@rowRanges@elementMetadata@listData[["Symbol"]]
    ### save in new directory
    filteredDir = paste0(sampleName,'/filtered_feature_bc_matrix')
    write10xCounts(path = filteredDir, x = countsFiltered, gene.id = geneID, gene.symbol = geneSymbol, version = '3')
    
    # Get cell cycle from filtered output
    ## generate filtered single cell experiment object
    sceFilt = SingleCellExperiment(assays = list(counts = countsFiltered))
    rownames(sceFilt) = geneID
    ## analyse cell cycle phase
    cellPhase = getCellCycle(sceFilt, param$refBuild)
    write_tsv(cellPhase,
              path=file.path(filteredDir,
                "CellCyclePhase.txt"))
    
    return("Success")
}
# create STARsolo shell command
makeSTARsoloCmd = function(param, refDir, sampleName, sampleDirs){
    require('tools')
    # define binary full path
    STARbin = "/usr/local/ngseq/packages/Aligner/STAR/2.7.3a/bin/STAR"

    ## decide which chemistry whitelist to take
    soloCBwhitelist = list(
        SC3Pv1 = '/usr/local/ngseq/opt/cellranger-3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/737K-april-2014_rc.txt',
        SC3Pv2 = '/usr/local/ngseq/opt/cellranger-3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/737K-august-2016.txt',
        SC3Pv3 = '/usr/local/ngseq/opt/cellranger-3.1.0/cellranger-cs/3.1.0/lib/python/cellranger/barcodes/3M-february-2018.txt'
    )
    
    ## decide soloUMIlen
    if(param[['soloUMIlen']]=='auto'){
        if(param[['soloCBwhitelist']]=='SC3Pv3'){
            soloUMIlen = '12'
        }else{
            soloUMIlen = '10'
        }
    }else{
        soloUMIlen = param[['soloUMIlen']]
    }
    
    ## create readFilesIn
    ### list files: the names are always the same for standard runs. Only the presence of indexes is optional.
    fastqfiles = sort(list.files(sampleDirs,full.names = TRUE,pattern = '.fastq'))
    cDNAfastqs = paste(grep('R2',fastqfiles, value = TRUE), collapse = ',')
    barcodesfastqs = paste(grep('R1',fastqfiles, value = TRUE), collapse = ',')
    indexfastqs = paste(grep('_I',fastqfiles,value = TRUE), collapse = ',')
    
    readFilesIn = trimws(paste(cDNAfastqs, indexfastqs, barcodesfastqs, collapse = ' ')) # remove last space, incase indexfastqs is empty.
    
    ## create readFilesCommand
    if(unique(file_ext(fastqfiles))=='gz'){
        readFilesCommand = 'zcat'
    }else{
        readFilesCommand = 'cat'
    }
    
    # create output folder for the sample (in CellRanger this was automated)
    if(!dir.exists(sampleName)){dir.create(sampleName)}
    
    # create full STARsolo command
    cmd = paste(STARbin,
                ## STAR general parameters
                paste0('--runThreadN ',param[['cores']]),
                paste0('--outFileNamePrefix ',sampleName,'/'),
                
                paste0('--readFilesIn ',readFilesIn),
                paste0('--readFilesCommand ',readFilesCommand),
                paste0('--genomeDir ',refDir),
                
                ## STARsolo parameters
                paste0('--soloType ',param[['soloType']]),
                paste0('--soloCBwhitelist ',soloCBwhitelist[[param[['soloCBwhitelist']]]]),
                paste0('--soloCBstart ',param[['soloCBstart']]),
                paste0('--soloCBlen ',param[['soloCBlen']]),
                paste0('--soloUMIstart ',param[['soloUMIstart']]),
                paste0('--soloUMIlen ',soloUMIlen),

                paste0('--soloUMIfiltering ',param[['soloUMIfiltering']]),
                paste0('--soloCBmatchWLtype ',param[['soloCBmatchWLtype']]),
                paste0('--soloCellFilter ',param[['soloCellFilter']])
    )
    
    if(ezIsSpecified(param$cmdOptions)){
        cmd = paste(cmd, param$cmdOptions)
    }
    
    return(cmd)
}















