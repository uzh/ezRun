###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

# notes for implementation

## STARsolo run:
"
/path/to/STAR --genomeDir /path/to/genome/dir/ \
              --readFilesIn cDNAfragmentSequence.fastq.gz CellBarcodeUMIsequence.fastq.gz \
              [...other parameters...] \
              --soloType ...  \
              --soloCBwhitelist /path/to/cell/barcode/whitelist
"

### extra parameters:

# barcode lengths: 
#   --soloUMIlen 12
# reads:
#   Importantly, in the --readFilesIn option, the 1st file has to be cDNA read, and the 2nd file has to be the barcode (cell+UMI) read
# For multiple lanes, use commas separated lists for Read2 and Read1:
#   --readFilesIn Read2_Lane1.fastq.gz,Read2_Lane2.fastq.gz,Read2_Lane3.fastq.gz  Read1_Lane1.fastq.gz,Read1_Lane2.fastq.gz,Read1_Lane3.fastq.gz
# barcode geometry:
#   --soloCBstart, --soloCBlen, --soloUMIstart, --soloUMIlen
# type of SCRNA-seq:
#   string(s): type of single-cell RNA-seq
#     CB_UMI_Simple   ... (a.k.a. Droplet) one UMI and one Cell Barcode of fixed length in read2, e.g. Drop-seq and 10X Chromium
#     CB_UMI_Complex  ... one UMI of fixed length, but multiple Cell Barcodes of varying length, as well as adapters sequences are allowed in read2 only, e.g. inDrop.



ezMethodSTARsolo = function(input=NA, output=NA, param=NA){
    sampleName = input$getNames()
    sampleDirs = strsplit(input$getColumn("RawDataDir"), ",")[[sampleName]]
    sampleDirs <- file.path(input$dataRoot, sampleDirs)
    sampleDir <- paste(sampleDirs, collapse=",")
    STARsoloFolder = paste0(sampleName, "-STARsolo")
    
    cmd = paste(STARsoloPath,
                # STARsolo parameters
                paste0('--genomeDir ',refDir),
                paste0('--readFilesIn ',cDNAfragmentSequence.fastq.gz,' ',CellBarcodeUMIsequence.fastq.gz),
                paste0('--soloType',)
                
                # general parameters
                paste0('--')  
      
    )



}
















