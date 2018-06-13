## Single cell
library(Rsubread)
setwd("/scratch/gtan/debug/p2488-scRNA/HiSeq4000_20180209_RUN424_o3954")
localBamFile <- "/srv/gstore/projects/p2488/SCCountsApp_24663_2018-06-10--18-09-52/20180209.B-single_cell_hippo_neurons_E1_1_unmapped.bam"

localBamFile <- "20180209.B-single_cell_hippo_neurons_E1_1_unmapped.bam"
gtfFile <- "/srv/GT/reference/Mus_musculus/Ensembl/GRCm38.p5/Annotation/Release_91-2018-02-26/Genes/genes.gtf"

countResult = featureCounts(localBamFile, annot.inbuilt=NULL,
                            annot.ext=gtfFile, isGTFAnnotationFile=TRUE,
                            GTF.featureType="exon",
                            GTF.attrType="gene_id",
                            strandSpecific=0,
                            nthreads=4,
                            byReadGroup=TRUE)
require(Rsamtools)
bamHeaders <- scanBamHeader(localBamFile)
