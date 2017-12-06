###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


ezMethodMothurPhyloSeqAnalysis = function(input=NA, output=NA, param=NA, 
                          htmlFile="00index.html"){
  require(rmarkdown)
  require(ShortRead)
  require(phyloseq)
  require(plyr)
  require(ape)
  require(ggplot2)
  require(DESeq2)
  dataset = input$meta

### Analyzes results with phyloseq: preparing objects to be processed in the Rmd file

### OTUs
copyOTUTablePacbioCmd <- paste("cp", input$getFullPaths("OTU_pacbio"),"./", sep = " ")
ezSystem(copyOTUTablePacbioCmd)
otuObjectPacBio <- phyloSeqOTU(basename(input$getFullPaths("OTU_pacbio")))

copyOTUTableIllCmd <- paste("cp", input$getFullPaths("OTU_Illumina"),"./", sep = " ")
ezSystem(copyOTUTableIllCmd)
otuObjectIllumina <- phyloSeqOTU(basename(input$getFullPaths("OTU_Illumina")))


### taxonomy
copytaxaTablePacbioCmd <- paste("cp", input$getFullPaths("Taxonomy_pacbio"),"./", sep = " ")
ezSystem(copytaxaTablePacbioCmd)
taxaObjectPacBio <- phyloSeqTaxa(basename(input$getFullPaths("Taxonomy_pacbio")))
copytaxaTableIllCmd <- paste("cp", input$getFullPaths("Taxonomy_Illumina"),"./", sep = " ")
ezSystem(copytaxaTableIllCmd)
taxaObjectIllumina <- phyloSeqTaxa(basename(input$getFullPaths("Taxonomy_Illumina")))
### pruning level
prunePB <- param$pacbioRepresentativeOTUs
pruneIll <- param$illuminaRepresentativeOTUs
### Samples
designMatrix <- ezRead.table("/srv/GT/analysis/course_sushi/public/projects/p2000/MetagenomicsCourseTestData/designMatrix.tsv")
sampleObjectIllumina <- designMatrix[designMatrix$`Technology [Factor]` == "Illumina",]
sampleObjectIllumina <- phyloSeqSample(sampleObjectIllumina)
sampleObjectPacbio <- designMatrix[designMatrix$`Technology [Factor]` == "PacBio",]
sampleObjectPacbio <- phyloSeqSample(sampleObjectPacbio)

### create, add trees, preprocess and prune phyloseq object 
physeqIllNoTree = phyloseq(otuObjectIllumina, taxaObjectIllumina,sampleObjectIllumina)
treeObjectIll = rtree(ntaxa(physeqIllNoTree), rooted=TRUE, tip.label=taxa_names(physeqIllNoTree))
physeqIll <- merge_phyloseq(physeqIllNoTree,treeObjectIll)
physeqIll <- phyloSeqPreprocess(physeqIll)
physeqIll <- prune_taxa(taxa_names(physeqIll)[1:pruneIll], physeqIll)

physeqPacBioNoTree = phyloseq(otuObjectPacBio, taxaObjectPacBio, sampleObjectPacbio)
treeObjectPacbio = rtree(ntaxa(physeqPacBioNoTree), rooted=TRUE, tip.label=taxa_names(physeqPacBioNoTree))
physeqPacBio <-  merge_phyloseq(physeqPacBioNoTree,treeObjectPacbio)
physeqPacBio <- phyloSeqPreprocess(physeqPacBio)
physeqPacBio <- prune_taxa(taxa_names(physeqPacBio)[1:prunePB], physeqPacBio)

physeqMergedNoTree <- merge_phyloseq(physeqIllNoTree,physeqPacBioNoTree)
treeMerged <- rtree(ntaxa(physeqMergedNoTree), rooted=TRUE, tip.label=taxa_names(physeqMergedNoTree))
physeqMerged <- merge_phyloseq(physeqMergedNoTree,treeMerged)
physeqMergedPruned <- prune_taxa(taxa_names(physeqMerged)[1:50], physeqMerged)

  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "MothurPhyloseqAnalysisMerged.Rmd", 
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  rmarkdown::render(input="MothurPhyloseqAnalysisMerged.Rmd", envir = new.env(),
                    output_dir=".", output_file=htmlFile, quiet=TRUE,fig.width=18, fig.height=12)
}
##' @template app-template
##' @templateVar method ezMethodMothurPhyloSeqAnalysis()
##' @templateVar htmlArg )
##' @description Use this reference class to run 
EzAppMothurPhyloSeqAnalysis <-
  setRefClass("EzAppMothurPhyloSeqAnalysis",
              contains = "EzApp",
              methods = list(
                initialize = function()
                {
                  "Initializes the application using its specific defaults."
                  runMethod <<- ezMethodMothurPhyloSeqAnalysis
                  name <<- "EzAppMothurPhyloSeqAnalysis"
                  appDefaults <<- rbind(pacbioRepresentativeOTUs = ezFrame(Type="numeric",  DefaultValue="",Description="Number of core OTUs for Pacbio samples."),
                                        illuminaRepresentativeOTUs = ezFrame(Type="character",  DefaultValue="",Description="Number of core OTUs for Illumina samples.")
                  )
                }
              )
  )
