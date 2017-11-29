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
  setwdNew(basename(output$getColumn("Static Report")))
  dataset = input$meta

### analyze results with phyloseq
### load OTUs, taxa and sample file 
otuFileName = "testMothur.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.shared"
taxaFileName = "testMothur.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.unique_list.0.03.cons.taxonomy"
sampleFileName = "mouse.time.design"
### OTUs
otuObjectPacBio <- phyloSeqOTU(input$OTU_pacbio)
otuObjectIllumina <- phyloSeqOTU(input$OTU_Illumina)
### taxonomy
taxaObjectPacBio <- phyloSeqTaxa(input$Taxonomy_pacbio)
taxaObjectIllumina <- phyloSeqTaxa(input$Taxonomy_Illumina)
### Samples
designMatrix <- param$designMatrix 
sampleObjectIllumina <- designMatrix[designMatrix$`Technology [Factor]` == "Illumina",]
sampleObjectIllumina <- phyloSeqSample(sampleObjectIllumina)
### create plots: 1. abundance
physeqIll = phyloseq(otuObjectIllumina, taxaObjectIllumina,sampleObjectIllumina)
physeqPacBio = phyloseq(otuObjectPacBio, taxaObjectPacBio)



plot_bar(physeq, fill = "Phylum")

### create plots: 2. tree
random_tree = rtree(ntaxa(physeq), rooted=TRUE, tip.label=taxa_names(physeq))
plot(random_tree)

### merge data
physeq2 = phyloseq(otuObject, taxaObject, sampleObject, random_tree)
physeqBacteroidetes <- subset_taxa(physeq2, Phylum == "Bacteroidetes")
plot_bar(physeqBacteroidetes, x="time", fill = "Genus")
plot_bar(physeqBacteroidetes, fill = "Genus", facet_grid=~time)
### create plots: 3. tree with sample info
plot_tree(physeq2, color="time", label.tips="taxa_names", ladderize="left", plot.margin=0.3)

### create plots: 3. heatmap
plot_heatmap(physeq2, taxa.label="Phylum")

  ## Copy the style files and templates
  styleFiles <- file.path(system.file("templates", package="ezRun"),
                          c("fgcz.css", "FastQC.Rmd", "FastQC_overview.Rmd",
                            "fgcz_header.html", "banner.png"))
  file.copy(from=styleFiles, to=".", overwrite=TRUE)
  
  plots = c("Summary of sequences in each sample.png",
            "Per sequence quality scores"="per_sequence_quality.png",
            "Per tile sequence quality"="per_tile_quality.png",
            "Per base sequence content"="per_base_sequence_content.png",
            "Per sequence GC content"="per_sequence_gc_content.png",
            "Per base N content"="per_base_n_content.png",
            "Sequence Length Distribution"="sequence_length_distribution.png",
            "Sequence Duplication Levels"="duplication_levels.png",
            "Adapter Content"="adapter_content.png",
            "Kmer Content"="kmer_profiles.png")
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
                  appDefaults <<- rbind(cutOff = ezFrame(Type="numeric",  DefaultValue="0,03",Description="Cut-off for OTU clustering.")
                  )
                }
              )
  )
