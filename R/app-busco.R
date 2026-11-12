###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

#busco -i /srv/gstore/projects/p25396/o30319_Spades_2023-02-06--14-53-26/DAMREC0027.fasta -o busco_test -m geno -l gammaproteobacteria_odb10 -c 4
#generate_plot.py -wd busco_test/
#/srv/GT/analysis/qiwei/p25396/o30319/o30319_Prodigal_2023-02-07--12-29-31

ezMethodBusco = function(
  input = NA,
  output = NA,
  param = NA,
  htmlFile = "00index.html"
) {
  opt = param$cmdOptions
  sampleName = input$getNames()
  draft = input$getFullPaths("Draft")
  cmd = paste(
    "busco -i ",
    draft,
    "-o",
    sampleName,
    "-l",
    param$lineage,
    "-m",
    param$mode,
    "-c",
    ezThreads(),
    opt,
    "1>",
    paste0(sampleName, "_busco.log")
  )
  ezSystem(cmd)
  cmd = paste("generate_plot.py -wd", sampleName)
  ezSystem(cmd)
  return("Success")
}


##' @template app-template
##' @templateVar method ezMethodProkka
##' @description Use this reference class to run
##' @seealso \code{\link{getPbmm2Reference}}
EzAppBusco <-
  setRefClass(
    "EzAppBusco",
    contains = "EzApp",
    methods = list(
      ## BUSCO: always runs (ezMethodBusco); Manni 2021 is the BUSCO v5 paper.
      ## OrthoDB v10: the lineage datasets of BUSCO 5 (a lineage without suffix gets _odb10);
      ## left out only when param$lineage names another OrthoDB release. OrthoDB v12 (Tegenfeldt
      ## 2025, the OrthoDB/BUSCO update that introduced the _odb12 sets) for an _odb12 lineage:
      ## BUSCO 5.8.0 accepts it and busco-data.ezlab.org serves _odb10 and _odb12 sets only (no
      ## _odb11 lineages exist, so none is cited).
      citation = function(param = list()) {
        c(
          "Manni, M., Berkeley, M.R., Seppey, M., Simão, F.A. & Zdobnov, E.M. BUSCO update: novel and streamlined workflows along with broader and deeper phylogenetic coverage for scoring of eukaryotic, prokaryotic, and viral genomes. Molecular Biology and Evolution 38(10), 4647-4654 (2021). https://doi.org/10.1093/molbev/msab199",
          if (!isTRUE(grepl("_odb(?!10)[0-9]+", param$lineage, perl = TRUE))) "Kriventseva, E.V. et al. OrthoDB v10: sampling the diversity of animal, plant, fungal, protist, bacterial and viral genomes for evolutionary and functional annotations of orthologs. Nucleic Acids Research 47(D1), D807-D811 (2019). https://doi.org/10.1093/nar/gky1053",
          if (isTRUE(grepl("_odb12", param$lineage))) "Tegenfeldt, F., Kuznetsov, D., Manni, M., Berkeley, M., Zdobnov, E.M. & Kriventseva, E.V. OrthoDB and BUSCO update: annotation of orthologs with wider sampling of genomes. Nucleic Acids Research 53(D1), D516-D522 (2025). https://doi.org/10.1093/nar/gkae987"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodBusco
        name <<- "EzAppBusco"
        appDefaults <<- rbind(
          kingdom = ezFrame(
            Type = "character",
            DefaultValue = "Bacteria",
            Description = "annotation mode and genetic code. Default is Bacteria"
          )
        )
      }
    )
  )
