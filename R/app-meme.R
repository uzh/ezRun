###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

ezMethodMEME = function(input = NA, output = NA, param = NA) {
  Sys.setenv(
    PERL5LIB = paste(
      "/usr/local/ngseq/packages/Dev/Perl/5.40.0/lib/perl5/site_perl/5.40.0",
      Sys.getenv("PERL5LIB"),
      sep = ":"
    )
  )
  opt = param$cmdOptions
  db = param$motifDB
  sampleName = input$getNames()
  if (param$filterPeaks) {
    param$relPos = c('upstream', 'overlapStart', 'inside', 'includeFeature')
    peaks = readxl::read_xlsx(input$getFullPaths("CalledPeaks"))
    peaks = peaks[abs(peaks$distancetoFeature) <= param$distance, ]
    peaks = peaks[peaks$insideFeature %in% param$relPos, ]
    peaks = peaks[abs(peaks$fold_enrichment) >= param$minFold, ]
    toRemove = which(duplicated(peaks$name))

    if (length(toRemove) > 0) {
      peaks = peaks[-toRemove, ]
    }
    ezWrite.table(
      peaks,
      basename(input$getFullPaths("CalledPeaks")),
      row.names = F
    )

    #CreateBedFile:
    bed = ezRead.table(input$getFullPaths("BED"), row.names = NULL, header = F)
    bed = bed[bed$V4 %in% peaks$name, ]
    bed = bed[order(bed$V4), ]
    peaks = peaks[order(peaks$name), ]
    bed$V4 = paste0(
      bed$V4,
      '__',
      paste(
        peaks$chr,
        peaks$start,
        peaks$end,
        peaks$gene_name,
        peaks$distancetoFeature,
        sep = '_'
      )
    )
    bedFileName = basename(input$getFullPaths("BED"))
    ezWrite.table(bed, bedFileName, row.names = F, col.names = F)

    refFile = param$ezRef["refFastaFile"]
    peakSeqFile = paste0(bedFileName, "_peaks.fa")
    cmd = paste(
      "bedtools",
      " getfasta -fi",
      refFile,
      " -bed ",
      bedFileName,
      " -name -fo ",
      peakSeqFile
    )
    ezSystem(cmd)
    cmd = paste(
      "meme-chip -oc",
      sampleName,
      "-time 300 -order 1",
      db,
      opt,
      peakSeqFile
    )
  } else {
    cmd = paste(
      "meme-chip -oc",
      sampleName,
      "-time 300 -order 1",
      db,
      opt,
      input$getFullPaths("PeakSequences")
    )
  }

  ezSystem(cmd)
  setwd(sampleName)
  ezSystem(paste('mv meme-chip.html ', paste0(sampleName, "_meme-chip.html")))
  setwd('..')
  return("Success")
}

##' @author Opitz, Lennart
##' @template app-template
##' @templateVar method ezMethodMEME(input=NA, output=NA, param=NA)
##' @description Use this reference class to run
EzAppMEME <-
  setRefClass(
    "EzAppMEME",
    contains = "EzApp",
    methods = list(
      ## meme-chip (ezMethodMEME) with motifDB and cmdOptions pasted on its command line:
      ## - MEME-ChIP, CentriMo and Tomtom: always.
      ## - MEME / STREME: de novo discovery unless cmdOptions sets -meme-nmotifs 0 / -streme-nmotifs 0.
      ## - SpaMo / FIMO: unless cmdOptions sets -spamo-skip / -fimo-skip.
      ## - JASPAR, HOCOMOCO and UniPROBE releases: when their file is given with -db in motifDB.
      citation = function(param = list()) {
        opt <- if (is.null(param$cmdOptions)) "" else param$cmdOptions
        db <- if (is.null(param$motifDB)) "" else param$motifDB
        c(
          "Machanick, P. & Bailey, T.L. MEME-ChIP: motif analysis of large DNA datasets. Bioinformatics 27(12), 1696-1697 (2011). https://doi.org/10.1093/bioinformatics/btr189",
          if (!grepl("-meme-nmotifs\\s+0(\\s|$)", opt)) "Bailey, T.L. & Elkan, C. Fitting a mixture model by expectation maximization to discover motifs in biopolymers. Proceedings of the Second International Conference on Intelligent Systems for Molecular Biology 2, 28-36 (1994).",
          if (!grepl("-streme-nmotifs\\s+0(\\s|$)", opt)) "Bailey, T.L. STREME: accurate and versatile sequence motif discovery. Bioinformatics 37(18), 2834-2840 (2021). https://doi.org/10.1093/bioinformatics/btab203",
          "Bailey, T.L. & Machanick, P. Inferring direct DNA binding from ChIP-seq. Nucleic Acids Research 40(17), e128 (2012). https://doi.org/10.1093/nar/gks433",
          "Gupta, S., Stamatoyannopoulos, J.A., Bailey, T.L. & Noble, W.S. Quantifying similarity between motifs. Genome Biology 8, R24 (2007). https://doi.org/10.1186/gb-2007-8-2-r24",
          if (!grepl("-spamo-skip", opt)) "Whitington, T., Frith, M.C., Johnson, J. & Bailey, T.L. Inferring transcription factor complexes from ChIP-seq data. Nucleic Acids Research 39(15), e98 (2011). https://doi.org/10.1093/nar/gkr341",
          if (!grepl("-fimo-skip", opt)) "Grant, C.E., Bailey, T.L. & Noble, W.S. FIMO: scanning for occurrences of a given motif. Bioinformatics 27(7), 1017-1018 (2011). https://doi.org/10.1093/bioinformatics/btr064",
          if (grepl("JASPAR[^ ]*2014", db)) "Mathelier, A. et al. JASPAR 2014: an extensively expanded and updated open-access database of transcription factor binding profiles. Nucleic Acids Research 42(D1), D142-D147 (2014). https://doi.org/10.1093/nar/gkt997",
          if (grepl("JASPAR[^ ]*2016", db)) "Mathelier, A. et al. JASPAR 2016: a major expansion and update of the open-access database of transcription factor binding profiles. Nucleic Acids Research 44(D1), D110-D115 (2016). https://doi.org/10.1093/nar/gkv1176",
          if (grepl("JASPAR2018", db)) "Khan, A. et al. JASPAR 2018: update of the open-access database of transcription factor binding profiles and its web framework. Nucleic Acids Research 46(D1), D260-D266 (2018). https://doi.org/10.1093/nar/gkx1126",
          if (grepl("JASPAR2020", db)) "Fornes, O. et al. JASPAR 2020: update of the open-access database of transcription factor binding profiles. Nucleic Acids Research 48(D1), D87-D92 (2020). https://doi.org/10.1093/nar/gkz1001",
          if (grepl("JASPAR2022", db)) "Castro-Mondragon, J.A. et al. JASPAR 2022: the 9th release of the open-access database of transcription factor binding profiles. Nucleic Acids Research 50(D1), D165-D173 (2022). https://doi.org/10.1093/nar/gkab1113",
          if (grepl("HOCOMOCOv9", db)) "Kulakovskiy, I.V. et al. HOCOMOCO: a comprehensive collection of human transcription factor binding sites models. Nucleic Acids Research 41(D1), D195-D202 (2013). https://doi.org/10.1093/nar/gks1089",
          if (grepl("HOCOMOCOv10", db)) "Kulakovskiy, I.V. et al. HOCOMOCO: expansion and enhancement of the collection of transcription factor binding sites models. Nucleic Acids Research 44(D1), D116-D125 (2016). https://doi.org/10.1093/nar/gkv1249",
          if (grepl("HOCOMOCOv11", db)) "Kulakovskiy, I.V. et al. HOCOMOCO: towards a complete collection of transcription factor binding models for human and mouse via large-scale ChIP-Seq analysis. Nucleic Acids Research 46(D1), D252-D259 (2018). https://doi.org/10.1093/nar/gkx1106",
          if (grepl("H12(CORE|INVITRO|INVIVO|RSNP)", db)) "Vorontsov, I.E. et al. HOCOMOCO in 2024: a rebuild of the curated collection of binding models for human and mouse transcription factors. Nucleic Acids Research 52(D1), D154-D163 (2024). https://doi.org/10.1093/nar/gkad1077",
          if (grepl("uniprobe", db, ignore.case = TRUE)) "Newburger, D.E. & Bulyk, M.L. UniPROBE: an online database of protein binding microarray data on protein-DNA interactions. Nucleic Acids Research 37, D77-D82 (2009). https://doi.org/10.1093/nar/gkn660"
        )
      },
      initialize = function() {
        "Initializes the application using its specific defaults."
        runMethod <<- ezMethodMEME
        name <<- "EzAppMEME"
      }
    )
  )
