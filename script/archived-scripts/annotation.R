##' @describeIn ezFeatureAnnotation Gets the annotation from a .gtf file and transforms it into a tab-separated tabular .txt file.
writeAnnotationFromGtf = function(param, featureFile=param$ezRef@refFeatureFile, featAnnoFile=param$ezRef@refAnnotationFile){
  gtf = ezLoadFeatures(param, featureFile=featureFile)
  gtf = gtf[gtf$type == "exon", ]
  seqAnno = transcriptAnnoFromGtf(gtf, attributes=c("gene_name"))
  trProps = getTranscriptGcAndWidth(param)
  seqAnno$gc = trProps$gc[rownames(seqAnno)]
  ezWrite.table(seqAnno, file=featAnnoFile)
  return(seqAnno)
}