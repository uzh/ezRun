.mergeBams =function(bams, output="merged.bam", prependSetId=FALSE){
  
  if (prependSetId){
    modBams = character()
    for (i in 1:length(bams)){
      modBam = basename(bams[i])
      stopifnot(!file.exists(modBam))
      modSam = sub(".bam$", ".sam", modBam)
      cmd = paste("samtools", "view -H", bams[i], ">", modSam)
      ezSystem(cmd)
      cmd = paste0("samtools", " view ", bams[i], " | awk '{print \"set_", i, "_\" $0}'", " >> ", modSam)
      ezSystem(cmd)
      cmd = paste("samtools", "view -S -b -h", modSam, ">", modBam)
      ezSystem(cmd)
      cmd = paste("samtools", "index", modBam)
      ezSystem(cmd)
      modBams[i] = modBam
    }
    cmd = paste("samtools", "merge -r", output, paste(modBams, collapse=" "))
    ezSystem(cmd)
    file.remove(modBams)
    file.remove(paste0(modBams, ".bai"))
  } else {
    cmd = paste("samtools", "merge -r", output, paste(bams, collapse=" "))
    ezSystem(cmd)
  }
  cmd = paste("samtools", "index", output)
  ezSystem(cmd)
  return(output)
}

# .samToBam = function(sam, bam, samtools ="samtools", maxMem="1000000000", removeSam=TRUE, sortIndexBam=TRUE){
#   tmpBam = sub("sam$", "tmp.bam", sam)
#   cmd = paste(samtools, "view -S", sam, "-b -o", tmpBam)
#   ezSystem(cmd)
#   if (removeSam){
#     file.remove(sam)
#   }
#   if (sortIndex){
#     ezSortIndexBam(tmpBam, bam, samtools=samtoosl, maxmem=maxMem, removeBam=TRUE)
#   }
# }