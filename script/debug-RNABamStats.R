
## posSpecErrorBam function
param = ScanBamParam(what=c("qname", "seq"),
                     flag=scanBamFlag(isUnmappedQuery=FALSE, isSecondaryAlignment=FALSE, isDuplicate=FALSE))
bamGA <- readGAlignments("/srv/gstore/projects/p2179/STAR_13239_2018-03-27--13-42-28/ILC3_D2_S16.bam", param=param)
genome <- readDNAStringSet("/srv/GT/reference/Homo_sapiens/Ensembl/GRCh38.p10/Sequence/WholeGenomeFasta/genome.fa")

errShort <- posSpecErrorBam(bamGA, genome)
