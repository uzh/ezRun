# The parameter secondRef

The parameter let's you specify additional reference sequences beyond the reference for a given species. These could be

* trans-genes
* virus sequences
* second species

The parameter needs to point to a fasta file that is readable at runtime. Suggested locations are

* /srv/gstore/projects/p<current project>/<extraSequences or any other informative name>
* /srv/GT/databases/extra_references/<your informative folder name>

Preferred is the location on gstore because that provides persistent storage.

The interpretation of the fasta file is as follows

* If `file.exists(sub(".fa", ".gtf", param$secondRef))` then the assumption is that the fasta file contains DNA sequences (chromosomes or contigs) and the gtf file indicates the gene coordinates
* If there is no gtf file, then the assumption is that these additional sequences are used as is and in RNA applications are considered as spliced transcripts. Each transcript representing a different gene

## Behavior

### STAR

If there is a gtf a new temporary ref is built from the concatenation of the reference genome and the secondRef.

If there is now gtf, STAR adds the sequences just ad-hoc in the alignment step. No new reference is built.

In both cases the secondRef parameter also must be specified at subsequent `featureCounts`.


### CellRanger

Always a temporary ref is built.


### Kallisto

Not yeat supported.

