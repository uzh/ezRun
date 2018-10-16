package ch.uzh.fgcz.queue

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.extensions.picard._
import org.broadinstitute.gatk.queue.util.Logging
import org.broadinstitute.gatk.tools.walkers.haplotypecaller.ReferenceConfidenceMode
import org.broadinstitute.gatk.utils.ValidationExclusion

import scala.io.Source


class VariantCallingRna extends QScript with Logging {

	qs =>

	@Input(doc = "The reference file in FASTA format",
			shortName = "R", fullName = "reference", required = true)
	var pathRef: File = _

	@Input(doc = "File describing libraries. The file should contain a header and six columns: " +
			"Sample, Library, Platform, RGID, FileBAM",
			shortName = "libs", fullName = "pathLibs", required = true)
	var pathLibs: File = _

	@Input(doc = "Output VCF file",
			shortName = "pathOut", fullName = "pathOut", required = true)
	var pathOut: File = "07_variants_filt.vcf.gz"

	@Argument(doc="File with filters. It should contain a header and two columns: filter name and " +
			"filter expression",	shortName = "filters", fullName = "pathFilters", required = false)
	var pathFilters: Option[File] = None

	@Argument(doc="File with genotype filters. It should contain a header and two columns: filter " +
			"name and filter expression",	shortName = "gtfilters", fullName = "pathGtFilters",
			required = false)
	var pathGtFilters: Option[File] = None

	@Argument(doc="Heterozygosity level for joint genotyping.",
			shortName = "het", fullName = "heterozygosity", required = false)
	var heterozygosity: Double = 0.001

	@Argument(doc="Number of data threads to use for jobs that support -nt. This does not apply to " +
			"jobs that support scattering.", shortName = "t", fullName = "threadN", required = false)
	var threadN: Int = 8

	@Argument(doc = "Number of CPU threads to use for jobs that support -nct. Scatter jobs will get " +
			"ct / sj threads each", shortName = "ct", fullName = "cpuThreadN", required = false)
	var cpuThreadN: Int = 8

	@Argument(doc = "Number of scatter jobs to generate",
			shortName = "sj", fullName = "scatterN", required = false)
	var scatterN: Int = 8

	@Argument(doc = "Maximum amount of memory to use in GB",
			shortName = "maxMem", fullName = "maxMemory", required = false)
	var maxMem: Int = 32

	@Argument(doc = "If specified, soft clipped bases in the reads will not be analysed.",
			fullName = "dontUseSoftClippedBases", required = false)
	var dontUseSoftClippedBases: Boolean = false

	trait commonArgs extends CommandLineGATK {
		this.reference_sequence = qs.pathRef
	}


	def readLibMeta(path: File) : List[LibMeta] = {
		val bfs = Source.fromFile(path)
		var libMetaBuf = scala.collection.mutable.ListBuffer[LibMeta]()
		try {
			// Drop the header
			for (line <- bfs.getLines.drop(1)) {
				libMetaBuf += LibMeta.fromString(line)
			}
		} finally {
			bfs.close()
		}

		val fileErrMsg = "File '%s' in library '%s' does not exist or cannot be found."
		for (s <- libMetaBuf) {
			if (! s.fileBam.isFile) {
				throw new IllegalArgumentException(fileErrMsg.format(s.fileBam.toString, s.libraryId))
			}
		}
		libMetaBuf.toList
	}


	def readFilters(path: File) : Map[String, String] = {
		val bfs = Source.fromFile(path)
		var filtBuf = scala.collection.mutable.Map[String, String]()
		try {
			for (line <- bfs.getLines().drop(1)) {
				val arr = line.split("\t", -1)
				if (arr.length != 2) {
					val msg = s"Invalid number of columns in the filter file. Expected 2 but found " +
							s"${arr.length} in '$line'."
					throw new IllegalArgumentException(msg)
				}
				filtBuf += arr(0) -> arr(1)
			}
		}
		filtBuf.toMap
	}


	def script() {
		val libraries = readLibMeta(qs.pathLibs)
		val filters = qs.pathFilters.map {
			readFilters(_)
		}.getOrElse(Map[String, String]())
		val gtFilters = qs.pathGtFilters.map {
			readFilters(_)
		}.getOrElse(Map[String, String]())
		val adjCpuThreadN = Math.max(1, qs.cpuThreadN / qs.scatterN)
		val adjMemoryLimit = qs.maxMem / qs.threadN * adjCpuThreadN

		for (lib <- libraries) {
			// Add read groups
			val addRG = new AddOrReplaceReadGroups {
				this.analysisName = "addRG"
				this.input :+= lib.fileBam
				this.output = lib.fileBamRg
				this.RGID = lib.rgId
				this.RGLB = lib.libraryId
				this.RGPL = lib.platform
				this.RGPU = lib.libraryId
				this.RGSM = lib.sampleId
				this.isIntermediate = true
				this.javaMemoryLimit = adjMemoryLimit
			}
			add(addRG)

			val markDuplicates = new MarkDuplicates {
				this.isIntermediate = true
				this.input :+= lib.fileBamRg
				this.output = lib.fileDedup
				this.javaMemoryLimit = adjMemoryLimit
			}
			add(markDuplicates)

			val splitNCigarReads = new SplitNCigarReads with commonArgs {
				this.isIntermediate = true
				this.input_file :+= lib.fileDedup
				this.out = lib.fileSplit
				this.unsafe = ValidationExclusion.TYPE.ALLOW_N_CIGAR_READS
				this.javaMemoryLimit = adjMemoryLimit
			}
			add(splitNCigarReads)

			// The filter class is generated on the fly. Make sure you use full argument names rather
			// than public field names, e.g. reassign_mapping_quality_from and not
			// reassignMappingQualityFrom
			val mqFix = new PrintReads with ReassignOneMappingQuality with commonArgs {
				this.isIntermediate = false
				this.input_file :+= lib.fileSplit
				this.out = lib.fileMqFix
				this.reassign_mapping_quality_from = 255
				this.reassign_mapping_quality_to = 60
				this.javaMemoryLimit = adjMemoryLimit
			}
			add(mqFix)
		}

		// Base Recalibration 
		// https://software.broadinstitute.org/gatk/documentation/article?id=2801
		/*
		for (lib <- libraries) {
			// Analyse the patterns of covariation
			val baseRecalibrator = new BaseRecalibrator with commonArgs {
				this.isIntermediate = true
				this.input_file :+= lib.fileSplit
				this.out = lib.recalTable
				this.knownSites :+= qs.knownSites
				this.scatterCount = qs.scatterN
				this.num_cpu_threads_per_data_thread = adjCpuThreadN
				this.memoryLimit = adjMemoryLimit
			}

			// Analyse the covariance after recalibration
			val baseRecalibratorCov = new BaseRecalibrator with commonArgs {
				this.isIntermediate = true
				this.input_file :+= lib.fileSplit
				this.out = lib.postRecalTable
				this.knownSites :+= qs.knownSites
				this.scatterCount = qs.scatterN
				this.num_cpu_threads_per_data_thread = adjCpuThreadN
				this.memoryLimit = adjMemoryLimit
			}
			
			// Plot base qualities before and after recalibration
			val plotMaker = new AnalyzeCovariates with commonArgs {
				this.isIntermediate = false
				this.beforeReportFile = lib.recalTable
				this.afterReportFile = lib.postRecalTable
				this.plotReportFile = lib.recalPlots
			}

			// Apply recalibration to the data
			val readPrinter = new PrintReads with commonArgs {
				this.isIntermediate = false
				this.input_file = lib.fileSplit
				this.BQSR = lib.recalTable
				this.out = lib.fileRecal
				this.num_cpu_threads_per_data_thread = qs.threadN
			}

			add(baseRecalibrator, baseRecalibratorCov, plotMaker, readPrinter)
		}
		*/

		// Call variants per genotype
		val samples = libraries.groupBy(x => x.sampleId)
		val sampleVarFiles: Map[String, File] = (for {
			s <- samples.keys
		} yield {
			s -> new File(s"06_variants_$s.g.vcf.gz")
		}) (collection.breakOut)
		for (s <- samples) {
			val bams: Seq[File] = for (lib <- s._2) yield lib.fileMqFix
			val haplotypeCaller = new HaplotypeCaller with commonArgs {
				this.isIntermediate = false
				this.emitRefConfidence = ReferenceConfidenceMode.GVCF
				this.scatterCount = qs.scatterN
				this.num_cpu_threads_per_data_thread = adjCpuThreadN
				this.memoryLimit = adjMemoryLimit
				this.input_file ++= bams
				this.out = sampleVarFiles(s._1)
				this.dontUseSoftClippedBases = qs.dontUseSoftClippedBases
			}
			add(haplotypeCaller)
		}

		// Joint genotyping
		val genotypeGVCFs = new GenotypeGVCFs with commonArgs {
			this.isIntermediate = false
			this.memoryLimit = qs.maxMem
			this.num_threads = qs.threadN
			this.heterozygosity = qs.heterozygosity
			this.variant = sampleVarFiles.values.toSeq
			this.out = "07_variants.vcf.gz"
		}
		add(genotypeGVCFs)

		// Filter variants
		val filterNames = filters.keys.toSeq
		val filterExprs = for {k <- filterNames} yield filters(k)
		val gtFilterNames = gtFilters.keys.toSeq
		val gtFilterExprs = for {k <- gtFilterNames} yield gtFilters(k)
		val variantFilter = new VariantFiltration with commonArgs {
			this.isIntermediate = false
			// Sometimes it gets very hungry for memory, e.g. 32 GB may not be enough
			this.memoryLimit = qs.maxMem
			this.variant = genotypeGVCFs.out
			this.out = qs.pathOut
			this.filterName = filterNames
			this.filterExpression = filterExprs
			this.genotypeFilterName = gtFilterNames
			this.genotypeFilterExpression = gtFilterExprs
			// The following is specific to RNA-seq
			// https://software.broadinstitute.org/gatk/guide/article?id=3891
			this.clusterSize = 3
			this.clusterWindowSize = 35
		}
		add(variantFilter)

	}



	object LibMeta {
		def fromString(s: String): LibMeta = {
			assert(s.length > 0)
			// The negative number is required to make sure the empty trailing fields are returned
			val arr = s.split("\t", -1)
			val argN = 5
			if (arr.length != argN) {
				throw new IllegalArgumentException(s"Incorrect number of fields (${arr.length} != $argN) in '$s'")
			}
			val Array(
					sampleId,
					libraryId,
					platform,
					rgId,
					pathBam
			) = arr
			new LibMeta(
					sampleId,
					libraryId,
					platform,
					rgId,
					new File(pathBam)
			)
		}
	}

	class LibMeta(
			val sampleId: String,
			val libraryId: String,
			val platform: String,
			val rgId: String,
			val fileBam: File
	) {
		val fileBamRg = new File(List("01_withRG", libraryId).mkString("_") + ".bam")
		val fileDedup = new File(List("02_dedup", libraryId).mkString("_") + ".bam")
		val fileSplit = new File(List("03_split", libraryId).mkString("_") + ".bam")
		// Apparently, we can't add a filter to a walker easily
		val fileMqFix = new File(List("04_splitMqFix", libraryId).mkString("_") + ".bam")
		//val recalTable = new File(List("05_recalAnte", libraryId).mkString("_") + ".table")
		//val postRecalTable = new File(List("05_recalPost", libraryId).mkString("_") + ".table")
		//val recalPlots = new File(List("05_recalPlots", libraryId).mkString("_") + ".pdf")
		//val fileRecal = new File(List("06_recal", libraryId).mkString("_") + ".bam")
		val fileGvcf = new File(s"05_$sampleId.g.vcf.gz")

		override def toString = List(sampleId, libraryId).mkString("_")

		def mkReadGroup : String = {
			List(s"ID:$rgId", s"PL:$platform", s"LB:$libraryId", s"SM:$sampleId").mkString("\t")
		}
	}

}

