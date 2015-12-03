###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


## this works only if the final URL of the session xml is known
## it will fail if the directory is moved
##' @title Writes a .jnlp file
##' @description Writes a .jnlp file and adds the \code{projectId} and the \code{sessionUrl} to the template.
##' @param jnlpFile the connection to write the IGV. 
##' @param projectId a character representing the project ID.
##' @param sessionUrl a character representing the URL of the session.
##' @template roxygen-template
##' @return Returns lines written to the specified connection.
writeIgvJnlp = function(jnlpFile, projectId, sessionUrl){
  jnlpLines = readLines(ezIgvTemplateFile())
  jnlpLines = sub("PROJECT_ID", projectId, jnlpLines)
  jnlpLines = sub("IGV_SESSION", sessionUrl, jnlpLines)
  writeLines(jnlpLines, con=jnlpFile)  
}

##' @describeIn writeIgvJnlp Gets the IGV template file.
ezIgvTemplateFile = function(){
  if (IGV_TEMPLATE_FILE == ""){
    return(system.file("extdata/igvTemplate.jnlp", package="ezRun", mustWork = TRUE))
  } else {
    if (file.exists(IGV_TEMPLATE_FILE)){
      return(IGV_TEMPLATE_FILE)
    } else {
      stop("file does not exist: ", IGV_TEMPLATE_FILE)
    }
  }
}

## this adds a jnlp link to a static html; the html can be moved as long as the session xml stays in the same directory
##' @title Writes a jnlp link to a static html
##' @description Writes a jnlp link to a static html. The html can be moved as long as the session xml stays in the same directory.
##' @template htmlFile-template
##' @param projectId a character representing the project ID.
##' @param html a character representing an html address.
##' @template roxygen-template
writeJavaScriptIgvStarter = function(htmlFile, projectId, html){
  jnlpLines1 = paste('<jnlp spec="6.0+" codebase="http://www.broadinstitute.org/igv/projects/current">',
                    '<information>',
                    '<title>IGV 2.3</title>',
                    '<vendor>The Broad Institute</vendor>',
                    '<homepage href="http://www.broadinstitute.org/igv"/>',
                    '<description>IGV Software</description>',
                    '<description kind="short">IGV</description>',
                    '<icon href="IGV_64.png"/>',
                    '<icon kind="splash" href="IGV_64.png"/>',
                    '<offline-allowed/>',
                    '</information>',
                    '<security>',
                    '<all-permissions/>',
                    '</security>',
                    '<update check="background" />', #check="never" policy="never"/>',
                    '<resources>',
                    '<java version="1.6+" initial-heap-size="256m" max-heap-size="1100m" />',
                    '<jar href="igv.jar" download="eager" main="true"/>',
                    '<jar href="batik-codec__V1.7.jar" download="eager"/>',
                    '<jar href="goby-io-igv__V1.0.jar" download="eager"/>',  
                    '<property name="apple.laf.useScreenMenuBar" value="true"/>',
                    '<property name="com.apple.mrj.application.growbox.intrudes" value="false"/>',
                    '<property name="com.apple.mrj.application.live-resize" value="true"/>',
                    '<property name="com.apple.macos.smallTabs" value="true"/>',
                    '<property name="http.agent" value="IGV"/>',
                    '<property name="development" value="false"/>',
                    '</resources>',
                    '<application-desc  main-class="org.broad.igv.ui.Main">',
                    '<argument>--genomeServerURL=http://fgcz-gstore.uzh.ch/reference/igv_genomes.txt</argument>',
                    paste0('<argument>--dataServerURL=', "http://fgcz-gstore.uzh.ch/list_registries/", projectId, '</argument>'),
                    '<argument>',
                    sep="\n")
  jnlpLines2 = paste("</argument>",
                    '</application-desc>',
                    '</jnlp>',
                    sep="\n")
  ezWrite("<script>", con=html)
  ezWrite("function startIgvFromJnlp(label, locus){", con=html)
  ezWrite("var theSession = document.location.href.replace('", htmlFile, "', 'igvSession.xml');", con=html)
  ezWrite("var igvLink = 'data:application/x-java-jnlp-file;charset=utf-8,';", con=html)
  ezWrite("igvLink += '", RCurl::curlEscape(jnlpLines1), "';", con=html)
  ezWrite("igvLink += theSession;", con=html)
  ezWrite("igvLink += '", RCurl::curlEscape(jnlpLines2), "';", con=html)
  ezWrite("document.write(label.link(igvLink))", con=html)
  ezWrite("}", con=html)
  ezWrite("</script>", con=html)
}

##' @title Writes an IGV session
##' @description Writes an IGV session into a separate .xml file
##' @param genome a character representing the build name.
##' @param refBuild a character representing a file path to the reference build.
##' @param file a character representing the name of the file to write the IGV session in.
##' @param bamUrls a character vector containing BAM file links.
##' @param vcfUrls a character vector containing VCF file links.
##' @param locus a character describing which loci were used.
##' @template roxygen-template
##' @return Returns a character containing the name of the newly written .xml file.
##' @seealso \code{\link[XML]{newXMLNode}}
##' @seealso \code{\link[XML]{saveXML}}
##' @seealso \code{\link[XML]{addChildren}}
##' @examples
##' param = ezParam()
##' genome = getIgvGenome(param)
##' writeIgvSession(genome, param$ezRef["refBuild"])
writeIgvSession = function(genome, refBuild, file="igvSession.xml", bamUrls=NULL, vcfUrls=NULL,  locus="All"){

  requireNamespace("XML", warn.conflicts=WARN_CONFLICTS, quietly=!WARN_CONFLICTS)
  
  session=newXMLNode("Session", attrs =c(genome=genome, locus=locus, version="4"))
  resources=newXMLNode("Resources", parent=session)
  for (du in c(bamUrls, vcfUrls)){
    resources = addChildren(resources, newXMLNode("Resource", attrs=c(path=du)))
  }
  gtfFile=paste(REF_HOST, refBuild, "Genes", "genes.sorted.gtf", sep="/")
  resources = addChildren(resources, newXMLNode("Resource", attrs=c(path=gtfFile)))
  
  dataPanel = newXMLNode("Panel",attrs=c(name="DataPanel"), parent=session)
  for (du in bamUrls){
    track = newXMLNode("Track", attrs=c(altColor="0,0,178", autoscale="false", color="175,175,175",
                                        colorScale="ContinuousColorScale;0.0;156.0;255,255,255;175,175,175",
                                        displayMode="COLLAPSED", featureVisibilityWindow="-1" , fontSize="10",
                                        id=paste0(du, "_coverage"), name=paste(basename(du), "Coverage"),
                                        showDataRange="true", visible="true"))
    newXMLNode("DataRange",attrs=c(baseline="0.0", drawBaseline="true", flipAxis="false", maximum="1000.0" ,
                                   minimum="0.0", type="LOG"), parent=track)
    dataPanel = addChildren(dataPanel, track)
    track = newXMLNode("Track", attrs=c(altColor="0,0,178", color="0,0,178", colorOption="READ_STRAND",
                                        displayMode="EXPANDED", featureVisibilityWindow="100000", fontSize="10",
                                        id=du, name=basename(du), showDataRange="true", sortByTag="", visible="true"))
    dataPanel = addChildren(dataPanel, track)
  }
  ## VCF template
  #   <Panel height="915" name="Panel1417429920553" width="1681">
  #     <Track SQUISHED_ROW_HEIGHT="8" altColor="0,0,178" autoScale="false" clazz="org.broad.igv.track.FeatureTrack" color="0,0,178" colorMode="GENOTYPE" 
  #                 displayMode="EXPANDED" featureVisibilityWindow="100000000" fontSize="10" id="/Users/hubert/tmp/vcfFilt/all-20141117-haplo-filt.vcf.gz" name="all-20141117-haplo-filt.vcf.gz"
  #                 renderer="BASIC_FEATURE" sortable="false" visible="true" windowFunction="count"/>
  #     </Panel>
  for (du in vcfUrls){
    track = newXMLNode("Track", attrs=c(SQUISHED_ROW_HEIGHT="8", altColor="0,0,178", autoscale="false", color="0,0,178",
                                        clazz="org.broad.igv.track.FeatureTrack",
                                        colorMode="GENOTYPE",
                                        displayMode="COLLAPSED", featureVisibilityWindow="100000" , fontSize="10" , height="60",
                                        id=du, name=basename(du),
                                        renderer="BASIC_FEATURE", sortable="false", visible="true", windowFunction="count"))
    dataPanel = addChildren(dataPanel, track)  
  }
  
  
  featurePanel = newXMLNode("Panel",attrs=c(name="FeaturePanel"), parent=session)
  track = newXMLNode("Track", attrs=c(altColor="0,0,178", autoscale="false", color="0,0,178",
                                      colorScale="ContinuousColorScale;0.0;190.0;255,255,255;0,0,178",
                                      displayMode="COLLAPSED", featureVisibilityWindow="-1" , fontSize="10",
                                      id=gtfFile, name="Gene structure",
                                      showDataRange="true", visible="true"))
  featurePanel = addChildren(featurePanel, track)
  newXMLNode("PanelLayout", attrs=c(dividerFractions="0.7"), parent=session)  
  saveXML(session, file=file, prefix='<?xml version="1.0"? >', encoding="UTF-8", standalone="no")
  return(file)  
}

##' @describeIn writeIgvSession Writes an IGV session link.
writeIgvSessionLink = function(genome, refBuild, bamFiles, html, locus="All", label="Open Integrative Genomics Viewer", baseUrl=PROJECT_BASE_URL){
  #TODO stopifnot(grepl("^p", bamFiles))
  #bamFiles = sub("^/.*?/p", "p", bamFiles) ## TODO: replace only data root item
  urls = paste(baseUrl, bamFiles, sep="/")
  writeIgvSession(genome, refBuild, bamUrls=urls, locus=locus)
  #ezWrite("<p><a onClick='startIgv(\"", locus, "\")'>", label, "</a></p>", con=html)
  ezWrite('<p><script>startIgvFromJnlp("', label, '", "', locus, '"); </script></p>', con=html)
  return()
}

##' @describeIn writeIgvSession Adds an IGV session link to an object of the class bsdoc.
addIgvSessionLink = function(genome, refBuild, bamFiles, doc, locus="All", label="Open Integrative Genomics Viewer", baseUrl=PROJECT_BASE_URL){
  urls = paste(baseUrl, bamFiles, sep="/")
  writeIgvSession(genome, refBuild, bamUrls=urls, locus=locus)
  for (url in urls){
    addParagraph(doc, pot(url, hyperlink = url))
  }
  addJavascript(doc, text=paste0('startIgvFromJnlp("', label, '", "', locus, '")'))
}

##' @describeIn writeIgvSession Gets the IGV genome if specified or otherwise tries to get the build name from the parameters.
getIgvGenome = function(param){
  ifelse(ezIsSpecified(param$igvGenome),
         param$igvGenome,
         param$ezRef["refBuildName"])
}

## REFAC, but function is currently unused.
##' @title Gets an IGV locus link in .html format
##' @description Gets an IGV locus link in html format using chromosome, start and end information.
##' @param chrom the name of the chromosome.
##' @param start the start position of the locus.
##' @param end the end position of the locus.
##' @template roxygen-template
##' @return Returns a character containing a link in .html format.
##' @examples
##' getIgvLocusLink(17,31,465)
getIgvLocusLink = function(chrom, start, end){
  
  locus = paste0(chrom, ":", start, "-", end)
  link = paste0("http://www.broadinstitute.org/igv/projects/current/launch.php?locus=", locus)
  return(paste0("<a href=", link, ">", locus, "</a>"))
}


#http://www.broadinstitute.org/igv/projects/current/launch.php?locus=chr16:36762879-36778737
#http://www.broadinstitute.org/igv/projects/current/launch.php?sessionURL=/Users/hubert/work/IGV/igvLaunch-mm9.xml&locus=chr16:36762879-36778737

# writeIgvSessionScript = function(htmlFile, html){
#   ezWrite("<script>", con=html)
#   ezWrite("function startIgv(label, locus){", con=html)
#   ezWrite("var theSession = document.location.href.replace('", htmlFile, "', 'igvSession.xml');", con=html)
# 	ezWrite("var igvlink = 'http://www.broadinstitute.org/igv/projects/current/launch.php?sessionURL=' + theSession + '&locus=' + locus;", con=html)
#   ezWrite("document.write(label.link(igvlink))", con=html)
# 	ezWrite("}", con=html)
#   ezWrite("</script>", con=html)
# }
