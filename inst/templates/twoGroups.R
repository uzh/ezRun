setwd("/srv/GT/analysis/gtan/p1997-edgeR/LR/EdgeR--AME--over--NS_Rmd")
library(ezRun)
param = list()
param[['cores']] = '1'
param[['ram']] = '2'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.PatchesSkipped/Annotation/Version-2015-07-05'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['testMethod']] = 'glm'
param[['grouping']] = 'Condition'
param[['sampleGroup']] = 'AME'
param[['refGroup']] = 'NS'
param[['normMethod']] = 'TMM'
param[['runGO']] = 'false'
param[['grouping2']] = ''
param[['backgroundExpression']] = '10'
param[['specialOptions']] = ''
param[['expressionName']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['comparison']] = 'AME--over--NS'
param[['name']] = 'AME--over--NS'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p1997/EdgeR_17729_AME--over--NS_2017-05-03--13-16-23'
param[['deTest']] = 'LR'
output = list()
output[['Name']] = 'AME--over--NS'
output[['Species']] = ''
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.PatchesSkipped/Annotation/Version-2015-07-05'
output[['Static Report [Link]']] = 'p1997/EdgeR_17729_AME--over--NS_2017-05-03--13-16-23/EdgeR--AME--over--NS/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-shiny.uzh.ch/fgcz_exploreDEG_app/?data=p1997/EdgeR_17729_AME--over--NS_2017-05-03--13-16-23/EdgeR--AME--over--NS/result-AME--over--NS-gbunszhvhmxz-EzResult.RData'
output[['Report [File]']] = 'p1997/EdgeR_17729_AME--over--NS_2017-05-03--13-16-23/EdgeR--AME--over--NS'
input = '/srv/gstore/projects/p1997/EdgeR_17729_AME--over--NS_2017-05-03--13-16-23/input_dataset.tsv'
param[['runGfold']] = 'false'
param[['runGO']] = 'true'
param[['doPrecomputeEnrichr']] = 'false'
EzAppEdger$new()$run(input=input, output=output, param=param)

# debug
output = EzDataset$new(meta=output, dataRoot=param$dataRoot)
load("/srv/GT/analysis/gtan/p1997-edgeR/LR/EdgeR--AME--over--NS_withSE_rightdeResult/result-AME--over--NS-gbunszhvhmxz-EzResult.RData")
deResult = EzResult(param=param, rawData=rawData, se=se)
param <- deResult$param
file.copy(from="/Users/gtan/Repos/FGCZ/ezRun/inst/templates/twoGroups.Rmd",
          to="twoGroups.Rmd", overwrite = TRUE)
rmarkdown::render(input="twoGroups.Rmd", envir = new.env(),
                  output_dir=".", output_file="00index.html")

# p2444
setwd("/srv/GT/analysis/gtan/p2444-edgeR")
library(ezRun)
param = list()
param[['cores']] = '1'
param[['ram']] = '2'
param[['scratch']] = '10'
param[['node']] = ''
param[['process_mode']] = 'DATASET'
param[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
param[['refFeatureFile']] = 'genes.gtf'
param[['featureLevel']] = 'gene'
param[['testMethod']] = 'glm'
param[['deTest']] = 'QL'
param[['grouping']] = 'Time'
param[['sampleGroup']] = 'd3'
param[['refGroup']] = 'BL'
param[['normMethod']] = 'TMM'
param[['runGO']] = 'true'
param[['grouping2']] = 'Subject'
param[['backgroundExpression']] = '10'
param[['transcriptTypes']] = 'protein_coding'
param[['specialOptions']] = ''
param[['expressionName']] = ''
param[['mail']] = 'ge.tan@fgcz.ethz.ch'
param[['comparison']] = 'd3--over--BL'
param[['name']] = 'd3--over--BL'
param[['dataRoot']] = '/srv/gstore/projects'
param[['resultDir']] = 'p2444/EdgeR_19220_healthy_pair_d3--over--BL_2017-07-03--22-42-08'
param[['runGfold']] = 'false'
param[['doPrecomputeEnrichr']] = 'false'
output = list()
output[['Name']] = 'd3--over--BL'
output[['Species']] = ''
output[['refBuild']] = 'Homo_sapiens/Ensembl/GRCh38.p10/Annotation/Release_89-2017-05-31'
output[['Static Report [Link]']] = 'p2444/EdgeR_19220_healthy_pair_d3--over--BL_2017-07-03--22-42-08/EdgeR--d3--over--BL/00index.html'
output[['Live Report [Link]']] = 'http://fgcz-shiny.uzh.ch/fgcz_exploreDEG_app/?data=p2444/EdgeR_19220_healthy_pair_d3--over--BL_2017-07-03--22-42-08/EdgeR--d3--over--BL/result-d3--over--BL-goeykvyjjsmb-EzResult.RData'
output[['Report [File]']] = 'p2444/EdgeR_19220_healthy_pair_d3--over--BL_2017-07-03--22-42-08/EdgeR--d3--over--BL'
input = '/srv/gstore/projects/p2444/EdgeR_19220_healthy_pair_d3--over--BL_2017-07-03--22-42-08/input_dataset.tsv'
EzAppEdger$new()$run(input=input, output=output, param=param)

## debug
output = EzDataset$new(meta=output, dataRoot=param$dataRoot)
load("/srv/gstore/projects/p2444/EdgeR_19220_healthy_pair_d3--over--BL_2017-07-03--22-42-08/EdgeR--d3--over--BL/result-d3--over--BL-goeykvyjjsmb-EzResult.RData")
deResult = EzResult(param=param, rawData=rawData, se=se)
setwd("/srv/GT/analysis/gtan/p2444-edgeR/EdgeR--d3--over--BL")
file.copy(from="/Users/gtan/Repos/FGCZ/ezRun/inst/templates/twoGroups.Rmd",
          to="twoGroups.Rmd", overwrite = TRUE)
rmarkdown::render(input="twoGroups.Rmd", envir = new.env(),
                  output_dir=".", output_file="00index.html")
