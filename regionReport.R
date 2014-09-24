
## ----vignetteSetup, echo = FALSE, message = FALSE, warning = FALSE-------
## Track time spent on making the vignette
startTimeVignette <- Sys.time()

## Bib setup
library('knitcitations')

## Load knitcitations with a clean bibliography
cleanbib()
cite_options(hyperlink = 'to.doc', citation_format = 'text', style = 'html')
# Note links won't show for now due to the following issue
# https://github.com/cboettig/knitcitations/issues/63

## Write bibliography information
write.bibtex(c(knitcitations = citation('knitcitations'), 
    derfinder = citation('derfinder'), 
    regionReport = citation('regionReport'),
    knitrBootstrap = citation('knitrBootstrap'),
    ggbio = citation('ggbio'),
    ggplot2 = citation('ggplot2'),
    rCharts = citation('rCharts'), 
    knitr = citation('knitr')[3],
    rmarkdown = citation('rmarkdown')),
    file = 'derfinderPlotRef.bib')
bib <- read.bibtex('derfinderPlotRef.bib')

## Fix some names to work with CRAN and GitHub versions
names(bib)[names(bib) == 'hester2013knitrbootstrap'] <- 'hester2014knitrbootstrap'


## ----loadDerfinder, bootstrap.show.code=TRUE-----------------------------
## Load derfinder
library('derfinder')

## The output will be saved in the 'report' directory
dir.create('report', showWarnings = FALSE, recursive = TRUE)


## ----runDerfinderFake, eval=FALSE, bootstrap.show.code=TRUE--------------
## ## Save the current path
## initialPath <- getwd()
## setwd(file.path(initialPath, 'report'))
## 
## ## Generate output from derfinder
## 
## ## Collapse the coverage information
## collapsedFull <- collapseFullCoverage(list(genomeData$coverage),
## verbose=TRUE)
## 
## ## Calculate library size adjustments
## sampleDepths <- sampleDepth(collapsedFull, probs=c(0.5), nonzero=TRUE,
## verbose=TRUE)
## 
## ## Build the models
## group <- genomeInfo$pop
## adjustvars <- data.frame(genomeInfo$gender)
## models <- makeModels(sampleDepths, testvars=group, adjustvars=adjustvars)
## 
## ## Analyze chromosome 21
## analysis <- analyzeChr(chr='21', coverageInfo=genomeData, models=models,
## cutoffFstat=1, cutoffType='manual', seeds=20140330, groupInfo=group,
## mc.cores=1, writeOutput=TRUE, returnOutput=TRUE)
## 
## ## Save the stats options for later
## optionsStats <- analysis$optionsStats
## 
## ## Change the directory back to the original one
## setwd(initialPath)


## ----runDerfinderReal, bootstrap.show.code=TRUE--------------------------
## Copy previous results
file.copy(system.file(file.path('extdata', 'chr21'), package='derfinder', 
mustWork=TRUE), 'report', recursive=TRUE)


## ----mergeResults, bootstrap.show.code=TRUE, bootstrap.show.message=FALSE----
## Merge the results from the different chromosomes. In this case, there's 
## only one: chr21
mergeResults(chrs = 'chr21', prefix = 'report',
    genomicState = genomicState$fullGenome)


## ----loadLib, message=FALSE, bootstrap.show.code=TRUE--------------------
## Load derfindeReport
library('regionReport')


## ----createReportFake, eval=FALSE, bootstrap.show.code=TRUE--------------
## ## Generate the HTML report
## report <- derfinderReport(prefix='report', browse=FALSE,
##     nBestRegions=15, makeBestClusters=TRUE, outdir='html',
##     fullCov=list('21'=genomeDataRaw$coverage), optionsStats=optionsStats)


## ----createReportReal, echo=FALSE, message=FALSE, bootstrap.show.code=FALSE----
## Generate the HTML report in a clean environment
library('devtools')
cat("## Generate the report in an isolated environment
## This helps avoids conflicts with generating the vignette
library(derfinder)
library(regionReport)

## Load optionsStats
load(file.path('report', 'chr21', 'optionsStats.Rdata'))

## Create report
report <- derfinderReport(prefix='report', browse=FALSE,
    nBestRegions=15, makeBestClusters=TRUE, outdir='html',
    fullCov=list('21'=genomeDataRaw$coverage), optionsStats=optionsStats)
    
## Clean up
file.remove('derfinderReport-isolated.R')
", file='derfinderReport-isolated.R')
clean_source('derfinderReport-isolated.R', quiet=TRUE)


## ----vignetteBrowse, eval=FALSE, bootstrap.show.code=TRUE----------------
## ## Browse the report
## browseURL(report)


## ----openIncludedReport, eval=FALSE--------------------------------------
## browseURL(system.file(file.path('basicExploration', 'basicExploration.html'),
##     package = 'regionReport', mustWork = TRUE))


## ----'advancedArg'-------------------------------------------------------
## URLs to advanced arguemtns
derfinder::advancedArg('derfinderReport', package = 'regionReport',
    browse = FALSE)
## Set browse = TRUE if you want to open them in your browser


## ----createVignette, eval=FALSE, bootstrap.show.code=TRUE----------------
## ## Create the vignette
## library('knitrBootstrap')
## 
## knitrBootstrapFlag <- packageVersion('knitrBootstrap') < '1.0.0'
## if(knitrBootstrapFlag) {
##     ## CRAN version
##     system.time(knit_bootstrap('regionReport.Rmd', chooser=c('boot', 'code'), show_code = TRUE))
##     unlink('regionReport.md')
## } else {
##     ## GitHub version
##     library('rmarkdown')
##     system.time(render('regionReport.Rmd',
##         'knitrBootstrap::bootstrap_document'))
## }
## ## Note: if you prefer the knitr version use:
## # library('rmarkdown')
## # system.time(render('regionReport.Rmd', 'html_document'))
## 
## ## Extract the R code
## library('knitr')
## knit('regionReport.Rmd', tangle = TRUE)
## 
## ## Copy report output to be distributed with the package for comparison
## ## purposes
## file.copy(file.path('report', 'html', 'basicExploration.html'),
##     '../inst/basicExploration/basicExploration.html', overwrite=TRUE)
## 
## ## Clean up
## file.remove('derfinderPlotRef.bib')
## #unlink('regionReport_files', recursive=TRUE)
## unlink('report', recursive = TRUE)


## ----vignetteReproducibility1, echo=FALSE--------------------------------
## Date the report was generated
Sys.time()


## ----vignetteReproducibility2, echo=FALSE, bootstrap.show.code=FALSE-----
## Processing time in seconds
totalTimeVignette <- diff(c(startTimeVignette, Sys.time()))
round(totalTimeVignette, digits=3)


## ----vignetteReproducibility3, echo=FALSE--------------------------------
## Session info
session_info()


## ----vignetteBiblio, results='asis', echo=FALSE--------------------------
## Print bibliography
bibliography()

