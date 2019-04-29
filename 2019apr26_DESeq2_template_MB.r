################################################################################
### R scripts BarSeq
###
### April 2019
### designed to be executed with SARTools 1.5.0
################################################################################

################################################################################
###                parameters: to be modified by the user                    ###
################################################################################
################################################################################

rm(list = ls())                                          # remove all the objects from the R session

workDir <-
  getwd()                                   # working directory for the R session

projectName <-
  "2019April27_HOP_mar29_MB"              # name of the project

author <-
  "gg"                                         # author of the statistical analysis/report

featuresToRemove <-
  NULL                               # names of the features to be removed, NULL if no feature to remove

cpmCutoff <- 1

varInt <-
  "drug"                                       # factor of interest

condRef <-
  "ctrl"                                      # reference biological condition

batch <-
  NULL                                          # blocking factor: NULL (default) or "batch" for example

fitType <-
  "local"                                     # mean-variance relationship: "parametric" (default),"local" or "mean"

cooksCutoff <-
  TRUE                                    # TRUE/FALSE to perform the outliers detection (default is TRUE)

independentFiltering <-
  TRUE                           # TRUE/FALSE to perform independent filtering (default is TRUE)

alpha <-
  0.05                                          # threshold of statistical significance

pAdjustMethod <-
  "BH"                                  # p-value adjustment method: "BH" (default) or "BY"

typeTrans <-
  "VST"                                     # transformation for PCA/clustering: "VST" or "rlog"

locfunc <-
  "median"                                  # "median" (default) or "shorth" to estimate the size factors


# #
workDir <- getwd()                                           # working directory for the R session


projectName <- "2019Apr26_barseqMJ"                # name of the project
author <- "gg"                                               # author of the statistical analysis/report


forceCairoGraph <- FALSE

source('~/RProjects/2019_April_barseq_pipeline/april25_2019_marjan_barseq_funcs.R')
source('~/RProjects/2019_April_barseq_pipeline/2019_apr25_source_functions.R')

# vector of colors of each biological condition on the plots

geneDir =  "/home/common/barseq_gg_pdata_exp_info/summary_barseq"

fdata = read.delim(
  paste(geneDir, "2019_apr21_fdata.txt", sep = '/'),
  stringsAsFactors = F,
  check.names = F
)

pool = "noness"

mycolors = c(
  "darkorange",
  "dodgerblue",
  "limegreen",
  "navy",
  "mediumpurple"  ,
  "royalblue3",
  "darkolivegreen4",
  "firebrick",
  "cyan",
  "hotpink",
  "plum4",
  "blue",
  "magenta2",
  "skyblue",
  "green",
  "red",
  "steelblue",
  "tomato",
  "purple",
  "yellow"
)

palette(mycolors)
################################################################################
#if using my runDESeq2 = runDESeq3, these parameters need to be set
#

Shrink = T             # yes or no to shrink after running DESeq

lfcThreshold = 1       # log2 fold change default for DESeq results if Shrink = F, or for                                lfcShrink if Shrink = T

type = "apeglm"        # shrinking method for lfcShrink, apeglm is now (2019 April) recommended


# Defaults are same as above, to return the equivant of runDESeq2, Shrink = F, lfcThreshold = 0
# with the exception that not all combinatorial pairwise drug combinations will be run
################################################################################
###                       LOADING DATA RUN SPECIFIC DATA                     ###
################################################################################
### get count matrix and associated target info
###

rawDir <-
  "/home/common/barseq_output"                 # path to the directory containing raw counts files

# specific count matrix of interest
countMat =  "190329_NB501711_0151_AHT77HBGX9-BS19-RGTBS10-56099-adapter-notrim-results-2019-04-10-26305-counts.rds"

# newpool
# "190219_NB501711_0141_AHWFM2BGX9-BS18-BS18-18060-adapter-mask-results-2019-02-20-20860-counts.rds"

# correct dox count mx
# "190329_NB501711_0151_AHT77HBGX9-BS19-RGTBS10-56099-adapter-notrim-results-2019-04-10-26305-counts.rds"
# had to be fixed as the original names with the screwups were duplicated
#190329_NB501711_0151_AHT77HBGX9-BS19-RGTBS10-56099-adapter-notrim-results-2019-04-10-26305-HOPfixedcounts.rds"

# gilead
#"190417_NB501711_0153_AHT55CBGX9-bs21-RGTBS21-155163-adapter-notrim-results-2019-04-18-4193-counts.rds"

countFiles = paste(rawDir, countMat, sep = "/")

expDir   = "/home/common/barseq_gg_pdata_exp_info"

# specific target file for count matrix of interest
expFiles = "2019april20_mar29_dox_ndea_pdata.txt"

targetFile = paste(expDir, expFiles, sep = "/")

library(SARTools)

# loading target file
#
#### this a screwy function; i prefer just to read the file in with read.delim
# target <-
#   loadTargetFile(
#     targetFile = targetFile,
#     varInt = varInt,
#     condRef = condRef,
#     batch = batch
#   )

#### this a screwy function; i prefer just to read the file in with read.delim
target = read.delim(targetFile,header=T,stringsAsFactors = F)


group = relevel(factor(target[, varInt]), ref = condRef)

target[,varInt] = group

rownames(target) = target$name

# loading counts

counts = mybarseqrows(readRDS(countFiles))

pool = "noness"

library(dplyr)

xpool = fdata %>% filter(essential == pool)

counts <- counts[xpool$strain, target$orig]

colnames(counts) = target$name

counts = counts[-which(rowSums(counts) == 0),]

counts = mysumtags(counts)

counts = counts[, target$name]

#target$name=as.character(target$name)
rownames(target) = target$name


# checking parameters

checkParameters.DESeq2(
  projectName = projectName,
  author = author,
  targetFile = targetFile,
  rawDir = rawDir,
  featuresToRemove = featuresToRemove,
  varInt = varInt,
  condRef = condRef,
  batch = batch,
  fitType = fitType,
  cooksCutoff = cooksCutoff,
  independentFiltering = independentFiltering,
  alpha = alpha,
  pAdjustMethod = pAdjustMethod,
  typeTrans = typeTrans,
  locfunc = locfunc,
  colors = mycolors
)

################################################################################
###                             running script                               ###
################################################################################

setwd(workDir)

library(SARTools)

if (forceCairoGraph)
  options(bitmapType = "cairo")


# description plots
majSequences <-
  descriptionPlots(counts = counts,
    group = target[, varInt],
    col = mycolors)

# analysis with DESeq2
out.DESeq2 <-
  run.DESeq2(
    counts = counts,
    target = target,
    varInt = varInt,
    batch = batch,
    locfunc = locfunc,
    fitType = fitType,
    pAdjustMethod = pAdjustMethod,
    cooksCutoff = cooksCutoff,
    independentFiltering = independentFiltering,
    alpha = alpha
  )

out.DESeq2$results= out.DESeq2$results[grep(condRef,out.DESeq2$results)]
#### alternatively for limiting all pairwise comparisons, for lfcShrink and setting lfcThresholds, try this (satisfaction not guaranteed) -- these parameters are definitely an art; so it would be good to systematically test and fine tune the parameters

out.DESeq3 <-
  myrunDESeq3(
    counts = counts,
    target = target,
    varInt = varInt,
    batch = batch,
    locfunc = locfunc,
    fitType = fitType,
    pAdjustMethod = pAdjustMethod,
    cooksCutoff = cooksCutoff,
    independentFiltering = independentFiltering,
    alpha = alpha,
    lfcThreshold = 0,
    Shrink = T, type = "apeglm",
    condRef = "ctrl"
  )


### MDS cluster PCA
exploreCounts(
  object = out.DESeq3$dds,
  group = target[, varInt],
  typeTrans = typeTrans,
  col = mycolors
)

mySAR_QCplots(counts = counts,group=group)

mySARplots(obj = out.DESeq3,condRef = condRef, limit = 60, add = F)

##### i don't like running this -- its a lot of data and i never use it -- breaks too
#####
# summaryResults <-
#   summarizeResults.DESeq2(
#     out.DESeq2,
#     group = target[, varInt],
#     col = mycolors,
#     independentFiltering = independentFiltering,
#     cooksCutoff = cooksCutoff,
#     alpha = alpha
#   )

##### instead i just use these
#####

diagSizeFactorsPlots(out.DESeq3$dds,group,col=mycolors)
dispersionsPlot(out.DESeq3$dds)

complete = out.DESeq3$results
complete = lapply(complete,function(x) {
  colnames(x) = gsub("svalue","pvalue",colnames(x))
  x
}
)
rawpHist(complete)
countsBoxplots(out.DESeq3$dds,group,col = mycolors)
myclusterPlot(counts,group,col=mycolors)

save.image(file = paste0(projectName, ".RData"))
#savehistory(file = paste0(projectName, ".Rhistory"))

check =  list.files("figures")

figs = c(
  'barplotTotal.png',
  'barplotNull.png',
  'densplot.png',
  'dens_norm.png',
  'majSeq.png',
  'pairwiseScatter.png',
  'cluster.png',
  'PCA.png',
  'MDS.png',
  'diagSizeFactorsHist.png',
  'diagSizeFactorsTC.png',
  'countsBoxplots.png'             ,
  'dispersionsPlot.png',
  'rawpHist.png',
  'MAPlot.png',
  'volcanoPlot.png',
  'FDPlots.png'
)

print(setdiff(figs, check))

mywritereport(target, counts, out.DESeq2, summaryResults, majSequences,
  workDir, projectName, author, targetFile, rawDir, featuresToRemove,
  varInt, condRef, batch, fitType, cooksCutoff, independentFiltering,
  alpha, pAdjustMethod, typeTrans, locfunc, colors = mycolors)

