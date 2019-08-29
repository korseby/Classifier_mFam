
# Script for the compilation of classifiers translating MS/MS spectrum information to metabolit family information

## Install

Install dependencies.
```
## matrix
install.packages("SparseM")
install.packages("slam")

## ??
install.packages("plotrix")
BiocManager::install("Biobase")
install.packages("MASS")
library("RMassBank")

## util
install.packages("stringr")
install.packages("stringi")
install.packages("ROCR")
install.packages("PRROC")

## classifier
install.packages("klaR")
install.packages("e1071")
install.packages("kohonen")
install.packages("nnet")
install.packages("rda")
install.packages("caret")
install.packages("caretEnsemble")

## MetFamily project parsing
install.packages("squash")
```



## Run classifier

<!-- Download MetFamily source. -->
<!--`git clone  https://github.com/ipb-halle/MetFamily.git` -->

Source the MetFamily source.
```
library("RCurl")
eval(parse(text = getURL("https://github.com/ipb-halle/MetFamily/blob/master/R/FragmentMatrixFunctions.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://github.com/ipb-halle/MetFamily/blob/master/R/Annotation.R", ssl.verifypeer = FALSE)))
eval(parse(text = getURL("https://github.com/ipb-halle/MetFamily/blob/master/R/DataProcessing.R", ssl.verifypeer = FALSE)))
#source("/home/htreutle/Code/Java/MetFamily/R/FragmentMatrixFunctions.R")
#source("/home/htreutle/Code/Java/MetFamily/R/Annotation.R")
#source("/home/htreutle/Code/Java/MetFamily/R/DataProcessing.R")
```

Source the classifier source files.
```
source("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam_code/Classifier_mFam/SubstanceClassClassifier.R")
source("/mnt/ifs/data/IPB/Projects/2017_005_MS-databases/mFam_code/Classifier_mFam/SubstanceClassClassifier_classifier.R")
```

```
##############################################################################
## general parameters

## cache
reParseMsp              <- FALSE
reProcessAnno           <- FALSE
reProcessFragmentMatrix <- FALSE
reProcessLibrary        <- FALSE

## I/O
resultFolderForClassifiers       <- "/home/htreutle/Data/SubstanceClasses/Classifier_ROC_Analysis"
resultFolderForMetfamilyProjects <- "/home/htreutle/Data/SubstanceClasses/MetFamily_class_projects"

progress    <- FALSE
outDetail    <- FALSE
maximumNumberOfScores <- 10000
removeRareFragments <- FALSE

builtClassifiers       <- TRUE
writeMetFamilyProjects <- TRUE
calculateRMDstatistics <- FALSE
computeFragmentFisherTest <- FALSE

writeClassifiers <- TRUE
writeResults     <- TRUE

outSuffix <- ""

## replicates
mergeDuplicatedSpectra              <- TRUE
takeSpectraWithMaximumNumberOfPeaks <- FALSE

## classes
classOfClass <- c(
  "ChemOnt|MainClass",   # 1
  "ChemOnt|AllClasses",  # 2
  "ChemOnt|Substituent", # 3
  "ChEBI|Identified",    # 4
  "ChEBI|Predicted",     # 5
  "Scaffold"             # 6
)[[2]]
#thisClass <- "Organic compounds; Lipids and lipid-like molecules; Prenol lipids; Terpene lactones; Sesquiterpene lactones"
thisClass <- NULL

## repeated random subsampling validation
minimumNumberOfPosSpectraPerClass <- 10
minimumNumberOfNegSpectraPerClass <- 10
numberOfDataSetDecompositions <- 10
proportionTraining <- 0.7
#fragmentColumnSelectionMethod <- "AbsoluteProportion"
#fragmentColumnSelectionMethod <- "ProportionOfSumOfFragmentPresence"
fragmentColumnSelectionMethod <- "ProportionOfHighestFragmentPresence"
minimumProportionOfPositiveFragments <- 0.05

## postprocessing
unitResolution <- FALSE

#######################################################################
## constraints

## library one
#thisLibrary    <- NULL
thisLibrary    <- "180912_MoNA-export-LC-MS-MS_Negative_Mode_processed.msp"
annoFile <- "/home/htreutle/Downloads/MetSWATH/MONA/190523_MSMS_HR_someScaffolds_completed.tsv"

splitClasses <- FALSE
splitClassesParameterSet <- list(
  minimumNumberOfChildren = 5,
  maximumNumberOfPrecursors = 1000,
  distanceMeasure = "Jaccard (intensity-weighted)",
  #distanceMeasure = "Jaccard",
  clusterMethod = "ward.D"
)

thisMethod <- "method=ColSumsPos; smoothIntensities=FALSE"
#thisMethod <- "method=caret; smoothIntensities=FALSE, modelName=binda"

thisParameterSet <- list(
  ## parse msp
  minimumIntensityOfMaximalMS2peak                  = 000,
  minimumProportionOfMS2peaks                       = 0.05,
  neutralLossesPrecursorToFragments                 = TRUE,
  neutralLossesFragmentsToFragments                 = FALSE,
  ## built fragment matrix
  mzDeviationAbsolute_grouping                      = 0.01,
  mzDeviationInPPM_grouping                         = 10,
  #doMs2PeakGroupDeisotoping                         = TRUE,
  doMs2PeakGroupDeisotoping                         = FALSE,
  mzDeviationAbsolute_ms2PeakGroupDeisotoping       = 0.01,
  mzDeviationInPPM_ms2PeakGroupDeisotoping          = 10,
  proportionOfMatchingPeaks_ms2PeakGroupDeisotoping = 0.9
)

##############################################################################
## box parameters
parameterSetAll <- list(
  ##############################################################################
  ## general parameters
  
  ## I/O
  progress              = progress,
  outDetail             = outDetail,
  resultFolderForClassifiers          = resultFolderForClassifiers,
  resultFolderForMetfamilyProjects = resultFolderForMetfamilyProjects,
  writeClassifiers      = writeClassifiers,
  writeResults          = writeResults,
  outSuffix             = outSuffix,
  maximumNumberOfScores = maximumNumberOfScores,
  removeRareFragments   = removeRareFragments,
  builtClassifiers      = builtClassifiers,
  writeMetFamilyProjects = writeMetFamilyProjects,
  calculateRMDstatistics = calculateRMDstatistics,
  ## classes
  mergeDuplicatedSpectra              = mergeDuplicatedSpectra,
  takeSpectraWithMaximumNumberOfPeaks = takeSpectraWithMaximumNumberOfPeaks,
  classOfClass = classOfClass,
  thisClass = thisClass,
  splitClasses = splitClasses,
  splitClassesParameterSet = splitClassesParameterSet,
  computeFragmentFisherTest         = computeFragmentFisherTest,
  ## repeated random subsampling validation
  minimumNumberOfPosSpectraPerClass = minimumNumberOfPosSpectraPerClass,
  minimumNumberOfNegSpectraPerClass = minimumNumberOfNegSpectraPerClass,
  numberOfDataSetDecompositions     = numberOfDataSetDecompositions,
  proportionTraining                = proportionTraining,
  fragmentColumnSelectionMethod     = fragmentColumnSelectionMethod,
  minimumProportionOfPositiveFragments = minimumProportionOfPositiveFragments,
  ## postprocessing
  unitResolution                    = unitResolution,
  #######################################################################
  ## cache
  reParseMsp = reParseMsp,
  reProcessAnno = reProcessAnno,
  reProcessLibrary = reProcessLibrary,
  reProcessFragmentMatrix = reProcessFragmentMatrix,
  
  #######################################################################
  ## constraints
  
  ## library one
  annoFile = annoFile,
  thisLibrary = thisLibrary,
  thisParameterSet = thisParameterSet,
  thisMethod  = thisMethod
)

runTest(parameterSetAll = parameterSetAll)
```

