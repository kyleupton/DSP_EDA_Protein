
###############################################################################################
###############################      Author : Kyle R. Upton      ##############################
###############################################################################################

### R version 4.3.2 (2023-10-31)

###############################################################################################
################################### Required R packages #######################################
# BiocManager::install("NormqPCR")
# BiocManager::install("ComplexHeatmap")
# BiocManager::install("vsn")
# BiocManager::install("edgeR")
# # NanoStringNorm_1.2.1.1.tar.gz is required but must be manually downloaded and installed. May require XQuartz on mac.

### Ensure all libraries are installed
Libraries <- c('matrixStats', 'ruv', 'NanoStringNorm', 'data.table','NormqPCR', 'optparse',
               'ComplexHeatmap','dplyr', 'plotrix', 'tuple', 'scales', 'ggplot2', 'stringr')
lapply(Libraries, require, character.only = TRUE)
rm(Libraries)

###############################################################################################
##########################  Read-in arguments for normalisation run ###########################
option_list = list(
  make_option(c("-d", "--datadir"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-s", "--subdir"), type="character", default="NSNorm", 
              help="dataset file name", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

###############################################################################################
####################################  Read-in probe counts ####################################
#***** Read in Nanostring raw/QC data
dataDir = opt$datadir
dataFile = opt$file

setwd(dataDir)
Nano_ExpressionMatrix <- read.delim(dataFile, 
                                    sep=',',
                                    stringsAsFactors = FALSE, 
                                    header = TRUE, 
                                    as.is = TRUE)
dim(Nano_ExpressionMatrix)

# #############################################################################################
# ###############  All 84 different normalization - using NanoStringNorm package ##############
exportDir = file.path(dataDir, opt$subdir)
dir.create(exportDir, showWarnings = FALSE)
setwd(exportDir)

#***** All different Nanostring normalisation options using NanoStringNorm.
### Create a matrix of normalisation options
CodeCount_Option <- c('none', 'sum', 'geo.mean')
BackGround_Option <- c('none', 'mean', 'mean.2sd', 'max')
HK_Option <- c('none', 'housekeeping.sum', 'housekeeping.geo.mean', 'total.sum', 'low.cv.geo.mean', 'top.mean', 'top.geo.mean')
AllOptions <- expand.grid(CodeCount_Option, BackGround_Option, HK_Option)
dim(AllOptions)

### Run normalisation and export data
for(i in 1:84){
  NanoString_mRNA_norm <- NanoStringNorm(x = Nano_ExpressionMatrix,
                               CodeCount = as.character(AllOptions[i , 1]),
                               Background = as.character(AllOptions[i , 2]),
                               SampleContent = as.character(AllOptions[i , 3]),
                               round.values = FALSE,
                               take.log = FALSE ,
                               return.matrix.of.endogenous.probes = TRUE)
  NSN <- paste('NanoStringNorm_', 
               str_pad(i, 2, pad = "0"),'_', 
               as.character(AllOptions[i , 1]), 
               '_',as.character(AllOptions[i , 2]), 
               '_',as.character(AllOptions[i , 3]), '.csv', 
               sep='')
  write.csv(NanoString_mRNA_norm,NSN, row.names=TRUE)
}