###############################################################################################
###############################      Author : Kyle R. Upton      ##############################
###############################################################################################

### R version 4.3.2 (2023-10-31)

### To Do.

# Review section 3 (Specific Experimental Designs) in User guide
# Run through with basic group-group comparisons (using exact test)

# Build up glm approach

###############################################################################################
################################### Required R packages #######################################
### Ensure all libraries are installed
Libraries <- c('edgeR', 'optparse', 'stringr')
print('Checking libraries:')
lapply(Libraries, require, character.only = TRUE)
rm(Libraries)
print('If any values above are FALSE, library availability must be checked')
print(' ')

###############################################################################################
###############################  Read-in arguments for DGE run ###############################

option_list = list(
  make_option(c("-c", "--configPath"), type="character", default='/Users/upton6/Library/CloudStorage/OneDrive-QueenslandUniversityofTechnology/Documents/notebooks/Nanostring/Adams_Bray/DSP_EDA_Protein', 
              help="dataset file name", metavar="character"),
  make_option(c("-d", "--rootdir"), type="character", default='/Users/upton6/Documents/Nanostring/projects/Adams/', 
              help="dataset file name", metavar="character"),
  make_option(c("-n", "--normpath"), type="character", default='Normalisation/NSNormDropped', 
              help="dataset file name", metavar="character"),
  make_option(c("-f", "--file"), type="character", default='NanoStringNorm_49_none_none_low.cv.geo.mean.csv', 
              help="dataset file name", metavar="character"),
  make_option(c("-e", "--exportdir"), type="character", default='EdgeR', 
              help="dataset file name", metavar="character"),
  make_option(c("-r", "--runname"), type="character", default='Default', 
              help="dataset file name", metavar="character"),
  make_option(c("-i", "--sampleinfo"), type="character", default='sampleInfo_with_Wells.csv', 
              help="dataset file name", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

# opt$runname

###############################################################################################
#####################################  Set up Functions  ######################################

setup_y <- function(counts,group){
  y <- DGEList(counts, group = group)
  # print(y$samples)
  return(y)
}

plot_First <- function(y){
  #### Plot first sample
  pdf( paste("plotMD_", opt$runname, ".pdf", sep = "") , width = 4 , height = 4 ) # in inches
  plotMD(cpm(y, log=TRUE), column=2)
  abline(h=0, col="red", lty=2, lwd=2)
  dev.off()
}

plot_MDS <- function(y, group){
  #### Plot Multi Dimensional Scaling
  #### Points and Colours need work
  pdf( paste("MDS_Plot_", opt$runname, ".pdf", sep = "") , width = 8 , height = 8 ) # in inches
  # points <- c(0,1,2,3,4,7,8)
  points <- rep(c(0,1,2,3,4,7,8,9,10,11,15,16,17,18,19,20,21),2)
  # colors <- c("blue", "blue", "blue", "blue", "blue", "darkgreen", "darkgreen", "darkgreen", "darkgreen", "darkgreen", "red", "red", "red", "red", "red")
  plotMDS(y, pch=points[group])
  legend("topleft", legend=levels(group), pch=points, ncol=1)
  dev.off()
}

make_Design <- function(group){
  #### Is this the best design set-up? Double check the use of 0 vs other options
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  return(design)
}

make_Disp <- function(y, design){
  pdf( paste(opt$runname, "_plotBCV.pdf", sep = "") , width = 4 , height = 4 ) # in inches
  y <- estimateDisp(y, design, robust=TRUE)
  y$common.dispersion
  plotBCV(y)
  dev.off()
  return(y)
}

make_Fit <- function(y, design){
  pdf( paste(opt$runname, "_plotQLDisp.pdf", sep = "") , width = 4 , height = 4 ) # in inches
  fit <- glmQLFit(y, design, robust=TRUE)
  head(fit$coefficients)
  plotQLDisp(fit)
  dev.off()
  return(fit)
}

get_Results_QLF <- function(fit, con, pdfName){
  qlf <- glmQLFTest(fit, contrast=con)
  print(topTags(qlf))
  print(summary(decideTests(qlf)))
  pdf(pdfName , width = 4 , height = 4 ) # in inches
  par( mfrow=c(1 ,1))
  plotMD(qlf)
  dev.off()
  return(qlf)
}

get_Results_TR <- function(fit, con, pdfName){
  tr <- glmTreat(fit, contrast=con, lfc=log2(1.2))
  print(topTags(tr))
  print(summary(decideTests(tr)))
  pdf(pdfName , width = 4 , height = 4 ) # in inches
  par( mfrow=c(1 ,1))
  plotMD(tr)
  dev.off()
  return(tr)
}


print_Results <- function(qlf, y, filename){
  resultsbyP <- topTags(qlf, n = nrow(qlf$table))$table
  resultsbyP
  wh.rows.glm <- match( rownames( resultsbyP ) , rownames( y$counts ) )
  results2.tbl <- cbind (resultsbyP, "Tgw.Disp"=y$tagwise.dispersion[wh.rows.glm], "UpDown" = decideTestsDGE(qlf)[wh.rows.glm,], y$counts[wh.rows.glm,] )
  # head (results2.tbl)
  write.table(results2.tbl, file = filename, sep = ",", row.names = TRUE)
}

parse_input <- function(vector){
  result = list()
  for (i in seq_along(vector)){
    c <- as.vector(strsplit(vector[i],':')[[1]])[2]
    c <- as.vector(strsplit(c,',')[[1]])
    result[[i]] <- c
  }
  return(result)
}


###############################################################################################
###############################  Read-in config ###############################

configDir = opt$configPath
setwd(configDir)
config = readLines("EdgeR_Config.txt")

groupHead = as.vector(config[grepl("GROUP:", config)])
# groupHead
groups = list()
for (i in seq_along(groupHead)){
  # print(groupHead[i])
  g <- as.vector(strsplit(groupHead[i],":")[[1]])[2]
  g <- as.vector(strsplit(g,"[.]")[[1]])
  groups[[i]] <- g
}
print('groups')
print(groups)

compRaw = as.vector(config[grepl("COMPARISON:", config)])
comps = parse_input(compRaw)
# comps

compNameRaw = as.vector(config[grepl("COMP_NAME:", config)])
compNames = parse_input(compNameRaw)
# compNames


###############################################################################################
###############################  Read-in probe counts ###############################
rootDir = opt$rootdir

normFile = file.path(rootDir, opt$normpath, opt$file)
rootDir
normFile
opt$normpath

setwd(rootDir)
# getwd()
# print(normFile)
raw.data <- read.table(file = normFile,
                       header = TRUE,
                       sep=",")
# head( raw.data )
# dim(raw.data)
counts <- raw.data[ , c(2:dim(raw.data)[2]) ]
# dim(counts)
# head(counts)
rownames( counts ) <- raw.data[ , 1 ] # gene names

keeps <- c(colnames( counts ))


###############################################################################################
###############################  Read-in sample annotations ###############################
sampleInfo = opt$sampleinfo 
infoRaw <- read.delim(sampleInfo, 
                      sep=",", 
                      header=TRUE)

info <- infoRaw[keeps]


targets <- info[ , c(1:dim(info)[2]) ]
info[ , 1 ]
rownames( targets ) <- infoRaw[ , 1 ]
targets = as.data.frame(t(targets))
# targets


###############################################################################################
########################################  Set up model ########################################

exportDir = file.path(rootDir, opt$exportdir)
dir.create(exportDir, showWarnings = FALSE)
setwd(exportDir)

groups
for (g in seq_along(groups)){
  # print(groups[g][[1]])
  test2 <- mapply(function(x,y) targets[x][,y], groups[g][[1]], 1)
  group = vector()
  for (x in suppressWarnings(1:range(dim(test2)[1]))){
    group <- append(group, str_c(test2[x,],collapse="."))
  }
  group <- factor(group)
  y <- setup_y(counts,group)
  plot_First(y)
  plot_MDS(y,group)
  
  #### Is this the best design set-up? Double check the use of 0 vs other options
  design <- make_Design(group)
  y <- make_Disp(y, design)
  fit <- make_Fit(y, design)
  
  
  ###############################################################################################
  ###############################  Run model and export results  ###############################
  
  for (c in 1:length(comps[[g]])){
    thisComp = comps[[g]][c]
    compName = compNames[[g]][c]                      
    con <- makeContrasts(paste(thisComp), levels=design)

    pdfName = paste("MD_plot_", compName, ".pdf", sep="")
    filename = paste("MD_plot_", compName, ".csv", sep="")
    pdfName2 = paste("MD_plot_", compName, "_tr.pdf", sep="")
    filename2 = paste("MD_plot_", compName, "_tr.csv", sep="")
    
    qlf <- get_Results_QLF(fit, con, pdfName)
    print_Results(qlf, y, filename)
    tr <- get_Results_TR(fit, con, pdfName2)
    print_Results(tr, y, filename2)
  }
}  
