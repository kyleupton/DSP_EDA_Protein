

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
Libraries <- c('edgeR', 'optparse')
lapply(Libraries, require, character.only = TRUE)
rm(Libraries)

###############################################################################################
###############################  Read-in arguments for DGE run ###############################
option_list = list(
  make_option(c("-d", "--rootdir"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-n", "--normpath"), type="character", default='Normalisation/NSNormDropped/', 
              help="dataset file name", metavar="character"),
  make_option(c("-f", "--file"), type="character", default='NanoStringNorm_49_none_none_low.cv.geo.mean.csv', 
              help="dataset file name", metavar="character"),
  make_option(c("-e", "--exportdir"), type="character", default="EdgeR", 
              help="dataset file name", metavar="character"),
  make_option(c("-r", "--runname"), type="character", default="Default", 
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
  print(y$samples)
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
  points <- c(0,1,2,3,4,7,8)
  # points <- rep(c(0,1,2,3,4,7,8,9,10,11,15,16,17,18,19,20,21),2)
  # colors <- c("blue", "blue", "blue", "blue", "blue", "darkgreen", "darkgreen", "darkgreen", "darkgreen", "darkgreen", "red", "red", "red", "red", "red")
  # colors <- c("blue", "blue", "blue", "darkgreen", "darkgreen", "darkgreen", "red", "red", "red", "red", "red")
  
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






###############################################################################################
###############################  Read-in config ###############################

configDir = '/Users/upton6/Documents/notebooks/Nanostring/Larisa_Spheroids/DSP_EDA_Protein'
setwd(configDir)
config = readLines("EdgeR_Config.txt")
# data[1]
# as.list(strsplit(data[1], ','))
groupHead = as.list(strsplit(config[1], ',')[[1]])[1]
comp1 = as.list(strsplit(config[2], ',')[[1]])[1]
compName = as.list(strsplit(config[3], ',')[[1]])[1]


groupHead = as.list(strsplit(config[grepl("GROUP:", config)]))#[1], ',')[[1]])[1]
groupHead

data = readLines("file.txt")
data[grepl("dfw: ", data)]


###############################################################################################
###############################  Read-in probe counts ###############################
# rootDir = opt$rootdir
rootDir = '/Users/upton6/Documents/Nanostring/projects/Larisa/2312_Run/DSP_Protein_Data/'

normFile = file.path(rootDir, opt$normpath, opt$file)
# normFile = file.path(rootDir, 'Normalisation/NSNormDropped/NanoStringNorm_49_none_none_low.cv.geo.mean.csv')
rootDir
normFile
opt$normpath

setwd(rootDir)
# getwd()
# print(normFile)
raw.data <- read.table(file = normFile,
                       header = TRUE,
                       sep=",")
head( raw.data )
dim(raw.data)
counts <- raw.data[ , c(2:dim(raw.data)[2]) ]
dim(counts)
head(counts)
rownames( counts ) <- raw.data[ , 1 ] # gene names


###############################################################################################
###############################  Read-in sample annotations ###############################
sampleInfo = opt$sampleinfo 
info <- read.delim(sampleInfo, 
                      sep=",", 
                      header=TRUE)
targets <- info[ , c(2:dim(info)[2]) ]
rownames( targets ) <- info[ , 1 ]
targets = as.data.frame(t(targets))
targets


###############################################################################################
########################################  Set up model ########################################

exportDir = file.path(rootDir, opt$exportdir)
# exportDir = opt$exportdir
# exportDir
dir.create(exportDir, showWarnings = FALSE)
# getwd()
setwd(exportDir)

#### Make sure pos and neg controls have been dropped!!!

groupHead = paste(groupHead)
# targets[groupHead][,1]



group <- factor(paste0(targets[groupHead][,1]))




# group <- factor(paste0(targets$"Substrate"))
# group <- factor(paste0("group", ".", targets$Substrate, ".", targets$Diff))
# group <- factor(paste0("group", ".", targets$Substrate, ".", targets$Diff, ".", targets$Layer))
group


# y <- DGEList(counts, group = group)
# # y <- calcNormFactors(y)
# y$samples
# # y$norm.factors


y <- setup_y(counts,group)



plot_First(y)
# pdf( paste("plotMD_", opt$runname, ".pdf", sep = "") , width = 4 , height = 4 ) # in inches
# #### Plot first sample
# plotMD(cpm(y, log=TRUE), column=1)
# abline(h=0, col="red", lty=2, lwd=2)
# dev.off() # this tells [R] to close and stop writing to the pdf.


#### Plot Multi Dimensional Scaling
#### Points and Colours need work
plot_MDS(y,group)
# pdf( paste("MDS_Plot_", opt$runname, ".pdf", sep = "") , width = 8 , height = 8 ) # in inches
# #set up to use colours for grade and symbols for label
# # points <- c(0,1,2,15,16,17)
# # points <- c(0,1,2,3,4,7,8,9,10,11,15,16,17,18,19,20,21)
# # points <- rep(c(0,1,2,3,4,5,6),2)
# points <- c(0,1,2,3,4,7,8)
# # points <- rep(c(0,1,2,3,4,7,8,9,10,11,15,16,17,18,19,20,21),2)
# # colors <- c("blue", "blue", "blue", "blue", "blue", "darkgreen", "darkgreen", "darkgreen", "darkgreen", "darkgreen", "red", "red", "red", "red", "red")
# colors <- c("blue", "blue", "blue", "darkgreen", "darkgreen", "darkgreen", "red", "red", "red", "red", "red")
# plotMDS(y, col=colors[group], pch=points[group])
# legend("topleft", legend=levels(group), col=colors[group], pch=points, ncol=1)
# dev.off() # this tells [R] to close and stop writing to the pdf.


#### Is this the best design set-up? Double check the use of 0 vs other options
design <- make_Design(group)
# design <- model.matrix(~0 + group)
# colnames(design) <- levels(group)
design

y <- make_Disp(y, design)
fit <- make_Fit(y, design)

# pdf( paste(opt$runname, "_plotBCV.pdf", sep = "") , width = 4 , height = 4 ) # in inches
# y <- estimateDisp(y, design, robust=TRUE)
# y$common.dispersion
# plotBCV(y)
# dev.off()

# 
# pdf( paste("plotQLDisp_", opt$runname, ".pdf", sep = "") , width = 4 , height = 4 ) # in inches
# fit <- glmQLFit(y, design, robust=TRUE)
# head(fit$coefficients)
# plotQLDisp(fit)
# dev.off()


###############################################################################################
###############################  Run model and export results  ###############################
paste(comp1)


con <- makeContrasts(paste(comp1), levels=design)
# con <- makeContrasts(Hydrogel - Spheroid, levels=design)

pdfName = paste("MD_plot", compName, ".pdf", sep="")
filename = paste("MD_plot", compName, ".csv", sep="")
pdfName2 = paste("MD_plot", compName, "_tr.pdf", sep="")
filename2 = paste("MD_plot", compName, "_tr.csv", sep="")


pdfName <- "MD_plot_Hydrogel_vs_Spheroid.pdf" #, width = 4 , height = 4 ) # in inches
filename <- "MD_plot_Hydrogel_vs_Spheroid.csv"
pdfName2 <- "MD_plot_Hydrogel_vs_Spheroid_tr.pdf"
filename2 <- "MD_plot_Hydrogel_vs_Spheroid_tr.csv"
qlf <- get_Results_QLF(fit, con, pdfName)
print_Results(qlf, y, filename)
tr <- get_Results_TR(fit, con, pdfName2)
print_Results(tr, y, filename2)




# 
# 
# 
# 
# con <- makeContrasts(group.Hydrogel.BDNF - group.Hydrogel.SD, levels=design)
# pdf( "MD_plot_Hydrogel_BDNF_vs_Hydrogel_SD.pdf" , width = 4 , height = 4 ) # in inches
# filename <- "MD_plot_Hydrogel_BDNF_vs_Hydrogel_SD.csv"
# 
# 
# con <- makeContrasts(group.Hydrogel.Dream - group.Hydrogel.SD, levels=design)
# pdf( "MD_plot_Hydrogel_Dream_vs_Hydrogel_SD.pdf" , width = 4 , height = 4 ) # in inches
# filename <- "MD_plot_Hydrogel_Dream_vs_Hydrogel_SD.csv"
# 
# 
# con <- makeContrasts(group.Hydrogel.PDGF - group.Hydrogel.SD, levels=design)
# pdf( "MD_plot_Hydrogel_PDGF_vs_Hydrogel_SD.pdf" , width = 4 , height = 4 ) # in inches
# filename <- "MD_plot_Hydrogel_PDGF_vs_Hydrogel_SD.csv"
# 
# 
# con <- makeContrasts(group.Hydrogel.Nil - group.Hydrogel.SD, levels=design)
# pdf( "MD_plot_Hydrogel_PDGF_vs_Hydrogel_SD.pdf" , width = 4 , height = 4 ) # in inches
# filename <- "MD_plot_Hydrogel_PDGF_vs_Hydrogel_SD.csv"
# 
# 
# con <- makeContrasts(group.Spheroid.Hep - group.Spheroid.SD, levels=design)
# pdf( "MD_plot_Spheroid.Hep_vs_Spheroid_SD.pdf" , width = 4 , height = 4 ) # in inches
# filename <- "MD_plot_Spheroid.Hep_vs_Spheroid_SD.csv"
# 
# 
# 
# con <- makeContrasts(group.Spheroid.Hep.Inner - group.Spheroid.SD.Inner, levels=design)
# pdf( "MD_plot_Spheroid.Hep.Inner_vs_Spheroid.SD.Inner.pdf" , width = 4 , height = 4 ) # in inches
# filename <- "MD_plot_Spheroid.Hep.Inner_vs_Spheroid.SD.Inner.csv"
# 
# 
# con <- makeContrasts(group.Spheroid.Hep.Outer - group.Spheroid.SD.Outer, levels=design)
# pdf( "MD_plot_Spheroid.Hep.Outer_vs_Spheroid.SD.Outer.pdf" , width = 4 , height = 4 ) # in inches
# filename <- "MD_plot_Spheroid.Hep.Outer_vs_Spheroid.SD.Outer.csv"
# 
# 
# con <- makeContrasts(group.Spheroid.SD.Inner - group.Spheroid.SD.Outer, levels=design)
# pdf( "MD_plot_Spheroid.SD.Inner_vs_Spheroid.SD.Outer.pdf" , width = 4 , height = 4 ) # in inches
# filename <- "MD_plot_Spheroid.SD.Inner_vs_Spheroid.SD.Outer.csv"
# 
# con <- makeContrasts(group.Spheroid.Hep.Inner - group.Spheroid.Hep.Outer, levels=design)
# pdf( "MD_plot_Spheroid.Hep.Inner_vs_Spheroid.Hep.Outer.pdf" , width = 4 , height = 4 ) # in inches
# filename <- "MD_plot_Spheroid.Hep.Inner_vs_Spheroid.Hep.Outer.csv"
# 
# 
# 
# qlf <- glmQLFTest(fit, contrast=con)
# topTags(qlf, 15)
# summary(decideTests(qlf))
# 
# ## Save plots and write results to file?
# # pdf( "MD_plot_Tumour_1_Tumour_0.pdf" , width = 4 , height = 4 ) # in inches
# par( mfrow=c(1 ,1) )
# plotMD(qlf)
# dev.off() # this tells [R] to close and stop writing to the pdf.
# 
# resultsbyP <- topTags(qlf, n = nrow(qlf$table))$table
# wh.rows.glm <- match( rownames( resultsbyP ) , rownames( y$counts ) )
# results2.tbl <- cbind (resultsbyP, "Tgw.Disp"=y$tagwise.dispersion[wh.rows.glm], "UpDown" = decideTestsDGE(qlf)[wh.rows.glm,], y$counts[wh.rows.glm,] )
# head (results2.tbl)
# write.table(results2.tbl, file = filename, sep = ",", row.names = TRUE)














