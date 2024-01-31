

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
  make_option(c("-d", "--datadir"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-n", "--normpath"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-f", "--file"), type="character", default=NULL, 
              help="dataset file name", metavar="character"),
  make_option(c("-e", "--exportdir"), type="character", default="NSNorm", 
              help="dataset file name", metavar="character"),
  make_option(c("-i", "--sampleinfo"), type="character", default=NULL, 
              help="dataset file name", metavar="character")
); 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);


# '/Users/upton6/Documents/Nanostring/projects/Larisa/2312_Run/DSP_Protein_Data/'

###############################################################################################
###############################  Read-in probe counts ###############################
rootDir = opt$datadir

normFile = file.path(rootDir, opt$normpath)
# normFile = opt$normpath

setwd(rootDir)
getwd()
print(normFile)

raw.data <- read.table(file = normFile,
                       header = TRUE,
                       sep=",")

# raw.data <- read.delim(dataFile, 
#                        sep=',',
#                        stringsAsFactors = FALSE, 
#                        header = TRUE, 
#                        as.is = TRUE)

head( raw.data )
dim(raw.data)
counts <- raw.data[ , c(2:120) ]
dim(counts)
head(counts)
rownames( counts ) <- raw.data[ , 1 ] # gene names

###############################################################################################
###############################  Read-in sample annotations ###############################

# exportDir = file.path(dataDir, opt$subdir)
exportDir = opt$exportdir
setwd(exportDir)
targets <- read.delim(opt$sampleinfo, 
                      sep=",", 
                      header=TRUE)
targets

###############################################################################################
########################################  Set up model ########################################

setwd(exportDir)

#### Make sure pos and neg controls have been dropped!!!

# group <- factor(paste0(targets$Grade_b, ".", targets$Grade_c))
group <- factor(paste0("group", ".", targets$Grade_b))
group <- factor(paste0("group", ".", targets$Grade_b, ".", targets$Grade_c))
group

dim(counts)
dim(targets)

y <- DGEList(counts, group = group)
# y <- calcNormFactors(y)
y$samples
# y$norm.factors


pdf( "plotMD.pdf" , width = 4 , height = 4 ) # in inches

#### Plot first sample
plotMD(cpm(y, log=TRUE), column=1)
abline(h=0, col="red", lty=2, lwd=2)

dev.off() # this tells [R] to close and stop writing to the pdf.


#### Plot Multi Dimensional Scaling
#### Points and Colours need work
pdf( "MDS_Plot.pdf" , width = 4 , height = 4 ) # in inches
#set up to use colours for grade and symbols for label
points <- c(0,1,2,15,16,17)
points <- c(0,1,2,3,4,7,8,9,10,11,15,16,17,18,19,20,21)
colors <- c("blue", "blue", "blue", "blue", "blue", "darkgreen", "darkgreen", "darkgreen", "darkgreen", "darkgreen", "red", "red", "red", "red", "red")
points <- rep(c(0,1,2,3,4,5,6),3)
plotMDS(y, col=colors[group], pch=points[group])
legend("topleft", legend=levels(group), col=colors[group], pch=points, ncol=3)
dev.off() # this tells [R] to close and stop writing to the pdf.


#### Is this the best design set-up? Double check the use of 0 vs other options
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
design


pdf( "plotBCV.pdf" , width = 4 , height = 4 ) # in inches
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
plotBCV(y)
dev.off() # this tells [R] to close and stop writing to the pdf.


pdf( "plotQLDisp.pdf" , width = 4 , height = 4 ) # in inches
fit <- glmQLFit(y, design, robust=TRUE)
head(fit$coefficients)
plotQLDisp(fit)
dev.off() # this tells [R] to close and stop writing to the pdf.


### Set up a list of different contrasts to make here

# # Tumour - TME
# con <- makeContrasts(Tumour - TME, levels=design)
# con <- makeContrasts(Tumour - Full_ROI, levels=design)
# con <- makeContrasts(TME - Full_ROI, levels=design)

# help(make.names)





###############################################################################################
###############################  Run model and export results  ###############################





# Grouped Comparisons
con <- makeContrasts(group.3 - group.1, levels=design)
pdf( "MD_plot_GroupedGrade3_GroupedGrade1.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade3_GroupedGrade1.csv"

con <- makeContrasts(group.2 - group.1, levels=design)
pdf( "MD_plot_GroupedGrade2_GroupedGrade1.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade2_GroupedGrade1.csv"

con <- makeContrasts(group.3 - group.2, levels=design)
pdf( "MD_plot_GroupedGrade3_GroupedGrade2.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade3_GroupedGrade2.csv"








con <- makeContrasts(group.1 - group.7, levels=design)
pdf( "MD_plot_GroupedGrade1_GroupedNorm.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade1_GroupedNorm.csv"

con <- makeContrasts(group.1 - group.5, levels=design)
pdf( "MD_plot_GroupedGrade1_GroupedNAT.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade1_GroupedNAT.csv"



con <- makeContrasts(group.2 - group.7, levels=design)
pdf( "MD_plot_GroupedGrade2_GroupedNorm.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade2_GroupedNorm.csv"

con <- makeContrasts(group.2 - group.5, levels=design)
pdf( "MD_plot_GroupedGrade2_GroupedNAT.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade2_GroupedNAT.csv"



con <- makeContrasts(group.3 - group.7, levels=design)
pdf( "MD_plot_GroupedGrade3_GroupedNorm.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade3_GroupedNorm.csv"

con <- makeContrasts(group.3 - group.5, levels=design)
pdf( "MD_plot_GroupedGrade3_GroupedNAT.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade3_GroupedNAT.csv"



con <- makeContrasts(group.3 - group.1, levels=design)
pdf( "MD_plot_GroupedGrade3_GroupedGrade1.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade3_GroupedGrade1.csv"





con <- makeContrasts(group.1.1 - group.5.1, levels=design)
pdf( "MD_plot_GroupedGrade1-1_GroupedNAT-1.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade1-1_GroupedNAT-1.csv"

con <- makeContrasts(group.2.2 - group.5.2, levels=design)
pdf( "MD_plot_GroupedGrade2-2_GroupedNAT-2.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade2-2_GroupedNAT-2.csv"

con <- makeContrasts(group.3.3 - group.5.3, levels=design)
pdf( "MD_plot_GroupedGrade3-3_GroupedNAT-3.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade3-3_GroupedNAT-3.csv"


# con <- makeContrasts(group.1.1 - group.6.1, levels=design)
# pdf( "MD_plot_GroupedGrade1-1_GroupedDAT-1.pdf" , width = 4 , height = 4 ) # in inches
# filename <- "MD_plot_GroupedGrade1-1_GroupedDAT-1.csv"

con <- makeContrasts(group.2.2 - group.6.2, levels=design)
pdf( "MD_plot_GroupedGrade2-2_GroupedDAT-2.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade2-2_GroupedDAT-2.csv"

con <- makeContrasts(group.3.3 - group.6.3, levels=design)
pdf( "MD_plot_GroupedGrade3-3_GroupedDAT-3.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade3-3_GroupedDAT-3.csv"


con <- makeContrasts(group.6.2 - group.5.2, levels=design)
pdf( "MD_plot_GroupedDAT-2_GroupedNAT-2.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedDAT-2_GroupedNAT-2.csv"

con <- makeContrasts(group.6.3 - group.5.3, levels=design)
pdf( "MD_plot_GroupedDAT-3_GroupedNAT-3.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedDAT-3_GroupedNAT-3.csv"





con <- makeContrasts(group.2.2 - group.1.1, levels=design)
pdf( "MD_plot_GroupedGrade2-2_GroupedGrade1-1.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade2-2_GroupedGrade1-1.csv"

con <- makeContrasts(group.3.3 - group.1.1, levels=design)
pdf( "MD_plot_GroupedGrade3-3_GroupedGrade1-1.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade3-3_GroupedGrade1-1.csv"

con <- makeContrasts(group.3.3 - group.2.2, levels=design)
pdf( "MD_plot_GroupedGrade3-3_GroupedGrade2-2.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedGrade3-3_GroupedGrade2-2.csv"



con <- makeContrasts(group.5.1 - group.7.7, levels=design)
pdf( "MD_plot_GroupedNAT-1_Norm.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedNAT-1_Norm.csv"

con <- makeContrasts(group.5.2 - group.7.7, levels=design)
pdf( "MD_plot_GroupedNAT-2_Norm.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedNAT-2_Norm.csv"

con <- makeContrasts(group.5.3 - group.7.7, levels=design)
pdf( "MD_plot_GroupedNAT-3_Norm.pdf" , width = 4 , height = 4 ) # in inches
filename <- "MD_plot_GroupedNAT-3_Norm.csv"






qlf <- glmQLFTest(fit, contrast=con)
topTags(qlf, 15)
summary(decideTests(qlf))

## Save plots and write results to file?
# pdf( "MD_plot_Tumour_1_Tumour_0.pdf" , width = 4 , height = 4 ) # in inches
par( mfrow=c(1 ,1) )
plotMD(qlf)
dev.off() # this tells [R] to close and stop writing to the pdf.

resultsbyP <- topTags(qlf, n = nrow(qlf$table))$table
wh.rows.glm <- match( rownames( resultsbyP ) , rownames( y$counts ) )
results2.tbl <- cbind (resultsbyP, "Tgw.Disp"=y$tagwise.dispersion[wh.rows.glm], "UpDown" = decideTestsDGE(qlf)[wh.rows.glm,], y$counts[wh.rows.glm,] )
head (results2.tbl)
write.table(results2.tbl, file = filename, sep = ",", row.names = TRUE)














