## Preparing the Testing Data ##
## Importing the Histone Data ##

h3k27ac <- read.csv("./data/H1_Cell_Line/H3K27ac/ENCFF663SAM.bw", sep = '\t', header = FALSE)
colnames(h3k27ac) <- c("chrom","start","end","peaks")
h3k4me1 <- read.csv("./data/H1_Cell_Line/H3K4me1/ENCFF441KOL.bw", sep = '\t', header = FALSE)
colnames(h3k4me1) <- c("chrom","start","end","peaks")
h3k4me2 <- read.csv("./data/H1_Cell_Line/H3K4me2/ENCFF799BDH.bw", sep = '\t', header = FALSE)
colnames(h3k4me2) <- c("chrom","start","end","peaks")
h3k4me3 <- read.csv("./data/H1_Cell_Line/H3K4me3/ENCFF340UJK.bw", sep = '\t', header = FALSE)
colnames(h3k4me3) <- c("chrom","start","end","peaks")


## Importing the merged BEDGRAPH/BW file

merged_bw <- read.csv("./data/H1_Cell_Line/bedtools_Merge.bw", sep = '\t', header = FALSE)
merged_bw <- merged_bw[-1,]
colnames(merged_bw) <- c("chrom", "start", "end", "peaks_h3k4me3", "peaks_h3k4me2", "peaks_h3k4me1", "peaks_h3k27ac")
head(merged_bw)
merged_bw$class <- "enhancer"
head(merged_bw)


## Identifying examples of standard chromosomes only and filtering the residuals.

chromosomes <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20","chr21", "chr22", "chrX", "chrY")
merged_bw<- as.data.frame(merged_bw[merged_bw$chrom %in% chromosomes, ])


## Deriving test data as input to the deep learning model.

test <- merged_bw[,c(4:8)]
test$peaks_h3k4me3 <- as.double(as.character(test$peaks_h3k4me3))
test$peaks_h3k4me2 <- as.double(as.character(test$peaks_h3k4me2))
test$peaks_h3k4me1 <- as.double(as.character(test$peaks_h3k4me1))
test$peaks_h3k27ac <- as.double(as.character(test$peaks_h3k27ac))
