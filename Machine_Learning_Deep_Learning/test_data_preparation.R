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

## Assigning class labels
labs <- c("Enhancer","Non-Enhancer")
merged_bw$class <- sample(labs,nrow(merged_bw),replace = TRUE)
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

## For fair results, let's make sure that no "NA" values exist.
na.omit(test)

## Min-Max normalization

test$peaks_h3k4me3 <- (test$peaks_h3k4me3-min(test$peaks_h3k4me3, na.rm = T))/(max(test$peaks_h3k4me3, na.rm = T)-min(test$peaks_h3k4me3, na.rm = T))
test$peaks_h3k4me2 <- (test$peaks_h3k4me2-min(test$peaks_h3k4me2, na.rm = T))/(max(test$peaks_h3k4me2, na.rm = T)-min(test$peaks_h3k4me2, na.rm = T))
test$peaks_h3k4me1 <- (test$peaks_h3k4me1-min(test$peaks_h3k4me1, na.rm = T))/(max(test$peaks_h3k4me1, na.rm = T)-min(test$peaks_h3k4me1, na.rm = T))
test$peaks_h3k27ac <- (test$peaks_h3k27ac-min(test$peaks_h3k27ac, na.rm = T))/(max(test$peaks_h3k27ac, na.rm = T)-min(test$peaks_h3k27ac, na.rm = T))

# Data Partition
set.seed(108)
ind <- sample(2, nrow(test), replace = TRUE, prob = c(0.8, 0.2))
training <- test[ind==1,]
testing <- test[ind==2,]


# Neural Networks
set.seed(007)
nn <- neuralnet(class~peaks_h3k4me3+peaks_h3k4me2+peaks_h3k4me1+peaks_h3k27ac,
                data = training,
                hidden = 5,
                err.fct = "sse",
                act.fct = "logistic",
                linear.output = FALSE)
plot(nn)