## Preparing the Testing Data ##
## Importing the Histone Data ##

h3k27ac <- read.csv("./data/H1_Cell_Line/H3K27ac/ENCFF663SAM2000_30.bw", sep = '\t', header = FALSE)
colnames(h3k27ac) <- c("chrom","start","end","reads")   
h3k4me3 <- read.csv("./data/H1_Cell_Line/H3K4me3/ENCFF340UJK2000_36.bw", sep = '\t', header = FALSE)
colnames(h3k4me3) <- c("chrom","start","end","reads")
h3k4me2 <- read.csv("./data/H1_Cell_Line/H3K4me2/ENCFF799BDH2000_36.bw", sep = '\t', header = FALSE)
colnames(h3k4me2) <- c("chrom","start","end","reads")
h3k4me1 <- read.csv("./data/H1_Cell_Line/H3K4me1/ENCFF441KOL2000_36.bw", sep = '\t', header = FALSE)
colnames(h3k4me1) <- c("chrom","start","end","reads")


## Convert the dataframes to Granges object

library(GenomicRanges)
h3k27ac_gr <- as(h3k27ac,"GRanges")
h3k4me3_gr <- as(h3k4me3,"GRanges")
h3k4me2_gr <- as(h3k4me2,"GRanges")
h3k4me1_gr <- as(h3k4me1,"GRanges")


## Importing the merged BEDGRAPH/BW file
bedtools_unionbedg_all.bw <- system ("bedtools unionbedg -i ./data/H1_Cell_Line/H3K27ac/ENCFF663SAM2000_30.bw ./data/H1_Cell_Line/H3K4me3/ENCFF340UJK2000_36.bw ./data/H1_Cell_Line/H3K4me2/ENCFF799BDH2000_36.bw ./data/H1_Cell_Line/H3K4me1/ENCFF441KOL2000_36.bw -header -names H3K27ac H3K4me3 H3K4me2 H3K4me1")
merged_bw <- read.csv("bedtools_unionbedg_all.bw", sep = '\t', header = TRUE)
#merged_bw <- merged_bw[-1,]
head(merged_bw)
merged_bw$class <- "enhancer"
head(merged_bw)

## Identifying examples of standard chromosomes only and filtering the residuals.

chromosomes <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20","chr21", "chr22", "chrX", "chrY")
merged_bw<- as.data.frame(merged_bw[merged_bw$chrom %in% chromosomes, ])


## Deriving test data as input to the deep learning model.

test <- merged_bw[,c(4:8)]
test$H3K27ac <- (test$H3K27ac-min(test$H3K27ac, na.rm = T))/(max(test$H3K27ac, na.rm = T)-min(test$H3K27ac, na.rm = T))
test$H3K4me3 <- (test$H3K4me3-min(test$H3K4me3, na.rm = T))/(max(test$H3K4me3, na.rm = T)-min(test$H3K4me3, na.rm = T))
test$H3K4me2 <- (test$H3K4me2-min(test$H3K4me2, na.rm = T))/(max(test$H3K4me2, na.rm = T)-min(test$H3K4me2, na.rm = T))
test$H3K4me1 <- (test$H3K4me1-min(test$H3K4me1, na.rm = T))/(max(test$H3K4me1, na.rm = T)-min(test$H3K4me1, na.rm = T))


## Merging Data on the basis of overlapping intervals

## Positive class data (labels)
positive_class

## Negative class (Labels)
negative_class

## Score Matrix (Input data)
input_score_data <- merged_bw
input_score_data<- as.data.frame(input_score_data[input_score_data$chrom %in% chromosomes, ])


## Replacing NAs with 0s(zeros) in the peaks' columns. Since the score data is not available, imputing empty cells with
## zero entries engenders mathematical convenience.

input_score_data$V4[is.na(input_score_data$V4)] <- 0 
input_score_data$V5[is.na(input_score_data$V5)] <- 0 
input_score_data$V6[is.na(input_score_data$V6)] <- 0 
input_score_data$V7[is.na(input_score_data$V7)] <- 0 

## Converting data to GRanges objects
library(GenomicRanges)
positive_class_labels <- GRanges(seqnames = positive_class$V1, ranges = IRanges(start = positive_class$V2, 
                                                                       end = positive_class$V3))
mcols(positive_class_labels) <- DataFrame(class= "Enhancer")


negative_class_labels <- GRanges(seqnames = negative_class$V1, ranges = IRanges(start = negative_class$V2, 
                                                                       end = negative_class$V3))
mcols(negative_class_labels) <- DataFrame(class= "Non-Enhancer")


input_score <- GRanges(seqnames = input_score_data$V1, ranges = IRanges(start = input_score_data$V2,
                                                                        end = input_score_data$V3))
mcols(input_score) <- DataFrame(peaks_h3k27ac = input_score_data$V4, peaks_h3k4me3 = input_score_data$V5,
                                peaks_h3k4me2 = input_score_data$V6, peaks_h3k4me1 = input_score_data$V7)

## Performing merge to figure out the score and class matrix.

intermatrix <- merge(as.data.frame(positive_class_labels), as.data.frame(negative_class_labels), all= TRUE) ## positive and negative classes ##


## Exporting 'intermatrix1' and 'input_score' as bed files to merge.
write.table(intermatrix, "./data/class_labels.bed", sep ='\t', quote = FALSE, row.names = FALSE)
write.table(input_score, "./data/score.bed", sep ='\t', quote = FALSE, row.names = FALSE)

## Importing the merged file
score_labels <- read.table("./data/score_labels.bed", sep = "\t", header = FALSE, stringsAsFactors=FALSE)

## Picking relevant columns
final_data <- score_labels[,c(6:9,15)]
colnames(final_data)=c("peaks_h3k27ac","peaks_h3k4me3","peaks_h3k4me2","peaks_h3k4me1","class")

## The matrix is a sparse matrix, having many non-zero entries. Let us convert it into a sparse matrix object.
library(Matrix)
final_sparse_data <- final_data
final_sparse_data$class <- as.numeric(as.factor(final_sparse_data$class))-1
final_sparse_data$peaks_h3k27ac <- as.numeric(final_sparse_data$peaks_h3k27ac)
final_sparse_data$peaks_h3k4me3 <- as.numeric(final_sparse_data$peaks_h3k4me3)
final_sparse_data$peaks_h3k4me2 <- as.numeric(final_sparse_data$peaks_h3k4me2)
final_sparse_data$peaks_h3k4me1 <- as.numeric(final_sparse_data$peaks_h3k4me1)

final_sparse_data <- final_sparse_data[complete.cases(final_sparse_data), ]

## Saving this data
saveRDS(final_sparse_data,"./data/ep_data.rds")

## Sample data
final_sparse_data_sample <- final_sparse_data[sample(nrow(final_sparse_data), 10000), ]
saveRDS(final_sparse_data_sample,"./data/ep_data_sample.rds")

## Converting data frame to matrix and then to sparse matrix.
final_sparse_data_matrix <- data.matrix(final_sparse_data)
final_sparse_data_matrix <- Matrix(final_sparse_data_matrix, sparse=TRUE)

## Let's visualize.
head(final_sparse_data_matrix)
print(object.size(final_sparse_data_matrix),units="auto")