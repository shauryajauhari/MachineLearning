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


## Replacing NAs with 0s(zeros) in the reads' columns. Since the score data is not available, imputing empty cells with
## zero entries engenders mathematical convenience.

input_score_data$H3K27ac[is.na(input_score_data$H3K27ac)] <- 0 
input_score_data$H3K4me3[is.na(input_score_data$H3K4me3)] <- 0 
input_score_data$H3K4me2[is.na(input_score_data$H3K4me2)] <- 0 
input_score_data$H3K4me1[is.na(input_score_data$H3K4me1)] <- 0 

## Converting data to GRanges objects
positive_class <- positive_class[-1,]
positive_class_labels <- GRanges(seqnames = positive_class$V1, ranges = IRanges(start = positive_class$V2, 
                                                                       end = positive_class$V3))
mcols(positive_class_labels) <- DataFrame(class= "enhancer")


negative_class_labels <- GRanges(seqnames = negative_class$V1, ranges = IRanges(start = negative_class$V2, 
                                                                       end = negative_class$V3))
mcols(negative_class_labels) <- DataFrame(class= "non-enhancer")


input_score <- GRanges(seqnames = input_score_data$chrom, ranges = IRanges(start = input_score_data$start,
                                                                        end = input_score_data$end))
mcols(input_score) <- DataFrame(reads_h3k27ac = input_score_data$H3K27ac, reads_h3k4me3 = input_score_data$H3K4me3,
                                reads_h3k4me2 = input_score_data$H3K4me2, reads_h3k4me1 = input_score_data$H3K4me1)

## Performing merge to figure out the score and class matrix.
## Appending negative and positive class intervals together.

positive_and_negative_class_intervals <- merge(as.data.frame(positive_class_labels), as.data.frame(negative_class_labels), all=TRUE) ## positive and negative classes

## Exporting 'positive_and_negative_class_intervals' and 'input_score' as bed files to merge.
write.table(positive_and_negative_class_intervals, "./data/class_labels.bed", sep ='\t', quote = FALSE, row.names = FALSE)
write.table(input_score, "./data/score.bed", sep ='\t', quote = FALSE, row.names = FALSE)

## Processing for syntax and merging

system("awk '{if (NR!=1) {print}}' ./data/class_labels.bed > ./data/class_labels_header_removed.bed")
system("awk '{if (NR!=1) {print}}' ./data/score.bed > ./data/score_header_removed.bed")

system("bedtools intersect -wa -wb -a ./data/score_header_removed.bed -b ./data/class_labels_header_removed.bed > ./data/score_labels.bed")

## While merging the intervals for class labels (positive and negative) and RPKM normalized scores, there are occassions where one interval of the latter may
## encompass several intervals from the former. This is reflected in the 'bedtools intersect' operation.

## Importing the merged file
score_labels <- read.table("./data/score_labels.bed", sep = "\t", header = FALSE, stringsAsFactors=FALSE)

## Picking relevant columns
final_data <- score_labels[,c(6:9,15)]
colnames(final_data)=c("reads_h3k27ac","reads_h3k4me3","reads_h3k4me2","reads_h3k4me1","class")

## Saving this data
saveRDS(final_data,"./data/ep_data.rds")

## Sample data
final_data_sample <- final_data[sample(nrow(final_data), 10000), ]
saveRDS(final_data_sample,"./data/ep_data_sample.rds")

## Converting data frame to matrix and then to sparse matrix.
final_sparse_data_matrix <- data.matrix(final_sparse_data)
final_sparse_data_matrix <- Matrix(final_sparse_data_matrix, sparse=TRUE)

## Let's visualize.
head(final_sparse_data_matrix)
print(object.size(final_sparse_data_matrix),units="auto")


## Saving relevant files
save(positive_class_labels, tss_rev, p300_rev_finale, input_score, hg38_random, dhs_min, file = "relevant_files.RData")
