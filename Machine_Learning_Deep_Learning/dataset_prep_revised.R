## Preparing the Testing Data ##
## Importing the Histone Data ##

h3k27ac <- read.csv("./data/H1_Cell_Line/H3K27ac/ENCFF663SAM.bw", sep = '\t', header = FALSE)
colnames(h3k27ac) <- c("chrom","start","end","peaks")
h3k4me3 <- read.csv("./data/H1_Cell_Line/H3K4me3/ENCFF340UJK.bw", sep = '\t', header = FALSE)
colnames(h3k4me3) <- c("chrom","start","end","peaks")
h3k4me2 <- read.csv("./data/H1_Cell_Line/H3K4me2/ENCFF799BDH.bw", sep = '\t', header = FALSE)
colnames(h3k4me2) <- c("chrom","start","end","peaks")
h3k4me1 <- read.csv("./data/H1_Cell_Line/H3K4me1/ENCFF441KOL.bw", sep = '\t', header = FALSE)
colnames(h3k4me1) <- c("chrom","start","end","peaks")

## Importing the merged BEDGRAPH/BW file

merged_bw <- read.csv("./data/H1_Cell_Line/bedtools_Merge_2000.bedgraph", sep = '\t', header = FALSE)
#merged_bw <- merged_bw[-1,]
colnames(merged_bw) <- c("chrom", "start", "end", "peaks_h3k27ac", "peaks_h3k4me3", "peaks_h3k4me2", "peaks_h3k4me1")
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


## Merging Data on the basis of overlapping intervals

## Positive class data (labels)
positive_class

## Negative class (Labels)
negative_class

## Score Matrix (Input data)
input_score_data <- read.table("./data/H1_Cell_Line/bedtools_Merge_2000.bedgraph", sep = "\t", header = FALSE)
input_score_data<- as.data.frame(input_score_data[input_score_data$V1 %in% chromosomes, ])


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

intermatrix1 <- merge(positive_class_labels, negative_class_labels, all= TRUE) ## positive and negative classes ##
intermatrix2 <- merge(intermatrix1, input_score, all= TRUE) ## scores and classes ##

## Pulling back the data frame from the GenomicRanges format for (i) sorting (ii) removing NAs.

final_data <- data.frame(intermatrix2)
str(final_data)

## Sorting on the basis of first two columns, viz. seqnames, start.
final_data <- final_data[order(final_data$seqnames, final_data$start),]
buffer_final_data <- final_data

## Pruning rows involving NA terms.

for(i in 1:nrow(buffer_final_data)) ## all rows
{
  for(j in 1:length(buffer_final_data)) ## all columns
  {
    if (is.na(buffer_final_data[i,j])) ## if a cell has 'NA'
    {
      buffer_final_data <- buffer_final_data[-i,] ## remove that particular row
    }
  }
}
