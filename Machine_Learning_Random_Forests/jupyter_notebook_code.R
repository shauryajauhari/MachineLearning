mydata_enhancers <- read.csv("./GSM929090_uniq.bed.gz", sep = '\t', header = FALSE)
mydata_nonenhancers <- read.csv("./GSM468792_uniq.bed.gz", sep = '\t', header = FALSE)

mydata_enhancers <- mydata_enhancers[-1,]
mydata_nonenhancers <- mydata_nonenhancers[-1,]

mydata_enhancers <- mydata_enhancers[,c(1,2,3)]
mydata_enhancers$Class <- as.factor("Enhancer")
colnames(mydata_enhancers) <- c("Chrom", "Start", "End", "Class")

set.seed(1)
nonenhancer_indexes <- sample(1:nrow(mydata_nonenhancers),(3/2*nrow(mydata_enhancers))) 
mydata_nonenhancers <- mydata_nonenhancers[nonenhancer_indexes,]

#mydata_nonenhancers <- mydata_nonenhancers[mydata_nonenhancers$V5 > 0,]
mydata_nonenhancers <- mydata_nonenhancers[,c(1,2,3)]
mydata_nonenhancers$Class <- "Non-Enhancer"
colnames(mydata_nonenhancers) <- c("Chrom", "Start", "End", "Class")

mydata_enhancers_sorted_chrom_names <- mydata_enhancers[with(mydata_enhancers, order(Chrom, Start)), ]
mydata_nonenhancers_sorted_chrom_names <- mydata_nonenhancers[with(mydata_nonenhancers, order(Chrom, Start)), ]


my_consolidated_dataset <- rbind(mydata_enhancers_sorted_chrom_names, mydata_nonenhancers_sorted_chrom_names)
my_consolidated_dataset <- my_consolidated_dataset[with(my_consolidated_dataset, order(Chrom, Start)), ]

chromosomes <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20","chr21", "chr22", "chrX", "chrY")
my_consolidated_dataset<- as.data.frame(my_consolidated_dataset[my_consolidated_dataset$Chrom %in% chromosomes, ])
head(my_consolidated_dataset)

set.seed(2)
my_consolidated_dataset$Class <- factor(my_consolidated_dataset$Class)
chr1_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr1", ]
chr1_data_sample_indexes <- sample(1:nrow(chr1_data), 10000) 
chr1_data_sample <- chr1_data[chr1_data_sample_indexes,]


## Creating bins

create_bins <- function(chr1_data_sample){
  chr1_data_sample$bin1 <- cbind(bin1 <- chr1_data_sample$Start+100)
}
