## Load the data files of respecive enhancers and non-enhancers.
##p300 binding sites
mydata_enhancers <- read.csv("./props/GSM929090_uniq.bed.gz", sep = '\t', header = FALSE)

## H3K27ac binding sites
mydata_nonenhancers <- read.csv("./props/GSM2341643.tar.xz", sep = '\t', header = FALSE)


## Removing the first row in each data frame. It is the file name and is not required.

mydata_enhancers <- mydata_enhancers[-1,]
mydata_nonenhancers <- mydata_nonenhancers[-1,]


## Scaling down the enhancers relative to the non-enhancers.
set.seed(1)
enhancer_indexes <- sample(1:nrow(mydata_enhancers),(1/2*nrow(mydata_nonenhancers))) 
mydata_enhancers <- mydata_enhancers[enhancer_indexes,]

## Renaming the columns and adding another column to represent the class- enhancer/ non-enhancer.

mydata_enhancers <- mydata_enhancers[,c(1,2,3)]
mydata_enhancers$Class <- as.factor("Enhancer")
colnames(mydata_enhancers) <- c("Chrom", "Start", "End", "Class")

mydata_nonenhancers <- mydata_nonenhancers[,c(1,2,3)]
mydata_nonenhancers$Class <- "Non-Enhancer"
colnames(mydata_nonenhancers) <- c("Chrom", "Start", "End", "Class")

## Sorting on the basis of chromosome number and start location.

mydata_enhancers_sorted_chrom_names <- mydata_enhancers[with(mydata_enhancers, order(Chrom, Start)), ]
mydata_nonenhancers_sorted_chrom_names <- mydata_nonenhancers[with(mydata_nonenhancers, order(Chrom, Start)), ]

## Merging both datasets to a consolidated dataset.

my_consolidated_dataset <- rbind(mydata_enhancers_sorted_chrom_names, mydata_nonenhancers_sorted_chrom_names)
my_consolidated_dataset <- my_consolidated_dataset[with(my_consolidated_dataset, order(Chrom, Start)), ]

## Identifying examples of standard chromosomes only and filtering the residuals.

chromosomes <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20","chr21", "chr22", "chrX", "chrY")
my_consolidated_dataset<- as.data.frame(my_consolidated_dataset[my_consolidated_dataset$Chrom %in% chromosomes, ])

## Visualization and output class structure.
head(my_consolidated_dataset)
table(my_consolidated_dataset$Class)

## Introducing new columns representing enhancer classes as p300 elements and H3K27ac sites as non-enhancers. 
my_consolidated_dataset$p300 <- as.factor(ifelse(my_consolidated_dataset$Class =="Enhancer", "Yes", "No"))
my_consolidated_dataset$H3K27ac <- as.factor(ifelse(my_consolidated_dataset$Class =="Non-Enhancer", "Yes", "No"))


## Since neural networks work on continuous/ integer data, we must convert factor variables to numeric.
my_consolidated_dataset$p300 <- as.integer(my_consolidated_dataset$p300)
my_consolidated_dataset$H3K27ac <- as.numeric(my_consolidated_dataset$H3K27ac)

## Creating a surrogate feature/ variable for chrom, start, and end features.
my_consolidated_dataset$Bin <- as.character(1:nrow(my_consolidated_dataset))




### Superimposing Data ###

library(GenomicRanges)
dhs <- read.csv("./data/41009_peaks.bed", sep = '\t', header = FALSE)
ep300 <- read.csv("./data/41271_peaks.bed", sep = '\t', header = FALSE)
histm <- read.csv("./data/41276_peaks.bed", sep = '\t', header = FALSE)

## Selecting the first 3 columns ##

dhs <- dhs[,c(1,2,3)]
colnames(dhs) <- c("chrom","start","end")
ep300 <- ep300[,c(1,2,3)]
colnames(ep300) <- c("chrom","start","end")
histm <- histm[,c(1,2,3)]
colnames(histm) <- c("chrom","start","end")

## Filtering out irrelavant/ random chromosome entries ##

chromosomes <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20","chr21", "chr22", "chrX", "chrY")
dhs <- as.data.frame(dhs[dhs$chrom %in% chromosomes, ])
ep300 <- as.data.frame(ep300[ep300$chrom %in% chromosomes, ])
histm <- as.data.frame(histm[histm$chrom %in% chromosomes, ])


## Ordering on the basis of chromosomes ##

dhs <- dhs[with(dhs, order(chrom, start)), ]
ep300 <- ep300[with(ep300, order(chrom, start)), ]
histm <- histm[with(histm, order(chrom, start)), ]

## Converting the dataframes as GRanges object ##

dhs <- GRanges(dhs)
ep300 <- GRanges(ep300)
histm <- GRanges(histm)

## Finding commonalities ##

in1 <- merge(dhs,ep300)
in2 <- merge(in1,histm)
