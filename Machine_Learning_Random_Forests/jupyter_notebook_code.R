mydata_enhancers <- read.csv("./GSM929090.bed", sep = '\t', header = FALSE)
mydata_nonenhancers <- read.csv("./GSM468792.bed", sep = '\t', header = FALSE)

mydata_enhancers <- mydata_enhancers[,c(1,2,3)]
mydata_enhancers$Class <- as.factor("Enhancer")
colnames(mydata_enhancers) <- c("Chrom", "Start", "End", "Class")

set.seed(007)
nonenhancer_indexes <- sample(1:nrow(mydata_nonenhancers),nrow(mydata_enhancers)) 
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

my_consolidated_dataset$Class <- as.factor(my_consolidated_dataset$Class)

chr1_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr1", ]
chr2_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr2", ]
chr3_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr3", ]
chr4_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr4", ]
chr5_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr5", ]
chr6_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr6", ]
chr7_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr7", ]
chr8_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr8", ]
chr9_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr9", ]
chr10_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr10", ]
chr11_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr11", ]
chr12_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr12", ]
chr13_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr13", ]
chr14_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr14", ]
chr15_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr15", ]
chr16_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr16", ]
chr17_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr17", ]
chr18_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr18", ]
chr19_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr19", ]
chr20_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr20", ]
chr21_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr21", ]
chr22_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chr22", ]
chrX_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chrX", ]
chrY_data<- my_consolidated_dataset[my_consolidated_dataset$Chrom == "chrY", ]

chr2_data$Start <- max(chr1_data$End)+as.double(chr2_data$Start)
chr2_data$End <- max(chr1_data$End)+as.double(chr2_data$End)
chr3_data$Start <- max(chr2_data$End)+as.double(chr3_data$Start)
chr3_data$End <- max(chr2_data$End)+as.double(chr3_data$End)
chr4_data$Start <- max(chr3_data$End)+as.double(chr4_data$Start)
chr4_data$End <- max(chr3_data$End)+as.double(chr4_data$End)
chr5_data$Start <- max(chr4_data$End)+as.double(chr5_data$Start)
chr5_data$End <- max(chr4_data$End)+as.double(chr5_data$End)
chr6_data$Start <- max(chr5_data$End)+as.double(chr6_data$Start)
chr6_data$End <- max(chr5_data$End)+as.double(chr6_data$End)
chr7_data$Start <- max(chr6_data$End)+as.double(chr7_data$Start)
chr7_data$End <- max(chr6_data$End)+as.double(chr7_data$End)
chr8_data$Start <- max(chr7_data$End)+as.double(chr8_data$Start)
chr8_data$End <- max(chr7_data$End)+as.double(chr8_data$End)
chr9_data$Start <- max(chr8_data$End)+as.double(chr9_data$Start)
chr9_data$End <- max(chr8_data$End)+as.double(chr9_data$End)
chr10_data$Start <- max(chr9_data$End)+as.double(chr10_data$Start)
chr10_data$End <- max(chr9_data$End)+as.double(chr10_data$End)
chr11_data$Start <- max(chr10_data$End)+as.double(chr11_data$Start)
chr11_data$End <- max(chr10_data$End)+as.double(chr11_data$End)
chr12_data$Start <- max(chr11_data$End)+as.double(chr12_data$Start)
chr12_data$End <- max(chr11_data$End)+as.double(chr12_data$End)
chr13_data$Start <- max(chr12_data$End)+as.double(chr13_data$Start)
chr13_data$End <- max(chr12_data$End)+as.double(chr13_data$End)
chr14_data$Start <- max(chr13_data$End)+as.double(chr14_data$Start)
chr14_data$End <- max(chr13_data$End)+as.double(chr14_data$End)
chr15_data$Start <- max(chr14_data$End)+as.double(chr15_data$Start)
chr15_data$End <- max(chr14_data$End)+as.double(chr15_data$End)
chr16_data$Start <- max(chr15_data$End)+as.double(chr16_data$Start)
chr16_data$End <- max(chr15_data$End)+as.double(chr16_data$End)
chr17_data$Start <- max(chr16_data$End)+as.double(chr17_data$Start)
chr17_data$End <- max(chr16_data$End)+as.double(chr17_data$End)
chr18_data$Start <- max(chr17_data$End)+as.double(chr18_data$Start)
chr18_data$End <- max(chr17_data$End)+as.double(chr18_data$End)
chr19_data$Start <- max(chr18_data$End)+as.double(chr19_data$Start)
chr19_data$End <- max(chr18_data$End)+as.double(chr19_data$End)
chr20_data$Start <- max(chr19_data$End)+as.double(chr20_data$Start)
chr20_data$End <- max(chr19_data$End)+as.double(chr20_data$End)
chr21_data$Start <- max(chr20_data$End)+as.double(chr21_data$Start)
chr21_data$End <- max(chr20_data$End)+as.double(chr21_data$End)
chr22_data$Start <- max(chr21_data$End)+as.double(chr22_data$Start)
chr22_data$End <- max(chr21_data$End)+as.double(chr22_data$End)
chrX_data$Start <- max(chr22_data$End)+as.double(chrX_data$Start)
chrX_data$End <- max(chr22_data$End)+as.double(chrX_data$End)
chrY_data$Start <- max(chrX_data$End)+as.double(chrY_data$Start)
chrY_data$End <- max(chrX_data$End)+as.double(chrY_data$End)

# interval_streamlining <- function(i){
#   for(i in chromosomes)
#   {
#     # eval(paste(eval(parse(text=i)),"_data"))$Start <- max(chr1_data$End)+eval(paste(eval(parse(text=i)),"_data"))$Start
#     # eval(paste(eval(parse(text=i)),"_data"))$End <- max(chr1_data$End)+eval(paste(eval(parse(text=i)),"_data"))$End
#   }
# }

my_consolidated_dataset <- my_consolidated_dataset[,-1]

