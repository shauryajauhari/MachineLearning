## Importing relevant data ##

dhs <- read.csv("./data/H1_Cell_Line/GSM878621_H1_DNase.bed.gz", sep = '\t', header = FALSE)
p300 <- read.csv("./data/H1_Cell_Line/GSM831036_H1_P300.bed.gz", sep = '\t', header = FALSE)

## We see that the chromosome names in p300 data are just numbers. Let's add "chr" in the beginning for consistency. 
p300$V1 <- paste0("chr",p300$V1)

## Selecting useful columns: chrom, start, end, peaks.
p300 <- p300[,c(1,2,3,5)]
dhs <- dhs[,c(1,2,3,5)]

## Sifting standard chromosomes.
colnames(dhs) <- c("chrom","start","end","peaks")
colnames(p300) <- c("chrom","start","end","peaks")

dhs <- as.data.frame(dhs[dhs$chrom %in% chromosomes, ])
p300 <- as.data.frame(p300[p300$chrom %in% chromosomes, ])

