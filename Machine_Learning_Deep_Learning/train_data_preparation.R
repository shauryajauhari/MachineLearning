## Importing relevant data ##
## Positive Class Labels ##

dhs <- read.csv("./data/H1_Cell_Line/GSM878621_H1_DNase_sorted.bed", sep = '\t', header = FALSE)
p300 <- read.csv("./data/H1_Cell_Line/GSM831036_H1_P300_sorted.bed", sep = '\t', header = FALSE)

## We see that the chromosome names in p300 data are just numbers. Let's add "chr" in the beginning for consistency. 
p300$V1 <- paste0("chr",p300$V1)

## Selecting useful columns: chrom, start, end.
p300 <- p300[,c(1,2,3)]
dhs <- dhs[,c(1,2,3)]

## Sifting standard chromosomes.
colnames(dhs) <- c("chrom","start","end")
colnames(p300) <- c("chrom","start","end")

chromosomes <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20","chr21", "chr22", "chrX", "chrY")
dhs <- as.data.frame(dhs[dhs$chrom %in% chromosomes, ])
p300 <- as.data.frame(p300[p300$chrom %in% chromosomes, ])



## Saving files

write.table(p300,"./data/H1_Cell_Line/h1_p300.bed", sep="\t", row.names=FALSE, quote = FALSE)
write.table(dhs,"./data/H1_Cell_Line/h1_dhs.bed", sep="\t", row.names=FALSE, quote = FALSE)

## Add class to the data: "enhancer"

h1_p300_dhs_intersect <- read.table("./data/H1_Cell_Line/h1_p300_dhs_intersect.bed", sep = "\t")
h1_p300_dhs_intersect$V4 <- "enhancer"
write.table(h1_p300_dhs_intersect,"./data/H1_Cell_Line/h1_p300_dhs_intersect_class.bed",sep="\t",row.names=FALSE, quote = FALSE)

#################################################################################################################

## Negative Class Labels ##

tss_sites <- read.table("./data/H1_Cell_Line/TSS_Indices_Human_Genome.txt", sep = "\t", header = TRUE)
tss_sites$Chromosome.scaffold.name <- paste0("chr", tss_sites$Chromosome.scaffold.name)
tss_sites <- as.data.frame(tss_sites[tss_sites$Chromosome.scaffold.name  %in% chromosomes, ])
tss_sites <- tss_sites[order(tss_sites$Chromosome.scaffold.name),]

## Export TSS sites to create an overlap with DHS sites

write.table(tss_sites,"./data/H1_Cell_Line/tss_sites.bed", sep="\t", row.names=FALSE, quote = FALSE)


## Shifting intervals forward by 3000 bp
tss_sites_test <- tss_sites
tss_sites_test$Transcript.start..bp. <- tss_sites_test$Transcript.start..bp.+3000
tss_sites_test$Transcript.end..bp. <- tss_sites_test$Transcript.end..bp.+3000

## Sourcing chromosomal lengths in the human genome for generating random tracks

hg19_chrom_sizes <- read.table(url("https://genome.ucsc.edu/goldenPath/help/hg19.chrom.sizes"), sep = "\t", header = FALSE, col.names = c("chrom", "size"))
hg19_chrom_sizes <- as.data.frame(hg19_chrom_sizes[hg19_chrom_sizes$chrom %in% chromosomes, ])

## Saving file for generating random tracks via 'bedtools random' function ##
write.table(hg19_chrom_sizes,"./data/H1_Cell_Line/hg19.genome", sep="\t", row.names=FALSE, quote = FALSE)

## Recalling

hg19_random_tracks <- read.table("./data/H1_Cell_Line/hg19_random_tracks.bed", sep = "\t", header = FALSE)
hg19_random_tracks <- hg19_random_tracks[,c(1,2,3)]
hg19_random_tracks_ordered <- hg19_random_tracks[order(hg19_random_tracks[,1],hg19_random_tracks[,2]),]

write.table(hg19_random_tracks_ordered,"./data/H1_Cell_Line/hg19_random_tracks_sorted_required.bed", sep="\t", row.names=FALSE, quote = FALSE)

## Choosing random sites distal to TSS

## The strategy is to combine random tracks distal to the TSS or p300 binding sites. Now, let us create a 
## combination of TSS and p300 binding sites and then subtract the random sites from these, thus giving us
## the residuals.

tss_final <- read.table("./data/H1_Cell_Line/tss_sites_header_removed.bed", sep = "\t", header = FALSE)
p300_final <- read.table("./data/H1_Cell_Line/h1_p300_merged.bed", sep = "\t", header = FALSE)

## Combine the intervals and not 'merge' them
tss_or_p300 <- rbind(tss_final,p300_final)

## sort on the basis of chromosome names.
tss_or_p300 <- tss_or_p300[order(tss_or_p300[,1]),]

## Output file.
write.table(tss_or_p300,"./data/H1_Cell_Line/tss_or_p300.bed", sep="\t", row.names=FALSE, quote = FALSE)


## Import resultant files from intersection.
negative_class <- read.table("./data/H1_Cell_Line/negative_class.bed", sep = "\t", header = FALSE)
negative_class$V4 <- "Non-Enhancer"
write.table(negative_class,"./data/H1_Cell_Line/negative_class.bed", sep="\t", row.names=FALSE, quote = FALSE)

positive_class <- read.table("./data/H1_Cell_Line/h1_p300_dhs_intersect_class_header_removed.bed", sep = "\t", header = FALSE)
