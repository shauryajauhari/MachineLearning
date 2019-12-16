## Importing relevant data ##
## Positive Class Labels ##

################################################ REVISION ##########################################################################
# INSTEAD OF DOWNLOADING THE BED FILES AND USING THEM DIRECTLY, LETS CONVERT THEM TO BAM AND MODIFY THEIR INTERVALS  ###############
# OR ######  USE THE FOLLOWING CHUNK ###############################################################################################
####################################################################################################################################

# Downloading the relevant files #
# DHS data #
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM878nnn/GSM878621/suppl/GSM878621%5FUW%2EH1%2EChromatinAccessibility%2EDS19100%2Ebed%2Egz", "./data/H1_Cell_Line/GSM878621_H1_DNase.bed.gz")
system("gunzip ./data/H1_Cell_Line/GSM878621_H1_DNase.bed.gz")

# P300 data #
download.file("ftp://ftp.ncbi.nlm.nih.gov/geo/samples/GSM831nnn/GSM831036/suppl/GSM831036%5F635EFAAXX%2E4%2EH1%5FP300%2Ealn%2Ebed%2Egz", "./data/H1_Cell_Line/GSM831036_H1_P300.bed.gz")
system("gunzip ./data/H1_Cell_Line/GSM831036_H1_P300.bed.gz")


# Sorting the data #
system("sortBed -i ./data/H1_Cell_Line/GSM878621_H1_DNase.bed > ./data/H1_Cell_Line/GSM878621_H1_DNase_sorted.bed")
system("sortBed -i ./data/H1_Cell_Line/GSM831036_H1_P300.bed > ./data/H1_Cell_Line/GSM831036_H1_P300_sorted.bed")

# Working on DHS data #

library(rtracklayer)
dhs_rev <- import.bed("./data/H1_Cell_Line/GSM878621_H1_DNase_sorted.bed")
chromosomes <- c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14",
                 "chr15", "chr16", "chr17", "chr18", "chr19", "chr20","chr21", "chr22", "chrX", "chrY")
dhs_rev_chr <- dhs_rev[seqnames(dhs_rev) %in% chromosomes]

library(GenomicRanges)
dhs_rev_chr_resize <- resize(dhs_rev_chr, width = 2000)
dhs_final <- as.data.frame(dhs_rev_chr_resize)
export.bed(dhs_final[,1:3], "dhs_final.bed")

# for data integrity
system("awk '$2>0 && $3>$2 {print $1 "\t" $2 "\t" $3}' dhs_final.bed > dhs_finale.bed")
system("sortBed -i dhs_finale.bed > dhs_finale_sorted.bed")

# Working on P300 data #

p300_rev <- import.bed("./data/H1_Cell_Line/GSM831036_H1_P300_sorted.bed")
BiocManager::install("diffloop")
library(diffloop)
p300_rev <- addchr(p300_rev) # Add prefix "chr" to the seqnames #
p300_rev_chr <- p300_rev[seqnames(p300_rev) %in% chromosomes]

p300_rev_chr_resize <- resize(p300_rev_chr, width = 2000)
p300_rev_final <- as.data.frame(p300_rev_chr_resize)
p300_rev_final <- unique(p300_rev_final[,c("seqnames","start","end")]) # removing redundant entries #
export.bed(p300_rev_final, "p300_rev_final.bed")

system("sortBed -i p300_rev_final.bed > p300_rev_final_sorted.bed")
system("awk '{print $1 "\t" $2 "\t" $3}' p300_rev_final_sorted.bed > p300_rev_finale.bed")

## Now that we have the p300 bindings genomewide, we shall filter them for being distant to the promoter/TSS regions. Note that we are 
## solely interested in the ones that are distal to the TSS; distal as in non-proximal! 
## GRanges objects are workabkle for set operations.

p300_rev_finale <- import.bed("p300_rev_finale.bed")

## Let's pull the TSS sites

tss_sites <- read.table("./data/H1_Cell_Line/TSS_Indices_Human_Genome.txt", sep = "\t", header = TRUE)
tss_sites$Chromosome.scaffold.name <- paste0("chr", tss_sites$Chromosome.scaffold.name)
tss_sites <- as.data.frame(tss_sites[tss_sites$Chromosome.scaffold.name  %in% chromosomes, ])
tss_sites <- tss_sites[order(tss_sites$Chromosome.scaffold.name),]
colnames(tss_sites) <- c("chrom", "start", "end")
tss_rev <- GRanges(tss_sites)

## We are going to shove off TSS sites that intersect with p300 binding sites. This will give us the distal p300 sites.

p300_nonTSS <- setdiff(p300_rev_finale,tss_rev)
export.bed(p300_nonTSS, "p300_nonTSS.bed")
system("awk '{print $1 "\t" $2 "\t" $3}' p300_nonTSS.bed > p300_nonTSS_final.bed")

## Finally, saving as a intersection of distal p300 and DHS sites

system("intersectBed -a dhs_finale_sorted.bed -b p300_nonTSS_final.bed > h1_distalp300_dhs_intersect.bed")


## Add class to the data: "enhancer"

h1_distalp300_dhs_intersect <- read.table("h1_distalp300_dhs_intersect.bed", sep = "\t")
h1_distalp300_dhs_intersect$V4 <- "Enhancer"
write.table(h1_distalp300_dhs_intersect,"h1_distalp300_dhs_intersect_class.bed",sep="\t",row.names=FALSE, quote = FALSE)

#################################################################################################################

## Negative Class Labels ##

# Here, we intend to find the intersection of TSS with DHS sites.
dhs_min <- import.bed("dhs_finale_sorted.bed")
tss_dhs <- intersect(tss_rev, dhs_min)


## The other aspect for mapping negative class labels (non-enhancers) entails sourcing chromosomal lengths in the human genome for generating random tracks

hg38_chrom_sizes <- read.table(url("https://genome.ucsc.edu/goldenPath/help/hg38.chrom.sizes"), sep = "\t", header = FALSE, col.names = c("chrom", "size"))
hg38_chrom_sizes <- as.data.frame(hg38_chrom_sizes[hg38_chrom_sizes$chrom %in% chromosomes, ])

## Saving file for generating random tracks via 'bedtools random' function ##
write.table(hg38_chrom_sizes,"./data/H1_Cell_Line/hg38.genome", sep="\t", row.names=FALSE, quote = FALSE)

## Generating random tracks

system("bedtools random -l 2000 -g ./data/H1_Cell_Line/hg38.genome > ./data/H1_Cell_Line/hg38_random_tracks.bed")

## Recalling

hg38_random_tracks <- read.table("./data/H1_Cell_Line/hg38_random_tracks.bed", sep = "\t", header = FALSE)
hg38_random_tracks <- hg38_random_tracks[,c(1,2,3)]
hg38_random_tracks_ordered <- hg38_random_tracks[order(hg38_random_tracks[,1],hg38_random_tracks[,2]),]

write.table(hg38_random_tracks_ordered,"./data/H1_Cell_Line/hg38_random_tracks_sorted_required.bed", sep="\t", row.names=FALSE, quote = FALSE)
colnames(hg38_random_tracks_ordered) <- c("chrom", "start", "end")
hg38_random <- GRanges(hg38_random_tracks_ordered)

## Choosing random sites distal to TSS and p300 bindings

## The strategy is to combine random tracks distal to the TSS or p300 binding sites. Now, let us create a 
## combination of TSS and p300 binding sites and then subtract the random sites from these, thus giving us
## the residuals.


## Combine the intervals and not 'merge' them
p300_tss <- union(p300_rev_finale, tss_rev)


## Inferring random tracks distal to p300 bindings and TSS sites.
random_nonp300TSS <- setdiff(hg38_random,p300_tss)


## Merging the aforementioned random tracks and the TSS, DHS intersection.
random_tss_dhs <- union(random_nonp300TSS, tss_dhs)

## Output file.
write.table(random_tss_dhs,"./data/H1_Cell_Line/random_tss_dhs.bed", sep="\t", row.names=FALSE, quote = FALSE)
system("awk '{if (NR!=1) {print}}' ./data/H1_Cell_Line/random_tss_dhs.bed > ./data/H1_Cell_Line/random_tss_dhs_header_removed.bed")

## Import resultant files from intersection.
negative_class <- read.table("./data/H1_Cell_Line/random_tss_dhs_header_removed.bed", sep = "\t", header = FALSE)
negative_class$V4 <- "Non-Enhancer"
negative_class$V5 <- c() # removing strand information
write.table(negative_class,"./data/H1_Cell_Line/negative_class.bed", sep="\t", row.names=FALSE, quote = FALSE)

positive_class <- read.table("h1_distalp300_dhs_intersect_class.bed", sep = "\t", header = FALSE)

## Saving auxilliary informaiton
save(p300_rev, dhs_rev, file="aux.RData")

