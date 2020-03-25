multiBamSummary bins --bamfiles \
/Users/soumyajauhari/Desktop/Machine_Learning/Machine_Learning_Deep_Learning/data/H1_Cell_Line/H3K27ac/ENCFF663SAM.bam \
/Users/soumyajauhari/Desktop/Machine_Learning/Machine_Learning_Deep_Learning/data/H1_Cell_Line/H3K4me1/ENCFF441KOL.bam \
/Users/soumyajauhari/Desktop/Machine_Learning/Machine_Learning_Deep_Learning/data/H1_Cell_Line/H3K4me2/ENCFF799BDH.bam \
/Users/soumyajauhari/Desktop/Machine_Learning/Machine_Learning_Deep_Learning/data/H1_Cell_Line/H3K4me3/ENCFF340UJK.bam \
--binSize 2000 \
--labels H3K27ac H3K4me1 H3K4me2 H3K4me3 \
-out readCountsFeatures.npz \
--outRawCounts readCountsFeatures.tab
