Code by Lauren Sanders

Purpose: Design guide RNAs for selected exons. 

gRNAs will be at three positions: mid-exon, 5' splice site, 3' splice site.
The 2 splice sites will have three gRNAs associated with them: 
1) gRNA with cutsite closest to splice site
2) gRNA with highest score near splice site
3) gRNA with next highest score near the splice site
The mid-exon site will only have one gRNA, the one closest to the mid-exon point.

Usage: 

Part 1: Crispr Files (ONE-TIME USE)
a. download the bigBed file crispr.bb from http://hgdownload.cse.ucsc.edu/gbdb/hg19/crispr/
b. perform the following Unix commands to download the bigBedToBed conversion tool, conver
    the bigBed file to a Bed File, and remove all gRNAs whose sequence is not unique in the genome: 
  > wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed
  > chmod +x bigBedToBed
  > ./bigBedToBed crispr.bb crispr.bed
  > sed '/MIT Spec. Score: -1/d' ./crispr.bed | sed '/Sequence is not unique in genome/d' | cut -f 1,2,3,6,12,14,15,16 > crispr_parsed.bed
c. run this command: python3 parse_crispr.py <crispr_parsed.bed (runtime ~ 3 hrs)
d. once you have the individual chromosome crispr files, you do not need to repeat this part.

Part 2: gRNA Design (Every time you need to design gRNAs for a new exon set)
a. Retrieve only the Associated_Exon_Coordinates column from your JuncBase file, and remove the header.
b. Remove the header using this command: sed -i '1d' infile
c. use this command the run the gRNA script: 
  > python3 guideRNAselection.py -f infile

Output:  3 files 
1) infile_5PrimeGuideRNAs.csv 
2) infile_3PrimeGuideRNAs.csv
3) infile_MidExonGuideRNAs.csv

The 5' and 3' files each have 12 columns as follows:

1) chromosome = exon chromosome number
2) exonStart = exon start coordinate
3) exonEnd = exon end coordinate
4) nearestGuide = guide nearest exon splice site regardless of score
5) nearestMIT = MIT score for the nearestGuide
6) nearestDoench = Doench score (percentage) for the nearestGuide
7) greenGuide = guide nearest splice site with Doench>60 and MIT>50
8) greenMIT = MIT score for greenGuide
9) greenDoench = Doench score (percentage) for the greenGuide
10) yellowGuide = guide nearest splice site with 30<Doench<60 and MIT>50
11) yellowMIT = MIT score for yellowGuide
12) yellowDoench = Doench score (percentage) for the yellowGuide

The MidExon file has 6 columns: chromosome, exonStart, exonEnd, MidExonGuide, guideMIT and guideDoench.
