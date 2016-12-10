####Purpose: Design guide RNAs for selected exons.<br />
####Author: Lauren Sanders
####Version 1: 12/9/16
<br /> 
This proram uses Max Haussler's guide RNA custom track on the UCSC Genome Browser to design guide RNAs for CRISPR screens targeting exons. <br /> 
The program produces guide RNAs targeting three positions on each exon: mid-exon, 5' splice site, 3' splice site.<br />
<br /> 
In order to provide flexibility regarding guide RNA quality score and guide RNA position, the 5' and 3'  splice sites will have three gRNAs associated with them: <br />
  - gRNA with cutsite closest to splice site<br />
  - gRNA with highest score within 200 bp of splice site<br />
  - gRNA with next highest score within 200 bp of splice site<br />

<br /> 
The mid-exon site will only have one gRNA, the one closest to the mid-exon point.<br />
<br />
Usage: <br />
<br />
Part 1: Crispr Files (ONE-TIME USE)<br />
a. download the bigBed file crispr.bb from http://hgdownload.cse.ucsc.edu/gbdb/hg19/crispr/<br />
b. perform the following Unix commands to download the bigBedToBed conversion tool, convert<br />
    the bigBed file to a Bed File, and remove all gRNAs whose sequence is not unique in the genome:<br /> 
  > wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed

  > chmod +x bigBedToBed
  
  > ./bigBedToBed crispr.bb crispr.bed
  
  > sed '/MIT Spec. Score: -1/d' ./crispr.bed | sed '/Sequence is not unique in genome/d' | cut -f 1,2,3,6,12,14,15,16 >crispr_parsed.bed
  
c. run this command (runtime ~ 3 hrs) to parse the crispr file into individual chromosome files: <br /> 
  > python3 parse_crispr.py <crispr_parsed.bed
d. once you have the individual chromosome crispr files, you do not need to repeat this part.<br />

Part 2: gRNA Design (Every time you need to design gRNAs for a new exon set)<br />
a. Retrieve only the Associated_Exon_Coordinates column from your JuncBase file, and remove the header.<br />
b. Remove the header using this command: 
    >sed -i '1d' infile
c. use this command the run the gRNA script: <br />
  > python3 guideRNAselection.py -f infile

Output:  3 files <br />
1) infile_5PrimeGuideRNAs.csv <br />
2) infile_3PrimeGuideRNAs.csv<br />
3) infile_MidExonGuideRNAs.csv<br />

The 5' and 3' files each have 12 columns as follows:<br />

1) chromosome = exon chromosome number<br />
2) exonStart = exon start coordinate<br />
3) exonEnd = exon end coordinate<br />
4) nearestGuide = guide nearest exon splice site regardless of score<br />
5) nearestMIT = MIT score for the nearestGuide<br />
6) nearestDoench = Doench score (percentage) for the nearestGuide<br />
7) greenGuide = guide nearest splice site with Doench>60 and MIT>50<br />
8) greenMIT = MIT score for greenGuide<br />
9) greenDoench = Doench score (percentage) for the greenGuide<br />
10) yellowGuide = guide nearest splice site with 30<Doench<60 and MIT>50<br />
11) yellowMIT = MIT score for yellowGuide<br />
12) yellowDoench = Doench score (percentage) for the yellowGuide<br />

The MidExon file has 6 columns: chromosome, exonStart, exonEnd, MidExonGuide, guideMIT and guideDoench.
