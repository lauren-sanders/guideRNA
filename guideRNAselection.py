#!/usr/bin/env python3
########################################################################
# File: guideRNAselection.py
# executable: 
# Purpose: Design guide RNAs for selected exons. 
#           gRNAs will be at three positions: mid-exon, 5' splice site, 3' splice site.
#           The 2 splice sites will have three gRNAs associated with them: 
#           1) gRNA with cutsite closest to splice site
#           2) gRNA with highest score near splice site
#           3) gRNA with next highest score near the splice site
#
#
#  Required to first use parse_crispr.py to parse a crispr file from the genome browser.
#  Input file: a file containing only the Associated_Exon_Coordinates column from a JuncBase file.
#  Output files: infile_3PrimeGuideRNAs.csv
#                infile_5PrimeGuideRNAs.csv
#                infile_MidExonGuideRNAs.csv
#          
# Author: Lauren M Sanders
# History: LMS 11/01/16 Created
########################################################################

import sys 
import argparse
import csv
import os
from argparse import RawTextHelpFormatter 

class CommandLine() :
    '''Implements a help option with program information, and an input file option.'''
    def __init__(self) :
        self.parser = argparse.ArgumentParser(description = "Gets guide RNA sequences for exons in an exon coordinate file.\n"
                                                            "Designed to take an input file containing the Associated_Exon_Coordinates column from a JuncBase file.\n"
                                                            "\n"
                                                            "Input file line format should be e.g. 1:1000-1100 with 1=chromosome, 1000=exonStart, 1100=exonEnd.\n"
                                                            "\n"
                                                            "Outputs 3 CSV files per exon coord infile; 1 file for 5' guides, 1 file for 3' guides, and 1 file for mid-exon guides.\n"
                                                            "\n"
                                                            "The 5' and 3' files each have 12 columns as follows:\n"
                                                            "\n"
                                                            "chromosome = exon chromosome number\n"
                                                            "exonStart = exon start coordinate\n"
                                                            "exonEnd = exon end coordinate\n"
                                                            "nearestGuide = guide nearest exon splice site regardless of score\n"
                                                            "nearestMIT = MIT score for the nearestGuide\n"
                                                            "nearestDoench = Doench score (percentage) for the nearestGuide\n"
                                                            "greenGuide = guide nearest splice site with Doench>60 and MIT>50\n"
                                                            "greenMIT = MIT score for greenGuide\n"
                                                            "greenDoench = Doench score (percentage) for the greenGuide\n"
                                                            "yellowGuide = guide nearest splice site with 30<Doench<60 and MIT>50\n"
                                                            "yellowMIT = MIT score for yellowGuide\n"
                                                            "yellowDoench = Doench score (percentage) for the yellowGuide\n"
                                                            "\n"
                                                            "The MidExon file has 6 columns: chromosome, exonStart, exonEnd, MidExonGuide, guideMIT and guideDoench.\n",
                                              formatter_class=RawTextHelpFormatter,
                                              add_help = True, 
                                              usage = 'guideRNAselection.py -f exonCoordinateFileExample.txt' ## fix this!!
                                                )
       
        self.parser.add_argument('-f', '--filename', action='store', help='File with exon coordinates')
        
        self.args = self.parser.parse_args()

class ExonFile : 
    '''
    Read a file containing one or more sets of exon coordinates for 
    which guide RNA design is desired.
            
    methods: 
    parseFile() : returns a list of all exon coordinates and chromosomes.
    
    '''
    def __init__(self, commandLine) : 
        self.cL = commandLine
        
    def parseFile(self):
        '''
        Add 'chr' before each chromosome number to match with crispr file.
        Convert each line to CSV format. Return a list with each exon as a separate element.
        '''
        allExons = []
        
        with open (self.cL.args.filename) as f:
        
            for exonLine in f:
                exonLine = 'chr'+exonLine.rstrip()
                exonLine = exonLine.replace(':',',')
                exonLine = exonLine.replace('-',',')
                
                allExons.append(exonLine)
            
        return allExons

class GuideRna : 
    '''
    Object is instantiated once for each chromosome represented in the exons file.
    Produces lists for each guide RNA location (5' of exon, 3' of exon, or mid-exon)
    and each subcategory within the location (closest to splice site, highest score,
    second highest score). 7 lists total.

    methods: 
    rangeLists() : Returns a list of exons in range of each guide RNA location (3 lists total)
    spliceSiteGuides() : Returns 3 lists each for 5' and 3' with guide RNAs closest to splice sites
    given the scoring constraints of each subcategory.
    midExonGuides() : Returns 1 list of guides nearest mid-exon location.
    
    '''
    
    def __init__(self, allExons, chr) :
        self.allExons = allExons
        self.chr = chr

    def rangeLists(self):
        '''
        Returns 3 lists: one for each guide RNA location. 
        for exon 5' splice site: all gRNAs with cutsite in range (200 nt 5' of exon start - midpoint of exon).
        for exon 3' splice site: all gRNAs with cutsite in range (midpoint of exon - 200 nt 3' of exon end)
        for mid-exon: all gRNAs with cutsite in range (exon start - exon end)
        '''
        exons_of_interest = []
        
        for exon in self.allExons:
            if exon.startswith(self.chr):    
                exons_of_interest.append(exon)
        
        guidesInFivePrimeRange = []
        guidesInThreePrimeRange = []
        guidesInMidRange = []
        
        with open(self.chr + 'crispr.txt', 'r') as h:
            
            # Iterate through each crispr file only once
            for crisprLine in h.readlines():
                crisprLine = crisprLine.strip()
                
                directionality = crisprLine.split(',')[3]
                guideStart = int(crisprLine.split(',')[1])
                guideSeq = crisprLine.split(',')[4]
                MITscore = crisprLine.split(',')[5]
                DoenchScore = crisprLine.split(',')[6].split(' ')[0].replace('%','')
                
                if directionality == '-' :
                        cutsite = int(guideStart+2)
                elif directionality == '+' :
                    cutsite = int(guideStart+16)
                
                # Iterate through exon file once per crispr line 
                for exon in exons_of_interest:
                    chromosome = exon.split(',')[0]
                    exonStart = int(exon.split(',')[1])
                    exonEnd = int(exon.split(',')[2])
                    
                    # Create a list of all gRNAs in the five prime cutsite range 
                    try:
                        if cutsite in range(int(exonStart-200),int((exonEnd-exonStart)/2)+exonStart) :
                            guidesInFivePrimeRange.append((cutsite,chromosome,exonStart,exonEnd,guideSeq,MITscore,DoenchScore))
                    except ValueError:
                        pass
                    
                    # Create a list of all gRNAs in the three prime cutsite range
                    try:
                        if cutsite in range(int(((exonEnd-exonStart)/2)+exonStart),int(exonEnd+200)) :
                            guidesInThreePrimeRange.append((cutsite,chromosome,exonStart,exonEnd,guideSeq,MITscore,DoenchScore))
                    except ValueError:
                        pass
                        
                    # Create a list of all gRNAs that fall within the exon 
                    try:
                        if cutsite in range(exonStart,exonEnd) :
                            guidesInMidRange.append((cutsite,chromosome,exonStart,exonEnd,guideSeq,MITscore,DoenchScore))
                    except ValueError:
                        pass
                        
            h.close()
                        
        return guidesInFivePrimeRange, guidesInThreePrimeRange, guidesInMidRange
        
    def spliceSiteGuides(self, guidesInFivePrimeRange, guidesInThreePrimeRange) :
        '''
        Returns 3 lists for each splice site location: 
        "nearest" = gRNA with splice site nearest splice site
        "green" = gRNA with MIT>50 and Doench>60 closest to splice site
        "yellow" = gRNA with MIT>50 and Doench>30 closest to splice site
        '''
        
        # Now we have three lists, holding ANY gRNA that exists in the 5', 3', or midexon range
        # The next step is to divide the 5' and 3' gRNAs up into the following categories: 
        # 1. "green" scores (MIT>50 and Doench>60)
        # 2. "yellow" scores (MIT>50 and Doench>30
        # 3. the gRNA that exists closest to the splice site itself, regardless of score
       
        # Create a list of five prime gRNAs that score in the "green" range
        greenFiveCandidates = []
        for guide in guidesInFivePrimeRange:
            if int(guide[5])>50 and int(guide[6])>60:
                greenFiveCandidates.append(guide)
        
        # Create a list of five prime gRNAs that score in the "yellow" range
        yellowFiveCandidates = []
        for guide in guidesInFivePrimeRange:
            if int(guide[5])>50 and int(guide[6])>30:
                yellowFiveCandidates.append(guide)
        
        # Create a list of three prime gRNAs that score in the "green" range
        greenThreeCandidates = []
        for guide in guidesInFivePrimeRange:
            if int(guide[5])>50 and int(guide[6])>60:
                greenThreeCandidates.append(guide)
        
        # Create a list of three prime gRNAs that score in the "yellow" range
        yellowThreeCandidates = []
        for guide in guidesInFivePrimeRange:
            if int(guide[5])>50 and int(guide[6])>30:
                yellowThreeCandidates.append(guide)
        
        def distance(guideRNA) :
            '''Returns distance between gRNA cutsite and exonStart'''
            
            cutsite = guideRNA[0]
            exonStart = guideRNA[2]
            
            distance = int(abs(exonStart-cutsite))
            return distance
            
        def nearestGuides(guideRangeList) :
            '''Returns guideRNAs nearest splice site in a given list of gRNAs'''
            
            guides_to_remove = []
            nearestGuides = []
            guides_left = True 
            guideRangeListLength = len(guideRangeList)
            while guides_left == True:
            
                same_exon = []
                distances = []
                
                try:
                    same_exon.append(guideRangeList[0])
                except IndexError: 
                    pass
                
                for j in range(1,len(guideRangeList)): 
                    if guideRangeList[0][2] == guideRangeList[j][2] : 
                        same_exon.append(guideRangeList[j])
                        
                for guide in same_exon: 
                    guideRangeList.remove(guide)
                    distances.append(distance(guide))
                try: 
                    nearestGuides.append(same_exon[distances.index(max(distances))])
                except ValueError:
                    pass
                    
                if len(guideRangeList) == 0 : 
                    guides_left = False
                    break
            
            return nearestGuides
            
        nearestFivePrime = nearestGuides(guidesInFivePrimeRange)
        greenFivePrime = nearestGuides(greenFiveCandidates)
        yellowFivePrime = nearestGuides(yellowFiveCandidates)
       
        nearestThreePrime = nearestGuides(guidesInThreePrimeRange)
        greenThreePrime = nearestGuides(greenThreeCandidates)
        yellowThreePrime = nearestGuides(yellowThreeCandidates)
       
       
        return nearestFivePrime, greenFivePrime, yellowFivePrime, nearestThreePrime, greenThreePrime, yellowThreePrime
        

    def midExonGuides(self, guideRangeList) :
        '''Returns a list of gRNAs, one nearest the midway point of each exon'''
        
        def distance(guideRNA) :
            '''Calculates the distance between a gRNA cutsite and the midway point of an exon'''
            
            cutsite = guideRNA[0]
            exonStart = guideRNA[2]
            exonEnd = guideRNA[3]
            
            distance = int(abs((((exonEnd-exonStart)/2)+exonStart)-cutsite))
            return distance
        
        
        guides_to_remove = []
        nearestGuides = []
        guides_left = True 
        while guides_left == True:
        
            same_exon = []
            distances = []
            
            try:
                same_exon.append(guideRangeList[0])
            except IndexError: 
                pass
            
            for j in range(1,len(guideRangeList)): 
                if guideRangeList[0][2] == guideRangeList[j][2] : 
                    same_exon.append(guideRangeList[j])
                    
            for guide in same_exon: 
                guideRangeList.remove(guide)
                distances.append(distance(guide))
            try: 
                nearestGuides.append(same_exon[distances.index(max(distances))])
            except ValueError:
                pass
                
            if len(guideRangeList) == 0 : 
                guides_left = False
                break
        
        return nearestGuides        

def main(cL=None):
    '''
    Instantiates all classes. 
    Passes each chromosome number to the GuideRNA class. 
    Retrieves lists of gRNAs for each subcategory of each position, 
    all chromosomes included. 
    Writes all gRNAs to 3 files, one for each position.
    '''
    cL = CommandLine ()
    
    allExons = ExonFile(cL).parseFile()
    
    chromList = ['chr1,','chr2,','chr3,','chr4,','chr5,','chr6,','chr7,','chr8,', 'chr9,','chr10,','chr11,','chr12,','chr13,','chr14,','chr15,','chr16,','chr17,','chr18,','chr19,','chr20,','chr21,','chr22,','chrX,','chrY,']

    # Lists for each of the 7 types; with ALL chromosomes in one list
    allnearestFivePrime = []
    allgreenFivePrime = []
    allyellowFivePrime = []
    allnearestThreePrime = []
    allgreenThreePrime = []
    allyellowThreePrime = []
    allMidExon = []
    
    for chrom in chromList:
        print("Finished guide RNAs for", chrom)
        
        newGuides = GuideRna(allExons, chrom)
    
        guidesInFivePrimeRange, guidesInThreePrimeRange, guidesInMidRange = newGuides.rangeLists()
    
        nearestFivePrime, greenFivePrime, yellowFivePrime, nearestThreePrime, greenThreePrime, yellowThreePrime = newGuides.spliceSiteGuides(guidesInFivePrimeRange, guidesInThreePrimeRange)
        
        midExonGuides = newGuides.midExonGuides(guidesInMidRange)
        
        # Combine each indiv chromosome list into a master list for each category
        for item in nearestFivePrime:
            allnearestFivePrime.append(item)
        for item in greenFivePrime:
            allgreenFivePrime.append(item)
        for item in yellowFivePrime:
            allyellowFivePrime.append(item)
        
        for item in nearestThreePrime:
            allnearestThreePrime.append(item)
        for item in greenThreePrime:
            allgreenThreePrime.append(item)
        for item in yellowThreePrime:
            allyellowThreePrime.append(item)
        
        for item in midExonGuides:
            allMidExon.append(item)
    
    header = "chromosome,exonStart,exonEnd,nearestGuide,nearestMIT,nearestDoench,greenGuide,greenMIT,greenDoench,yellowGuide,yellowMIT,yellowDoench"
    gap3 = ['-','-','-']
    gap6 = ['-','-','-','-','-','-']
    
    #########################
    # Write Five Prime File #
    #########################
    
    with open('allnearestFivePrime.csv', 'a') as f:
        for item in allnearestFivePrime:
            item = ','.join(map(str, item[1:]))
            f.writelines(item + '\n')
        f.close()
    
    with open('allnearestFivePrime.csv', 'r') as infile:
        with open('nearestAndgreenFivePrime.csv', 'a') as outfile:
            for item in allgreenFivePrime:
                written = False 
                for line in infile:
                    if int(line.split(',')[1]) == int(item[2]):
                        newline = line.rstrip('\n') + ',' + ','.join(map(str, item[4:]))
                        outfile.writelines(newline + '\n')
                        written = True 
                        break
                # if no previous entry exists for this exon, create one
                if written == False :
                    item = ','.join(map(str, item[1:]))
                    item = item.split(',')[:3] + gap3 + item.split(',')[3:]
                    outfile.writelines(','.join(map(str, item)) + '\n')
            outfile.close()
        infile.close()
                        
    with open('nearestAndgreenFivePrime.csv', 'r') as infile:
        with open(cL.args.filename.split('.')[0] + '_5PrimeGuideRNAs.csv', 'a') as outfile:
            outfile.writelines(header + "\n")
            for item in allyellowFivePrime:
                written = False
                for line in infile:
                    if int(line.split(',')[1]) == int(item[2]):
                        newline = line.rstrip('\n') + ',' + ','.join(map(str, item[4:]))
                        outfile.writelines(newline + '\n')
                        written = True
                        break
                if written == False : 
                    item = ','.join(map(str, item[1:]))
                    item = item.split(',')[:3] + gap6 + item.split(',')[3:]
                    outfile.writelines(','.join(map(str, item)) + '\n')
            outfile.close()
        infile.close()
    
    os.remove('allnearestFivePrime.csv')
    os.remove('nearestAndgreenFivePrime.csv')
    
    ##########################
    # Write Three Prime File #
    ##########################
    
    with open('allnearestThreePrime.csv', 'a') as f:
        for item in allnearestFivePrime:
            item = ','.join(map(str, item[1:]))
            f.writelines(item + '\n')
        f.close()
    
    with open('allnearestThreePrime.csv', 'r') as infile:
        with open('nearestAndgreenThreePrime.csv', 'a') as outfile:
            for item in allgreenThreePrime:
                written = False
                for line in infile:
                    if int(line.split(',')[1]) == int(item[2]):
                        newline = line.rstrip('\n') + ',' + ','.join(map(str, item[4:]))
                        outfile.writelines(newline + '\n')
                        written = True 
                        break
                # if no previous entry exists for this exon, create one
                if written == False :
                    item = ','.join(map(str, item[1:]))
                    item = item.split(',')[:3] + gap3 + item.split(',')[3:]
                    outfile.writelines(','.join(map(str, item)) + '\n')
            outfile.close()
        infile.close()
                        
    with open('nearestAndgreenThreePrime.csv', 'r') as infile:
        with open(cL.args.filename.split('.')[0] + '_3PrimeGuideRNAs.csv', 'a') as outfile:
            outfile.writelines(header + "\n")
            for item in allyellowThreePrime:
                written = False
                for line in infile:
                    if int(line.split(',')[1]) == int(item[2]):
                        newline = line.rstrip('\n') + ',' + ','.join(map(str, item[4:]))
                        outfile.writelines(newline + '\n')
                        written = True
                        break
                if written == False : 
                    item = ','.join(map(str, item[1:]))
                    item = item.split(',')[:3] + gap6 + item.split(',')[3:]
                    outfile.writelines(','.join(map(str, item)) + '\n')
            outfile.close()
        infile.close()
    
    os.remove('allnearestThreePrime.csv')
    os.remove('nearestAndgreenThreePrime.csv')
    
    #######################
    # Write Mid Exon File #
    #######################
    
    with open (cL.args.filename.split('.')[0] + '_MidExonGuideRNAs.csv', 'a') as m:
        m.writelines("chromosome,exonStart,exonEnd,midExonGuide,MIT,Doench")
        for item in allMidExon:
            item = ','.join(map(str, item[1:]))
            m.writelines(item + '\n')
        m.close()
        
        
    print("Finished designing guide RNAs for the given exons.")
    print("5' splice site, 3' splice site, and mid-exon gRNAs can be found in the following files:")
    print(cL.args.filename.split('.')[0] + '_5PrimeGuideRNAs.csv')
    print(cL.args.filename.split('.')[0] + '_3PrimeGuideRNAs.csv')
    print(cL.args.filename.split('.')[0] + '_MidExonGuideRNAs.csv')
    #######################
    # Write Skipped Exons #
    #######################
    
    with open(cL.args.filename, 'r') as f:
        exons = [part for part in [entry.strip().replace(':','-') for entry in f.readlines()[1:]]]
        exonstart = [part for part in [entry.split('-')[1] for entry in exons]]
    
    # Exons with no 5' guides
    missed_exons_five = []
    with open(cL.args.filename.split('.')[0] + '_5PrimeGuideRNAs.csv', 'r') as five:
        guides = [part for part in [entry.strip().split(',')[1] for entry in five.readlines()[1:]]]
        for exon in exonstart:
            if exon not in guides:
                missed_exons_five.append(exon)
                
    print("No 5' splice site guide RNAs were found for the following exons:")
    for miss in missed_exons_five:
        for exon in exons:
            if miss in exon:
                print(exon)
                
    # Exons with no 3' guides
    missed_exons_three = []
    with open (cL.args.filename.split('.')[0] + '_3PrimeGuideRNAs.csv', 'r') as three:
        guides = [part for part in [entry.strip().split(',')[1] for entry in three.readlines()[1:]]]
        for exon in exonstart:
            if exon not in guides:
                missed_exons_three.append(exon)
                
    print("No 3' splice site guide RNAs were found for the following exons:")
    for miss in missed_exons_three:
        for exon in exons:
            if miss in exon:
                print(exon)
                
    # Exons with no mid-exon guides
    missed_exons_mid = []
    with open (cL.args.filename.split('.')[0] + '_MidExonGuideRNAs.csv', 'r') as m:
        guides = [part for part in [entry.strip().split(',')[1] for entry in m.readlines()[1:]]]
        for exon in exonstart:
            if exon not in guides:
                missed_exons_mid.append(exon)
                
    print("No mid-exon guide RNAs were found for the following exons:")
    for miss in missed_exons_mid:
        for exon in exons:
            if miss in exon:
                print(exon)
    f.close()
    five.close()
    three.close()
    m.close()
    
if __name__ == "__main__":
    main();
    raise SystemExit
        
