#!/usr/bin/env python3
########################################################################
# File: parse_crispr.py
# Purpose: Convery a bigBed file from Max Haussler's custom Crispr guideRNA
#          genome browser track into subfiles corresponding to chromosome number.
#
#
#
#
#

# Author: Lauren M Sanders
# History: LMS 11/01/16 Created
########################################################################


import sys 
import argparse
import csv
from argparse import RawTextHelpFormatter 

class CommandLine() :
    '''Implements a help option with program information.'''
    def __init__(self) :
        self.parser = argparse.ArgumentParser(description = "Converts a bigBed file from Max Haeussler's custom Crispr guideRNA genome browser track\n"
                                                            "into subfiles corresponding to chromosome number, for downstream useability.\n"
                                                            "\n"
                                                            "First, download crispr.bb from http://hgdownload.cse.ucsc.edu/gbdb/hg19/crispr/ \n"
                                                            "Then perform the following linux commands:\n"
                                                            " \n"
                                                            "$ wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/bigBedToBed \n"
                                                            "$ chmod +x bigBedToBed \n"
                                                            "$ ./bigBedToBed crispr.bb crispr.bed\n"
                                                            "$ sed '/MIT Spec. Score: -1/d' ./crispr.bed | sed '/Sequence is not unique in genome/d' | cut -f 1,2,3,6,12,14,15,16 > crispr_parsed.bed \n"
                                                            " \n"
                                                            "Then execute this script on crispr_parsed.bed (runtime ~ 3 hrs) \n",
                                              formatter_class=RawTextHelpFormatter,
                                              add_help = True, 
                                              usage = 'crisprParse.py <crispr_parsed.bed'
                                                )
               
        self.args = self.parser.parse_args()
        
class CrisprReader : 
    '''
    to do
    '''
    
    def __init__ (self, infile):
        self.infile = infile 
        
    def readcrispr(self):
        '''Reads input file and writes each chromosome to a different file.'''
                        
        for line in self.infile:
            line = line.replace('\t',',').strip()
            with open(line.split(',')[0] + ',' + 'crispr.txt', 'a') as h:
                h.write(line+'\n')
                
                
        h.close()
        
        return 'Done parsing CRISPR guide RNA file'
                
def main(cL=None):
    '''Reads in a file from standard input. Calls the readcrispr function.'''
    cL = CommandLine ()
    
    crisprFiles = CrisprReader(sys.stdin)        

    makeFiles = crisprFiles.readcrispr()
    print(makeFiles)
    
if __name__ == "__main__":
    main();
    raise SystemExit
        
