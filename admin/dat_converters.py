#!/usr/bin/python

'''
Converts binary file formats used by maracluster
'''

import sys
import os
import csv

import struct
import string

def main(argv):
    usage = "Usage: python converters.py pvalBinToTsv|pvalTsvToBin|pvecBinToTsv|specBinToTsv|scanBinToTsv"
    if argv[0] == '-h':
        print(usage)
        return 0
    elif len(argv) < 3:
        print(usage)
        sys.exit("Not enough input parameters")
    
    inFile = argv[1]
    outFile = argv[2]
    if outFile.endswith("csv"):
        outDelimiter = ','
    else:
        outDelimiter = '\t'
    reader = csv.reader(open(inFile,'rb'), delimiter='\t')
    writer = csv.writer(open(outFile,'w'), delimiter=outDelimiter)
    
    if argv[0] == "pvalBinToTsv": # pvalue files, e.g. 601.dat.pvalues.dat
        convertPvalTripletsBinToTsv(inFile, writer)
    elif argv[0] == "pvalTsvToBin":
        convertPvalTripletsTsvToBin(reader, outFile)
    elif argv[0] == "specPrecMapToTsv":
        convertSpecPrecMapBinToTsv(inFile, writer)
    elif argv[0] == "pvecBinToTsv": # pvalue vector files, e.g. 601.dat.pvalue_vectors.dat
        convertPvecBinToTsv(inFile, writer)
    elif argv[0] == "specBinToTsv": # spectrum files, e.g. 358.dat
        convertSpecBinToTsv(inFile, writer)
    elif argv[0] == "scanBinToTsv": # scan info files, e.g. MaRaCluster.scan_info.dat
        convertScanBinToTsv(inFile, writer)
    else:
        print("Unknown mode:", argv[0])
        print(usage)
        
#########################################################
##  unsigned int fileIdx1, scannr1, fileIdx2, scannr2; ##
##  float pval;                                        ##
#########################################################

# python converters.py /media/storage/mergespec/data/batchcluster/Linfeng/601.dat.pvalues.dat /media/storage/mergespec/data/batchcluster/Linfeng/601.dat.pvalues.tsv

def convertPvalTripletsBinToTsv(inFile, writer):
    with open(inFile,'rb') as f:
        while True:
            a = f.read(20)
            if not a:
                break
            b = struct.unpack('=IIIIf', a)
            writer.writerow(b)


def convertPvalTripletsTsvToBin(reader, outFile):
    with open(outFile,'w') as w:
        for row in reader:
            row = [int(row[0]), int(row[1]), int(row[2]), int(row[3]), float(row[4])]
            b = struct.pack('=IIIIf', *row)
            w.write(b)

######################################
##  unsigned int fileIdx, scannr;   ##
##  float precMz, rTime, intensity; ##
##  float rtStart, rtEnd;           ##
##  int charge;                     ##
##  int fileIdx, featureIdx;        ##
######################################

# python converters.py specPrecMapToTsv /media/hdd/matthew/mergespec/data/iPRG2016/maraquant/iPRG2016/maracluster_extra_features/iPRG2016.spectrum_to_precursor_map.dat /media/hdd/matthew/mergespec/data/iPRG2016/maraquant/iPRG2016/maracluster_extra_features/iPRG2016.spectrum_to_precursor_map.tsv

def convertSpecPrecMapBinToTsv(inFile, writer):
    with open(inFile,'rb') as f:
        i = 0
        while True:
            a = f.read(40)
            if not a:
                break
            b = struct.unpack('=IIfffffiii', a)
            writer.writerow(b)
            i += 1
            
#####################################
##  double precMass;               ##
##  double precMz;                 ##
##  double retentionTime;          ##
##  double polyfit[6];             ##
##  short peakBins[40];            ##
##  short peakScores[40];          ##
##  int charge, queryCharge;       ##
##  ScanId scannr;                 ##
#####################################

# python converters.py pvecBinToTsv /media/storage/mergespec/data/batchcluster/Linfeng/601.dat.pvalue_vectors.dat /media/storage/mergespec/data/batchcluster/Linfeng/601.dat.pvalue_vectors.tsv

def convertPvecBinToTsv(inFile, writer):
    with open(inFile,'rb') as f:
        i = 0
        while True:
            a = f.read(248)
            if not a:
                break
            b = struct.unpack('=ddd6d40h40hiiII', a)
            writer.writerow(b)
            i += 1
            
#####################################
##  ScanId scannr;                 ##
##  unsigned int charge;           ##
##  float precMass, retentionTime; ##
##  short fragBins[100];           ##
#####################################

# python converters.py /media/storage/mergespec/data/batchcluster/Linfeng/601.dat /media/storage/mergespec/data/batchcluster/Linfeng/601.tsv

def convertSpecBinToTsv(inFile, writer):
    with open(inFile,'rb') as f:
        i = 0
        while True:
            a = f.read(100)
            if not a:
                break
            b = struct.unpack('=IIIff40h', a)
            writer.writerow(b)
            i += 1

def convertScanBinToTsv(inFile, writer):
    with open(inFile,'rb') as f:
        i = 0
        while True:
            a = f.read(16)
            if not a:
                break
            b = struct.unpack('=IIff', a)
            writer.writerow(b)
            i += 1

def countNumSpectraUnderNPeaks(inFile, N = 40):
    with open(inFile,'rb') as f:
        i = 0
        count = 0
        lastScannr = 0
        while True:
            a = f.read(220)
            if not a:
                break
            b = struct.unpack('=IIIff100h', a)
            firstZero = next((i for i, x in enumerate(b[4:]) if x == 0), 100)
            if firstZero < N and lastScannr != b[0]:
                lastScannr = b[0]
                #print b[0], firstZero
                count += 1
            #print b[val-1], b[val], val
            i += 1
        print(count)


def precMzFromPrecMass(pmass, z):
    return (float(pmass) + 1.00727646677 * (int(z) - 1)) / int(z)


def precMassFromPrecMz(pmz, z):
    return pmz * z - 1.00727646677 * (z - 1)


if __name__ == "__main__":
     main(sys.argv[1:])
