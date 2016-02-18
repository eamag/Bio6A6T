#! /usr/bin/env python
# -*- coding: utf-8 -*-
from Bio import SeqIO
from Levenshtein import distance
import os

NumOfDif = 2  # input("Enter Number of difference from 6A6T: ")+1
way = '/home/shared/STEMLOOPS/hg19/S15-30_L0-10_M5/'  # 'C:\моё\универ\прогаем\Bio6A6T\Bio\human'  # '/home/eamag/Bio'
cout = open('humanOut.txt', 'w')
path2 = '/home/magas/Bio/human/human_hg19_refseq'  # 'C:\моё\универ\прогаем\Bio6A6T\Bio\human\human_hg19_refseq'
x = "AAAAAATTTTTT"
cout.write('Name\tMotif\tCoordinate\tStrand\tName\tUpstream\tUpstream dist\tDownstream\tDownstream dist\n')
temp = ['name', 'chrom', '+', '0', '100000000000']
for d, dirs, files in os.walk(way):
    for names in files:
        if (names[-4:-1] + names[-1]) == '.fna':  # find all .fna files, parse it and find motif via levenshtein
            path = os.path.join(d, names)
            if os.path.exists(path):
                handle = open(path)
                for record in SeqIO.parse(handle, "fasta"):
                    w = record.seq
                    for i in range(0, len(str(w)) - 12):
                        if distance(x, str(w)[i:(i + 12)]) < NumOfDif:
                            chrom = names[0:-4]
                            handle2 = open(path2)
                            handle2.readline()
                            for line in handle2:
                                parts = line.rstrip().split('\t')
                                loc = [parts[1], parts[2], parts[3], parts[4], parts[5]]
                                name = loc[0]
                                chrNo = loc[1]
                                strand = loc[2]
                                start = loc[3]
                                end = loc[4]
                                if int(start) < i < int(end) and chrom == chrNo:  # in gene
                                    cout.write(chrNo + '\t' + str(w)[i:(i + 12)] + '\t' + str(i) + '\t' + strand +
                                               '\t' + name + '\tnone\tnone\tnone\tnone\n')
                                if int(temp[4]) < i < int(start) and chrom == chrNo:
                                    if temp[2] == '-' and strand == '-':  # strand is '-' both
                                        cout.write(chrNo + '\t' + str(w)[i:(i + 12)] + '\t' + str(i) + '\tnone\tnone' +
                                                   '\t' + temp[0] + '\t' + str(i - int(temp[4])) + '\t' + name + '\t' +
                                                   str(int(start) - i) + '\n')
                                    elif temp[2] == '+' and strand == '+':  # strand is '+' both
                                        cout.write(chrNo + '\t' + str(w)[i:(i + 12)] + '\t' + str(i) + '\tnone\tnone' +
                                                   '\t' + name + '\t' + str(int(start) - i) + '\t' + temp[0] + '\t' +
                                                   str(i - int(temp[4])) + '\n')
                                    elif temp[2] == '+' and strand == '-':  # first strand + and second - "down/down"
                                        cout.write(chrNo + '\t' + str(w)[i:(i + 12)] + '\t' + str(i) + '\tnone\tnone' +
                                                   '\tnone\tnone\t' + name + '/' + temp[0] + '\t' +
                                                   str(int(start) - i) + '/' + str(i - int(temp[4])) + '\n')
                                    elif temp[2] == '-' and strand == '+':  # first strand - and second + "up/up"
                                        cout.write(chrNo + '\t' + str(w)[i:(i + 12)] + '\t' + str(i) + '\tnone\tnone' +
                                                   '\tnone\tnone\t' + name + '/' + temp[0] + '\t' +
                                                   str(int(start) - i) + '/' + str(i - int(temp[4])) + '\tnone\tnone\n')
                                    else:
                                        print('error\n')
                                temp = loc
                            handle2.close()
                handle.close()
cout.close()
