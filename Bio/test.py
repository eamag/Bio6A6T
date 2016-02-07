#! /usr/bin/env python
# -*- coding: utf-8 -*-
from Bio import SeqIO
from Levenshtein import distance
import os

NumOfDif = 2  # input("Enter Number of difference from 6A6T: ")+1
way = 'C:\моё\универ\прогаем\Bio6A6T\Bio\human'  # '/home/eamag/Bio'
cout = open('humanOut.txt', 'w')
path2 = 'C:\моё\универ\прогаем\Bio6A6T\Bio\human\human_hg19_refseq'
x = "AAAAAATTTTTT"

handle2 = open(path2)
cout.write('Name\tMotif\tCoordinate\tStrand\tName\tUpstream\tUpstream dist\tDownstream\tDownstream dist\n')
#
#
# def SplitForLoc(handle):
#     line = handle.readline()
#
#     line = handle.readline()
#     parts = line.rstrip().split()
#     if len(parts) != 2 or parts[1].lower() != "proteins":
#         raise SyntaxError("Second line not recognised as an NCBI Protein Table (PTT file)")
#     line = handle.readline().strip()
#     if line.rstrip() != "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct":
#         raise SyntaxError("Third line not recognised as an NCBI Protein Table (PTT file)")
#     location = 0
#     strand = 1
#     product = 8
#
#     for line in handle:
#         parts = line.rstrip().split('\t')
#         loc = parts[location].split('..')
#         loc.append(parts[product])
#         loc.append(parts[strand])
#         yield loc


def split(handle):
    line = handle.readline()  # first line contents column names
    # if line.rstrip() != '#bin\tname\tchrom\tstrand\ttxStart\ttxEnd\tcdsStart\tcdsEnd\texonCount\texonStarts\texonEnds'
    #  \
    #                     '\tscore\tname2\tcdsStartStat\tcdsEndStat\texonFrames':
    #     raise SyntaxError('That is not a human_hg19_refseq')
    for line in handle:
        parts = line.rstrip().split('\t')
        loc = [parts[1], parts[2], parts[3], parts[4], parts[5]]
        yield loc

temp = ['name', 'chrom', '+', '0', '1000']
for d, dirs, files in os.walk(way):
    for names in files:
        if (names[-4:-1] + names[-1]) == '.fna':  # find all .fna files, parse it and find motif via levenshtein
            path = os.path.join(d, names)
            if os.path.exists(path):
                handle = open(path)
                for record in SeqIO.parse(handle, "fasta"):
                    w = record.seq
                    for i in range(0, len(str(w)) - 12):
                        # if str(w)[i] == 'N':
                        #     continue
                        if distance(x, str(w)[i:(i + 12)]) < NumOfDif:
                            chrNo = names[0:-4]

                            for loc in split(handle2):
                                name = loc[0]
                                chrom = loc[1]
                                strand = loc[2]
                                start = loc[3]
                                end = loc[4]  # temp is the same
                                if int(start) < i < int(end) and chrom == chrNo:  # in gene
                                    cout.write(chrNo + '\t' + str(w)[i:(i + 12)] + '\t' + str(i) + '\t' + strand +
                                               '/t' + name + '\tnone\tnone\tnone\tnone\n')
                                    print('hi')
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
                                                   '\t' + name + '/' + temp[0] + '\t' + str(int(start) - i) + '/' +
                                                   str(i - int(temp[4])) + '\tnone\tnone\n')
                                    else:
                                        print('error\n')
                                temp = loc

                            print(temp, i)
                            # for d, dirs, files in os.walk(way):  # find compared .ptt file and take info from there
                            #     for names2 in files:
                            #         if (names[0:-4] + '.ptt') == names2:
                            #             path2 = os.path.join(d, names2)
                            #
                            #
                            #
                            #                         elif loclist[k - 2][3] == '-' and loclist[k - 1][3] == '+':
                            #                             # first strand - and second + "upstream/upstream"
                            #
                            #                             cout.write(path3[-2] + '\t' + str(w)[i:(i + 12)] + '\t' +
                            #                                        str(i) + '\t' + 'none' + '\t' + 'none' + '\t' +
                            #                                       loclist[k - 1][2] + '/' + loclist[k - 2][2] + '\t' +
                            #                                        str(int(loclist[k - 1][0]) - i) + '/' +
                            #                                        str(i - int(loclist[k - 2][1])) + '\t' + 'none' +
                            #                                        '\t' + 'none' + '\n')
                            #
                            #                         else:
                            #                             print('ERROR.ERROR.ERROR.ERROR' + '\n')
                            #                 handle2.close()

                handle.close()
handle2.close()
cout.close()
