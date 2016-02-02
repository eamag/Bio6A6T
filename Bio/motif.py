#! /usr/bin/env python
# -*- coding: utf-8 -*- 
from Bio import SeqIO
from Levenshtein import distance
import os

NumOfDif = 2  # input("Enter Number of difference from 6A6T: ")+1
way = '/home/shared'  # '/home/eamag/Bio'
cout = open('out.txt', 'w')
cout.write(
    'Name\tMotif\tCoordinate\tStrand\tName\tUpstream\tUpstream dist\tDownstream\tDownstream dist\n')


def SplitForLoc(handle):
    line = handle.readline()

    line = handle.readline()
    parts = line.rstrip().split()
    if len(parts) != 2 or parts[1].lower() != "proteins":
        raise SyntaxError("Second line not recognised as an NCBI Protein Table (PTT file)")
    line = handle.readline().strip()
    if line.rstrip() != "Location\tStrand\tLength\tPID\tGene\tSynonym\tCode\tCOG\tProduct":
        raise SyntaxError("Third line not recognised as an NCBI Protein Table (PTT file)")
    location = 0
    strand = 1
    product = 8

    mapping = {
        location: 'location',
        strand: 'strand',
        product: 'product',
    }

    for line in handle:
        parts = line.rstrip().split('\t')
        loc = parts[location].split('..')
        loc.append(parts[product])
        loc.append(parts[strand])
        yield loc


for d, dirs, files in os.walk(way):
    for names in files:
        if (names[-4:-1] + names[-1]) == '.fna':  # find all .fna files, parse it and find motif via levenshtein

            path = os.path.join(d, names)

            if os.path.exists(path):
                handle = open(path)
                for record in SeqIO.parse(handle, "fasta"):
                    x = "AAAAAATTTTTT"
                    w = record.seq

                    for i in range(0, len(str(w)) - 12):
                        if distance(x, str(w)[i:(i + 12)]) < NumOfDif:
                            for d, dirs, files in os.walk(way):  # find compared .ptt file and take info from there
                                for names2 in files:
                                    if (names[0:-4] + '.ptt') == names2:
                                        path2 = os.path.join(d, names2)
                                        if os.path.exists(path2):
                                            path3 = path.split('/')  # first column (name of file)
                                            handle2 = open(path2)
                                            loclist = []
                                            k = 0
                                            for loc in SplitForLoc(handle2):
                                                loclist.append(loc)
                                                k += 1  # mistake, should do it in the end, now k-1 everywhere

                                                if int(loc[0]) < i < int(loc[1]):  # in gene
                                                    cout.write(
                                                        path3[-2] + '\t' + str(w)[i:(i + 12)] + '\t' + str(i) + '\t' +
                                                        loc[3] + '\t' + loc[
                                                            2] + '\t' + 'none' + '\t' + 'none' + '\t' + 'none' + '\t' +
                                                        'none' + '\n')

                                                if int(loclist[k - 2][1]) < i < int(loclist[k - 1][0]):

                                                    if loclist[k - 2][3] == '-' and loclist[k - 1][3] == '-':
                                                        # strand is '-' both
                                                        cout.write(path3[-2] + '\t' + str(w)[i:(i + 12)] + '\t' +
                                                                   str(i) + '\t' + 'none' + '\t' + 'none' + '\t' +
                                                                   loclist[k - 2][2] + '\t' +
                                                                   str(i - int(loclist[k - 2][1])) + '\t' +
                                                                   loclist[k - 1][2] + '\t' +
                                                                   str(int(loclist[k - 1][0]) - i) + '\n')

                                                    elif (loclist[k - 2][3] == '+' and loclist[k - 1][
                                                        3] == '+'):  # strand is '+' both
                                                        cout.write(path3[-2] + '\t' + str(w)[i:(i + 12)] + '\t' +
                                                                   str(i) + '\t' + 'none' + '\t' + 'none' + '\t' +
                                                                   loclist[k - 1][2] + '\t' +
                                                                   str(int(loclist[k - 1][0]) - i) + '\t' +
                                                                   loclist[k - 2][2] + '\t' +
                                                                   str(i - int(loclist[k - 2][1])) + '\n')

                                                    elif (loclist[k - 2][3] == '+' and loclist[k - 1][
                                                        3] == '-'):  # first strand + and second - "down/down"
                                                        cout.write(path3[-2] + '\t' + str(w)[i:(i + 12)] + '\t' + str(
                                                            i) + '\t' + 'none' + '\t' + 'none' + '\t' + 'none' + '\t' +
                                                                   'none' + '\t' +
                                                                   loclist[k - 1][2] + '/' + loclist[k - 2][2] +
                                                                   '\t' + str(int(loclist[k - 1][0]) - i) +
                                                                   '/' + str(i - int(loclist[k - 2][1])) + '\n')

                                                    elif loclist[k - 2][3] == '-' and loclist[k - 1][3] == '+':
                                                        # first strand - and second + "upstream/upstream"

                                                        cout.write(path3[-2] + '\t' + str(w)[i:(i + 12)] + '\t' +
                                                                   str(i) + '\t' + 'none' + '\t' + 'none' + '\t' +
                                                                   loclist[k - 1][2] + '/' + loclist[k - 2][2] + '\t' +
                                                                   str(int(loclist[k - 1][0]) - i) + '/' +
                                                                   str(i - int(loclist[k - 2][1])) + '\t' + 'none' +
                                                                   '\t' + 'none' + '\n')

                                                    else:
                                                        print('ERROR.ERROR.ERROR.ERROR' + '\n')
                                            handle2.close()
                            i += 1

                handle.close()

cout.close()
