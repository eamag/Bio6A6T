#! /usr/bin/env python
# -*- coding: utf-8 -*- 
from Bio import SeqIO
from Levenshtein import distance
import os


def splitout(handle):
    line = handle.readline().strip()
    if line.rstrip() <> ('Name\tMotif'+'\t'+'Coordinate'+'\t'+'Strand'+'\t'+'Name'+'\t'+'Upstream'+'\t'+
                             'Upstream dist'+'\t'+'Downstream'+'\t'+'Downstream dist') :
        raise SyntaxError("Third line not recognised as OUT")
    # NAME=0
    # MOTIF=1
    # COORDINATE=2
    # STRAND=3
    # NAME=4
    # UPSTREAM=5
    # UPSTREAMD=6
    # DOWNSTREAM=7
    # DOWNSTREAMD=8
    for line in handle:
        parts=line.rstrip().split('\t')
        yield parts


def analysis(lst,dct):

    for i in lst:
        if i in dct:
            dct[i] += 1
        else:
            dct[i] = 1

InGene = open('InGene.txt', 'w')
Upstream = open('Upstream.txt','w')
Downstream = open('Downstream.txt','w')
InGene.write('Name\tMotif'+'\t'+'Coordinate'+'\t'+'Strand'+'\t'+'Name'+'\t'+'Upstream'+
             '\t'+'Upstream dist'+'\t'+'Downstream'+'\t'+'Downstream dist'+'\n')
Upstream.write('Name\tMotif'+'\t'+'Coordinate'+'\t'+'Strand'+'\t'+'Name'+'\t'+'Upstream'+
               '\t'+'Upstream dist'+'\t'+'Downstream'+'\t'+'Downstream dist'+'\n')
Downstream.write('Name\tMotif'+'\t'+'Coordinate'+'\t'+'Strand'+'\t'+'Name'+'\t'+'Upstream'+
                 '\t'+'Upstream dist'+'\t'+'Downstream'+'\t'+'Downstream dist'+'\n')

gene = []
gened={}
upstream = []
upsd={}
downstream = []
downd={}
path='/home/eamag/Bio/motif_out.txt'
handle = open(path)
for parts in (splitout(handle)):
    if parts[4] != 'none':
        for item in parts:
            InGene.write(item+'\t')
        InGene.write('\n')
    if parts[5] != 'none':
        for item in parts:
            Upstream.write(item+'\t')
        Upstream.write('\n')
    if parts[7] != 'none':
        for item in parts:
            Downstream.write(item+'\t')
        Downstream.write('\n')
    if (parts[4]!= 'none' and parts[4] != 'hypothetical protein'):
        gene.append(parts[4])
    up=parts[5].split('/')
    up.append('hypothetical protein')
    down=parts[7].split('/')
    down.append('hypothetical protein')
    if (parts[5]!= 'none' and up[0]!= 'hypothetical protein'):
        upstream.append(up[0])
    if (parts[5]!= 'none' and up[1]!= 'hypothetical protein'):
        upstream.append(up[1])
    if (parts[7]!= 'none' and down[0]!= 'hypothetical protein'):
        downstream.append(down[0])
    if (parts[7]!= 'none' and down[1]!= 'hypothetical protein'):
        downstream.append(down[1])


handle.close()
InGene.close()
Upstream.close()
Downstream.close()
analysis(gene,gened)
analysis(upstream,upsd)
analysis(downstream,downd)
AGene=open('AnalysisGene.txt','w')
AUps=open('AnalysisUpstream.txt','w')
ADowns=open('AnalysisDownstream.txt','w')

for i in (sorted(gened.items(), key=lambda (k, v): v, reverse=True)):
    AGene.write(str(i[0])+'\t'+str(i[1])+'\n')
for i in (sorted(upsd.items(), key=lambda (k, v): v, reverse=True)):
    AUps.write(str(i[0])+'\t'+str(i[1])+'\n')
for i in (sorted(downd.items(), key=lambda (k, v): v, reverse=True)):
    ADowns.write(str(i[0])+'\t'+str(i[1])+'\n')

AGene.close()
AUps.close()
ADowns.close()