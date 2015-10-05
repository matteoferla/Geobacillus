__author__ = 'Matteo'
__doc__ = ''''''

N = "\n"
T = "\t"
# N="<br/>"

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
import csv
import math
import re
from difflib import Differ
import itertools

def get_genome(fp):
    return  list(SeqIO.parse(open(fp, "rU"), "genbank"))

def genelooper(genome):
    return [gene for chr in genome for gene in chr.features if gene.type =="CDS"]


def check4():
    fieldnames="id id_raw element_raw element length query_coverage plus query_start query_end sbj_start sbj_end sense al_sbj al_query mutationsT2D mutations".split()
    w=csv.DictWriter(open('TMorf_onContigs.csv','w'),fieldnames=fieldnames)
    w.writeheader()
    for rec in NCBIXML.parse(open("align4.xml")):
        id=rec.query
        if rec.alignments:
            align=rec.alignments[0]
            hit=align.hsps[0]
            if hit.sbjct_start<hit.sbjct_end:
                strand="same"
            else:
                strand="opposite"
                #print(rec.query+" is antisense)
            ps=Seq(hit.sbjct.replace('-','')).translate()
            pq=Seq(hit.query.replace('-','')).translate()
            mutations=[]
            for i in range(len(ps)):
                if len(pq)-1<i:
                    mutations.append(str(ps[i])+str(i+1)+'?')
                elif ps[i] != pq[i]:
                    mutations.append(str(ps[i])+str(i+1)+str(pq[i]))
            best={"id":id,
                    "id_raw": rec.query,
                    "element_raw":align.hit_def,
                    "element":'NA',
                    "length":rec.query_length,
                    "query_coverage":(hit.query_end-hit.query_start+1)/rec.query_length,
                    "plus":hit.positives,
                    "query_start":hit.query_start,
                    "query_end":hit.query_end,
                    "sbj_start":hit.sbjct_start,
                    "sbj_end":hit.sbjct_end,
                    "sense":strand,
                    "al_sbj":hit.sbjct,
                    "al_query":hit.query,
                    "mutationsT2D":"+".join(mutations),
                    "mutations":len(mutations)}
        else:
            best={"id":id,
                    "id_raw": rec.query,
                    "element_raw":"",
                    "element":"",
                    "length":rec.query_length,
                    "query_coverage":0,
                    "plus":0,
                    "query_start":0,
                    "query_end":0,
                    "sbj_start":0,
                    "sbj_end":0,
                    "sense":0,
                    "al_sbj":0,
                    "al_query":0,
                    "mutationsT2D":"NA",
                    "mutations":"NA"}
        w.writerow(best)

if __name__ == "__main__":
    check4()