__author__ = 'Matteo'
__doc__='''Comparison of two genomes. One whole, the other in contigs.
As I am currently working on a windows, the two blast commands were done manually'''

from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SeqIO
import csv

#cmd=NcbiblastnCommandline(out='align.xml', outfmt=5, query='LAKX01.1.fsa_nt', db='Gthg.fa', evalue=0.001)()
#Fecking doesn;t work. it is not the cmd='blastn' or the quotes

#Blasted manually in cmd

def get_genome(fp):
    return  list(SeqIO.parse(open(fp, "rU"), "genbank"))

def genelooper(genome):
    return [gene for chr in genome for gene in chr.features if gene.type =="CDS"]

def get_diff():
    __doc__='''A slightly misguided module. It spits the list of positions of the contigs against the reference genome.
    It relyes on a function called picker whose working is a tad silly as blast results are sorted already.'''
    genome = get_genome("Gthg_TM242_v3.0.gb")

    def picker_old(rec):
        best={"id": rec.query_id,"len":rec.query_length,"cov":0,"plus":0,"gaps":0}
        for a,align in enumerate(rec.alignments):
            for hit in align.hsps:
                cov=(hit.query_end-hit.query_start+1)/rec.query_length
                if cov>best['cov']:
                    best={"id": rec.query_id,"len":rec.query_length,"cov":cov,"plus":hit.positives,"gaps":hit.gaps}
        return best


    def picker(rec):
        best={"id": rec.query_id,"element":"","len":rec.query_length,"cov":0,"plus":0,"diff":0,"start":-1,"end":-1}
        old_pair=[0,0]
        old_sub=[]
        for a,align in enumerate(rec.alignments):
            for hit in align.hsps:
                pair=[hit.query_start,hit.query_end]
                if hit.sbjct_start<hit.sbjct_end:
                    sub=[hit.sbjct_start,hit.sbjct_end,1]
                else:
                    sub=[hit.sbjct_end,hit.sbjct_start,-1] #so weird sub's not reserved.
                cov=(hit.query_end-hit.query_start+1)/rec.query_length
                if pair[0]>old_pair[1]:
                    best['element']=align.hit_def
                    best["cov"]+=cov
                    best["plus"]+=hit.positives
                    old_pair=pair
                    if old_sub:
                        best["diff"]+=old_sub[1]-sub[0]
                        if sub[2]==1:
                            best['end']=sub[1]
                        else:
                            best['start']=sub[0]
                    else:
                        best["start"]=sub[0]
                        best["end"]=sub[1]
                    old_sub=sub
        return best

    #hit attr
    #'align_length', 'bits', 'expect', 'frame', 'gaps', 'identities', 'match', 'num_alignments', 'positives',
    # 'query', 'query_end', 'query_start', 'sbjct', 'sbjct_end', 'sbjct_start', 'score', 'strand']
    #"LAKX01.1.gb"
    w=csv.DictWriter(open("synteny2.csv","w"),"id len cov plus diff start end element".split())
    w.writeheader()
    for rec in NCBIXML.parse(open("align2.xml")):
        w.writerow(picker(rec))

def get_prot():
    __doc__='''Parse the blast results of prot vs. prot and give a table '''

    def picker(rec):
        if rec.alignments:
            hit=rec.alignments[0].hsps[0]
            mut=[]
            mut2=""
            for i,m in enumerate(hit.match):
                if m==" " or m=="+":
                    if i==0:
                        mut2=hit.sbjct[i]+str(i+int(hit.sbjct_start))+hit.query[i]
                    else:
                        mut.append(hit.sbjct[i]+str(i+int(hit.sbjct_start))+hit.query[i])
            top={"qid": rec.query,"sid":rec.alignments[0].hit_def,"len":rec.query_length,
                 "qst":hit.query_start,"qnd":hit.query_end,"sst":hit.sbjct_start,"snd":hit.sbjct_end,"ssq":hit.sbjct,"qsq":hit.query,"m":hit.match,"mut":mut,"V1M":mut2,"nm":len(mut)}
        else:
            top={"qid": rec.query,"sid":"#UNMATCHED","len":rec.query_length,
                 "qst":0,"qnd":0,"sst":0,"snd":0,"ssq":"","qsq":"","m":"","nm":"","V1M":""}
        return top

    w=csv.DictWriter(open("prot_cons.csv","w"),"qid sid len qst qnd sst snd qsq ssq m mut V1M nm".split())
    w.writeheader()
    for rec in NCBIXML.parse(open("pralign.xml")):
        w.writerow(picker(rec))

if __name__=="__main__":
    pass