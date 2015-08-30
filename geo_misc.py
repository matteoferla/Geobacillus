__author__ = 'Matteo'
__doc__=''' Misc. scripts'''

N="\n"

from Bio import SeqIO
from Bio import Seq
import csv


def get_genome(fp):
    return  list(SeqIO.parse(open(fp, "rU"), "genbank"))
def genelooper(genome):
    return [gene for chr in genome for gene in chr.features if gene.type =="CDS"]


def lister(genome):
    table=[]
    for gene in [x for chr in genome for x in chr.features]:
        if gene.type=='CDS':
            table.append({"locus":gene.qualifiers['locus_tag'][0],"product":gene.qualifiers['product'][0],"start":gene.location.start.position, "end":gene.location.end.position})
    return table

def protein_fasta(genome,name):
    n=0
    w=open(name+"_prot.fa","w")
    for chr in genome:
        for gene in chr.features:
            if gene.type=='CDS':
                n+=1
                if 'translation' in gene.qualifiers:
                    w.writelines([">"+gene.qualifiers['locus_tag'][0]+N,gene.qualifiers['translation'][0]+N])
                else:
                    w.writelines([">"+gene.qualifiers['locus_tag'][0]+N,str(gene.extract(chr).seq.translate(to_stop=1))+N])
    print(n)


def rRNA_fastamaker(genome):
    __doc__='''Make a fasta of a single 5S, 6S and 23S gene, plus whatever tRNAs you can find.'''
    seq=["","",""]
    count=[0,0,0]
    ribo="5S 16S 23S".split()
    trna=[]
    for chr in genome:
        for gene in chr.features:
            if gene.type =="rRNA":
                for xi,x in enumerate(ribo):
                    if gene.qualifiers['note'][0].find(x)>-1:
                        seq[xi]=gene.extract(chr).seq
                        count[xi]+=1
                        break
                else:
                    print("ERROR")
                    print(gene)
            if gene.type =="tRNA":
                seq.append(gene.extract(chr).seq)
                trna.append(gene.qualifiers['note'][0])


    for xi,x in enumerate(ribo+trna):
        print(">"+x+"\n"+seq[xi])
    print(count)




if __name__== "__main__":
    pass