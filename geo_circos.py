__author__ = 'Matteo'
__doc__ = '''Circos may be pretty, but it is really unfriendly. This is my code to make the input files from a genbank'''

N = "\n"
T = "\t"
# N="<br/>"
from Bio import SeqIO
from Bio import Seq
import csv

def get_genome(fp):
    return  list(SeqIO.parse(open(fp, "rU"), "genbank"))
def genelooper(genome):
    return [gene for chr in genome for gene in chr.features if gene.type =="CDS"]


def highlight_circos(fp, high,chr_names):
    __doc__='''tile file to make all genes grey bar the ones with product name in "high" string.
    The pretty chromosome names are third element'''
    genome=get_genome(fp)
    for ci, chr in enumerate(genome):
        #prefix="band chr1 1.1 1.1 "
        prefix="chr1 "
        open("chr_trans"+str(ci)+".txt","w").writelines(["chr - chr1 "+chr_names[ci]+" 0 "+str(len(chr))+" black"])
        w=open("tiles_trans"+str(ci)+".txt","w")
        for gene in chr.features:
            if gene.type=="CDS":
                if gene.qualifiers['product'][0].lower().find(high) >-1:
                    suffix=" color=red"
                else:
                    suffix=" color=grey"  ###GREY not GRAY
                w.write(prefix+str(gene.location.start)+" "+str(gene.location.end)+suffix+N)
        w.close()
        #[print() for gene in chr.features if gene.type=="CDS" and gene.qualifiers['product'][0].lower().find("transposase") >-1]


if __name__ == "__main__":
    pass