__author__ = 'Matteo'
__doc__ = '''Circos may be pretty, but it is really unfriendly. This is my code to make the input files from a genbank'''

N = "\n"
T = "\t"
# N="<br/>"
from Bio import SeqIO
from Bio import Seq
import csv, math

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

def highlight_distr(fp, high,span):
    __doc__='''Get the distribution of transposase. I am sure there is a witty way of doing it with vlookup and pivot tables,
    but I want to get dinner and this won't take a second'''
    genome=get_genome(fp)
    for ci, chr in enumerate(genome):
        g=[0]*(math.floor(len(chr)/span)+2) #there must me a more gangsta way of zeroes().
        for gene in chr.features:
            if gene.type=="CDS":
                if gene.qualifiers['product'][0].lower().find(high) >-1:
                    up=gene.location.start
                    down=gene.location.end
                    for x in range(math.floor(up/span),math.ceil(down/span)):
                            g[x]+=1  #I did g[x]++ again. Sob.

        open(high+"_dist_"+str(ci)+".txt","w").writelines([str(j)+N for j in g])

def unoverlap_genes():
    __doc__='''For Circos overlapping genes is a bad thing'''
    genome=get_genome("Gthg_TM242_v3.0.gb")
    table=[]
    for chr in genome:
        genes=[[gene.qualifiers['locus_tag'][0],gene.location.start,gene.location.end] for gene in chr.features if gene.type.lower() !='source']
        n=2 #triple overlap.
        while n>0:
            for i,g in enumerate(genes):
                if i == len(genes)-1:
                    pass
                elif genes[i]==None: #gene has been killed
                    pass
                elif genes[i][2]>=genes[i+1][2]: #kill the gene wholly overlapped
                    print(genes[i+1][0])
                    genes[i+1]=None
                elif genes[i][2]>=genes[i+1][1]: #overlap
                    if genes[i][2]-genes[i][1] > genes[i+1][2]-genes[i+1][1]: #trim next if smaller
                        genes[i+1][1]=genes[i][2]+1
                    else:   #trim this
                        genes[i][2]=genes[i+1][1]-1
                else:  #no overlap
                    pass
            n=n-1
            print(len(genes))
            genes=[g for g in genes if g]
        table.extend(genes)
    return table

def csvwrapper(fn,data,headers):
    file=open(fn,'w')
    out=csv.DictWriter(file,headers,lineterminator='\n')
    out.writeheader()
    out.writerows(data)
    file.close()




if __name__ == "__main__":
    data=unoverlap_genes()
    file=open('unoverlap_list.csv','w')
    out=csv.writer(file,lineterminator='\n')
    out.writerows(data)
    file.close()
    #highlight_distr("Gthg_TM242_v3.0.gb","transposase",10000)