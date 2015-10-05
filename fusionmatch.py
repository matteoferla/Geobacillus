__author__ = 'Matteo'
__doc__ = ''''''

N = "\n"
T = "\t"
# N="<br/>"

from Bio import SeqIO
from Bio import Seq
import csv
import math
from Bio.Blast import NCBIXML

def get_genome(fp):
    return  list(SeqIO.parse(open(fp, "rU"), "genbank"))
def genelooper(genome):
    return [gene for chr in genome for gene in chr.features if gene.type =="CDS"]

def depr_joiner(chr, a,b,c,reversed=False,to_stop=False):  #this assumes that the locations of the genes are set right.
    if reversed:
        fore=chr[a:b-1-((c-b+1) % 3)].seq
        print((c-b+1) % 3 +1)
        aft=chr[b:c].seq
        new=fore+aft
        return  new.reverse_complement().translate(to_stop=to_stop)
    else:
        fore=chr[a:b+((a-b+1) % 3)].seq
        aft=chr[b:c].seq
        new=fore+aft
        return new.translate(to_stop=to_stop)

def joiner1(chr, a,b,c,reversed=False,to_stop=True):  #junked
    prot=[]
    target=(c-a+1)/3
    if reversed:
        prot.append(str(chr[a:c].seq.reverse_complement().translate(to_stop=to_stop)))
        for n in range(60):
            new=chr[a:b].seq+chr[b-n:c].seq
            prot.append(str(new.reverse_complement().translate(to_stop=to_stop)))
            new=chr[a:b-n].seq+chr[b:c].seq
            prot.append(str(new.reverse_complement().translate(to_stop=to_stop)))
    else:
        prot.append(str(chr[a:c].seq.translate(to_stop=to_stop)))
        for n in range(60):
            new=chr[a:b-n].seq+chr[b:c].seq
            prot.append(str(new.translate(to_stop=to_stop)))
            new=chr[a:b].seq+chr[b-n:c].seq
            prot.append(str(new.translate(to_stop=to_stop)))
    p2=sorted(prot,key=lambda p: len(p))
    print(target,len(p2[-1]))
    #print(N.join(p2))
    return p2[-1]

def joiner(chr, a,b,c,reversed=False,to_stop=True): #shoddy
    target=(c-a+1)/3
    if reversed:
        trans=chr[a:c].seq.reverse_complement().translate(to_stop=to_stop)
        print(len(trans),trans)
        n=1
        while len(trans)< target-2-(n/1.5):
            seqs=[]
            print(len(trans),target-2-(n/1.5))
            new=chr[a:b].seq+chr[b-n:c].seq
            seqs.append(new.reverse_complement().translate(to_stop=to_stop))
            new=chr[a:b-n].seq+chr[b:c].seq
            seqs.append(new.reverse_complement().translate(to_stop=to_stop))
            new=chr[a:b-n].seq+chr[b-n:c].seq
            seqs.append(new.reverse_complement().translate(to_stop=to_stop))
            trans=max(seqs,key=lambda s:len(s))
            n+=1
    else:
        prot.append(str(chr[a:c].seq.translate(to_stop=to_stop)))
        for n in range(60):
            new=chr[a:b-n].seq+chr[b:c].seq
            prot.append(str(new.translate(to_stop=to_stop)))
            new=chr[a:b].seq+chr[b-n:c].seq
            prot.append(str(new.translate(to_stop=to_stop)))
    p2=sorted(prot,key=lambda p: len(p))
    print(target,len(p2[-1]))
    #print(N.join(p2))
    return p2[-1]

def make_slipcandidates(genome):
    fasta=open('slip_candidates.fa','w')
    threshhold=50
    for chr in genome:
        previous=None
        for gene in chr.features:
            if gene.type=='CDS':
                if previous:
                    if gene.location.strand==previous.location.strand and (gene.qualifiers['product'][0]==previous.qualifiers['product'][0] or gene.location.start-previous.location.end<threshhold):
                        if gene.location.strand<0:
                            gene.qualifiers['locus_tag'][0]=gene.qualifiers['locus_tag'][0]+'-'+previous.qualifiers['locus_tag'][0]
                            gene.qualifiers['translation'][0]=gene.qualifiers['translation'][0]+previous.qualifiers['translation'][0]
                        else:
                            gene.qualifiers['locus_tag'][0]=previous.qualifiers['locus_tag'][0]+'-'+gene.qualifiers['locus_tag'][0]
                            gene.qualifiers['translation'][0]=previous.qualifiers['translation'][0]+gene.qualifiers['translation'][0]
                    else:
                        fasta.write('>'+previous.qualifiers['locus_tag'][0]+N+previous.qualifiers['translation'][0]+N)
                previous=gene

def locus2prod(genome,infile,outfile):
    text=open(infile,'r').read()
    text=text.replace('-','---')
    for gene in genelooper(genome):
        text=text.replace(gene.qualifiers['locus_tag'][0],gene.qualifiers['product'][0])
    open(outfile,'w').write(text)


def get_prot(infile,outfile):
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

    w=csv.DictWriter(open(outfile,"w"),"qid sid qlen qst qnd sst snd qsq ssq m mut V1M nm".split())
    w.writeheader()
    for rec in NCBIXML.parse(open(infile)):
        w.writerow(picker(rec))


if __name__ == "__main__":
    genome=get_genome("Gthg_TM242_v3.0.gb")
    #make_slipcandidates(genome)
    #get_prot("enriched.xml","prot_cons_fusion.csv")
    locus2prod(genome,'mutanda.txt','mutata.txt')