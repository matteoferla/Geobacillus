__author__ = 'Matteo'
__doc__ = '''Given a series of contigs and a reference genome. What are the position of the gene of the contigs acoording to the reference?
Get the two genomes and their alignment and
then get the position to make a map for the protein
try to match the protein with each other.
I have not fixed the issue with contig 4 and with the loss of everything beyond the origin. therefore so seq is lost

The main fnx to run the operatings in in true_place().

trans_ values are the transform of the position to the reference chromosome.
trans_ is trying to match the offset of spaces and such and is buggy
trans2_ uses a blastn alignment of the orf sequences onto the reference genome.
So trans_ may be off by a few bases, while trans2_ might mismatch theoretically



'''

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


def namepicker(query,choices):
    for name in choices:
        if str(query).find(name)>0:
            el= name
            break
    else:
        print("unmatched "+str(query))
        el=query
    return el

def offsetter(query,subj, start, end):
    off=query[start:end].count("-")-subj[start:end].count("-")
    if off>0:
        off+=offsetter(query,subj, end, end+off)
    elif off<0:
        off-=offsetter(query,subj, end, end+off)
    return off

def picker(rec,regex=None):
    if regex:
        idtemp=re.search(regex,rec.query)
        if idtemp:
            id=idtemp.group(0)
        else:
            id=rec.query
            print(rec.query+" not matched with "+regex)
    else:
        id=rec.query
    if rec.alignments:
        align=rec.alignments[0]
        hit=align.hsps[0]
        if hit.sbjct_start<hit.sbjct_end:
            strand="same"
        else:
            strand="opposite"
            #print(rec.query+" is antisense")
        el=namepicker(align.hit_def,'chromosome large small'.split())
        best={"id":id,
                "id_raw": rec.query,
                "element_raw":align.hit_def,
                "element":el,
                "length":rec.query_length,
                "query_coverage":(hit.query_end-hit.query_start+1)/rec.query_length,
                "plus":hit.positives,
                "query_start":hit.query_start,
                "query_end":hit.query_end,
                "sbj_start":hit.sbjct_start,
                "sbj_end":hit.sbjct_end,
                "sense":strand,
                "al_sbj":hit.sbjct,
                "al_query":hit.query}
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
                "al_query":0}
    return best

def aligntabler(af,txt=None):
    return [picker(rec,txt) for rec in NCBIXML.parse(open(af))]

def enrich(genome,genref, alignball,orfball):
    for c,contig in enumerate(genome):
        if alignball:
            #print(genome[c].name) #NZ_ ?!
            a=["NZ_"+x['id'] for x in alignball].index(genome[c].name)
            if alignball[a]['sense']=='opposite':
                alignball[a]['sbj_start'],alignball[a]['sbj_end']=alignball[a]['sbj_end'],alignball[a]['sbj_start']
                genome[c]=contig.reverse_complement()
                alignball[a]['al_query']=str(Seq(alignball[a]['al_query']).reverse_complement())
                alignball[a]['al_sbj']=str(Seq(alignball[a]['al_sbj']).reverse_complement())
                #the sub are already swapped
            setattr(genome[c],'blast',alignball[a])
            contig=genome[c] #now that an in place change has happned I can
            for gene in contig.features:
                offset=contig.blast['sbj_start']-(contig.blast['query_start']-1)
                start_offset=offset+offsetter(alignball[a]['al_query'],alignball[a]['al_sbj'],1,gene.location.start.position)
                end_offset=offset+offsetter(alignball[a]['al_query'],alignball[a]['al_sbj'],1,gene.location.end.position)
                #the maths is so mad I have not checked it out...
                setattr(gene.location,'trans_start',start_offset+gene.location.start.position)
                setattr(gene.location,'trans_end',end_offset+gene.location.end.position)
            #via blastn of seq.
        for gene in contig.features:
            if 'locus_tag' in gene.qualifiers and gene.type=="CDS":
                o=[x['id_raw'] for x in orfball].index(gene.qualifiers['locus_tag'][0])
                if orfball[o]['sbj_start'] < orfball[o]['sbj_end']:
                    setattr(gene.location,'trans2_start',orfball[o]['sbj_start']-1) #off by one. empirical. no idea why.
                    setattr(gene.location,'trans2_end',orfball[o]['sbj_end']-1)
                else:
                    setattr(gene.location,'trans2_end',orfball[o]['sbj_start']-1)
                    setattr(gene.location,'trans2_start',orfball[o]['sbj_end']-1)
                setattr(gene.location,'trans2_place',orfball[o]['element'])
                setattr(gene.location,'trans2_accuracy',orfball[o]['query_coverage'])




def find_gene(genome,where,name):
    if not name:
        return []
    c='chromosome large small'.split().index(name) #lol. I am so sick of this
    chr=genome[c]
    matched=[]
    if where[0]>where[1]:
            where=[where[1],where[0]]
    for gene in chr.features:
        if gene.type=="CDS": #Remember to check if they actually madly overide type()
            start=gene.location.start.position
            end=gene.location.end.position
            if start>end:
                start,end=end,start
            if start<=where[0]<=end or start<=where[1]<=end:
                matched.append(gene)
    return matched

def dratio(a,b):
    awords=set(a.split())
    bwords=set(b.split())
    k=0
    for c in awords.intersection(bwords):
        k+=len(c)+1
    if len(bwords)>len(awords):
        fixed, shuffled=bwords.difference(awords), awords.difference(bwords)
    else:
        shuffled, fixed=bwords.difference(awords), awords.difference(bwords)
    diff=[]
    do=Differ()
    fword=' '.join(fixed)
    for draw in itertools.permutations(shuffled):
        sword=' '.join(draw)
        diff.append(([x[0:1] for x in do.compare(fword,sword)].count(" ")+k)/(k+len(sword)))
    return sorted(diff)[-1]

def überlister(genome,genref):
    table=[]
    for chr in genome:
        for gene in chr.features:
            if gene.type=='CDS':
                print(gene.qualifiers['locus_tag'][0],gene.qualifiers['product'][0])
                gmatch=find_gene(genref,[gene.location.trans2_start,gene.location.trans2_end],gene.location.trans2_place)
                if gmatch:
                    for m in gmatch:
                        print('...comparison with ',m.qualifiers['locus_tag'][0],m.qualifiers['product'][0])
                        setattr(m,'diff',dratio(gene.qualifiers['product'][0].lower(),m.qualifiers['product'][0].lower()))
                    dgene=sorted(gmatch,key=lambda m: m.diff)[-1]
                   # mgene=sorted(gmatch,key=lambda m: len(m))[0]
                    mainid,mainprod=dgene.qualifiers['locus_tag'][0],dgene.qualifiers['product'][0]
                    maindiff=dratio(gene.qualifiers['product'][0].lower(),dgene.qualifiers['product'][0].lower())
                    mainstart,mainend=dgene.location.start.position,dgene.location.end.position
                else:
                    mainid,mainprod,mainstart,maindiff, mainend="","",0,0,0
                table.append({"locus":gene.qualifiers['locus_tag'][0],
                              "product":gene.qualifiers['product'][0],
                              "start":gene.location.start.position,
                              "end":gene.location.end.position,
                              "contig":chr.blast['id'],
                              "trans_start":gene.location.trans_start,
                              "trans_end":gene.location.trans_end,
                              "trans2_start":gene.location.trans2_start,
                              "trans2_end":gene.location.trans2_end,
                              "trans2_place":gene.location.trans2_place,
                              "trans2_accuracy":gene.location.trans2_accuracy,
                              "trans_place":chr.blast['element'],
                              "occupancy":len(gmatch),
                              "occupancy_detail":[gene.qualifiers['locus_tag'][0]+" ("+str(gene.location.start.position)+","+str(gene.location.end.position)+")" for gene in gmatch],
                              "main_occupant":mainid,
                              "main_occupant_prod":mainprod,
                              "main_occupant_start":mainstart,
                              "main_occupant_end":mainend,
                              "difference":maindiff})
    return table


def überlister_control(genome,genref):
    table=[]
    for chr in genome:
        for gene in chr.features:
            if gene.type=='CDS':
                print(gene.qualifiers['locus_tag'][0],gene.qualifiers['product'][0])
                gmatch=find_gene(genref,[gene.location.trans2_start,gene.location.trans2_end],gene.location.trans2_place)
                if gmatch:
                    for m in gmatch:
                        print('...comparison with ',m.qualifiers['locus_tag'][0],m.qualifiers['product'][0])
                        setattr(m,'diff',dratio(gene.qualifiers['product'][0].lower(),m.qualifiers['product'][0].lower()))
                    dgene=sorted(gmatch,key=lambda m: m.diff)[-1]
                   # mgene=sorted(gmatch,key=lambda m: len(m))[0]
                    mainid,mainprod=dgene.qualifiers['locus_tag'][0],dgene.qualifiers['product'][0]
                    maindiff=dratio(gene.qualifiers['product'][0].lower(),dgene.qualifiers['product'][0].lower())
                    mainstart,mainend=dgene.location.start.position,dgene.location.end.position
                else:
                    mainid,mainprod,mainstart,maindiff, mainend="","",0,0,0
                table.append({"locus":gene.qualifiers['locus_tag'][0],
                              "product":gene.qualifiers['product'][0],
                              "start":gene.location.start.position,
                              "end":gene.location.end.position,
                              "contig":"NA",
                              "trans_start":"NA",
                              "trans_end":"NA",
                              "trans2_start":gene.location.trans2_start,
                              "trans2_end":gene.location.trans2_end,
                              "trans2_place":gene.location.trans2_place,
                              "trans2_accuracy":gene.location.trans2_accuracy,
                              "trans_place":"NA",
                              "occupancy":len(gmatch),
                              "occupancy_detail":[gene.qualifiers['locus_tag'][0]+" ("+str(gene.location.start.position)+","+str(gene.location.end.position)+")" for gene in gmatch],
                              "main_occupant":mainid,
                              "main_occupant_prod":mainprod,
                              "main_occupant_start":mainstart,
                              "main_occupant_end":mainend,
                              "difference":maindiff})
    return table


def lister(genome):
    table=[]
    for chr in genome:
        for gene in chr.features:
            if gene.type=='CDS':
                el=namepicker(chr.description,'chromosome large small'.split())
                table.append({"locus":gene.qualifiers['locus_tag'][0],"product":gene.qualifiers['product'][0],"place":el,"start":gene.location.start.position, "end":gene.location.end.position})
    return table

def csvwrapper(fn,data,headers):
    file=open(fn,'w')
    out=csv.DictWriter(file,headers,lineterminator='\n')
    out.writeheader()
    out.writerows(data)
    file.close()

def true_place():
        #load genomes
    dsm=get_genome("LAKX01.1.gb")
    tm=get_genome("Gthg_TM242_v3.0.gb")
    #get blast data and add it to dsm as .blast
    al=aligntabler("align.xml",'LAKX010\d+')
    al2=aligntabler("align_orf2.xml")
    enrich(dsm,tm,al,al2)
    fieldnames="locus product start end contig "+\
               "trans_start trans_end trans_place "+\
               "trans2_start trans2_end trans2_place trans2_accuracy "+\
               "occupancy occupancy_detail main_occupant main_occupant_start main_occupant_end main_occupant_prod difference"
    csvwrapper("truesite.csv",überlister(dsm,tm),fieldnames.split())
    #csvwrapper("ref.csv",lister(tm),"locus product place start end".split())


def true_place_control():
        #load genomes
    tm=get_genome("Gthg_TM242_v3.0.gb")
    #get blast data and add it to dsm as .blast
    al2=aligntabler("align_orf_ref.xml")
    enrich(tm,tm,[],al2)
    fieldnames="locus product start end contig "+\
               "trans_start trans_end trans_place "+\
               "trans2_start trans2_end trans2_place trans2_accuracy "+\
               "occupancy occupancy_detail main_occupant main_occupant_start main_occupant_end main_occupant_prod difference"
    csvwrapper("truesite.csv",überlister_control(tm,tm),fieldnames.split())
    #csvwrapper("ref.csv",lister(tm),"locus product place start end".split())

def orfome_fasta(genome,name):
    n=0
    w=open(name+"_orf.fa","w")
    for chr in genome:
        for gene in chr.features:
            if gene.type=='CDS':
                n+=1
                w.writelines([">"+gene.qualifiers['locus_tag'][0]+N,str(gene.extract(chr).seq)+N])
    print(n)


def get_prot():
    __doc__='''Parse the blast results of prot vs. prot and give a table '''
    genome = get_genome("Gthg_TM242_v3.0.gb")
    w=csv.DictWriter(open("prot_cons.csv","w"),"qid sid len qst qnd sst snd qsq ssq m mut V1M nm".split())
    w.writeheader()
    for rec in NCBIXML.parse(open("align_orf.xml")):
        w.writerow(picker(rec))

def check_two(genome):
    for chr in genome:
        for gene in chr.features:
            if 'translation' in gene.qualifiers:
                emp=str(gene.extract(chr).seq.translate(to_stop=1))
                exp=gene.qualifiers['translation'][0]
                d=[x[0:1] for x in Differ().compare(exp,emp)].count(" ")/max(len(emp),len(exp))
                print(d,gene.qualifiers['locus_tag'][0])

def check_three(): #UNfINISHED
            #load genomes
    dsm=get_genome("LAKX01.1.gb")
    tm=get_genome("Gthg_TM242_v3.0.gb")
    #get blast data and add it to dsm as .blast
    al=aligntabler("align.xml",'LAKX010\d+')
    al2=aligntabler("align_orf2.xml")
    enrich(dsm,tm,al,al2)
    fieldnames="locus length_np f1 f2 f3 r1 r2 r3"
    lenball=[]
    for chr in dsm:
        for gene in chr.features:
            if gene.type=='CDS':
                ref=tm['chromosome large small'.split().index(gene.location.trans2_place)]
                frames={'locus':gene.qualifiers['locus_tag'][0], 'length_np':gene.location.trans2_end-gene.location.trans2_start}
                print(frames)
                for x in (1,2,3):
                    frames['f'+str(x)]=len(ref[gene.location.trans2_start+x:gene.location.trans2_end+x].seq.translate(to_stop=1))
                    frames['r'+str(x)]=len(ref[gene.location.trans2_start+x:gene.location.trans2_end+x].seq.reverse_complement().translate(to_stop=1))
                lenball.append(frames)
    csvwrapper("check3.csv",lenball,fieldnames.split())

if __name__ == "__main__":
    #dsm=get_genome("LAKX01.1.gb")
    #orfome_fasta(dsm,"DSM")
    #SeqIO.convert("Gthg_TM242_v3.0.gb", "genbank", "Gthg.fasta", "fasta")
    #orfome_fasta(get_genome("Gthg_TM242_v3.0.gb"),"ref_orfome.fa")
    #true_place()
    #true_place_control()
    #check_two(get_genome("Gthg_TM242_v3.0.gb"))
    check_three()
