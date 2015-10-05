__author__ = 'Matteo'
__doc__=''' Misc. scripts'''

N="\n"

from Bio import SeqIO
from Bio import Seq
import csv
import math
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm


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

def protein_length(genome):
    for chr in genome:
        for gene in chr.features:
            if gene.type=='CDS':
                print(len(gene))


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

def in_silico_frag(genome,size):
    new=[]
    for chr in genome:
        chunks =math.ceil(len(chr)/size)
        for i in range(chunks):
            new.append(chr[i*size:(i+1)*size])
    return new

def locusmap(genome):
    for gene in genelooper(genome):
        if 'locus_tag' in gene.qualifiers and 'old_locus_tag' in gene.qualifiers:
            print(gene.qualifiers['locus_tag'][0]+"\t"+gene.qualifiers['old_locus_tag'][0])

def original_draw():
    g=get_genome("Gthg_TM242_v3.0.gb")[0][1077101:1137446]
    #SeqIO.write([g[3720000:3727000]], "clipping.gb", "genbank")
    gd_diagram = GenomeDiagram.Diagram(g.id)
    gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
    gd_feature_set = gd_track_for_features.new_set()
    for feature in g.features:
        if feature.qualifiers['product'][0].lower().find('hypo')>-1:
            color = colors.lightblue
            gd_feature_set.add_feature(feature, sigil="ARROW",
                                   color=color, label=False)
        elif feature.qualifiers['product'][0].lower().find('transposase')>-1:
            color = colors.lightgrey
            gd_feature_set.add_feature(feature, sigil="ARROW",
                                   color=color, label=False)
        else:
            color = colors.blue
            shortened=feature.qualifiers['product'][0]
            if len(shortened)>20:
                shortened=shortened[0:20]+'...'
            gd_feature_set.add_feature(feature, sigil="ARROW",
                                   color=color, label=True,
                                   name=shortened,
                                   label_size = 8, label_angle=30)
    gd_diagram.draw(format="linear", pagesize='LEGAL', orientation='landscape', fragments=4,
            start=0, end=len(g))
    gd_diagram.write("C:\\Users\\Cass\\Desktop\\Geo 2\\hemicellulose\\hemi.pdf", "PDF")
    gd_diagram.write("C:\\Users\\Cass\\Desktop\\Geo 2\\hemicellulose\\hemi.eps", "EPS")
    gd_diagram.write("C:\\Users\\Cass\\Desktop\\Geo 2\\hemicellulose\\hemi.svg", "SVG")
    gd_diagram.write("C:\\Users\\Cass\\Desktop\\Geo 2\\hemicellulose\\hemi.png", "PNG")


def draw(genoslice, dataset,op,intensity=False):
    #SeqIO.write([g[3720000:3727000]], "clipping.gb", "genbank")
    gd_diagram = GenomeDiagram.Diagram(genoslice.id)
    gd_track_for_features = gd_diagram.new_track(1, name="Annotated Features")
    gd_feature_set = gd_track_for_features.new_set()
    for feature in genoslice.features:
        if feature.type =="CDS":
        # feature.qualifiers['product'][0].lower().find('hypo')>-1:
        # feature.qualifiers['product'][0].lower().find('transposase')>-1:
            locus=feature.qualifiers['locus_tag'][0]
            red,green, blue=(1,1,1)
            border=colors.gainsboro
            if locus in dataset:
                if intensity:
                    intense=float(dataset[locus][1])
                    red,green=(1-intense/10,1-intense/10)
                else:
                    change=float(dataset[locus][1])
                    sig=float(dataset[locus][2])
                    if change >0:
                        red,blue=(1-change/10,1-change/10)
                    else:
                        green,blue=(1-abs(change/10),1-abs(change/10))
                    if sig<0.05:
                        border=colors.black
                    else:
                        border=colors.gray
            color = colors.Color(red,green,blue)
            shortened=feature.qualifiers['locus_tag'][0]+' '+feature.qualifiers['product'][0]
            if len(shortened)>20:
                shortened=shortened[0:20]+'...'
            gd_feature_set.add_feature(feature, sigil="ARROW",
                                   color=color, label=True,
                                   name=shortened,
                                   border=border,
                                   label_size = 8, label_angle=30)
    gd_diagram.draw(format="linear", pagesize='LEGAL', orientation='landscape', fragments=4,
            start=0, end=len(g))
    gd_diagram.write(op+".pdf", "PDF")
    gd_diagram.write(op+".eps", "EPS")
    gd_diagram.write(op+".svg", "SVG")
    gd_diagram.write(op+".png", "PNG")

def ptt(chr):
    ptt=csv.DictWriter(open('g.ptt','w'),fieldnames='Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product'.split(),delimiter='\t',lineterminator = '\n')
    ptt.writeheader()
    fasta=open('g_prot.fasta','w')
    for gene in chr.features:
        if gene.type=='CDS':
           #>gi|38505774|ref|NP_942288.1| exodeoxyribonuclease V alpha chain [Synechocystis sp. PCC 6803]
            fasta.write('>gi|'+gene.qualifiers['locus_tag'][0]+'|ref|'+gene.qualifiers['locus_tag'][0]+'| '+gene.qualifiers['product'][0]+' [Geobacillus thermoglucosidasius TM242 complete chromosome]\n'+str(gene.location.extract(chr).seq.translate())+'\n')
            if gene.location.strand >0:
                strand='+'
            else:
                strand='-'
            ptt.writerow({'Location':str(gene.location.start)+'..'+str(gene.location.end),
                          'Strand':strand,
                          'Length':str(math.floor(gene.location.end/3-gene.location.start/3)),
                          'PID':gene.qualifiers['locus_tag'][0],
                          'Gene':gene.qualifiers['locus_tag'][0],
                          'Synonym':'-',
                          'Code':'-',
                          'COG':'-',
                          'Product':gene.qualifiers['product'][0]
            })
            #Location	Strand	Length	PID	Gene	Synonym	Code	COG	Product
            #3093..3851	-	252	-	-	-	Gthg00024	-	Hypothetical protein

def framejoin(chr, a,b,c,reversed=False):
    if reversed:
        print('whole section: ',chr[a:c].seq.reverse_complement().translate())
        new=chr[a:b-1].seq+chr[b:c].seq
        print(len(new)/3,'frame 1: ',new.reverse_complement().translate())
        new=chr[a:b-2].seq+chr[b:c].seq
        print(len(new)/3,' frame 2: ',new.reverse_complement().translate())
        new=chr[a:b-3].seq+chr[b:c].seq
        print(len(new)/3,'frame 3: ',new.reverse_complement().translate())
    else:
        print('whole section: ',chr[a:c].seq.translate())
        new=chr[a:b-1].seq+chr[b:c].seq
        print(len(new)/3,'frame 1: ',new.translate())
        new=chr[a:b-2].seq+chr[b:c].seq
        print(len(new)/3,' frame 2: ',new.translate())
        new=chr[a:b-3].seq+chr[b:c].seq
        print(len(new)/3,'frame 3: ',new.translate())


def alottable(genome):
    w=csv.DictWriter(open('detailed.csv','w'),fieldnames='locus old_locus type chr length start end strand operon product EC_number pfam clan COG pfam_notes pfam_data_notes other_notes'.split())
    w.writeheader()
    chrname='main_chromosome large_plasmid small_plasmid'.split()
    for ci,chr in enumerate(genome):
        for gene in chr.features:
            if 'locus_tag' in gene.qualifiers:
                data={'type':gene.type,
                    'locus':gene.qualifiers['locus_tag'][0],
                    'chr':chrname[ci],
                    'start':gene.location.start,
                    'end':gene.location.end,
                    'strand':gene.location.strand}
                if 'note' in gene.qualifiers:
                    data['product']=gene.qualifiers['note'][0]
                if gene.type=='CDS':
                    data['old_locus']=gene.qualifiers['old_locus_tag'][0]
                    data['product']=gene.qualifiers['product'][0]
                    data['length']=len(gene.qualifiers['translation'][0])
                    #assert not len(gene.qualifiers['translation'][0])-len(gene.extract(chr).seq.translate(to_stop=True)), 'length mismatch '+gene.qualifiers['locus_tag'][0]
                    if 'note' in gene.qualifiers:
                        #gene.qualifiers['note'].remove(gene.qualifiers['product'][0])
                        for n in gene.qualifiers['note']:
                            if n.find('Pfam data:')>-1:
                                data['pfam_data_notes']=n.replace('Pfam data:','')
                            elif n.find('COG: ')>-1:
                                data['COG']=n.replace('COG: ','')
                            elif n.find('Pfam: ')>-1:
                                data['pfam_notes']=n.replace('Pfam: ','')
                            else:
                                data['other_notes']=n
                    if 'EC_number' in gene.qualifiers:
                        data['EC_number']=gene.qualifiers['EC_number'][0]
                    if 'db_xref' in gene.qualifiers:
                        for d in gene.qualifiers['db_xref']:
                            if d.find('Pfam:PF')>-1:
                                data['pfam']=n.replace('Pfam:','')
                            elif d.find('Pfam:CL')>-1:
                                data['clan']=n.replace('Pfam:','')
                            else:
                                print('unknown? ',d)
                    if 'operon' in gene.qualifiers:
                        data['operon']=gene.qualifiers['operon'][0]
                w.writerow(data)




if __name__== "__main__":
    #draw()
    #"Gthg_TM242_v3.0.gb"
    #data=lister(get_genome("Gthg_TM242_v3.0.gb"))
    #w=csv.DictWriter(open("list_Gthg.csv","w"),"start locus end product".split())
    #w.writeheader()
    #w.writerows(data)

   # SeqIO.write(in_silico_frag(get_genome("LAKX01.1.gb"),1000),open("frag_LAKX01.1.fa","w"),"fasta")
    #lister(get_genome("Gthg_TM242_v3.0.gb"))
   #locusmap(get_genome("Gthg_TM242_v3.0.gb"))
    genome=get_genome("Gthg_TM242_v3.1.gb")
    alottable(genome)
    ## [2892463:2929663] brokenglu

    #Gthg03548 [2874505:2890826]

    g=genome[0][2125752:2134345]
    d={r['locus']:[r['locus'],math.log2(float(r['rpkm'])+1),0] for r in csv.DictReader(open('rpkm.csv','r'))}
    draw(g,d,"C:\\Users\\Cass\\Desktop\\Gthg02630-Gthg02637")
    #locusmap(get_genome('Gthg_TM242_v3.0.gb'))


    #framejoin(g,2135547,2136997,2138344,reversed=True)
    #print(joiner(g,2135547,2136997,2138344,reversed=True))

    #print(len(chr)) #3859975
    #SeqIO.write(chr[1:1031300],open('G_1.gb','w'),"genbank")
    #SeqIO.write(chr[1031300:1993000],open('G_2.gb','w'),"genbank")
    #SeqIO.write(chr[1993000:3859975],open('G_3.gb','w'),"genbank")