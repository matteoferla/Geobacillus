__author__ = 'Matteo'
__doc__ = '''operate adds operon RockHopper data to GEnBank file. '''

N = "\n"
T = "\t"
# N="<br/>"
import csv, re, math
from Bio import SeqIO
from Bio import Seq
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Graphics import GenomeDiagram
from reportlab.lib import colors
from reportlab.lib.units import cm

def get_genome(fp):
    return  list(SeqIO.parse(open(fp, "rU"), "genbank"))

def genelooper(genome):
    return [gene for chr in genome for gene in chr.features if gene.type =="CDS"]


def operate(g,ci,fp):
    #Start	Stop	Strand	Number of Genes	Genes
    for ri,row in enumerate(csv.DictReader(open(fp,'r'),delimiter='\t')):
        op='Operon_p'+str(ci)+'_'+str(ri+1)
        for gi, gene in enumerate(g[ci].features):
            if gene.location.start>int(row['Start']):
                g[ci].features.insert(gi,SeqFeature(FeatureLocation(int(row['Start']), int(row['Stop'])), type="operon", strand=int(row['Strand']+'1'),qualifiers={'note':['Determined by RockHopper with AnRi52 dataset.'],'operon':[op]}))
                ni=gi+1
                mx=len(g[ci].features)
                while ni<mx and g[ci].features[ni].location.end <= int(row['Stop']):
                    g[ci].features[ni].qualifiers['operon']=[op]
                    ni+=1
                break
        else:
            print('Not matched',row)

def clean(gene):
    gene.qualifiers['product'][0]=re.sub(' \(EC .*?\)','',gene.qualifiers['product'][0])
    if 'note' in gene.qualifiers and gene.qualifiers['product'][0] in gene.qualifiers['note']:
        gene.qualifiers['note'].remove(gene.qualifiers['product'][0])
    if 'function' in gene.qualifiers and gene.qualifiers['product'][0] in gene.qualifiers['function']:
        gene.qualifiers['function'].remove(gene.qualifiers['product'][0])

def geneannotate(genome):
    #operate(genome,2,'small_plasmid_operons.txt')
    #operate(genome,1,'large_plasmid_operons.txt')
    deathlist=[gene[0] for gene in csv.reader(open('deathlist.txt','r'),delimiter='\t')]
    fusedict={fusion.split('-')[0]: fusion for fusion in open('fuse.txt','r').readlines()}
    skiplist=[]
    #if gene.qualifiers['locus_tag'][0] in fusedict:
    #                chr.features.insert(gi,SeqFeature(FeatureLocation(7, 110, strand=1), type="gene"))
    for gene in genelooper(genome):
        clean(gene)

    for fusion in open('fuse.txt','r').read().split():
        splits=fusion.split('-')
        print('hunting: ',splits)
        min=99999999999999999  #pretty sure they are sorted by appearance
        max=0
        mi=0
        cii=0
        strand=0
        prod=''
        geneball={}
        for ci,chr in enumerate(genome):
            for gi,gene in enumerate(chr.features):
                if gene.type == 'CDS':
                    if gene.qualifiers['locus_tag'][0] in splits:
                        #pseudotag
                        gene.qualifiers['pseudo']=None
                        #annotate
                        if not 'note' in gene.qualifiers:
                            gene.qualifiers['note']=[]
                        gene.qualifiers['note'].append('The gene appears to have been pseudogenised going from NCIB 11955 to TM242')
                        #add to skiplist
                        skiplist.append(gene.qualifiers['locus_tag'][0])
                        geneball[gene.qualifiers['locus_tag'][0]]=gene
                        #find ends
                        if gene.location.start < min:
                            min=gene.location.start
                            mi=gi
                        if gene.location.end > max:
                            max=gene.location.end
                        strand=gene.location.strand #by def they should all have the same
                        cii=ci
                        prod=gene.qualifiers['product'][0]
        genome[cii].features.insert(mi,SeqFeature(FeatureLocation(min, max, strand=strand), type="gene", qualifiers={'locus_tag':[fusion],'product':[prod], 'pseudo':None, 'note':['The gene appears to have been pseudogenised going from NCIB 11955 to TM242']}))


    for ci,chr in enumerate(genome):
        upshift=0 #coding with the flu is the worst
        for gi,gene in enumerate(list(chr.features)):
            if gene.type == 'CDS':
                print(gene.qualifiers['locus_tag'][0])
                if gene.qualifiers['locus_tag'][0] in skiplist:
                    print(gi, 'skiplist ',gene.qualifiers['locus_tag'][0])
                elif gene.qualifiers['locus_tag'][0] in deathlist:
                    print(gi, 'pseudogene ',gene.qualifiers['locus_tag'][0])
                    gene.qualifiers['pseudo']=None
                    if not 'note' in gene.qualifiers:
                        gene.qualifiers['note']=[]
                    gene.qualifiers['note'].append('No evidence of transcription in RNASeq data')
                    chr.features.insert(gi+upshift,SeqFeature(gene.location, type="gene",qualifiers={'locus_tag':gene.qualifiers['locus_tag'][0],'pseudo':None, 'Note':'No evidence of transcription in RNASeq data'}))
                    upshift+=1
                else:
                    print(gi, 'normal ',gene.qualifiers['locus_tag'][0])
                    chr.features.insert(gi+upshift,SeqFeature(gene.location, type="gene",qualifiers={'locus_tag':gene.qualifiers['locus_tag'][0]}))
                    upshift+=1

def list_possible_splits(genome): #variant of fusionmatch
    threshhold=20
    w=csv.writer(open('possibles.txt','w'),delimiter='\t')
    for chr in genome:
        previous=None
        for gene in chr.features:
            if gene.type=='CDS':
                if previous:
                    cond=False
                    if gene.location.strand==previous.location.strand:
                        if gene.qualifiers['product'][0]==previous.qualifiers['product'][0]:
                            cond=True
                            #print('match by product')
                        #if gene.location.start-previous.location.end<threshhold: #this time I need to be stringent.
                            #if gene.qualifiers['product'][0].find('ypothetic') >0 or gene.qualifiers['product'][0].find('ypothetic') >0:
                         #   cond=True
                                #print('match by hypo')
                        if cond:
                            if gene.location.strand>0:
                                w.writerow([previous.qualifiers['locus_tag'][0]+'-'+gene.qualifiers['locus_tag'][0],previous.location.end-previous.location.start,gene.location.start-previous.location.end,gene.location.end-gene.location.start, previous.qualifiers['product'][0],gene.qualifiers['product'][0]])
                            else:
                                w.writerow([gene.qualifiers['locus_tag'][0]+'-'+previous.qualifiers['locus_tag'][0],gene.location.end-gene.location.start,gene.location.start-previous.location.end,previous.location.end-previous.location.start, previous.qualifiers['product'][0],gene.qualifiers['product'][0]])
                previous=gene

                
def operate(genome,dataset,label, intensity=False): #go through and take snaps of each operon.
    #draw(g,d,"C:\\Users\\Cass\\Desktop\\Gthg03545-Gthg03544")
    for chr in genome:
        step=50000
        print(len(chr),' in ',step, ' is ', len(chr)/step)
        i=0
        for x in range(step+1,len(chr)-100,step):
            print(x-step,x+100)
            gs=chr[x-step:x+100]
            hyperdraw(chr[x-step:x+100],dataset,'C:\\Users\\Cass\\Desktop\\test\\'+label+gs.features[1].qualifiers['locus_tag'][0]+'-'+gs.features[-1].qualifiers['locus_tag'][0],intensity)
            i+=1


def hyperdraw(genome,dataset,label, intensity=False):
    for genoslice in genome:
        step=50000
        print(len(genoslice),' in ',step, ' is ', len(genoslice)/step)

        gd_diagram = GenomeDiagram.Diagram(genoslice.id)
        gd_track_for_features = gd_diagram.new_track(1,
                                                     name="Annotated Features",
                                                     scale_ticks=1,
                                                     scale_largetick_interval=1000,
                                                     scale_smalltick_interval=100,
                                                     scale_smallticks=0.05,
                                                     scale_largeticks=0.2,
                                                     scale_smalltick_labels=0

        )
        gd_feature_set = gd_track_for_features.new_set()
        for feature in genoslice.features:
            #if feature.type == "operon":
            #    gd_feature_set.add_feature(feature, sigil="BOX", ##pointy boxes
            #                           color=colors.grey, label=False)
            if feature.type =="rRNA" or feature.type =="tRNA":
                gd_feature_set.add_feature(feature, sigil="OCTO", ##pointy boxes
                                       color=colors.grey)
            elif feature.type =="CDS":
            # feature.qualifiers['product'][0].lower().find('hypo')>-1:
            # feature.qualifiers['product'][0].lower().find('transposase')>-1:
                locus=feature.qualifiers['locus_tag'][0]
                red,green, blue=(1,1,1)
                border=colors.gainsboro
                ah=0
                if locus in dataset:
                    if intensity:
                        intense=float(dataset[locus][1])
                        red,green=(1-intense/10,1-intense/10)
                    else:
                        ah=dataset[locus][3]
                        if ah<0:
                            ah=0
                        elif ah>1:
                            ah=1
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
                    labellable=True
                else:
                    labellable=False
                if len(feature)<200:
                    labellable=False
                color = colors.Color(red,green,blue)

                #shortened=feature.qualifiers['locus_tag'][0]+' '+feature.qualifiers['product'][0]
                #if len(shortened)>20:
                #    shortened=shortened[0:20]+'...'
                shortened=feature.qualifiers['product'][0]
                for badword in 'hypothetical putative probable family domain unknown possible partial'.split():
                    if feature.qualifiers['product'][0].lower().find(badword)>-1:
                        shortened=feature.qualifiers['locus_tag'][0]
                allocated=math.floor(len(feature)/1000*30)
                if len(shortened)>allocated:
                    shortened=shortened[0:allocated]+'...'

                gd_feature_set.add_feature(feature, sigil="ARROW",arrowshaft_height=ah, ##pointy boxes=1
                                       color=color, label=labellable,
                                       height=0.8, #not the actually setting...
                                       name=shortened,
                                       label_position="start",
                                       border=border,
                                       #label_strand=1,
                                       label_size = 8,
                                       label_angle=0
                                       )

        i=0
        for x in range(step+1,len(genoslice)-100,step):
            print(x-step,x+100)
            gs=genoslice[x-step:x+100]
            gd_diagram.draw(format="linear", pagesize='A4', orientation='landscape', fragments=10,
                    start=x-step, end=x+100)
            print(gs.features[1].qualifiers['locus_tag'][0],gs.features[-1].qualifiers['locus_tag'][0])
            op='C:\\Users\\Cass\\Desktop\\test2\\'+label+gs.features[1].qualifiers['locus_tag'][0]+'-'+gs.features[-1].qualifiers['locus_tag'][0]
            gd_diagram.write(op+".pdf", "PDF")
            gd_diagram.write(op+".eps", "EPS")
            gd_diagram.write(op+".svg", "SVG")
            gd_diagram.write(op+".png", "PNG")

if __name__ == "__main__":
    genome=get_genome('Gthg_TM242_v3.1.gb')
    #geneannotate(genome)
    #SeqIO.write(genome, open('Gthg_TM242_v3.1_test.gb', "w"), "genbank")
    #list_possible_splits(genome)
    #math.log2(float(r['RvM_log2FoldChange'])+1)
    #operate(genome,{r['locus']:[r['locus'],float(r['RvM_log2FoldChange']),float(r['RvM_padj'])] for r in csv.DictReader(open('change.csv','r'))},'Rich')
    #operate(genome,{r['locus']:[r['locus'],float(r['AvU_log2FoldChange']),float(r['AvU_padj'])] for r in csv.DictReader(open('change.csv','r'))},'Aero')
    #operate(genome,{r['locus']:[r['locus'],math.log2(float(r['rpkm'])+1),0] for r in csv.DictReader(open('rpkm.csv','r'))},'Int',intensity=True)
    #operate(genome,{r['locus']:[r['locus'],float(r['RvM_log2FoldChange']),float(r['RvM_padj']),math.log2(float(r['RvM_baseMean'])+1)/15] for r in csv.DictReader(open('change.csv','r'))},'Rich Int ')
    #operate(genome,{r['locus']:[r['locus'],float(r['AvU_log2FoldChange']),float(r['AvU_padj']),math.log2(float(r['AvU_baseMean'])+1)/15] for r in csv.DictReader(open('change.csv','r'))},'Aero Int ')
    #hyperdraw(genome,{r['locus']:[r['locus'],float(r['AvU_log2FoldChange']),float(r['AvU_padj']),math.log2(float(r['AvU_baseMean'])+1)/15] for r in csv.DictReader(open('change.csv','r'))},'Aero Int ')
    #hyperdraw(genome,{r['locus']:[r['locus'],float(r['RvM_log2FoldChange']),float(r['RvM_padj']),math.log2(float(r['RvM_baseMean'])+1)/15] for r in csv.DictReader(open('change.csv','r'))},'Rich Int ')




