__author__ = 'Matteo'
__doc__ = ''' There seems to be no standard locus names for Geobacillus thermoglucosidasius. So I went with gthm_
C56-YS93  goes GEOTH_RS
TNO-09.020 goes GT20_RS
NBRC 107763 goes GT2_01_
DSM 2542 goes WH82_

Before I went ahead I did some checks.
contig_check. How many bits? 1 chr and two plasmids
ori_check. Where does GC skew say the Ori is? base 1.
roundtrip_check. GenBank unusual coding that will be lost? no.
sequentiality_check. Is everything aligned? well. 499 overlaps! Eco has 345.
shifter. shift and clean.

'''
from Bio import SeqIO

def roundtrip_check(filepath):
    import os
    input_handle = open("TMO.gbk", "rU")
    SeqIO.convert(filepath, "genbank", "test.gbk", "genbank")
    print("went from " + str(os.stat('TMO.gbk').st_size) + " to " + str(os.stat("test.gbk").st_size))
    os.remove("test.gbk")


def contig_check():
    __doc__ = '''
Contig0003 9256
    genes: 4569
Contig0001 227
    genes: 113
Contig0002 119
    genes: 59
total genes: 4741
while http://www.ncbi.nlm.nih.gov/genome/2405 says:
Type	Name	RefSeq	INSDC	Size (Mb)	GC%	Protein	rRNA	tRNA	Other RNA	Gene	Pseudogene
Chr	-	NC_015660.1	CP002835.1	3.89	44.0	3,630	27	90	1	3,895	118
Plsm	pGEOTH01	NC_015665.1	CP002836.1	0.080849	44.0	67	-	-	-	77	6
Plsm	pGEOTH02	NC_015661.1	CP002837.1	0.019638	40.4	20	-	-	-	20	-
Oh. three bits. Doh.


So Contig 3 > Chr
Contig 1 > bigger plasmid
Contig 2 > smaller plasmid
'''
    from collections import Counter

    os.chdir("genome_shift")
    input_handle = open("TMO.gbk", "rU")
    x = 0
    for record in SeqIO.parse(input_handle, "genbank"):
        print(record.id + " " + str(len(record.features)))
        tally = Counter([gene.type for gene in record.features])
        x += tally['gene']
        print("    genes: " + str(tally['gene']))
    print("total genes: " + str(x))
    input_handle.close()


def ori_check(filepath):
    # from Bio.SeqUtils import xGC_skew
    input_handle = open(filepath, "rU")
    from itertools import tee

    genome = SeqIO.parse(input_handle, "genbank")
    g1, g2 = tee(genome)
    x = 1
    element = 0
    for chr in g1:
        chr.id = ["Chromosome", "pGTH1", "pGTH2"][element]
        element += 1
    for chr in g2:
        #print(chr.annotations)
        for i in range(10):
            #print(chr.features[i])
            pass
        #xGC_skew(chr.seq,50)
        #crashes?
        xGC_skew_mod(chr.id, chr.seq, str(chr.id) + " (GC skew plot).ps")
        element += 1


def xGC_skew_mod(name, seq, outpath, window=None, zoom=100, r=300, px=100, py=100):
    """Calculates and plots normal and accumulated GC skew (GRAPHICS !!!)."""
    # Copy of from Bio.SeqUtils import xGC_skew. to make it print
    if not window:
        window = round(len(seq) / 720)
    from Bio.SeqUtils import GC, GC_skew
    from math import pi, sin, cos

    try:
        import Tkinter as tkinter  # Python 2
    except ImportError:
        import tkinter  # Python 3

    yscroll = tkinter.Scrollbar(orient=tkinter.VERTICAL)
    xscroll = tkinter.Scrollbar(orient=tkinter.HORIZONTAL)
    canvas = tkinter.Canvas(yscrollcommand=yscroll.set, xscrollcommand=xscroll.set, background='white')
    win = canvas.winfo_toplevel()
    win.geometry('900x900')

    yscroll.config(command=canvas.yview)
    xscroll.config(command=canvas.xview)
    yscroll.pack(side=tkinter.RIGHT, fill=tkinter.Y)
    xscroll.pack(side=tkinter.BOTTOM, fill=tkinter.X)
    canvas.pack(fill=tkinter.BOTH, side=tkinter.LEFT, expand=1)
    canvas.update()

    X0, Y0 = r + px, r + py
    x1, x2, y1, y2 = X0 - r, X0 + r, Y0 - r, Y0 + r

    ty = Y0
    canvas.create_text(X0, ty, text=name)
    ty += 20
    canvas.create_text(X0, ty, text='GC %3.2f%%' % (GC(seq)))
    ty += 20
    canvas.create_text(X0, ty, text='GC Skew', fill='blue')
    ty += 20
    canvas.create_text(X0, ty, text='Accumulated GC Skew', fill='magenta')
    ty += 20
    canvas.create_oval(x1, y1, x2, y2)

    acc = 0
    start = 0
    for gc in GC_skew(seq, window):
        r1 = r
        acc += gc
        # GC skew
        alpha = pi - (2 * pi * start) / len(seq)
        r2 = r1 - gc * zoom
        x1 = X0 + r1 * sin(alpha)
        y1 = Y0 + r1 * cos(alpha)
        x2 = X0 + r2 * sin(alpha)
        y2 = Y0 + r2 * cos(alpha)
        canvas.create_line(x1, y1, x2, y2, fill='blue')
        # accumulated GC skew
        r1 = r - 50
        r2 = r1 - acc
        x1 = X0 + r1 * sin(alpha)
        y1 = Y0 + r1 * cos(alpha)
        x2 = X0 + r2 * sin(alpha)
        y2 = Y0 + r2 * cos(alpha)
        canvas.create_line(x1, y1, x2, y2, fill='magenta')

        canvas.update()
        start += window

    canvas.configure(scrollregion=canvas.bbox(tkinter.ALL))
    canvas.postscript(file=outpath)
    canvas.destroy()


def sequentiality_check(filepath):
    __doc__ = '''ERROR: ['RTMO01197'] [12893:13148](+) comes after ['RTMO01196'] [12658:12919](+)offset: -26
ERROR: ['RTMO01188'] [23286:24321](+) comes after ['RTMO01185'] [21436:23290](+)offset: -4
ERROR: ['RTMO01189'] [24317:24866](+) comes after ['RTMO01188'] [23286:24321](+)offset: -4
ERROR: ['RTMO01175'] [30384:31023](+) comes after ['RTMO01174'] [30205:30388](+)offset: -4
ERROR: ['RTMO01176'] [31250:32321](+) comes after ['RTMO05742'] [31102:31282](+)offset: -32'''

    # from Bio.SeqUtils import xGC_skew
    input_handle = open(filepath, "rU")
    from itertools import tee

    genome = SeqIO.parse(input_handle, "genbank")
    g1, g2 = tee(genome)
    x = 1
    element = 0
    for chr in g1:
        chr.id = ["Chromosome", "pGTH1", "pGTH2"][element]
        element += 1
    for chr in g2:
        print("*" * 20)
        print(chr.id)
        min = 0
        max = 0
        historic = "initialised"
        for gene in chr.features:
            if 'locus_tag' in gene.qualifiers:
                #print(str(gene.location)+" "+str(gene.qualifiers['locus_tag']))
                new_min = gene.location.start.position
                new_max = gene.location.end.position

                if new_min < max and new_min != min:
                    x += 1
                    print("ERROR: " + str(gene.qualifiers['locus_tag']) + " " + str(
                        gene.location) + " comes after " + historic + "offset: " + str(new_min - max))
                min = new_min
                max = new_max
                historic = str(gene.qualifiers['locus_tag']) + " " + str(gene.location)
    print("=" * 20)
    print("number of overlaps: " + str(x))


def shifter(filepath, outpath, affix):
    bug = "Geobacillus thermoglucosidasius NCIMB 11955"

    genome = list(SeqIO.parse(open(filepath, "rU"), "genbank"))
    x = 1
    element = 0
    genome[0].features[0].qualifiers['chromosome']="primary"
    genome[1].features[0].qualifiers['plasmid']="p"+affix+"1"
    genome[2].features[0].qualifiers['plasmid']="p"+affix+"2"
    for chr in genome:
        chr.name = ["Chromosome", "p"+affix+"1", "p"+affix+"2"][element]
        chr.id = "NC_xxxxxx"
        chr.description = \
        [bug + " complete " + z for z in ["chromosome", "large plasmid (of two)", "small plasmid (of two)"]][element]
        chr.annotations['references'] = []
        chr.annotations['source'] = bug
        chr.annotations['organism'] = bug
        chr.annotations['taxonomy'] = ['Bacteria', 'Firmicutes', 'Bacilli', 'Bacillales', 'Bacillaceae', 'Geobacillus']
        chr.annotations['sequence_version'] = "v.2"
        element += 1
    for chr in genome:
        for gene in chr.features:
            if gene.type == "source":
                gene.qualifiers['strain'] = "NCIMB 11955 (KP 1006 = R-35637 = ATCC 43742 = BGSC 95A1 = CCUG 28887 = CIP 106930 = DSM 2542 = LMG 7137 = NCIMB 11955 = NRRL B-14516)"
                gene.qualifiers['organism'] = "Geobacillus thermoglucosidasius"
                gene.qualifiers['mol_type']="genomic DNA"
                gene.qualifiers['db_xref']="taxon:1426"
                gene.qualifiers['sub_clone']="laboratory collection of Leak, DJ, University of Bath"
            elif gene.type == "tRNA" or gene.type == "rRNA":
                new = "Gtg" + "{:0>5d}".format(x)
                gene.qualifiers['locus_tag'] = [new]
                x += 1
            elif gene.type == "gene":
                new = affix + "{:0>5d}".format(x)
                gene.qualifiers['old_locus_tag'] = gene.qualifiers['locus_tag']
                gene.qualifiers['locus_tag'] = [new]
                # 'protein_id': ['igergo:RTMO05622'], should protein_id="NP_414551.1"
                # x++ feck
                x += 1
            elif gene.type == "CDS":
                gene.qualifiers['old_locus_tag'] = gene.qualifiers['locus_tag']
                gene.qualifiers['locus_tag'] = [new]
                gene.qualifiers['protein_id'] = ["NP_xxxxxx"]
            else:
                print("UNKNOWN FEATURE TYPE: " + str(gene.type))

    SeqIO.write(genome, open(outpath, "w"), "genbank")


def botch(filepath, v):
    f = open(filepath, "rU")
    txt = f.read()
    txt=txt.replace("VERSION     NC_xxxxxx", "VERSION     " + str(v)).replace("DNA              ", "DNA     circular ")
    f.close()
    open(filepath, "w").write(txt)

def fastidiate(filepath): #How punderful! I am fastidious about the fasta header.
    N="\n"
    T="\t"
    protball=[]
    for chr in SeqIO.parse(filepath,"genbank"):
        protball+=[">"+gene.qualifiers["locus_tag"][0]+N+gene.qualifiers['translation'][0] for gene in chr.features if gene.type=="CDS"]
        #Nice that there is a translation qualifier
    open(filepath.replace(".gbk",".fasta").replace(".gb",".fasta").replace(".gen",".fasta"),"w").write(N.join(protball)+N)
    #Feck. Python, I hate your backwards join (and your lack of autoincrementor and conditional assignment).
    #Annoyme again and I will switch to whitespace

def pfaming(filepath, pfilepath, outpath):
    N="\n"
    T="\t"
    #print([line.split() for line in open(pfilepath,"r").readlines()][1])
    #['<seq', 'id>', '<alignment', 'start>', '<alignment', 'end>', '<envelope', 'start>', '<envelope', 'end>', '<hmm', 'acc>',
    # '<hmm', 'name>', '<type>', '<hmm', 'start>', '<hmm', 'end>', '<hmm', 'length>', '<bit', 'score>', '<E-value>',
    # '<significance>', '<clan>']
    fields=['id','al_s','al_e','en_s','en_e','hmm_id','hmm_name','type','hmm_s','hmm_e','hmm_l','bit','evalue','sig','clan']
    #pfam=[dict(zip(fields,line.split())) for line in open(pfilepath,"r").readlines()][2:]
    #Pythonic bliss. However I try it never gets satanic as Perl.
    #damn. I had change approach
    pfam={}
    f=open(pfilepath,"r")
    next(f)
    next(f)
    for line in f:
        l=dict(zip(fields,line.split()))
        pfam[l['id']]=l
    #Regret in not making a lambda, coucld a dict comprehension accept a tuple returned by a lambda?
    genome = list(SeqIO.parse(open(filepath, "rU"), "genbank"))
    for chr in genome:
        for gene in chr.features:
            if gene.type=="CDS":
                locus=gene.qualifiers['locus_tag'][0]
                if locus in pfam.keys():
                    #db_xref="<database>:<identifier>"
                    if not 'db_xref' in gene.qualifiers:
                        gene.qualifiers['db_xref']=[]
                    gene.qualifiers['db_xref'].append("Pfam:"+pfam[locus]['hmm_id'])
                    if pfam[locus]['clan'] !='No_clan':
                        gene.qualifiers['db_xref'].append("Pfam:"+pfam[locus]['clan'])
                    if not 'note' in gene.qualifiers:
                        gene.qualifiers['note']=[]
                    gene.qualifiers['note'].append("Pfam data:"+", ".join(pfam[locus].values()))
                else:
                    print("no "+locus)
    SeqIO.write(genome, open(outpath, "w"), "genbank")

def more_tweaks():
    pass

def nominate(filepath,outpath):
    N="\n"
    T="\t"
    nameball=[gene.qualifiers["locus_tag"][0]+T+gene.qualifiers['product'][0] for chr in SeqIO.parse(filepath,"genbank") for gene in chr.features if gene.type=="CDS"]
    #incomprehensible list comprehensions make me happy
    open(outpath,"w").write(N.join(nameball))


def length_table(filepath):
    #copy paste from sequentiality
    space=[]
    genome = list(SeqIO.parse(open(filepath, "rU"), "genbank"))
    for gene in genome[0].features:
        if gene.type=="CDS" or gene.type=="tRNA" or gene.type=="rRNA":
            continue
        if 'locus_tag' in gene.qualifiers:
            #print(str(gene.location)+" "+str(gene.qualifiers['locus_tag']))
            space.append([str(gene.qualifiers['locus_tag'][0]),gene.location.start.position, gene.location.end.position])
    print(len(genome[0]))
    return space


if __name__ == "__main__":
    ref = "Eco NC_000913.gb"
    geo = "TMO.gbk"
    fp = geo
    # sequentiality_check(fp)
    #precheck(fp)
    nfp="Gth_NCIMB_11955_v2.3.gb"
    #shifter(fp, nfp,"Gthg")
    #botch(nfp,2.2)
    #fastidiate(nfp)
    #nominate(nfp,"name_list 2.2.txt")
    geo="Gth_NCIMB_11955_v2.4.gb"
    #pfaming(nfp, "pfam2.2.txt",n2fp)

    #Eco NC_000913
    sequentiality_check(geo)
    import csv
   # csv.writer(open("Geo_out.txt","w")).writerows([['locus', 'start','end']]+length_table(geo))