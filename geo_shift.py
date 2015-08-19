__author__ = 'Matteo'
__doc__=''' There seems to be no standard locus names for Geobacillus thermoglucosidasius. So I went with gthm_
C56-YS93  goes GEOTH_RS
TNO-09.020 goes GT20_RS
NBRC 107763 goes GT2_01_
DSM 2542 goes WH82_

Before I went ahead I did some checks.
contig_check. How many bits? 1 chr and two plasmids
ori_check. Where does GC skew say the Ori is? base 1.
roundtrip_check. GenBank unusual coding that will be lost? no.
sequenciality_check. Is everything aligned? well. 499 overlaps! Eco has 345.
shifter. shift and clean.

'''



def roundtrip_check(filepath):
    from Bio import SeqIO
    import os
    input_handle = open("TMO.gbk", "rU")
    SeqIO.convert(filepath, "genbank", "test.gbk", "genbank")
    print("went from "+str(os.stat('TMO.gbk').st_size)+" to "+str(os.stat("test.gbk").st_size))
    os.remove("test.gbk")

def contig_check():
    __doc__='''
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
    from Bio import SeqIO
    from collections import Counter
    os.chdir("genome_shift")
    input_handle = open("TMO.gbk", "rU")
    x=0
    for record in SeqIO.parse(input_handle, "genbank") :
        print(record.id+" "+str(len(record.features)))
        tally=Counter([gene.type for gene in record.features])
        x+=tally['gene']
        print("    genes: "+str(tally['gene']))
    print("total genes: "+str(x))
    input_handle.close()



def ori_check(filepath):
    from Bio import SeqIO
    #from Bio.SeqUtils import xGC_skew
    input_handle = open(filepath, "rU")
    from itertools import tee
    genome=SeqIO.parse(input_handle, "genbank")
    g1, g2=tee(genome)
    x=1
    element=0
    for chr in g1:
        chr.id=["Chromosome","pGTH1","pGTH2"][element]
        element+=1
    for chr in g2:
        #print(chr.annotations)
        for i in range(10):
            #print(chr.features[i])
            pass
        #xGC_skew(chr.seq,50)
        #crashes?
        xGC_skew_mod(chr.id, chr.seq,str(chr.id)+" (GC skew plot).ps")
        element+=1

def xGC_skew_mod(name, seq, outpath, window=None, zoom=100,r=300, px=100, py=100):
    """Calculates and plots normal and accumulated GC skew (GRAPHICS !!!)."""
    #Copy of from Bio.SeqUtils import xGC_skew. to make it print
    if not window:
        window=round(len(seq)/720)
    from Bio.SeqUtils import GC, GC_skew
    from math import pi, sin, cos
    try:
        import Tkinter as tkinter # Python 2
    except ImportError:
        import tkinter # Python 3

    yscroll = tkinter.Scrollbar(orient=tkinter.VERTICAL)
    xscroll = tkinter.Scrollbar(orient=tkinter.HORIZONTAL)
    canvas = tkinter.Canvas(yscrollcommand=yscroll.set,xscrollcommand=xscroll.set, background='white')
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
        alpha = pi - (2*pi*start)/len(seq)
        r2 = r1 - gc*zoom
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
    canvas.postscript(file= outpath)
    canvas.destroy()

def sequenciality_check(filepath):
    __doc__='''ERROR: ['RTMO01197'] [12893:13148](+) comes after ['RTMO01196'] [12658:12919](+)offset: -26
ERROR: ['RTMO01188'] [23286:24321](+) comes after ['RTMO01185'] [21436:23290](+)offset: -4
ERROR: ['RTMO01189'] [24317:24866](+) comes after ['RTMO01188'] [23286:24321](+)offset: -4
ERROR: ['RTMO01175'] [30384:31023](+) comes after ['RTMO01174'] [30205:30388](+)offset: -4
ERROR: ['RTMO01176'] [31250:32321](+) comes after ['RTMO05742'] [31102:31282](+)offset: -32'''
    from Bio import SeqIO
    #from Bio.SeqUtils import xGC_skew
    input_handle = open(filepath, "rU")
    from itertools import tee
    genome=SeqIO.parse(input_handle, "genbank")
    g1, g2=tee(genome)
    x=1
    element=0
    for chr in g1:
        chr.id=["Chromosome","pGTH1","pGTH2"][element]
        element+=1
    for chr in g2:
        print("*"*20)
        print(chr.id)
        min=0
        max=0
        historic="initialised"
        for gene in chr.features:
            if 'locus_tag' in gene.qualifiers:
                #print(str(gene.location)+" "+str(gene.qualifiers['locus_tag']))
                if gene.location.strand ==1:
                    new_min=gene.location.start.position
                    new_max=gene.location.end.position
                else:
                    new_max=gene.location.start.position
                    new_min=gene.location.end.position

                if new_min < max and new_min !=min:
                    x+=1
                    print("ERROR: "+str(gene.qualifiers['locus_tag'])+" "+str(gene.location)+" comes after "+historic+  "offset: "+str(new_min-max))
                min=new_min
                max=new_max
                historic=str(gene.qualifiers['locus_tag'])+" "+str(gene.location)
    print("="*20)
    print("number of overlaps: "+str(x))


def shifter(filepath,outpath):
    bug="Geobacillus thermoglucosidasius NCIMB 11955"
    from Bio import SeqIO
    genome=list(SeqIO.parse(open(filepath, "rU"), "genbank"))
    x=1
    element=0
    for chr in genome:
        chr.name=["Chromosome","pGTH1","pGTH2"][element]
        chr.id="NC_xxxxxx"
        chr.description=[bug+" complete "+z for z in ["chromosome","large plasmid (of two)","small plasmid (of two)"]][element]
        chr.annotations['references']=[]
        chr.annotations['source']=bug
        chr.annotations['organism']=bug
        chr.annotations['taxonomy']=['Bacteria','Firmicutes','Bacilli','Bacillales','Bacillaceae','Geobacillus']
        chr.annotations['sequence_version']="v.2"
        element+=1
    for chr in genome:
        for gene in chr.features:
            if gene.type =="source":
                gene.qualifiers['strain']="NCIMB 11955"
                gene.qualifiers['organism']=bug
            if 'locus_tag' in gene.qualifiers:
                new="Gthm"+"{:0>5d}".format(x)
                #x++ feck
                x+=1
                gene.qualifiers['old_locus_tag']=gene.qualifiers['locus_tag']
                gene.qualifiers['locus_tag']=[new]
                # 'protein_id': ['igergo:RTMO05622'], should protein_id="NP_414551.1"
                gene.qualifiers['protein_id']=["NP_xxxxxx"]
    SeqIO.write(genome,open(outpath, "w"),"genbank")

if __name__=="__main__":
    ref="Eco NC_000913.gb"
    geo="TMO.gbk"
    fp=geo
    sequenciality_check(fp)
    #precheck(fp)
    #shifter(fp,"Gth_NCIMB_11955_v2.gb")

    #Eco NC_000913
