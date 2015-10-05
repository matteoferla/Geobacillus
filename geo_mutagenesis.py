__author__ = 'Matteo'
__doc__='''This could be made into a handy mutagenesis library if I had time.'''

from Bio.Seq import Seq,MutableSeq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from difflib import Differ

def Gthg01471():
    ori=Seq("ATGAGCATAAGTTTATCGGTTCCAAAATGGTTATTAACAGTTTTATCAATTTTATCTTTAGTCGTAGCATTTATTTTCGGTACCGTTTCCAATGCATCAGCAACAATTAACTATGGGGAGGAAGTCGCGGCAGTAGCAAATGACTATGTAGGAAGCCCATATAAATATGGAGGTACAACGCCAAAAGGATTTGATGCGAGTGGCTTTACTCAGTATGTGTATAAAAATGCTGCAACCAAATTGGCTATTCCGCGAACGAGTGCCGCACAGTATAAAGTCGGTAAATTTGTTAAACAAAGTGCGTTACAAAGAGGCGATTTAGTGTTTTATGCAACAGGAGCAAAAGGAAAGGTATCCTTTGTGGGAATTTATAATGGAAATGGTACGTTTATTGGTGCCACATCAAAAGGGGTAAAAGTGGTTAAAATGAGTGATAAATATTGGAAAGACCGGTATATAGGGGCTAAGCGAGTCATTAAGTAA", IUPAC.unambiguous_dna)
    mut=MutableSeq("ATGAGCATAAGTTTATCGGTTCCAAAATGGTTATTAACAGTTTTATCAATTTTATCTTTAGTCGTAGCATTTATTTTCGGTACCGTTTCCAATGCATCAGCAACAATTAACTATGGGGAGGAAGTCGCGGCAGTAGCAAATGACTATGTAGGAAGCCCATATAAATATGGAGGTACAACGCCAAAAGGATTTGATGCGAGTGGCTTTACTCAGTATGTGTATAAAAATGCTGCAACCAAATTGGCTATTCCGCGAACGAGTGCCGCACAGTATAAAGTCGGTAAATTTGTTAAACAAAGTGCGTTACAAAGAGGCGATTTAGTGTTTTATGCAACAGGAGCAAAAGGAAAGGTATCCTTTGTGGGAATTTATAATGGAAATGGTACGTTTATTGGTGCCACATCAAAAGGGGTAAAAGTGGTTAAAATGAGTGATAAATATTGGAAAGACCGGTATATAGGGGCTAAGCGAGTCATTAAGTAA", IUPAC.unambiguous_dna)

    a="AGTCGA"
    b="GACTAG"
    for i,v in enumerate([259,277,282,295,299,306]):
        print(mut[v-1]+a[i])
        mut[v-1]=b[i]
    print(ori.translate())
    print(mut.toseq().translate())

def Gthg04369():
    filepath="Gthg_from_embl_pfamed.gb"
    genome = list(SeqIO.parse(open(filepath, "rU"), "genbank"))
    z=genome[0].seq[3583975:3585290].translate(to_stop=1)
    x=genome[0].seq[3583975:3585290].tomutable()
    print(x.pop(895-1))
    y=x.toseq().translate(to_stop=1)
    print(z)
    print(y)
    print(list(Differ().compare(str(z),str(y))))
    print(len(z),len(y))

def Gthg01115():
    filepath="Gthg_from_embl_pfamed.gb"
    genome = list(SeqIO.parse(open(filepath, "rU"), "genbank"))
    z=genome[0].seq[891404:892205].reverse_complement().translate(to_stop=1)
    x=genome[0].seq[891404:892205].reverse_complement().tomutable()
    print(x.pop(421-1))
    y=x.toseq().translate(to_stop=1)
    print(z)
    print(y)
    print(list(Differ().compare(str(z),str(y))))
    print(len(z),len(y))

def Gthg03544():
    filepath="Gthg_from_embl_pfamed.gb"
    genome = list(SeqIO.parse(open(filepath, "rU"), "genbank"))
    z=genome[0].seq[2885410:2887572].reverse_complement().translate(to_stop=1)
    x=genome[0].seq[2885410:2887572].reverse_complement().tomutable()
    print(x.pop(1748-1))
    y=x.toseq().translate(to_stop=1)
    print(z)
    print(y)
    print(list(Differ().compare(str(z),str(y))))
    print(len(z),len(y))


if __name__ == "main":
    pass