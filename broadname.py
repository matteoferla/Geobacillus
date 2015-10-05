__author__ = 'Matteo'
__doc__ = ''''''

N = "\n"
T = "\t"
# N="<br/>"

from collections import OrderedDict
import csv

def make_symname(fp):
        symball={}
        for el in csv.reader(open(fp),delimiter=T):
            if len(el)<7:
                pass
            elif el[2] not in symball:
                symball[el[2]]=el[9]
        return symball




if __name__ == "__main__":
    sym=make_symname('gene_association.metacyc')
    keywords=OrderedDict([('transpos',0),
                          ('viral',0),
                          ('phage',0),
                          ('ase',1),
                          ('enzyme',1),
                          ('ribosomal',5),
                          ('biosynthesis',1),
                          ('export',3),
                          ('import',3),
                          ('symporter',3),
                          ('transport',3),
                          ('channel',3),
                          ('regulat',4),
                          ('factor',4),
                          ('repress',4),
                          ('induc',4),
                          ('translation',4),
                          ('hypothetical',2),
                          ('putative',2),
                          ('family',2)])
    terms=['transposon','enzyme','hypothetical','transporter','TF','ribosomal']
    for line in open('name_list 2.2.txt','r').readlines():
        if line:
            (locus,name)=line.rstrip().split(T)
            if len(name)<5:
                if name in sym:
                    name=sym[name]
                elif name[0]=='y':
                    name+='hypothetical'
            name=name.lower()
            for k in keywords:
                if name.find(k)>-1:
                    mapped=terms[keywords[k]]
                    break
            else:
                mapped='other'
                #print('--------',name)
        print(locus,mapped)
