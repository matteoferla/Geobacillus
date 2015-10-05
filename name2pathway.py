__author__ = 'Matteo'
__doc__ = '''Given a list of gene product names match to a reaction ID and therefore to a pathway from BioCyc.
KEGG is not free and GO's flatfiles lack names. So BioCyc was chosen, which I like more than KEGG TBH and certainly GO.
It uses the enzymes.col and gene_association.metacyc flatfiles of metacyc. It has is mostly accurate.
It matches by words where if the match has longer than half matched words it is chosen.
I did not use Differ as I don't think that would have made too many false positives.'''

import csv, re
import numpy
#from collections import OrderedDict
N = "\n"
T = "\t"
# N="<br/>"

def wordmatching(mylist,refp,ref2p,fout):
    def make_symname(fp):
        symball={}
        for el in csv.reader(open(fp),delimiter=T):
            if len(el)<7:
                pass
            elif el[2] not in symball:
                symball[el[2]]=el[9]
        return symball

    def make_namepath(fp):
        nameball={}
        for line in open(fp).readlines():
            if line[0]=='#':
                continue
            #UNIQUE-ID	NAME	REACTION-EQUATION	PATHWAYS	PATHWAYS	PATHWAYS	PATHWAYS
            parts=line.split(T)
            name=' '.join(munger(parts[1]))
            paths=set(filter(None,parts[3:6]))
            if name not in nameball: #new key
                nameball[name]=paths
            else:  #seen key
                nameball[name].update(paths)
        return {k:sorted(nameball[k]) for k in nameball}


    def munger(txt):
        blacklist='transport binding family'.split()
        blacktype='transport tranport hypo'.split()
        for black in blacklist:
            if txt.find(black+' protein')>0:
                words=['CENSORED', blacktype[blacklist.index(black)]]
                break
        else:
            w=txt.lower().replace('protein',' ').replace('polypeptide',' ').replace('peptide',' ').replace('subunit',' ') #not family
            w=w.replace('oxo','keto').replace('yl','yl ') #off for now. Numbers?
            w=re.sub('\-[0-9]\-','',w)
            #w=w.replace('-','')
            words=re.sub('[^0-9a-zA-Z]+', ' ',w).split()
            blacklist='hypothetical putative family symporter transporter channel regulator factor repressor inducer'.split()
            blacktype='hypo hypo hypo transport  transport  transport TF TF TF TF'.split()
            for black in blacklist:#hypos ect will not be matched to a pathway.
                if black in words:
                    words=['CENSORED', blacktype[blacklist.index(black)]]
        return set(words)

    def matcher(queries,refs,sym):
        kwords=sorted(refs.keys(),key=lambda x: len(refs[x]))
        rpage=[set(k.split()) for k in kwords]
        score=[]
        censor=0
        hit=[['error'] for i in queries]  #???
        for qi,q in enumerate(queries):
            if len(q)<5:
                if q[1] in sym:
                    #print(q[1]+' is '+sym[q[1]])
                    q[1]=sym[q[1]]
            qwords=munger(q[1])
            if 'CENSORED' in qwords:
                censor+=censor
                hit[qi]=sorted(qwords-{'CENSORED'})
                continue
            max=0
            mi=None
            ms=999
            for ri,rwords in enumerate(rpage):
                inter=qwords.intersection(rwords)
                n=len(inter)
                s=len(' '.join(inter))
                if n>max or (n==max and s>ms):
                    max=n
                    mi=ri
                    ms=s
            score.append(max/len(qwords))
            if not mi:
                pass
                hit[qi]=['unmatched'] #dodgy
               # print(qwords,'unmatched')  #unmatched
            elif len(qwords) == len(rpage[mi]):
                hit[qi]=refs[kwords[mi]] #perferct
            elif max/len(qwords)>0.5:
                hit[qi]=refs[kwords[mi]] #dodgy
            else:
                hit[qi]=['unmatched'] #issue
        print(numpy.mean(score),score.count(0)/(len(queries)-censor),score.count(1)/(len(queries)-censor))
        return hit
        #0.295555237493 0.5176123180763552 0.15228854672010125 #hyphen=''
        #0.323181296193 0.47880194051887787 0.16199114110947058 #hyphen=' '
        #0.509206555418 0.11347816916262392 0.16199114110947058  skipping censored
        #0.54410711962 0.10693946424804894 0.17085003163889476
        #0.578467839238 0.08437038599451592 0.17274836532377136

    ref=make_namepath(refp)
    sym=make_symname(ref2p)
    qry=list(csv.reader(open(mylist),delimiter=T))
    ht=matcher(qry,ref,sym)
    w=csv.writer(open(fout,'w'))
    for i, qe in enumerate(qry):
        w.writerow(qe+ht[i])

if __name__ == "__main__":
    wordmatching('name_list 2.2.txt','enzymes.col','gene_association.metacyc','pairings.csv')
