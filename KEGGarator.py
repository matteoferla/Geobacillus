__author__ = 'Matteo'
__doc__ = '''Get data off SMBL and get data from KEGG.'''

N = "\n"
T = "\t"
# N="<br/>"
from urllib.request import urlopen
from urllib.error import HTTPError
import re
import pickle
import csv
from collections import OrderedDict
from itertools import chain

def keggmaster(rxn,main=True):
    try:
        url='http://rest.kegg.jp/get/'+rxn
        rep=urlopen(url).read()
        entry={'ID':{rxn}}
        last=0
        for line in rep.decode('UTF-8').split(N):
            parts=re.split(r'\s{2,}', line)
            print(parts)
            if len(parts)==1:
                pass
            elif parts[0] and parts[0] not in entry:
                entry[parts[0]]={' '.join(parts[1:])}
                last=parts[0]
            else:
                entry[last].add(' '.join(parts[1:]))
        if 'REMARK' in entry and main:
            for r in entry['REMARK']:
                if r.find('Same as: ')>-1:
                    target=r.replace('Same as: ','')
                    print('NESTED REQUEST...')
                    new=keggmaster(target,False)
                    for k in new:
                        if k=='REMARK':
                            pass
                        elif k in entry:
                            print(k,entry[k])
                            entry[k].update(new[k])
                        else:
                            entry[k]=new[k]
        return entry
    except HTTPError as e:
        print('*****')
        print(rxn)
        print(e)
        return {'ID':{rxn},'REMARK':str(e)}


def keggarate(todo):
    enzymeball=pickle.load(open('enzymeball.p','rb'))+[keggmaster(rxn) for rxn in todo]
    pickle.dump(enzymeball,open('enzymeball2.p','wb'))
    w=csv.DictWriter(open('enzymeball2.tsv','w'),fieldnames='ID PATHWAY NAME ENZYME COMMENT DEFINITION ENTRY MODULE ORTHOLOGY EQUATION REFERENCE RPAIR REMARK'.split(),delimiter=T )
    w.writeheader()
    w.writerows(enzymeball)

def keggexpander(todo):
    enzymeball=pickle.load(open('enzymeball.p','rb'))
    ids=[]
    for e in enzymeball:
        if 'ID' in e:
            ids.extend(e['ID'])
    print(ids)
    for t in todo:
        if t not in ids:
            print(t,'fetching')
            enzymeball.append(keggmaster(t))
            #pickle.dump(enzymeball,open('enzymeball_temp.p','wb'))
        else:
            print(t,'already got it')
    pickle.dump(enzymeball,open('enzymeball.p','wb'))
    w=csv.DictWriter(open('enzymeball.tsv','w'),fieldnames='ID PATHWAY NAME ENZYME COMMENT DEFINITION ENTRY MODULE ORTHOLOGY EQUATION REFERENCE RPAIR REMARK'.split(),delimiter=T )
    w.writeheader()
    w.writerows(enzymeball)


def counter():
    enzymeball=pickle.load(open('enzymeball3.p','rb'))
    paths=set()
    for e in enzymeball:
        if 'PATHWAY' in e:
            paths.update(e['PATHWAY'])
            e['PATHWAY']=sorted(e['PATHWAY'])
    for p in sorted(paths):
        print(p)
    print(len(paths))

def rxn2path(todo):
    d={}
    for e in pickle.load(open('enzymeball3.p','rb')):
        if e=={}:
            continue
        for e2 in e['ID']:
            if 'PATHWAY' in e:
                p=sorted(e['PATHWAY'])[0]
            else:
                p='NA'
            d[e2]=p
    for e in todo:
        if e in d:
            x=d[e]
        else:
            x='ERROR'
            print(e)
        print(x)

def locus2pathrix():
    __doc__ = '''Collage. I need to read new2old. Get the old key.
Get the sets of rxn keys associated with each old key.
Get the set of path keys associated with each rxn key.
Make a boolian matrix of path vs. new key.
This uses data from KEGGarator.py'''
    #make path dict
    path=dict()
    for p in pickle.load(open('enzymeball.p','rb')):
        if type(p) is not dict:
            print('!!!!!')
            print(p)
            print(type(p))
            continue
        if 'ID' not in p or 'PATHWAY' not in p:
            continue
        for i in p['ID']:
            if i in path:
                path[i].update(p['PATHWAY'])
            elif type(p['PATHWAY']) is set:
                path[i]=p['PATHWAY']
            else:
                path[i]={p['PATHWAY']}
    #Make pathlist
    plist=set()
    for e in path.values():
        plist.update(e)
    print(plist)
    #Make rxn dict
    rxn={}
    uber={}
    for r in csv.reader(open('fromSBML.csv','r')):
        if len(r)<2:
            continue
        if r[0] in rxn:
            rxn[r[0]].add(r[1])
        else:
            rxn[r[0]]={r[1]}
        if r[1] in path:
            if r[0] in uber:
                uber[r[0]].update(path[r[1]])
            else:
                uber[r[0]]=path[r[1]]
        else:
            uber[r[0]]=set()
    print(uber)
    w=csv.DictWriter(open('pathrix.csv','w'),fieldnames=['locus','old']+sorted(plist))
    w.writeheader()
    for line in open('new2old.txt','r').readlines():
        prot={}
        if not line:
            continue
        prot['locus'],prot['old']=line.split()
        if prot['old'] in uber:
            for f in uber[prot['old']]:
                prot[f]=1
        w.writerow(prot)

def SBML_looter(fp):
    import re,csv
    import xml.etree.ElementTree as ET
    geneball=[]
    tree = ET.parse(fp)
    root = tree.getroot()[0][3]
    for child in root:
        for tag in child[0]:
            #print(tag)
            if hasattr(tag,'text'):
                matches=re.search('(RTMO\d+)',tag.text)
                if matches:
                    for m in matches.groups():
                        geneball.append([m,child.attrib['name']])
    csv.writer(open('fromSBML.csv','w')).writerows(geneball)

def getR(ids):
    ball=dict()
    for p in pickle.load(open('enzymeball.p','rb')): #list of dict
        if type(p) is not dict:
            print('!!!!!')
            print(p)
            print(type(p))
            continue
        if 'ID' not in p:
            print('no id?',p)
            continue
        for i in p['ID']:
            ball[i]=p
            ball[i]['UID']=i
    w=csv.DictWriter(open('holes.csv','w'),dialect='excel',fieldnames='UID ID PATHWAY NAME ENZYME COMMENT DEFINITION ENTRY MODULE ORTHOLOGY EQUATION REFERENCE RPAIR REMARK'.split())
    w.writeheader()
    issue=[]
    for id in ids:
        if id in ball:
            w.writerow(ball[id])
        else:
            print(id,' invalid')
            issue.append(id)
    return issue

def allR():
    ball=dict()
    for p in pickle.load(open('enzymeball.p','rb')): #list of dict
        if type(p) is not dict:
            print('!!!!!')
            print(p)
            print(type(p))
            continue
        if 'ID' not in p:
            print('no id?',p)
            continue
        for i in p['ID']:
            ball[i]=p.copy()
            ball[i]['UID']=i
            if i=='R02421':
                print('YES',ball[i])
    w=csv.DictWriter(open('allRxn.csv','w'),dialect='excel',fieldnames='UID ID PATHWAY NAME ENZYME COMMENT DEFINITION ENTRY MODULE ORTHOLOGY EQUATION REFERENCE RPAIR REMARK'.split())
    w.writeheader()
    w.writerows(ball.values())

if __name__ == "__main__":
    #SBML_looter('Gthg_2.1.xml')
    #r=[x for x in open('rxn.txt','r').read().split() if x.find('R')>-1]
    #print('requests',r)
    #keggexpander(['R06087'])
    #i=getR(r)
    #allR()

    import json
    prot={}
    for line in open('new2old.txt','r').readlines():
        print(line)
        new,old=line.split()
        prot[old]=new
    print(prot)
    print(json.dumps(prot))