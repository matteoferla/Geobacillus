__author__ = 'Matteo'
__doc__ = '''I need to change the locus tags on the SBML file.'''

N = "\n"
T = "\t"
# N="<br/>"

def update(fp,fp2):
    text=open(fp,'r').read()
    for line in open('new2old.txt','r').readlines():
        new,old=line.split()
        text=text.replace(old,new)
    import re
    text=re.sub('RTMO\d+','Gthg00000',text)
    open(fp2,'w').write(text)

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
                matches=re.search('(Gthg\d+)',tag.text)
                if matches:
                    for m in matches.groups():
                        geneball.append([m,child.attrib['name']])
    csv.writer(open('fromSBML.csv','w')).writerows(geneball)

def parse(fp):
    from libsbml import SBMLReader
    document = SBMLReader().readSBML(fp)
    print(document.getNumErrors())



if __name__ == "__main__":
    #update('depracated model.xml','depracated Gthg_model.xml')
    #SBML_looter('Depracated Gthg_model.xml')
    parse('Gthg_2.2 (manual).xml')

