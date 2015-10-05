#! python
__author__ = 'Matteo'
__doc__ = '''
I would like to add KEGG data to the notes.
I would like to add reactions.

Unfortunately... Adding species is easy, except for a few dodgy cases. The code is over tampered with. It works fine for most bar them.


'''

N = "\n"
T = "\t"
# N="<br/>"

import libsbml, pickle, re, csv
from collections import Counter, defaultdict
from urllib.request import urlopen
from urllib.error import HTTPError

def SBML_looter(model):  #rewritten
    geneball=[]
    for rxn in model.getListOfReactions():
        matches=re.search('(Gthg\d+)',rxn.getNotesString())
        if matches:
            for m in matches.groups():
                geneball.append([m, rxn.getName()])
    csv.writer(open('fromSBML.csv','w')).writerows([['locus','rxn']]+geneball)


def generow(model):
    nameball=defaultdict(list)
    pathball=defaultdict(list)
    chemball=defaultdict(list)
    for rxn in model.getListOfReactions():
        pathex=re.search('PATHWAY: (.*?)<',rxn.getNotesString())
        if pathex:
            path=pathex.groups(1)
        else:
            path=''
        chemex=re.search('DEFINITION: (.*?)<',rxn.getNotesString())
        if chemex:
            chem=chemex.groups(1)
        else:
            chem=''
        genes=re.search('(Gthg\d+)',rxn.getNotesString())
        if genes:
            for g in genes.groups():
                nameball[g].append(rxn.getName())
                pathball[g].append(path)
                chemball[g].append(chem)

    csv.writer(open('genesfromSBML.csv','w',encoding='utf32'),delimiter='\t').writerows([['locus','rxns','pathways','defs']]+[[gene,nameball[gene],pathball[gene],chemball[gene]] for gene in nameball])

def check(value, message): #COPY PASTE FROM MANUAL with the except that I added libsbml. as I did not import all.
   """If 'value' is None, prints an error message constructed using
   'message' and then exits with status code 1.  If 'value' is an integer,
   it assumes it is a libSBML return status code.  If the code value is
   LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
   prints an error message constructed using 'message' along with text from
   libSBML explaining the meaning of the code, and exits with status code 1.
   """
   if value == None:
     raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
   elif type(value) is int:
     if value == libsbml.LIBSBML_OPERATION_SUCCESS:
       return
     else:
       err_msg = 'Error encountered trying to ' + message + '.' \
                 + 'LibSBML returned error code ' + str(value) + ': "' \
                 + libsbml.OperationReturnValue_toString(value).strip() + '"'
       raise SystemExit(err_msg)
   else:
     return

def get_kegg_data():
    keggtionary={}
    for k in pickle.load(open('enzymeball.p','rb')):
        if 'ID' in k:
            for i in k['ID']:
                    keggtionary[i]=k.copy()
    return keggtionary

def keggmaster(rxn,main=True):
    try:
        url='http://rest.kegg.jp/get/'+rxn
        rep=urlopen(url).read()
        entry={'ID':{rxn}}
        last=0
        for line in rep.decode('UTF-8').split(N):
            parts=re.split(r'\s{2,}', line)
            #print(parts)
            if len(parts)==1:
                pass
            elif parts[0] and parts[0] not in entry:
                entry[parts[0]]={' '.join(parts[1:])}
                last=parts[0]
            else:
                entry[last].add(' '.join(parts[1:]))
        print(entry)
        return entry
    except HTTPError as e:
        print('*****')
        print(rxn)
        print(e)
        return {'ID':{rxn},'REMARK':str(e)}

def hyperwinner(hypercount,bij=defaultdict(str),n=[999999]):
#def hyperwinner(hypercount,bij=defaultdict(set),n=[999999]):
    __doc__='This is a mess, because I tried to automate something that really could not. It does 90% of the job.' \
            'You\'ll have to figure out the rest.'
    print(len(bij))
    if n[0]==len(hypercount):
        #salvage... repeat through each iteration...
        #for k in list(hypercount.keys()):
            #for cand in hypercount[k]:
                #print('cadn ',cand,'k ',hypercount[k],bij[cand[0][0]])
                #if type(bij[cand[0][0]]) is set:
                    #bij[cand[0][0]].add(k) #die if a string!
                #print(hypercount[k])
        return bij
    else:
        n[0]=len(hypercount)
    for k in list(hypercount.keys()):
        del hypercount[k]['Σ']
        #del hypercount[k]['Σ']
        c=hypercount[k].most_common(2)
        #assert c[0][1] == Σ,'Max less than sigma for '+c[0][0]

        #assert c[0][1] > c[1][1],'DRAW!!'+str(c)
        if len(c)==0:
            pass
            #print('>>>',k,hypercount[k])
        elif len(c)==1 or c[0][1] > c[1][1]:
            #assert c[0][0] not in bij, 'Double match for '+c[0][0]
            bij[c[0][0]]=k
            del hypercount[k]
            for k2 in hypercount.keys():
                del hypercount[k2][c[0][0]]
        else:
            pass
    if hypercount:
        hyperwinner(hypercount)
    return bij

def rxn_kegg2sbml():
    kict=get_kegg_data()
    sbml=libsbml.SBMLReader().readSBML("Gthg_2.2 (manual).xml")
    print('The model has ',sbml.getNumErrors(),' errors')
    model=sbml.getModel()
    hypercount=defaultdict(Counter)
    for rxn in model.getListOfReactions():
        if rxn.getName() in kict:
            #Also should add: to RDF
            #http://www.kegg.jp/entry/R01090
            details=kict[rxn.getName()]
            for f in sorted(details):
                if f !='ENTRY':
                    #How do I encode <=>?
                    #Unicode? HTML ascii?
                    #https://en.wikipedia.org/wiki/Arrow_(symbol)
                    #encode('ascii','xmlcharrefreplace'))
                    rxn.appendNotes('<html:p>'+f+': '+'; '.join(details[f]).replace("<=>","\u21cc").replace('=>','\u2192').replace('<=','\u2190')+'</html:p>')
            #if 'DEFINITION' in details and 'EQUATION' in details:
            if 'EQUATION' in details:
                #Match the species.
                #prods=[model.getElementBySId(rxn.getProduct(n).getSpecies()) for n in range(rxn.getNumProducts())]
                #reacs=[model.getElementBySId(rxn.getReactant(n).getSpecies()) for n in range(rxn.getNumReactants())]
                #Split the def and eq lines and zip then
                #terms=list(zip([[re.search("([CG]\d+)",s2).group(0) for s2 in re.sub('\(.*?\)','',s).split('+')] for s in re.split('<?=>?',list(details['EQUATION'])[0])],[s.split('+') for s in re.split('<?=>?',list(details['DEFINITION'])[0])]))
                for k in [re.search("([CG]\d+)",s2).group(0) for s in re.split('<?=>?',list(details['EQUATION'])[0]) for s2 in re.sub('\(.*?\)','',s).split('+')]:
                    hypercount[k].update(['Σ']+[rxn.getProduct(n).getSpecies() for n in range(rxn.getNumProducts())]+[rxn.getReactant(n).getSpecies() for n in range(rxn.getNumReactants())])
        else:
            print('Unknown ',rxn.getName())
    #print(hypercount)
    mapID=hyperwinner(hypercount)
    print('data:',mapID)
    for sp in model.getListOfSpecies():
        tag=sp.getId().replace('_e','_c')
        if tag in mapID:
            if type(mapID[tag]) is str:
                sp.appendNotes('<html:p>KEGG: '+mapID[tag]+'</html:p>')
            #else:
               # sp.appendNotes('<html:p>KEGG: either '+' or '.join(mapID[tag])+'</html:p>')
        else:
            sp.appendNotes('<html:p>KEGG: !!!!GOGGLE IT</html:p>')
            print(tag,' not an id element')
    libsbml.writeSBMLToFile(sbml,"test enrich.xml")

def error_check(sbml):
    print('The model has ',sbml.getNumErrors(),' errors')
    for n in range(sbml.getNumErrors()):
        e=sbml.getError(n)
        print('cat',e.getCategory())
        print('issue',e.getMessage())
        print('line',e.getLine())

def notify(cid, keggball={}):
    newotes='            <body xmlns="http://www.w3.org/1999/xhtml">\n            <p>KEGG: '+cid+'</p>\n'
    if cid in keggball:
        kict=keggball[cid]
    else:
        kict=keggmaster(cid)
    like='NAME MODULE FORMULA EXACT_MASS'.split()
    for l in like:
        if l in kict:
            newotes+='            <p>'+l+': '+"; ".join(kict[l])+'</p>\n'
    if 'COMMENT' in kict:
        newotes+='            <p>KEGG_COMMENT: '+"; ".join(kict['COMMENT'])+'</p>\n'
    if 'DBLINKS' in kict:
        for link in kict['DBLINKS']:
            newotes+='            <p>'+link+'</p>\n'
    #print(newotes)
    return newotes+'            </body>'

def sbml_addkeggchem():
    kict=get_kegg_data()
    sbml=libsbml.SBMLReader().readSBML("Gthg_2.3.xml")
    error_check(sbml)
    model=sbml.getModel()
    for sp in model.getListOfSpecies():
        notes=sp.getNotesString()
        regexmatch=re.search('KEGG:\s+(\w+)',notes)
        if regexmatch:
            cid=regexmatch.group(1)

            check(sp.setNotes(notify(cid)),'species notes')
        else:
            print('KEGG-less:',notes)
    libsbml.writeSBMLToFile(sbml,"test enrich.xml")

def sbml_name_fixer():
    __doc__='Whereas cobra does not care of names, SBML2Latex does...'
    sbml=libsbml.SBMLReader().readSBML("Gthg_2.4.xml")
    error_check(sbml)
    model=sbml.getModel()
    log=open('Log of names SBML Gthg2.4.txt','w',encoding='utf-8')
    for sp in model.getListOfSpecies():
        notes=sp.getNotesString()
        regexmatch=re.search('NAME:\s+(.*?)[<;]',notes)
        old=sp.getName()
        if regexmatch:
            name=regexmatch.group(1)
            sp.setName(name)
        else:
            name='UNMATCHED!'
        log.write(old+T+name+N)

def nominate(sp): #serious botch
    notes=sp.getNotesString()
    regexmatch=re.search('NAME:\s+(.*?)[<;]',notes)
    old=sp.getName()
    if regexmatch:
        name=regexmatch.group(1)
        return name
    else:
        raise Exception('UNMATCHED!')


def sbml_rdf_fixer():
    __doc__='Too hard to workout. Abbanded.'
    sbml=libsbml.SBMLReader().readSBML("Gthg_2.4.xml")
    error_check(sbml)
    model=sbml.getModel()
    for sp in model.getListOfSpecies():
        if not sp.getCVTerms():
            print(sp.addCVTerm(libsbml.CVTerm(libsbml.BIOLOGICAL_QUALIFIER), True))
        bag=sp.getCVTerm(0)
        notes=sp.getNotesString()
        regexmatch=re.search('CHEBI:\s+(.*?)[<;]',notes)
        if regexmatch:
            name=regexmatch.group(1)
            bag.addResource("urn:miriam:chebi:CHEBI%3A"+name)
        regexmatch=re.search('KEGG:\s+(.*?)[<;]',notes)
        if regexmatch:
            name=regexmatch.group(1)
            bag.addResource("urn:miriam:kegg.compound:"+name)
    libsbml.writeSBMLToFile(sbml,"test enrich.xml")

def sbml_body():
    __doc__='SBML2Latex does not like the notes'
    sbml=libsbml.SBMLReader().readSBML("Gthg_2.4.xml")
    error_check(sbml)
    model=sbml.getModel()
    for sp in model.getListOfSpecies():
        notes='<body xmlns="http://www.w3.org/1999/xhtml">\n'+sp.getNotesString()+'\n</body>'
        sp.setNotes(notes.replace('<notes>\n','').replace('</notes>\n','').replace('html:',''))
    for rxn in model.getListOfReactions():
        notes='<body xmlns="http://www.w3.org/1999/xhtml">\n'+rxn.getNotesString()+'\n</body>'
        rxn.setNotes(notes.replace('<notes>\n','').replace('</notes>\n','').replace('html:',''))
    libsbml.writeSBMLToFile(sbml,"Gthg_2.4_html_tweak.xml")

def pathway_changes():
    pass

def get_vacantSpeciesIDs(model): #I could use a regex. but it is nice to be safe.
    taken={int(sp.getId().replace('M_','').replace('_c','').replace('_e','')) for sp in model.getListOfSpecies() if sp.getId() != 'M_Biomass_c' and sp.getId() != 'M_Biomass_e'}
    return [i for i in range(200,max(taken)) if i not in taken] #starting from 200... no reason.

def get_vacantReactionIDs(model):
    taken={int(sp.getId().replace('R_','').replace('M_','').replace('_c','').replace('_e','').replace('_out','').replace('biomass','0').replace('Biomass','0')) for sp in model.getListOfReactions()}
    return [i for i in range(200,max(taken)) if i not in taken] #starting from 200... no reason.


def add_KEGGspecies(model, cid): #C01516 morphine. test.
    #species id="M_24_c" name="Oxoglutaric acid" compartment="Cytosol" initialConcentration="1" boundaryCondition="false"
    mid='M_'+str(get_vacantSpeciesIDs(model)[0])+'_c'
    s1 = model.createSpecies()
    check(s1,                                 'create species s1')
    check(s1.setId(mid),       'set species s1 id')
    check(s1.setCompartment('Cytosol'),            'set species s1 compartment')
    check(s1.setConstant(False),              'set "constant" attribute on s1')
    check(s1.setInitialAmount(1),             'set initial amount for s1')
    #check(s1.setSubstanceUnits('mole'),       'set substance units for s1')
    check(s1.setBoundaryCondition(True),     'set "boundaryCondition" on s1')
    #check(s1.setHasOnlySubstanceUnits(False), 'set "hasOnlySubstanceUnits" on s1')
    check(s1.setNotes(notify(cid)),'adding notes')
    check(s1.setName(nominate(s1)),'name')
    return mid

def makespecies(model,cid,kegg={}):
    #getSpeciesbyKEGG
    if kegg =={}:
        for sp in model.getListOfSpecies():
            regexmatch=re.search('KEGG:\s+(.*?)[<;]',sp.getNotesString())
            if regexmatch:
                kegg[regexmatch.group(1)]=sp
    if cid in kegg:
        print(cid+' is '+kegg[cid].getId())
        return kegg[cid].getId()
    else:
        new=add_KEGGspecies(model, cid)
        kegg[cid]=model.getElementBySId(new)
        assert kegg[cid], 'Impossible!'
        print(cid+' added '+kegg[cid].getId()+' ('+new+')')
        return new


def parse_add(ceq,species_ref1):
    stex=re.match('\W*(\d+)',ceq)
    if stex:
        stoi=int(stex.group(1))
    else:
        stoi=1
    cex=re.search('([CGD]\d*)',ceq)
    if cex:
        cid=cex.group(1)
    else:
        raise Exception(ceq)
    print('working on: ',cid,stoi)
    sid=makespecies(model,cid)
    assert sid,'no sid'
    check(species_ref1.setSpecies(sid),      'assign reactant species')
    check(species_ref1.setStoichiometry(stoi),     'set stoi')

def add_KEGGreaction(model,kid, subsys, locus, reversable=False, backwards=False, override=None):  #R03591 morphine oxidoreductase. 2 new, 2 old
    rid='R_'+str(get_vacantReactionIDs(model)[0])
    rxn=model.createReaction()
    check(rxn,                                 'create reaction')
    check(rxn.setId(rid),                     'set reaction id')
    check(rxn.setReversible(reversable),            'set reaction reversibility flag')
    check(rxn.setName(kid),            'set name')
    check(rxn.setMetaId('metaid_'+rid),'set metaid')
    #check(rxn.setFast(False),                  'set reaction "fast" attribute')
    details=keggmaster(kid)
    notes='            <body xmlns="http://www.w3.org/1999/xhtml">\n'
    notes+='<p>SUBSYSTEM: '+subsys+' </p>\n'
    notes+='<p>GENE_ASSOCIATION: '+locus+'</p>\n'
    for f in sorted(details):
        if f =='ENTRY':
            pass
        if f =='EQUATION':
            notes+='<p>'+f+': '+'; '.join(details[f]).replace("<=>","\u21cc").replace('=>','\u2192').replace('<=','\u2190')+'</p>\n'
            equa=list(details[f])[0]
            if override:
                for o in override:
                    print('change')
                    print(equa)
                    equa=equa.replace(o,override[o])
                    print(equa)
            notes+='<p>EQUATION_USED_HEREIN: '+'; '.join(details[f]).replace("<=>","\u21cc").replace('=>','\u2192').replace('<=','\u2190')+'</p>\n'
        else:
            notes+='<p>'+f+': '+'; '.join(details[f]).replace("<=>","\u21cc").replace('=>','\u2192').replace('<=','\u2190')+'</p>\n'
    notes+='            </body>'
    print(notes)
    check(rxn.setNotes(notes),'adding notes...')

    if 'EQUATION' in details:
        equa=list(details['EQUATION'])[0]
        r,p=equa.replace('<','').split('=>') #I'll decide manually if reversible.
        print('parsing ',details['EQUATION'])
        assert p, 'no product'
        if backwards:
            r,p=p,r


    for sp in r.split('+'):
        species_ref1 = rxn.createReactant()
        check(species_ref1,                       'create reactant')
        parse_add(sp,species_ref1)
    for sp in p.split('+'):
        species_ref1 = rxn.createProduct()
        check(species_ref1,                       'create prod')
        parse_add(sp,species_ref1)

    kinetic_law = rxn.createKineticLaw()
    check(kinetic_law,                        'create kinetic law')
    print('remember to search and replace <kineticLaw/> with <kineticLaw><math xmlns="http://www.w3.org/1998/Math/MathML"> <ci> FLUX_VALUE </ci> </math><listOfParameters> <parameter id="LOWER_BOUND" value="-1000" units="mmol_per_gDW_per_hr"/><parameter id="UPPER_BOUND" value="1000" units="mmol_per_gDW_per_hr"/> <parameter id="OBJECTIVE_COEFFICIENT" value="0" units="mmol_per_gDW_per_hr"/><parameter id="FLUX_VALUE" value="0" units="mmol_per_gDW_per_hr"/></listOfParameters></kineticLaw>')


if __name__ == "__main__":
    #sbml_addkeggchem()
    #sbml_name_fixer()
    #sbml_rdf_fixer()
    sbml=libsbml.SBMLReader().readSBML("Gthg_2.5.xml")
    error_check(sbml)
    model=sbml.getModel()
    generow(model)
    #SBML_looter(model)
    #uncommitted...
    #add_KEGGreaction(model,'R03026', 'Butanoate metabolism', '( Gthg02239 or Gthg01064 or Gthg02240 or Gthg03320 )', reversable=True)




    #error_check(sbml)
    #libsbml.writeSBMLToFile(sbml,"Gthg_2.6.xml")


'''
    add_KEGGreaction(model,'R00717', 'Glycoxylate metabolism', 'Gthg03271', reversable=True)
    add_KEGGreaction(model,'R02946', 'Butanoate metabolism', 'Gthg04107', reversable=True)
    add_KEGGreaction(model,'R02855', 'Butanoate metabolism', 'Gthg04106', reversable=False)
    add_KEGGreaction(model,'R09524', 'Butanoate metabolism', '( Gthg01052 or Gthg01053 or Gthg01054 )', reversable=False, override={'C00466':'C00810'})
    add_KEGGreaction(model,'R00226', 'Butanoate metabolism', '( Gthg02903 or Gthg02904 )', reversable=False, backwards=True)
    add_KEGGreaction(model,'R10506', 'Butanoate metabolism', '( spontaneous )', reversable=False, override={'C00028':'C00016','C00030':'C01352'})
    add_KEGGreaction(model,'R01394', 'Glycoxylate metabolism', 'Gthg02058', reversable=True)
    add_KEGGreaction(model,'R01745', 'Glycoxylate metabolism', '( Gthg00526 or Gthg02057 )', reversable=True)
    add_KEGGreaction(model,'R01171', 'Butanoate metabolism', 'Gthg00577', reversable=True)

    '''