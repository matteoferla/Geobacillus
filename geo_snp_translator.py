__author__ = 'Matteo'
__doc__='''
R3 is different that R4. Changed:
=IFERROR(LEFT(B2,FIND("^",B2)-1),B2)
=IFERROR(RIGHT(B2,LEN(B2)-FIND("^",B2))-LEFT(B2,FIND("^",B2)-1),1)
header standardised.
'''
from Bio import SeqIO
import csv
#import XlsxWriter

def translate(listfp,genfp,outfp):
    for id in open(listfp):
        #RTMO00107-1358A-
        #I forgot that hypens were a deletions. Oh well.
        #Damn. Need to start from scratch
        pass

def get_common(fp1,fp2,genfp,outfp):
    def parse(fp):
        f=csv.reader(open(fp))
        header=next(f)
        data=[{header[i]: line[i] for i,v in enumerate(line)} for line in f]
        id=["|".join([row['Chromosome'],row['Region'],row['length'],row['Reference'],row['Allele']]) for row in data]
        return (data,id)

    def mutate(originalis,varietas,positio,longinquitas,ex,ad): #Modded Copypasted method.
    # I dont remember why it is in latin, bar for the mutanda joke
        verbositas=0
        positio=int(positio)
        longinquitas=int(longinquitas)
        mutanda=originalis.tomutable()
        if varietas == "SNV" or varietas == "SNP":
            if verbositas: print("SNP")
            if mutanda[positio-1] !=ex:
                print("SNP mismatch",ex,mutanda[positio-1])
            else:
                mutanda[positio-1]=ad
        elif varietas == "Deletion":
            if verbositas: print("Del")
            if ad !="-":
                print("Not a deletion:",ad)
            else:
                for n in range(longinquitas):
                    mutanda.pop(positio-1)
        elif varietas == "Insertion":
            if verbositas: print("Ins")
            if ex !="-":
                print("Not a insertion:",ex)
            else:
                for n in range(longinquitas):
                    mutanda.insert(positio-1+n,ad[n])
        else:
            print("ERROR")
        return mutanda.toseq()

    data1,id1=parse(fp1)
    data2,id2=parse(fp2)
    genome=list(SeqIO.parse(open(genfp, "rU"), "genbank"))
    idc=set(id1).intersection(set(id2))
    pairball=[[idc_j, data1[id1.index(idc_j)],data2[id2.index(idc_j)]] for idc_j in idc]
    comball=[]
    for pair in pairball:
        #Chromosome	Region	Type	length	Reference	Allele	Reference allele	Zygosity	Count	Coverage	Frequency	Forward/reverse balance	Average quality
        mutation={"name":pair[0],
                        "gene":pair[1]['Chromosome'],
                        "position":int(pair[1]['Region']),
                        "type":pair[1]['Type'],
                        "length":int(pair[1]['length']),
                        "from":pair[1]['Reference'],
                        "to":pair[1]['Allele'],
                        "count":[pair[1]['Count'],pair[2]['Count']],
                        "coverage":[pair[1]['Coverage'],pair[2]['Coverage']],
                        "frequency":[pair[1]['Frequency'],pair[2]['Frequency']],
                        "feq_avg":sum([float(pair[1]['Frequency']),float(pair[2]['Frequency'])])/2,
                        "count_avg":sum([float(pair[1]['Count']),float(pair[2]['Count'])])/2,
                        "cov_avg":sum([float(pair[1]['Coverage']),float(pair[2]['Coverage'])])/2,
                        "quality":[pair[1]['Average quality'],pair[2]['Average quality']]}
        for chr in genome:
            for gene in chr.features:
                if gene.type =="CDS":
                    #print(gene.qualifiers["locus_tag"],mutation["gene"])
                    #I forgot about the old tag
                    if gene.qualifiers["old_locus_tag"][0] == mutation["gene"]:
                        mutation['locus_tag']=gene.qualifiers["locus_tag"][0]
                        mutation['feature']=str(gene).replace("\n","   ").replace("\r","")
                        mutation['ori_seq']=gene.location.extract(chr).seq
                        mutation['mut_seq']=mutate(mutation['ori_seq'],mutation["type"],mutation["position"],mutation['length'],mutation['from'],mutation['to'])
                        mutation['ori_trans']=mutation['ori_seq'].translate()
                        mutation['mut_trans']=mutation['mut_seq'].translate()
                        mutation['Pfam']=""
                        if "db_xref" in gene.qualifiers:
                            for i in range(len(gene.qualifiers["db_xref"])):
                                if gene.qualifiers["db_xref"][i].find("Pfam:")>-1:
                                   mutation['Pfam']=gene.qualifiers["db_xref"][i].replace("Pfam:","")
                                elif gene.qualifiers["db_xref"][i].find("CL:")>-1:
                                   mutation['CL']=gene.qualifiers["db_xref"][i].replace("CL:","")
                        if "EC_number" in gene.qualifiers:
                            mutation['EC']=gene.qualifiers["EC_number"][0]
                        else:
                            mutation['EC']=""
                        #print(gene.qualifiers.keys())
                        if "function" in gene.qualifiers and gene.qualifiers["function"][0]:
                            #print(gene.qualifiers["function"][0])
                            mutation["function"]=gene.qualifiers["function"][0]
                        elif gene.qualifiers['product'][0]=="Hypothetical protein" or gene.qualifiers['product'][0]=="hypothetical proteins":
                            mutation["function"]="unknown"
                        else:
                            mutation["function"]=""
                        if "product" in gene.qualifiers:
                            mutation["product"]=gene.qualifiers["product"][0]
                        else:
                            mutation["product"]="unknown"
                        if "MUTATED" in gene.qualifiers:
                            gene.qualifiers["MUTATED"]+=1
                        else:
                            gene.qualifiers["MUTATED"]=1




                        #from difflib import Differ
                        #d=list(Differ().compare(str(mutation['ori_trans']),str(mutation['mut_trans'])))
                        #start=[(i,l[-1]) for i,l in enumerate(d) if l[0]=="-"]
                        #stop=[(i,l[-1]) for i,l in enumerate(d) if l[0]=="+"]
                        #print(start,stop)
                        import math
                        n=0
                        p=math.floor((mutation["position"]-1)/3)
                        timebomb=0
                        mutation['AA_diff']=[p,'','']
                        mutation['AA_result']="missense" #default
                        while 1==1:
                            if len(mutation['mut_trans'])-1 <n+p:
                                break
                            a=mutation['ori_trans'][p+n]
                            b=mutation['mut_trans'][p+n]
                            mutation['AA_diff'][1]+=str(a)
                            mutation['AA_diff'][2]+=str(b)
                            if b=="*":
                                mutation['AA_result']="nonsense"
                                break
                            if a == b:
                                if timebomb==0:
                                    timebomb=1
                                else:
                                    if n==1:
                                        mutation['AA_result']="silent"
                                    else:
                                        mutation['AA_result']="missense"
                                    break
                            n+=1
        comball.append(mutation)

    for chr in genome:
            for gene in chr.features:
                if gene.type =="CDS":
                    if 'MUTATED' not in gene.qualifiers:
                        gene.qualifiers['MUTATED']=0
                    if 'function' not in gene.qualifiers:
                        gene.qualifiers['function']=[""]
                    print(str(gene.qualifiers['locus_tag'][0])+"\t"+str(gene.qualifiers['MUTATED'])+"\t"+str(gene.qualifiers['function'][0])+"\t"+str(gene.qualifiers['product'][0])+"\t"+"\n")
    k=list(comball[1].keys())
    done=csv.DictWriter(open(outfp,"w"),k)
    done.writeheader()
    done.writerows(comball)
    return comball



if __name__=="__main__":
    databall=get_common("R2A4.csv","R3A3.csv","Gthg_from_embl_pfamed.gb","out2.csv")



    #translate("common.txt","Gth_NCIMB_11955_v2.4","out.csv")