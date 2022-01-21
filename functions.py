#########################################################
# functions                                             #
# University of Rome Torvergata                         #
# Author: Bruno Ariano (bruno.ariano.87@gmail.com)      #
#########################################################
import random
import os
from Bio import SeqIO,AlignIO
from Bio.Align import MultipleSeqAlignment
import subprocess
import re
import sys
sys.path.append(os.getcwd()+"/PACKAGES/biopython-1.68")
dir_path = os.getcwd()
def determinants(f1):
    determinanti=open(f1,"r").read().split("\n")
    deter=[]
    for i in range(len(determinanti)-1):
        determinanti[i]=determinanti[i].split("\t")
    del determinanti[-1]
    del determinanti[0]
    del determinanti[0]
    for i in range(len(determinanti)):
        deter.append(int(determinanti[i][0]))
    deter=set(deter)
    deter=sorted(deter)
    return deter


def mapping_position(arg1,arg2,arg3):
    gapsChar = '-'
    class KinaseMSA:
        def __init__(self, query_kinases, QueryMSA, queryformat):
            tempalign = AlignIO.read(open(QueryMSA, 'rU'), queryformat)
            tempalign.sort()
            RecordsList = []
            query_kin_records = []
            for i in range(len(query_kinases)):
                query_kinases[i]=query_kinases[i][1:]
            self.nquerykin = len(query_kinases)
            for k in tempalign:
                if not k.id in query_kinases:
                    RecordsList.append(k)
                else:
                    query_kin_records.append(k)
            self.queryMSA = MultipleSeqAlignment(RecordsList)
            self.queryMSA.sort()
            self.queryMSA.extend(query_kin_records)
            self.queryMSAncols = self.queryMSA.get_alignment_length()
            self.queryMSAnrows = len(self.queryMSA)
            self.queryStart = self.queryMSAnrows-self.nquerykin
            self.TrainGapRE = re.compile('['+gapsChar+']'+'{'+str(self.queryStart)+'}')
        def map_new_alignment_to_trainMSA(self, TrainMSA, trainformat, TrainExcludedPos = []):
            self.trainMSA = AlignIO.read(open(TrainMSA, 'rU'), trainformat)
            self.trainMSA.sort()
            self.trainMSAncols = self.trainMSA.get_alignment_length()
            trainMSApos = 1
            self.queryMSAposToTrainMSApos = {}
            for queryMSApos in range(self.queryMSAncols):
                queryTrainCol = self.queryMSA[:self.queryStart,queryMSApos-1].upper()
                oldTrainCol = self.trainMSA[:,trainMSApos-1].upper()
                if queryTrainCol == oldTrainCol:
                    self.queryMSAposToTrainMSApos[queryMSApos] = trainMSApos
                    trainMSApos+=1
                    if trainMSApos == self.trainMSAncols:
                        break
            return self.queryMSAposToTrainMSApos
    mapping =KinaseMSA([arg1], arg2, 'stockholm')
    dictionary_mapping=mapping.map_new_alignment_to_trainMSA(dir_path+"/db/domain_sequence_align",'stockholm')
    return dictionary_mapping
    

def write_domain_from_table(fasta_sequence, table, destination):
    file_input = open(table).read().split("\n")
    del(file_input[-1])
    dictionary_table = {}
    dictionary_sequence = {}
    fasta=[]
    domain=[]
    for x in file_input:
        x = x.split("\t")
        fasta.append(x[0])
        ID_table = x[0]
        if "|" in ID_table:
            ID_table=ID_table.split("|")[1]
        
        if ID_table in domain:
            position = x[5:7]
            domain_count += 1
            dictionary_table[ID_table+"_domain"+str(domain_count)]=position
        else:
            position = x[5:7]
            domain_count = 1
            dictionary_table[ID_table+"_domain"+str(domain_count)]=position
        domain.append(ID_table)
    for seq_record in SeqIO.parse(fasta_sequence, "fasta"):
            if "|" in seq_record.id:
                ID_fasta=seq_record.id.split("|")[1]
            else:
                ID_fasta = seq_record.id
            dictionary_sequence[ID_fasta]=seq_record.seq

    dictionary_domain = {}
    for i in dictionary_table:
        if "|" in i:
            i=i.split("|")[1]+ "_" + i.split("_")[-1]
        dictionary_domain[i] = dictionary_sequence[i.split("_domain")[0]][int(dictionary_table[i][0])-1:int(dictionary_table[i][1])]
       
    h = open(destination, "w")
    for i in dictionary_domain:
        h.write(">" + i + "\n" + str(dictionary_domain[i])[1:-1] + "\n")
    return dictionary_domain.keys()

def find_domain(directory,arg1):
    dict_domain={}
    dir_path = directory
    os.system("hmmscan --domtblout " + dir_path + "/tmp/domain_kinase_t " + dir_path + "/db/hmmdb/CompleteKinomeAlignmentManuallyRefined_v3.hmm "+arg1 + " > output_hmm")
    os.system("rm " + dir_path +"/output_hmm")
    #subprocess.call(["hmmscan","--domtblout " + dir_path + "/output/tmp/domain_kinase_t " + dir_path + "/db/hmmdb/CompleteKinomeAlignmentManuallyRefined_v3.hmm "+arg1])
    os.system("python "+dir_path+"/script/parse_domain_table.py -f " + dir_path + "/tmp/domain_kinase_t")
    dict_domain=write_domain_from_table(arg1, dir_path + "/tmp/domain_kinase_t.out",  dir_path + "/tmp/domain_found.fasta")
    if len(dict_domain)==0:
        return 1
    else:
        return dict_domain


def code_amino(amino):
    cod={"A":0,"R":0,"N":0,"D":0,"C":0,"E":0,"Q":0,"G":0,"H":0,"I":0,"L":0,"K":0,"M":0,"F":0,"P":0,"S":0,"T":0,"W":0,"Y":0,"V":0,"-":0}
    if amino.upper() in cod.keys():
        cod[amino.upper()]+=1
        return cod.values()
    else:
        cod["-"]+=1
    return cod.values()


def code_fasta_file(fasta_file):
	code_seq={}
	f1=open(fasta_file,"r").read().split(">")
	del(f1[0])
	for i in f1:
		i=i.split("\n")
		seq="".join(i[1:])
		code_seq[i[0]]=[]
        	if i[0]!="":
        		for j in seq:
				code_seq[i[0]]+=code_amino(j)
	return code_seq


def select_peptide_window(position,amino_seq):
    aminos=set(["A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V","-"])
    peptide=""
    for res in range(position-7,position+8):
        if res<0 or res>=len(amino_seq) or amino_seq[res] not in aminos:
            peptide+="-"
        else:
            peptide+=amino_seq[res]
    return peptide
