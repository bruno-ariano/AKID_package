#########################################################
# find_determinant				                        #
# University of Rome Torvergata                         #
# Author: Bruno Ariano (bruno.ariano.87@gmail.com)      #
#########################################################

import os,sys,time,sys
from functions import mapping_position, determinants,find_domain,select_peptide_window
import argparse
from Bio import AlignIO, SeqIO
directory=os.getcwd()
usage = """%(prog)s reads fasta sequence file and return the determinants of the kinase domains)."""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-k", "--k", "-kinase", dest="kinase",
                  help="name of the kinase file in the 'INPUT' directory")
p.add_argument("-p", "--p", "-peptide", "--peptide", dest="peptide",
                  help="name of the peptide file in the 'INPUT' directory")
p.add_argument("-o", "--o", "-output", "--output", dest="output",
                  help="name of the output files")
p.add_argument("-t", "--t", "-type", "--type", dest="type",
                  help="type of prediction, single (default) or discovery")

args = p.parse_args()
if len(sys.argv)==1: print "No valid arguments in input\n\n";p.print_help();sys.exit(1)

arg1 =os.getcwd()+"/input/"+args.kinase
arg2 =os.getcwd()+"/input/"+args.peptide
arg3 =args.output

# check if the input files exist
try: open(arg1)
except IOError: print ("No kinase file with name \"" + args.kinase + "\"");sys.exit(1)

try: open(arg2)
except IOError: print ("No peptide file with name \"" + args.peptide + "\"");sys.exit(1)

# creating a directory to store the temporary files
dir_path = os.getcwd()


os.system("mkdir " + dir_path + "/tmp")
os.system("mkdir " + dir_path + "/tmp/kinase_sequences")
os.system("mkdir " + dir_path + "/tmp/alignments")


#this function check if kinase domains were found
dict_domain=find_domain(dir_path,arg1)
if dict_domain==1: print "No kinase domains found in input sequences"; exit(1)

#this part is used to map the determinants in the kinase domains
domains=open(dir_path+"/db/domain_sequence","r").read()
deter=determinants( dir_path + "/db/determinant")
deter_2=open(dir_path + "/kinase_domain_determinants/"+ arg3 +"_determinants.fasta","w")
domain={}
for seq_record in SeqIO.parse(dir_path + "/tmp/domain_found.fasta", "fasta"):
	if "|" in seq_record.id:
		seq_record.id=seq_record.id.split("")[1]
	if seq_record.id in dict_domain and "X" not in seq_record.seq:
		domain[seq_record.id]=">" + domains + ">" + seq_record.id + "\n" + seq_record.seq
		sequences_output=open(dir_path + "/tmp/kinase_sequences/" + seq_record.id,"w")
		sequences_output.write(str(domain[seq_record.id])[1:-1])
		sequences_output.close()
		os.system(dir_path + "/PACKAGES/hmmer-3.1b2-linux-intel-x86_64/binaries/hmmalign -o " + dir_path + "/tmp/alignments/align_"+seq_record.id+" "+ dir_path+"/db/hmmdb/CompleteKinomeAlignmentManuallyRefined_v3.hmm" + " " + dir_path + "/tmp/kinase_sequences/"+seq_record.id)
		print "mapping " + seq_record.id
		diz1=mapping_position(">" + seq_record.id, dir_path + "/tmp/alignments/align_" + seq_record.id,  dir_path + "/tmp/dictionary/diz_"+seq_record.id)
		a1=AlignIO.read(open(dir_path + "/tmp/alignments/align_" + seq_record.id,"r"),"stockholm")
        for record in a1:
        	if record.id==seq_record.id:
				k_deter=[]
				for j in deter:
					k_deter+=record.seq[diz1.keys()[diz1.values().index(int(j)-1)]]
				deter_2.write(">"+seq_record.id+"\n"+"".join(k_deter)+"\n")


#Here we divide the substrate in 15mers sequences

peptides_in=open(arg2).read().split(">")
peptides_out=open(dir_path+"/peptide/"+ arg3 +"_peptides.fasta","w")
del(peptides_in[0])
for record in peptides_in:
	record=record.split("\n")
	amino_seq="".join(record[1:])
	if args.type=='discovery':
		for pos,amino in enumerate(amino_seq):
			if amino in ["S","T","Y"]:
				pep_seq=select_peptide_window(pos,amino_seq)
				peptides_out.write(">"+record[0]+"_"+amino+str(pos+1)+"\n"+pep_seq+"\n")
	elif args.type=='single':
		peptides_out.write(">"+record[0]+"\n"+amino_seq+"\n")
os.system("rm -r " + dir_path + "/tmp")
sys.exit()
############
