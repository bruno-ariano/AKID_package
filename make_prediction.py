#########################################################
# make_prediction                       				#
# University of Rome Torvergata                         #
# Author: Bruno Ariano (bruno.ariano.87@gmail.com)      #
#########################################################

import os
import sys
import glob
from Bio import AlignIO
import argparse
from functions import code_fasta_file
dir_path = os.getcwd()
usage = """%(prog)s reads determinants and peptide files and returns the probability of interaction between kinases and peptides."""
p = argparse.ArgumentParser(description=usage)
p.add_argument("-dos", "--dos", dest="dos",
                   help="determinant of specificity residues file")
p.add_argument("-p", "--p", "--peptides", "-peptides", dest="peptides",
                   help="peptides file")
p.add_argument("-o", "--o", "--output", "-output", dest="output",
                   help="output prediction file")
p.add_argument("--temp", "-temp", dest="temporary",
                   help="temporary folder")

args = p.parse_args()
arg1 = dir_path + "/kinase_domain_determinants/" + args.dos
arg2 = dir_path + "/peptide/" + args.peptides
arg3 = args.output

tmp_dest=dir_path
if os.path.exists(tmp_dest+"/tmp")==False:os.system("mkdir " + tmp_dest+"/tmp")
if os.path.exists(tmp_dest+"/tmp/predictions/")==False:os.system("mkdir " + tmp_dest+"/tmp/predictions/")

peptide=code_fasta_file(arg2)
k_deter_code=code_fasta_file(arg1)
print "creating the files for the prediction"

okk={}
for i in k_deter_code:
	print 'Preparing kinase %d of %d\r'%(k_deter_code.keys().index(i)+1,len(k_deter_code)),;sys.stdout.flush()
	f1=open(tmp_dest+"/tmp/predictions/"+i,"w")
	f2=open(tmp_dest+"/tmp/predictions/"+i+"_interaction","w")
        for j in peptide:
		f1.write(str(k_deter_code[i])[1:-1]+","+str(peptide[j])[1:-1]+"\n")
		f2.write(i+"\t"+j+"\n")
	f1.close()
	f2.close()
	okk[i]=''

prediction_out=open(dir_path + "/prediction/"+ arg3 +"_predictions","w")
prediction_out.write("kinase_domain\tpeptide\tscore\n")
prediction_out.close()
print "\nmaking the predictions\n"

for i in okk:
	print "processing Kinase %d"%(okk.keys().index(i)+1) +" of %d (%s)"%(len(okk),i)
	os.system("Rscript "+dir_path+"/script/network.r " + tmp_dest+"/tmp/predictions/"+str(i)+ " " + tmp_dest+"/tmp/predictions/"+str(i)+"_interaction > "+ tmp_dest+"/tmp/output_predictions")
	os.system("cat " + tmp_dest+"/tmp/predictions/"+str(i)+ "_pred >> " + dir_path + "/prediction/"+ arg3 +"_predictions")

#prediction_count=len(glob.glob(tmp_dest+"/tmp/predictions/*"))/2
#m=0
#for i in range(1,t+1):
#	m+=1
#	print "processing "+ str(m) +" file of "+ str(prediction_count)
#	os.system("Rscript "+dir_path+"/script/network.r " + tmp_dest+"/tmp/predictions/"+str(i)+ " " + tmp_dest+"/tmp/predictions/"+str(i)+"_interaction > "+ tmp_dest+"/tmp/output_predictions")
#	os.system("cat " + tmp_dest+"/tmp/predictions/"+str(i)+ "_pred >> " + dir_path + "/prediction/"+ arg3 +"_predictions")
os.system("rm -r " + tmp_dest+"/tmp")
