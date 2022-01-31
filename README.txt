# Prediction of kinase peptide interaction in Eukaryotic organisms using aminoacid sequences information (Parca L.*, Ariano B*, et al. 2019)

###########################
#                         #
# R E Q U I R E M E N T S #
#                         #
###########################


# This tool require: an R-cran version >=3.1, Biopython >=1.68, hmmer>=3.1b2, python>=3.0

# To install the network:

unzip AKID_package.zip
cd AKID_package
Rscript install_AKID.R

#############################
#                           #
# P R E P A R E   I N P U T #
#                           #
#############################

# To find the determinant sequences and generate the peptides, type from the AKID directory the following command.
# The determinant sequences and the peptides will be stored respectively in kinase_domain_determinants and peptide directory.

python create_DoS-peptides.py -k example_kinases.fasta -p example_peptides.fasta -o example -t single

# Note that both the kinase sequences and the peptides must be stored in the input/ directory.
# If you already have a list of on specific 15-residues-long peptides (centered on S/T/Y residues) run the script with the "-t single" option.
# In this case you can put the peptides in a single fasta-formatted file in the peptide/ directory (see the example file).
# PLEASE USE UNIQUE IDENTIFIER (E.G. UNIPROT ID PLUS THE RESIDUE NUMBER AT THE CENTER OF THE QUERY PEPTIDE) IN THE FASTA HEADER OF EACH QUERY PEPTIDE.

# If you want to analyze every S/T/Y residues in whole protein sequences run the script with the "-t discovery" option.
# In this case you can put the protein sequences in a single fasta-formatted file in the input/ directory. The script will prepare the query 15-residues-long peptides in the peptide/ directory.
# You can see the determinants calculated for all the kinases in a fasta-formatted file in the kinase_domain_determinants/ directory.


#######################
#                     #
# P R E D I C T I O N #
#                     #
#######################

# In order to make the predictions run the following command.
# The determinant sequences must be in the kinase_domain_determinants/ direectory, the query peptides must be in the peptide/ directory.
# The prediction results will be reported in the prediction directory.
# PLEASE note that for large input AKID will write rather large files (approximately 100MB for the prediction of 1 Kinase against 20000 query peptides). These files are written in a temporary folder 
# (tmp/ folder in the main AKID directory) that will be deleted at the end of the prediction. Make sure you have enough space for the prediction. In case you want to redirect the temporary folder 
# somewhere else you can edit the make_prediction.py script, changing the "tmp_dest" variable (31st line of the code) to the desired directory: this will be the directory where the tmp/ folder will be created and deleted.

python make_prediction.py -dos example_determinants.fasta -p example_peptides.fasta -o example
