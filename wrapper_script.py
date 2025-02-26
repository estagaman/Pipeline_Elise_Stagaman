#input files needed:
    #metadata file with columns Sample, Donor, and Condition - can have additional columns, but not required
    #a directory of paired-end fastq files 
    #the sample ID used for the fastq files must match sample ID in the metadata

#argparse to take in input files
import sys
import argparse
#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="wrapper_script.py") #specify the name of the script
    parser.add_argument("-m", "--metadata", #add metadata file argument
    help="metadata csv file",
    required=True)
    parser.add_argument("-d", "--directory", #add directory with test data argument
    help="fastq file directory",
    required=True)
    return parser.parse_args(args)

#retrieve command line arguments 
arguments = check_arg(sys.argv[1:]) #get the arguments that I inputted from the terminal
metadata_csv = arguments.metadata #load metadata
data_dir = arguments.directory #load fastq directory

import os #for running command line tools from script
import pandas as pd #for data frames
from Bio import SeqIO #for parsing fasta files and interacting with NCBI
import glob #for finding fastq files in our fastq directory

##STEP 1: DIRECTORY SETUP
os.system("mkdir PipelineProject_Elise_Stagaman") #make new directory for my output

info_by_sample = pd.read_csv(metadata_csv) #load metadata into pandas data frame

samples = info_by_sample["Sample"].to_list() #use "Sample" column as sample names

fastq_names = glob.glob(data_dir + "/*.fastq") #find the path to each individual fastq file

for name in fastq_names: #for each fastq file
    os.system(f"cp {name} PipelineProject_Elise_Stagaman") #copy it into our output directory so I can refer to it easily later

os.chdir("PipelineProject_Elise_Stagaman") #move into my output directory

##STEP 2: build a transcriptome index for HCMV NC_006273.2
os.system("efetch -db nuccore -id NC_006273.2 -format genbank > HCMV.gb") #get the GenBank record for HCMV

cds_count = 0
with open("HCMV_cds.fasta", 'w') as out_fasta: #open up a fasta file to write the CDS records to
    for record in SeqIO.parse("HCMV.gb", "genbank"): #parse the GenBank record
        for feature in record.features: #check each feature
            if feature.type == "CDS": #find CDS
                cds_count += 1 #add 1 to our count of CDS sequences

                #write the protein ID and sequence for that CDS to the fasta file
                out_fasta.write(">" + feature.qualifiers["protein_id"][0] + "\n" + str(feature.location.extract(record).seq) + "\n")

#add the number of CDS sequences to the log file
with open("PipelineProject.log", "w") as f:
    f.write(f'The HCMV genome (NC_006273.2) has {cds_count} CDS.\n')

##STEP 3: KALLISTO
os.system("kallisto index -i HCMV_index.idx HCMV_cds.fasta") #index the HCMV fasta

os.system("mkdir results") #make a directory to hold the kallisto esults 

info_by_sample = info_by_sample.set_index("Sample") #make "Sample" the index of our metadata for easy access later

#do the quantification
for samplename in samples: #for each sample
    forward = samplename + "_1.fastq" #find the forward and reverse filenames
    reverse = samplename + "_2.fastq"
    output = "results/" + samplename #make a results folder for each sample
    os.system(f"time kallisto quant -i HCMV_index.idx -o {output} -b 10 -t 2 {forward} {reverse}") #quantify using kallisto

#start list to hold tpm statistics from kallisto
tpm_stats = [["sample", "condition", "min_tpm", "med_tpm", "mean_tpm", "max_tpm"]]

#find statistics for each sample and add to our list
for samplename in samples: #for each sample
    abund_file = "results/" + samplename + "/abundance.tsv" #access the abundance file
    abundances = pd.read_csv(abund_file, sep="\t") #read it in as a pandas data frame
    tpm = abundances["tpm"] #extract the tpms

    #save the tpm information I want as a list, append it to tpm_stats
    tpm_stats.append([samplename, str(info_by_sample.loc[samplename, "Condition"]), str(tpm.min()), str(tpm.median()), str(tpm.mean()), str(tpm.max())])

#output tpm table to log file 
with open("PipelineProject.log", "a") as file:
    for row in tpm_stats: #add each list of the tpm stats as a line
        file.write("\t".join(row) + "\n")

#STEP 4: SLEUTH

#make txt to input into sleuth
sleuth_data = info_by_sample.reset_index() #make "Sample" its own column again
sleuth_data = sleuth_data[["Sample", "Condition"]] #take just the Sample and Condition columns

#create the results paths for each sample
results_paths = [] 
for sample in samples:
    results_path = "results/" + sample 
    results_paths.append(results_path)

sleuth_data["Path"] = results_paths #add results paths to our dataframe

sleuth_data.columns = ["sample", "condition", "path"] #rename columns so sleuth can recognize them
sleuth_data.to_csv("sleuth_input.txt", sep="\t", index=False) #save as a txt to our current directory

os.system("Rscript ../sleuth_script.R") #run the Sleuth R script

##STEP 5: BOWTIE2, create genome index, map strains to patient samples

#use NCBI datasets to retrieve HCMV genome and unzip, build the index as HCMV Bowtie
os.system("datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report")
os.system("unzip ncbi_dataset.zip")
os.system("bowtie2-build ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna HCMVBowtie")

for sample in samples: #for each sample, find forward and reverse read paths
    forward = sample + "_1.fastq"
    reverse = sample + "_2.fastq"
    mapped_sample = sample + "_mapped_%.fq.gz" #name the mapped sample

    #map reads to the HCMV genome, keeping only the mapped reads in a new fastq file
    os.system(f"bowtie2 --quiet -x HCMVBowtie -1 {forward} -2 {reverse} -S HCMVmap.sam --al-conc-gz {mapped_sample}")

#counting the number of mapped reads:
os.system("mkdir bowtie_counts") #make a new directory to hold count files in

for sample in samples: #for each sample
    forward = sample + "_1.fastq" #find forward reads before mapping
    for_mapped = sample + "_mapped_2.fq.gz" #find forward reads after mapping

    outfile = "bowtie_counts/BT_out_" + sample #format output file for counts
    for file in [forward, for_mapped]: #for before and after mapping
        os.system(f"wc -l < {file} >> {outfile}") #count the number of lines in a file and save it to the output location

    with open(outfile, "r") as f: #open each output file again
        counts = [int(line.strip()) for line in f] #make a list for with the counts "before" and counts "after"
        before = str(round(counts[0]/4)) #divide counts by 4 because every 4th line in a fastq is the sequence
        after = str(round(counts[1]/4))
        with open("PipelineProject.log", "a") as file:
            #save output string with counts to our log file
            out_str = "Donor " + str(info_by_sample.loc[sample, "Donor"]) + " (" + str(info_by_sample.loc[sample, "Condition"]) + ") had " + before + " read pairs before Bowtie2 filtering and " + after + " read pairs after. \n"
            file.write(out_str)

##STEP 6: SPAdes - one assembly for each patient 
for donor in set(info_by_sample["Donor"].to_list()): #for each donor
    one_donor = info_by_sample.index[info_by_sample["Donor"] == donor].tolist() #find the samples that match that donor
    for_1 = one_donor[0] + "_mapped_1.fq.gz" #forward for 1st donor
    rev_1 = one_donor[0] + "_mapped_2.fq.gz" #reverse for 1st donor
    for_2 = one_donor[1] + "_mapped_1.fq.gz" #forward for 2nd donor
    rev_2 = one_donor[1] + "_mapped_2.fq.gz" #reverse for 2nd donor
    output_location = "Donor" + str(donor) + "assembly/" #output file for the assembly

    #make the spades command for all 4 sequences
    spades_command = f"spades.py -k 77 -t 2 --only-assembler --pe-1 1 {for_1} --pe-2 1 {rev_1} --pe-1 2 {for_2} --pe-2 2 {rev_2} -o {output_location}"
    os.system(spades_command) #run the spades command
    with open("PipelineProject.log", "a") as file:
        file.write(spades_command + "\n") #save the command to log file

#STEP 7: BLAST
os.system("datasets download virus genome taxon Betaherpesvirinae --refseq --include genome") #extract Betaherpesvirinae RefSeq entries using NCBI datasets
os.system("unzip ncbi_dataset.zip") #unzip the datasets retrieval

file_name = "ncbi_dataset/data/genomic.fna" #find the fasta file
db_name = "betaherpesvirinae" #name the database
makeblast_command='makeblastdb -in '+file_name+' -out '+db_name+' -title '+db_name+' -dbtype nucl' #build blast command

os.system(makeblast_command) #make a blast database

output_list = []

donors = list(set(info_by_sample["Donor"].to_list())) #find unique donors in our data
for i in range(len(donors)):
    donors[i] = "Donor" + str(donors[i]) #make a list ["Donor1", "Donor3"]

for assembly in donors: #for each donor's assembly
    contig_path = assembly + "assembly/contigs.fasta" #construct a path to the contigs

    #find the longest contig - make into it's own fasta and save
    with open(contig_path, 'r') as f:
        longest = next(SeqIO.parse(f, "fasta"))

        out_fasta = assembly + "_longest.fasta"

        with open(out_fasta, "w") as out_fas:
            SeqIO.write(longest, out_fas, "fasta")
        
    query_seqfile = out_fasta #call that fasta again

    output_file = "BLAST_" + assembly + ".csv" #save BLAST output as csv

    output_list.append(output_file) #hold the name of that csv in a list for later

    #run blast on the longest contig, saving results to specified output file
    blast_command = "blastn -query " + query_seqfile + " -db betaherpesvirinae -out " + output_file + " -max_hsps 1 -outfmt '10 sacc pident length qstart qend sstart send bitscore evalue stitle'"
    os.system(blast_command)

donors = list(set(info_by_sample["Donor"].to_list())) #find unique donors again

for i in range(2): #for each donor
    blast_result = output_list[i] #look at the BLAST output, load it in as a pandas dataframe
    results = pd.read_csv(blast_result, header = None, usecols = range(11))

    #name columns of this dataframe
    results.columns = ['sacc', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'evalue', 'stitle', 'genome']

    #add "complete genome" to the stitle instead of having it separate
    results['stitle'] = results['stitle'] + ", " + results['genome']
    results = results.drop("genome", axis=1)

    if len(results) > 10: #if there are more than 10 sequences in the results, keep only the top 10
        results = results.head(10)

    #add BLAST results for each donor to the log file
    with open("PipelineProject.log", "a") as file:
        file.write("Donor " + str(donors[i]) + ": \n")
        results.to_csv(file, sep='\t', index=False, mode='a')
