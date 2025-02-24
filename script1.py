#input files needed:
    #metadata file with samplenames titled Sample, Condition column, and Donor column
    #a directory of paired-end fastq files - titled whatever you would like 
    #the sample ID used for the fastq files must match sample ID in the metadata

#argparse 
import sys
import argparse
#function to parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="wrapper_script.py") #specify the name of the script
    parser.add_argument("-m", "--metadata", #add input file argument
    help="metadata csv file",
    required=True)
    parser.add_argument("-d", "--directory", #add output file argument
    help="fastq file directory",
    required=True)
    return parser.parse_args(args)
#retrieve command line arguments 
arguments = check_arg(sys.argv[1:]) #get the arguments that I inputted from the terminal
metadata_csv = arguments.metadata 
data_dir = arguments.directory 

import os 
import pandas as pd
from Bio import SeqIO
import glob

os.system("mkdir PipelineProject_Elise_Stagaman")

info_by_sample = pd.read_csv(metadata_csv)

samples = info_by_sample["Sample"].to_list()

fastq_names = glob.glob(data_dir + "/*.fastq")

for name in fastq_names:
    copy_command = f"cp {name} PipelineProject_Elise_Stagaman"
    os.system(copy_command)

os.chdir("PipelineProject_Elise_Stagaman")

#STEP 2: build a transcriptome index for HCMV NC_006273.2
os.system("efetch -db nuccore -id NC_006273.2 -format genbank > HCMV.gb") #get the fasta format of HCMV

cds_count = 0
with open("HCMV_cds.fasta", 'w') as out_fasta:
    for record in SeqIO.parse("HCMV.gb", "genbank"):
        if record.features:
            for feature in record.features:
                if feature.type == "CDS":
                    cds_count += 1
                    out_fasta.write(">" + feature.qualifiers["protein_id"][0] + "\n" + str(feature.location.extract(record).seq) + "\n")

with open("PipelineProject.log", "w") as f:
    f.write(f'The HCMV genome (NC_006273.2) has {cds_count} CDS.\n')

#index with kallisto
os.system("kallisto index -i HCMV_index.idx HCMV_cds.fasta")

#STEP 3: quantify the TPM of each CDS in each transcriptome using kallisto
    #each sample and time period
    #min, med, mean, max TPM --> abundance.tsv
    #write calculations to log file

os.system("mkdir results") #make a directory to hold the results 

#create a dataframe of metadata by sample
info_by_sample = info_by_sample.set_index("Sample")

#do the quantification
for samplename in samples: 
    forward = samplename + "_1.fastq"
    reverse = samplename + "_2.fastq"
    output = "results/" + samplename
    os.system(f"time kallisto quant -i HCMV_index.idx -o {output} -b 10 -t 2 {forward} {reverse}")

#start table to hold tpm statistics 
tpm_stats = [["sample", "condition", "min_tpm", "med_tpm", "mean_tpm", "max_tpm"]]

#find statistics for each sample and add to our list
for samplename in samples:
    abund_file = "results/" + samplename + "/abundance.tsv"
    abundances = pd.read_csv(abund_file, sep="\t")
    tpm = abundances["tpm"]
    tpm_stats.append([samplename, str(info_by_sample.loc[samplename, "Condition"]), str(tpm.min()), str(tpm.median()), str(tpm.mean()), str(tpm.max())])

#output tpm table to log file 
with open("PipelineProject.log", "a") as file:
    for row in tpm_stats: 
        file.write("\t".join(row) + "\n")

#STEP 4: use kallisto output for sleuth
    #find DEGs between timepoints
    #write the results for each significant transcript to log file

#make txt to input into sleuth
sleuth_data = info_by_sample.reset_index()
sleuth_data = sleuth_data[["Sample", "Condition"]]
results_paths = []

for sample in samples:
    results_path = "results/" + sample
    results_paths.append(results_path)

sleuth_data["Path"] = results_paths
sleuth_data.columns = ["sample", "condition", "path"]
sleuth_data.to_csv("sleuth_input.txt", sep="\t", index=False)

os.system("Rscript ../sleuth_script.R")

#bowtie2: create genome index, map strains to patient samples

os.system("datasets download genome accession GCF_000845245.1 --include gff3,rna,cds,protein,genome,seq-report")
os.system("unzip ncbi_dataset.zip")
os.system("bowtie2-build ncbi_dataset/data/GCF_000845245.1/GCF_000845245.1_ViralProj14559_genomic.fna HCMVBowtie")

for sample in samples: 
    forward = sample + "_1.fastq"
    reverse = sample + "_2.fastq"
    mapped_sample = sample + "_mapped_%.fq.gz"
    os.system(f"bowtie2 --quiet -x HCMVBowtie -1 {forward} -2 {reverse} -S HCMVmap.sam --al-conc-gz {mapped_sample}")

#counting the number of mapped reads:
os.system("mkdir bowtie_counts")

for sample in samples:
    forward = sample + "_1.fastq"
    for_mapped = sample + "_mapped_2.fq.gz"

    outfile = "bowtie_counts/BT_out_" + sample
    for file in [forward, for_mapped]:
        os.system(f"wc -l < {file} >> {outfile}")

    with open(outfile, "r") as f:
        counts = [int(line.strip()) for line in f]
        before = str(round(counts[0]/4))
        after = str(round(counts[1]/4))
        with open("PipelineProject.log", "a") as file:
            write_str = "Donor " + str(info_by_sample.loc[sample, "Donor"]) + " (" + str(info_by_sample.loc[sample, "Condition"]) + ") had " + before + " read pairs before Bowtie2 filtering and " + after + " read pairs after. \n"
            file.write(write_str)

#SPAdes - one assembly for each patient 
for donor in set(info_by_sample["Donor"].to_list()):
    one_donor = info_by_sample.index[info_by_sample["Donor"] == donor].tolist()
    for_1 = one_donor[0] + "_mapped_1.fq.gz"
    rev_1 = one_donor[0] + "_mapped_2.fq.gz"
    for_2 = one_donor[1] + "_mapped_1.fq.gz"
    rev_2 = one_donor[1] + "_mapped_2.fq.gz"
    output_location = "Donor" + str(donor) + "assembly/"

    spades_command = f"spades.py -k 77 -t 2 --only-assembler --pe-1 1 {for_1} --pe-2 1 {rev_1} --pe-1 2 {for_2} --pe-2 2 {rev_2} -o {output_location}"
    os.system(spades_command)
    with open("PipelineProject.log", "a") as file:
        file.write(spades_command + "\n")

#do BLAST
os.system("datasets download virus genome taxon Betaherpesvirinae --refseq --include genome")
os.system("unzip ncbi_dataset.zip")

file_name = "ncbi_dataset/data/genomic.fna" #make sure this will overwrite the other dataset - idk how this works
db_name = "betaherpesvirinae"
makeblast_command='makeblastdb -in '+file_name+' -out '+db_name+' -title '+db_name+' -dbtype nucl' 

os.system(makeblast_command)

output_list = []

donors = list(set(info_by_sample["Donor"].to_list()))
for i in range(len(donors)):
    donors[i] = "Donor" + str(donors[i])

for assembly in donors: #add the filenames here
    contig_path = assembly + "assembly/contigs.fasta"
    #find the longest contig - make into it's own fasta and save
    with open(contig_path, 'r') as f:
        longest = next(SeqIO.parse(f, "fasta"))

        out_fasta = assembly + "_longest.fasta"

        with open(out_fasta, "w") as out_fas:
            SeqIO.write(longest, out_fas, "fasta")
        
    query_seqfile = out_fasta #call that fasta again

    output_file = "BLAST_" + assembly + ".csv"

    output_list.append(output_file)

    blast_command = "blastn -query " + query_seqfile + " -db betaherpesvirinae -out " + output_file + " -max_hsps 1 -outfmt '10 sacc pident length qstart qend sstart send bitscore evalue stitle'"
    os.system(blast_command)

donors = list(set(info_by_sample["Donor"].to_list()))

for i in range(2):
    blast_result = output_list[i]
    results = pd.read_csv(blast_result, header = None)
    results.columns = ['sacc', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'evalue', 'stitle', 'genome']
    results['stitle'] = results['stitle'] + ", " + results['genome']
    results = results.drop("genome", axis=1)
    with open("PipelineProject.log", "a") as file: 
        file.write("Donor " + str(donors[i]) + ": \n")
        results.to_csv(file, sep='\t', index=False, mode='a')
