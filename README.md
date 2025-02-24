# A Pipeline to Analyze Paired-End Reads from HCMV Transcriptomes 

# Getting Started With the Pipeline:
## Step 1: Download all files in this repository to a directory on your local/remote machine
### Where to Find Test Data: 
Test data is in the folder sample_data. It includes paired-end reads from samples SRR5660030, SRR5660033, SRR5660044, and SRR5660045. Keep it together in the directory sample_data. You will use this directory as input for the wrapper script.

Metadata for these samples is found in sample_metadata.csv. This includes the donor number for each sample and the condition (2 or 6 days post-infection). Metadata must include column "Sample", which specifies the sample name attached to our corresponding paired-end fastq files. If paired end fastq files are SRR5660044_1.fastq and SRR5660044_2.fastq, the "Sample" must be SRR5660044. Metadata also must include columns "Donor" and "Condition."

## Step 2: Make sure all dependencies are installed

Command Line: <br>
  -- Python <br>
  -- R <br>
  -- BLAST+ <br>
  -- SPades <br>
  -- Bowtie2 <br>
  -- Kallisto <br>
  -- NCBI datasets tools <br>

Python Packages:  <br> 
  -- Biopython (version 1.83)  <br> 
  -- Pandas (version 2.1.4)  <br> 
  -- sys  <br> 
  -- argparse  <br> 
  -- os <br>
  -- glob <br>

R Packages: <br>
  -- Sleuth <br>

## Step 3: Use the command line to run the wrapper script as follows: 

Before running the script, clone this repository and cd into it:
```
git clone https://github.com/estagaman/Pipeline_Elise_Stagaman.git

cd Pipeline_Elise Stagaman
```

### Running the Wrapper Script:
Input flags:
  -m path to the metadata csv
  -d path to the directory containing paired-end fastq files

```
#for running test data: 
python3 wrapper_script.py -m <<metadata>> -d <<sample data>>

#for running your own data:
python3 wrapper_script.py -m sample_metadata.csv -d sample_data &
```

If you want to be able to exit the server while the script is running, use command: 

```
nohup python3 wrapper_script.py -m sample_metadata.csv -d sample_data &
```
