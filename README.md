# A Pipeline to Analyze Paired-End Reads from HCMV Transcriptomes 

## Where to Find Test Data: 
Test data is in the folder test_data. It includes paired-end reads from samples SRR5660030, SRR5660033, SRR5660044, and SRR5660045. 

Metadata for these samples is found in info_by_sample.txt. This tells you which samples are from which donors/patients. It also includes information on whether the sample was taken 2 or 6 days post-infection. 

# Getting Started With the Pipeline:
## Step 1: Download all data in this repository to a directory on your local/remote machine

## Step 2: Use the command line to run the wrapper script as follows: 

```
python3 wrapper_script.py
```

If you want to be able to exit the server while the script is running, use command: 

```
nohup python3 wrapper_script.py &
```
