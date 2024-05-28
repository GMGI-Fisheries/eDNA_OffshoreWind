# Environmental DNA (eDNA) workflow for vertebrate targets with 12S

See https://github.com/emmastrand/GMGI_Notebook/blob/main/posts/2024-02-27%20eDNA%20workflow%2012S.md for more details on the programs used.

### Workflow steps 

0. FastQC (`00-fastqc.sh`)  
1. Nf-core Ampliseq: Implements CutAdapt and DADA2 pipelines (`01-ampliseq.sh`).     
2. Taxonomic identification with blastn and taxonkit (`02-tax_ID.sh`).  
3. Taxonomic assignment and data table preparation (`03-datatable_prep_report_generation.Rmd`)

## 00. FastQC and MultiqQC 

I previously set up conda environments with multiqc installed. I should reorganize these environments, but for now use `haddock_methylation` by:  
- `source ~/../../work/gmgi/miniconda3/bin/activate`  
- `conda activate haddock_methylation`  

To run a slurm array, create a rawdata list `ls -d /work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/*.gz > rawdata`. I did August and September later:  
- `ls -d /work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/August/*.gz > rawdata_August` (sbatch --array=1-92 00-fastqc_August.sh)
- `ls -d /work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/September/*.gz > rawdata_September` (sbatch --array=1-98 00-fastqc_September.sh)

`00-fastqc.sh`:  

```
#!/bin/bash
#SBATCH --error=fastqc_output/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=fastqc_output/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=23:00:00
#SBATCH --job-name=fastqc
#SBATCH --mem=5GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

module load OpenJDK/19.0.1 ## dependency on NU Discovery cluster 
module load fastqc/0.11.9

raw_path="/work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data"
dir="/work/gmgi/Fisheries/eDNA/offshore_wind2023/QC/raw_fastqc"

## File name based on rawdata list
mapfile -t FILENAMES < ${raw_path}/rawdata
i=${FILENAMES[$SLURM_ARRAY_TASK_ID]}

## FastQC program
fastqc ${i} --outdir ${dir}
```

To run slurm array = `sbatch --array=0-290 00-fastqc.sh`.

Once complete, `cat *output.* > ../fastqc_output.txt` to create one file with all the output. The length of this file should be 290. 

`00-mutliqc.sh`: 

```
#!/bin/bash
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=23:00:00
#SBATCH --job-name=multiqc
#SBATCH --mem=5GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

dir="/work/gmgi/Fisheries/eDNA/offshore_wind2023/QC/raw_fastqc"
multiqc_dir="/work/gmgi/Fisheries/eDNA/offshore_wind2023/QC/multiqc/"

multiqc --interactive ${dir} -o ${multiqc_dir} --filename multiqc_raw.html
```

Check this output before proceeding. 

## 01. Nf-core ampliseq pipeline 

Testing 5 months of data for July - November of offshore wind eDNA 12S metabarcoding data. 

### Metadata 

#### 12S primer sequences (required)

Below is what we used for vertebrate testing between Riaz and Degenerate. Ampliseq will automatically calculate the reverse compliment and include this for us.

MiFish 12S amplicon F: ACTGGGATTAGATACCCY (if adding RC, ...CTAGAGGAGCCTGTTCTA)      
MiFish 12S amplicon R: TAGAACAGGCTCCTCTAG (if adding RC, ...RGGGTATCTAATCCCAGT)     

#### Metadata sheet (optional) 

The metadata file has to follow the QIIME2 specifications (https://docs.qiime2.org/2021.2/tutorials/metadata/). Below is a preview of the sample sheet used for this test. Keep the column headers the same for future use. The first column needs to be "ID" and can only contain numbers, letters, or "-". This is different than the sample sheet. NAs should be empty cells rather than "NA". 

#### Samplesheet information (required)

This file indicates the sample ID and the path to R1 and R2 files. Below is a preview of the sample sheet used in this test. File created on RStudio Interactive on Discovery Cluster using (`create_metadatasheets.R`).  

- sampleID (required): Unique sample IDs, must start with a letter, and can only contain letters, numbers or underscores (no hyphons!).  
- forwardReads (required): Paths to (forward) reads zipped FastQ files  
- reverseReads (optional): Paths to reverse reads zipped FastQ files, required if the data is paired-end  
- run (optional): If the data was produced by multiple sequencing runs, any string  

| July_501_1_B   | /work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/July-501-1BottomDegen_R1.fastq.gz  | /work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/July-501-1BottomDegen_R2.fastq.gz  |
|----------------|-----------------------------------------------------------------------------------------|-----------------------------------------------------------------------------------------|
| July_501_1_S   | /work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/July-501-1SurfaceDegen_R1.fastq.gz | /work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/July-501-1SurfaceDegen_R2.fastq.gz |
| July_501_2_B   | /work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/July-501-2BottomDegen_R1.fastq.gz  | /work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/July-501-2BottomDegen_R2.fastq.gz  |
| July_501_2_S   | /work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/July-501-2SurfaceDegen_R1.fastq.gz | /work/gmgi/Fisheries/eDNA/offshore_wind2023/raw_data/July-501-2SurfaceDegen_R2.fastq.gz |


### Workflow Overview 

![](https://raw.githubusercontent.com/nf-core/ampliseq/2.8.0//docs/images/ampliseq_workflow.png)

#### Pipeline summary

nf-core's ampliseq pipeline uses the following steps and programs. All programs are loaded by ampliseq. 
- Quality control (FastQC)   
- Trimming (Cutadapt)  
- Infer Amplicon Sequence Variants (DADA2)  
- Predict ribosomal RNA sequences from ASVs (Barnap) 
- Taxonomic classification (DADA2)  
- Pipeline QC summaries (MultiQC)  
- Pipeline summary report (R Markdown)  

#### Container information 

Singularity is the container loaded onto NU's cluster: https://sylabs.io/docs/. 

#### Pipeline updates

> When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you’re running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```
module load nextflow/23.10.1
nextflow pull nf-core/ampliseq
```

Notes:   
- Checked for most recent version 2024-03-29.   
- 2024-04-01 tried while short partition was down `srun --partition=short --nodes=1 --cpus-per-task=2 --ntasks=24 --pty /bin/bash`. Nodes were available soon after and ran as a slurm script. 

### Slurm script 

Slurm script to run: 

`01-ampliseq.sh`:

```
#!/bin/bash
#SBATCH --error="%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output="%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=ampliseq_OW
#SBATCH --mem=20GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

# singularity module version 3.10.3
module load singularity/3.10.3

# nextflow module loaded on NU cluster is v23.10.1
module load nextflow/23.10.1

#set paths 
metadata="/work/gmgi/Fisheries/eDNA/offshore_wind2023/metadata" 

cd /work/gmgi/Fisheries/eDNA/offshore_wind2023

nextflow run nf-core/ampliseq -resume \
   -profile singularity \
   --input ${metadata}/samplesheet.csv \
   --FW_primer "ACTGGGATTAGATACCCY" \
   --RV_primer "TAGAACAGGCTCCTCTAG" \
   --outdir ./results \
   --trunclenf 100 \
   --trunclenr 100 \
   --trunc_qmin 25 \
   --max_len 200 \
   --max_ee 2 \
   --min_len_asv 100 \
   --max_len_asv 115 \
   --sample_inference pseudo \
   --skip_taxonomy \
   --ignore_failed_trimming
```

Spaces are not allowed after each \ otherwise nf-core will not read the parameter. 

**Notes**  
- A couple samples have very little sequences... Will have to come back to this? I needed to add `--ignore_failed_trimming` in order for the script to continue running after a failed sample. 
- July - November took ~50 minutes for `01-ampliseq.sh` to run with all flags above. 


# 2. Taxonomic Identification

Taxonkit has already been downloaded and databases have already been made into ncbi search-able databases.  

Until the -remote option is fixed within scripts (waiting on NU response), I'm running Mito and GMGI within slurm script and blastn separately. 

Checking Mito db is the most updated version: https://mitofish.aori.u-tokyo.ac.jp/download/. 

### Slurm script

The following script takes the `ASV_seqs.len.fasta` file from DADA2 output (`asv_length_filter/`) and uses `blastn` to compare sequences to three databases: Blast nt, Mitofish, and our in-house GMGI database.

Make a directory called `BLASToutput` within the project directory.

`02-taxonomicID.sh`: 

```
#!/bin/bash
#SBATCH --error=output_messages/"%x_error.%j" #if your job fails, the error report will be put in this file
#SBATCH --output=output_messages/"%x_output.%j" #once your job is completed, any final job report comments will be put in this file
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --time=20:00:00
#SBATCH --job-name=tax_ID
#SBATCH --mem=30GB
#SBATCH --ntasks=24
#SBATCH --cpus-per-task=2

## load modules needed
module load ncbi-blast+/2.13.0

# set path to ASV_seqs.fasta
## needs to be changed for every project
fasta="/work/gmgi/Fisheries/eDNA/offshore_wind2023/results/asv_length_filter"
out="/work/gmgi/Fisheries/eDNA/offshore_wind2023/BLASToutput"
db="/work/gmgi/Fisheries/databases/12S/reference_fasta"
taxonkit="/work/gmgi/Fisheries/databases/taxonkit"

#### DATABASE QUERY ####
### NCBI database 
blastn -remote -db nt \
   -query ${fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_NCBI.txt \
   -max_target_seqs 10 -perc_identity 100 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid   sscinames   staxid pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

## Mitofish database 
blastn -db ${db}/Mitofish.fasta \
   -query ${fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_Mito.txt \
   -max_target_seqs 10 -perc_identity 100 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid  pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

## GMGI database 
blastn -db ${db}/GMGIVertRef.fasta \
   -query ${fasta}/ASV_seqs.len.fasta \
   -out ${out}/BLASTResults_GMGI.txt \
   -max_target_seqs 10 -perc_identity 98 -qcov_hsp_perc 95 \
   -outfmt '6  qseqid   sseqid   pident   length   mismatch gapopen  qstart   qend  sstart   send  evalue   bitscore'

############################

#### TAXONOMIC CLASSIFICATION #### 
## creating list of staxids from all three files 
awk -F $'\t' '{ print $4}' ${out}/BLASTResults_NCBI.txt | sort -u > ${out}/NCBI_sp.txt

## annotating taxid with full taxonomic classification
cat ${out}/NCBI_sp.txt | ${taxonkit}/taxonkit reformat -I 1 -r "Unassigned" > NCBI_taxassigned.txt
```

### Output 

- `BLASTResults_GMGI.txt`: GMGI's curated database results    
- `BLASTResults_Mito.txt`: MitoFish database results    
- `BLASTResults_NCBI.txt`: NCBI results    
- `NCBI_taxassigned.txt`: full taxonomic classification for all staxids found in NCBI output. 

Files to use within `/work/gmgi/Fisheries/databases/12S`: `taxonomic_classification_fishbase.csv` and `taxonomic_classification_mitofish.csv`. 

# 3. Datatable preparation and report generation 

I'm giving each database a rank - 1.) GMGI, 2.) MitoFish, 3.) NCBI. This way in the R script, we can give priority to the taxonomic match according to rank. 

RStudio on my desktop (so the output will go into our Box folder): Open Offshore Wind R Project then `03-datatable_prep_report_generation.Rmd`. 

Sections:  
1. Data input: Load libraries, metadata, raw ASV_table, NCBI results from GMGI, Mito, NCBI, speices_lookup tab, and NCBI_taxassigned.   
2. BLAST output data table preparation: Editing BLAST output dataframes so each df has a 'Species_name' column and rank (1-3).  
3. Taxonomic assignment: Assigning the Species_name to the ASV_table based on rank order above (1:GMGI, 2:Mito, 3:NCBI).    
4. Editing taxonomic assignment if necessary: If an ASV has multiple hits or NCBI/Mito assignment includes multiple hits, then edit if needed.  
5. Filtering: First remove any reads that make up less than 0.1% of ASV total per sample, THEN collapsing species name together (5 cod groups into 1 cod group).  
6. Reports: Creating graphics of filtering and ASV statistics.   

### Output 

- `SampleReport_taxonomic_mismatch.xlsx`: List of taxonomic mismatch between GMGI, Mito, and NCBI list. Review this to confirm OK with GMGI assignment over Mito and NCBI. Big question here is: Are there any entries in Mito or NCBI that should be in our GMGI database?  
- `SampleReport_FilteringStats.png`: Graphics of filtering information (e.g., number of reads passed through each cutadapt and DADA2 step)  
- `SampleReport_unassigned.png`: Graphics of assigned ASVs to taxonomic units and information on unknown detections.   
- `SampleReport_Speciesname_abundance.png`: Mean read count per common name bin for each Month of offshore wind as first glance at data. 
- `Results.xlsx`: Metadata with SampleID, species name and common name with number of reads  
- `Results_matrix.xlsx`: No metadata, each row is a species group with a column for each sample ID and number of reads population the df. 

Input needed from Fisheries team member: Compare the list of mismatching scientific names from GMGI vs. Mito vs. NCBI:  
- Does this list make sense?  
- GMGI's annotation is given preference, is this OK after seeing this list?  
- Are there any spp from Mito or NCBI that should be added to GMGI's? 

