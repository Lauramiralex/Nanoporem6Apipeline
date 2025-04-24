# Pipeline for Nanopore Data Analysis

This pipeline is designed for the analysis of Nanopore sequencing data, specifically POD5 files.  
It processes experimental and control group data by performing basecalling and genome alignment 
followed by differentially methylated region (DMR) analysis.  

# Introduction
## Biological Motivation

High-resolution DNA sequencing and analysis are fundamental to modern biological research and clinical diagnostics. By extracting DNA material and subjecting it to high-throughput sequencing, we gain detailed insights into the genetic and epigenetic composition of a sample. This provides a foundation for studying gene regulation, disease mechanisms, evolutionary processes and  developmental biology.

Sequencing DNA allows researchers to determine the nucleotide sequence of specific genomic regions or even entire genomes. This information is essential for identifying genetic variations such as insertions or deletions, characterizing genes and their regulatory elements and investigating genome-wide changes in different biological conditions or individuals.

## Nanopore

Nanopore sequencing is a third-generation sequencing technology that enables the direct reading of DNA or RNA strands as they pass through a biological or synthetic nanopore embedded in a membrane. 
As nucleotides pass through the nanometer-scale hole which is embedded in a membrane (pore), they cause characteristic changes in an ionic current. These changes are detected in real-time and translated into nucleotide sequences through basecalling algorithms.

Nanopore sequencing offers several distinct advantages over traditional sequencing technologies, particularly in its ability to generate long reads and directly detect native DNA or RNA molecules. This does not only simplify experimental workflows but also expands the range of biological insights that can be obtained from a single sequencing run.
A distinct advatage of nanopore sequencing is its capacity to detect base modifications directly from the ionic current signal eliminating the need for chemical conversion protocols which are traditionally used in DNA methylation analysis.
This simplifies the experimental workflow and enables to preserve the biological context, including strand orientation and base modification. Moreover, by eliminating PCR amplification, nanopore sequencing reduces bias and enables the analysis of native DNA and RNA molecules.

## m6A

m6A is a chemical modification in which a methyl group is added to the nitrogen at position 6 of adenosine. It is the most abundant internat modification found in messenger RNA and long non-coding RNA in eukaryotic cells. m6A is typically found in conserved sequence motifs and is frequently enriched near stop codons, within long exons and in 3' untranslated regions. Its writer mark is deposited by a methyltransferase complex and it can be removed by demethylases such as FTO. Specialized m6A-binding proteins recognize and interprete this mark to influence downstream RNA metabolism.
Functinally, m6A plays a central role in post-transcriptional gene regulation influencing critical processes such as mRNA stability, translation efficiency and RNA splicing. Changes in its patterns have been associated with cancer, viral infections and neurological disorders.


# Pipeline Overview / Features
Basecalling: Converts raw electrical signal data from Nanopore sequencing into nucleotide sequences. This happens with Dorado which is an open-source basecaller 
developed by Oxford Nanopore Technologies (ONT) for fast and accurate basecalling of POD5 files.  

Genome Alignment: Maps the sequences of basecalled reads efficiently to a reference genome while preserving information about modified bases (e.g., methylation). 
This happens with Modkit.  

Pileup: Operation that summarizes aligned sequencing reads at each position of a reference genome.  

Differentially Methylated Region (DMR) Analysis: DMR analysis helps identify regions of the genome where methylation patterns significantly differ between experimental and control groups.  

All applications of the pipeline can also be used individually to rerun specific sections.  

# Installation and Dependencies
To successfully run your Nanopore data analysis pipeline, the following tools must be installed.   
Follow the steps below to install the necessary tools on a Linux system.  

## Installing Conda
Go to the Anaconda [download page](https://www.anaconda.com/download/)  
Download the Anaconda3 (Python 3 version) 64-bit installer for Linux.  
1. Open a terminal and run the following command to download and install Anaconda:  
   ``` wget https://repo.anaconda.com/archive/Anaconda3-2025.03-Linux-x86_64.sh ```  
2. After the download completes, run the following command to start the Anaconda installation:  
   ``` bash Anaconda3-2025.03-Linux-x86_64.sh ```  
3. Follow onscreen instructions such as:
   -Accept the license agreement.  
   -Choose the installation directory.  
   -Allow the installer to initialize Anaconda by adding it to your shell startup files (this is recommended).  
4. Verify installation:  
   Check if Anaconda (and Conda) was successfully installed by running:  
   ```conda --version```  
   You should now see the conda version which confirms the successfull installation

For an in-depth installation guide please look into the [Anaconda Documentation](https://www.anaconda.com/docs/getting-started/anaconda/install)
### Creating Conda Environment
After installing Conda one should create a dedicated environment for the project in which one cann install all further necessary packages.
```conda create --name nanoporepipe python=3.8```  
Activate the environment with:  
```conda activate nanoporepipe```  
The activation step needs to be done everytime one wants to work on the project.  

## Installing Bioconda
Bioconda is a collection of bioinformatics software packages available through Conda. Configure your Conda channels:
``` #Add the required channels
conda config --add channels conda-forge
conda config --add channels bioconda
conda config --add channels defaults
 ```   
For an in-depth installation guide please look into the [Bioconda Documentation](https://bioconda.github.io)
## Installing Snakemake
Snakemake is used for workflow management. Install Snakemake in terminal:  
```conda install -c conda-forge snakemake ```  
For an in-depth installation guide please look into the [Snakemake Documentation](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)

## Installing the Tools
### Modkit
Modkit is a bioinformatics tool available via Bioconda. Install Modkit in terminal:  
```conda install -c bioconda modkit```  
### Dorado
Dorado is a tool for nanopore data analysis provided by Oxford Nanopore Technologies.Install Dorado in terminal:  
```conda install -c bioconda dorado```  
For an in-depth installation guide please look into the [Dorado guide by Nanopore](https://github.com/nanoporetech/dorado)
### Samtools
Samtools is a widely used tool for processing SAM/BAM files. Install Samtools in terminal:  
```conda install -c bioconda samtools```   
For an in-depth installation guide please look into the [Samtools Documentation](https://www.htslib.org)
### Bgzip and Tabix
Bgzip is a version of gzip optimized for handling large files.  
Tabix is a tool for indexing compressed files, especially useful with VCF and GFF files.  
Install Bgzip and Tabix in terminal:  
```
conda install -c bioconda bgzip
conda install -c bioconda tabix
```  
## Verify Installation
To check the installed tools and verify they work correctly run following commands in terminal:  
```
snakemake --version
samtools --version
bgzip --version
tabix --version
```  
If these commands return version numbers, the tools have been installed successfully.  

# Setup Instructions for Usage
To successfully run the scripts in this repository, follow the steps below to set up the required folder structure and input files.  

## 1. Setup Working Structure
Before running the script, the following structure should be in place:  
   - The script and the config folder(download from repository) have to be in the same directory   
   
## 2. Organizing Input Files
The script requires POD5 files to function correctly. Follow these rules when organizing them:      
   - Ensure that the folder containing the POD5 files is uniquely named according to the experiment      name.   
## 3. Customizing Config File
In your config file you will find different contents which need customization:     
   - config.yaml :
        - provide the **workdirectory** under "workdirection", from the workdirectory it is                       critical to be able to access all further needed folders from chosen directory.   
        **This will be the start of all further provided paths**   
        - provide the desired **resultdirectory** in which all results will be accessible.   
          The result directorys path needs to be provided starting by workdirectory.   
        **example: workdirectory : "/vol/spacex/nanopore" which lies before result directory
              resultdirectory : "test_pipeline/final"**  
        - provide the path for the **genome fasta and genome.fai**  
        - provide the path and modelname of used **basecall_model**  
            ```basecall_model: 'your_model_name_here' ```    
  
## 4. Configuring TSV Files
Two TSV files inside the `config` folder need to be correctly filled:  
### allfiles.tsv
In the `allfiles` column, list all experiment names and control group names and their individual path starting from the workdirectory.  
### sample_sheet.tsv
In the `experiment` column, list all experiment names and their individual path starting from the workdirectory.  
In the `control` column, list all control group names and their individual path starting from the workdirectory.  

#### ! important: The listed names of experiments and control groups must be identical to correlating individual pod5 files !

With this setup, the script should run correctly.

# How to run the script
Befor running the script make sure to follow the Setup instructions.  
After activating the Environment ( see step 'Creating Conda Environment') one is ready to let the script run.  
While being in the `main folder` the command  
   ```snakemake --cores all -s scriptname.py --jobs 1 ```  
runs the desired script on the previous determined pod5  files.  

In case one wants to only run a certain part, make sure that the required data is correctly named in the folders from the step bevor.  
Run the same command with the name of the script for the required part.

### completerun.py
For running `completerun.py` one needs the pod5 files in a file which path is provided in allfiles.tsv and sample_sheet.tsv.  
The script produces data from basecalling, alignment, pileup and dmr while also providing summarys after each step.  
The dataformats produced by this script include .bam files, .txt files, .tsv files, .bed files and .log files.

### basecaller.py
For running `basecaller.py` one needs the pod5 files in a file which path is provided in allfiles.tsv and sample_sheet.tsv.  
The script produces data from basecalling while also providing a summary.  
The dataformats produced by this script include .bam files, .txt files, .tsv files and .log files.

### aligner.py
For running `aligner.py` one needs the pod5 files in a file which path is provided in allfiles.tsv and sample_sheet.tsv.  
The script produces data from basecalling and alignment while also providing summarys after each step.  
The dataformats produced by this script include .bam files,.txt files, .tsv files and .log files.

### pileup.py
For running `pileup.py` one needs the pod5 files in a file which path is provided in allfiles.tsv and sample_sheet.tsv`.  
The script produces data from basecalling, alignment and modkit pileup while also providing summarys after each step.  
The dataformats produced by this script include .bam files,.txt files, .bed files, .bai files, .tsv files and .log files.  

### only_dmr.py
For running `only_dmr.py` one needs the sorted.bam files in a file `{path to pod5 files}/{name of experiment|control}/results/pileup/pileup.bed.gz`. This is achievable with the `aligner.py` script. 
The script produces data from modkit extract and dmr analysis, while also providing summarys after each step.  
The dataformats produced by this script include .txt files, .tsv files and .log files.  

### only_align.py
For running `only_align.py` one needs the basecalled/{}.bam files in a file called `{path to pod5 files}/{name of experiment|control}/results/basecalled/{name of experiment/control}.bam`. This is achievable with the `basecaller.py` script.  
The script produces data from dorado aligner  while also providing a summary.  
The dataformats produced by this script include .txt files, .bam files, .tsv files and .log files.


### only_pileup.py
For running `only_pileup.py` one needs the results/aligned/indexed.bam.bai and results/aligned/sorted.bam files in a file called `{path to pod5 files}/{name of experiment|control}/results/aligned`. This is achievable with the `aligner.py` script or the `onlyalign.py`script.   
The script produces data from modkit pileup  while also providing a summary.  
The dataformats produced by this script include .txt files, .bai files, .bed files, .bam files and .log files.






