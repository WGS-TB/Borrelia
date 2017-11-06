# Pre-requisites

1) R
2) VelvetOptimiser: Please download it here: https://github.com/tseemann/VelvetOptimiser
3) Cookicutter    : Please download it here: https://github.com/ad3002/Cookiecutter



# Inputs

1) A main folder which contains 
1) For every sample i, we assign a folder S_i which consists pair end reads as two seprate FASTQ files. All Folders S_i's are stored in a main folder M.  

2) All sequences of the reference alleles as a nexus file.

# How to run the pipeline

1) The required scripts and files are in PIPELINE folder.
3) Make the script PIPELINE.sh executable by doing as:
```
chmod +x PIPELINE.sh
```
2) Run the following script in the PIPELINE folder as follows:
```
./PIPELINE.sh   dir    path
```
dir  : The directory of the main folder M.

path : The path of the nexus file which contains all reference alleles. 

# Output
1) The pipeline create a folder named RESULT in PIPELINE folder. 
2) The folder RESULT consists a folder S_i for every sample. Every folder S_i consists two fastq files and a folder named leftOver.
The fastq files are reads likely to come from the gene of interest from all the reads in sample S_i. The folder leftOver consists two fastq files and a folder named Results. The fastq files are the unassinged reads. The folder Results contains two folders: Fraction, assembly_result. The folder Fraction consists the fraction of alleles which being present in the sample in a csv file.  The folder assembly_result consists the output of velvetompimiser over unassgined reads; Moreover, it contains Longestcontigs.fa which is the longest conting assembled from the unassgined reads, and AlignScoreRef.fa, AlignScoreRef_RC.fa which are the alignment scores of the longest contig with the reference alleles and their reverse complement respectively.
