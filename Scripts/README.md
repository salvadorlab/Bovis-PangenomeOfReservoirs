# Welcome to the Scripts folder of the Bovis-PangenomeOfReservoirs project.
In here, you will find the scripts used to generate the necessary data for this analysis.
Every script will be listed with a short description of what each program does, alongside the inputs and expected outputs.

## download_fastq.py
python script

Run for loop that will submit a submission job for each accession in an accession list

Input: A text file with 1 accession on each line

Output: The working directory will contain the downloaded accessions

Last Update: May 4th 2020

## read_trimming.py
python script

Trim all the isolates using trimmomatic

Input: no input, must run this program in a directory with untrimmed FASTQ files

Output: The working directory will contain trimmed FASTQ files

Last Update: May 11th 2020

## assemble_isolates.py
python script

Run for loop that will submit a submission job for each assembly

Input: no input, must run this program in a directory with trimmed FASTQ files

Output: The working directory will contain genome assembles

Last Update: May 10th 2020

## assembly_stats.py
python script

Input: no input, must run this program in a directory that contains scaffold FASTA files 

Generate Statistics for Genome Assemblies using QUAST

Last Update: May 14th 2020

## pangenome_figures.R
R script

Generate ggplot figures associated with outputs from the analysis

Input: Multiple inputs, inspect the comments inside the script to learn more about what is needed to produce each figure.

Output: Figures crafted using ggplot library 

Last Update: May 15th 2020
