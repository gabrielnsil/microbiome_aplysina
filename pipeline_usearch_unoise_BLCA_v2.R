
# Brazil 16S Aplysina fulva seasonal samples
# author: Cristiane Hardoim and Marwan
# date: "30/05/2021"
# Modified by: Gabriel Nascimento-Silva
# date_mod: 11/10/2021
# output: html_document

# This script includes the use of the new database DTGB and BLCA method. 
# Also includes chloroplast removal and deseq normalisation method.

#-------------------------------------------------------------------------------


## set the working directory

setwd("C:/Users/gabri/Desktop/microbiota_seasonal/aplysina_data")

## Requirements
library(seqinr)
library(dplyr)
library(tidyverse)
library(reshape2)

#system("ln -s 00_DataNeeded/usearch/usearch11.0.667_i86osx64 usearch11.0.667 ; chmod +x usearch11.0.667") # MacOS
#system("cp 00_DataNeeded/usearch/usearch11.0.667_i86linux64 usearch11.0.667 ; chmod +x usearch11.0.667") # Linux
shell("copy 00_DataNeeded\\usearch\\usearch11.0.667_win64.exe usearch11.0.667") # Windows
system('./usearch11.0.667') # Test

## Quality filtering

dir.create("02_PreProcessedData")
dir.create("03_MergedData")
dir.create("04_FilteredData")
dir.create("05_TrimmedData")
dir.create("06_DereplicatedData")

files = dir(path = '01_RawData/', pattern = '_R1_')
files
for (file1 in files)
{
  file2 = gsub(x = file1, pattern = '_R1_', replacement = '_R2_')
  #filename = gsub(x = file1, pattern = "_L001_R1_001.fastq.*", replacement = '', perl = T)
  filename = gsub(x = file1, pattern = "_R1_001.fastq.*", replacement = '', perl = T)
  
  #	Quality trimming
#  command <- paste('java -jar 00_DataNeeded/Trimmomatic-0.38/trimmomatic-0.38.jar PE 01_RawData/', file1, " 01_RawData/", file2, " 02_PreProcessedData/", filename, "_pF.fastq 02_PreProcessedData/", filename, "_upF.fastq 02_PreProcessedData/",filename, "_pR.fastq 02_PreProcessedData/",filename, "_upR.fastq  SLIDINGWINDOW:4:20 MINLEN:100", sep = "")
#   command <- paste('java -jar 00_DataNeeded/Trimmomatic-0.38/trimmomatic-0.38.jar PE 01_RawData/', file1, " 01_RawData/", file2, " 02_PreProcessedData/", filename, "_pF.fastq 02_PreProcessedData/", filename, "_upF.fastq 02_PreProcessedData/",filename, "_pR.fastq 02_PreProcessedData/",filename, "_upR.fastq  SLIDINGWINDOW:4:22 MINLEN:100", sep = "")
   command <- paste('java -jar 00_DataNeeded/Trimmomatic-0.38/trimmomatic-0.38.jar PE 01_RawData/', file1, " 01_RawData/", file2, " 02_PreProcessedData/", filename, "_pF.fastq 02_PreProcessedData/", filename, "_upF.fastq 02_PreProcessedData/",filename, "_pR.fastq 02_PreProcessedData/",filename, "_upR.fastq  SLIDINGWINDOW:4:25 MINLEN:100", sep = "") # Checagem de 4 nucleotídeos e mensura o score. Se não atigir 25 é excluído. MINLEN, tudo que for menor que 100 será removido.
   system(command)
  
  #	Merge paired reads (range needs to be adjusted)
  command = paste("./usearch11.0.667 -fastq_mergepairs 02_PreProcessedData/", filename, "_pF.fastq -reverse 02_PreProcessedData/", filename, "_pR.fastq -fastqout 03_MergedData/",filename, ".fastq -relabel @ -fastq_maxdiffs 5 -fastq_pctid 80 -fastq_minmergelen 230 -fastq_maxmergelen 300 -sample ", filename, sep = "")
  system(command)
  
  # #	Quality filtering, descarta qualquer sequência que tenha nucleotídeos 'n', e qualquer sequência que tem o número de erros igual ou maior 1
  command = paste("./usearch11.0.667 -fastq_filter 03_MergedData/", filename, ".fastq -fastaout 04_FilteredData/", filename, ".fasta -fastq_maxns 1 -fastq_maxee 1", sep = "")
  system(command)

  # #	Check and remove primers
  #FQ = read.fasta(file = paste0("04_FilteredData/", filename, ".fasta"), as.string = T)
  #for_primer_present = grep(pattern = "^[ATGC]{0,2}CCTACGGG[ATGC]GGC[AT]GCAG", x = FQ, ignore.case = T)
  #rev_primer_present = grep(pattern = "GGATTAGATACCC[CGT][AGT]GTAGTC[ATGC]{0,2}$", x = FQ, ignore.case = T)
  #both_present = intersect(for_primer_present, rev_primer_present)
  
  #FWD:GTGYCAGCMGCCGCGGTAA (Parada et al., 2016)-> GTG[CT]CAGC[AC]GCCGCGGTAA (PRIMERS 515 e 806)
  #REV:GGACTACNVGGGTWTCTAAT (Apprill et al., 2015)-> reverse complemente ATTAGA[AT]ACCC[CGT][ATCG]GTAGTCC[ATGC]
  #(PRIMERS 515 e 806)
  
  FQ = read.fasta(file = paste0("04_FilteredData/", filename, ".fasta"), as.string = T)
  for_primer_present = grep(pattern = "^[ATGC]{0,2}GTG[CT]CAGC[AC]GCCGCGGTAA", x = FQ, ignore.case = T)
  rev_primer_present = grep(pattern = "ATTAGA[AT]ACCC[CGT][ATCG]GTAGTCC[ATGC]{0,2}$", x = FQ, ignore.case = T)
  both_present = intersect(for_primer_present, rev_primer_present)
  
  
  #FQ1 = FQ[both_present]
  #IDs = names(FQ1)
  #SEQs = as.character(FQ1)
  #SEQs = gsub(pattern = ""^[ATGC]{0,2}CCTACGGG[ATGC]GGC[AT]GCAG", replacement = "", SEQs, ignore.case = T)
  #SEQs = gsub(pattern = "GGATTAGATACCC[CGT][AGT]GTAGTC[ATGC]{0,2}$", replacement = "", SEQs, ignore.case = T)
  #OUT = file(description = paste0("05_TrimmedData/", filename, ".fasta"), open = "w")
  #for(i in 1:length(IDs)) write(x = paste0(">", IDs[i], "\n", SEQs[i]), file = OUT)
  #close(OUT)
  
  FQ1 = FQ[both_present]
  IDs = names(FQ1)
  SEQs = as.character(FQ1)
  SEQs = gsub(pattern = "^[ATGC]{0,2}GTG[CT]CAGC[AC]GCCGCGGTAA", replacement = "", SEQs, ignore.case = T)
  SEQs = gsub(pattern = "ATTAGA[AT]ACCC[CGT][ATCG]GTAGTCC[ATGC]{0,2}$", replacement = "", SEQs, ignore.case = T)
  OUT = file(description = paste0("05_TrimmedData/", filename, ".fasta"), open = "w")
  for(i in 1:length(IDs)) write(x = paste0(">", IDs[i], "\n", SEQs[i]), file = OUT)
  close(OUT)
  
  
  #	Dereplication
  command = paste("./usearch11.0.667  -fastx_uniques 05_TrimmedData/", filename, ".fasta -fastaout 06_DereplicatedData/", filename, ".fasta -sizeout", sep = "")
  system(command)
}


#below is an example of what the data looks like when merging and trimming.


rm(list=ls())

#Check merging range - checked and changed to 250-550 (this was the most appropriate for the current study)

seq_lengths = NULL
for(f in dir(path = '05_TrimmedData/', full.names = T))
{
  print(f)
  fasta_file = read.fasta(f)
  seq_lengths = c(seq_lengths, as.numeric(lengths(fasta_file)))
}
hist(seq_lengths) # USe in line 27/28 to adjust range
quantile(seq_lengths)

#> Gabriel quantile(seq_lengths)
# 0%  25%  50%  75% 100% 
# 191  253  253  253  261 

#	Join
dir.create("07_JoinedData")

# Please decomment
#system("cat 06_DereplicatedData/*.fasta > 07_JoinedData/AllSamples.fasta") # Linux and MacOS
shell("type 06_DereplicatedData\\*.fasta > 07_JoinedData\\AllSamples.fasta") # Windows

FNA = readLines("07_JoinedData/AllSamples.fasta")
FNA[grep(pattern = ">", x = FNA, invert = T)] = toupper(FNA[grep(pattern = ">", x = FNA, invert = T)])
write(x = FNA, file = "07_JoinedData/AllSamples2.fasta")

#	Dereplication

dir.create("08_UniqueSequences")
system("./usearch11.0.667 -fastx_uniques 07_JoinedData/AllSamples2.fasta -fastaout 08_UniqueSequences/AllSamples_uniques.fasta -sizein -sizeout -strand both")


#	Generating unique sequences using UNOISE
dir.create("09_DenoisedSequences")

system("./usearch11.0.667 -unoise3 08_UniqueSequences/AllSamples_uniques.fasta -zotus 09_DenoisedSequences/AllSamples_denoised.fasta")

#	Chimera Removal
dir.create("10_UchimeReference")
system("./usearch11.0.667 -uchime2_ref 09_DenoisedSequences/AllSamples_denoised.fasta -db 00_DataNeeded/BLCA/db_GTDB_SSU/ssu_all_r95_BLCAparsed.fasta -strand plus -mode high_confidence -notmatched 10_UchimeReference/AllSamples_unoise_nc.fasta")

#	OTU table generation
# An OTU table is made by the otutab command
dir.create("11_OtuTable")
system("./usearch11.0.667 -otutab 07_JoinedData/AllSamples2.fasta -otus 10_UchimeReference/AllSamples_unoise_nc.fasta -id 0.97 -otutabout 11_OtuTable/AllSamples_unoise_otu_table1.txt") # id de 97 serve para permitir uma certa variabilidade genética entre as espécies biológicas


#instructions on how to install Blast can be accessed here (https://github.com/mhahsler/rBLAST)
#	Mapping of OTUs on Reference Database



##############################################################################################################
#################################################### BLCA ####################################################
##############################################################################################################

# add clustalo and muscle to system path (required by BLCA)
#system("cp 00_DataNeeded/clustalo/clustal-omega-1.2.3-macosx 00_DataNeeded/clustalo/clustalo ; chmod +x 00_DataNeeded/clustalo/clustalo") # MacOS
#system("cp 00_DataNeeded/muscle/muscle3.8.31_i86darwin64 00_DataNeeded/muscle/muscle ; chmod +x 00_DataNeeded/muscle/muscle") # MacOS
#Sys.setenv(PATH = paste(Sys.getenv("PATH"), "00_DataNeeded/clustalo", sep = ":"))
#Sys.setenv(PATH = paste(Sys.getenv("PATH"), "00_DataNeeded/muscle", sep = ":"))

dir.create("12_TaxAssignmentGTDB_BLCA")

# run BLCA agsinst GTDB SSU database

# If you saw a "BioPython is not detected!" error, decomment the following line to install biopython
#system(paste(system("which python3", intern = TRUE), '-m pip install biopython --user'))

# r89
# system('python3 00_DataNeeded/BLCA/2.blca_main.py -r 00_DataNeeded/BLCA/db_GTDB_SSU/GTDB_bac120_ar122_ssu_r89_BLCAparsed.taxonomy -q 00_DataNeeded/BLCA/db_GTDB_SSU/GTDB_bac120_ar122_ssu_r89_BLCAparsed.fasta -i 10_UchimeReferenceGTDB/zOTU_nc.fasta -o 12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.1.txt')

# If you saw a "BioPython is not detected!" error, decomment the following line to install biopython
system('python -m pip install biopython --user')
system('python -m pip install --upgrade pip --user')

# First step is to build a taxonomy database from the NCBI 16S microbial database
#To compile, subset the 16S Microbial database
#system('python3 00_DataNeeded/BLCA/1.subset_db_acc.py --dir 00_DataNeeded/BLCA') #it is not downloading the taxonomy file properly
# make index of your own dataset using terminal
#makeblastdb -in YourDatabase.fasta -dbtype nucl -parse_seqids -out YourDatabase


# r95
# Second step: Run your analysis with the compiled database
system('python3 00_DataNeeded/BLCA/2.blca_main.py -r 00_DataNeeded/BLCA/db_GTDB_SSU/ssu_all_r95_BLCAparsed.taxonomy -q 00_DataNeeded/BLCA/db_GTDB_SSU/ssu_all_r95_BLCAparsed -i 10_UchimeReference/AllSamples_unoise_nc.fasta -o 12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.1.txt')

  # remove taxonomic ranks from BLCA's classification and put number in parenthesis (make it easy to read)
m = 1
for (each_line in readLines(file('12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.1.txt',open="r")) ){
  each_line_split = strsplit(each_line, '\t')
  OTU_ID = each_line_split[[1]][1]
  taxonomy = each_line_split[[1]][2]
  taxonomy_split = strsplit(taxonomy, ';')
  taxonomy_no_rank = ''
  n = 1
  for (taxon in taxonomy_split[[1]]){
    if (n%%2 == 1){
      taxon_split = strsplit(taxon, ':')
      if (length(taxon_split[[1]]) ==2)
      {taxon_no_rank = taxon_split[[1]][2]} 
      else 
      {taxon_no_rank = taxon_split[[1]][1]}
      taxonomy_no_rank = paste(taxonomy_no_rank, taxon_no_rank, sep = ");")} 
    else 
    {taxonomy_no_rank = paste(taxonomy_no_rank, taxon, sep = "(")}
    n = n + 1
  }
  taxonomy_no_rank = paste(taxonomy_no_rank, ')', sep = "")
  taxonomy_no_rank = substr(taxonomy_no_rank, 3, nchar(taxonomy_no_rank))
  if (taxonomy_no_rank == "Unclassified)"){taxonomy_no_rank = "Unclassified"}
  
  taxonomy_no_rank_with_OTU = paste(OTU_ID, taxonomy_no_rank, sep = "\t")
  
  if (m == 1)
  {cat(taxonomy_no_rank_with_OTU,file='12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.2.txt',sep="\n", append=FALSE)} 
  else 
  {cat(taxonomy_no_rank_with_OTU,file='12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.2.txt',sep="\n", append=TRUE)}
  m = m +1
}


# # Merge Table and Taxonomy
dir.create("13_FinalOtuTableGTDB_BLCA")
OTU = read.delim("11_OtuTable/AllSamples_unoise_otu_table1.txt", header = T)
#I manually removed the numbers in parentesis after the taxonomy e.g. (100.0000). Then I saved the file as "AllSamples_unoise_BLCA_out.2filtered2_split.txt". here the taxonomy is split.
TAX = read.delim("12_TaxAssignmentGTDB_BLCA/AllSamples_unoise_BLCA_out.1.txt", header = F)
names(TAX) = c("X.OTU.ID", "Taxonomy")
# 
OTU_TAX = merge(OTU, TAX, by = "X.OTU.ID")
write.table(OTU_TAX, "13_FinalOtuTableGTDB_BLCA/AllSamples_unoise_otu_table_BLCA_filtered2.txt", sep = "\t", row.names = F, col.names = T, quote = F)
# 
rm(list=ls(all=TRUE))

#######################################################################
################GREENGENES TO DETECT ORGANELLS#########################
#######################################################################
dir.create("14_Unclassified_zotus")#prepare a fasta file with all unclassified zOTUs and move to the folder
dir.create("15_Unclassified_gg")

## For creation of fasta file with all unclassified zOTUS we can use awk (Linux)

# PT-BR: Primeiro vamos imprimir todas as linhas que contêm a palavra 'Unclassified'
# e enviar o resultado para um arquivo chamado unclassified_otus.txt

#awk '/Unclassified/ {print $1}' 13_FinalOtuTableGTDB_BLCA/AllSamples_unoise_otu_table_BLCA_filtered2.txt > ./unclassified_otus.txt

# PT-BR: Utilizando grep, usamos as flags -w (whole words) e -A 4 (para pular as 4 linhas de nucleotídeos)
# a opção '--no-group-separator' serve para não aparecer hífens entre os matches.
# a flag -f indica que estamos selecionando todas as zOTUS 'Unclassified' do arquivo ./unclassified_otus
# O resultado é salvo em um arquivo chamado unclassified_seqs.fasta dentro da pasta 14_unclassified_zotus
# observação: estamos usando o grep no arquivo fasta '10_UchimeReference/AllSamples_unoise_nc.fasta' pós a remoção das quimeras e identificação taxonômica no BLCA

# grep -w -A 4 -f ./unclassified_otus 10_UchimeReference/AllSamples_unoise_nc.fasta --no-group-separator > 14_unclassified_zotus/unclassified_seqs.fasta

#Install the greengenes db
#system('python3 00_DataNeeded/greengenes/1.subset_db_gg.py --dir 00_DataNeeded/greengenes')

#Run your analysis with the compiled database
system('python3 00_DataNeeded/BLCA/2.blca_main.py -r 00_DataNeeded/BLCA/gg_13_5_taxonomy.taxonomy -q 00_DataNeeded/BLCA/gg_13_5 -i 14_unclassified_zotus/unclassified_seqs.fasta -o 15_Unclassified_gg/unclassified_greengenes_taxonomy.txt')


#######################################################################
####################NORMALIZED FOR ALPHA-DIVERSITY#####################
#######################################################################
# https://drive5.com/usearch/manual/cmd_otutab_rare.html

# What is Alpha-Diversity? It is the mean diversity of species in different sites or habitats within a local scale.
# Subsamples ("rarefies") an OTU table to a fixed number of reads per sample using random subsampling without replacement. 
# This is likely the best strategy because preserves the shape of the abundance distribution in each sample more accurately 
# than systematic rounding as used in the obsolete otutab_norm command.
# The author of USEARCH recommend the use of this command to normalize samples to the same number of reads so that they are comparable to each other.

# Samples which have < sample_size reads are discarded.
# Code example: usearch -otutab_rare otutab.txt -sample_size 5000 -output otutab_5k.txt

dir.create("16_Normalized")

# PT-BR: Para normalizar a diversidade alfa entre as amostras para torná-las comparáveis é necessário
# remover todas as zOTUs que são classificadas como 'Mitochondria' e 'Chloroplast' 
# Além disso, devemos remover a última coluna que contém a classificação taxonômica



# awk '/[Cc]hloroplast/ || /[Mm]itochond/ {print $1}' 15_Unclassified_gg/unclassified_greengenes_taxonomy.txt > ./15_Unclassified_gg/zOTUs_removal_file.txt

# Você pode visualizar no terminal como ficará seu arquivo da seguinte forma:

# grep -w -v -f ./15_Unclassified_gg/zOTUs_removal_file.txt 13_FinalOtuTableGTDB_BLCA/AllSamples_unoise_otu_table_BLCA_filtered2.txt | awk -F '\t' -v OFS='\t' '{$NF=""; print $0}' | column -t -s $'\t' | less -S

# Salvar o arquivo final direto na pasta 16_Normalized:
# grep -w -v -f ./15_Unclassified_gg/zOTUs_removal_file.txt 13_FinalOtuTableGTDB_BLCA/AllSamples_unoise_otu_table_BLCA_filtered2.txt | awk -F '\t' -v OFS='\t' '{$NF=""; print $0}' > ./16_Normalized/finalzOTUsTable_not_normalized.txt

# Agora só editar o primeiro campo do arquivo para ele ser aceito no USEARCH
# sed 's/X.OTU.ID/#OTU.ID/' ./16_Normalized/finalzOTUsTable_not_normalized.txt | sed 's/\t$//g' > ./16_Normalized/finalzOTUsTable_not_normalized_formated.txt

install.packages('janitor')
library(janitor)
finalTable = read_tsv('./16_Normalized/finalzOTUsTable_not_normalized_formated.txt')

finalTable2 = finalTable %>%
  adorn_totals("row")

min_size = apply(tail(finalTable2, 1), 1, FUN = min)

typeof(min_size)

# system("./usearch11.0.667 -otutab_rare 16_Normalized/finalzOTUsTable_not_normalized.txt -sample_size 500 -output 16_Normalized/normalized.txt")

system(paste("./usearch11.0.667 -otutab_rare 16_Normalized/finalzOTUsTable_not_normalized_formated.txt -sample_size ", min_size, "-output 16_Normalized/normalized.txt"))


# Confirmação da amostragem:

normalized = read_tsv('./16_Normalized/normalized.txt')
normalized = normalized %>%
  adorn_totals("row")

tail(normalized, 1)

rm(list = ls())

