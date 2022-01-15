
setwd("C:/Users/gabri/Desktop/microbiota_seasonal/aplysina_data")
setwd("C:/Users/gabri/Desktop/microbiota_seasonal/Analyses_after_pipeline/05_bubble_chart")

# install.packages("svglite")
# install.packages("Cairo")

library(janitor)
library(tidyverse)
library(ggplot2)
library(Cairo)
library(reshape2)
library(viridis)


#citation(package = "ggplot2")
#citation(package = "reshape2")

grep -w -v -f ./15_Unclassified_gg/zOTUs_removal_file.txt 13_FinalOtuTableGTDB_BLCA/AllSamples_unoise_otu_table_BLCA_filtered2.txt > ../Analyses_after_pipeline/05_bubble_chart/raw_zotuTable

sed 's/^Zotu/ASV/g' ../Analyses_after_pipeline/05_bubble_chart/raw_zotuTable.txt | sed 's/X\.OTU\.ID/ASV_ID/' > ../Analyses_after_pipeline/05_bubble_chart/ASV_table.txt

awk -F '\t' '{print $NF}' ../Analyses_after_pipeline/05_bubble_chart/ASV_table.txt | awk -F ';' '{print $3}' | awk -F ':' '{print $2}' | sed 's/^$/Unclassified/g' | sed '1,1 s/Unclassified/Phylum/' > ../Analyses_after_pipeline/05_bubble_chart/phyla_for_appending.txt

awk -F '\t' '{$1="";$NF =""; print}' ../Analyses_after_pipeline/05_bubble_chart/ASV_table.txt | sed 's:\(^ \| $\)::g' | tr " " '\t' > ../Analyses_after_pipeline/05_bubble_chart/samples.txt

paste ../Analyses_after_pipeline/05_bubble_chart/phyla_for_appending.txt ../Analyses_after_pipeline/05_bubble_chart/samples.txt > ../Analyses_after_pipeline/05_bubble_chart/ASV_phylumTable.txt



phylumTable <- read_tsv("./ASV_phylumTable.txt")


phylumTable2 <- phylumTable %>%
  group_by(Phylum) %>%
  summarise(across(where(is.numeric), sum))

totals = colSums(phylumTable2[,-1])

relAbund <- mapply('/', phylumTable2[,-1]*100, totals)

relAbund <- as.data.frame(round(relAbund, digits=6))

relAbund <- tibble(phylumTable2[,1], relAbund)

relAbund_analysis <- relAbund %>%
  rowwise() %>% 
  mutate( mean = mean(c_across(where(is.numeric))))


write.table(relAbund_analysis, file='relAbund_for_formatting.tsv', quote=FALSE, sep='\t', row.names = FALSE)


## formatação básica feita no Excel ou Bloco de Notas

input <- read_tsv("./relAbund_formatted.txt")

input <- input %>%
  rowwise() %>%
  mutate(Autumn = mean(c_across(starts_with('Autumn')))) %>%
  mutate(Winter = mean(c_across(starts_with('Winter')))) %>%
  mutate(Spring = mean(c_across(starts_with('Spring')))) %>%
  mutate(Summer = mean(c_across(starts_with('Summer')))) 


input <- input %>%
  group_by(Phylum) %>%
  filter(mean > 1) %>%
  select(tail(names(.), 4))



input2 <- as.data.frame(t(input)) %>%
  rownames_to_column() %>%
  row_to_names(row_number = 1) %>%
  rename("Sample" = "Phylum") %>%
  ungroup() 

cols <- names(input2)[-1]

input2[cols] <- lapply(input2[cols], as.numeric)
sapply(input2, class)



phylum = input2

rm(list = ls())



########################################################################################################
#################################################PHYLUM#################################################
########################################################################################################
#upload your data to R - exchange "Your_csv_file.csv" with the name of 
#your csv file
phylum = read.csv("Pasta1.csv", header = TRUE)

#convert data frame from a "wide" format to a "long" format
pcm_phylum = melt(phylum, id = c("Sample"))

#pick colours
colours = viridis_pal(option = "D")(13)



#Keep the order of samples from your excel sheet:
pcm_phylum$Sample <- factor(pcm_phylum$Sample,levels=unique(pcm_phylum$Sample))

#Plot it! For a bubble plot, you are using geom_point and scaling the size to your value (relative abundance) column.
xx = ggplot(pcm_phylum, aes(x = Sample, y = variable)) +
  geom_point(aes(size = value, fill = variable), alpha = 0.75, shape = 21) +
  scale_size_continuous(limits = c(0.000001, 100), range = c(1,17), breaks = c(1,10,50,75)) +
  labs( x= "", y = "", size = "Relative Abundance (%)", fill = "") +
  theme(legend.key=element_blank(),
        axis.text.x = element_text(colour = "black", size = 12, face = "bold", angle = 90, vjust = 0.3, hjust = 1),
        axis.text.y = element_text(colour = "black", face = "bold", size = 11),
        legend.text = element_text(size = 10, face ="bold", colour ="black"),
        legend.title = element_text(size = 12, face = "bold"),
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.position = "right") +
  scale_fill_manual(values = colours, guide = "none") +
  scale_y_discrete(limits = rev(levels(pcm_phylum$variable)))
xx
#Warning message:
#It is deprecated to specify `guide = FALSE` to remove a guide. Please use `guide = "none"` instead.

ggsave("Bubble_phylum.svg")
ggsave("Bubble_phylum2.png", type = 'cairo-png')

