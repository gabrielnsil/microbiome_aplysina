#Bioinformatics amplicon analysis short course
#Post-Mothur analysis in R
#The goal of this tutorial is to take you from the Mothur output (.shared file) 
#to meaningful results from your dataset.  We will use the statistical program R 
#to analyze and visualize our data.  

#1.	Open R and check to see the location of your working file directory.  
#All files that we output in this exercise will be located in this directory.
#getwd()
#Should you want to change your working directory enter  
#setwd(enter directory path here) in the R console.

setwd("C:/Users/gabri/Desktop/microbiota_seasonal/Analyses_after_pipeline/02_nMDS")
rm(list = ls())

#install.packages("vegan")
library(vegan)
library(picante)
library(tidyverse)
library(janitor)
library(Cairo)


aplysina <- read_tsv('01_zotus_table_for_nmds_formatted.txt')

aplysina1 <- data.frame(t(aplysina)) %>%
  row_to_names(row_number = 1) %>%
  mutate(across(everything(), as.numeric))

sapply(aplysina1, class)


#4.	Confirm your file is correct
class(aplysina1)     # checks the file type (should be âdata.frame)
dim(aplysina1)   # gives you the dimensions of your table
rownames(aplysina1) # prints the row names
head(colnames(aplysina1))  # prints the first 5 column names

#5.	Check for the total number of reads in eaplysina1h sample
apply(aplysina1, 1,sum)
#Are the row sums equal?  If not we need to transform the data.

#6.	Transform the data to relative abundance. Turn raw number of reads into relative 
#abundance by dividing each value by sample total abundance.
aplysina2 <- decostand(aplysina1, method = "total")

#7.	Confirm transformation (using proka3 and in step 5).  Total abundance of each sample should now equal 1.
apply(aplysina2, 1, sum)

#8.	Import your metadata file into the R console
metadata2 <- read.table("02_design_seasons.txt", header = TRUE, row.names = 1)

#9.	Take a look at the metadata file
head(metadata2)
tail(metadata2)
class(metadata2)     # checks the file type (should be data.frame)
dim(metadata2)   # gives you the dimensions of your table
rownames(metadata2) # prints the row names
head(colnames(metadata2))  # prints the first 5 column names

#10.	Sort the rows in âproka3 to match the row order of âmetadata2
aplysina3 <- aplysina2[rownames(metadata2), ]

#11.	Check to make sure the row names for eaplysina1h file match.
all.equal(rownames(aplysina3), rownames(metadata2))

#12.	Multivariate analysis: Calculate the Bray-Curtis distance among samples
aplysina2.bc.dist <- vegdist(aplysina2, method = "bray")

#13.	Cluster the communities using average-linkage algorithm
aplysina3.bc.clust <- hclust(aplysina2.bc.dist, method = "average")

#14.	Visualize community dissimilarity with a cluster dendrogram
plot(aplysina3.bc.clust, cex = 0.5, ylab = "Bray-Curtis dissimilarity")

#15.	Visualize the Bray-Curtis distance with an MDS plot
aplysina3.bc.mds <- metaMDS(aplysina3, dist = "bray")
ordiplot(aplysina3.bc.mds, display = "sites", type = "text")

#16.	Customize your  MDS visualization


stress <- aplysina3.bc.mds$stress
number_asv = nrow(aplysina)
stressplot(aplysina3.bc.mds)


mds.fig <- ordiplot(aplysina3.bc.mds, type = "none")
points(mds.fig, "sites", pch = 19, col = "magenta", select = metadata2$Season == "Autumn")
points(mds.fig, "sites", pch = 19, col = "green4", select = metadata2$Season == "Winter")
points(mds.fig, "sites", pch = 19, col = "tan1", select = metadata2$Season == "Spring")
points(mds.fig, "sites", pch = 19, col = "mediumblue", select = metadata2$Season == "Summer")

ordiellipse(aplysina3.bc.mds, metadata2$Season, conf = 0.95, label = TRUE)

legend("bottomleft", 
       legend = c(paste("2D Stress = ", round(stress, 5)), paste("Number of zOTUs = ", number_asv)),
       cex = 1, bty = "n")

       
Cairo(file="Cairo_PNG_72_dpi.png", 
      type="png",
      units="in", 
      width=5, 
      height=4, 
      pointsize=12, 
      dpi=300)
dev.off()

#ordicluster(proka2.bc.mds, proka2.bc.clust, col = "gray")

#17.	Test for differences in Bray-Curtis dissimilarity among host species
adonis1 <- adonis(aplysina2.bc.dist ~ Season, data = metadata2, perm = 9999)


#18. Anosim
anosim1 <- anosim(aplysina2.bc.dist, metadata2$Season, perm = 9999)

print(adonis1) 
print(anosim1)





#A list of class "simper" with following items:
#species=The species names.
#average=Average contribution to overall dissimilarity.
#overall=The overall between-group dissimilarity.
#sd=Standard deviation of contribution.
#ratio=Average to sd ratio.
#ava, avb=Average abundances per group.
#ord=An index vector to order vectors by their contribution or order cusum baplysina1k to the original data order.
#cusum=Ordered cumulative contribution.
#p=Permutation p-value. Probability of getting a larger or equal average contribution in random permutation of the group faplysina1tor.




##############
#https://rdrr.io/rforge/vegan/man/adonis.html
#adonis/PERMANOVA also with the crossed design of habitat and site 
adonis(aplysina2.bc.dist ~ Habitat * Site, data = metadata2, perm = 10000)




########################################
library(vegan)
data(varespec)
vare_dist <- vegdist(varespec, method="bray")
vare_mds <- metaMDS(comm = vare_dist, autotransform = FALSE) 
# aplysina1tually autotransform = FALSE doesn't seem to change the results
plot(vare_mds, type = "t")

vare_mds_2 <- metaMDS(comm = varespec, distance = "bray", k =2)
plot(vare_mds_2, display = "sites", type = "t")

# plots above are different and the stress values below as well
vare_mds$stress; vare_mds_2$stress
# [1] 0.1000211
# [1] 0.1843196













class(aplysina1)     # checks the file type (should be a data.frame)
dim(aplysina1)   # gives you the dimensions of your table
rownames(aplysina1) # prints the row names
head(colnames(aplysina1))  # prints the first 5 column names



set.seed(2)
community_matrix=matrix(
  sample(1:100,300,replace=T),nrow=10,
  dimnames=list(paste("community",1:10,sep=""),paste("sp",1:30,sep="")))

example_NMDS=metaMDS(community_matrix, # Our community-by-species matrix
                     k=2) # The number of reduced dimensions

example_NMDS=metaMDS(community_matrix,k=2,trymax=100)

stressplot(example_NMDS)

plot(example_NMDS)
ordiplot(example_NMDS,type="n")
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",cex=1,air=0.01)


treat=c(rep("Treatment1",5),rep("Treatment2",5))
ordiplot(example_NMDS,type="n")
ordihull(example_NMDS,groups=treat,draw="polygon",col="grey90",label=F)
orditorp(example_NMDS,display="species",col="red",air=0.01)
orditorp(example_NMDS,display="sites",col=c(rep("green",5),rep("blue",5)),
         air=0.01,cex=1.25)
