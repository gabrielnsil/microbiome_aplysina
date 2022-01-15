####################################################
#################ALL SAMPLES########################
####################################################
rm(list = ls())
library(vegan)
library(tidyverse)
library(janitor)
library(Cairo)

setwd("C:/Users/gabri/Desktop/microbiota_seasonal/Analyses_after_pipeline/03_rarefaction_curve")

aplysina_notNorm <- read_tsv("finalzOTUsTable_formatted.txt")
aplysina_norm <- read_tsv("normalized_aplysina_formatted.txt")


test1 <- adorn_totals(aplysina_notNorm, where = c("row", "col"))
test2 <- adorn_totals(aplysina_norm, where = c("row", "col"))

inputVegan <- aplysina_notNorm %>%
  rowwise() %>%
  mutate(Autumn = mean(c_across(starts_with('Autumn')))) %>%
  mutate(Winter = mean(c_across(starts_with('Winter')))) %>%
  mutate(Spring = mean(c_across(starts_with('Spring')))) %>%
  mutate(Summer = mean(c_across(starts_with('Summer')))) %>%
  select(Autumn, Winter, Spring, Summer) %>%
  mutate_if(is.numeric, round)
 

row.names(inputVegan) <- aplysina_notNorm$Species
inputVegan <- as.data.frame(t(inputVegan))
str(inputVegan)

diversity(inputVegan,index = "simpson")
simpson <- diversity(inputVegan, "simpson") # or assign to var.
simpson 

shannon <- diversity(inputVegan) # note that Shannon's is default
shannon #Typically ranges from 1.5 - 3.4, higher = more diverse 


# lets compare the two
par(mfrow = c(1, 2))  # use par to generate panels with 1 row of 2 graphs
hist(simpson)
hist(shannon)

######Rarefaction Species Richness, rarefy(vegan)


inputVegan
S <- specnumber(inputVegan) # observed number of species
S

raremax <- min(rowSums(inputVegan)) # is the minimum sample count achieved over the 4 samples (seasons)
raremax

Srare <- rarefy(inputVegan, raremax)
Srare

par(mfrow = c(1,1))
col <- c("black", "darkred", "darkblue", "hotpink")
lty <- c("solid", "dotdash")
lwd <- 1.5
pars <- expand.grid(col = col, lty = lty, lwd = lwd, stringsAsFactors = FALSE)
head(pars)

plot(S, Srare, xlab = "Observed No. of Species", ylab = "Rarefied No. of Species")
abline(0, 1)

rarecurve <- with(pars[1:4, ], rarecurve(inputVegan, step = 20, sample = raremax, 
                                         col = col, lty = lty, lwd = lwd, cex = 0.6))

Cairo(file="rarefactionCurve_seasons.png", 
      type="png",
      dpi=72)
rarecurve <- with(pars[1:4, ], rarecurve(inputVegan, step = 20, sample = raremax, 
                                         col = col, lty = lty, lwd = lwd, cex = 0.6))
dev.off()

