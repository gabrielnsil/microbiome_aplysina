#install.packages("mvabund")
library(mvabund)
library(tidyverse)


####################################################
#################ALL SAMPLES########################
####################################################



setwd("C:/Users/gabri/Desktop/microbiota_seasonal/Analyses_after_pipeline/01_mvabund")

#import the data (reminder: normalized data)

aplysina_microbiome <- read_tsv("normalized_aplysina_formatted.txt") %>%
  column_to_rownames(., var = "Species")

t_aplysina = data.frame(t(aplysina_microbiome))
samples <- read_tsv("samples.txt") # rows with the name of the seasons that the sample was collected

class(t_aplysina)

aplysina <- cbind.data.frame(samples, t_aplysina)

# rename the column "Samples" to "season" and remove rownames
aplysina <- rename(aplysina, "season" = "Samples")

row.names(aplysina) <- NULL


#check if the data are correct

aplysina %>% 
  select(1:3) %>% 
  head(., 10)




#convert it to an mvabund object

seasons_aplysina <- mvabund(aplysina[,2:ncol(aplysina)])


#Having a quick look at the spread of the data using the boxplot function.
par(mar=c(1,1,1,1)) # adjusts the margins par(mar = c(bottom, left, top, right))
boxplot(seasons_aplysina[,1:10], horizontal = TRUE,las=2, main="Abundance")

#Check the mean-variance relationship among samples
par(mar = c(5,5,3,1))
meanvar.plot(seasons_aplysina)

#To contrast transformed abundances to the predictor variables selected,
# To work, we must transform the predictor variable into a factor

season_factor <- as.factor(aplysina$season)

plot(seasons_aplysina ~ season_factor, cex.axis=0.8, cex=1)


##The model syntax that will fit the response variable to the predictor variable 

mod1 <- manyglm(seasons_aplysina ~ aplysina$season, family="poisson")
mod2 <- manyglm(seasons_aplysina ~ aplysina$season, family="negative_binomial")

#Check the model assumptions.
plot(mod1)
plot(mod2)


# If a random scatter of points occurs, the model's assumption is correct. 
# In the case of sponge, it did not happen using "poisson".
# Use the family argument to choose a distribution which is better suited to the data. 
# For count data which does not fit the 'poisson distribution, we can use the negative_binomial distribution.


####Interpreting the results
#Test the multivariate hypothesis of whether species composition varied across the seasons by using the anova function. 
anova_mod2 <- anova(mod2)
capture.output(print(anova_mod2), file="results_anova.txt")

# The results show that the species composition of bacteria and archea of the
# microbiome of Aplysina fulva varies significantly between seasons.



#Run univariate tests for each species separately.This is done by using the p.uni="adjusted" argument in the 
#anova function. 
#mod2_an <-anova(mod2, p.uni="adjusted")
#write.table(mod2_an, "mod2_ANOVA.txt", sep = "\t", quote = F)


