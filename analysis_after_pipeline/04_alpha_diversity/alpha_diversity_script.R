
# [Community Ecology Package](https://cran.r-project.org/web/packages/vegan/index.html)
library(vegan)

# Data Analysis Libraries
library(tidyverse)
library(plyr)
library(agricolae)
library(picante)
library(multcomp)
library(reshape2)
library(gridExtra)

rm(list = ls())

setwd("C:/Users/gabri/Desktop/microbiota_seasonal/Analyses_after_pipeline/04_alpha_diversity")


#import data
aplysina_microbiome <- read_tsv("normalized_aplysina_formatted.txt") %>%
  column_to_rownames(., var = "Species")


t_aplysina = data.frame(t(aplysina_microbiome))

class(t_aplysina)
row.names(t_aplysina)
rowSums(t_aplysina)
rowMeans(t_aplysina)

samples <- read_tsv("Samples.txt")

## Calculating diversity: Shannon, Simpson and Inverted Simpson
shannon <- diversity(t_aplysina, index = "shannon", MARGIN = 1, base = exp(1))
simpson <- diversity(t_aplysina, "simpson", MARGIN = 1)
inverted_simpson <- diversity(t_aplysina, "inv", MARGIN = 1)

## Calculating expected species richness 
rarefied <- rarefy(t_aplysina, 2)

## Calculating fisher.alpha (Genuine count of individuals)
fisher_alpha <- fisher.alpha(t_aplysina)

## Richness and Evenness
observed_richness <- specnumber(t_aplysina, MARGIN = 1) ## Observed richness
pielous_evenness <- shannon/log(observed_richness) #(Pielous Evenness)

## estimating the number of unseen species
##Extrapolated Species Richness in a Species Pool
## CHAO
## Quantitative model
CHAO.P <- t(estimateR(t_aplysina))

## Make a table with the metrics
metrics.P <- data.frame(
  Samples = Samples,
  H_Shannon = shannon,
  H_Simp = simpson,
  D_Inv_Simp = inverted_simpson,
  rarefy = rarefied,
  alpha_fisher = fisher_alpha,
  richness = observed_richness,
  evenness = pielous_evenness,
  CHAO = CHAO.P[, 2],
  ACE = CHAO.P[, 4]) %>%
  remove_rownames()
  
  
write.table(metrics.P, "alpha_diversity.txt", quote=FALSE, row.names = FALSE, sep="\t")


####################################################
#############ANOVA for Shannon (H)##################
####################################################
# We can use analysis of variance (ANOVA) to tell 
# if at least one of the diversity means is different from the rest.

aovShannon.P <- aov(shannon ~ Samples, data=metrics.P)
summary(aovShannon.P)

## Post-hoc Tukey Test for Shannon
# A Tukey's Honest Significant Difference (HSD) test 
# can do pairwise comparisons of the means to find this out.
TukeyHSD(aovShannon.P,ordered=T)

tukey_Shannon.P <- HSD.test(aovShannon.P, "Samples", group = TRUE)
print(tukey_Shannon.P)

#Tukey with multicomp
#tukey_Shannon_multcomp <- summary(glht(aovShannon.P, linfct = mcp(Samples = "Tukey")))
#tukey_Shannon_multcomp.P # Shannon multicomp



## to summarize your data (with mean, standard deviation), broken down by group.
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/


# Run the functions length, mean, and sd on the value of "change" for each group, 
# broken down by sex + condition
summary_Shannon.P <- ddply(metrics.P, c("Samples"), summarise,
                     N    = length(H_Shannon),
                     mean = mean(H_Shannon),
                     sd   = sd(H_Shannon),
                     se   = sd / sqrt(N)
)
summary_Shannon.P

write.table(summary_Shannon.P, file = "HSDTkey_shannon_proka.txt", quote = FALSE, 
            sep = "\t")

group_data_H.P <- tukey_Shannon.P$groups[order(rownames(tukey_Shannon.P$groups)),]


ggplot(metrics.P, aes(x = Samples, y = H_Shannon)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_H.P), y = max(metrics.P$H_Shannon)+0.05 
                , label = group_data_H.P$groups),
            col = 'black',
            size = 8) +
  geom_boxplot() + 
  scale_color_brewer(palette="Dark2")+
  ggtitle("Alpha diversity of Proka Shannon") +
  xlab("") +
  ylab("Shannon diversity index")
  





####################################################
##################ANOVA for Simpson#################
####################################################
# We can use analysis of variance (ANOVA) to tell 
# if at least one of the diversity means is different from the rest.
aovSimpson.P <- aov(H_Simp ~ Samples, data=metrics.P)
summary(aovSimpson.P)
# yes they are

## Post-hoc Tukey Test for Simpson
# A Tukey's Honest Significant Difference (HSD) test 
# can do pairwise comparisons of the means to find this out.
TukeyHSD(aovSimpson.P, ordered=T)

tukey_sim.P <- HSD.test(aovSimpson.P, "Samples", group = TRUE)
print(tukey_sim.P)

#Tukey with multicomp
tukey_sim.P <- summary(glht(aovSimpson.P, linfct = mcp(Samples = "Tukey")))
tukey_sim.P # Simpson multicomp

## Plot a barchart of Simpson with means
group_data_sim.P <- tukey_sim.P$groups[order(rownames(tukey_sim.P$groups)),]

ggplot(metrics.P, aes(x = Samples, y = H_Simp)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_sim.P), y = max(metrics.P$H_Simp)+ 0.001, 
                label = group_data_sim.P$groups),
            col = 'black',
            size = 8) +
  geom_boxplot() +
  ggtitle("Alpha diversity of Proka Simpson") +
  xlab("Samples") +
  ylab("Simpson diversity index")

ggsave("Simpson.svg")

## to summarize your data (with mean, standard deviation), broken down by group.
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/


# Run the functions length, mean, and sd on the value of "change" for each group, 
# broken down by sex + condition
summary_sim.P <- ddply(metrics.P, c("Samples"), summarise,
                       N    = length(H_Simp),
                       mean = mean(H_Simp),
                       sd   = sd(H_Simp),
                       se   = sd / sqrt(N)
)
summary_sim.P

write.table(summary_sim.P, file = "HSDTkey_Simpson_proka.txt", quote = FALSE, 
            sep = "\t")


####################################################
############ANOVA for Inverse_Simpson##############
####################################################
# We can use analysis of variance (ANOVA) to tell 
# if at least one of the diversity means is different from the rest.
aovInvSimp.P <- aov(D_Inv_Simp ~ Samples, data=metrics.P)
summary(aovInvSimp.P)
# yes they are

## Post-hoc Tukey Test for Inverted_Simpson
# A Tukey's Honest Significant Difference (HSD) test 
# can do pairwise comparisons of the means to find this out.
TukeyHSD(aovInvSimp.P,ordered=T)

tukey_D.P <- HSD.test(aovInvSimp.P, "Samples", group = TRUE)
print(tukey_D.P)

#Tukey with multicomp
tukey_D.P <- summary(glht(aovInvSimp.P, linfct = mcp(Samples = "Tukey")))
tukey_D.P # Inverted_Simpson multicomp


## Plot a barchart of Inverted_Simpson with means
group_data_D.P <- tukey_D.P$groups[order(rownames(tukey_D.P$groups)),]


ggplot(metrics.P, aes(x = Samples, y = D_Inv_Simp)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_D.P), y = max(metrics.P$D_Inv_Simp) 
                + 1, label = group_data_D.P$groups),
            col = 'black',
            size = 8) +
  geom_boxplot() +
  ggtitle("Alpha diversity of Proka Inverted_Simpson") +
  xlab("Samples") +
  ylab("Inverted_Simpson diversity index")

ggsave("Inverted_Simpson.svg")

## to summarize your data (with mean, standard deviation), broken down by group.
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/


# Run the functions length, mean, and sd on the value of "change" for each group, 
# broken down by sex + condition
summary_D.P <- ddply(metrics.P, c("Samples"), summarise,
                     N    = length(D_Inv_Simp),
                     mean = mean(D_Inv_Simp),
                     sd   = sd(D_Inv_Simp),
                     se   = sd / sqrt(N)
)
summary_D.P

write.table(summary_D.P, file = "HSDTkey_Inverted_Simpson_proka.txt", quote = FALSE, 
            sep = "\t")


####################################################
#############ANOVA for Expected species#############
####################################################
# We can use analysis of variance (ANOVA) to tell 
# if at least one of the richness means is different from the rest.
aovr.P <- aov(rarefy ~ Samples, data=metrics.P)
summary(aovr.P)
# yes they are

## Post-hoc Tukey Test for Expected species
# A Tukey's Honest Significant Difference (HSD) test 
# can do pairwise comparisons of the means to find this out.
TukeyHSD(aovr.P,ordered=T)

tukey_r.P <- HSD.test(aovr.P, "Samples", group = TRUE)
print(tukey_r.P)

#Tukey with multicomp
tukey_r.P <- summary(glht(aovr.P, linfct = mcp(Samples = "Tukey")))
tukey_r.P # Expected species multicomp


## Plot a barchart of Expected species with means
group_data_r.P <- tukey_r.P$groups[order(rownames(tukey_r.P$groups)),]


ggplot(metrics.P, aes(x = Samples, y = rarefy)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_r.P), y = max(metrics.P$rarefy) 
              , label = group_data_r.P$groups),
            col = 'black',
            size = 8) +
  geom_boxplot() +
  ggtitle("Richness of Proka Expected species") +
  xlab("Samples") +
  ylab("Expected species richness index")

ggsave("Expected_species_richness.svg")


## to summarize your data (with mean, standard deviation), broken down by group.
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/


# Run the functions length, mean, and sd on the value of "change" for each group, 
# broken down by sex + condition
summary_r.P <- ddply(metrics.P, c("Samples"), summarise,
                     N    = length(rarefy),
                     mean = mean(rarefy),
                     sd   = sd(rarefy),
                     se   = sd / sqrt(N)
)
summary_r.P

write.table(summary_r.P, file = "HSDTkey_Expected species_proka.txt", quote = FALSE, 
            sep = "\t")


####################################################
##############ANOVA for Genuine_count###############
####################################################
# We can use analysis of variance (ANOVA) to tell 
# if at least one of the richness means is different from the rest.
aovalpha.P <- aov(alpha_fisher ~ Samples, data=metrics.P)
summary(aovalpha.P)
# yes they are

## Post-hoc Tukey Test for Genuine_count
# A Tukey's Honest Significant Difference (HSD) test 
# can do pairwise comparisons of the means to find this out.
TukeyHSD(aovalpha.P,ordered=T)

tukey_alpha.P <- HSD.test(aovalpha.P, "Samples", group = TRUE)
print(tukey_alpha.P)

#Tukey with multicomp
tukey_alpha.P <- summary(glht(aovalpha.P, linfct = mcp(Samples = "Tukey")))
tukey_alpha.P # Observed OTUs multicomp

## Plot a barchart of Genuine_count with means
group_data_alpha.P <- tukey_alpha.P$groups[order(rownames(tukey_alpha.P$groups)),]


ggplot(metrics.P, aes(x = Samples, y = alpha_fisher)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_alpha.P), y = max(metrics.P$a) 
                + 1, label = group_data_alpha.P$groups),
            col = 'black',
            size = 8) +
  geom_boxplot() +
  ggtitle("Richness of Proka Genuine_count") +
  xlab("Samples") +
  ylab("Genuine_count richness index")

ggsave("Genuine_count_richness.svg")


## to summarize your data (with mean, standard deviation), broken down by group.
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/


# Run the functions length, mean, and sd on the value of "change" for each group, 
# broken down by sex + condition
summary_alpha.P <- ddply(metrics.P, c("Samples"), summarise,
                         N    = length(alpha_fisher),
                         mean = mean(alpha_fisher),
                         sd   = sd(alpha_fisher),
                         se   = sd / sqrt(N)
)
summary_alpha.P

write.table(summary_alpha.P, file = "HSDTkey_Genuine_count_proka.txt", quote = FALSE, 
            sep = "\t")


####################################################
############ANOVA for Observed_OTUs################
####################################################
# We can use analysis of variance (ANOVA) to tell 
# if at least one of the richness means is different from the rest.
aovS.P <- aov(richness ~ Samples, data=metrics.P)
summary(aovS.P)
# yes they are

## Post-hoc Tukey Test for Observed_OTUs
# A Tukey's Honest Significant Difference (HSD) test 
# can do pairwise comparisons of the means to find this out.
TukeyHSD(aovS.P,ordered=T)

tukey_S.P <- HSD.test(aovS.P, "Samples", group = TRUE)
print(tukey_S.P)

#Tukey with multicomp
tukey_S.P <- summary(glht(aovS.P, linfct = mcp(Samples = "Tukey")))
tukey_S.P # Observed OTUs multicomp

## Plot a barchart of Observed_OTUs with means
group_data_S.P <- tukey_S.P$groups[order(rownames(tukey_S.P$groups)),]

ggplot(metrics.P, aes(x = Samples, y = richness)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_S.P), y = max(metrics.P$richness) 
                + 5, label = group_data_S.P$groups),
            col = 'black',
            size = 8) +
  geom_boxplot() +
  ggtitle("Richness of Proka Observed_zOTUs") +
  xlab("Samples") +
  ylab("Observed_zOTUs richness index")

ggsave("Observed_zOTUs.svg")

## to summarize your data (with mean, standard deviation), broken down by group.
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/


# Run the functions length, mean, and sd on the value of "change" for each group, 
# broken down by sex + condition
summary_S.P <- ddply(metrics.P, c("Samples"), summarise,
                     N    = length(richness),
                     mean = mean(richness),
                     sd   = sd(richness),
                     se   = sd / sqrt(N)
)
summary_S.P

write.table(summary_S.P, file = "HSDTkey_Observed_OTUs_proka.txt", quote = FALSE, 
            sep = "\t")


####################################################
#############ANOVA for Pielous_evenness#############
####################################################
# We can use analysis of variance (ANOVA) to tell 
# if at least one of the evenness means is different from the rest.
aovJ.P <- aov(evenness ~ Samples, data=metrics.P)
summary(aovJ.P)
# yes they are

## Post-hoc Tukey Test for Pielous_evenness
# A Tukey's Honest Significant Difference (HSD) test 
# can do pairwise comparisons of the means to find this out.
TukeyHSD(aovJ.P,ordered=T)

tukey_J.P <- HSD.test(aovJ.P, "Samples", group = TRUE)
print(tukey_J.P)

#Tukey with multicomp
tukey_J.P <- summary(glht(aovJ.P, linfct = mcp(Samples = "Tukey")))
tukey_J.P # Pielous_evenness multicomp


## Plot a barchart of Pielous_evenness with means
group_data_J.P <- tukey_J.P$groups[order(rownames(tukey_J.P$groups)),]

ggplot(metrics.P, aes(x = Samples, y = evenness)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_J.P), y = max(metrics.P$evenness) 
                + 0.01, label = group_data_J.P$groups),
            col = 'black',
            size = 8) +
  geom_boxplot() +
  ggtitle("Richness of Proka Pielous_evenness") +
  xlab("Samples") +
  ylab("Pielous_evenness evenness index")

ggsave("Pielous_evenness.svg")

## to summarize your data (with mean, standard deviation), broken down by group.
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/


# Run the functions length, mean, and sd on the value of "change" for each group, 
# broken down by sex + condition
summary_J.P <- ddply(metrics.P, c("Samples"), summarise,
                     N    = length(evenness),
                     mean = mean(evenness),
                     sd   = sd(evenness),
                     se   = sd / sqrt(N)
)
summary_J.P

write.table(summary_J.P, file = "HSDTkey_Pielous_evenness_proka.txt", quote = FALSE, 
            sep = "\t")



####################################################
###################ANOVA for CHAO##################
###################################################
# We can use analysis of variance (ANOVA) to tell 
# if at least one of the estimated richness means is different from the rest.
aovC.P <- aov(CHAO ~ Samples, data=metrics.P)
summary(aovC.P)
# yes they are

## Post-hoc Tukey Test for CHAO
# A Tukey's Honest Significant Difference (HSD) test 
# can do pairwise comparisons of the means to find this out.
TukeyHSD(aovC.P,ordered=T)

tukey_C.P <- HSD.test(aovC.P, "Samples", group = TRUE)
print(tukey_C.P)

#Tukey with multicomp
tukey_C.P <- summary(glht(aovC.P, linfct = mcp(Samples = "Tukey")))
tukey_C.P # CHAO multicomp

## Plot a barchart of CHAO with means
group_data_C.P <- tukey_C.P$groups[order(rownames(tukey_C.P$groups)),]

ggplot(metrics.P, aes(x = Samples, y = CHAO)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_C.P), y = max(metrics.P$CHAO) 
                + 1, label = group_data_C.P$groups),
            col = 'black',
            size = 8) +
  geom_boxplot() +
  ggtitle("Alpha diversity of Proka CHAO") +
  xlab("Samples") +
  ylab("CHAO diversity index")

ggsave("CHAO.svg")

## to summarize your data (with mean, standard deviation), broken down by group.
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/


# Run the functions length, mean, and sd on the value of "change" for each group, 
# broken down by sex + condition
summary_C.P <- ddply(metrics.P, c("Samples"), summarise,
                     N    = length(CHAO),
                     mean = mean(CHAO),
                     sd   = sd(CHAO),
                     se   = sd / sqrt(N)
)
summary_C.P

write.table(summary_C.P, file = "HSDTukey_CHAO_Proka.txt", quote = FALSE, 
            sep = "\t")


####################################################
####################ANOVA for ACE###################
####################################################
# We can use analysis of variance (ANOVA) to tell 
# if at least one of the estimated richness means is different from the rest.
aovA.P <- aov(ACE ~ Samples, data=metrics.P)
summary(aovA.P)
# yes they are

## Post-hoc Tukey Test for ACE
# A Tukey's Honest Significant Difference (HSD) test 
# can do pairwise comparisons of the means to find this out.
TukeyHSD(aovA.P,ordered=T)

tukey_A.P <- HSD.test(aovA.P, "Samples", group = TRUE)
print(tukey_A.P)

#Tukey with multicomp
tukey_A.P <- summary(glht(aovA.P, linfct = mcp(Samples = "Tukey")))
tukey_A.P # ACE multicomp

## Plot a barchart of ACE with means
group_data_A.P <- tukey_A.P$groups[order(rownames(tukey_A.P$groups)),]

ggplot(metrics.P, aes(x = Samples, y = ACE)) +
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_A.P), y = max(metrics.P$ACE) 
                + 1, label = group_data_A.P$groups),
            col = 'black',
            size = 8) +
  geom_boxplot() +
  ggtitle("Alpha diversity of Proka ACE") +
  xlab("Samples") +
  ylab("ACE diversity index")

ggsave("ACE.svg")

## to summarize your data (with mean, standard deviation), broken down by group.
# http://www.cookbook-r.com/Manipulating_data/Summarizing_data/


# Run the functions length, mean, and sd on the value of "change" for each group, 
# broken down by sex + condition
summary_A.P <- ddply(metrics.P, c("Samples"), summarise,
                     N    = length(ACE),
                     mean = mean(ACE),
                     sd   = sd(ACE),
                     se   = sd / sqrt(N)
)
summary_A.P

write.table(summary_A.P, file = "HSDTukey_ACE_Proka.txt", quote = FALSE, 
            sep = "\t")



####################################################
####################SAME GRID GRAPHS################
####################################################


observed_bp <- ggplot(metrics.P, aes(x = Samples, y = richness, color=Samples))+
  geom_boxplot() +  scale_color_brewer(palette="Dark2") + 
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_A.P), y = max(metrics.P$richness)+
            15, label = group_data_A.P$groups),
            col = 'black',
            size = 5) +
  xlab("") +
  ylab("Observed_zOTUs richness index") + theme(legend.position="none")

observed_bp

shannon_bp <- ggplot(metrics.P, aes(x = Samples, y = H_Shannon, color=Samples))+
  geom_boxplot() +  scale_color_brewer(palette="Dark2") + 
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_A.P), y = max(metrics.P$H_Shannon) 
                + 0.05, label = group_data_A.P$groups),
            col = 'black',
            size = 5) +
  xlab("") +
  ylab("Shannon diversity index") + theme(legend.position="none")

shannon_bp


chao_bp <- ggplot(metrics.P, aes(x = Samples, y = CHAO, color=Samples))+
  geom_boxplot() +  scale_color_brewer(palette="Dark2") + 
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_A.P), y = max(metrics.P$CHAO) 
                + 5, label = group_data_A.P$groups),
            col = 'black',
            size = 5) +
  xlab("") +
  ylab("CHAO diversity index") + theme(legend.position="none")

chao_bp


simpson <- ggplot(metrics.P, aes(x = Samples, y = ACE, color=Samples))+
  geom_boxplot() +  scale_color_brewer(palette="Dark2") + 
  geom_text(data = data.frame(),
            aes(x = rownames(group_data_A.P), y = max(metrics.P$ACE) 
                + 5, label = group_data_A.P$groups),
            col = 'black',
            size = 5) +
  ggtitle("Alpha diversity of Proka ACE") +
  xlab("Samples") +
  ylab("ACE diversity index") + theme(legend.position="none")








