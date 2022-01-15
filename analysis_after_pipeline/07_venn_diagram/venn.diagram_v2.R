

library(tidyverse)
library(janitor)
library(VennDiagram)
library(janitor)

# Resolve problema das fontes para Windows

library(extrafont)
library(remotes)
remotes::install_version("Rttf2pt1", version = "1.3.8")
font_import()




setwd("C:/Users/gabri/Desktop/microbiota_seasonal/Analyses_after_pipeline/07_venn_diagram")


###Load the table with the OTU_ID per sample
aplysina_notNorm <- read_tsv("finalzOTUsTable_formatted.txt")


inputVenn <- aplysina_notNorm %>%
  rowwise() %>%
  mutate(Autumn = mean(c_across(starts_with('Autumn')))) %>%
  mutate(Winter = mean(c_across(starts_with('Winter')))) %>%
  mutate(Spring = mean(c_across(starts_with('Spring')))) %>%
  mutate(Summer = mean(c_across(starts_with('Summer')))) %>%
  select(Autumn, Winter, Spring, Summer) %>%
  mutate_if(is.numeric, round)


#### Select the data from the table, in this case, diagram with only sponges
autumn <- inputVenn$Autumn
winter <- inputVenn$Winter
spring <- inputVenn$Spring
summer <- inputVenn$Summer


#### Samples' name outside the graph
venn.plot_pdf <- venn.diagram(
  x = list(
    autumn = autumn,
    winter = winter,
    spring = spring,
    summer = summer
  ), 
  filename = NULL,
  col = "black",
  fill = c("magenta", "green", "red", "blue"),
  na = "remove",
  cex = 1.5,
  cat.col = c("black", "black", "black", "black"),
  cat.cex = 1.5,
)
pdf(file="venn_pdf_seasons.pdf")
grid.draw(venn.plot_pdf)
dev.off()

#### Samples' name outside the graph and file in svg format
venn.plot_svg <- venn.diagram(
  x = list(
    autumn = autumn,
    winter = winter,
    spring = spring,
    summer = summer
  ), 
  filename = NULL,
  col = "black",
  fill = c("magenta", "green", "red", "blue"),
  na = "remove",
  cex = 1.5,
  fontfamily = "Arial",
  cat.col = c("black", "black", "black", "black"),
  cat.cex = 1.5,
)

svg("venn_svg_seasons.svg")
grid.draw(venn.plot_svg)
dev.off()


#### Samples's name outside the graph and file in png format

venn.plot_png <- venn.diagram(
  x = list(
    autumn = autumn,
    winter = winter,
    spring = spring,
    summer = summer
  ), 
  filename = NULL,
  col = "black",
  fill = c("magenta", "green", "red", "blue"),
  na = "remove",
  cex = 1.5,
  cat.col = c("black", "black", "black", "black"),
  cat.cex = 1.5,
)

png("venn_png_seasons.png", type = "cairo")
grid.draw(venn.plot_png)
dev.off()

#From Venn_diagram manual
# Save picture to non-TIFF file type
# currently working on adding this functionality directly into venn.diagram
#venn.plot <- venn.diagram(
 # venn.diagram
  #x = list (
    #A = 1:10,
    #B = 6:25)
  #filename = NULL
#jpeg("venn_jpeg.jpg");
#grid.draw(venn.plot);
#dev.off();
## End(Not run)
