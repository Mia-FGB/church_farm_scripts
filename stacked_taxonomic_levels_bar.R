setwd("/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm")

library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(lubridate)
library(scales)
library(patchwork)
library(viridis)
library(RColorBrewer)

#read in the data, this data has not been filtered but merged with metadata and lineage in python script
#Scripts/Initial_CF_2023_analysis.ipynb from the MARTi output
data <-  read.csv("Regular_Collections_2023/marti_all_data_summed.csv") #note this takes a while to read in

#sort the dates - Nice.Date is the Date collected but easier to read in
data$Date.collected = as.Date(data$Date.collected)

#filter to phylum level and on rel abundance
phylum_data <- data %>%
  filter(NCBI.Rank == "phylum") %>% #2301 entries, many with a 0 read count
  filter(rel_abundance > 0.1) #reduces entries can change the stringency

# Group by Date.collected and Phylum Name, and calculate the mean of rel_abundance
phylum_data_avg <- phylum_data %>%
  group_by(Date.collected, Name) %>%
  summarise(mean_abundance = mean(rel_abundance)) %>%
  mutate(Name = factor(Name, levels = Name[order(-mean_abundance)])) %>% #to order the bars nicer
  ungroup()

#Setting x-axis start and end
break.vec <- c(as.Date("2022-10-01"), seq(from=as.Date("2022-10-01"), to=as.Date("2023-08-31"),by="1 week")) 

#Plotting - spaced by when they were collected
p <- ggplot(
  phylum_data_avg, aes(x = Date.collected, y = mean_abundance, fill = Name)) + 
  geom_bar(stat = "identity", width = 5, color = "black", size = 0.3,) + #Using position = fill makes the bar fill the axis
  ylab("Relative Abundance (> 0.15%)") +
  xlab("Date Collected") +
  ggtitle("Phylum Composition over time") +
  scale_y_continuous(expand = c(0, 0)) +
  #scale_fill_manual(values = brewer.pal(n = 15, name = "Set3")) + #need a big enough palette
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
    axis.title.x = element_blank(), #remove x-axis title
    axis.line = element_line(size = 0.7)) + #thicker axis
  scale_x_date(breaks = break.vec, date_labels = "%d-%b")

# Save as SVG
ggsave(filename = "Regular_Collections_2023/graphs/stacked_phylum_graph.svg", plot = p, width = 8, height = 6)

#------------------
#If I wanted them to be alongside each other, I think this is misleading though --------
# # Convert Date.collected to a factor to treat dates as discrete categories
# phylum_data_avg$Date.collected <- as.factor(phylum_data_avg$Date.collected)
# 
# # Plotting with close spacing
# ggplot(phylum_data_avg, aes(x = Date.collected, y = mean_abundance, fill = Name)) + 
#   geom_bar(stat = "identity", width = 0.8, color = "black") +  # Add thin outline and dodge
#   scale_fill_viridis_d(option = "D") +  # Use viridis color palette
#   theme(axis.title.x = element_blank()) +   # Remove x-axis title
#   ylab("Average Relative Abundance") + 
#   ggtitle("Phylum Composition over Time") + 
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
#------------------

#Stacked bar chart of the pathogenic genera ------------
#read in the pathogenic data
pathogen_data <- read.csv("Regular_Collections_2023/just_pathogen_data.csv", header = TRUE)
pathogen_data$Date.collected = as.Date(pathogen_data$Date.collected)

# Group by Date.collected and Name, and calculate the mean of rel_abundance - they are low levels 
pathogen_data_avg <- pathogen_data %>%
  group_by(Date.collected, Name) %>%
  summarise(mean_abundance = mean(rel_abundance)) %>%
  mutate(Name = factor(Name, levels = sort(unique(Name)))) %>% #order alphabetically
  ungroup()

custom_colors <- c(
  "Zymoseptoria" = "#8DD3C7", 
  "Fusarium" = "#BC80BD", 
  "Pyrenophora" = "#BEBADA", 
  "Parastagonospora" = "#FB8072", 
  "Blumeria" = "#80B1D3", 
  "Phaeosphaeria" = "#FDB462",
  "Puccinia"= "#B3DE69",
  "Ustilago" = "#FCCDE5")

#Plotting pathogen data as a stacked bar
path <- ggplot(
  pathogen_data_avg, aes(x = Date.collected, y = mean_abundance, fill = Name)) + 
  geom_bar(stat = "identity", width = 5, color = "black", size = 0.3) + #Using position = "fill" makes the bar fill the axis
  ylab("Relative Abundance") +
  xlab("Date Collected") +
  ggtitle("Pathogen Composition over time") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
    axis.title.x = element_blank(), #remove x-axis title
    axis.line = element_line(size = 0.7)) + #thicker axis
  scale_x_date(breaks = break.vec, date_labels = "%d-%b") +
  scale_fill_manual(values = custom_colors)  # Use your custom color palette


#ggsave(filename = "Regular_Collections_2023/graphs/stacked_pathogen_graph_proportion.svg", plot = path, width = 8, height = 6)
 
