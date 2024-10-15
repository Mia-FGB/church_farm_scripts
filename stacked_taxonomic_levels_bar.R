setwd("/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm")

library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(lubridate)
library(scales)
library(patchwork)
library(viridis)

#read in the data, this data has not been filtered but merged with metadata and lineage in python script
#Scripts/Initial_CF_2023_analysis.ipynb from the MARTi output
data <-  read.csv("Regular_Collections_2023/marti_all_data_summed.csv") #note this takes a while to read in

#sort the dates - Nice.Date is the Date collected but easier to read in
data$Date.collected = as.Date(data$Date.collected)

#filter to phylum level and on rel abundance
phylum_data <- data %>%
  filter(NCBI.Rank == "phylum") %>% #2301 entries, many with a 0 read count
  filter(rel_abundance > 0.05) #reduces entries can change the stringency

# Group by Date.collected and Phylum Name, and calculate the mean of rel_abundance
phylum_data_avg <- phylum_data %>%
  group_by(Date.collected, Name) %>%
  summarise(mean_abundance = mean(rel_abundance))

#Setting x-axis start and end
break.vec <- c(as.Date("2022-10-01"), seq(from=as.Date("2022-10-01"), to=as.Date("2023-08-31"),by="1 month")) 

#Plotting - spaced by when they were collected
p <- ggplot(
  phylum_data_avg, aes(x = Date.collected, y = mean_abundance, fill = Name)) + 
  geom_bar(stat = "identity", width = 5, color = "black", size = 0.3, position = "fill") + #Using position = fill makes the bar fill the axis
  theme(axis.title.x = element_blank()) +   # Remove x axis title
  scale_fill_viridis_d(option = "D") +  # Use the discrete version of viridis
  ylab("Relative Abundance (> 0.05%)") +
  xlab("Date Collected") +
  ggtitle("Phylum Composition over time") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  scale_x_date(breaks = break.vec, date_labels = "%b") 

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
  summarise(mean_abundance = mean(rel_abundance))

#Plotting pathogen data as a stacked bar
path <- ggplot(
  pathogen_data_avg, aes(x = Date.collected, y = mean_abundance, fill = Name)) + 
  geom_bar(stat = "identity", width = 5, color = "black", size = 0.3, position = "fill") + #Using position = fill makes the bar fill the axis
  theme(axis.title.x = element_blank()) +   # Remove x axis title
  scale_fill_viridis_d(option = "D") +  # Use the discrete version of viridis
  ylab("Relative Abundance (> 0.05%)") +
  xlab("Date Collected") +
  ggtitle("Pathogen Composition over time") +
  theme_bw() + 
  theme(panel.grid.minor = element_blank()) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=0.5)) +
  scale_x_date(breaks = break.vec, date_labels = "%b") 

ggsave(filename = "Regular_Collections_2023/graphs/stacked_pathogen_graph.svg", plot = path, width = 8, height = 6)

