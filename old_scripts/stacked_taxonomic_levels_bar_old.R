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


# Graphs without unclassified proportion ----------------------------------


#sort the dates - Nice.Date is the Date collected but easier to read in
data$Date.collected = as.Date(data$Date.collected)

#filter to phylum level and on rel abundance
#want to retain those with Unassigned Phyla, but current filtering isn't working
#since these are summed reads can't just take the unassigned numbers as they are everything
#could use NumContigs from contig stats?

phylum_data <- data %>%
  filter(NCBI.Rank == "phylum") %>% #2301 entries, many with a 0 read count
  filter(rel_abundance > 0.01) #reduces entries can change the stringency

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
#ggsave(filename = "Regular_Collections_2023/graphs/stacked_phylum_graph.svg", plot = p, width = 8, height = 6)

#------------------

#If I wanted them to be alongside each other, I think this is misleading though 

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

#Pathogen only stacked bar chart, without unclassifieds  ----------------

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



# Including Unclassified -------------------------------------------------
#Phylum level first ----

#Create filtered datasets with the same columns so I can join them altogether for plotting
#note if I change the filtering I need to remake the phylum_data and classified_reads datasets 
#Normal phylum data, filtered to be rel_abundance > 0.1
phylum_data_uncl <- data %>%
  filter(NCBI.Rank == "phylum") %>% #2301 entries, many with a 0 read count
  filter(rel_abundance > 0.01) %>%  #reduces entries can change the stringency
  select(Sample, Date.collected, read.count, NumContigs, phylum)

# Group and aggregate the candidatus data which has two species classified as phylum?
phylum_data_uncl <- phylum_data_uncl %>%
  # Modify 'phylum' to group the selected rows into "Candidatus"
  mutate(phylum = ifelse(phylum %in% c("Candidatus Binatota", "Candidatus Krumholzibacteriota"), 
                         "Candidatus", phylum)) %>%
  # Group by Sample
  group_by(Sample, phylum, Date.collected, NumContigs) %>%
  # Sum read counts within each group
  summarise(read.count = sum(read.count), .groups = "drop")  # Keep other columns as needed


# Ensure all phyla are present in each sample, filling missing values with 0
# This is how I am creating the averages whilst only plotting those with rel_abundance >0.01
phylum_data_uncl <- phylum_data_uncl %>%
  complete(Sample, phylum, fill = list(read.count = 0))

#Fill missing Date.collected and NumContigs within each Sample, these were created using complete
phylum_data_uncl <- phylum_data_uncl %>%
  group_by(Sample) %>%
  fill(Date.collected, NumContigs, .direction = "downup") %>%
  ungroup()

# "classified" reads (sum of reads that did not pass the filter)
classified_data <- data %>%
  filter(NCBI.Rank == "phylum", rel_abundance <= 0.01) %>%  # Entries that did not pass the filter
  group_by(Sample, Date.collected, NumContigs) %>%  # Group by sample
  summarise(read.count = sum(read.count)) %>%  # Sum the reads not passing the filter
  mutate(phylum = "Classified") 

# unclassified reads - calculated by deducting the total of the phylum level classified reads from NumContigs
unclassified_data <- data %>% 
  group_by(Sample, Date.collected, NumContigs) %>% 
  filter(NCBI.Rank == "phylum") %>% 
  summarise(total_class_reads = sum(read.count)) %>% #sum up the total phylum level classified reads
  mutate(read.count = NumContigs - total_class_reads, phylum = "Unclassified") %>% 
  select(Sample, Date.collected, read.count, NumContigs, phylum)

# concatenate the datasets together
full_data <- bind_rows(phylum_data_uncl, classified_data, unclassified_data)

# add a column for relative abundance 
full_data <-  full_data %>% 
  mutate(rel_abundance = (read.count/NumContigs)*100, HPM = (read.count/NumContigs)*1000000)

#calculate average from 2 replicates on the same day
full_data_avg <-  full_data %>% 
  group_by(Date.collected, phylum) %>%
  summarise(mean_abundance = mean(rel_abundance)) %>% 
  mutate(phylum = factor(phylum, levels = phylum[order(-mean_abundance)])) %>% #to order the bars nicer
  ungroup()

full_data_avg$Date.collected = as.Date(full_data_avg$Date.collected)

# Set custom order for phylum factor
full_data_avg <- full_data_avg %>%
  mutate(phylum = factor(phylum, 
                         levels = c(sort(setdiff(unique(phylum), c("Classified", "Unclassified"))), 
                                    "Classified", "Unclassified")))
#To get the colour palette --

# Identify unique phyla and exclude "classified" and "unclassified"
phyla <- sort(setdiff(unique(full_data_avg$phylum), c("Classified", "Unclassified")))
# Define color palette for the remaining phyla
n_phyla <- length(phyla)
main_palette <- colorRampPalette(brewer.pal(12, "Set3"))(n_phyla)
# Combine main palette with black and grey for "classified" and "unclassified"
custom_colors <- c(setNames(main_palette, phyla), 
                   "Classified" = "#D3D3D3", 
                   "Unclassified" = "#4F4F4F")

#Plotting phylum level data with classified & unclassified ---
up <- ggplot(
  full_data_avg, aes(x = Date.collected, y = mean_abundance, fill = phylum)) + 
  geom_bar(stat = "identity", color = "black", size = 0.3,) + #Using position = fill makes the bar fill the axis
  scale_fill_manual(values = custom_colors) + 
  ylab("Relative Abundance (> 0.01%)") +
  xlab("Date Collected") +
  ggtitle("Phylum Composition over time") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
    axis.title.x = element_blank(), #remove x-axis title
    axis.line = element_line(size = 0.7)) + #thicker axis
  scale_x_date(breaks = break.vec, date_labels = "%d-%b")

#ggsave(filename = "Regular_Collections_2023/graphs/stacked_phylum_graph_unclassified.svg", plot = up, width = 15, height = 6)

#Monthly phylum plot ---------
#also selecting samples from the weekly data that are closest to a month on
#will all have to change once I get the actual dates
month_dates <- as.Date(
  c("2022-10-03", "2022-11-14", "2022-12-14", "	2023-01-14",
    "2023-02-14", "2023-03-14", "2023-04-14", "2023-05-25",
    "2023-06-14", "2023-07-12", "2023-08-17"))

# Filter to include only rows with dates in `month_dates`
month_data <- full_data_avg %>%
  filter(Date.collected %in% month_dates)  %>%
  mutate(Date.collected = factor(Date.collected, levels = month_dates))

#Plot
month <- ggplot(
  month_data, aes(x = Date.collected, y = mean_abundance, fill = phylum)) + 
  geom_bar(stat = "identity", color = "black", size = 0.3,) + #Using position = fill makes the bar fill the axis
  scale_fill_manual(values = custom_colors) + 
  ylab("Relative Abundance (> 0.01%)") +
  xlab("Date Collected") +
  ggtitle("Monthly Phylum Composition") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
    axis.title.x = element_blank(), #remove x-axis title
    axis.line = element_line(size = 0.7)) #thicker axis

#ggsave(filename = "Regular_Collections_2023/graphs/month_stacked_phylum_graph_unclassified.svg", plot = month, width = 8, height = 6)

#Weekly Phylum-------
week_dates <- as.Date(
  c("2023-05-25", "2023-06-02", "2023-06-06", "2023-06-14",
    "2023-06-21","2023-06-30", "2023-07-07", "2023-07-12",
    "2023-07-20", "2023-07-27", "2023-08-03", "2023-08-09",
    "2023-08-17"))

# Filter to include only rows with dates in `month_dates`
week_data <- full_data_avg %>%
  filter(Date.collected %in% week_dates)  %>%
  mutate(Date.collected = factor(Date.collected, levels = week_dates))


#Plot
week <- ggplot(
  week_data, aes(x = Date.collected, y = mean_abundance, fill = phylum)) + 
  geom_bar(stat = "identity", color = "black", size = 0.3,) + #Using position = fill makes the bar fill the axis
  scale_fill_manual(values = custom_colors) + 
  ylab("Relative Abundance (> 0.01%)") +
  xlab("Date Collected") +
  ggtitle("Weekly Phylum Composition") +
  scale_y_continuous(expand = c(0, 0)) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
    axis.title.x = element_blank(), #remove x-axis title
    axis.line = element_line(size = 0.7)) #thicker axis

#ggsave(
#  filename = "Regular_Collections_2023/graphs/week_stacked_phylum_graph_unclassified.svg",
#       plot = week, width = 8, height = 6)
