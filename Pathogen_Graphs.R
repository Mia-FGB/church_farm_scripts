# Script to create pathogen plots from MARTi read count data

#Packages ------
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lubridate) # for date time
library(reshape2)  #to melt data
library(RColorBrewer)

# Set plotting theme
theme_set(theme_bw())

#Reading in & Procesisng data ----------

# Read counts ---
# Using summed data so I can work at Genus level 
data <-  read.delim("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/all_data_mar25_summed.tsv", check.names = FALSE)


pathogens <- c("Puccinia", "Blumeria", "Fusarium", "Zymoseptoria", "Ustilago", "Magnaporthe",
               "Claviceps", "Pyrenophora", "Parastagonospora", "Phaeosphaeria")

pathogen_data <- data %>% 
  filter(`NCBI Rank` == 'genus') %>% 
  filter(Name %in% pathogens)

# Metadata ---
meta <- read.delim("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/metadata/metadata_cleaned.tsv", check.names=FALSE)

#convert to date time
meta$DateTime_UTC <- ymd_hms(meta$DateTime_UTC)

#  Read numbers ---
read_no <- read.delim("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/metadata/read_numbers.tsv")
read_no <- read_no  %>%
  rename(Sample_ID = ID)

# Melting & merging the data ------

# Melt data to long 
pathogen_data <- melt(pathogen_data, id.vars = c("Name", "NCBI ID", "NCBI Rank"), 
                      variable.name = "Sample_ID", value.name = "Read_Count")

# Merge on metadata & read numbers
pathogen_meta <- left_join(pathogen_data, meta, by = "Sample_ID")
pathogen_meta <- left_join(pathogen_meta, read_no, by = "Sample_ID")

# Calculate HPM
pathogen_meta <- pathogen_meta %>%
  mutate(HPM = (Read_Count / ReadsPassBasecall) * 1000000)

# Removing low read count values, but keeping as 0 for average calcualtion 
pathogen_meta <- pathogen_meta %>% 
  mutate(HPM = ifelse(HPM < 100, 0, HPM)) %>% 
  filter(!Sample %in% c("Positive", "Negative"))


# Plotting --------
plot_pathogen_bar <- function(years, name){
  filtered_data <- pathogen_meta %>% filter(Year %in% years, Name == name)
  
  # Check if there is any data with HPM > 0
  if (sum(filtered_data$HPM > 0, na.rm = TRUE) == 0) {
    message("No data with HPM > 0 for ", name, " in " , years, ". Exiting function.")
    return()  # Exit the function if there is no HPM > 0 data
  }
  
  # Compute mean HPM for each Sample_Time
  summary_data <- filtered_data %>%
    group_by(Sample, DateTime_UTC) %>%
    summarise(
      Avg_HPM = mean(HPM),
      SD_HPM = sd(HPM),  # Standard deviation
      SE_HPM = SD_HPM / sqrt(n()),  # Standard error
      .groups = "drop"
    ) %>% arrange(DateTime_UTC)
  
  # Ensure Sample is uniquely ordered by Sample
  summary_data$Sample <- factor(summary_data$Sample, levels = unique(summary_data$Sample))
  
  pathogen_plot <- ggplot(summary_data, aes(x = Sample, y = Avg_HPM)) +
    geom_bar(stat = "identity", position = "identity") +
    geom_errorbar(aes(ymin = Avg_HPM - SE_HPM, ymax = Avg_HPM + SE_HPM), color = "black") +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", linewidth = 0.8),  # Make axis lines thicker and black
      panel.border = element_blank()  # Remove the panel border
    ) +
    scale_y_continuous(expand = c(0, 0)) +  # Set the y-axis to start at 0
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
    labs(title = paste(name, paste(years, collapse = ", ")),
         x = "Date Collected",
         y = "Average HPM")
  
  filename <- paste0("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/pathogen_graphs/", 
                     name, "_", paste(years, collapse = "_"), ".svg")
  ggsave(filename, plot = pathogen_plot, device = "svg", width = 10, height = 6)
}

#Example for one plot:
plot_pathogen_bar(2024, "Zymoseptoria")

# Loop through each pathogen name 2024
for (name in pathogens) {
  # Call the function with the current month, year, and pathogen name
  plot_pathogen_bar(years = 2024, 
                         name = name)
}

# Loop through each pathogen name 2022 - 2023
for (name in pathogens) {
  # Call the function with the current month, year, and pathogen name
  plot_pathogen_bar(years = c(2022, 2023), 
                    name = name)
}


# Multiple years on the same axis --------
# Is working but in order for the bars to have a nice width, have Sample on x-axis which is a bit messy
# Tried using the date-formtatitng approach but it looked worse

pathogen_meta$Date_No_Yr <- format(pathogen_meta$DateTime_UTC, "1900-%m-%d")
# Convert Date_No_Yr to a proper date format (assuming it's in MM-DD format)
pathogen_meta$Date_No_Yr <- as.Date(pathogen_meta$Date_No_Yr, format = "%Y-%m-%d")

plot_pathogen_bar_yrs <- function(name){
  filtered_data <- pathogen_meta %>% filter(Name == name)
  
  # Compute mean HPM for each Sample_Time
  summary_data <- filtered_data %>%
    group_by(Sample, Date_No_Yr, Year) %>%
    summarise(
      Avg_HPM = mean(HPM),
      SD_HPM = sd(HPM),  # Standard deviation
      SE_HPM = SD_HPM / sqrt(n()),  # Standard error
      .groups = "drop"
    ) %>% arrange(Date_No_Yr)
  
  # Ensure Sample is uniquely ordered by Sample
  summary_data$Sample <- factor(summary_data$Sample, levels = unique(summary_data$Sample))
  
  pathogen_plot <- ggplot(summary_data, aes(x = Sample, y = Avg_HPM, fill = as.factor(Year))) +
    geom_bar(stat = "identity", position = "identity", alpha = 0.7) +
    geom_errorbar(aes(ymin = Avg_HPM - SE_HPM, ymax = Avg_HPM + SE_HPM), color = "black") +
    theme(
      panel.grid.major = element_blank(),  # Remove major grid lines
      panel.grid.minor = element_blank(),  # Remove minor grid lines
      axis.line = element_line(color = "black", linewidth = 0.8),  # Make axis lines thicker and black
      panel.border = element_blank()  # Remove the panel border
    ) +
    scale_y_continuous(expand = c(0, 0)) +  # Set the y-axis to start at 0
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5)) +
    labs(title = paste(name),
         x = "Date Collected",
         y = "Average HPM",
         fill = "Year")
  
  # filename <- paste0("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/pathogen_graphs/", 
  #                   name, "_bar.svg")
  # ggsave(filename, plot = pathogen_plot, device = "svg", width = 10, height = 6)
  print(pathogen_plot)
}

plot_pathogen_bar_yrs("Zymoseptoria")

# Loop through each pathogen name
for (name in pathogens) {
  # Call the function with the current month, year, and pathogen name
  plot_pathogen_bar_yrs(name)
}

print(pathogens)


