# Script to create pathogen plots from MARTi read count data

#Packages ------
library(ggplot2)
library(dplyr)
library(tidyverse)
library(lubridate) # for date time
library(reshape2)  #to melt data
library(RColorBrewer)
library(patchwork)

# Set plotting theme
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

#Reading in & Procesisng data ----------

# Read counts ---
# Using summed data so I can work at Genus level 
data <-  read.delim("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/all_data_mar25_summed.tsv", check.names = FALSE)


pathogens <- c("Puccinia", "Blumeria", "Fusarium", "Zymoseptoria", "Ustilago", "Magnaporthe",
               "Claviceps", "Pyrenophora", "Phaeosphaeria", "Parastagonospora", )

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
                     name, "_", paste(years, collapse = "_"), ".pdf")
  ggsave(filename, plot = pathogen_plot, device = "svg", width = 10, height = 6)
  
  return(pathogen_plot)
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

# Single Year panel plots  ---------
plot_pathogen_panel <- function(data, years, panel_label) {
  filtered_data <- data %>%
    filter(Year %in% years, Name %in% pathogens) %>%
    group_by(Name, Sample, DateTime_UTC) %>%
    summarise(
      Avg_HPM = mean(HPM),
      SD_HPM = sd(HPM),
      SE_HPM = SD_HPM / sqrt(n()),
      .groups = "drop"
    ) %>%
    arrange(DateTime_UTC)
  
  # Ensure samples are chronologically ordered
  filtered_data$Sample <- factor(filtered_data$Sample, levels = unique(filtered_data$Sample))
  
  p <- ggplot(filtered_data, aes(x = Sample, y = Avg_HPM)) +
    geom_bar(stat = "identity", fill = "#BD98A2") +
    geom_errorbar(aes(ymin = Avg_HPM - SE_HPM, ymax = Avg_HPM + SE_HPM), width = 0.2) +
    facet_wrap(~ Name, scales = "free_y", ncol = 2) +  # Adjust ncol as needed
    custom_theme +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      strip.text = element_text(size = 10, face = "bold"),
      panel.grid = element_blank()
    ) +
    labs(title = paste("Pathogen HPM Panel -", panel_label),
         x = "Date",
         y = "Average HPM")

  ggsave(
    paste0("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/pathogen_graphs/Panel_", panel_label, ".pdf"),
    plot = p, width = 12, height = 10
  )
  
  return(p)
}

panel_2024 <- plot_pathogen_panel(pathogen_meta, years = 2024, panel_label = "2024")
panel_2022_2023 <- plot_pathogen_panel(pathogen_meta, years = c(2022, 2023), panel_label = "2022_2023")

# Panel Plot ------------

## Create a function to return a dual-year plot per pathogen
plot_pathogen_patch <- function(name, show_strip = FALSE) {
  filtered_data <- pathogen_meta %>%
    filter(Name == name, Year %in% c(2022, 2023, 2024)) %>%
    mutate(Year_Group = ifelse(Year == 2024, "2024", "2022–2023"))
  
  if (sum(filtered_data$HPM > 0, na.rm = TRUE) == 0) {
    message("Skipping: no data for ", name)
    return(NULL)
  }
  
  summary_data <- filtered_data %>%
    group_by(Year_Group, Sample, DateTime_UTC) %>%
    summarise(
      Avg_HPM = mean(HPM),
      SD_HPM = sd(HPM),
      SE_HPM = SD_HPM / sqrt(n()),
      .groups = "drop"
    ) %>%
    arrange(Year_Group, DateTime_UTC)
  
  summary_data$Sample <- factor(summary_data$Sample, levels = unique(summary_data$Sample))
  
  base_plot <- ggplot(summary_data, aes(x = Sample, y = Avg_HPM)) +
    geom_bar(stat = "identity", fill = "#B2909F") +
    geom_errorbar(aes(ymin = Avg_HPM - SE_HPM, ymax = Avg_HPM + SE_HPM), width = 0.2) +
    facet_wrap(~ Year_Group, scales = "free_x", ncol = 2) +
    labs(title = name, x = "Date Collected", y = "Average HPM") +
    custom_theme +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 6),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(size = 10, face = "bold"),
      panel.grid = element_blank()
    )
  
  # Optionally hide facet strip
  if (!show_strip) {
    base_plot <- base_plot + theme(strip.text = element_blank())
  }
  
  return(base_plot)
}

#Generate plots for each pathogen
plot_list <- lapply(seq_along(pathogens), function(i) {
  name <- pathogens[i]
  show_strip <- (i == 1)                             # Only top plot shows year strip
  plot_pathogen_patch(name, show_strip = show_strip)
})


plot_list <- Filter(Negate(is.null), plot_list)

# Strip x-axis text and label from all but the last
for (i in seq_len(length(plot_list) - 1)) {
  plot_list[[i]] <- plot_list[[i]] +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank()
    )
}

# Combine plots
combined_panel <- wrap_plots(plot_list, ncol = 1)

# Save
ggsave("../Graphs/pathogen_graphs/combined_pathogen_patchwork.pdf", combined_panel, width = 12, height = length(plot_list) * 2)



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


