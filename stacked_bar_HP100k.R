# Stacked Phylum level pltos from MARTi output 

#Library -----------------------------------------------------------------------
library(ggplot2)    # Generate visualization plots 
library(ggsignif)   # Visualize comparisons and the significance between two groups
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(forcats)
library(lubridate)
library(patchwork)

# Using this script as inspo 
# https://github.com/Mia-FGB/Sampler_Comparison_Scripts/blob/main/top_taxa_graphs.R

# Read in data ----------
marti <- read.delim("../all_data_mar25_summed.tsv", check.names = FALSE) #2023 & 2024 samples
meta <- read.delim("../metadata/metadata_cleaned.tsv")
read_numbers <- read.delim("../metadata/read_numbers.tsv")

# Data wrangling -----------

## MARTi -----
marti_long <- marti %>%
  select(-`NCBI ID`) %>%                                   # drop ID
  pivot_longer(
    cols = -c(Name, `NCBI Rank`),                          # pivot everything except these
    names_to  = "Sample_ID",
    values_to = "Count"
  )

marti_long$Count <- as.numeric(marti_long$Count)

## Merge on metadata ----------

# Merge on read numbers and calculate

marti_num <- marti_long %>%
  left_join(read_numbers, by = c("Sample_ID" = "ID"))

marti_num <- marti_num %>% 
  mutate(
    percent_classified_read = Count / ReadsPassBasecall * 100,
    HPM = Count / ReadsPassBasecall * 1e6,  #hits per million
    HP100k = Count / ReadsPassBasecall * 100000
  )

# metadata
marti_meta <- marti_num %>%
  left_join(meta, by = "Sample_ID")


marti_meta <- marti_meta %>%
  mutate(DateTime_UTC = as.POSIXct(DateTime_UTC,
                                   format = "%Y-%m-%d %H:%M:%S",
                                   tz = "UTC"))

# Filter the data on read count and HPM
filtered_marti <- marti_meta  %>% 
  filter(HP100k > 1000) # to match previous graphs

# Create filtered df for months where weekly samples were collected in both years (May - Aug)
fil_marti_date <- filtered_marti %>% 
  filter(month(DateTime_UTC ) >= 5 & month(DateTime_UTC ) <= 8)


# Phylum -----------

# This is missing higehr taxa, also can deal with that later if this graph looks better
phylum <- fil_marti_date  %>% 
  filter(`NCBI Rank` == "phylum")

# Filtering has produced different results to the other graph, less phyla now?
sort(unique(phylum$Name))

# But keep colours the same for now 
phy_colours <- c(
  Acidobacteriota = "#7BCED9",
  Actinomycetota = "#BEE6BE",
  Arthropoda = "#FBF3A1",
  Ascomycota = "#E49599",
  Bacillota = "#CECBD0",
  Bacteroidota = "#B0A7DA",
  Basidiomycota = "#FFD8C4",
  Chordata = "#D96C4E",
  Cressdnaviricota = "#BD98A2",
  Cyanobacteriota = "#87ADCC",
  Deinococcota = "#AEB2A8",
  Mollusca = "#E5B377",
  Myxococcota = "#EABA63",
  Oomycota = "#CAD066",
  Planctomycetota = "#ACDB78",
  Pseudomonadota = "#DCD4AE",
  Streptophyta = "#FCCDE5",
  Higher_Taxa = "darkgrey"
)

plot_phylum_by_year <- function(data,
                                years,
                                colours = phy_colours,
                                save_path,
                                title,
                                drop_unused_phyla = FALSE) {
  # Filter data by the selected year
  filtered_data <- data %>% filter(Year %in% years)
  
  # Order data by DateTime_UTC
  filtered_data <- filtered_data %>% arrange(DateTime_UTC, Repeat)
  
  # Ensure Nice_ID is uniquely ordered by Nice_ID
  filtered_data$Nice_ID <- factor(filtered_data$Nice_ID, levels = unique(filtered_data$Nice_ID))
  
  # Create folder if it doesn't exist
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  # Create the stacked bar plot
  p <- ggplot(filtered_data, aes(x = Nice_ID, y = HP100k, fill = Name)) +  
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colours, drop = drop_unused_phyla) + # for combining, legend will have all phyla
    labs(fill = "Phylum") +
    ggtitle(paste(title)) +
    ylab("Hits per 100k") +
    xlab("Date Collected and Replicate") +
    custom_theme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
  
  
  # Save the plot
  file_name <- paste0(save_path, "Phylum_HP100k", title, ".pdf")
  # ggsave(file_name, p, width = 12, height = 6)
  
  return(p)  # Return the plot
}

# Set plotting theme
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

# Ensure phylum column is a factor with ordered levels - not for now
# fil_samp_phylum$Name<- factor(fil_samp_phylum$Name, levels = ordered_phyla)

# get the phyla in the filtered data
present_phyla <- unique(phylum$Name)

# Subset the colour palette to only those phyla
filtered_colours <- phy_colours[names(phy_colours) %in% present_phyla]

# Run the function
p2024 <- plot_phylum_by_year(data = phylum, 
                             years = 2024,
                             colours = filtered_colours,
                             save_path = "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                             title = "2024",
                             drop_unused_phyla = FALSE) # so the Key can be shared across both 

p2023 <- plot_phylum_by_year(data = phylum,
                             years = 2023,
                             colours = filtered_colours,
                             save_path = "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                             title = "2023",
                             drop_unused_phyla = FALSE) + 
  theme(axis.title.x = element_blank()) # so it isn't repeated in both facets

# Combine them vertically with shared legend
combined_plot <- p2023 / p2024 + plot_layout(guides = "collect") & theme(legend.position = "right")
combined_plot

ggsave("../Graphs/Phylum_Combined_2023_2024_HP100k.pdf", combined_plot, width = 12, height = 10)


# Fungal and oomycete ---------

# Filter for these 
fungal_phyla <- c( "Basidiomycota", "Ascomycota", "Oomycota")

# total reads per sample
totals <- phylum %>%
  distinct(Sample_ID, total_reads = ReadsPassBasecall)

# reads per phylum per sample
reads_by <- phylum %>%
  filter(Name %in% fungal_phyla) %>%
  group_by(Sample_ID, Name) %>%
  summarise(reads = sum(Count, na.rm = TRUE), .groups = "drop")

# add explicit zeros for missing phyla within each sample, then compute relative abundance
abundance_by_sample <- reads_by %>%
  complete(Sample_ID, Name = fungal_phyla, fill = list(reads = 0)) %>%
  left_join(totals, by = "Sample_ID") %>%
  mutate(abundance = (reads / total_reads) * 100)

# now summarise across samples
phylum_stats_zero_included <- abundance_by_sample %>%
  group_by(Name) %>%
  summarise(
    n_samples = n(),
    mean_abundance = mean(abundance, na.rm = TRUE),
    se = sd(abundance, na.rm = TRUE) / sqrt(n()),
  )
