#Analysis of the weekly CF data - building on 24hr script
# Want to plot stacked bar charts

#Library -----------------------------------------------------------------------
library(phyloseq)   # Facilitate the import, storage, analysis, and graphical display of microbiome census data.
library(vegan)      # Analysis of variance using distance matrices using the adonis2 function
library(Maaslin2)   # Determine multivariable association between metadata and microbial meta-omics features; differential abundance analysis
library(ggplot2)    # Generate visualization plots 
library(ggsignif)   # Visualize comparisons and the significance between two groups
library(dplyr)
library(RColorBrewer)
library(tidyr)
library(forcats)
library(lubridate)
library(patchwork)

# Set plotting theme
custom_theme <- theme_minimal(base_size = 12) +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank()
  )

# Load in data ------------------------------------------------------------

#OTU table --
#generated the assigned only read tsv from MARTi output using Scripts/split_marti_taxa.py 
# Use assigned reads for phyloseq as builds the objevt using lineage data
otu <- read.delim("../all_data_mar25_assigned.tsv", check.names = FALSE) #2023 & 2024 samples
otu <- otu[,-c(1,3)]          # remove Name & Rank
rownames(otu) <- otu[,1]      # Making the rownames the NCBI ID
otu <- otu[,-1]             # Then drop the repeated col 

#Taxa table --
#generated the lineages using - /Scripts/get_lineage_from_marti.sh from the marti output
taxa <- read.csv("../metadata/all_data_marti_assignments_lca_0.1_all_levels_2025-MAR-31_9-49-19_taxaID_lineage.csv")
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)
#Remove first row (NCBI ID which is empty)
taxa <- taxa[-c(1), ]

rownames(taxa) <- rownames(otu) #same rownames as the otu table
taxa <- taxa[,-1]    


#Metadata --
#Edited in python to deal with missing dates and BST / GMT timezones 
meta <- read.delim("../metadata/metadata_cleaned.tsv")
#Drop Empty SampleID rows
meta <- meta %>% filter(Sample_ID != "" & !is.na(Sample_ID))
rownames(meta) <- meta[,1]       #Make the row names the Sample ID

#tables need to be matrices
otu_mat<- as.matrix(otu)
tax_mat<- as.matrix(taxa)

# Phyloseq ----------------------------------------------------------------

#Make phyloseq object
phylo_OTU<- otu_table(otu_mat, taxa_are_rows = TRUE)
phylo_TAX<- tax_table(tax_mat)
phylo_samples <- sample_data(meta)
phylo_object <- phyloseq(phylo_OTU, phylo_TAX, phylo_samples) #Bring them together

# Functions ------------------------------------------------
#Function to plot ----


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
  p <- ggplot(filtered_data, aes(x = Nice_ID, y = Abundance, fill = phylum)) +  
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colours, drop = drop_unused_phyla) + # for combining, legend will have all phyla
    labs(fill = "Phylum") +
    ylab("Relative Abundance (Phyla > 0.01%) \n") +
    xlab("Date Collected and Replicate") +
    ggtitle(paste(title)) +
    custom_theme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))

  
  # Save the plot
  file_name <- paste0(save_path, "Phylum_", title, ".pdf")
  ggsave(file_name, p, width = 12, height = 6)
  
  return(p)  # Return the plot
}

#Function to plot --
plot_fungi_by_year <- function(data,
                               years,
                               colours = fung_colours,
                               save_path,
                               title,
                               drop_unused_fung = FALSE) {
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
  p <- ggplot(filtered_data, aes(x = Nice_ID, y = Abundance, fill = genus)) +  
    geom_bar(stat = "identity") +
    scale_fill_manual(values = colours, drop = drop_unused_fung) +
    labs(fill = "Genus") +
    ylab("Relative Abundance (Genus > 0.02%) \n") +
    xlab("Date Collected and Replicate") +
    ggtitle(paste(title)) +
    custom_theme +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
  
  # Save the plots 
  file_name <- paste0(save_path, "fungi_", title, ".pdf")
  ggsave(file_name, p, width = 12, height = 6)
  
  return(p)  # Display the plot
}


# Phylum -----------------------------------------------------------------


#Collapse to phylum level and filter on abundance  ----
samp_phylum <- phylo_object %>%
  tax_glom(taxrank = "phylum") %>%                     # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.01) %>%                         # Filter out low abundance taxa <0.01%
  arrange(phylum)                                      # Sort data frame alphabetically by phylum

# Ensure Start_Time is in a proper date format (handle NAs if needed)
samp_phylum$DateTime_UTC <- as.POSIXct(samp_phylum$DateTime_UTC, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

#Get unique phylum in this dataset 
uniq_phyla <- unique(samp_phylum$phylum)
# Get unique phylum names, ensuring 'Higher_Taxa' is last
ordered_phyla <- sort(setdiff(uniq_phyla, "Higher_Taxa"))  # Alphabetical order
ordered_phyla <- c(ordered_phyla, "Higher_Taxa")  # Append 'Higher_Taxa' at the end

# # Assign colors, keeping the existing ones but setting 'Higher_Taxa' to grey
# phy_colours <- colorRampPalette(brewer.pal(8, "Set3"))(length(ordered_phyla) - 1)  # Exclude 'Higher_Taxa'
# phy_colours <- c(phy_colours, "darkgrey")  # Add grey for 'Higher_Taxa'
# phy_colours <- setNames(phy_colours, ordered_phyla)

# Created from above, but decided to edit to increase differences
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

# Ensure phylum column is a factor with ordered levels
samp_phylum$phylum <- factor(samp_phylum$phylum, levels = ordered_phyla)



## Call phylum function on full data ----
plot_phylum_by_year(data = samp_phylum, 
                    years = 2024,
                    colours = phy_colours,
                    save_path = "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                    title = "All_samples_2024",
                    drop_unused_phyla = TRUE)

plot_phylum_by_year(data = samp_phylum, 
                    years = c(2022, 2023),
                    colours = phy_colours,
                    save_path = "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                    title = "2022 - 2023",
                    drop_unused_phyla = TRUE)

## Filtered phylum data ----
# Create filtered df for months where weekly samples were collected in both years (May - Aug)
fil_samp_phylum <- samp_phylum %>% 
  filter(month(DateTime_UTC ) >= 5 & month(DateTime_UTC ) <= 8)

# Ensure Start_Time is in a proper date format (handle NAs if needed)
fil_samp_phylum$DateTime_UTC <- as.POSIXct(fil_samp_phylum$DateTime_UTC, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")
# Ensure phylum column is a factor with ordered levels
fil_samp_phylum$phylum <- factor(fil_samp_phylum$phylum, levels = ordered_phyla)

# get the phyla in the filtered data
present_phyla <- unique(fil_samp_phylum$phylum)

# Subset the colour palette to only those phyla
filtered_colours <- phy_colours[names(phy_colours) %in% present_phyla]

# Run the function
p2024 <- plot_phylum_by_year(data = fil_samp_phylum, 
                             years = 2024,
                             colours = filtered_colours,
                             save_path = "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                             title = "2024",
                             drop_unused_phyla = FALSE) # so the Key can be shared across both 

p2023 <- plot_phylum_by_year(data = fil_samp_phylum,
                             years = 2023,
                             colours = filtered_colours,
                             save_path = "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                             title = "2023",
                             drop_unused_phyla = FALSE) + 
  theme(axis.title.x = element_blank()) # so it isn't repeated in both facets

# Combine them vertically with shared legend
combined_plot <- p2023 / p2024 + plot_layout(guides = "collect") & theme(legend.position = "right")
combined_plot

ggsave("../Graphs/Phylum_Combined_2023_2024.pdf", combined_plot, width = 12, height = 10)


# Fungi -------------------------------------------------------------------


#Filter the data to just be fungi
#First want to filter to genus level - takes a little while to run
samp_genus <- phylo_object %>%
  tax_glom(taxrank = "genus") %>%                      # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa diff threshold to before
  arrange(genus)   

# Get all sample IDs and DateTime_UTC from meta to add back in 
meta_data <- meta %>% select(Nice_ID, DateTime_UTC, Year)

all_samples <- meta_data$Nice_ID

# See which phyla there are & which are fungal
# unique(samp_genus$phylum)

# Filter for these - leaves 25 genera ( maybe too many)
fungal_phyla <- c( "Basidiomycota", "Ascomycota", "Oomycota")
fungi_data <- samp_genus %>% 
  filter(phylum %in% fungal_phyla)

# Ensure all sample IDs are included, even if they have no fungal data
# Complete missing samples keeping DateTime_UTC
fungi_data <- fungi_data %>%
  complete(Nice_ID = all_samples, fill = list(Abundance = 0)) %>%  # Keeps genus values
  mutate(genus = ifelse(is.na(genus), "Higher_Taxa", genus)) %>%  # Fill only missing genus
  left_join(meta_data, by = "Nice_ID") %>%   # Add back DateTime_UTC
  mutate(DateTime_UTC = coalesce(DateTime_UTC.x, DateTime_UTC.y),
         Year = coalesce(Year.x, Year.y)) %>%  # Fill NA values
  select(-DateTime_UTC.x, -DateTime_UTC.y, -Year.x, -Year.y)  # Remove duplicate column

#Get unique genera in this dataset 
uniq_genera <- unique(fungi_data$genus)
# Get unique genera names, ensuring 'Higher_Taxa' is last
ordered_genera <- sort(setdiff(uniq_genera, "Higher_Taxa"))  # Alphabetical order
ordered_genera <- c(ordered_genera, "Higher_Taxa")  # Append 'Higher_Taxa' at the end

# # Assign colors, keeping the existing ones but setting 'Higher_Taxa' to grey
# fung_colours <- colorRampPalette(brewer.pal(8, "Set3"))(length(ordered_genera) - 1)  # Exclude 'Higher_Taxa'
# fung_colours <- c(fung_colours, "darkgrey")  # Add grey for 'Higher_Taxa'
# fung_colours <- setNames(fung_colours, ordered_genera)

fung_colours <- c(Alternaria ="#8DD3C7",
                  Ascochyta =  '#F6A6A4',                                      
                  Botrytis = "#D2EDBA",
                  Cladosporium = "#FBF3A1",
                  Coprinellus = "#F0F0BB",
                  Diaporthe = "#C5B2CC",
                  Dioszegia = "#31A465",
                  Exserohilum = "#DDDBC7",                                           
                  Fulvia = "#D8A0AC",
                  Fusarium = "#AFE0C0",
                  Hyaloperonospora =  "#F58276",
                  Malassezia = "#D09193",
                  Neohortaea = "#9D7BB0",
                  Parastagonospora = "#85AECE",
                  Pestalotiopsis = "#C6B293",                                   
                  Pleurotus =  "#C3D467",
                  Puccinia = "#D7C865",
                  Pyrenophora = '#D5B030',
                  Ramularia = "#F5A469",
                  Septoria =  "#8FCF3C" ,
                  Peronospora = '#4B6587',
                  Trametes ="#CFD799",
                  Ustilago = '#A9746E',
                  Zymoseptoria =  "#FCCDE5",
                  Higher_Taxa = "darkgrey"
                  )


# Ensure genus is a factor with correct ordering
fungi_data$genus <- factor(fungi_data$genus, levels = ordered_genera)

# Ensure Start_Time is in a proper date format (handle NAs if needed)
fungi_data$DateTime_UTC <- as.POSIXct(fungi_data$DateTime_UTC, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")


## Plot fungi data ----
plot_fungi_by_year(data = fungi_data,
                   years = 2024,
                   colours = fung_colours,
                   save_path = "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                   title = "All 2024",
                   drop_unused_fung = TRUE)

plot_fungi_by_year(data = fungi_data, 
                   years = c(2022, 2023),
                   colours = fung_colours,
                   save_path =  "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                   title = "2022 - 2023",
                   drop_unused_fung = TRUE)


## Filtered fungi genus data ----
# Create filtered df for months where weekly samples were collected in both years (May - Aug)
fil_samp_genus <- fungi_data %>% 
  filter(month(DateTime_UTC ) >= 5 & month(DateTime_UTC ) <= 8)

# Ensure genus column is a factor with ordered levels
fil_samp_genus$genus <- factor(fil_samp_genus$genus, levels = ordered_genera)

# get the phyla in the filtered data
present_genera <- unique(fil_samp_genus$genus)

# Subset the colour palette to only those phyla
filter_fung_colours <- fung_colours[names(fung_colours) %in% present_genera]

# Run the function
f2024 <- plot_fungi_by_year(data = fil_samp_genus, 
                             years = 2024,
                             colours = filter_fung_colours,
                             save_path = "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                             title = "2024",
                             drop_unused_fung = FALSE) # so the Key can be shared across both 

f2023 <- plot_fungi_by_year(data = fil_samp_genus,
                             years = 2023,
                             colours = filter_fung_colours,
                             save_path = "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                             title = "2023",
                             drop_unused_fung = FALSE) + 
  theme(axis.title.x = element_blank()) # so it isn't repeated in both facets

# Combine them vertically with shared legend
combined_plot <- f2023 / f2024 + plot_layout(guides = "collect") & theme(legend.position = "right")
combined_plot

ggsave("../Graphs/fung_genus_Combined_2023_2024.pdf", combined_plot, width = 12, height = 10)



