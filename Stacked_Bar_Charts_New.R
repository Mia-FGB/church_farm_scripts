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

# Set plotting theme
theme_set(theme_bw())

# Load in data ------------------------------------------------------------

#OTU table --
#generated the assigned only read tsv from MARTi output using Scripts/split_marti_taxa.py 
# Use assigned reads for phyloseq as builds the objevt using lineage data
otu <- read.delim("all_data_mar25_assigned.tsv", check.names = FALSE) #2023 & 2024 samples
otu <- otu[,-c(1,3)]          # remove Name & Rank
rownames(otu) <- otu[,1]      # Making the rownames the NCBI ID
otu <- otu[,-1]             # Then drop the repeated col 

#Taxa table --
#generated the lineages using - /Scripts/get_lineage_from_marti.sh from the marti output
taxa <- read.csv("metadata/all_data_marti_assignments_lca_0.1_all_levels_2025-MAR-31_9-49-19_taxaID_lineage.csv")
taxa<- taxa %>% 
  mutate_if(is.character, as.factor)
#Remove first row (NCBI ID which is empty)
taxa <- taxa[-c(1), ]

rownames(taxa) <- rownames(otu) #same rownames as the otu table
taxa <- taxa[,-1]    


#Metadata --
#Edited in python to deal with missing dates and BST / GMT timezones 
meta <- read.delim("metadata/metadata_cleaned.tsv")
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

# Plotting -----------------------------------------------------------------

# Phylum --

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

# Assign colors, keeping the existing ones but setting 'Higher_Taxa' to grey
phy_colours <- colorRampPalette(brewer.pal(8, "Set3"))(length(ordered_phyla) - 1)  # Exclude 'Higher_Taxa'
phy_colours <- c(phy_colours, "darkgrey")  # Add grey for 'Higher_Taxa'
phy_colours <- setNames(phy_colours, ordered_phyla)

# Ensure phylum column is a factor with ordered levels
samp_phylum$phylum <- factor(samp_phylum$phylum, levels = ordered_phyla)

#Function to plot --
plot_phylum_by_year <- function(data, years, save_path, title) {
  # Filter data by the selected year
  filtered_data <- data %>% filter(Year %in% years)
  
  # Order data by DateTime_UTC
  filtered_data <- filtered_data %>% arrange(DateTime_UTC)
  
  # Ensure Nice_ID is uniquely ordered by Nice_ID
  filtered_data$Nice_ID <- factor(filtered_data$Nice_ID, levels = unique(filtered_data$Nice_ID))
  
  # Create folder if it doesn't exist
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  # Create the stacked bar plot
  p <- ggplot(filtered_data, aes(x = Nice_ID, y = Abundance, fill = phylum)) +  
    geom_bar(stat = "identity") +
    scale_fill_manual(values = phy_colours) +
    theme(axis.title.x = element_blank()) +
    labs(fill = "Phylum") +
    ylab("Relative Abundance (Phyla > 0.01%) \n") +
    ggtitle(paste("Phylum Composition -", title)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
  
  # Save the plot
  file_name <- paste0(save_path, "Phylum_", title, ".svg")
  ggsave(file_name, p, width = 12, height = 6)
  
  print(p)  # Display the plot
}

#Run function:
plot_phylum_by_year(samp_phylum, 2024,
                    "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                    "2024")
plot_phylum_by_year(samp_phylum, c(2022, 2023),
                    "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                    "2022 - 2023")



# Fungi -------------------------------------------------------------------


#Filter the data to just be fungi
#First want to filter to genus level - takes a little while to run
samp_genus <- phylo_object %>%
  tax_glom(taxrank = "genus") %>%                      # agglomerate at phylum level
  transform_sample_counts(function(x) {x/sum(x)} ) %>% # Transform to rel. abundance - compared to pass reads
  psmelt() %>%                                         # Melt to long format
  filter(Abundance > 0.02) %>%                         # Filter out low abundance taxa diff threshold to before
  arrange(genus)                                      

#See which phyla there are & which are fungal
unique(samp_genus$phylum)
# Filter for these - leaves 25 genera ( maybe too many)
fungal_phyla <- c( "Basidiomycota", "Ascomycota", "Oomycota")
fungi_data <- samp_genus %>% 
  filter(phylum %in% fungal_phyla)

#Get unique genera in this dataset 
uniq_genera <- unique(fungi_data$genus)
# Get unique genera names, ensuring 'Higher_Taxa' is last
ordered_genera <- sort(setdiff(uniq_genera, "Higher_Taxa"))  # Alphabetical order
ordered_genera <- c(ordered_genera, "Higher_Taxa")  # Append 'Higher_Taxa' at the end

# Assign colors, keeping the existing ones but setting 'Higher_Taxa' to grey
fung_colours <- colorRampPalette(brewer.pal(8, "Set3"))(length(ordered_genera) - 1)  # Exclude 'Higher_Taxa'
fung_colours <- c(fung_colours, "darkgrey")  # Add grey for 'Higher_Taxa'
fung_colours <- setNames(fung_colours, ordered_genera)

# Ensure genus is a factor with correct ordering
fungi_data$genus <- factor(fungi_data$genus, levels = ordered_genera)

# Ensure Start_Time is in a proper date format (handle NAs if needed)
fungi_data$DateTime_UTC <- as.POSIXct(fungi_data$DateTime_UTC, format = "%Y-%m-%d %H:%M:%S", tz = "UTC")

#Function to plot --
plot_fungi_by_year <- function(data, years, save_path, title) {
  # Filter data by the selected year
  filtered_data <- data %>% filter(Year %in% years)
  
  # Order data by DateTime_UTC
  filtered_data <- filtered_data %>% arrange(DateTime_UTC)
  
  # Ensure Nice_ID is uniquely ordered by Nice_ID
  filtered_data$Nice_ID <- factor(filtered_data$Nice_ID, levels = unique(filtered_data$Nice_ID))
  
  # Create folder if it doesn't exist
  if (!dir.exists(save_path)) {
    dir.create(save_path, recursive = TRUE)
  }
  
  # Create the stacked bar plot
  p <- ggplot(filtered_data, aes(x = Nice_ID, y = Abundance, fill = genus)) +  
    geom_bar(stat = "identity") +
    scale_fill_manual(values = fung_colours) +
    theme(axis.title.x = element_blank()) +
    labs(fill = "Genus") +
    ylab("Relative Abundance (Genus > 0.02%) \n") +
    ggtitle(paste("Fungal Genera Composition -", title)) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0.5))
  
  # Save the plot
  file_name <- paste0(save_path, "fungi_", title, ".svg")
  ggsave(file_name, p, width = 12, height = 6)
  
  print(p)  # Display the plot
}

#Run function:
plot_fungi_by_year(fungi_data, 2024,
                    "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                    "2024")
plot_fungi_by_year(fungi_data, c(2022, 2023),
                    "/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm/Graphs/",
                    "2022 - 2023")

