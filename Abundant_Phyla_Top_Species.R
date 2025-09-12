# Look specifically at the Streptophyta reads

#Library -----------------------------------------------------------------------
library(ggplot2)    # Generate visualization plots 
library(dplyr)
library(RColorBrewer)
library(tidyr)

# Data ------
# season monitroing assigned data all samples
data <- read.delim("../all_data_mar25_assigned.tsv", check.names = FALSE) 
sum_data <- read.delim("../all_data_mar25_summed.tsv", check.names = FALSE) 

# Lineage
taxa <- read.csv("../metadata/all_data_marti_assignments_lca_0.1_all_levels_2025-MAR-31_9-49-19_taxaID_lineage.csv")

# Top phyla ------
phyla <- sum_data %>% 
  filter(`NCBI Rank` == "phylum")

# Pivot data 
phyla  <- phyla %>%
  pivot_longer(-c(Name, `NCBI ID`, `NCBI Rank`), names_to = "Sample_ID", values_to = "reads")

# Look at the most abundant across samples  ---------
phyla_abund <- phyla  %>%
  group_by(Name) %>%
  summarise(
    total_reads = sum(reads, na.rm = TRUE),     # overall abundance
    n_samples   = n_distinct(Sample_ID),        # in how many samples it occurs
    mean_reads  = mean(reads, na.rm = TRUE)
  ) %>%
  arrange(desc(total_reads))

# Filtering ---------
# Filter the lineage data 

summarise_marti_by_phyla <- function(taxa, data, phyla, level){ 
  filtered_data <- taxa %>% 
    filter(phylum == phyla)
  
  # Unique level
  unique_level <- filtered_data %>%
    distinct(.data[[level]]) %>%
    pull(1)
  
  # Use these to filter the MARTi data to see which are abundant -------
  marti_filter_level <- data %>%
    filter(Name %in% unique_level)
  
  # Pivot data 
  marti_filter_level <- marti_filter_level %>%
    pivot_longer(-c(Name, `NCBI ID`, `NCBI Rank`), names_to = "Sample_ID", values_to = "reads")
  
  # Look at the most abundant across samples  ---------
  level_abund <- marti_filter_level %>%
    group_by(Name) %>%
    summarise(
      total_reads = sum(reads, na.rm = TRUE),     # overall abundance
      n_samples   = n_distinct(Sample_ID),        # in how many samples it occurs
      mean_reads  = mean(reads, na.rm = TRUE)
    ) %>%
    arrange(desc(total_reads))
  
  return(level_abund)
}

# Call for top 3 phyla at species
strepto_species <- summarise_marti_by_phyla(taxa, sum_data, "Streptophyta", "species")
pseudo_species  <- summarise_marti_by_phyla(taxa, sum_data, "Pseudomonadota", "species")
asco_species    <- summarise_marti_by_phyla(taxa, sum_data, "Ascomycota", "species")

# And for genus  
strepto_genus <- summarise_marti_by_phyla(taxa, sum_data, "Streptophyta", "genus")
pseudo_genus  <- summarise_marti_by_phyla(taxa, sum_data, "Pseudomonadota", "genus")
asco_genus    <- summarise_marti_by_phyla(taxa, sum_data, "Ascomycota", "genus")


