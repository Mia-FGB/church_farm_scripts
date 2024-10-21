setwd("/Users/berelsom/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Air_Samples/Church_farm")

library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(dplyr)
library(lubridate)
library(scales)
library(patchwork)

#Read in the dataframe, this was generated in Scripts/Initial_CF_2023_analysis.ipynb from the MARTi output
#This data has not been filtered
data <- read.csv("Regular_Collections_2023/just_pathogen_data.csv", header = TRUE)

#sort the dates - Nice.Date is the Date collected but easier to read in
data$Nice.Date = as.Date(data$Nice.Date, format="%d/%m/%y")

# Read number table grouped by collection date ---------------------------
read_num <- data %>% 
  group_by(Nice.Date, Name) %>% #Name is the diff genera
  summarise(mean_read_num = mean(read.count), mean_rel_abundance = mean(rel_abundance))


#Graphs ---------------------------------------------

#mean, min & max data for each genus for errorbar just HPM
# data_min_max <- data %>%
#   group_by(Nice.Date, Name) %>%
#   summarise(mean = mean(HPM), min = min(HPM), max = max(HPM))

data$ln_HPM[is.infinite(data$ln_HPM)] <- 0

data_min_max <- data %>%
  group_by(Nice.Date, Name) %>%
  summarise(mean = mean(ln_HPM), min = min(ln_HPM), max = max(ln_HPM))

#Setting x-axis start and end
break.vec <- c(as.Date("2022-10-01"), seq(from=as.Date("2022-10-01"), to=as.Date("2023-08-31"),by="1 month")) 


# Function to create plots for a specific genus and save them to a PDF - need to imrpove the aesthetics

# create_genus_plot <- function(genus_name) {
#   p <- ggplot(data_min_max[data_min_max$Name == genus_name, ], aes(x = Nice.Date, y = mean)) +
#    # facet_grid(cols = vars(location), scales ="free_x") + 
#     geom_line(size = 1.2) + 
#     geom_point(size = 4) +
#     geom_errorbar(aes(ymin = min, ymax = max), width = 10, colour = "grey39") +
#     theme_bw() + 
#     theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
#     xlab("Date collected (Week)") + 
#     ylab("Hits per Million") +
#     scale_x_date(breaks = break.vec, date_labels = "%d/%m/%y", limits = range(break.vec)) +
#     ggtitle(paste("Hits per Million - Genus:", genus_name))
#   return(p)
# }
# 
# # Creating plots for each genus, pull out all the unique genera in the dataframe
# genus_list <- unique(data_min_max$Name)
# 
# for (genus in genus_list) {
#   plot <- create_genus_plot(genus)
#   filename <- paste0("Regular_Collections_2023/graphs/R_plots", genus, ".pdf")
#   ggsave(filename, plot = plot, width = 24, height = 5, units = "in")
# }


#Facet plots for all genera
p <- ggplot(data_min_max, aes(x = Nice.Date, y = mean)) + 
  facet_wrap(~ Name, ncol = 2) +  
  geom_line(size = 0.6, color = "#014d4e") + 
  geom_point(size = 1, color = "#014d4e") +
  geom_errorbar(aes(ymin = min, ymax = max), width = 10, size = 0.6, colour = "grey39") +
  theme_bw() + 
  theme(
    panel.grid.minor = element_blank(), 
    panel.grid.major = element_blank(),
  theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) +
  xlab("Date collected") + 
  ylab("Hits per Million (Ln)") +
  scale_x_date(breaks = break.vec, date_labels = "%b", limits = range(break.vec))

# Save as SVG
ggsave(filename = "Regular_Collections_2023/graphs/sep_pathogens_line_error.svg", plot = p, width = 10, height = 9)
