#setting working directory
setwd("~/Library/CloudStorage/OneDrive-NorwichBioScienceInstitutes/Church farm")
#loading libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(stringr)
library(stringi)
library(tidyr)
#Reading in the data file
data <- read.csv("cfarm_read_count.csv", header=TRUE)


#--Normalising the data---------------

#Need to normalise the reads for each location, e.g. hits per million reads
#Therefore we need the total number of reads for each sample, not just the classified reads
data <- 
  data %>% 
  #only include species woth more than 5 reads
  filter(count >=5) %>% 
  mutate(hits_per_100.000 = (count * 100000)/(total_read)) %>% 
  #log scale-Mutate the data to add a new row
  mutate(log_hits = log10(hits_per_100.000))

#Graph
cfarm.heatmap.phi <- ggplot(data = data, mapping = aes(x = label,
                                                        y = species,
                                                        fill = hits_per_100.000)) +
  geom_tile() +
  xlab(label = "Location") + ylab(label = "Species") +
  scale_fill_distiller(name = "Hits per 100,000", palette = "RdBu") +
  theme_bw() +
  ggtitle(label = "Church Farm Air-Seq Species Abundance (hits per 100,000 >5)") +
  #Editing axis title size
  theme(axis.text.y = element_text(vjust = 0.5, hjust=1, size = 15,)) +
  theme(axis.title.y = element_text(size = 20, vjust = 0.5)) +
  theme(axis.title.x = element_text(size = 20, hjust = 0.5)) +
  theme(axis.text.x = element_text(vjust=0.5, size = 15,)) +
  # Change legend key size and key width
  theme(legend.key.size = unit(6, "cm"),legend.key.width = unit(1.5,"cm"))+
  #Change legend font size
  theme(legend.title = element_text(size = 15)) +
  theme(legend.text = element_text(size =15)) +
  #Changing the title size
  theme(plot.title = element_text(size = 20, hjust = 0.5))

cfarm.heatmap.phi
#Saving the plot to right dimensions, much taller y axis
ggsave(filename= "cfarm_initial_sampling_heatmap_phibase.png", plot = cfarm.heatmap.phi, width=20, height=20)




