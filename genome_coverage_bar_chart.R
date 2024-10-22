library(dplyr)
library(ggplot2)
library(readr)
library(stringr)
library(lubridate)
library(RColorBrewer)

#Reading in the data ---------
z_tritici <- read_tsv("Regular_Collections_2023/genome_coverage/z_tritici_coverage.tsv")
p_nodorum <- read_tsv("Regular_Collections_2023/genome_coverage/p_nodorum_coverage.tsv")
p_avenae <- read_tsv("Regular_Collections_2023/genome_coverage/p_avenae_coverage.tsv")

metadata <- read_csv("CF_2023_metadata.csv")


#Add reference genome length --------
z_tritici$ref_length_bp <-39686251
p_nodorum$ref_length_bp <-37213987
p_avenae$ref_length_bp  <-36888127

#Merge the datasets -------
combined_data <- bind_rows(z_tritici, p_nodorum, p_avenae)

#Calculate coverage -------
combined_data <- combined_data %>% 
  mutate(coverage = (`Number of Mapped Reads`/ref_length_bp)*100) %>% 
# Extract only the number from the 'barcode' column
  mutate(Barcode = str_extract(Barcode, "\\d+"))

#Add relevant metadata ---------
merged_data <- merge(combined_data, metadata, by = "Barcode")

# Fixing the date column
merged_data  <- merged_data  %>%
  mutate(`Date collected` = ifelse(grepl("^[A-Za-z]{3}-\\d{2}$", `Date collected`),  # Check for format like 'Apr-23'
                                   paste0("15-", `Date collected`),                   # Append '-15' to create '15-Apr-23'
                                   `Date collected`),                                  # Otherwise leave the date unchanged
         # Now, try to parse both types of dates using parse_date_time, which handles multiple formats
         `Date collected` = parse_date_time(`Date collected`, orders = c("dmy", "mdy")))

merged_data$`Date collected` <- as.Date(merged_data$`Date collected`)

#Plot coverage over time--------
#facet by reference genome 

# Group by Date collected and reference genome, and calculate the mean coverage
merged_data_avg <- merged_data %>%
  group_by(`Reference Genome`, `Date collected`) %>%
  summarise(mean_coverage = mean(coverage), min_coverage = min(coverage),
            max_coverage = max(coverage),
            mean_mapped_reads = mean(`Number of Mapped Reads`),
            min_mapped_reads = min(`Number of Mapped Reads`),
            max_mapped_reads = max(`Number of Mapped Reads`))

#transform the coverage data for better visualisation
merged_data_avg$sqrt_coverage <- sqrt(merged_data_avg$mean_coverage + 0.001)
merged_data_avg$sqrt_min_coverage <- sqrt(merged_data_avg$min_coverage + 0.001)
merged_data_avg$sqrt_max_coverage <- sqrt(merged_data_avg$max_coverage + 0.001)

#Setting x-axis start and end
break.vec <- c(as.Date("2022-10-01"), seq(from=as.Date("2022-10-01"), to=as.Date("2023-08-31"),by="1 week")) 

#Bar plots ----------

#Genome coverage plot!
gc <-  ggplot(
  merged_data_avg, aes(x = `Date collected`, y = mean_coverage, fill = `Reference Genome`)) + 
  #facet_wrap(~ `Reference Genome`, ncol = 1) +  
  geom_bar(stat = "identity",  position = position_dodge()) +
  geom_errorbar(
    aes(ymin = min_coverage, ymax = max_coverage),
    width = 4, size = 0.4, colour = "grey39", position = position_dodge()) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Set3")) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
    axis.title.x = element_blank(), #remove x-axis title
    axis.line = element_line(size = 0.7)) + #thicker axis
  scale_x_date(breaks = break.vec, date_labels = "%d-%b") +
  ylab("Average Genome coverage %") +
  xlab("Date Collected") +
  ggtitle("Genome coverage over time") 

# Save as SVG
ggsave(filename = "Regular_Collections_2023/graphs/genome_coverage.svg", plot = gc, width = 10, height = 6)


#Genome coverage sqrt transformed plot!
#I don't think this provides any better info than the standard plot
tf <- ggplot(
  merged_data_avg, aes(x = `Date collected`, y = sqrt_coverage, fill = `Reference Genome`)) + 
  #facet_wrap(~ `Reference Genome`, ncol = 1) +  
  geom_bar(stat = "identity",  position = position_dodge()) +
  geom_errorbar(
    aes(ymin = sqrt_min_coverage, ymax = sqrt_max_coverage),
    width = 4, size = 0.4, colour = "grey39", position = position_dodge()) +
  ylab("Average Genome coverage % (sqrt transformed)") +
  xlab("Date Collected") +
  ggtitle("Genome coverage over time") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Set3")) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
    axis.title.x = element_blank(), #remove x-axis title
    axis.line = element_line(size = 0.7)) + #thicker axis
  scale_x_date(breaks = break.vec, date_labels = "%d-%b")

# Save as SVG
ggsave(filename = "Regular_Collections_2023/graphs/genome_coverage_sqrt.svg", plot = tf, width = 10, height = 6)
 
#Number mapped read plot!
mr <-  ggplot(
   merged_data_avg, aes(x = `Date collected`, y = mean_mapped_reads, fill = `Reference Genome`)) + 
   #facet_wrap(~ `Reference Genome`, ncol = 1) +  
   geom_bar(stat = "identity",  position = position_dodge()) +
   geom_errorbar(
     aes(ymin = min_mapped_reads, ymax = max_mapped_reads),
     width = 4, size = 0.4, colour = "grey39", position = position_dodge()) +
   ylab("Average Number of Mapped Reads") +
   xlab("Date Collected") +
   ggtitle("Number of Mapped Reads over time") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_fill_manual(values = brewer.pal(n = 3, name = "Set3")) +
  theme_bw() + 
  theme(
    panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.border = element_blank(), 
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=0.5),
    axis.title.x = element_blank(), #remove x-axis title
    axis.line = element_line(size = 0.7)) + #thicker axis
  scale_x_date(breaks = break.vec, date_labels = "%d-%b")

 # Save as SVG
 ggsave(filename = "Regular_Collections_2023/graphs/mapped_reads.svg", plot = mr, width = 10, height = 6)
