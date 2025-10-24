library(dplyr)

# Read metadata
metadata <- read.delim("experiment_report_2025_120.tsv",
                       header = TRUE, sep = "\t", skip = 1, stringsAsFactors = FALSE)

# NOTE: ONLY ONE SINGLE SAMPLE has a replicate 
# why idk 
# but this affects the distribution & needs to be considered 

# only keep distinct (unique) Accession IDs
unique_samples <- df %>%
  distinct(Accession, .keep_all = TRUE)

# get sample count per diagnosis group
sample_counts <- unique_samples %>%
  group_by(diagnosis) %>%
  summarise(Count = n())

# visualize distribution 
library(ggplot2)
ggplot(sample_counts, aes(x = diagnosis, y = Count, fill = diagnosis)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Count), vjust = -0.5, fontface = "bold") +
  theme_minimal() +
  labs(title = "Unique Biological Sample Count by Diagnosis Group",
       x = "Diagnosis Group", y = "Sample Count")


ggplot(unique_samples, aes(x = diagnosis, fill = Sex)) +
  geom_bar(position = "stack") +
  geom_text(
    stat = 'count',
    aes(label = ..count..), 
    position = position_stack(vjust = 0.5),
    color = "white",
    fontface = "bold"
  ) +
  labs(title = "Sex Distribution by Diagnosis (Unique Samples)",
       x = "Diagnosis Group", y = "Sample Count") +
  theme_minimal()


