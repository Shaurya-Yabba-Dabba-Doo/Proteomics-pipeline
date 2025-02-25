library(dplyr)
library(stringr)
library(tidyr)
library(tibble)
library(readr)

# Load data
full_data <- read_tsv(
    system.file("extdata/proteinGroups.txt", 
                package = "proDA", mustWork = TRUE),
    col_types = cols(Reverse = col_character())
)

print(head(full_data))

# Process data
tidy_data <- full_data %>%
    dplyr::select(ProteinID = Protein.IDs, starts_with("LFQ.intensity.")) %>%
    gather(Sample, Intensity, starts_with("LFQ.intensity.")) %>%
    mutate(Condition = str_match(Sample, "LFQ\\.intensity\\.([[:alnum:]]+)\\.\\d+")[, 2]) %>%
    mutate(Replicate = as.numeric(str_match(Sample, "LFQ\\.intensity\\.[[:alnum:]]+\\.(\\d+)")[, 2])) %>%
    mutate(SampleName = paste0(Condition, ".", Replicate))

print(head(tidy_data))

# Transform intensity values
data <- tidy_data %>%
    mutate(Intensity = ifelse(Intensity == 0, NA, log2(as.numeric(Intensity)))) %>%
    dplyr::select(ProteinID, SampleName, Intensity) %>%
    spread(SampleName, Intensity) %>%
    column_to_rownames("ProteinID") %>%
    as.matrix()

# Display a subset of the data matrix
print(data[1:4, 1:7])

# Create annotation df
annotation_df <- tidy_data %>%
    dplyr::select(SampleName, Condition, Replicate) %>%
    distinct() %>%
    arrange(Condition, Replicate) %>%
    as.data.frame() %>%
    column_to_rownames("SampleName")

# Display the annotation dataframe
print(annotation_df)
