
### Generate manifest file for FragPipe

library(here)
library(tidyverse)
library(here)
library(magrittr)


# Get list of mzML files with full path names
ms_data_dir <- here("data", "ms_data", "mzml")
ms_files <- list.files(ms_data_dir, pattern = "mzML", full.names = TRUE)

# create data frame for fragpipe manifest
# extract experiments (condition_time) from file names 
fp_manifest <- data.frame(file = ms_files, stringsAsFactors = FALSE) %>%
  mutate(experiment = str_extract(file, "N_[A-Z]_t[0-5]"),
         # extract bioreplicates from file names
         bioreplicate = gsub(".mzML", "", str_extract(file, "[1-3].mzML$")),
         # add column with data type "DDA"
         data_type = "DDA") %>%
  # select rows where experiment is not NA
  filter(!is.na(experiment)) %>%
  set_colnames(c("file", "experiment", "bioreplicate", "data_type"))

# export tab separated file with the name "nepenthes.fp-manifest"
write.table(fp_manifest, file = here("data" ,"nepenthes.fp-manifest"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# repeat for sarracenia
fp_manifest <- data.frame(file = ms_files, stringsAsFactors = FALSE) %>%
  mutate(experiment = str_extract(file, "S_[A-Z]_t[0-5]"),
         bioreplicate = gsub(".mzML", "", str_extract(file, "[1-3].mzML$")),
         data_type = "DDA") %>%
  filter(!is.na(experiment)) %>%
  set_colnames(c("file", "experiment", "bioreplicate", "data_type"))

write.table(fp_manifest, file = here("data", "sarracenia.fp-manifest"),
            sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

