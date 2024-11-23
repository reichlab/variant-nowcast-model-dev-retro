#' Target_data_maker
#'
#' A script to create data files with data available on a given date.

require(tidyverse)
require(arrow)
# reading in the data; if you want to use the 100k sample to test instead of the full metadata
# change the file path and use the commended #full_data_file <- read_tsv(file_path) line,
# because of the different file types.
file_path <- "./input-data/full_metadata_2024-11-22.tsv.zst"
#full_data_file <- read_tsv(file_path)
full_data_file <- read_delim_arrow(file_path, delim = "\t")

# removing any data before 2022, non-humans, and not in the 52 regions we care about
full_data_file <- filter(full_data_file, date_submitted >= as.Date("2022-01-01"),
                         host == "Homo sapiens", division %in% c(state.name, "Washington DC","Puerto Rico"))
# removing all the columns we don't need
full_data_file <- full_data_file[c("division", "Nextstrain_clade", "date", "date_submitted")]
# removing any NA's
full_data_file <- filter(full_data_file, !is.na(date))
# reducing each clade name to there 3 character designation except recombinant
full_data_file$Nextstrain_clade <- substr(full_data_file$Nextstrain_clade,1,3)
full_data_file$Nextstrain_clade[full_data_file$Nextstrain_clade == "rec"] <- rep("recombinant", length(full_data_file$Nextstrain_clade[full_data_file$Nextstrain_clade == "rec"]))
# creating a list of target dates to create versions of, starting for the first Monday in Aug 2022
# and ending on the first Monday in Aug 2024
target_days <- seq(from = as.Date("2022-08-01"), to = as.Date("2024-08-05"), by = 7)
# for each date subsetting by date_submitted, grouping and creating a parquet file.
for(date in target_days[1]){
  date <- as.Date(date)
  data_as_of_date <- full_data_file[full_data_file$date_submitted <= date, ]
  data_to_save <- data_as_of_date %>%
    group_by(division, Nextstrain_clade, date) %>%
    summarize(observation = n(), .groups = "drop")
  data_to_save <- data_to_save %>%
    rename(location = division, clade = Nextstrain_clade)
  file_name <- paste0("./target-data/data_as_of ",date,".parquet" )
  write_parquet(data_to_save, file_name)
}
