#' Target_data_maker
#'
#' A script to create data files with data available on a given date.

require(tidyverse)
require(arrow)
# reading in the data
file_path <- "./auxiliary-data/raw-target-data/metadata.csv"
full_data_file <- read_csv(file_path)

# reducing each clade name to there 3 character designation except recombinant
full_data_file <- full_data_file %>%
   mutate(clade = ifelse(clade == "recombinant", "recombinant", substr(clade, 1, 3)))
# creating a list of target dates to create versions of, starting for the first Monday in Aug 2022
# and ending on the first Monday in Aug 2024
target_days <- seq(from = as.Date("2022-08-01"), to = as.Date("2024-08-05"), by = 7)
# for each date subsetting by date_submitted, grouping and creating a parquet file.
# matching the structure of write_dataset
file_name <- "part-0.parquet"
for(date in target_days[1]) {
  date <- as.Date(date)
  data_as_of_date <- full_data_file[full_data_file$date_submitted <= date, ]
  data_to_save <- data_as_of_date %>%
    group_by(location, clade, date) %>%
    summarize(observation = n(), .groups = "drop")
  #creating the partition folders
  folder_name <- paste0("./target-data/as_of=", date)
  if(!file.exists(folder_name)) {
    dir.create(folder_name)
  }
  write_parquet(data_to_save, file.path(folder_name, file_name))
}
