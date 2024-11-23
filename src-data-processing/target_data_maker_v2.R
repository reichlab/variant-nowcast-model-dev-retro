require(tidyverse)
require(arrow)
# reading in the data; if you want to use the 100k sample to test instead of the full metadata
# change the file path and use the commended #full_data_file <- read_tsv(file_path) line,
# because of the different file types.
file_path <- "./input-data/full_metadata_2024-11-22.tsv.zst"
#full_data_file <- read_tsv(file_path)
full_data_file <- read_delim_arrow(file_path, delim = "\t")

# removing any data submitted before 2022, after 2024-08-05, non-humans, and not in the 52 regions we care about
full_data_file <- filter(full_data_file, date_submitted >= as.Date("2022-01-01") & date_submitted <= as.Date("2024-08-05"),
                         host == "Homo sapiens", division %in% c(state.name, "Washington DC","Puerto Rico"))
# removing all the columns we don't need
full_data_file <- full_data_file[c("division", "Nextstrain_clade", "date", "date_submitted")]
# removing any NA's from the data
full_data_file <- filter(full_data_file, !is.na(date))
# reducing each clade name to their 3 character designation, except recombinant
full_data_file$Nextstrain_clade <- substr(full_data_file$Nextstrain_clade,1,3)
full_data_file$Nextstrain_clade[full_data_file$Nextstrain_clade == "rec"] <- rep("recombinant", length(full_data_file$Nextstrain_clade[full_data_file$Nextstrain_clade == "rec"]))
# creating a list of target dates to create versions of, starting for the first Monday in Aug 2022
# and ending on the first Monday in Aug 2024
target_days <- seq(from = as.Date("2022-08-01"), to = as.Date("2024-08-05"), by = 7)
# creating a as_of column that denotes what week the seq was added
full_data_file$as_of <- rep(as.Date("2022-08-01"), length(full_data_file$Nextstrain_clade))

for(t_date in target_days){
  t_date <- as.Date(t_date)
  full_data_file$as_of <- if_else(full_data_file$date_submitted >= as.Date(t_date - 6), t_date,full_data_file$as_of)
}
# grouping the data and renaming some of the columns
data_to_save <-  full_data_file %>%
  group_by(division, Nextstrain_clade, date, as_of) %>%
  summarize(observation = n(), .groups = "drop")
data_to_save <- data_to_save %>%
  rename(location = division, clade = Nextstrain_clade)
# writing the datasets, here each dataset except the first contains all sequences submitted
# in the week before it's date; i.e, 2024-07-08 contains all the sequences submitted between
# 2024-07-01 and 2024-07-08. The first dataset 2022-08-01 contains all sequences submitted
# in 2022 before 2022-08-01.
write_dataset(dataset = data_to_save, path = "./target-data/",
              format = "parquet", partitioning = "as_of" )
