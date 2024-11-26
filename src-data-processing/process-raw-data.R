## script to filter and save metadata on git lfs

## the full open `metadata.tsv` file will be downloaded, made smaller,
## and placed into the raw-target-data folder

## some code adapted from : https://tutorials.inbo.be/tutorials/r_large_data_files_handling/
tsv.name <- "./auxiliary-data/raw-target-data/metadata.tsv"
tsv.url <- paste("https://data.nextstrain.org/files/ncov/open/metadata.tsv.gz",
                 sep = "")
if (!file.exists(tsv.name)) {
  download.file(tsv.url, destfile = paste0(tsv.name, ".gz"))
  R.utils::gunzip(paste0(tsv.name, ".gz"))
}

## read data, with only the necessary columns
dat <- data.table::fread(tsv.name,
                         select = c(division = "character",
                                    Nextstrain_clade = "character",
                                    date = "character",
                                    date_submitted = "Date",
                                    host = "character"))
dat <- dat[, date_new := as.Date(date, format = "%Y-%m-%d")]

## the thinking here is that we only include sequences within 150 days of the first forecast date
dat <- dat[date_new >= as.Date("2022-08-01") - 150,]

## exclude any rows that have date_submitted as NA
dat <- dat[!is.na(date_submitted),]

## only include human specimens
dat <- dat[host == "Homo sapiens",]

## final subset and also renaming columns
## now that we have filtered by host, excluding that column
dat <- dat[division %in% c(state.name, "Washington DC","Puerto Rico"),
           .(location = division, clade = Nextstrain_clade, date=date_new, date_submitted)]

readr::write_csv(dat, "./auxiliary-data/raw-target-data/metadata.csv")
