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
                                    date = "character", ## enough "?"s that it can't be read directly as date
                                    date_submitted = "Date",
                                    host = "character"))

## filter by location
dat <- dat[division %in% c(state.name, "Washington DC","Puerto Rico"),]

## exclude any rows that have NA dates
dat <- dat[, date := as.Date(date, format = "%Y-%m-%d")]
dat <- dat[!is.na(date),]
dat <- dat[!is.na(date_submitted),]

## only include human specimens
dat <- dat[host == "Homo sapiens",]

## rename columns, and now that we have filtered by host, excluding that column
dat <- dat[, .(location = division, clade = Nextstrain_clade, date, date_submitted)]

readr::write_csv(dat, "./auxiliary-data/raw-target-data/metadata.csv")
