#' the script for creating the submission parquet files

# loading in the required functions and libraries
source("./src-model/func_for_model_creation.R")
source("./src-model/func_for_model_submission.R")
require(yaml)
require(arrow)
require(rjson)

#' submission_creator, creates the submission dataframes for a given model
#'
#' @param stan_file_name, the name of the stan file to be used, a string should have .Stan on the end
#'
#' @param model_name the name of the model to be used; needs to match the model in the model metadata section.
#'
#' @returns nothing; creates the submission files
submission_creator <- function(stan_file_name, model_name){
  #the mondays before the dates we want to model
  target_dates <- seq(from = as.Date("2022-08-01"), to = as.Date("2024-08-05"), by = 7)
  # the paths to the files we will need
  path_to_json  <- "./auxiliary-data/modeled-clades"
  path_to_data <- "./target-data/time-series"
  path_to_stan <- "./src-model"
  path_for_output <- "./model-output"
  # creating the folder for output if necessary
  if(!file.exists(path_for_output)) {
    dir.create(path_for_output)
  }
  for(date in target_dates){
    date <- as.Date(date)
    json_file_name <- file.path(path_to_json, paste0(as.Date(date + 2), ".json"))
    # reading in the clades for the week
    clades <- fromJSON(file = json_file_name)
    data_file_name <- file.path(path_to_data, paste0("as_of=", date), "part-0.parquet")
    # reading in the data for the week
    data <- read_parquet(data_file_name)
    # only keeping the clades we want to model
    trimed_data <- trim_clades(data, clades = clades$clades)
    stan_file <- file.path(path_to_stan, stan_file_name)
    # the model for the week
    fitted_model <- stan_maker(trimed_data, stan_file = stan_file, target_date = as.Date(date + 2), num_seq = 1, num_days = 150,
                               interations = 15000, warmup = 7000 )
    set.seed(1)
    # the submission for the week
    submission_df <- prediction_sampler(fitted_model, given_date = as.Date(date + 2))
    file_name <- paste0(paste(as.Date(date + 2), "UMass", model_name, sep = "-"), ".parquet")
    # writing out the submission
    write_parquet(submission_df, file.path(path_for_output, file_name))
  }
}
