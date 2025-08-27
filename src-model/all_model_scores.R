library("tidyverse")
library("dplyr")
library("ggplot2")
library("arrow")
library("scoringRules")
require(rjson)
source("./src-model/func_for_model_creation.R")
source("./src-model/func_for_model_submission.R")
set.seed(1989)
get_model_scores <- function(df_model_output, df_validation){
  # getting the order of the clades
  clades <- unique(df_model_output$clade)
  locs <- sort(unique(df_model_output$location))
  # finding the dates and subsetting
  date_min <- min(df_model_output$target_date)
  date_max <- max(df_model_output$target_date)
  dates <- seq.Date(from = date_min, to = date_max, by = "day")

  targets <- df_validation

  targets <- targets |> subset((date %in% dates) &
                                 (location %in% locs) &
                                 (clade %in% clades))
  #getting all the columns needed
  columns <- rep("", 4 + 2*length(clades))
  columns[1:4] <- c("es_score","brier_score", "location", "target_date")
  PI_levels <- c("50", "95")
  for(j in 1:length(PI_levels)){
    for(i in 1:length(clades)){
      columns[4 + i + (j-1)*length(clades)] <- paste0(clades[i],"_",PI_levels[j],"_CI")
    }
  }
  # creating the scores data frame
  df_scores <- data.frame(matrix(nrow = 0, ncol = length(columns)))
  colnames(df_scores) <- columns
  for(loc in locs){
    print(loc)
    for(day in dates){
      #  if there is no data on the day, if not skip to next
      if(sum(filter(df_validation, location == loc, date == as.Date(day))$observation) == 0){
        x = list(NA,NA, loc, as.Date(day))
        for(w in 5:(4 + 2*length(clades))){
          x[w] <- NA
        }
        df_temp <- as.data.frame(x, col.names = columns)
        df_scores <- rbind(df_scores, df_temp)

        next
      }

      es <- as.numeric()

      # Observed proportions
      df_obs <- filter(targets, location == loc, date == as.Date(day))
      obs_count <- rep(NA, length(clades))
      for(i in 1:length(clades)){
        obs_count[i] <- sum(filter(df_obs, clade == clades[i])$observation)
      }
      # MCMC sample of modeled COUNTS
      df_samp <- subset(df_model_output, target_date == as.Date(day) & location == loc & output_type == "sample") |>
        group_by(clade)

      # Pivot wider to get to MCMC format for scoring
      df_samp_wide <- pivot_wider(df_samp, names_from = output_type_id, values_from = value) |>
        group_by(clade)
      df_samp_wide <- subset(df_samp_wide,
                             select = -c(nowcast_date, target_date,
                                         clade, location, output_type))
      # get df for means to calculate brier scores
      df_mean <- subset(df_model_output, target_date == as.Date(day) & location == loc & output_type == "mean") |>
        group_by(clade)
      df_mean <- subset(df_mean,
                        select = -c(nowcast_date, target_date,
                                    location, output_type))
      # convert samples to matrix for scoringRules
      samp_matrix <- as.matrix(df_samp_wide)

      # create matrix to store the multinomial sample
      samp_counts <- matrix(, nrow = length(clades), ncol = 10000)
      j = 1
      ## Generate 100 multinomial observations for samp_props
      # Need the N for each loc/day from the validation data
      N <- sum(subset(targets,
                      location == loc &  date == as.Date(day))$observation)
      cov_mat_50 <- rep(NA,length(clades))
      cov_mat_95 <- rep(NA,length(clades))
      for(col in 1:dim(samp_matrix)[2]){

        # Get sample clade proportions from predictive distr
        samp_props <- samp_matrix[,col]

        # 100 multinomial samples
        samp_counts[, (100*(j-1) + 1):(100*j)] <- rmultinom(n = 100, size = N, prob = samp_props)
        j = j + 1
      }
      # calculating whether prediction intervals cover
      for(m in 1:length(clades)){
        cont <- quantile(samp_counts[m,], probs = c(0.25, 0.75))
        if( cont[1] <= obs_count[m] && cont[2] >= obs_count[m] ){
          cov_mat_50[m] <- 1
        } else{
          cov_mat_50[m] <- 0
        }
        cont <- quantile(samp_counts[m], probs = c(0.025, 0.975))
        if( cont[1] <= obs_count[m] && cont[2] >= obs_count[m] ){
          cov_mat_95[m] <- 1
        } else{
          cov_mat_95[m] <- 0
        }
      }
      # Energy score for the 10000 multinomial samples for day/loc
      es <- es_sample(y = obs_count, dat = samp_counts)

      #brier scores for mean predictions
      brier_score <- 0
      for(k in 1:length(clades)){
        if(clades[k] != df_mean$clade[k]){
          print("failed")
        }
        brier_score <- brier_score + obs_count[k]*(df_mean$value[k] - 1)^2 + (N- obs_count[k])*(df_mean$value[k])^2
      }
      brier_score <- brier_score/(N)
      #standardizing to normal brier score
      brier_score <- brier_score/2
      # Store scores and coverage to data frame
      CIcov <- c(cov_mat_50,cov_mat_95)
      x = list(es,brier_score, loc, as.Date(day))
      for(w in 5:(4 + 2*length(clades))){
        x[w] <- CIcov[w- 4]
      }
      df_temp <- as.data.frame(x, col.names = columns)
      df_scores <- rbind(df_scores, df_temp)

    }
  }
  return(df_scores)
}
target_dates <- seq(from = as.Date("2022-08-01"), to = as.Date("2024-08-05"), by = 7)
#models <- c("FDLMN","TDLDN","TDLMN","FDLDN", "FDSDN", "FDSMN", "TDSDN", "TDSMN", "baseline")
models <- c("FDLMN")
path_to_json  <- "./auxiliary-data/modeled-clades"
path_to_data <- "./target-data/time-series/as_of="
path_to_model <- "./model-output"
path_for_output <- "./model-scores"
int_number <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
date <- target_dates[int_number]
print(date)
for(model in models){
    json_file_name <- file.path(path_to_json, paste0(as.Date(date + 2), ".json"))
    clades <- fromJSON(file = json_file_name)
    if(model == "baseline"){
      submission_file <-read_parquet(paste0(path_to_model,"/Hub-",model,"/",as.Date(date + 2), "-Hub-", model,".parquet"))
    } else{
     submission_file <-read_parquet(paste0(path_to_model,"/UMass-HMLR_",model,"/",as.Date(date + 2), "-UMass-HMLR_", model,".parquet"))
    }
    data_to_score_on <- read_parquet(paste0(path_to_data,as.Date(date + 91),"/part-0.parquet"))
    data_to_score_on <- trim_clades(data_to_score_on, clades$clades)
    data_to_score_on$location <- convert_to_abbreviation(data_to_score_on$location)
    df_scores <- get_model_scores(submission_file, data_to_score_on)
    if( model == "baseline"){
      output_direct = paste0(path_for_output,"/Hub-", model)
    } else{
      output_direct <- paste0(path_for_output,"/UMass-HMLR_", model)
    }
    if(!file.exists(output_direct)) {
      dir.create(output_direct, recursive = T)
    }
    if( model == "baseline"){
      write.csv(df_scores, file = paste0(output_direct,"/",as.Date(date + 2),"-Hub-",model,".csv"))
    } else{
      write.csv(df_scores, file = paste0(output_direct,"/",as.Date(date + 2),"-UMass-HMLR_",model,".csv"))
    }
    print(model)
}
