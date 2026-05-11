require(dplyr)
require(rjson)
require(arrow)
require(rstan)
# helper function to get the mean probabilities
mean_probs <- function(stan, num_days, num_clades){
  means <- extract(stan, pars = c("raw_alpha", "raw_beta")) # the alpha and beta
  days <- rep(0, num_clades) # the probabilities for each day
  coef <- means$raw_beta # the beta's
  dates <- c(1:num_days) # the days we will calculate the probabilities for
  intercepts <- means$raw_alpha # the alphas
  probs <- array( dim = c(num_clades, num_days, length(intercepts[1,]))) # the matrix of probabilities for each day for each of the samples
  mean_probs <- matrix( nrow = num_clades, ncol = num_days) # the mean probabilities over all the samples
      for(j in 1:length(intercepts[1, ])){
        for( i in 1:num_days){ # getting the probabilities
          days[1:(num_clades-1)] <- exp(intercepts[j, ] + coef[j,]*dates[i])/(sum(exp(intercepts[j, ] + coef[j, ]*dates[i]))+1) # calculating the probablities
          #for all but the reference
          days[num_clades] <- 1 - sum(days[1:(num_clades-1)]) # getting the probability for the reference
          probs[  , i, j] <- days
        }
      }
    for(i in 1:num_days){
      for(j in 1:num_clades){
        mean_probs[j, i] <- mean(probs[j, i, ]) # finding the mean over the sampled probabilities
      }
    }
  return(mean_probs)
}
create_baseline_submission <- function(target_date, num_days){
  if(target_date == as.Date("2024-10-09")){
    data_lag = 9
  } else{
    data_lag = 2
  }
  # reading in the data and summarizing it the way we need to
  data <- read_parquet(paste("https://covid-clade-counts.s3.amazonaws.com/",
                             as.Date(target_date - data_lag), "_covid_clade_counts.parquet", sep = ""))
  json_path <- "../auxiliary-data/modeled-clades/"
  clades <- fromJSON(file = paste0(json_path, target_date, ".json"))$clades
  data <- mutate(data, clade = ifelse(clade %in% clades, clade, "other"))
  data <- filter(data, date >= as.Date(target_date - num_days))
  data$days <- as.numeric(data$date - as.Date(target_date - num_days))
  data <- aggregate(count ~ days + clade, data = data, sum)
  data <- tidyr::complete(data = data, days, clade, fill = list(count=0))
  clade_totals <- aggregate(count ~ clade, data = data, sum)
  max_clade <- clade_totals$clade[which.max(clade_totals$count)]
  clades <- clades[-which(clades == max_clade)]
  clades[length(clades) + 1] <- max_clade
  data$clade <- as.numeric(factor(data$clade, levels = clades ))
  num_clades <- length(clades)
  mlr_data <- list(
    weights = data$count,
    y = data$clade,
    x = data$days, # the days the cases happened (first day of dataset treated as 0)
    N = length(data$clade), # how many cases we had
    K = num_clades # number of different clades
  )
  model  <- stan( # fitting the model
    file = "./Baseline.stan",
    data = mlr_data,
    chains = 1,
    refresh = 500
  )
  model_draws <- extract(model, pars = c("raw_alpha", "raw_beta"))
  values <- rep(NA, 42*100*num_clades)
  samples <- rep("sample", 42*100*num_clades*52)
  nowcast_date <- rep(target_date,42*100*num_clades*52 )
  clade_id <- rep(clades,42*100*52 )
  target_dates <- rep(NA,num_clades*42*100)
  state_ab <- c(state.abb, "PR", "DC")
  location <- rep(NA, 42*100*num_clades*52)
  mean_locations <- c(rep("0",52*num_clades*42 ))
  ids <- rep(NA, num_clades*100)
  date_seq <- as.Date((target_date - 31):(target_date + 10))
  output_type_id <- rep(NA, 42*100*num_clades*52)
  for(i in 1:52){
    location[(42*100*num_clades*(i-1) + 1):(i*(42*100*num_clades))] <- rep(state_ab[i],42*100*num_clades)
    mean_locations[(42*(i -1)*num_clades + 1):(i*(42*num_clades))] <- rep(state_ab[i],42*num_clades)
    for(w in 1:100){
      if(w < 10){
        ids[((w - 1)*num_clades + 1):(w*num_clades)] <- rep(paste0(state_ab[i], "0",(w - 1) ), num_clades)
      } else{
        ids[((w - 1)*num_clades + 1):(w*num_clades)] <- rep(paste0(state_ab[i],(w - 1) ), num_clades)
      }
    }
    output_type_id[((i-1)*num_clades*100*42 + 1):((i)*num_clades*100*42)] <- rep(ids, 42)
  }
  for(j in 1:42){
    target_dates[ (1 + 100*num_clades*(j-1)):(j*100*num_clades)] <- rep(date_seq[j], 100*num_clades)
  }
  all_target_dates <- rep(target_dates, 52)
  total_draws <- length(model_draws$raw_alpha[, 1])
  num <- runif(1, min = 1, max = total_draws/100)
  alphas <- matrix(nrow = num_clades - 1, ncol = 100)
  betas <- matrix(nrow = num_clades - 1, ncol = 100)
  for(m in 1:(num_clades - 1)){
    alphas[m, ] <- model_draws$raw_alpha[seq(num, total_draws - 20 + num, length.out = 100), m]
    betas[m, ] <- model_draws$raw_beta[seq(num, total_draws - 20 + num, length.out = 100 ), m]
  }
    for(j in 1:42){
      for(k in 1:100){
        sample_prop <- rep(NA, length(clades))
        sample_prop[1:(num_clades - 1)] <- exp(alphas[, k] + betas[, k]*(num_days - 32 + j))/(sum(exp(alphas[, k] + betas[, k]*(num_days - 32 + j))) + 1)
        sample_prop[num_clades] <- 1 - sum(sample_prop[1:(num_clades - 1)])
        values[( (j-1)*100*(length(clades)) + (k-1)*(length(clades)) +1):( (j-1)*100*(length(clades)) + (k)*(length(clades)))] <- sample_prop
      }
    }
  all_values <- rep(values, 52)
  #computing the means
  means <- mean_probs(stan = model, num_days = num_days + 10, num_clades = num_clades)
  mean_values <- c(rep(0, num_clades*42))
  mean_sample_ids <- c(rep(NA,52*num_clades*42 ))
  mean_output_type <- c(rep("mean",52*num_clades*42 ))
  mean_clade_ids <- c(rep(clades, 52*42))
  mean_origin_date <- c(rep(target_date,num_clades*42*52))
  mean_horizon <- rep(0,num_clades*42)
    for(i in 1:42){
      mean_values[(1 + num_clades*(i-1)):(i*num_clades)] <- means[, num_days - 32 + i]
      mean_horizon[(1 + num_clades*(i-1)):(i*num_clades) ] <- rep(date_seq[i], num_clades)
    }
  #converting output types
  mean_horizon <- as.Date(mean_horizon)
  all_mean_horizon <- rep(mean_horizon, 52)
  all_mean_values <- rep(mean_values, 52)
  # the full submission data frame
  submission_df <- data_frame(nowcast_date = c(nowcast_date,mean_origin_date), target_date = c(as.Date(all_target_dates),as.Date(all_mean_horizon)), clade = c(clade_id, mean_clade_ids), location = c(location, mean_locations), output_type = c(samples, mean_output_type), output_type_id = c(output_type_id, mean_sample_ids), value = c(all_values, all_mean_values))
  return(submission_df)
}
set.seed(787076)
target_date <- c(as.Date("2025-04-16"))
submission <- create_baseline_submission(target_date, num_days = 50)
write_parquet(submission, paste0("./model-output/Hub-baseline/", as.Date(target_date), "-Hub-baseline.parquet"))
hubValidations::validate_model_data(hub_path=".", file_path="UMass-HMLR/2025-04-16-UMass-HMLR.parquet")
plot_summary_graphs(model_output_file = "Hub-baseline/2025-04-16-Hub-baseline.parquet",
                    s3_data_date = "2025-04-14", hub_path = "./", save_path = "./")
set.seed(787076)
target_dates <- seq(as.Date("2024-10-16"),as.Date("2025-03-26"), by = 7 )
target_dates <- c(as.Date("2024-12-18"))
for(date in target_dates){
  submission <- create_baseline_submission(as.Date(date), num_days = 50)
  write_parquet(submission, paste0("../model-output/Variant_hub-baseline/", as.Date(date), "-Variant_hub-baseline.parquet"))
  print(as.Date(date))
}
