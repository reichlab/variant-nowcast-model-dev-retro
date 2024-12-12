#' All the functions required for creating the model submissions
require(splines)

#' creates the data_frame that will be turned into the submission
#' @param stan the list returned by stan_maker or stan_maker_splines
#' @param given_date what day was the model fit to, a date object
#' @param N The number of different trajectories to draw, default is 1000
#' @param dates which days to model, defaults to 31 days before to 10 after the given date
#' @param splines Is the model the spline model, default false
#' @returns Returns a DF containing  all the information required required for the hub submission
prediction_sampler <- function(stan, given_date, N = 1000, dates = c(119:160), splines = FALSE){
  K <- stan$K # the number of clades
  L <- stan$L # the number of locations modleed
  target_lo <- convert_to_abbreviation(stan$target_lo) # mapping states to there two-letter form
  sample_ids <- c(rep("0", N*length(target_lo)*length(dates)*length(stan$clades))) # getting the right form for the ids
  clade_ids <- c(rep(stan$clades, N*length(target_lo)*length(dates))) # the ids for each clade
  origin_date <- c(rep(given_date,N*length(target_lo)*length(dates)*length(stan$clades))) # the date the forecast was created
  location <- c(rep("0",N*length(target_lo)*length(dates)*length(stan$clades) )) # the locations
  mean_locations <- c(rep("0",L*K*length(dates) ))
  horizon <- c(rep(0, N*length(target_lo)*length(dates)*length(stan$clades))) # the day we are forecasting or nowcasting
  output_type <- c(rep("sample",N*length(target_lo)*length(dates)*length(stan$clades) ))
  values <- c(rep(0, N*length(target_lo)*length(dates)*length(stan$clades))) # the samples
  temp <- c(rep(0, length(stan$clades))) # the samples for a given day
  draws <- extract(stan$mlr_fit) # extracting the MCMC samples
  if(splines){
    random_draws <- array(dim = c(K-1,1 + stan$B,N))
    spline <- bs(1:160, degree =  stan$B)
  } else{
    random_draws <- array(dim = c(K-1,2,N))
  }
  # indexing the right location dates
  for(i in 1:length(target_lo)){
    location[(1 + (i-1)*(N*length(stan$clades)*length(dates))):( (i)*(N*length(stan$clades)*length(dates)))] <- rep(target_lo[i],N*length(stan$clades)*length(dates))
    mean_locations[(1 + (i-1)*K*length(dates)):((i)*K*length(dates))] <- rep(target_lo[i],length(stan$clades)*length(dates))
  }
  # creating the sample ids for each location
  for(l in 1:L){
    for(i in 1:(length(dates))){
      for(m in 1:N){
        if( m-1 < 10){
          word <- rep(paste0(target_lo[l],"0",m-1), K)
        } else{
          word <- rep(paste0(target_lo[l],m-1), K)
        }
        sample_ids[(1 + (m-1)*K + (i-1)*N*(K) + (l-1)*(N)*(K)*(length(dates))):((m)*K + (i-1)*N*K + (l-1)*(N)*(K)*(length(dates)))] <- word
      }
    }
  }
  # getting the sample predictions
  if(splines){
    for(l in 1:L){
      for(n in 1:N){
        c <- ceiling(runif(1, min = 0, max = length(draws$raw_alpha[ ,1, 1])/N))
        for(q in 1:(K-1)){
          random_draws[q, 1,n] <- draws$raw_alpha[c + (n-1)*length(draws$raw_alpha[ ,1, 1])/N,l,q] # getting the random draws
          for(b in 1:stan$B){
            random_draws[q, b + 1,n] <- draws$raw_beta[c + (n-1)*length(draws$raw_alpha[ ,1, 1])/N,l,q, b]
          }
        }
      }
      for(i in 1:length(dates)){
        for(m in 1:N){
          for(k in 1:(K-1)){
            temp[k] <- exp(random_draws[k, 1, m] + sum(random_draws[k, 2:(stan$B + 1), m]*spline[dates[i],]))
          }
          temp[1:(K-1)] <- temp[1:(K-1)]/(sum(temp[1:(K-1)]) + 1)
          temp[K] <- 1 - sum(temp[1:(K-1)])
          values[ (1 + (m-1)*K + (i-1)*N*(K) + (l-1)*(N)*(K)*(length(dates))):((m)*K + (i-1)*N*K + (l-1)*(N)*(K)*(length(dates)))  ] <- temp
          horizon[ (1 + (m-1)*K + (i-1)*N*(K) + (l-1)*(N)*(K)*(length(dates))):((m)*K + (i-1)*N*K + (l-1)*(N)*(K)*(length(dates)))  ] <- rep(given_date + i - length(dates) + 10, K)
        }
      }
    }
  } else{
    for(l in 1:L){
      for(n in 1:N){
        c <- ceiling(runif(1, min = 0, max = length(draws$raw_alpha[ ,1, 1])/N))
        for(q in 1:(K-1)){
          random_draws[q, ,n] <- c(draws$raw_alpha[c + (n-1)*length(draws$raw_alpha[ ,1, 1])/N,l,q], draws$raw_beta[c + (n-1)*length(draws$raw_alpha[ ,1, 1])/N,l,q] ) # getting the random draws
        }
      }
      for(i in 1:length(dates)){
        for(m in 1:N){
          temp[1:(K-1)] <- exp(random_draws[, 1, m] + random_draws[, 2, m]*dates[i])/(sum(exp(random_draws[, 1, m] + random_draws[, 2, m]*dates[i]))+1)
          temp[K] <- 1 - sum(temp[1:(K-1)])
          values[ (1 + (m-1)*K + (i-1)*N*(K) + (l-1)*(N)*(K)*(length(dates))):((m)*K + (i-1)*N*K + (l-1)*(N)*(K)*(length(dates)))  ] <- temp
          horizon[ (1 + (m-1)*K + (i-1)*N*(K) + (l-1)*(N)*(K)*(length(dates))):((m)*K + (i-1)*N*K + (l-1)*(N)*(K)*(length(dates)))  ] <- rep(given_date + i - length(dates) + 10, K)
        }
      }
    }
  }
  horizon <- as.Date(horizon)
  # getting the mean probabilities
  if(splines){
    means <- mlr_probs_splines(stan = stan, num_days = max(dates))
  } else{
    means <- mlr_probs(stan = stan, num_days = max(dates))
  }
  mean_values <- c(rep(0, L*K*(length(dates))))
  mean_sample_ids <- c(rep(NA,L*K*length(dates) ))
  mean_output_type <- c(rep("mean",L*K*length(dates) ))
  mean_clade_ids <- c(rep(stan$clades, length(target_lo)*length(dates)))
  mean_origin_date <- c(rep(given_date,K*length(dates)*L))
  mean_horizon <- rep(0 , K*L*length(dates))
  for(l in 1:L){
    for(i in 1:length(dates)){
      mean_values[(1 + (i-1)*K + (l-1)*(length(dates))*K):( (i)*K + (l-1)*(length(dates))*K) ] <- means[[stan$target_lo[l]]][, dates[i]]
      mean_horizon[(1 + (i-1)*K + (l-1)*(length(dates))*K):( (i)*K + (l-1)*(length(dates))*K) ] <- rep(given_date + i - length(dates) + 10, K)
    }
  }
  #converting output types
  mean_horizon <- as.Date(mean_horizon)
  clade_ids <- as.character(clade_ids)
  mean_clade_ids <- as.character(mean_clade_ids)
  # creating the dataframe
  df <- data.frame(nowcast_date = c(origin_date, mean_origin_date), target_date = c(horizon, mean_horizon), clade = c(clade_ids, mean_clade_ids), location = c(location, mean_locations), output_type = c(output_type, mean_output_type) , output_type_id = c(sample_ids,mean_sample_ids), value = c(values, mean_values))
  return(df)
}

#' returns a list of probabilities for all locations, called in prediction_sampler
#'
#' @param stan the list returned by stan_maker
#' @param num_days the number of days of probability wanted, counts up from the first day that the model was fit to
#' @param shifted if T uses then uses time since the variant was introduced instead of days from beginning of dataset
#' should be set to F at present.
#' @returns a named list containing mean probablities for each location, indexed by location

mlr_probs <- function(stan, num_days, shifted = F){
  full_probs <- list()
  means <- extract(stan$mlr_fit, pars = c("raw_alpha", "raw_beta")) # the alpha and beta
  # the number of location
  L <- stan$L
  # the number of clades
  if(!is.null(stan$K)){
    K <- stan$K
  } else{
    K <- stan$V
  }
  days <- rep(0, K) # the probablities for each day
  for(l in 1:L){
    intercepts <- means$raw_alpha[, l, ] # the alphas
    coef <- means$raw_beta[ , l , ] # the beta's
    probs <- array( dim = c(K, num_days, length(intercepts[1,]))) # the matrix of probabilities for each day for each of the samples
    mean_probs <- matrix( nrow = K, ncol = num_days) # the mean probabilities over all the samples
    dates <- c(1:num_days) # the days we will calculate the probabilities for
    if(!(shifted)){
      for(j in 1:length(intercepts[1, ])){
        for( i in dates){ # getting the probabilities
          days[1:(K-1)] <- exp(intercepts[j, ] + coef[j,]*dates[i])/(sum(exp(intercepts[j, ] + coef[j, ]*dates[i]))+1) # calculating the probablities
          #for all but the reference
          days[K] <- 1 - sum(days[1:(K-1)]) # getting the probability for the reference
          probs[  , i, j] <- days
        }
      }
    } else {
      for( j in 1:length(intercepts[1, ])){
        for( i in dates){ # getting the probabilities
          days[1:(K-1)] <- exp(intercepts[j, ] + coef[j, ]*(dates[i] - stan$start_times))
          for( m in 1:(K-1)){
            if(stan$start_times[m] >= i ){
              days[m] <- 0.001
            }
          }
          days[1:(K-1)] <- days[1:(K-1)]/(sum(days[1:(K-1)]) + 1)
          days[K] <- 1 - sum(days[1:(K-1)])
          probs[, i, j] <- days
        }
      }
    }
    for(i in 1:num_days){
      for(j in 1:K){
        mean_probs[j, i] <- mean(probs[j, i, ]) # finding the mean over the sampled probabilities
      }
    }
    full_probs[[stan$target_lo[l]]] <-mean_probs # saving the probabilities for the location
  }
  for(lo in stan$target_lo){
    row.names(full_probs[[lo]]) <- stan$clades # indexing the probablities by clade
  }
  return(full_probs)
}

#' returns a list of probabilities for all locations used for spline model, called in prediction_sampler
#'
#' @param stan the list returned by stan_maker_splines
#' @param num_days the number of days of probability wanted, counts up from the first day that the model was fit to
#' #' @returns a named list containing mean probablities for each location, indexed by location
mlr_probs_splines <- function(stan, num_days){
  full_probs <- list()
  means <- extract(stan$mlr_fit, pars = c("raw_alpha", "raw_beta")) # the alpha and beta
  L = stan$L
  B = stan$B
  if(!is.null(stan$K)){
    K = stan$K
  } else{
    K = stan$V
  }
  for(l in 1:L){
    intercepts <- means$raw_alpha[, l, ] # the alphas
    coef <- means$raw_beta[ , l , , ] # the beta's
    probs <- array( dim = c(K, num_days, length(intercepts[1,]))) # the matrix of probabilities for each day
    mean_probs <- matrix( nrow = K, ncol = num_days)
    dates <- bs(c(1:num_days), degree =  B)
    days <- c(rep(0, K))
    for(j in 1:length(intercepts[1, ])){
      for( i in 1:num_days){ # getting the probabilities
        for(k in 1:K-1){
          days[k] <- exp(intercepts[j, k] + sum(coef[j,k,]*dates[i, ]))
        }
        days[1:K-1] <- days[1:K-1]/(sum(days[1:K-1]) +1)
        days[K] <- 1 - sum(days[1:(K-1)])
        probs[  , i, j] <- days
      }
    }
    for(i in 1:num_days){
      for(j in 1:K){
        mean_probs[j, i] <- mean(probs[j, i, ])
      }
    }
    full_probs[[stan$target_lo[l]]] <- mean_probs
  }
  for(lo in stan$target_lo){
    row.names(full_probs[[lo]]) <- stan$clades
  }
  return(full_probs)
}

#' Function to convert state names to abbreviations, Chat GPT created function, called in prediction_sampler
#'
#' @param states the state names we went to convert to abbreviations
#' @returns the abbreviations of the state names
convert_to_abbreviation <- function(states) {
  state_abbreviation_map <- setNames(c(state.abb,"PR","DC"), c(state.name, "Puerto Rico", "Washington DC"))
  sapply(states, function(state) state_abbreviation_map[[state]])
}
