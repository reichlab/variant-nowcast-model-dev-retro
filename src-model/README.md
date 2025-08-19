# src-model

This directory should contain code for models. Every model represented in model-output and model-metadata should have a corresponding source file.

To create a model using this repo and these scripts, follow the following steps:

1. Clone the repo and set the root of the project as the working directory.
2. Source `func_for_model_creation` in the `src-model` folder. This will load the functions and libraries needed for the model creation. 
3. Choose what dataset you want to use; the datasets are in the target-data folder and are indexed by Mondays between 2022-08-01 and 2024-08-05.
4. Find the clades modeled for the week. These are in `auxiliary-data/modeled-clades`. These are indexed by Wednesday, so if you choose 2022-08-01 for the dataset you want 2022-08-03 for the clades.
5. Use `trim_clades` to make the clades in the dataset match the clades modeled. If you have loaded in the data as data and clades as clades you would use the following syntax:
   ```r
   data <- trim_clades(data = data, clades = clades)
   ```
6. Fit the model using `stan_maker` or `stan_maker_dirichlet`; the following example would fit a TDLMN model to 2022-08-03 (you want to use the Wednesday from above) to all states with at least one sequence in the last 60 days, using the last 150 days of data, for 3000 total iterations and 1000 warmup iterations:
   ```r
   fitted_model <- stan_maker( data = data,
                               stan_file = "./src-model/HMLR_TDLMN.stan",
                               num_seq = 1,
                               target_date = as.Date("2022-08-03"),
                               num_days = 150,
                               iterations = 3000,
                               warmup = 1000)
   ```
   The function will return a list with the Stan draws along with some information that will differ depending on the exact Stan model fit. See the description of `stan_maker` in `func_for_model_creation` for more details.
7. You can then source `func_for_model_submission` to use functions that create the submission data files and find the mean probabilities. 
