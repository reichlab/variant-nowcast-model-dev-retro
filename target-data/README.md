# target-data folder

This directory should contain target data.

One sub-directory entitled `timeseries` that will have a set of nested set of subdirectories.
Each subdirectory should be named something like `as_of=2024-11-18` and then inside that will live a parquet file that will have the target data. This would be the output of a call to `arrow::write_dataset()`, where we are partitioning on a new variable that we will create called "as_of".

We will download once the full-open metadata file at [the Nextstrain site](https://docs.nextstrain.org/projects/ncov/en/latest/reference/remote_inputs.html#remote-inputs-open-files) and create a new `as_of` variable based on the `date_submitted` column. 
