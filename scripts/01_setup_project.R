# Create the directory structure expected by the analysis scripts.
dirs <- c(
  "data_raw",
  "data_processed",
  "metadata",
  "results",
  "results/validation",
  "figures",
  "supplementary"
)
for (d in dirs) dir.create(d, recursive = TRUE, showWarnings = FALSE)
message("Project directories created.")
