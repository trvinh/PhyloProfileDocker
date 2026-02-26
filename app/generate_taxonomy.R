library(PhyloProfile)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

# If an argument is provided, use it as data_dir; otherwise fallback to default
data_dir <- if (length(args) >= 1) args[1] else "/srv/shiny-server/data"

cat("Writing preProcessedTaxonomy.txt to:", data_dir, "\n")

# Generate taxonomy
preProcessedTaxonomy <- processNcbiTaxonomy()

# Write file into package data directory
write.table(
  preProcessedTaxonomy,
  file = file.path(data_dir, "preProcessedTaxonomy.txt"),
  col.names = TRUE,
  row.names = FALSE,
  quote = FALSE,
  sep = "\t"
)
