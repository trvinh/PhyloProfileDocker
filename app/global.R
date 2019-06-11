#' Startup script for PhyloProfile
#' 1) install and load packages
#' 2) start the PhyloProfile app

source("R/functions.R")

# List of dependent packages --------------------------------------------------
packages <- c("ape", "colourpicker", "data.table", "dendextend", "devtools",
              "plyr", "dplyr", "DT", "energy", "ggdendro", "ggplot2", "gplots",
              "grid", "gridExtra", "gtable", "RCurl", "RColorBrewer",
              "reshape2", "scales", "shiny", "shinyBS", "shinyjs", "stringr",
              "tidyr", "zoo")

# Set path for install packages while deploy into shiny server ----------------
# (from https://gist.github.com/wch/c3653fb39a00c63b33cf)

# By default, the directories in .libPaths() aren't writable on shinyapps.io
# so create a subdir where we'll install our package.
if (!file.exists("R-lib")) {
  dir.create("R-lib")
}

# Unfortunately, there's no way to get deployapp() to ignore this directory,
# so make sure to remove it locally before you call deployapp().
# This can be done with:
unlink("R-lib", recursive = TRUE)

# You may also need to restart R before calling deployapp(), because calling
# runApp() will modify your libpath (below), which can confuse deployapp().
# Add ./R-lib/ to the libPaths
.libPaths( c(normalizePath("R-lib/"), .libPaths()) )

# Find & install missing packages ---------------------------------------------
missing_packages <-
	packages[!(packages %in% installed.packages()[, "Package"])]
print(missing_packages)
# installPackages(packages)

# Load packages
lapply(packages, library, character.only = TRUE)

# Install packages from bioconductor ------------------------------------------
bioconductor_pkgs <- c("Biostrings", "bioDist")
# installPackagesBioconductor(bioconductor_pkgs)
lapply(bioconductor_pkgs, library, character.only = TRUE)

# Install OmaDB's dependencies "GO.db", "GenomeInfoDbData"
# and ExperimentHub
oma_pkgs <- c("GO.db", "GenomeInfoDbData", "ExperimentHub")
# installPackagesBioconductor(oma_pkgs)
if (packageVersion("ExperimentHub") < "1.11.1")
    BiocManager::install(pkgs = "ExperimentHub", version = "devel")
lapply(oma_pkgs, library, character.only = TRUE)

# Load demo data from PhyloProfileData package
eh = ExperimentHub()
myData <- query(eh, "PhyloProfileData")
