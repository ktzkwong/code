# Papers with multiple diagnostics

## Editing Author: Kevin Kwong
## Date: 4/3/2017
## Purpose: Find all sources using more than 1 diagnostic in the full GAHI dataset

#########################################################################################################
rm(list=ls())
dev.off()

pacman::p_load(data.table, ggplot2, magrittr, reshape2)

workdir <- "H:/code/"
setwd(workdir)

# load datasets
gahi <- fread("gahi_3.28.csv")
extraction <- fread("geospatial_LF_cserway_Mar162017.csv")

# get unique prevalence fields in GAHI dataset
names <- names(gahi)
prev_names <- grep("prev", names, value = T)
prev_names <- grep("restr", prev_names, value = T)
prev_names <- prev_names[! prev_names %in% "prev_dx_restr"]

# convert unique prevalence field to binary (is NA or not)
for (col in prev_names){
  gahi[[col]] <- !is.na(gahi[[col]])
}

# get count of diagnostics per row
num_diagnostics <- apply(gahi[, .SD, .SDcol = prev_names], 1, function(x) length(which(x == T)))
gahi <- cbind(gahi, num_diagnostics)

# subset sources that used multiple diagnostics
multi_sources <- gahi[num_diagnostics > 1, .(Author_restr, Reference_restr, year_pub_restr)]
multi_sources <- unique(multi_sources)

write.csv(multi_sources, "sources_multi_diag.csv")

# get sources that use multiple diagnostics in extraction sheet
multi <- as.character()

for (source in unique(extraction$field_citation_value)) {
  temp <- extraction[field_citation_value == source,]
  l <- length(unique(temp$case_diagnostics)) 
  if (l > 1) {
    multi <- c(multi, source)
  } 
}

# get corresponding file paths and copy files to target directory
multi_subset <- extraction[field_citation_value %in% multi, .(field_citation_value, file_path, case_diagnostics)]
multi_subset <- multi_subset[field_citation_value != "",]
filepaths <- unique(multi_subset$file_path)
targetdir <- "J:/Project/NTDS/data/LF/multiple_diagnostics/"
file.copy(from= filepaths, to = targetdir, overwrite = F, recursive = F, copy.mode = T)

write.csv(multi_subset, "sources_multi_diagnosis_extraction.csv")
