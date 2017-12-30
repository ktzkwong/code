# Multiple Diagnostics/ MF - ICT crosswalk

## Editing Author: Kevin Kwong
## Date: 4/3/2017
## Purpose: Identify all sources using more than 1 diagnostic in the full GAHI and lit extraction datasets. 
## Create a MF - ICT crosswalk dataset (same pop tested). Vet data quality issues. Preliminary exploration of crosswalk.

#########################################################################################################
rm(list=ls())
dev.off()

pacman::p_load(data.table, ggplot2, magrittr, reshape2, stringi, dplyr, ResourceSelection, glmmML, lme4, scales)

workdir <- "H:/code/"
setwd(workdir)

# load datasets
gahi_clean <- fread("gahi_3.28.csv") #KK- Liz said to not use this one right now
extraction <- fread("LF_Extraction_with_NIDs_cleaned_Diagnostics.csv")
#gahi_KK <- fread("gahi_KK_collapse.csv")
gahi_KK <- fread("gahi_KK_collapse_new.csv")
cw_analysis <- fread("cw_analysis.csv")

# source shared functions
central_dir <- 'J:/temp/central_comp/libraries/current/r/'
source(paste0(central_dir, 'get_location_metadata.R'))

#get locations for later
loc <- get_location_metadata(location_set_id = 2) %>% as.data.table
loc_short <- loc[,.(location_id, super_region_id, super_region_name, region_id, region_name, developed)]
loc_short$location_id <- as.integer(loc_short$location_id)

#get GAHI loc file
gahi_loc_dict <- fread("GAHI_location_dictionary.csv")
gahi_loc_dict <- gahi_loc_dict[,.(location_name, super_region_id, super_region_name, region_id, region_name, developed)]

psum <- function(..., na.rm=T) { 
  x <- list(...)
  rowSums(matrix(unlist(x), ncol=length(x)), na.rm=na.rm)
} 

### GAHI ######################################################################
if(FALSE){ #not using this per Liz
  gahi <- copy(gahi_clean)
  
  # fix missing restr prevalences for mf and ict
  gahi[is.na(prev_ict_restr) == T, prev_ict_restr := prev_ict]
  gahi[is.na(prev_mf_filtration_restr) == T, prev_mf_filtration_restr := prev_mf_filtration]
  gahi[is.na(prev_mf_restr) == T, prev_mf_restr := prev_mf]
  
  ## Get diagnostic count per row
  
  # get unique prevalence fields in GAHI dataset
  names <- names(gahi)
  prev_names_all <- grep("prev", names, value = T)
  prev_names_restr <- grep("restr", prev_names_all, value = T)
  prev_names_not_restr <- prev_names_all[!(prev_names_all %in% prev_names_restr)]
  prev_names_restr <- prev_names_restr[! prev_names_restr %in% "prev_dx_restr"]
  prev_names_not_restr <- prev_names_not_restr[! prev_names_not_restr %in% "prev_dx"]
  
  
  # convert unique prevalence field to binary (is NA or not)
  for (col in prev_names_all){
    gahi[[col]] <- !is.na(gahi[[col]])
  }
  
  # get count of diagnostics per row both for restr columns and not restr columns
  num_diagnostics_restr <- apply(gahi[, .SD, .SDcol = prev_names_restr], 1, function(x) length(which(x == T)))
  gahi <- cbind(gahi, num_diagnostics_restr)
  
  num_diagnostics_not_restr <- apply(gahi[, .SD, .SDcol = prev_names_not_restr], 1, function(x) length(which(x == T)))
  gahi <- cbind(gahi, num_diagnostics_not_restr)
  
  # identify data points where collapse code failed and restr does not match not restr
  # output if for if we need to crosswalk any other diagnostics
  problem_rows <- gahi[num_diagnostics_restr != num_diagnostics_not_restr]
  problem_rows <- problem_rows[num_diagnostics_not_restr != 0 & num_diagnostics_not_restr > num_diagnostics_restr &
                                 (num_diagnostics_not_restr > 1 | num_diagnostics_restr > 1)]
  
  problem_rows <- gahi_clean[uid_from_source_restr %in% problem_rows$uid_from_source_restr]
  
  #used previously to identify and check rows with problems in ict/mf prev fields
  #problem_rows <- problem_rows[(is.na(prev_ict_restr) == F | is.na(prev_ict) == F) & (is.na(prev_mf_restr) == F | is.na(prev_mf) == F | is.na(prev_mf_filtration_restr) == F | is.na(prev_mf_filtration_restr) == F)]
  #output csv to use for tracking down GAHI sources with multiple diagnostics
  #multi_sources <- multi_sources[num_diagnostics_restr > 1,]
  #multi_sources <- multi_sources[, .(Author_restr, Reference_restr, year_pub_restr)] %>% unique
  #write.csv(multi_sources, "sources_multi_diag.csv")
  
  # get data points that have both ICT and MF prev measurements
  
  gahi_ict_mf <- gahi_clean[, .(uid_from_source_restr, prev_mf_filtration_restr, prev_mf_restr, prev_ict_restr, 
                                pop_mf_filtration_restr, np_mf_filtration_restr, pop_mf_restr, np_mf_restr,
                                pop_ict_restr, np_ict_restr)]
  gahi_ict_mf <- melt(gahi_ict_mf, id.vars = "uid_from_source_restr", na.rm = T) %>% unique
  gahi_ict_mf <- dcast(gahi_ict_mf, uid_from_source_restr ~ variable, value.var = "value", fun.aggregrate = mean) %>% data.table
  
  gahi_ict_mf[is.na(prev_mf_filtration_restr) == F & is.na(prev_mf_restr) == F, 
              prev_mf_restr := (np_mf_filtration_restr + np_mf_restr)/(pop_mf_restr + pop_mf_filtration_restr)]
  
  gahi_ict_mf[is.na(prev_mf_filtration_restr) == F & is.na(prev_mf_restr) == T, 
              prev_mf_restr := prev_mf_filtration_restr]
  gahi_ict_mf <- gahi_ict_mf[is.na(prev_mf_restr) == F & is.na(prev_ict_restr) == F,]
  gahi_ict_mf[, sample_size_ict := pop_ict_restr]
  gahi_ict_mf[, sample_size_mf := rowSums(.SD, na.rm = T), .SDcols = c("pop_mf_restr", "pop_mf_filtration_restr")]
  gahi_ict_mf[, source := "GAHI"]
  gahi_ict_mf <- gahi_ict_mf[, .(uid_from_source_restr, prev_ict_restr, prev_mf_restr, source, sample_size_ict, sample_size_mf)]
  setnames(gahi_ict_mf, names(gahi_ict_mf), c("pop_id", "ict", "mf", "source", "sample_size_ict", "sample_size_mf"))
}

### Using KK collapsed GAHI dataset
# get data points with ict and mf measurements
ict_KK <- grep("ICT", gahi_KK$diagnostics)
mf_KK <- c(grep("Filtration", gahi_KK$diagnostics), grep("Blood", gahi_KK$diagnostics)) %>% unique
mf_ict_KK <- c(ict_KK, mf_KK) %>% table %>% as.data.table
mf_ict_KK <- mf_ict_KK[N > 1, .] %>% as.integer
gahi_KK_sub <- gahi_KK[mf_ict_KK, .(pop_id, Record_ID, np_ICT, np_Blood_smear, np_Filtration, N_rounds, ADM0, ADM1, ADM2, TypeSite, Site_name,
                                    Year_start, Year_end, Age_start, Age_end, Age_Range, Sex, pop_ICT, 
                                    pop_Blood_smear, pop_Filtration, prev_Blood_smear, prev_Filtration, prev_ICT)]

gahi_KK_sub[is.na(prev_Blood_smear) == F, c("prev_mf", "pop_mf", "np_mf", "filtration") := list(prev_Blood_smear, pop_Blood_smear, np_Blood_smear, 0)]
gahi_KK_sub[is.na(prev_Blood_smear) == T, c("pop_mf", "np_mf", "filtration") := list(pop_Filtration, np_Filtration, 1)]
gahi_KK_sub[,source := "GAHI"]
gahi_KK_sub[N_rounds == 0, mda := 0]
gahi_KK_sub[N_rounds > 0, mda := 1]
setnames(gahi_KK_sub, c("ADM0", "Year_start", "Year_end", "prev_mf", "prev_ICT", "pop_ICT", "pop_mf", "Sex", "Age_start", "Age_end"), 
         c("location_name", "year_start", "year_end", "mf", "ict", "sample_size_ict", "sample_size_mf", "sex", "age_start", "age_end"))
gahi_KK_sub[, c("np_ICT", "np_Blood_smear", "pop_Filtration", "pop_Blood_smear", "np_Filtration", "prev_Blood_smear", "prev_Filtration") := NULL]
gahi_KK_sub <- merge(gahi_KK_sub, gahi_loc_dict, by = "location_name")

### Lit Extraction from Chelsey ########################################################

# remove 2 rows missing file_path and citation; will probably be fixed by Chelsey. Check!!!
multi_subset <- extraction[file_path != "" & is.na(file_path) == F,]

# standardize diagnostics into buckets
ict <- grep("ICT", multi_subset$case_diagnostics)
mf <- grep("microf", multi_subset$case_diagnostics)
elisa <- grep("ELISA", multi_subset$case_diagnostics)
igg4 <- grep("IgG4", multi_subset$case_diagnostics)
clinical <- grep("clinical", multi_subset$case_diagnostics)
fts <- grep("FTS", multi_subset$case_diagnostics)
diagnostics <- c("ict", "mf", "elisa", "igg4", "clinical", "fts")

i <- 1
for (diag in diagnostics) {
  multi_subset[get(diag), diagnostic := i]
  i <- i + 1
}

diagnostic_n <- data.table(diagnostic = c(1:6), diagnostic_name = diagnostics)
multi_subset <- merge(multi_subset, diagnostic_n, by = "diagnostic")

# denote if mf is filtered
multi_subset[grep("filtered", case_diagnostics), filtration := 1]
multi_subset[is.na(filtration) == T, filtration := 0]

# get data points with ICT and MF measurements
ict_mf_extract <- c()
for (g in unique(multi_subset$group)) {
  diag <- multi_subset[group == g, diagnostic] %>% unique %>% sort
  if ((1 %in% diag) == T & (2 %in% diag) == T) {
    ict_mf_extract <- c(ict_mf_extract, g)
  }
}

ict_mf_extract <- multi_subset[group %in% ict_mf_extract & (diagnostic == 1 | diagnostic == 2)]

if(FALSE){# move pdfs for review
  multi_subset <- ict_mf_extract[, .(field_citation, file_path, case_diagnostics)]
  filepaths <- unique(multi_subset$file_path)
  targetdir <- "J:/Project/NTDS/data/LF/multiple_diagnostics/"
  file.copy(from= filepaths, to = targetdir, overwrite = F, recursive = F, copy.mode = T)
  
  # create reference csv
  write.csv(multi_subset, "sources_multi_diagnosis_extraction.csv")
}

### Collapse MF and ICT measurements in extraction data ####################################

# generate unique ids for same loc + time + pop; first split by specificity to prevent double counting and allow for choice to use age/sex/split data

# need to first split by specificity
i <- 1
ref_list <- list()
id_list <- list()
extract_list <- list()
for (spec in unique(ict_mf_extract$specificity)) {#discarding total specificity?
  ref_list[[i]] <- ict_mf_extract[specificity == spec,] %>% as.data.table
  
  id_list[[i]] <- ref_list[[i]][, .N, by = list(field_citation, file_path, location_name, lat, long, 
                                                poly_id, poly_reference, sex, year_start, year_end, 
                                                month_start, month_end, age_start, age_end, site_memo)] %>% as.data.table
  id_list[[i]] <- id_list[[i]][N > 1,]
  id_list[[i]]$N <- NULL
  pop_id <- 1:nrow(id_list[[i]])
  id_list[[i]] <- cbind(id_list[[i]], pop_id)
  
  
  ref_list[[i]] <- merge(ref_list[[i]], id_list[[i]], by = c("field_citation", "file_path", "location_name", "lat", "long", 
                                                             "poly_id", "poly_reference", "sex", "year_start", "year_end", 
                                                             "month_start", "month_end", "age_start", "age_end", "site_memo"), all.x = T)
  ref_list[[i]] <- ref_list[[i]][is.na(pop_id) == F] %>% unique
  
  keep <- c()
  
  for (id in unique(ref_list[[i]]$pop_id)) {
    diag <- ref_list[[i]][pop_id == id, diagnostic] %>% unique %>% sort
    if ((1 %in% diag) == T & (2 %in% diag) == T) {
      keep <- c(keep, id)
    }
    
  }
  
  extract_list[[i]] <- ref_list[[i]][pop_id %in% keep, .(pop_id, group, mean, upper, lower, cases, sample_size, diagnostic_name)]
  pop_group_ref <- ref_list[[i]][, .(pop_id, group, field_citation, file_path, location_name, lat, long, location_id,
                                     poly_id, poly_reference, sex, year_start, year_end, cv_MDAbase, cv_MDAend, cv_MDAschool,
                                     month_start, month_end, age_start, age_end, site_memo, filtration, n_rounds)] %>% unique
  names_numeric <- c("mean", "upper", "lower", "cases", "sample_size")
  for (col in names_numeric) set(extract_list[[i]], j=col, value=as.numeric(extract_list[[i]][[col]]))
  extract_list[[i]][is.na(cases) == F & is.na(sample_size) == F, mean := cases/sample_size]
  extract_list[[i]] <- data.table::dcast(extract_list[[i]], formula = pop_id ~ diagnostic_name, value.var = c("mean", "sample_size"), fun = mean, drop = F) %>% as.data.table
  extract_list[[i]][, source := "Lit Extraction"]
  extract_list[[i]][, specificity := spec]
  setnames(extract_list[[i]], c("mean_ict", "mean_mf"), c("ict", "mf"))
  
  extract_list[[i]] <- merge(extract_list[[i]], pop_group_ref, by = "pop_id")
  names_numeric <- c("year_start", "year_end", "month_start", "month_end", "age_start", "age_end")
  for (col in names_numeric) set(extract_list[[i]], j=col, value=as.numeric(extract_list[[i]][[col]]))
  
  extract_list[[i]]$location_id <- as.integer(extract_list[[i]]$location_id)
  extract_list[[i]][poly_id == 4197, location_id := as.integer(27)] #one row missing location_id from extraction, is American Samoa
  extract_list[[i]][location_id == 4870, location_id := 43937] #Tamil Nadu (total) location_id not in get_location, switch to Tamil Nadu_rural
  extract_list[[i]] <- merge(extract_list[[i]], loc_short, by = "location_id", all.x = T)
  
  extract_list[[i]][is.na(cv_MDAend) == T & is.na(cv_MDAschool) ==T, mda := 0]
  extract_list[[i]][is.na(mda) == T, mda := 1]
  
  i <- i + 1
}

# combine list
final_extract <- rbindlist(extract_list)
uid <- 1:nrow(final_extract)
final_extract <- cbind(final_extract, uid)
final_extract$pop_id <- NULL
keep <- c()
for (g in unique(final_extract$group)){
  t <- table(final_extract[group == g, specificity]) %>% as.data.table
  spec <- t[N == max(N), V1]
  keep <- c(keep, final_extract[group == g & specificity == spec[1], uid])
}

final_extract <- final_extract[uid %in% keep,]
pop_id <- 1:nrow(final_extract)
final_extract <- cbind(final_extract, pop_id)
final_extract$uid <- NULL

#get rid of filtered/not-filtered duplicates from extraction
exclude <- c()
for (i in unique(final_extract$group)){
  if (length(final_extract[group == i, filtration]) >1){
    exclude <- c(exclude, final_extract[group == i & filtration == 1, pop_id])
  }
}

final_extract <- final_extract[!(pop_id %in% exclude)]
setnames(final_extract, "site_memo", "Site_name")

### Check sources that were subsetted ##########################################################

if(FALSE){# get source information to check against Cano paper sources
  sources_gahi <- gahi_clean[, .(Reference_restr, Author_restr, Source_restr)] %>% unique
  sources_gahi_subset <- gahi_clean[uid_from_source_restr %in% grep("LF", gahi_ict_mf$pop_id, value = T), .(Reference_restr, Author_restr, Source_restr)]%>% unique
  sources_extraction <- extraction[, .(field_citation, file_path)] %>% unique
  sources_extraction_subset <- ict_mf_extract[, .(field_citation, file_path)] %>% unique
  sources <- c("sources_gahi", "sources_gahi_subset", "sources_extraction", "sources_extraction_subset")
  for (i in sources){
    write.csv(get(i), paste0(i, ".csv"))
  }
}

### Combine GAHI and extraction cleaned crosswalk datasets ################################################################
gahi_bind <- gahi_KK_sub[, .(pop_id, Record_ID, source, year_start, year_end, sex, ict, mf, sample_size_ict, sample_size_mf, location_name,
                             age_start, age_end, super_region_id, super_region_name, region_id, region_name, developed, mda, filtration, Site_name)]
extract_bind <- final_extract[, .(pop_id, source, year_start, year_end, sex, ict, mf, sample_size_ict, sample_size_mf, location_name,
                                  age_start, age_end, super_region_id, super_region_name, region_id, region_name, developed, mda, filtration, Site_name)]
extract_bind[, Record_ID := NA]
crosswalk_agg <- rbind(gahi_bind, extract_bind)

crosswalk_agg[mda == 1, mda_name := "post-MDA"]
crosswalk_agg[mda == 0, mda_name := "pre-MDA"]
crosswalk_agg$sample_size_ict <- as.integer(crosswalk_agg$sample_size_ict)
crosswalk_agg$sample_size_mf <- as.integer(crosswalk_agg$sample_size_mf)
crosswalk_agg$mda_name <- factor(crosswalk_agg$mda_name, levels = c("pre-MDA", "post-MDA"))
crosswalk_agg[,weight := (sample_size_ict + sample_size_mf)/2]
crosswalk_agg <- crosswalk_agg[is.na(mf) == F & is.na(ict) == F,]

### Data Quality Checks/Fixes ####################################################################

### Flag sample sizes issues
#zero denominators
crosswalk_agg[sample_size_ict == 0 | sample_size_mf == 0, flag_0_denom := 1]
crosswalk_agg[is.na(flag_0_denom) == T, flag_0_denom := 0]

#sample sizes that are greater than 50% different for review
crosswalk_agg[sample_size_ict/sample_size_mf >= 1.5 | sample_size_ict/sample_size_mf <= 0.5, flag_sample_size := 1]
crosswalk_agg[is.na(flag_sample_size) == T, flag_sample_size := 0]

#write.csv(gahi_KK[pop_id %in% crosswalk_agg[source == "GAHI" & (flag_0_denom == 1 | flag_sample_size == 1), pop_id],], "check_GAHI_sample_size_new.csv", row.names = F)

#read in csv from checking odd GAHI sample sizes and exclude data points found to be erroneous (68 dropped)
##NEED TO CHECK IF ANY DROPPED DATA CAN BE RECOVERED; ALSO CHECK FOR WHAT ICT TEST USED AND IF MORE DATA AVAILABLE
check_ss <- gahi_KK[pop_id %in% crosswalk_agg[source == "GAHI" & (flag_0_denom == 1 | flag_sample_size == 1), pop_id],]
checked_ss_gahi <- fread("check_GAHI_sample_size.csv")

check_ss <- merge(check_ss, checked_ss_gahi[,.(Record_ID, Comments)], by = "Record_ID", all.x = T)

#check remaining GAHI sources (37)
#write.csv(gahi_KK[pop_id %in% crosswalk_agg[source == "GAHI", pop_id] & !(pop_id %in% checked_ss_gahi[Comments == "valid", pop_id]),], "check_GAHI_all.csv", row.names = F)
check_gahi_all <- gahi_KK[pop_id %in% crosswalk_agg[source == "GAHI", pop_id] & !(pop_id %in% check_ss[Comments == "valid", pop_id]),]
checked_gahi_all <- fread("check_GAHI_all.csv")

check_gahi_all <- merge(check_gahi_all, checked_gahi_all[,.(Record_ID, Comments)], by = "Record_ID", all.x = T)

#check original extraction rows used in final dataset. 
extract_check <- final_extract[, .(specificity, file_path)] %>% unique
extract_check <- cbind(extract_check, 1:nrow(extract_check))
extract_check <- merge(extraction, extract_check, by = c("specificity", "file_path"))
#write.csv(extract_check, "extraction_crosswalk_review.csv", row.names = F)

if(TRUE){#adjust extraction data points that only tested mf in ICT+. Specificity from Irvine paper (98.4%)
  adjust_mf <- final_extract[file_path %in% (extract_check[grep("Only ICT+", extract_check$note_SR), file_path] %>% unique), pop_id]
  adjust_mf <- c(adjust_mf, check_gahi_all[grep("ICT positive", check_gahi_all$Comments), pop_id], 
                 check_ss[grep("only done", check_ss$Comments), pop_id])
  adjust_mf <- crosswalk_agg[pop_id %in% adjust_mf]
  
  #drop all that are missing sample size (cannot adjust without denominator) or has prev_ict of 0 (adjustment wouldn't do anything)
  adjust_mf <- adjust_mf[sample_size_ict != 0 & sample_size_mf != 0 & is.na(sample_size_ict) == F & is.na(sample_size_mf) == F & mf != 0,]
  adjust_mf[sample_size_ict != sample_size_mf, mf_denom_type := "only ICT+"]
  adjust_mf[sample_size_ict == sample_size_mf, mf_denom_type := "total"]
  adjust_mf[, mf_adj_sample_size := as.integer(ict * sample_size_ict/.984)]
  adjust_mf[mf_denom_type == "total", mf := (mf * sample_size_mf)/mf_adj_sample_size]
  adjust_mf[, mf_np := mf * mf_adj_sample_size]
  adjust_mf[, mf := mf_np/sample_size_ict]
  adjust_mf[, sample_size_mf := sample_size_ict]
  adjust_mf[, c("mf_adj_sample_size", "mf_np", "mf_denom_type") := NULL]
}

# drop gahi and extraction rows founds to be erroneous
crosswalk_agg <- crosswalk_agg[!(pop_id %in% check_ss[Comments != "valid" & is.na(Comments) == F, pop_id]),]
crosswalk_agg <- crosswalk_agg[!(pop_id %in% checked_gahi_all[Comments != "valid" & is.na(Comments) == F, pop_id]),]
crosswalk_agg <- crosswalk_agg[!(pop_id %in% final_extract[file_path %in% (extract_check[grep("Only ICT+", extract_check$note_SR), file_path] %>% unique), pop_id])]

if(exists("adjust_mf")){
  crosswalk_agg <- rbind(crosswalk_agg, adjust_mf)
}

if(FALSE){#move all sources used in crosswalk to a folder
  #move lit extraction sources
  filepaths <- extract_check[,file_path] %>% unique
  targetdir <- "J:/Project/NTDS/data/LF/multiple_diagnostics/crosswalk/"
  file.copy(from= filepaths, to = targetdir, overwrite = T, recursive = F, copy.mode = T)
  
  # create GAHI sources csv
  write.csv(rbind(checked_ss_gahi[, .(Reference, Author, Year)], checked_gahi_all[, .(Reference, Author, Year)]) %>% unique,
            paste0(targetdir, "sources_GAHI_crosswalk.csv"), row.names = F)
}

#drop GAHI extracted row of duplicated source (keep lit Extraction rows)
crosswalk_agg <- crosswalk_agg[!(pop_id %in% (gahi_KK[grep("Spatial clustering of filarial transmission", gahi_KK$Reference), pop_id] %>% unique))]

crosswalk_agg <- unique(crosswalk_agg)
crosswalk_agg$pop_id <- as.character(crosswalk_agg$pop_id)
crosswalk_agg[source == "Lit Extraction", pop_id := paste0(pop_id, "L")]

### Add age extraction data #########################################################################

age_extract <- fread("J:/Project/NTDS/data/LF/multiple_diagnostics/crosswalk/age_extractions_gahi.csv")
age_extract[Record_ID == "", Record_ID := Author]
age_extract <- age_extract[, .(pop_id, Record_ID, np_ICT, np_Blood_smear, np_Filtration, N_rounds, ADM0, ADM1, ADM2, TypeSite, Site_name,
                               Year_start, Year_end, Age_start, Age_end, Age_Range, Sex, pop_ICT, 
                               pop_Blood_smear, pop_Filtration, prev_Blood_smear, prev_Filtration, prev_ICT)]

age_extract[is.na(prev_Blood_smear) == F, c("prev_mf", "pop_mf", "np_mf", "filtration") := list(prev_Blood_smear, pop_Blood_smear, np_Blood_smear, 0)]
age_extract[is.na(prev_Blood_smear) == T, c("pop_mf", "np_mf", "filtration") := list(pop_Filtration, np_Filtration, 1)]
age_extract[,source := "GAHI"]
age_extract[N_rounds == 0, mda := 0]
age_extract[N_rounds > 0, mda := 1]
setnames(age_extract, c("ADM0", "Year_start", "Year_end", "prev_mf", "prev_ICT", "pop_ICT", "pop_mf", "Sex", "Age_start", "Age_end", "np_ICT"), 
         c("location_name", "year_start", "year_end", "mf", "ict", "sample_size_ict", "sample_size_mf", "sex", "age_start", "age_end", "np_ict"))
age_extract[, c("np_Blood_smear", "pop_Filtration", "pop_Blood_smear", "np_Filtration", "prev_Blood_smear", "prev_Filtration") := NULL]
age_extract <- merge(age_extract, gahi_loc_dict, by = "location_name", all.x = T)

age_extract <- age_extract[, .(pop_id, Record_ID, source, year_start, year_end, sex, ict, mf, sample_size_ict, sample_size_mf, location_name,
                               age_start, age_end, super_region_id, super_region_name, region_id, region_name, developed, mda, filtration, np_ict, np_mf, Site_name, N_rounds)]

age_extract[mda == 1, mda_name := "post-MDA"]
age_extract[mda == 0, mda_name := "pre-MDA"]
age_extract$sample_size_ict <- as.integer(age_extract$sample_size_ict)
age_extract$sample_size_mf <- as.integer(age_extract$sample_size_mf)
age_extract$mda_name <- factor(age_extract$mda_name, levels = c("pre-MDA", "post-MDA"))
age_extract[,weight := (sample_size_ict + sample_size_mf)/2]

crosswalk_agg[,np_ict := round(ict * sample_size_ict) %>% as.integer]
crosswalk_agg[,np_mf := round(mf * sample_size_mf) %>% as.integer]
crosswalk_agg[, c("flag_0_denom", "flag_sample_size") := NULL]
dropped <- crosswalk_agg[pop_id %in% age_extract$pop_id,]
crosswalk_agg <- rbind(crosswalk_agg[!(pop_id %in% age_extract$pop_id) & is.na(sample_size_ict) == F & is.na(sample_size_mf) == F,],
                       age_extract[sample_size_ict >= 20 & sample_size_mf >= 20])
### Covariates/grouping variables ###################################################
# age distribution
crosswalk_agg[(age_start == 0| is.na(age_start) == T) & (age_end == 0 | is.na(age_end) == T), age_group := "Unknown"]
crosswalk_agg[age_end >= 20 & age_start <= 20, age_group := "All ages"]
crosswalk_agg[age_start < 20 & age_end >= 10 & age_start >= 5 & age_end <= 20, age_group := "Adolescents"]
crosswalk_agg[age_start <= 8 & age_end <= 12 & (age_start != 0 & age_end != 0), age_group := "Children"]
crosswalk_agg[age_start >= 10 & age_end > 20, age_group := "Adults"]
crosswalk_agg[is.na(age_group) == T, age_group := "Adolescents"]

crosswalk_agg[pop_id %in% (gahi_KK_sub[pop_id %in% crosswalk_agg[age_group == "Unknown", pop_id]][Age_Range == "Adults/Children", pop_id]), age_group := "All ages"]
crosswalk_agg[pop_id %in% (gahi_KK_sub[pop_id %in% crosswalk_agg[age_group == "Unknown", pop_id]][Age_Range == "Adults", pop_id]), age_group := "Adults"]

#These data points were described in this paper (http://www.tandfonline.com/doi/pdf/10.1179/000349802125001735?needAccess=true)
crosswalk_agg[pop_id %in% (gahi_KK[pop_id %in% crosswalk_agg[age_group == "Unknown", pop_id] & Reference == "Programme to Eliminate Lymphatic Filariasis. Country Report WHO;GHA0020_LF2001", pop_id]), age_group  := "Adults"]

#These data points were described in this paper (http://journals.plos.org/plosntds/article?id=10.1371/journal.pntd.0001346)
crosswalk_agg[pop_id %in% (gahi_KK[pop_id %in% crosswalk_agg[age_group == "Unknown", pop_id] & Reference == "Programme to Eliminate Lymphatic Filariasis. Country Report WHO; NGA0013_LF2008", pop_id]), age_group  := "All ages"]

crosswalk_agg$age_group <- factor(crosswalk_agg$age_group, c("Children", "Adolescents", "Adults", "All ages"))

crosswalk_agg[weight < 50, ss_group := "less_50"]
crosswalk_agg[weight >= 50 & weight < 100, ss_group := "b_50_100"]
crosswalk_agg[weight >= 100 & weight < 500, ss_group := "b_100_500"]
crosswalk_agg[weight >= 500 & weight < 1000, ss_group := "b_500_1000"]
crosswalk_agg[weight >= 1000, ss_group := "more_1000"]
crosswalk_agg$ss_group <- factor(crosswalk_agg$ss_group, c("less_50", "b_50_100", "b_100_500", "b_500_1000", "more_1000"))

# precision estimates
crosswalk_agg[, c("sd_np_ict", "sd_np_mf") := 
                list(sqrt(sample_size_ict *ict * (1-ict)), sqrt(sample_size_mf *mf * (1-mf)))]
crosswalk_agg[, c("upper_ict", "upper_mf", "lower_ict", "lower_mf") := 
                list(((ict * sample_size_ict) + 1.96*sd_np_ict)/sample_size_ict, (((mf * sample_size_mf) + 1.96*sd_np_mf)/sample_size_mf),
                     (((ict * sample_size_ict) - 1.96*sd_np_ict)/sample_size_ict), (((mf * sample_size_mf) - 1.96*sd_np_mf)/sample_size_mf))]
crosswalk_agg$age_group <- factor(crosswalk_agg$age_group, c("Children","Adolescents", "Adults", "Adults & Children", "Unknown"))

# potential loa loa 
crosswalk_agg[unique(c(grep("Cameroon", location_name), grep("Nigeria", location_name))), loa := 1]
crosswalk_agg[is.na(loa) == T, loa := 0]


#get n_rounds for gahi data
n_rounds_gahi <- gahi_KK_sub[,.(pop_id, N_rounds)]
n_rounds_extract <- final_extract[, .(pop_id, n_rounds)]
n_rounds_extract[, pop_id := paste0(pop_id, "L")]
setnames(n_rounds_extract, "n_rounds", "N_rounds")
n_rounds_ageextract <- age_extract[grep("E", pop_id), .(pop_id, N_rounds)]
n_rounds_connect <- rbind(n_rounds_gahi, n_rounds_extract, n_rounds_ageextract)
crosswalk_agg <- merge(crosswalk_agg, n_rounds_connect, by = "pop_id", all.x= T)
crosswalk_agg[is.na(N_rounds) == T, N_rounds := 0]
crosswalk_agg$N_rounds <- as.integer(crosswalk_agg$N_rounds)

#year bins
crosswalk_agg[year_start < 2000, year_bin := "pre-2000"]
crosswalk_agg[year_start >= 2000 & year_start < 2005, year_bin := "2000-2004"]
crosswalk_agg[year_start >= 2005 & year_start < 2010, year_bin := "2005-2009"]
crosswalk_agg[year_start >= 2010, year_bin := "post-2010"]
crosswalk_agg$year_bin <- factor(crosswalk_agg$year_bin, levels = c("pre-2000", "2000-2004", "2005-2009", "post-2010"))

#N rounds bins
crosswalk_agg[N_rounds >= 3, N_rounds_bins := ">= 3"]
crosswalk_agg[N_rounds < 3, N_rounds_bins := as.character(N_rounds)]
crosswalk_agg$N_rounds_bins <- factor(crosswalk_agg$N_rounds_bins, levels = c("0", "1", "2", ">= 3"))

#ict prev bins
ict_groups <- cut(crosswalk_agg$ict, breaks = 3)
crosswalk_agg <- cbind(crosswalk_agg, ict_groups)

#ict^2
crosswalk_agg[, ict2 := ict^2]

#write.csv(crosswalk_agg, "cw_analysis.csv", row.names = F)

### Plots #########################################################################
simple_plot <- ggplot(data = crosswalk_agg, aes(x= ict, y= mf, weight = weight)) +
  geom_point(aes(shape = source, color = super_region_name)) +
  facet_wrap(~ mda_name) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "pICT and pMF in the same study group (GAHI and Lit Extraction)", x = "pICT (antigenemia)", y = "pMF (microfilariae)") +
  theme_minimal() +
  coord_fixed(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
  guides(shape = guide_legend(title="Source"), color = guide_legend(title="GBD Super Region"))
simple_plot

plot_mda_year <- ggplot(data = crosswalk_agg, aes(x = ict, y = mf, weight = weight)) +
  geom_point(aes(shape = source, color = mda_name)) +
  facet_wrap (~ year_start) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "pICT and pMF in the same study group pre-MDA (GAHI and Lit Extraction)", x = "pICT (antigenemia)", y = "pMF (microfilariae)") +
  theme_minimal() +
  coord_fixed(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
  guides(shape = guide_legend(title="Source"), color = guide_legend(title="MDA"))
plot_mda_year

region_plot <- ggplot(data = crosswalk_agg, aes(x= ict, y= mf, weight = weight)) +
  geom_point(aes(shape = source, color = mda_name)) +
  facet_wrap(~ region_name) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "pICT and pMF in the same study group (GAHI and Lit Extraction)", x = "pICT (antigenemia)", y = "pMF (microfilariae)") +
  theme_minimal() +
  coord_fixed(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
  guides(shape = guide_legend(title="Source"), color = guide_legend(title="Pre/Post MDA"))
region_plot

# year and sample size distribution of all data
year_plot <- ggplot(crosswalk_agg[year_start > 1980], aes(year_start, year_end)) +
  geom_point() +
  geom_jitter() +
  theme_minimal()
year_plot

sample_size_plot <- ggplot(data = crosswalk_agg, aes(x = sample_size_ict, y = sample_size_mf)) +
  geom_point(aes(shape = source, color = mda_name), position = "jitter") +
  facet_wrap(~ region_name) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Sample Size Comparison Between ICT and MF tests (same study)", x = "ICT sample size", y = "MF sample size") +
  theme_minimal() +
  coord_fixed(xlim = c(0, 1500), ylim = c(0, 1500)) +
  guides(shape = guide_legend(title="Source"), color = guide_legend(title="Pre/Post MDA"))
sample_size_plot

### Descriptive statistics #############################################################

age_table <- crosswalk_agg[,.(age_group, mda_name)] %>% table %>% as.data.table

age_viz <- ggplot(age_table, aes(x = mda_name, fill = mda_name)) +
  geom_col(aes(y = N)) +
  facet_wrap( ~ age_group) +
  theme_minimal() +
  scale_x_discrete(limits = c("pre-MDA", "post-MDA")) +
  scale_fill_discrete(name = "MDA", guide = guide_legend(reverse=TRUE)) +
  geom_text(aes(y = N, label= N), vjust= -0.5, size = 3) +
  labs(title = "Age Group Breakdown by MDA", x = "MDA", y = "Count") +
  guides(fill = F)

#year stats
year_range <- 1990:2015
n_points <- rep(0, length(year_range)) %>% as.integer
b_50_100 <- rep(0, length(year_range)) %>% as.integer
b_100_500 <- rep(0, length(year_range)) %>% as.integer
b_500_1000 <- rep(0, length(year_range)) %>% as.integer
more_1000 <- rep(0, length(year_range)) %>% as.integer

year_stats <- data.table(year_range, n_points, b_50_100, b_100_500, b_500_1000, more_1000)
setnames(year_stats, "year_range", "year")

for (i in year_range){
  temp <- crosswalk_agg[year_start == i & is.na(weight) == F,]
  year_stats[year == i, n_points := nrow(temp)]
  
  if(nrow(temp) > 0) {
    year_stats[year == i, less_50 := temp[weight < 50,] %>% nrow]
    year_stats[year == i, b_50_100 := temp[weight >= 50 &weight < 100,] %>% nrow]
    year_stats[year == i, b_100_500 := temp[weight >= 100 & weight < 500,] %>% nrow]
    year_stats[year == i, b_500_1000 := temp[weight >= 500 & weight < 1000,] %>% nrow]
    year_stats[year == i, more_1000 := temp[weight >= 1000,] %>% nrow]
    year_stats[year == i, adults := temp[age_group == "Adults",] %>% nrow]
    year_stats[year == i, children := temp[age_group == "Children",] %>% nrow]
    year_stats[year == i, all_ages := temp[age_group == "Adults & Children",] %>% nrow]
  }else{next}
}

year_viz <- ggplot(crosswalk_agg, aes(x = year_start)) +
  geom_histogram(binwidth=1, aes(fill = mda_name)) +
  theme_minimal() +
  labs(title = "Data by Year", x = "Year", y = "Count") +
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-0.5)

ss_table <- data.table::melt(year_stats[,.(year, less_50, b_50_100, b_100_500, b_500_1000, more_1000)], id = "year")
ss_viz <- ggplot(ss_table, aes(x = year)) +
  geom_col(aes(y = value, fill = variable)) +
  theme_minimal() +
  labs(title = "Sample Size breakdown by Year", x = "Year", y = "Count")

age_year_viz <- ggplot(crosswalk_agg, aes(x = year_start)) +
  geom_histogram(binwidth=1, aes(fill = mda_name)) +
  facet_wrap(~ age_group) +
  theme_minimal() +
  labs(title = "Data by Year & Age", x = "Year", y = "Count") +
  stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-0.5)

#with errors

simple_plot_err <- ggplot(data = crosswalk_agg, aes(x= ict, y= mf, color = super_region_name)) +
  geom_pointrange(aes(shape = source, ymin = lower_mf, ymax = upper_mf)) +
  geom_errorbarh(aes(xmin = lower_ict, xmax = upper_ict)) +
  facet_wrap(~ mda_name) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "pICT and pMF in the same study group (GAHI and Lit Extraction)", x = "pICT (antigenemia)", y = "pMF (microfilariae)") +
  theme_minimal() +
  coord_fixed(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
  guides(shape = guide_legend(title="Source"), color = guide_legend(title="GBD Super Region"))
simple_plot_err

## with errors and sample_size

simple_plot_err_ss <- ggplot(data = crosswalk_agg[is.na(ss_group) == F], aes(x= ict, y= mf, color = ss_group)) +
  geom_pointrange(aes(shape = source, ymin = lower_mf, ymax = upper_mf)) +
  geom_errorbarh(aes(xmin = lower_ict, xmax = upper_ict)) +
  facet_wrap(~ mda_name) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "pICT and pMF in the same study group (GAHI and Lit Extraction)", x = "pICT (antigenemia)", y = "pMF (microfilariae)") +
  theme_minimal() +
  coord_fixed(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
  guides(shape = guide_legend(title="Source"), color = guide_legend(title="Sample Size")) +
  scale_fill_discrete(limits=c("less_50", "b_50_100", "b_100_500", "b_500_1000", "NA"))
simple_plot_err_ss

## with age_group
simple_plot_err_ag <- ggplot(data = crosswalk_agg[is.na(ss_group) == F], aes(x= ict, y= mf, color = age_group)) +
  geom_pointrange(aes(shape = source, ymin = lower_mf, ymax = upper_mf)) +
  geom_errorbarh(aes(xmin = lower_ict, xmax = upper_ict)) +
  facet_wrap(~ mda_name) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "pICT and pMF in the same study group (GAHI and Lit Extraction)", x = "pICT (antigenemia)", y = "pMF (microfilariae)") +
  theme_minimal() +
  coord_fixed(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
  guides(shape = guide_legend(title="Source"), color = guide_legend(title="Sample Size"))
simple_plot_err_ag

# age group + mda boxplots
give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

crosswalk_agg_m <- data.table::melt(crosswalk_agg[,.(mda_name, ict, mf, age_group)], measure.vars = c("ict", "mf"))
crosswalk_agg_m$age_group <- factor(crosswalk_agg_m$age_group, c("Children", "Adults & Children", "Adults","NA"))
age_prev_bp <- ggplot(crosswalk_agg_m[is.na(age_group) == F], aes(y = value, x = interaction(age_group, variable))) +
  geom_boxplot(aes(fill = age_group)) +
  facet_grid( ~ mda_name)+
  theme_classic() +
  guides(fill = F) +
  labs(title = "Prevalence by Age Group and MDA", x = "Age Group", y = "Prevalence")+ 
  stat_summary(fun.data = give.n, geom = "text")

# age group + mda boxplots
crosswalk_agg_m <- data.table::melt(crosswalk_agg[,.(mda_name, ict, mf, ss_group)], measure.vars = c("ict", "mf"))
crosswalk_agg_m$ss_group <- factor(crosswalk_agg_m$ss_group, c("less_50", "b_50_100", "b_100_500", "b_500_1000", "more_1000"))
ss_prev_bp <- ggplot(crosswalk_agg_m[is.na(ss_group) == F], aes(y = value, x = interaction(ss_group, variable))) +
  geom_boxplot(aes(fill = ss_group)) +
  facet_grid( ~ mda_name)+
  theme_classic() +
  guides(fill = F) +
  labs(title = "Prevalence by Sample Size and MDA", x = "Sample Size", y = "Prevalence") + 
  scale_fill_brewer(palette="YlOrRd")+ 
  stat_summary(fun.data = give.n, geom = "text")

# sample size group and region
ss_region_viz <- ggplot(crosswalk_agg, aes(x = mda_name)) +
  geom_bar(aes(fill = ss_group)) +
  facet_grid( ~ region_name) +
  theme_classic() +
  labs(title = "Crosswalk Data by Sample Size, Region, and MDA", x = "Region", y = "Count")

age_extracted <- ggplot(data = crosswalk_agg[pop_id %in% age_extract$pop_id], aes(x= age_start, y= ict, color = location_name)) +
  geom_pointrange(aes(shape = source, ymin = lower_ict, ymax = upper_ict)) +
  facet_wrap(~ mda_name) +
  labs(title = "Extracted Age Split Age Trend (ICT)", x = "Age (age_start)", y = "pMF (microfilariae)") +
  theme_minimal() +
  guides(shape = guide_legend(title="Source"))

year_prev <- ggplot(data = crosswalk_agg, aes(x= year_start, y= mf*100, color = ict)) +
  geom_pointrange(aes(shape = source, ymin = lower_mf, ymax = upper_mf)) +
  geom_errorbarh(aes(xmin = lower_ict, xmax = upper_ict)) +
  facet_wrap(~ mda_name) +
  labs(title = "Relationship Between Survey Year and MF/ICT prevalences", x = "Survey Year (Start)", y = "pMF (% microfilariae)") +
  theme_minimal() +
  coord_fixed(xlim = c(1990, 2020), ylim = c(0, 80)) +
  geom_smooth(method = "loess") +
  guides(shape = guide_legend(title="Source"), color = guide_legend(title="ICT prevalence"))

#year bins

simple_plot_err_ybin <- ggplot(data = crosswalk_agg[mda == 1], aes(x= ict, y= mf)) +
  geom_pointrange(aes(shape = source, ymin = lower_mf, ymax = upper_mf, color = as.factor(N_rounds))) +
  geom_errorbarh(aes(xmin = lower_ict, xmax = upper_ict,color = as.factor(N_rounds))) +
  facet_wrap(~ year_bin) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "Post-MDA pICT and pMF (GAHI and Lit Extraction)", x = "pICT (antigenemia)", y = "pMF (microfilariae)") +
  theme_minimal() +
  coord_fixed(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
  guides(shape = guide_legend(title="Source"), color = guide_legend(title="N_rounds"))
simple_plot_err_ybin

simple_plot_err_rounds <- ggplot(data = crosswalk_agg, aes(x= ict, y= mf, weight = weight)) +
  geom_pointrange(aes(shape = source, ymin = lower_mf, ymax = upper_mf, color = super_region_name)) +
  geom_errorbarh(aes(xmin = lower_ict, xmax = upper_ict,color = super_region_name)) +
  facet_wrap(~ N_rounds_bins) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "pICT and pMF by # MDA Rounds (GAHI and Lit Extraction)", x = "pICT (antigenemia)", y = "pMF (microfilariae)") +
  theme_minimal() +
  coord_fixed(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
  geom_smooth(method = "loess", span = 3) +
  guides(shape = guide_legend(title="Source"), color = guide_legend(title="GBD Super Regions"))
simple_plot_err_rounds

age_stratified <- ggplot(data = crosswalk_agg[mda == 0], aes(x= ict, y= mf, weight = weight)) +
  geom_pointrange(aes(shape = source, ymin = lower_mf, ymax = upper_mf, color = region_name)) +
  geom_errorbarh(aes(xmin = lower_ict, xmax = upper_ict, color = region_name)) +
  facet_grid(~ age_group) +
  geom_abline(intercept = 0, slope = 1) +
  labs(title = "pICT and pMF (Pre-MDA)", x = "pICT (antigenemia)", y = "pMF (microfilariae)") +
  theme_minimal() +
  coord_fixed(xlim = c(0, 0.8), ylim = c(0, 0.8)) +
  geom_smooth(method = "lm") +
  guides(shape = guide_legend(title="Source"), color = guide_legend(title="GBD region"))
age_stratified