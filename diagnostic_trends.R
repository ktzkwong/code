# Multiple Diagnostics/ MF - ICT crosswalk

## Editing Author: Kevin Kwong
## Date: 4/3/2017
## Purpose: Identify all sources using more than 1 diagnostic in the full GAHI and lit extraction datasets. 
## Create a MF - ICT crosswalk dataset (same pop tested). Vet data quality issues. Preliminary exploration of crosswalk.

#########################################################################################################
rm(list=ls())
dev.off()

pacman::p_load(data.table, ggplot2, magrittr, reshape2, stringi, dplyr, ggrepel)

workdir <- "H:/code/"
setwd(workdir)

# load datasets
gahi <- fread("H:/code/gahi_KK_collapse_new.csv")
extraction <- fread("LF_Extraction_with_NIDs_cleaned_Diagnostics.csv")


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

#define some functions
psum <- function(..., na.rm=T) { 
  x <- list(...)
  rowSums(matrix(unlist(x), ncol=length(x)), na.rm=na.rm)
} 

names_same <- function(dt_new, dt_ref){
  names(dt_new)[names(dt_new) %in% names(dt_ref) == T]
}

names_diff <- function(dt_new, dt_ref){
  names(dt_new)[!(names(dt_new) %in% names(dt_ref) == T)]
}

### GAHI processing #######################################################################
gahi[, source := "GAHI"]

# diagnostics
dx_names <- unlist(strsplit(paste(gahi$diagnostics %>% unique, collapse = " "), split = " "))
dx_names <- dx_names[!(dx_names %in% "+")] %>% unique

for (name in dx_names){
  gahi[grep(name, gahi$diagnostics), c(name) := list(1)]
  gahi[!(grep(name, gahi$diagnostics)), c(name) := list(0)]
}

# mda
gahi[N_rounds == 0, c("mda", "mda_name") := list(0, "pre-MDA")]
gahi[N_rounds > 0, c("mda", "mda_name") := list(1, "post-MDA")]

#get GAHI loc file
gahi_loc_dict <- fread("GAHI_location_dictionary.csv")
gahi_loc_dict <- gahi_loc_dict[,.(location_name, super_region_id, super_region_name, region_id, region_name, developed)]
gahi <- merge(gahi, gahi_loc_dict, by.x = "ADM0", by.y = "location_name")
setnames(gahi, "ADM0", "location_name")
### Lit Extraction processing ################################################################

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

# need to first split by specificity
i <- 1
ref_list <- list()
id_list <- list()
extract_list <- list()
for (spec in unique(multi_subset$specificity)) {
  ref_list[[i]] <- multi_subset[specificity == spec,] %>% as.data.table
  
  id_list[[i]] <- ref_list[[i]][, .N, by = list(field_citation, file_path, location_name, lat, long, 
                                                poly_id, poly_reference, sex, year_start, year_end, 
                                                month_start, month_end, age_start, age_end, site_memo)] %>% as.data.table

  pop_id <- 1:nrow(id_list[[i]])
  id_list[[i]] <- cbind(id_list[[i]], pop_id)
  
  
  ref_list[[i]] <- merge(ref_list[[i]], id_list[[i]], by = c("field_citation", "file_path", "location_name", "lat", "long", 
                                                             "poly_id", "poly_reference", "sex", "year_start", "year_end", 
                                                             "month_start", "month_end", "age_start", "age_end", "site_memo"), all.x = T)
  
  extract_list[[i]] <- ref_list[[i]][, .(pop_id, group, mean, upper, lower, cases, sample_size, diagnostic_name)]
  pop_group_ref <- ref_list[[i]][, .(pop_id, group, field_citation, file_path, location_name, lat, long, location_id,
                                     poly_id, poly_reference, sex, year_start, year_end, cv_MDAbase, cv_MDAend, cv_MDAschool,
                                     month_start, month_end, age_start, age_end, site_memo, filtration)] %>% unique
  names_numeric <- c("mean", "upper", "lower", "cases", "sample_size")
  for (col in names_numeric) set(extract_list[[i]], j=col, value=as.numeric(extract_list[[i]][[col]]))
  extract_list[[i]][is.na(cases) == F & is.na(sample_size) == F, mean := cases/sample_size]
  extract_list[[i]] <- data.table::dcast(extract_list[[i]], formula = pop_id ~ diagnostic_name, value.var = c("mean", "sample_size", "cases"), 
                                         fun = mean, drop = F, fill = NA) %>% as.data.table
  extract_list[[i]][, source := "Lit Extraction"]
  extract_list[[i]][, specificity := spec]
  
  extract_list[[i]] <- merge(extract_list[[i]], pop_group_ref, by = "pop_id")
  names_numeric <- c("year_start", "year_end", "month_start", "month_end", "age_start", "age_end")
  for (col in names_numeric) set(extract_list[[i]], j=col, value=as.numeric(extract_list[[i]][[col]]))
  
  extract_list[[i]]$location_id <- as.integer(extract_list[[i]]$location_id)
  extract_list[[i]][poly_id == 4197, location_id := as.integer(27)] #one row missing location_id from extraction, is American Samoa
  extract_list[[i]][location_id == 4870, location_id := 43937] #Tamil Nadu (total) location_id not in get_location, switch to Tamil Nadu_rural
  extract_list[[i]] <- merge(extract_list[[i]], loc_short, by = "location_id", all.x = T)
  
  extract_list[[i]][is.na(cv_MDAend) == T & is.na(cv_MDAschool) ==T, c("mda", "mda_name") := list(0, "pre-MDA")]
  extract_list[[i]][is.na(mda) == T, c("mda", "mda_name") := list(1, "post-MDA")]
  
  i <- i + 1
}

# combine list
final_extract <- rbindlist(extract_list, fill = T)
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

setnames(final_extract, grep("mean", names(final_extract), value = T), gsub(".*_", "prev_", grep("mean", names(final_extract), value = T)))
setnames(final_extract, grep("cases", names(final_extract), value = T), gsub(".*_", "np_", grep("cases", names(final_extract), value = T)))
setnames(final_extract, grep("sample_size", names(final_extract), value = T), gsub(".*_", "pop_", grep("sample_size", names(final_extract), value = T)))
setnames(final_extract, names(final_extract), gsub("ict", "ICT", names(final_extract)))
setnames(final_extract, names(final_extract), gsub("mf", "Blood_smear", names(final_extract)))
setnames(final_extract, names(final_extract), gsub("fts", "New_rapid_test", names(final_extract)))
setnames(final_extract, names(final_extract), gsub("elisa", "ELISA", names(final_extract)))
setnames(final_extract, names(final_extract), gsub("igg4", "Others", names(final_extract)))
setnames(final_extract, names(final_extract), gsub("clinical", "Any_clinic", names(final_extract)))
setnames(final_extract, names(final_extract), gsub("age", "Age", names(final_extract)))
setnames(final_extract, c("field_citation","file_path", "filtration", "group", "lat", "long", "sex", "year_start", "year_end"),
                                   c("Library_ID", "Reference", "Filtration", "Record_ID", "Latitude", "Longitude", "Sex", "Year_start", "Year_end"))
final_extract[, names_diff(final_extract, gahi) := NULL]
final_extract[, source := "Lit Extraction"]

# get diagnostic columns and count for each row
dx_names <- c("Blood_smear", "ICT", "Others", "Any_clinic", "ELISA", "New_rapid_test")

for (name in dx_names){
  dx_cols <- grep(name, names(final_extract), value = T)
  for (col in dx_cols){
    final_extract[is.na(final_extract[, get(col)]) == F, c(name) := list(1)]
  }
  final_extract[is.na(get(name)) == T, c(name) := list(0)]
}

gahi[, source := "GAHI"]
total <- rbind(gahi, final_extract, fill = T)


# age distribution
total[Age_start > 18 & Age_end > 18, age_group := "Adults"]
total[Age_start <= 18 & Age_end <= 18 & Age_start != 0 & Age_end != 0, age_group := "Children"]
total[Age_start <= 18 & Age_end > 18, age_group := "Adults & Children"]
total[Age_start == 0 & Age_end == 0, age_group := "Unknown"]

#regenerate full dx_names
dx_names <- unlist(strsplit(paste(gahi$diagnostics %>% unique, collapse = " "), split = " "))
dx_names <- dx_names[!(dx_names %in% "+")] %>% unique

### total Year Stats ######################################################################

age_table <- total[,.(age_group, mda_name)] %>% table %>% as.data.table

levels(total$mda_name) <- c("pre-MDA", "post-MDA")
age_viz <- ggplot(age_table, aes(x = mda_name, fill = mda_name)) +
  geom_col(aes(y = N)) +
  facet_wrap( ~ age_group) +
  theme_minimal() +
  scale_fill_discrete(name = "MDA", guide = guide_legend(reverse=TRUE)) +
  geom_text(aes(y = N, label= N), vjust= -0.5, size = 3) +
  labs(title = "Age Group Breakdown by MDA", x = "MDA", y = "Count") +
  scale_x_discrete (limits = c("pre-MDA", "post-MDA")) +
  guides(fill = F)

# sample size groups by year
year_range <- 1990:2015
n_points <- rep(0, length(year_range)) %>% as.integer

year_stats <- data.table(year_range, n_points)
setnames(year_stats, "year_range", "year")

for (i in year_range){
  temp <- total[Year_start == i,]
  year_stats[year == i, n_points := nrow(temp)]
  year_stats[year == i, n_sources_tot := length(unique(temp$Library_ID))]
  year_stats[year == i , n_sources_site := length(unique(temp[is.na(object_ids) == T, Library_ID]))]
  year_stats[year == i , n_sources_gahi := length(unique(temp[is.na(object_ids) == F, Library_ID]))]
  
  if(nrow(temp) > 0) {
    year_stats[year == i, adults := temp[age_group == "Adults",] %>% nrow]
    year_stats[year == i, children := temp[age_group == "Children",] %>% nrow]
    year_stats[year == i, all_ages := temp[age_group == "Adults & Children",] %>% nrow]
    year_stats[year == i, pre_mda := temp[mda == 0,] %>% nrow]
    year_stats[year == i, post_mda := temp[mda == 1,] %>% nrow]
    
    for (name in dx_names){
      year_stats[year == i, c(name) := nrow(temp[get(name) == 1,])]
    }
    
    for (reg in unique(total$REGION)){
      year_stats[year == i, c(reg) := nrow(temp[REGION == reg,])]
    }
    
  }else{next}
}

adults <- rep(0, length(dx_names))
children <- rep(0, length(dx_names))
all_ages <- rep(0, length(dx_names))
n_points <- rep(0, length(dx_names))
dx_stats <- data.table(dx_names, n_points, adults, children, all_ages)

age_prev_dt <- data.table(pop_id = integer(), diagnostic = character(), prev = double(), age_group = character(), mda_name = character())

for (dx in dx_names){
  temp <- total[get(dx) == 1,]
  if(nrow(temp) > 0) {
    dx_stats[dx_names == dx, n_points := temp %>% nrow]
    dx_stats[dx_names == dx, only_dx := temp[num_diagnostics == 1,] %>% nrow]
    dx_stats[dx_names == dx, no_ict_mf := temp[ICT == 0 & Blood_smear == 0 & Filtration == 0,] %>% nrow]
    dx_stats[dx_names == dx, adults := temp[age_group == "Adults",] %>% nrow]
    dx_stats[dx_names == dx, children := temp[age_group == "Children",] %>% nrow]
    dx_stats[dx_names == dx, all_ages := temp[age_group == "Adults & Children",] %>% nrow]
    dx_stats[dx_names == dx, pre_mda := temp[mda == 0,] %>% nrow]
    dx_stats[dx_names == dx, post_mda := temp[mda == 1,] %>% nrow]
    for (i in c(0, 1)){
      if (i == 0){
        m <- "pre_mda"
      } else{ m <- "post_mda"}
      
      for (ag in c("Adults", "Children", "Adults & Children", "Unknown")){
        rows <- temp[mda == i & age_group == ag,] %>% nrow
        if (rows == 0){
          dx_stats[dx_names == dx, paste0(m, "_", ag) := 0]
        } else{
          dx_stats[dx_names == dx, paste0(m, "_", ag) := rows]
          
          names <- c("pop_id", "age_group", paste0("prev_", dx), "mda_name")
          t <- temp[mda == i & age_group == ag, names, with = F]
          diagnostic <- rep(dx, times = rows)
          t <- cbind(t, diagnostic)
          setnames(t, paste0("prev_", dx), "prev")
          age_prev_dt <- rbind(age_prev_dt, t, fill = T)
        }
      }
    }
  }
}

give.n <- function(x){
  return(c(y = mean(x), label = length(x)))
}

age_prev_dt$age_group <- factor(age_prev_dt$age_group, c("Children", "Adults & Children", "Adults", "Unknown"))
age_prev_pre_MDA <- ggplot(age_prev_dt[mda_name == "pre-MDA" & diagnostic %in% c("Blood_smear", "ICT") & prev != 0], aes(y = prev)) +
  geom_boxplot(aes(x = age_group, color = age_group)) +
  facet_grid( ~ diagnostic)+
  theme_minimal() +
  guides(fill = F) +
  labs(title = "Pre-MDA Prevalence by Age Group and Diagnostic", x = "Age Group", y = "Prevalence")

dx_cols <- c("dx_names", grep("pre_mda_", names(dx_stats), value = T), grep("post_mda_", names(dx_stats), value = T))
dx_tb <- dx_stats[, dx_cols, with = F]
dx_tb <- data.table::melt(dx_tb, id = "dx_names")
dx_tb[grep("pre_mda", dx_tb$variable), mda_name := "pre_MDA"]
dx_tb[grep("post_mda", dx_tb$variable), mda_name := "post_MDA"]
dx_tb[grep("Adults", dx_tb$variable), age_group := "Adults"]
dx_tb[grep("Children", dx_tb$variable), age_group := "Children"]
dx_tb[grep("Unknown", dx_tb$variable), age_group := "Unknown"]
dx_tb[grep("Adults & Children", dx_tb$variable), age_group := "Adults & Children"]

common_dx <- c("Any_clinic", "Blood_smear", "ELISA", "Hydrocele", "ICT", "Lymphedema")

dx_viz <- ggplot(dx_tb[dx_names %in% common_dx,], aes(x = mda_name)) +
  geom_col(aes(y = value, fill = age_group)) +
  facet_grid( ~ dx_names) +
  theme_classic() +
  labs(title = "Age Group Breakdown by Diagnostic and MDA (GAHI and Lit Extraction)", x = "Age Group", y = "Count") +
  guides(fill = guide_legend(title="Age Group"))

# KK - should look into why some prev are >1
# KK - why are there so many 0 ICT prev for pre-MDA? Because not truly pre-MDA?
