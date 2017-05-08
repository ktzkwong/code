# Collapse gahi dataset

## Author: Kevin Kwong
## Date: 4/17/2017
## Purpose: Collapse long format of gahi dataset so multiple diagnostics for the same location + year + population is swung wide.
##          Compare to what Erin Stearns did. 

#########################################################################################################
rm(list=ls())
dev.off()

pacman::p_load(data.table, ggplot2, magrittr, reshape2, stringi)

workdir <- "H:/code/"
setwd(workdir)

### load datasets ###############################################################################
gahi_raw <- fread("GAHI_LF.csv") #restricted and unrestricted GAHI data; from Rachel Pullen
site_data <- fread("J:/DATA/Incoming Data/thiswormyworld_data/Liz/LF/site_extra.csv") #pulled from GAHI dataset by Erin S.

# data quality datasets created by Liz; filtered from site_data
verify_polygon <- fread("J:/Project/NTDS/data/LF/verify_polygons/verify_polygons.csv") #KK - reviewed website data point/polygon
new_batch <- fread("J:/Project/NTDS/data/LF/verify_polygons/new_batch.csv") #KK - batch2 of reviewing point/polygon
verify_duplicates <- fread("J:/Project/NTDS/data/LF/verify_duplicates/verify_duplicates.csv") #KK - reviewed location + year duplicates

if(TRUE){### Fix data quality issues in Rachel's dataset #################################################
#Documentation of changes in "J:\Project\NTDS\data\LF\outstanding_gahi_issues.xlsx"

gahi_raw[Record_ID == "LF2148KE", ADM0_ID := 133]
gahi_raw <- gahi_raw[!(OBJECTID %in% c(5012, 5022, 5036, 5048, 5052, 5066, 9308, 9309, 9310, 7047, 7046, 5012, 5022, 5036, 5048, 5066,
                           6154, 6156, 6158, 6160, 6162, 6164, 6166, 9835, 9836, 9111, 9112))]
gahi_raw[Record_ID %in% c("LF10421TZ", "LF10422TZ", "LF10423TZ"), c("Age_start", "Age_end") := list(5, 16)]
gahi_raw[OBJECTID == 6377, pop := 159]
gahi_raw[OBJECTID == 9306, Record_ID := paste(Record_ID, "_1")]
gahi_raw[OBJECTID == 3406, c("pop", "np", "Method_2") := list(770, 127, "ELISA")]
gahi_raw[OBJECTID %in% c(15280, 15279), Record_ID := "LF10605GH"]
gahi_raw[OBJECTID %in% c(15283, 15284), Record_ID := "LF10608GH"]
gahi_raw[OBJECTID %in% c(15278, 15282), LF_species := "Brugia malayi"]
gahi_raw[OBJECTID == 9112, Record_ID := "LF7129PF"]

# Q!? - has 0 or super early year_start rows been checked for extraction error? Might be throwing away a lot of data needlessly (3473)
gahi_raw <- gahi_raw[Year_start > 1980 | Year_start == 0,]

#fix coordinates
coords_found <- verify_duplicates[is.na(Latitude) == F, .(Latitude, Longitude, Autonum_ref)]
for (num in coords_found$Autonum_ref){
  gahi_raw[Record_ID == num, Latitude := coords_found[Autonum_ref == num, Latitude]]
  gahi_raw[Record_ID == num, Longitude := coords_found[Autonum_ref == num, Longitude]]
}

#drop points
gahi_raw <- gahi_raw[!(Record_ID %in% verify_duplicates[delete == "yes", Autonum_ref])]

# fix years from verify_duplicates
for (num in verify_duplicates[is.na(Year_start) == F, Autonum_ref]){
  gahi_raw[Record_ID == num, Year_start := verify_duplicates[Autonum_ref == num, Year_start]]
}

for (num in verify_duplicates[is.na(Year_end) == F, Autonum_ref]){
  gahi_raw[Record_ID == num, Year_end := verify_duplicates[Autonum_ref == num, Year_end]]
}
}

if(TRUE){### clean up data pulled from the GAHI website (by Erin) ###############################
names_same <- function(dt_new, dt_ref){
  names(dt_new)[names(dt_new) %in% names(dt_ref) == T]
}

names_diff <- function(dt_new, dt_ref){
  names(dt_new)[!(names(dt_new) %in% names(dt_ref) == T)]
}

#check field similarities and differences 
shared_fields <- names_same(site_data, gahi_raw)
names_diff(site_data, gahi_raw)

gahi_name_list <- c("Record_ID","CONTINENT", "REGION", "WHO_region", "ISO_CTRY", "ADM0", "ADM0_ID", "ADM1", "ADM1_ID", "ADM2",
                    "ADM2_ID", "TypeSite", "Site_name", "Precision", "Latitude", "Longitude", "PoT", "Coverage", "Trial", "TimePeriod",
                    "Year_start", "Year_end", "N_rounds", "Method_0", "LF_species", "Sex", "Age_Range", "Age_start", "Age_end",
                    "Library_ID", "Author", "Source", "Year", "Approval", "Reference")

missing_id_fields <- gahi_name_list[!(gahi_name_list %in% shared_fields)]

# update site data with coordinates found and verified
coords_found <- rbind(verify_polygon[is.na(Latitude) == F, .(Latitude, Longitude, Autonum_ref)],
                      new_batch[is.na(Latitude) == F, .(Latitude, Longitude, Autonum_ref)])
for (num in coords_found$Autonum_ref){
  site_data[Autonum_ref == num, Latitude := coords_found[Autonum_ref == num, Latitude]]
  site_data[Autonum_ref == num, Longitude := coords_found[Autonum_ref == num, Longitude]]
}

# delete point site data rows which cannot be geolocated and/or whose sources are not accessible (Verify_polygon, new_batch)
# delete unrepresentative study and study from 1930s (new_batch)
unusable_data <- c(verify_polygon[is.na(Latitude) == T, Autonum_ref], 
                   new_batch[Source_data == "IDN0040LF" | Source_data == "TLS0004LF", Autonum_ref])
site_data <- site_data[!(Autonum_ref %in% unusable_data)]

# fix incorrect Record_IDs
fix_rec_ids <- c("LF4670HA", "LF4671HA", "LF4672HA", "LF4673HA", "LF4674HA", "LF4675HA", "LF4667HA", 
                 "LF4668HA", "LF4669HA", "LF4717HA")

site_data[Autonum_ref %in% fix_rec_ids, Autonum_ref := gsub('.{1}$', 'T', Autonum_ref)]

#rename fields in site_data
setnames(site_data, c("Autonum_ref", "Low_ADM_Level", "ADM0_NAME", "ADM0_CODE", "ADM1_NAME", "ADM1_CODE", 
                      "Site.name", "GEO_RELIABILITY", "Control_trial", "Period", 
                      "Source_data", "Source_data1", "Source_data2", "Source_data3"), 
         c("Record_ID", "TypeSite", "ADM0", "ADM0_ID", "ADM1", "ADM1_ID", "Site_name", "Precision", 
           "Trial", "TimePeriod", "Library_ID", "Author", "Source", "Year"))

#take care of ADM2 columns (given v. assigned)
site_data[is.na(ASSIGNED_ADM2_NAME) == F, ADM2 := ASSIGNED_ADM2_NAME]
site_data[is.na(ASSIGNED_ADM2_CODE) == F, ADM2_ID := ASSIGNED_ADM2_CODE]
site_data[is.na(ASSIGNED_ADM2_NAME) == T, ADM2 := GIVEN_ADM2_NAME]
site_data[is.na(ASSIGNED_ADM2_CODE) == T, ADM2_ID := GIVEN_ADM2_CODE]

#melt 
site_data[, c("V1", "X", "GIVEN_ADM2_NAME", "GIVEN_ADM2_CODE","ASSIGNED_ADM2_NAME", "ASSIGNED_ADM2_CODE", "Method_1",
              "Type_District", "Name_Positioning", "Urban.rural", "Cartography","Year_survey", "month", "Method_2",
              "Day.Night", "Status", "mfload_method", "Reviewed", grep("BS_", names(site_data), value = T)) := NULL]

site_data_m <- data.table::melt(site_data, measure.vars = patterns("pop_", "prev_", "np_"), value.name = c("pop", "prev", "np")) %>% unique
Method_1 <- c("Parasitological", "Serological", NA, "Clinical", "Clinical", "Clinical", "Parasitological")
Method_2 <- c("Blood smear", "ICT", "Others", "Hydrocele", "Lymphedema", "Any clinic", "Filtration")
variable <- c(1:7)
var_names <- data.table(variable, Method_1, Method_2)
site_data_m <- merge(site_data_m, var_names, by = "variable", all.x = T)
site_data_m[, names_diff(gahi_raw, site_data_m) := NA]
names_numeric <- c("pop", "np", "prev", "Blood_volu")
for (col in names_numeric) set(site_data_m, j=col, value=as.numeric(site_data_m[[col]]))
site_data_m <- site_data_m[is.na(pop) == F | is.na(np) == F | is.na(prev) == F,]
site_data_m[is.na(prev) == T & is.na(pop) == F & is.na(np) == F, prev := np/pop]
site_data_m[, variable := NULL]

# append site data to GAHI dataset provided by Rachel Pullen
gahi_raw <- rbind(gahi_raw, site_data_m)
}


# generate unique ids for same loc + time + pop
# any data points that are not identical to any other data point for all the referenced fields have an unique id
gahi_ref <- gahi_raw[, .N, by = list(Record_ID, CONTINENT, REGION, WHO_region, ISO_CTRY, ADM0, ADM0_ID, ADM1, ADM1_ID, ADM2,
                                     ADM2_ID, TypeSite, Site_name, Precision, Latitude, Longitude, PoT, Coverage, Trial, TimePeriod,
                                     Year_start, Year_end, N_rounds, Method_0, LF_species, Sex, Age_Range, Age_start, Age_end,
                                     Library_ID, Author, Source, Year, Approval, Reference)]
pop_id <- c(1: nrow(gahi_ref))
gahi_ref <- cbind(gahi_ref, pop_id)
gahi_ref$N <- NULL

gahi <- merge(gahi_raw, gahi_ref, by = gahi_name_list, all.x = T) %>% as.data.table

# cast measurement columns wide by unique pop_id
gahi <- gahi[,Method_2 := gsub(" ", "_", Method_2, fixed = TRUE)]

measurements_wide <- data.table::dcast(gahi, pop_id ~ Method_2, value.var = c("Blood_volu", "prev", "density", "pop", "np"), 
                                       fun.aggregate = mean, fill = NA)

# get count of diagnostics, vector of diagnostic names, and vector of record_ids collapsed into pop_id for each unique pop_id 
for (i in unique(gahi$pop_id)) {
  diag_vec <- gahi[pop_id == i, Method_2] %>% unique
  diag_vec_collapsed <- paste(diag_vec, collapse = " + ")
  
  o_ids <- gahi[pop_id == i, OBJECTID] %>% unique
  object_ids_collapsed <- paste(o_ids, collapse = " + ")
  
  measurements_wide[pop_id == i, diagnostics := diag_vec_collapsed %>% as.character]
  measurements_wide[pop_id == i, num_diagnostics := length(diag_vec) %>% as.integer]
  measurements_wide[pop_id == i, object_ids := object_ids_collapsed %>% as.character]
}

# create final collapsed, cleaned dataset and export
gahi_clean <- merge(gahi_ref, measurements_wide, by = "pop_id")

# double checking
if(TRUE){
  check <- gahi_clean$Library_ID %>% table %>% sort(decreasing = T) %>% as.data.table
  check_ids <- check[N > 1, .]
  
  # legitimate reasons for multiple pops for same Library_ID (different age group, pre/post treatment, MDA round, time or location of survey)
  correct_dups_reasons <- c("Age_start", "Age_end", "Age_Range", "ADM0", "ADM1", "ADM2", "Site_name", "TypeSite",
                            "Latitude", "Longitude", "Year_start", "Year_end", "N_rounds", "PoT") 
 
  double_check <- c()
  cor_dup <- c()
  dup_reason <- c()
  
  for (id in check_ids) {
    check_temp <- gahi_clean[Library_ID == id]
    check_temp[N_rounds == 0, N_rounds := 100] #for some reason setdiff doesn't work when value = 0
    check_temp[N_rounds == 1, N_rounds := 200] #also when value = 1
    cor_dup_test <- c()
    for (i in 1:(nrow(check_temp) - 1)){
      diff <- base::setdiff(check_temp[i], check_temp[i+1])
      cor_dup_test <- c(cor_dup_test, intersect(correct_dups_reasons, names(diff)))
    }
    cor_dup_test <- unique(cor_dup_test)
    if ( length(cor_dup_test) == 0) { #rows that differ only in diagnostics columns without being different for legitimate reasons are captured here
      double_check <- c(double_check, id)
    } else { #all others rows are sorted by the manner by which they differ from others with the same Library_ID
      cor_dup <- c(cor_dup, id)
      reason <- c()
      if (length(grep("Age", cor_dup_test)) > 0) {
        reason <- c(reason, "age")
      }
      if (length(grep("ADM", cor_dup_test)) > 0 | length(grep("tude", cor_dup_test)) > 0 | length(grep("Site", cor_dup_test)) > 0) {
        reason <- c(reason, "location")
      }
      if (length(grep("Year", cor_dup_test)) > 0) {
        reason <- c(reason, "year")
      }
      if (length(grep("rounds", cor_dup_test)) > 0) {
        reason <- c(reason, "MDA round")
      }
      if (length(grep("PoT", cor_dup_test)) > 0) {
        reason <- c(reason, "Pre/Post Intervention")
      }
      dup_reason <- c(dup_reason, paste(reason, collapse = " + "))
    }
  }
  
  double_check <- gahi_clean[Library_ID %in% double_check]
  correct_dups_ref <- cbind(cor_dup, dup_reason) %>% as.data.table
  setnames(correct_dups_ref, c("cor_dup", "dup_reason"), c("Library_ID", "specificity"))
  correct_dups <- merge(correct_dups_ref, gahi_clean, by = "Library_ID",all.x = T)
  
  if (nrow(double_check) > 1) {
    print("Warning! Potential collapse/data issue(s) detected, please double_check")
  }
  
  if (nrow(correct_dups) > 1) {
    print("There appears to be legitimate instances of multiple groups for the same Library_ID. Details in correct_dups and correct_dups_ref")
  }
}

gahi_clean <- merge(gahi_clean, correct_dups_ref, by = "Library_ID", all.x = T)
gahi_clean[is.na(specificity) == T, specificity := "Total"]

# fix aggregate prevalences
if(FALSE){ # since we're including record_id in determining pop_id, we don't need to do this
gahi_clean[grep(" +", gahi_clean$record_ids), prev_mean_Any_clinic := np_sum_Any_clinic/pop_sum_Any_clinic]
gahi_clean[grep(" +", gahi_clean$record_ids), prev_mean_Blood_smear := np_sum_Blood_smear/pop_sum_Blood_smear]
gahi_clean[grep(" +", gahi_clean$record_ids), prev_mean_Chamber := np_sum_Chamber/pop_sum_Chamber]
gahi_clean[grep(" +", gahi_clean$record_ids), prev_mean_ELISA := np_sum_ELISA/pop_sum_ELISA]
gahi_clean[grep(" +", gahi_clean$record_ids), prev_mean_Filtration := np_sum_Filtration/pop_sum_Filtration]
gahi_clean[grep(" +", gahi_clean$record_ids), prev_mean_Hydrocele := np_sum_Hydrocele/pop_sum_Hydrocele]
gahi_clean[grep(" +", gahi_clean$record_ids), prev_mean_ICT := np_sum_ICT/pop_sum_ICT]
gahi_clean[grep(" +", gahi_clean$record_ids), prev_mean_IFI := np_sum_IFI/pop_sum_IFI]
gahi_clean[grep(" +", gahi_clean$record_ids), prev_mean_Knott := np_sum_Knott/pop_sum_Knott]
gahi_clean[grep(" +", gahi_clean$record_ids), prev_mean_Lymphedema := np_sum_Lymphedema/pop_sum_Lymphedema]
gahi_clean[grep(" +", gahi_clean$record_ids), prev_mean_New_rapid_test := np_sum_New_rapid_test/pop_sum_New_rapid_test]
gahi_clean[grep(" +", gahi_clean$record_ids), prev_mean_PCR := np_sum_PCR/pop_sum_PCR]
}

write.csv(gahi_clean, "gahi_KK_collapse_new.csv", row.names = F)
