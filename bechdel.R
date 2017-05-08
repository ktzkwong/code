## Playing around - fivethirtyeight - Bechdel test
## Editing Author: ktzkwong
## Date: 3/22/2017
## Purpose: 
##
#########################################################################################################

### Script Prep #################################################################################

## Set directory and load libraries
rm(list=ls())

if(!require(pacman)) {
  install.packages("pacman"); require(pacman)}

pacman::p_load(dplyr, data.table, ggplot2, knitr, magrittr, broom, stringr, ggthemes, scales, fivethirtyeight)

## OS locals
os <- .Platform$OS.type
if (os == "windows") {
  h <- "H:/"
} else {
  h <- "/home/h/"
}

work_dir <- paste0(h, "code/")
setwd(work_dir)

### Playground #############################################################################
vignette("bechdel", package = "fivethirtyeight")
options(scipen = 99)
data("bechdel")

#restrict dataset between 1990 and 2013
bechdel90_13 <- bechdel %>% filter(between(year, 1990, 2013)) %>% as.data.table

#create descriptive stats columns regarding gross
bechdel90_13 %<>% 
  mutate(int_only = intgross_2013 - domgross_2013,
         roi_total = intgross_2013 / budget_2013,
         roi_dom = domgross_2013 / budget_2013,
         roi_int = int_only / budget_2013)

#generous test variable based on whether or not either test or clean_test returns the movie as "ok" or "dubious"
#regarding passing the Bechdel test. Binary is the more rigorous test

bechdel90_13 %<>%
  mutate(generous = ifelse(test = clean_test %in% c("ok", "dubious"),
                           yes = TRUE,
                           no = FALSE))

#calculate median ROI using Binary as grouping var. 
by_binary <- bechdel90_13[, list(median_ROI = median(roi_total, na.rm = TRUE), 
                                     mean_ROI = mean(roi_total, na.rm = TRUE), 
                                 median_budget = median(budget_2013, na.rm = TRUE), 
                                 mean_budget = mean(budget_2013, na.rm = TRUE)), by = list(binary)]
by_binary

#overall median and mean ROI
bechdel90_13 %>% 
  summarize(
    `Median Overall Return on Investment` = median(roi_total, na.rm = TRUE), `Mean` = mean(roi_total, na.rm = T))

#graphs
ggplot(data = bechdel90_13, mapping = aes(x = budget)) +
  geom_histogram(color = "white", bins = 20) +
  labs(title = "Histogram of budget")

ggplot(data = bechdel90_13, mapping = aes(x = log(budget))) +
  geom_histogram(color = "white", bins = 20) +
  labs(title = "Histogram of Logarithm of Budget")

ggplot(data = bechdel90_13, mapping = aes(x = intgross_2013)) +
  geom_histogram(color = "white", bins = 20) +
  labs(title = "Histogram of International Gross")

ggplot(data = bechdel90_13, mapping = aes(x = log(intgross_2013))) +
  geom_histogram(color = "white", bins = 20) +
  labs(title = "Histogram of Logarithm of International Gross")

ggplot(data = bechdel90_13, mapping = aes(x = roi_total)) +
  geom_histogram(color = "white", bins = 20) +
  labs(title = "Histogram of ROI")

ggplot(data = bechdel90_13, mapping = aes(x = roi_total)) +
  geom_histogram(color = "white", bins = 20) +
  labs(title = "Histogram of ROI") +
  xlim(0, 25)

ggplot(data = bechdel90_13, mapping = aes(x = log(roi_total))) +
  geom_histogram(color = "white", bins = 20) +
  labs(title = "Histogram of Logarithm of ROI")