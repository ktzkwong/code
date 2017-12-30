## Playing around - FedData
## Editing Author: ktzkwong
## Date: 3/22/2017
## Purpose: 
## Learn to code and work with spatial data. 
## http://zevross.com/blog/2016/03/15/using-the-new-r-package-feddata-to-access-federal-open-datasets-including-interactive-graphics/ 

#########################################################################################################

### Script Prep #################################################################################

## Set directory and load libraries
rm(list=ls())

if(!require(pacman)) {
  install.packages("pacman"); require(pacman)}

pacman::p_load(Rtools, data.table, ggplot2, magrittr, reshape2, dplyr, Hmisc., plyr, rgdal, 
               FedData, raster, leaflet, ggmap, rgeos, tidyr, xts, dygraphs, sp, tigris)

## OS locals
os <- .Platform$OS.type
if (os == "windows") {
  h <- "H:/"
} else {
  h <- "/home/h/"
}

work_dir <- paste0(h, "code/")
setwd(work_dir)

### Define extent #############################################################################
seattle <- get_map("Seattle, WA", zoom = 11)
ggmap(seattle)

#get bounding box coordinates from seattle map
bb <- attr(seattle, "bb")

extentB<-polygon_from_extent(raster::extent(bb$ll.lon, bb$ur.lon,
                                            bb$ll.lat, bb$ur.lat), 
                             proj4string="+proj=longlat +ellps=GRS80 +datum=NAD83 +no_defs")

ggmap(seattle) + geom_polygon(aes(x=long, y=lat), size=3, color='purple',
                             data=extentB, fill=NA)

# Download NED elevation data. This will return a raster
# of elevation data cropped to our extent area. Res = "1" 
# refers to 1 arc-second NED. This is the default. To 
# download 1/3 arc-second NED use res = "13". I'm adding 
# force.redo = F to both download functions. This will 
# check to see if the data already exists in which case 
# it will not re-download it.
ned_seattle<-get_ned(template=extentB, label="ned_seattle", 
                    res="1", force.redo = F)


# Download GHCN daily climate data. We're interested in 
# precipitation (PRCP) but many other climate-related
# elements exist (i.e. max and min temp, snowfall, etc).
# See the `FedData` documentation for more details. Also 
# of note is that a character vector of station IDs can 
# be used in place of our extent object
ghcn_seattle<-get_ghcn_daily(template=extentB, label="ghcn_seattle",
                            elements="prcp", force.redo = F)

sp::plot(ned_seattle)
sp::plot(ghcn_seattle$spatial, add=TRUE)

# Plot the data using ggmap and ggplot2. Note 
# convert to data table from list for plotting with ggmap.
ghcn_pts<-data.frame(ID = ghcn_seattle$spatial$ID,
                     lat = ghcn_seattle$spatial@coords[,2],
                     long = ghcn_seattle$spatial@coords[,1])

ghcn_pts<-mutate(ghcn_pts, ID = as.character(as.factor(ID)))


# Plot the points again.
ggmap(seattle) + geom_point(aes(x=long, y=lat),
                           size=3, color='red', data=ghcn_pts)

ghcn_dat<-ghcn_seattle$tabular


# Pull out the PRCP (precipitation) tables. If there were other
# climate elements we could access these in a 
# similar way (i.e. ghcn_dat, "[[", "TMIN")
prcp<-lapply(ghcn_dat, "[[", "PRCP")


# Dissolve the list into a single data.frame and add a station id
# from the row names (then replace row names with a number sequence)
prcp<-do.call("rbind", prcp)
prcp$station<-substring(row.names(prcp),1,11)
row.names(prcp)<-1:nrow(prcp)

# I'm going to convert to "long" format for easier calculations.
# The gather function is from the tidyr package
prcp<-gather(prcp, day, precip, -station, -YEAR, -MONTH)


# Add the date
prcp<-mutate(prcp,
             MONTH = stringr::str_pad(MONTH, 2, side="left", pad="0"),
             day = stringr::str_pad(gsub("D", "", day), 2,
                                    side="left", pad="0"),
             date = as.Date(paste(YEAR, MONTH, day, sep="-")))


# Exclude non NA values and limit to current years
prcp<-dplyr::filter(prcp, !is.na(date), !is.na(precip),
                    date > as.Date("2008-01-01"),
                    date < as.Date("2016-01-01"))

# Do a count by station using dplyr functions
cnt<-group_by(prcp, station) %>% 
  summarize(cnt = n()) %>% arrange(cnt) %>% 
  mutate(station = factor(station, levels = station))


ggplot(cnt, aes(station, cnt)) + 
  geom_bar(stat="identity", fill="cadetblue") + 
  coord_flip() +
  geom_hline(yintercept=2000) +
  ggtitle("Count of days with data at climate stations")