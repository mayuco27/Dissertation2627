# Load required libraries
library(sf)
library(rnaturalearth)
library(dplyr)

# 1. Download the Atlantic Ocean polygon data
ocean_features <- ne_download(scale = 110, type = "geography_marine_polys", category = "physical", returnclass = "sf")
atlantic_polygon <- ocean_features %>% 
  filter(grepl("Atlantic Ocean", name, ignore.case = TRUE))

# 2. Read your CSV file
locations <- read_csv("~/Dissertation2627/Data/2016 Traits/locations.csv")
View(locations)

# 2.1. FIX: Remove rows where longitude or latitude is missing (NA)
locations_clean <- locations %>%
  filter(!is.na(longitude) & !is.na(latitude)) 

# 3. Convert the CLEANED data into an R spatial object
my_points_sf <- st_as_sf(
  locations_clean,                          # Use the cleaned data frame here
  coords = c("longitude", "latitude"), 
  crs = 4326, 
  remove = FALSE  
)

# 4. Filter coordinates inside the Atlantic Ocean
atlantic_points_sf <- st_intersection(my_points_sf, atlantic_polygon)

# 5. Save as a new CSV
result_df <- atlantic_points_sf %>% select(1:note)

summarytable <- summary(result_df[, c("latitude", "longitude")])
print(summarytable)

write.csv(result_df, "atlantic_coordinates.csv", row.names = FALSE)
cat("Success! Filtered coordinates have been saved.\n")



