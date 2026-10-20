library(sf)
library(tidyverse)

# 1. Load the datasets Force R to read ALL columns as character strings first to eliminate parsin
octo_data_raw <- read_csv("OctocoralTraits_v2_2.csv", col_types = cols(.default = "c"))
atlantic_coords <- read_csv("/Users/may_urico_/Dissertation2627/Data/2016 Traits/atlantic_coordinates.csv")

# 2. 
# 2.a. First, convert text columns into true numeric numbers
octo_numeric <- octo_data_raw %>%
  mutate(
    latitude = as.numeric(latitude),
    longitude = as.numeric(longitude)
  )

# 2.b. Drop the rows that failed the numeric conversion (the text/empty anomalies)
octo_clean <- octo_numeric %>%
  drop_na(latitude, longitude)

# 2.bb. Check how many objects are dropped == 12079
nrow(octo_data_raw) - nrow(octo_clean)

# ---
# Vast majority of data dropped due to "global estimate" on locations column, also don't have latitude and longitude coordinates. 
# -> will shift to downloading OBIS species occurence data 



# 2.c. Now build your spatial object safely using WGS84 (CRS 4326 is standard for WGS84)
octo_sf <- octo_clean %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# 2.d. Convert Octocoral data into a spatial points object (WGS84 projection)
octo_sf <- octo_sf %>%
  drop_na(latitude, longitude) %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4001, remove = FALSE)

# 3. Convert your Atlantic coordinates into a spatial polygon or point grid
# Option A: If your CSV is a bounding grid of points, turn it into points:
atlantic_sf <- atlantic_coords %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4001)

# Option B: If your CSV marks boundary vertices, you can build a bounding box or polygon:
# atlantic_poly <- atlantic_sf %>% st_combine() %>% st_cast("POLYGON")

# 4. Filter Octocoral data points that intersect/overlap your Atlantic points
# (Using st_intersects or st_is_within_distance if you want to allow a tiny buffer tolerance)
matched_spatial_indices <- st_intersects(octo_sf, atlantic_sf, sparse = FALSE)

# Extract rows where any intersection evaluates to TRUE
atlantic_octo_final <- octo_data[rowSums(matched_spatial_indices) > 0, ]

# 5. Save the output
write_csv(atlantic_octo_final, "atlantic_octocoraltraits_spatial.csv")
