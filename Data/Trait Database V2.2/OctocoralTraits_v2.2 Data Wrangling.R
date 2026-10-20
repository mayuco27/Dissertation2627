library(readr)
library(dplyr)
library(stringr)


# 0. Structure of whole dataset
traits <- read_csv("OctocoralTraits_v2_2.csv")
names(traits) # Column names
glimpse(traits) # Data structure

# 0.1. Trait Counts 
traits %>% 
  count(trait_name, sort = TRUE) 

# 0.2. Species Counts
traits %>%
  summarise(
    n_species = n_distinct(specie_name),
    n_locations = n_distinct(location_id),
    n_traits = n_distinct(trait_name)
  )

# 0.3. Cleaning Dataset

cleaned_traits <- traits %>%
  # Remove rows that don't have an actual trait value
  filter(!is.na(trait_name)) %>% 
  # Standardize text columns
  mutate(
    trait_name = tolower(trait_name),
    specie_name = trimws(specie_name) # Removes accidental leading/trailing spaces
  )

# Check missing coordinates (only 27943 records out of 148493 records lacking either or both of the coordinates.)
# Global Estimates: Global estimate is special location reserved for species-level traits. Such traits include individual-level traits that do not vary among individuals of a species over ecological time-scales (e.g., sexual system), as well as emergent characteristics of individuals within a species (such as geographical range size and greatest depth observed). Species-level traits do not have contextual traits because, by definition, they apply to all individuals, globally.

cleaned_traits %>%
  summarise(
    total_rows = n(),
    missing_lat = sum(is.na(latitude)),
    missing_lon = sum(is.na(longitude)),
    missing_either = sum(
      is.na(latitude) | is.na(longitude)
    )
  )

# Check duplicate coordinates in atlantic_coords
atlantic_coords %>%
  count(latitude, longitude) %>%
  filter(n > 1) %>%
  arrange(desc(n))

# Check whether coordinates are valid
cleaned_traits %>%
  filter(
    !is.na(latitude),
    !is.na(longitude)
  ) %>%
  summarise(
    invalid_lat = sum(latitude < -90 | latitude > 90),
    invalid_lon = sum(longitude < -180 | longitude > 180)
  )

# Verify the cleaned data
summary(cleaned_traits)

# ----------------------------------------------
# 1. Spatial filtering (Isolate Atlantic ocean)
# Identify species-level records
global_traits <- cleaned_traits %>%
  filter(
    tolower(trimws(location_name)) == "global estimate"
  )

# Keep records with actual location information separately
location_traits <- cleaned_traits %>%
  filter(
    is.na(location_name) |
      tolower(trimws(location_name)) != "global estimate"
  )

# Compare the two groups
tibble(
  dataset = c("Global estimate", "Other locations"),
  n_records = c(nrow(global_traits), nrow(location_traits)),
  n_species = c(
    n_distinct(global_traits$specie_name, na.rm = TRUE),
    n_distinct(location_traits$specie_name, na.rm = TRUE)
  )
)

atlantic_traits <- cleaned_traits %>%
  inner_join(atlantic_coords, by = c("latitude", "longitude"))


