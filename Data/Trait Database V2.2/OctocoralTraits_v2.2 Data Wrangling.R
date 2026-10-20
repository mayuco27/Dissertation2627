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

# Verify the cleaned data
summary(cleaned_traits)

# ----------------------------------------------
# 1. Spatial filtering (Isolate Atlantic ocean)
atlantic_traits <- cleaned_traits %>%
  inner_join(atlantic_coords, by = c("latitude", "longitude"))


