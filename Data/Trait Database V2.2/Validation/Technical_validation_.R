###Code used to validate OctocoralTraits v2.2 data released in the manuscript:
##Gómez-Gras et al.The Octocoral Trait Database: a global database of trait information for octocoral species 
#The code includes:
#a) Detection of outliers and extreme outliers for quantitative traits, which were there manually 
#evaluated to check for correctness.
#b) Detection of wrong category names in categorical traits
#c) Detection of duplicated rows
#d) Ensuring that each observation_id is associated to only one location, resource and species

#We set the directory 
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

##We load the dataset
OTD<-read.csv("../OctocoralTraits_v2_2.csv")
df<-OTD

##We remove contextual traits
df<-OTD[OTD$trait_category!="Contextual",]

## a) Detection of outliers and Extreme outliers in quantitative (numerical) traits

#Libraries
library("dplyr")
library("ggplot2")
library("gridExtra")
library("grid")
library("ggpubr")

#Function to identify outliers
identify_outliers <- function(values) {
  q1 <- quantile(values, 0.25)  # First quartile
  q3 <- quantile(values, 0.75)  # Third quartile
  iqr <- IQR(values)  # Interquartile range
  lower_bound <- q1 - 1.5 * iqr  # Lower bound for outliers
  upper_bound <- q3 + 1.5 * iqr  # Upper bound for outliers
  return(values < lower_bound | values > upper_bound)  # Return a logical vector of outliers
}

#Function to identify extreme outliers
identify_extreme_outliers <- function(values) {
  q1 <- quantile(values, 0.25)  # First quartile
  q3 <- quantile(values, 0.75)  # Third quartile
  iqr <- IQR(values)  # Interquartile range
  lower_bound <- q1 - 3 * iqr  # Lower bound for extreme outliers
  upper_bound <- q3 + 3 * iqr  # Upper bound for extreme outliers
  return(values < lower_bound | values > upper_bound)  # Return a logical vector of extreme outliers
}

##With the following code we will check if each trait is associated to its correct variable type, and whether there are
# traits associated to multiple.

# How This Works:
# get_trait_type():
# - Mixed types check: The function first checks if there are both numeric and non-numeric values for a trait using is_numeric. If both are found, it flags the trait as "mixed (numeric and character/factor)".
# - Binomial detection: If a trait has exactly two unique values, it’s classified as "binomial (factor)".
# - Numeric detection: If all values are numeric, the trait is classified as "numeric".
# - Factor detection: If the trait has repeated non-numeric values (indicating a categorical variable with more than two levels), it's classified as a factor.
# - Character detection: If none of the conditions are met, the trait is classified as character.
# 
# tapply(): Applies this function to each trait by grouping values in the value column according to the trait names.
#
# result_df: The final result is stored in a dataframe with each trait and its corresponding value type.

#Actual method

# Define a function to check for mixed types and classify the trait
get_trait_type <- function(x) {
  unique_values <- unique(x)
  
  # Step 1: Check for mixed types
  is_numeric <- !is.na(suppressWarnings(as.numeric(x)))  # Detect numeric values
  if (any(is_numeric) && any(!is_numeric)) {
    return("mixed (numeric and character/factor)")
  }
  
  # Step 2: Handle consistent types
  # Binomial check: If there are exactly two unique values, it's binomial (factor)
  if (length(unique_values) == 2) {
    return("binomial (factor)")
  }
  
  # Numeric check: If all values are numeric
  if (all(!is.na(suppressWarnings(as.numeric(x))))) {
    return("numeric")
  }
  
  # Factor check: If the trait has repeated non-numeric values
  if (length(unique_values) < length(x)) {
    return("factor")
  }
  
  # If none of the above, classify as character
  return("character")
}

# Apply the function to each trait to get the type of values
types_per_trait <- tapply(df$value, df$trait_name, get_trait_type)

# Convert the result into a dataframe
result_df <- data.frame(
  trait = names(types_per_trait),
  value_type = as.character(types_per_trait)
)

##We replace the value of colony contractability, since it is indeed binomial
result_df <- result_df %>%
  mutate(value = ifelse(trait == "Colony contractability" & value_type == "numeric", 
                        "binomial (factor)", 
                        value_type))


#### Make subsets of numeric traits
numeric_traits<-result_df[result_df$value_type=="numeric",]
numeric_traits$trait<-as.factor(numeric_traits$trait)

trait_names<-numeric_traits$trait
for (trait in trait_names) {
  # Create a variable name for the subset
  subset_name <- make.names(trait)  # Create a valid R variable name
  
  # Create the subset and assign it to a variable
  assign(subset_name, df[df$trait_name == trait, ])
}


#### Make subsets of categorical traits
categorical_traits<-result_df[result_df$value_type=="factor",]
categorical_traits$trait<-as.factor(categorical_traits$trait)

trait_names<-categorical_traits$trait
for (trait in trait_names) {
  # Create a variable name for the subset
  subset_name <- make.names(trait)  # Create a valid R variable name
  # Create the subset and assign it to a variable
  assign(subset_name, df[df$trait_name == trait, ])}
  

### For each quantitative trait, we standardize values to common unit and detect outliers and 
### extreme outliers as follows:

#1) Age at maturity
Age.at.maturity$value=as.numeric(Age.at.maturity$value)

#We check standard units and standardized in case needed
#in this case, month to years
Age.at.maturity$standard_unit=as.factor(Age.at.maturity$standard_unit)
levels(Age.at.maturity$standard_unit)
Age.at.maturity$standardized_unit = "years"

# Create the standardized_value column
Age.at.maturity <- Age.at.maturity %>%
  mutate(standardized_value = ifelse(standard_unit == "months", value / 12, value))

# Identify outliers 
Age.at.maturity <- Age.at.maturity %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Age.at.maturity <- Age.at.maturity %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

#Merge both
Age.at.maturity <- Age.at.maturity %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Age.at.maturity <- Age.at.maturity %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))


# Create a table of outliers 
outlier_Age.at.maturity <- Age.at.maturity %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Age.at.maturity)

num_outliers1 <- nrow(Age.at.maturity %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Age.at.maturity %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Age.at.maturity %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p1 <- ggplot(Age.at.maturity, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Age.at.maturity$trait_name), "- N:", nrow(Age.at.maturity)),
       x = "Trait Name",
       y = "Years",
       color = "Outlier Status") +
  theme_minimal()
p1

#2) Axis radial growth rate

Axis.radial.growth.rate$value=as.numeric(Axis.radial.growth.rate$value)

#We check standard units and standardized in case needed
#in this case, month to years
Axis.radial.growth.rate$standard_unit=as.factor(Axis.radial.growth.rate$standard_unit)
levels(Axis.radial.growth.rate$standard_unit)
Axis.radial.growth.rate$standardized_unit = "mm yr^-1"

# Create the standardized_value column
Axis.radial.growth.rate <- Axis.radial.growth.rate %>%
  mutate(standardized_value = case_when(
    standard_unit == "mm month^-1" ~ value * 12,       # Convert to mm yr^-1
    standard_unit == "μm yr^-1" ~ value / 1000,        # Convert to mm yr^-1
    TRUE ~ value  # Already in mm yr^-1
  ))

# Identify outliers 
Axis.radial.growth.rate <- Axis.radial.growth.rate %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Axis.radial.growth.rate <- Axis.radial.growth.rate %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

#Merge both
Axis.radial.growth.rate <- Axis.radial.growth.rate %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Axis.radial.growth.rate <- Axis.radial.growth.rate %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))


# Create a table of outliers 
outlier_Axis.radial.growth.rate <- Axis.radial.growth.rate %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Axis.radial.growth.rate)

num_outliers1 <- nrow(Axis.radial.growth.rate %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Axis.radial.growth.rate %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Axis.radial.growth.rate %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p2 <- ggplot(Axis.radial.growth.rate, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Axis.radial.growth.rate$trait_name), "- N:", nrow(Axis.radial.growth.rate)),
       x = "Trait Name",
       y = "mm y^-1",
       color = "Outlier Status") +
  theme_minimal()
p2


#3) Branch diameter

Branch.diameter$value=as.numeric(Branch.diameter$value)

#We check standard units and standardized in case needed
#in this case, month to years
Branch.diameter$standard_unit=as.factor(Branch.diameter$standard_unit)
levels(Branch.diameter$standard_unit)
Branch.diameter$standardized_unit = "mm"

# No need to create the standardized_value since they are all mm
Branch.diameter <- Branch.diameter %>%
  mutate(standardized_value = case_when(
    TRUE ~ value  # Already in mm
  ))


# Identify outliers 
Branch.diameter <- Branch.diameter %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Branch.diameter <- Branch.diameter %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

#Merge both
Branch.diameter <- Branch.diameter %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Branch.diameter <- Branch.diameter %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))


# Create a table of outliers 
outlier_Branch.diameter <- Branch.diameter %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Branch.diameter)

num_outliers1 <- nrow(Branch.diameter %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Branch.diameter %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Branch.diameter %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p3 <- ggplot(Branch.diameter, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Branch.diameter$trait_name), "- N:", nrow(Branch.diameter)),
       x = "Trait Name",
       y = "mm",
       color = "Outlier Status") +
  theme_minimal()
p3



#4) Ca.Mg.ratio.in.sclerites

Ca.Mg.ratio.in.sclerites$value=as.numeric(Ca.Mg.ratio.in.sclerites$value)

#We check standard units and standardized in case needed
#in this case, month to years
Ca.Mg.ratio.in.sclerites$standard_unit=as.factor(Ca.Mg.ratio.in.sclerites$standard_unit)
levels(Ca.Mg.ratio.in.sclerites$standard_unit)
Ca.Mg.ratio.in.sclerites$standardized_unit = "mol:mol"

# No need to create the standardized_value since they are all mm
Ca.Mg.ratio.in.sclerites <- Ca.Mg.ratio.in.sclerites %>%
  mutate(standardized_value = case_when(
    TRUE ~ value  # Already in mol:mol
  ))


# Identify outliers 
Ca.Mg.ratio.in.sclerites <- Ca.Mg.ratio.in.sclerites %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Ca.Mg.ratio.in.sclerites <- Ca.Mg.ratio.in.sclerites %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

#Merge both
Ca.Mg.ratio.in.sclerites <- Ca.Mg.ratio.in.sclerites %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Ca.Mg.ratio.in.sclerites <- Ca.Mg.ratio.in.sclerites %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))


# Create a table of outliers 
outlier_Ca.Mg.ratio.in.sclerites <- Ca.Mg.ratio.in.sclerites %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Branch.diameter)

num_outliers1 <- nrow(Ca.Mg.ratio.in.sclerites %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Ca.Mg.ratio.in.sclerites %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Ca.Mg.ratio.in.sclerites %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p4 <- ggplot(Ca.Mg.ratio.in.sclerites, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Ca.Mg.ratio.in.sclerites$trait_name), "- N:", nrow(Ca.Mg.ratio.in.sclerites)),
       x = "Trait Name",
       y = "mol:mol",
       color = "Outlier Status") +
  theme_minimal()
p4



###5) Calcification.rate

Calcification.rate$value=as.numeric(Calcification.rate$value)

#We check standard units and standardized in case needed
#in this case, month to years
Calcification.rate$standard_unit=as.factor(Calcification.rate$standard_unit)
levels(Calcification.rate$standard_unit)

Calcification.rate_a<-Calcification.rate[Calcification.rate$standard_unit ==" µgmol CaCO3 cm^-2 d^-1",]
Calcification.rate_b<-Calcification.rate[Calcification.rate$standard_unit =="% wt month^-1",]
Calcification.rate_c<-Calcification.rate[Calcification.rate$standard_unit =="g cm ^-2 yr^-1",]
Calcification.rate_d<-Calcification.rate[Calcification.rate$standard_unit =="nmol mg_protein^-1 d^-1",]
Calcification.rate_e<-Calcification.rate[Calcification.rate$standard_unit =="µg Ca mg N^-1 h^-1",]

Calcification.rate_a$standardized_unit = " µgmol CaCO3 cm^-2 d^-1"

# No need to create the standardized_value since they are all mm
Calcification.rate_a <- Calcification.rate_a %>%
  mutate(standardized_value = case_when(
    TRUE ~ value  # Already in " µgmol CaCO3 cm^-2 d^-1"
  ))

# Identify outliers 
Calcification.rate_a <- Calcification.rate_a %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Calcification.rate_a <- Calcification.rate_a %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

#Merge both
Calcification.rate_a <- Calcification.rate_a %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Calcification.rate_a <- Calcification.rate_a %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))


# Create a table of outliers 
outlier_Calcification.rate_a <- Calcification.rate_a %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Calcification.rate_a)

num_outliers1 <- nrow(Calcification.rate_a %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Calcification.rate_a %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Calcification.rate_a %>% filter(outliers == "Not Outlier"))




# Create the scatter plot
p5 <- ggplot(Calcification.rate_a, aes(x = trait_name, y = as.numeric(standardized_value), color = outlier)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "red"), labels = c("Not Outlier", paste("Outlier (N =", nrow(outlier_Calcification.rate_a), ")"))) +
  labs(title = paste("Plot for", unique(Calcification.rate_a$trait_name),"- N:", nrow(Calcification.rate_a)),
       x = "Trait Name",
       y = " µgmol CaCO3 cm^-2 d^-1",
       color = "Outlier Status") +
  theme_minimal()
p5

# There were not outliers, so there is no need for table 

#6) Calyx height

Calyx.height$value=as.numeric(Calyx.height$value)

#We check standard units and standardized in case needed
#in this case, month to years
Calyx.height$standard_unit=as.factor(Calyx.height$standard_unit)
levels(Calyx.height$standard_unit)
Calyx.height$standardized_unit = "mm"

# Create the standardized_value column
Calyx.height <- Calyx.height %>%
  mutate(standardized_value = case_when(
    standard_unit == "cm" ~ value * 10,       # Convert to mm
    TRUE ~ value  # Already in mm
  ))

# Identify outliers 
Calyx.height <- Calyx.height %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Calyx.height <- Calyx.height %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Calyx.height <- Calyx.height %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Calyx.height <- Calyx.height %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Calyx.height <- Calyx.height %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Calyx.height)

num_outliers1 <- nrow(Calyx.height %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Calyx.height %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Calyx.height %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p6 <- ggplot(Calyx.height, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Calyx.height$trait_name), "- N:", nrow(Calyx.height)),
       x = "Trait Name",
       y = "mm",
       color = "Outlier Status") +
  theme_minimal()
p6

#7) Calyx width

Calyx.width$value=as.numeric(Calyx.width$value)

#We check standard units and standardized in case needed
#in this case, month to years
Calyx.width$standard_unit=as.factor(Calyx.width$standard_unit)
levels(Calyx.width$standard_unit)
Calyx.width$standardized_unit = "mm"

# Create the standardized_value column
Calyx.width <- Calyx.width %>%
  mutate(standardized_value = case_when(
    standard_unit == "µm" ~ value / 1000,       # Convert to mm
    TRUE ~ value  # Already in mm
  ))

# Identify outliers 
Calyx.width <- Calyx.width %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Calyx.width <- Calyx.width %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Calyx.width <- Calyx.width %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Calyx.width <- Calyx.width %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Calyx.width <- Calyx.width %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Calyx.width)

num_outliers1 <- nrow(Calyx.width %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Calyx.width %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Calyx.width %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p7 <- ggplot(Calyx.width, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Calyx.width$trait_name), "- N:", nrow(Calyx.width)),
       x = "Trait Name",
       y = "mm",
       color = "Outlier Status") +
  theme_minimal()
p7
##8) Capture.rate

Capture.rate$value=as.numeric(Capture.rate$value)

#We check standard units and standardized in case needed
#in this case, month to years
Capture.rate$standard_unit=as.factor(Capture.rate$standard_unit)
levels(Capture.rate$standard_unit)

##There are plenty of valid standards with few data so we will simply confirm in original sources that the data has correctly been extracted

#9) Colony age

Colony.age$value=as.numeric(Colony.age$value)

#We check standard units and standardized in case needed
#in this case, month to years
Colony.age$standard_unit=as.factor(Colony.age$standard_unit)
levels(Colony.age$standard_unit)
Colony.age$standardized_unit = "years"

# No need to standardize since all is in years. 
Colony.age <- Colony.age %>%
  mutate(standardized_value = case_when(
    TRUE ~ value  # Already in years
  ))

# Identify outliers 
Colony.age <- Colony.age %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Colony.age <- Colony.age %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Colony.age <- Colony.age %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Colony.age <- Colony.age %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Colony.age <- Colony.age %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Colony.age)

num_outliers1 <- nrow(Colony.age %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Colony.age %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Colony.age %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p9 <- ggplot(Colony.age, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Colony.age$trait_name), "- N:", nrow(Colony.age)),
       x = "Trait Name",
       y = "years",
       color = "Outlier Status") +
  theme_minimal()
p9

#10) Colony area

Colony.area$value=as.numeric(Colony.area$value)

#We check standard units and standardized in case needed
#in this case, month to years
Colony.area$standard_unit=as.factor(Colony.area$standard_unit)
levels(Colony.area$standard_unit)
Colony.area$standardized_unit = "mm^2"

# Create the standardized_value column
Colony.area <- Colony.area %>%
  mutate(standardized_value = case_when(
    standard_unit == "cm^2" ~ value * 100,       # Convert to mm^2
    TRUE ~ value  # Already in mm^2
  ))

# Identify outliers 
Colony.area <- Colony.area %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Colony.area <- Colony.area %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Colony.area <- Colony.area %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Colony.area <- Colony.area %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Colony.area <- Colony.area %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Colony.area)

num_outliers1 <- nrow(Colony.area %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Colony.area %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Colony.area %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p10 <- ggplot(Colony.area, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Colony.area$trait_name), "- N:", nrow(Colony.area)),
       x = "Trait Name",
       y = "cm^2",
       color = "Outlier Status") +
  theme_minimal()
p10


#11) Colony fecundity

Colony.fecundity$value=as.numeric(Colony.fecundity$value)

#We check standard units and standardized in case needed
#in this case, month to years
Colony.fecundity$standard_unit=as.factor(Colony.fecundity$standard_unit)
levels(Colony.fecundity$standard_unit)
#There is only one standard
Colony.fecundity$standardized_unit = "oocytes colony^-1"

# Create the standardized_value column
Colony.fecundity <- Colony.fecundity %>%
  mutate(standardized_value = case_when(
    TRUE ~ value  # Already in oocytes colony^-1
  ))
# Identify outliers 
Colony.fecundity <- Colony.fecundity %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Colony.fecundity <- Colony.fecundity %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Colony.fecundity <- Colony.fecundity %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Colony.fecundity <- Colony.fecundity %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Colony.fecundity <- Colony.fecundity %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Colony.fecundity)

num_outliers1 <- nrow(Colony.fecundity %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Colony.fecundity %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Colony.fecundity %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p11 <- ggplot(Colony.fecundity, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Colony.fecundity$trait_name), "- N:", nrow(Colony.fecundity)),
       x = "Trait Name",
       y = "oocytes colony^-1",
       color = "Outlier Status") +
  theme_minimal()
p11

#12) Colony height

Colony.height$value=as.numeric(Colony.height$value)

#We check standard units and standardized in case needed
#in this case, month to years
Colony.height$standard_unit=as.factor(Colony.height$standard_unit)
levels(Colony.height$standard_unit)
Colony.height$standardized_unit = "cm"

# Create the standardized_value column
Colony.height <- Colony.height %>%
  mutate(standardized_value = case_when(
    standard_unit == "cm" ~ value,                # Already in cm
    standard_unit == "ft" ~ value * 30.48,        # Convert feet to cm
    standard_unit == "in." ~ value * 2.54,        # Convert inches to cm
    standard_unit == "m" ~ value * 100,           # Convert meters to cm
    standard_unit == "mm" ~ value / 10,           # Convert mm to cm
    TRUE ~ NA_real_                                # Assign NA for unsupported units
  ))

# Identify outliers 
Colony.height <- Colony.height %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Colony.height <- Colony.height %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Colony.height <- Colony.height %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Colony.height <- Colony.height %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Colony.height <- Colony.height %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Colony.height)

num_outliers1 <- nrow(Colony.height %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Colony.height %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Colony.height %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p12 <- ggplot(Colony.height, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Colony.height$trait_name), "- N:", nrow(Colony.height)),
       x = "Trait Name",
       y = "cm",
       color = "Outlier Status") +
  theme_minimal()
p12


#13) Colony width

Colony.width$value=as.numeric(Colony.width$value)

#We check standard units and standardized in case needed
#in this case, month to years
Colony.width$standard_unit=as.factor(Colony.width$standard_unit)
levels(Colony.width$standard_unit)
Colony.width$standardized_unit = "cm"

# Create the standardized_value column
Colony.width <- Colony.width %>%
  mutate(standardized_value = case_when(
    standard_unit == "cm" ~ value,                # Already in cm
      standard_unit == "m" ~ value * 100,           # Convert meters to cm
    standard_unit == "mm" ~ value / 10,           # Convert mm to cm
    TRUE ~ NA_real_                                # Assign NA for unsupported units
  ))

# Identify outliers 
Colony.width <- Colony.width %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Colony.width <- Colony.width %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Colony.width <- Colony.width %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Colony.width <- Colony.width %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Colony.width <- Colony.width %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Colony.width)

num_outliers1 <- nrow(Colony.width %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Colony.width %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Colony.width %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p13 <- ggplot(Colony.width, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Colony.width$trait_name), "- N:", nrow(Colony.width)),
       x = "Trait Name",
       y = "cm",
       color = "Outlier Status") +
  theme_minimal()
p13


#14) Depth.lower

Depth.lower$value=as.numeric(Depth.lower$value)

#We check standard units and standardized in case needed
#in this case, month to years
Depth.lower$standard_unit=as.factor(Depth.lower$standard_unit)
levels(Depth.lower$standard_unit)
Depth.lower$standardized_unit = "m"

# Create the standardized_value column
Depth.lower <- Depth.lower %>%
  mutate(standardized_value = case_when(
    standard_unit == "m" ~ value,                   # Already in meters
    standard_unit == "ft" ~ value * 0.3048,         # Convert feet to meters
    standard_unit == "fm" ~ value * 0.000508,       # Convert fathoms to meters (1 fm = 1.8288 m)
    TRUE ~ NA_real_                                  # Assign NA for unsupported units
  ))

# Identify outliers 
Depth.lower <- Depth.lower %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Depth.lower <- Depth.lower %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Depth.lower <- Depth.lower %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Depth.lower <- Depth.lower %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Depth.lower <- Depth.lower %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Depth.lower)

num_outliers1 <- nrow(Depth.lower %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Depth.lower %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Depth.lower %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p14 <- ggplot(Depth.lower, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Depth.lower$trait_name), "- N:", nrow(Depth.lower)),
       x = "Trait Name",
       y = "m",
       color = "Outlier Status") +
  theme_minimal()

p14

#15) Depth.upper

Depth.upper$value=as.numeric(Depth.upper$value)

#We check standard units and standardized in case needed
#in this case, month to years
Depth.upper$standard_unit=as.factor(Depth.upper$standard_unit)
levels(Depth.upper$standard_unit)
Depth.upper$standardized_unit = "m"

# Create the standardized_value column
Depth.upper <- Depth.upper %>%
  mutate(standardized_value = case_when(
    standard_unit == "m" ~ value,                   # Already in meters
    standard_unit == "ft" ~ value * 0.3048,         # Convert feet to meters
    standard_unit == "fm" ~ value * 0.000508,       # Convert fathoms to meters (1 fm = 1.8288 m)
    TRUE ~ NA_real_                                  # Assign NA for unsupported units
  ))

# Identify outliers 
Depth.upper <- Depth.upper %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Depth.upper <- Depth.upper %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Depth.upper <- Depth.upper %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Depth.upper <- Depth.upper %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Depth.upper <- Depth.upper %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Depth.upper)

num_outliers1 <- nrow(Depth.upper %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Depth.upper %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Depth.upper %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p15 <- ggplot(Depth.upper, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Depth.upper$trait_name), "- N:", nrow(Depth.upper)),
       x = "Trait Name",
       y = "m",
       color = "Outlier Status") +
  theme_minimal()

p15

#16) Diameter of axis skeleton

Diameter.of.axis.skeleton$value=as.numeric(Diameter.of.axis.skeleton$value)

#We check standard units and standardized in case needed
#in this case, month to years
Diameter.of.axis.skeleton$standard_unit=as.factor(Diameter.of.axis.skeleton$standard_unit)
levels(Diameter.of.axis.skeleton$standard_unit)
Diameter.of.axis.skeleton$standardized_unit = "mm"

# Create the standardized_value column
Diameter.of.axis.skeleton <- Diameter.of.axis.skeleton %>%
  mutate(standardized_value = case_when(
    standard_unit == "cm" ~ value * 10,         # Convert cm to mm
    standard_unit == "mm" ~ value))

# Identify outliers 
Diameter.of.axis.skeleton <- Diameter.of.axis.skeleton %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Diameter.of.axis.skeleton <- Diameter.of.axis.skeleton %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Diameter.of.axis.skeleton <- Diameter.of.axis.skeleton %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Diameter.of.axis.skeleton <- Diameter.of.axis.skeleton %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Diameter.of.axis.skeleton <- Diameter.of.axis.skeleton %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Diameter.of.axis.skeleton)

num_outliers1 <- nrow(Diameter.of.axis.skeleton %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Diameter.of.axis.skeleton %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Diameter.of.axis.skeleton %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p16 <- ggplot(Diameter.of.axis.skeleton, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Diameter.of.axis.skeleton$trait_name), "- N:", nrow(Diameter.of.axis.skeleton)),
       x = "Trait Name",
       y = "mm",
       color = "Outlier Status") +
  theme_minimal()

p16

#17) Distance between polyps

Distance.between.polyps$value=as.numeric(Distance.between.polyps$value)

#We check standard units and standardized in case needed
#in this case, month to years
Distance.between.polyps$standard_unit=as.factor(Distance.between.polyps$standard_unit)
levels(Distance.between.polyps$standard_unit)
Distance.between.polyps$standardized_unit = "mm"

# Create the standardized_value column
Distance.between.polyps <- Distance.between.polyps %>%
  mutate(standardized_value = case_when(
    standard_unit == "mm" ~ value))

# Identify outliers 
Distance.between.polyps <- Distance.between.polyps %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Distance.between.polyps <- Distance.between.polyps %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Distance.between.polyps <- Distance.between.polyps %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Distance.between.polyps <- Distance.between.polyps %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Distance.between.polyps <- Distance.between.polyps %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

print(outlier_Distance.between.polyps)

num_outliers1 <- nrow(Distance.between.polyps %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Distance.between.polyps %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Distance.between.polyps %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p17 <- ggplot(Distance.between.polyps, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Distance.between.polyps$trait_name), "- N:", nrow(Distance.between.polyps)),
       x = "Trait Name",
       y = "mm",
       color = "Outlier Status") +
  theme_minimal()

p17

#18) Eggs per area

Eggs.per.area$value=as.numeric(Eggs.per.area$value)

#We check standard units and standardized in case needed
#in this case, month to years
Eggs.per.area$standard_unit=as.factor(Eggs.per.area$standard_unit)
levels(Eggs.per.area$standard_unit)

#There are few trait observations and with different standards with are not comparable. Thus
#we ensured the accuracy of the incorporated data by manually reviewing them in the sources.



#19) Fractal dimension
Fractal.dimension$value=as.numeric(Fractal.dimension$value)

#We check standard units and standardized in case needed
#in this case, month to years
Fractal.dimension$standard_unit=as.factor(Fractal.dimension$standard_unit)
levels(Fractal.dimension$standard_unit)
Fractal.dimension$standardized_unit = "x:y"

# Create the standardized_value column
Fractal.dimension <- Fractal.dimension %>%
  mutate(standardized_value = case_when(
    standard_unit == "x:y" ~ value))

# Identify outliers 
Fractal.dimension <- Fractal.dimension %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Fractal.dimension <- Fractal.dimension %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Fractal.dimension <- Fractal.dimension %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Fractal.dimension <- Fractal.dimension %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Fractal.dimension <- Fractal.dimension %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

##There are not outliers

num_outliers1 <- nrow(Fractal.dimension %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Fractal.dimension %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Fractal.dimension %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p19 <- ggplot(Fractal.dimension, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Fractal.dimension$trait_name), "- N:", nrow(Fractal.dimension)),
       x = "Trait Name",
       y = "x:y",
       color = "Outlier Status") +
  theme_minimal()

p19

#20) Gross photosynthesis

Gross.photosynthesis$value=as.numeric(Gross.photosynthesis$value)

#We check standard units and standardized in case needed
#in this case, month to years
Gross.photosynthesis$standard_unit=as.factor(Gross.photosynthesis$standard_unit)
levels(Gross.photosynthesis$standard_unit)
Gross.photosynthesis$standardized_unit = "mm"

##There are few observations and from two different standards that cannot be grouped. Thus, we ensured the accuracy of the data by simply
#looking at the original sources.


#21) Growth.rate

Growth.rate$value=as.numeric(Growth.rate$value)

#We check standard units and standardized in case needed
#in this case, month to years
Growth.rate$standard_unit=as.factor(Growth.rate$standard_unit)
levels(Growth.rate$standard_unit)
Growth.rate$standardized_unit = "mm year^-1"

#We first remove from the list standards tha are not in length over time units, since those cannot be compared
removed_units <- c("BFSA% yr^-1", "g year^-1", "polyps d^-1","cm^2 month^-1","mm^2 per 4 months^-1","polyps year^-1","cm^2 230d^-1","cm^2 yr^-1","mm^2 d^-1")
Growth.rate_b <- Growth.rate[!Growth.rate$standard_unit %in% removed_units,]

# Create the standardized_value column
Growth.rate_b <- Growth.rate_b %>%
  mutate(standardized_value = case_when(
    standard_unit == "cm yr^-1" ~ value * 10,                   # Convert cm to mm
    standard_unit == "cm month^-1" ~ value * 10 * 12,             # Convert cm to mm and then to year^-1
    standard_unit == "cm week^-1" ~ value * 10 * 52,              # Convert cm to mm and then to year^-1
    standard_unit == "cm 6months^-1" ~ value * 10 * 2,            # Convert cm to mm and then to year^-1
    standard_unit == "cm 11months^-1" ~ value * 10 * (11/12),     # Convert cm to mm and then to year^-1
    standard_unit == "cm 230d^-1" ~ value * 10 * (230/365),       # Convert cm to mm and then to year^-1
    standard_unit == "mm yr^-1" ~ value,                         # Already in mm year^-1
    standard_unit == "mm month^-1" ~ value * 12,                  # Convert mm month^-1 to mm year^-1
    standard_unit == "mm week^-1" ~ value * 52,                   # Convert mm week^-1 to mm year^-1
    standard_unit == "mm d^-1" ~ value * 365,                     # Convert mm d^-1 to mm year^-1
    TRUE ~ NA_real_                                                # Assign NA for unsupported units
  ))


# Identify outliers 
Growth.rate_b <- Growth.rate_b %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Growth.rate_b <- Growth.rate_b %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Growth.rate_b <- Growth.rate_b %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Growth.rate_b <- Growth.rate_b %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Growth_rate_b <- Growth.rate_b %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Growth_rate_b)

# Count the number of outlier categories
num_outliers1 <- nrow(Growth.rate_b %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Growth.rate_b %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Growth.rate_b %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p21 <- ggplot(Growth.rate_b, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Growth.rate_b$trait_name), "- N:", nrow(Growth.rate_b)),
       x = "Trait Name",
       y = "mm/year",
       color = "Outlier Status") +
  theme_minimal()

p21


#22) Larval length

Larval.length$value=as.numeric(Larval.length$value)

#We check standard units and standardized in case needed
#in this case, month to years
Larval.length$standard_unit=as.factor(Larval.length$standard_unit)
levels(Larval.length$standard_unit)
Larval.length$standardized_unit = "mm"

# Create the standardized_value column
Larval.length <- Larval.length %>%
  mutate(standardized_value = case_when(
    standard_unit == "µm" ~ value / 1000,                # Convert micras (µm) to mm
    standard_unit == "mm" ~ value,                         # Already in mm year^-1
      TRUE ~ NA_real_                                                # Assign NA for unsupported units
  ))

# Identify outliers 
Larval.length <- Larval.length %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Larval.length <- Larval.length %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Larval.length <- Larval.length %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Larval.length <- Larval.length %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Larval_length <- Larval.length %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Larval_length)

# Count the number of outlier categories
num_outliers1 <- nrow(Larval.length %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Larval.length %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Larval.length %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p22 <- ggplot(Larval.length, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Larval.length$trait_name), "- N:", nrow(Larval.length)),
       x = "Trait Name",
       y = "mm",
       color = "Outlier Status") +
  theme_minimal()

p22


#23) Larval width

Larval.width$value=as.numeric(Larval.width$value)

#We check standard units and standardized in case needed
#in this case, month to years
Larval.width$standard_unit=as.factor(Larval.width$standard_unit)
levels(Larval.width$standard_unit)
Larval.width$standardized_unit = "mm"

# Create the standardized_value column
Larval.width <- Larval.width %>%
  mutate(standardized_value = case_when(
    standard_unit == "µm" ~ value / 1000,                # Convert micras (µm) to mm
    standard_unit == "mm" ~ value,                         # Already in mm year^-1
    TRUE ~ NA_real_                                                # Assign NA for unsupported units
  ))

# Identify outliers 
Larval.width <- Larval.width %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Larval.width <- Larval.width %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Larval.width <- Larval.width %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Larval.width <- Larval.width %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Larval_width <- Larval.width %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Larval_width)

# Count the number of outlier categories
num_outliers1 <- nrow(Larval.width %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Larval.width %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Larval.width %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p23 <- ggplot(Larval.width, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Larval.width$trait_name), "- N:", nrow(Larval.width)),
       x = "Trait Name",
       y = "mm",
       color = "Outlier Status") +
  theme_minimal()

p23

#24) Longevity

Longevity$value=as.numeric(Longevity$value)

#We check standard units and standardized in case needed
#in this case, month to years
Longevity$standard_unit=as.factor(Longevity$standard_unit)
levels(Longevity$standard_unit)
Longevity$standardized_unit = "years"

# Create the standardized_value column
Longevity <- Longevity %>%
  mutate(standardized_value = case_when(
     standard_unit == "years" ~ value,                         # Already in years
    TRUE ~ NA_real_                                                # Assign NA for unsupported units
  ))

# Identify outliers 
Longevity <- Longevity %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Longevity <- Longevity %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Longevity <- Longevity %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Longevity <- Longevity %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Longevity <- Longevity %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Longevity)

# Count the number of outlier categories
num_outliers1 <- nrow(Longevity %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Longevity %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Longevity %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p24 <- ggplot(Longevity, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Longevity$trait_name), "- N:", nrow(Longevity)),
       x = "Trait Name",
       y = "Years",
       color = "Outlier Status") +
  theme_minimal()

p24

#25) Mg.Ca.ratio.in.sclerites

Mg.Ca.ratio.in.sclerites$value=as.numeric(Mg.Ca.ratio.in.sclerites$value)

#We check standard units and standardized in case needed
#in this case, month to years
Mg.Ca.ratio.in.sclerites$standard_unit=as.factor(Mg.Ca.ratio.in.sclerites$standard_unit)
levels(Mg.Ca.ratio.in.sclerites$standard_unit)
Mg.Ca.ratio.in.sclerites$standardized_unit = "mol:mol"

# We create the standardized_value since they are all mm
Mg.Ca.ratio.in.sclerites<-Mg.Ca.ratio.in.sclerites%>%
mutate(standardized_value = case_when(
  standard_unit == "mmol/mol" ~ value * 1e-3,  # Convert mmol:mol to mol:mol
  standard_unit == "mol:mol" ~ value,          # Already in mol:mol
  TRUE ~ NA_real_                              # Assign NA for unsupported units
))

# Identify outliers 
Mg.Ca.ratio.in.sclerites <- Mg.Ca.ratio.in.sclerites %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Mg.Ca.ratio.in.sclerites <- Mg.Ca.ratio.in.sclerites %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Mg.Ca.ratio.in.sclerites <- Mg.Ca.ratio.in.sclerites %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Mg.Ca.ratio.in.sclerites <- Mg.Ca.ratio.in.sclerites %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Mg_Ca_ratio <- Mg.Ca.ratio.in.sclerites %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Mg_Ca_ratio)

# Count the number of outlier categories
num_outliers1 <- nrow(Mg.Ca.ratio.in.sclerites %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Mg.Ca.ratio.in.sclerites %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Mg.Ca.ratio.in.sclerites %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p25 <- ggplot(Mg.Ca.ratio.in.sclerites, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Mg.Ca.ratio.in.sclerites$trait_name), "- N:", nrow(Mg.Ca.ratio.in.sclerites)),
       x = "Trait Name",
       y = "mol:mol",
       color = "Outlier Status") +
  theme_minimal()

p25

#26) Number.of.pinnules per tentacle

Number.of.pinnules.per.tentacle$value=as.numeric(Number.of.pinnules.per.tentacle$value)

#We check standard units and standardized in case needed
#in this case, month to years
Number.of.pinnules.per.tentacle$standard_unit=as.factor(Number.of.pinnules.per.tentacle$standard_unit)
levels(Number.of.pinnules.per.tentacle$standard_unit)
Number.of.pinnules.per.tentacle$standardized_unit = "units"

# Convert standard units
Number.of.pinnules.per.tentacle <- Number.of.pinnules.per.tentacle %>%
  mutate(standardized_value = case_when(
    standard_unit == "units" ~ value,          # Already in units
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))
# Identify outliers
Number.of.pinnules.per.tentacle <- Number.of.pinnules.per.tentacle %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers
Number.of.pinnules.per.tentacle <- Number.of.pinnules.per.tentacle %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Number.of.pinnules.per.tentacle <- Number.of.pinnules.per.tentacle %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Number.of.pinnules.per.tentacle <- Number.of.pinnules.per.tentacle %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers
outlier_Number_of_pinnules <- Number.of.pinnules.per.tentacle %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Count the number of outlier categories
num_outliers1 <- nrow(Number.of.pinnules.per.tentacle %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Number.of.pinnules.per.tentacle %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Number.of.pinnules.per.tentacle %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p26 <- ggplot(Number.of.pinnules.per.tentacle, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Number.of.pinnules.per.tentacle$trait_name), "- N:", nrow(Number.of.pinnules.per.tentacle)),
       x = "Trait Name",
       y = "units",  # Y label set to "units"
       color = "Outlier Status") +
  theme_minimal()

# Display the plot
p26



#27) Number.of.tentacles.per.polyp

Number.of.tentacles.per.polyp$value=as.numeric(Number.of.tentacles.per.polyp$value)

#We check standard units and standardized in case needed
#in this case, month to years
Number.of.tentacles.per.polyp$standard_unit=as.factor(Number.of.tentacles.per.polyp$standard_unit)
levels(Number.of.tentacles.per.polyp$standard_unit)
Number.of.tentacles.per.polyp$standardized_unit = "units"

# Convert standard units
Number.of.tentacles.per.polyp <- Number.of.tentacles.per.polyp %>%
  mutate(standardized_value = case_when(
    standard_unit == "units" ~ value,          # Already in units
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers 
Number.of.tentacles.per.polyp <- Number.of.tentacles.per.polyp %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Number.of.tentacles.per.polyp <- Number.of.tentacles.per.polyp %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Number.of.tentacles.per.polyp <- Number.of.tentacles.per.polyp %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Number.of.tentacles.per.polyp <- Number.of.tentacles.per.polyp %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Number_of_tentacles <- Number.of.tentacles.per.polyp %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Number_of_tentacles)

#All values are correctly selected as 8, so there is no need for further plots.


#28) Oocytes size at maturity

Oocyte.size.at.maturity$value=as.numeric(Oocyte.size.at.maturity$value)

#We check standard units and standardized in case needed
#in this case, month to years
Oocyte.size.at.maturity$standard_unit=as.factor(Oocyte.size.at.maturity$standard_unit)
levels(Oocyte.size.at.maturity$standard_unit)
Oocyte.size.at.maturity$standardized_unit = "µm"

# Convert standard units
Oocyte.size.at.maturity <- Oocyte.size.at.maturity %>%
  mutate(standardized_value = case_when(
    standard_unit == "mm" ~ value * 1000,      # Convert mm to micras
    standard_unit == "µm" ~ value,              # Already in micras
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers 
Oocyte.size.at.maturity <- Oocyte.size.at.maturity %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Oocyte.size.at.maturity <- Oocyte.size.at.maturity %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Oocyte.size.at.maturity <- Oocyte.size.at.maturity %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Oocyte.size.at.maturity <- Oocyte.size.at.maturity %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Oocyte_size <- Oocyte.size.at.maturity %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Oocyte_size)

# Count the number of outlier categories
num_outliers1 <- nrow(Oocyte.size.at.maturity %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Oocyte.size.at.maturity %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Oocyte.size.at.maturity %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p28 <- ggplot(Oocyte.size.at.maturity, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Oocyte.size.at.maturity$trait_name), "- N:", nrow(Oocyte.size.at.maturity)),
       x = "Trait Name",
       y = "µm",  # Y label changed to "µm"
       color = "Outlier Status") +
  theme_minimal()

p28

#29) Oogenic cycle duration

Oogenic.cycle.duration$value=as.numeric(Oogenic.cycle.duration$value)

#We check standard units and standardized in case needed
#in this case, month to years
Oogenic.cycle.duration$standard_unit=as.factor(Oogenic.cycle.duration$standard_unit)
levels(Oogenic.cycle.duration$standard_unit)
Oogenic.cycle.duration$standardized_unit = "months"

# Convert standard units
Oogenic.cycle.duration <- Oogenic.cycle.duration %>%
  mutate(standardized_value = case_when(
      standard_unit == "months" ~ value,              # Already in months
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers 
Oogenic.cycle.duration <- Oogenic.cycle.duration %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Oogenic.cycle.duration <- Oogenic.cycle.duration %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Oogenic.cycle.duration <- Oogenic.cycle.duration %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Oogenic.cycle.duration <- Oogenic.cycle.duration %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Oogenic_cycle <- Oogenic.cycle.duration %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

#There are not outliers

# Count the number of outlier categories
num_outliers1 <- nrow(Oogenic.cycle.duration %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Oogenic.cycle.duration %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Oogenic.cycle.duration %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p29 <- ggplot(Oogenic.cycle.duration, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Oogenic.cycle.duration$trait_name), "- N:", nrow(Oogenic.cycle.duration)),
       x = "Trait Name",
       y = "months",  # Y label changed to "months"
       color = "Outlier Status") +
  theme_minimal()

p29

#30) Peduncle length

Peduncle.length$value=as.numeric(Peduncle.length$value)

#We check standard units and standardized in case needed
#in this case, month to years
Peduncle.length$standard_unit=as.factor(Peduncle.length$standard_unit)
levels(Peduncle.length$standard_unit)
Peduncle.length$standardized_unit = "mm"

# Convert standard units
Peduncle.length <- Peduncle.length %>%
  mutate(standardized_value = case_when(
    standard_unit == "cm" ~ value * 10,      # Convert cm to mm
    standard_unit == "mm" ~ value,              # Already in mm
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers 
Peduncle.length <- Peduncle.length %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Peduncle.length <- Peduncle.length %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Peduncle.length <- Peduncle.length %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Peduncle.length <- Peduncle.length %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Peduncle_length <- Peduncle.length %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Peduncle_length)

# Count the number of outlier categories
num_outliers1 <- nrow(Peduncle.length %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Peduncle.length %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Peduncle.length %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p30 <- ggplot(Peduncle.length, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Peduncle.length$trait_name), "- N:", nrow(Peduncle.length)),
       x = "Trait Name",
       y = "mm",  # Y label changed to "mm"
       color = "Outlier Status") +
  theme_minimal()

p30

#31) Photosynthesis/respiration ratio

Photosynthesis.respiration.ratio..P.R.$value=as.numeric(Photosynthesis.respiration.ratio..P.R.$value)

#We check standard units and standardized in case needed
#in this case, month to years
Photosynthesis.respiration.ratio..P.R.$standard_unit=as.factor(Photosynthesis.respiration.ratio..P.R.$standard_unit)
levels(Photosynthesis.respiration.ratio..P.R.$standard_unit)
Photosynthesis.respiration.ratio..P.R.$standardized_unit = "x:y"

# Convert standard units
Photosynthesis.respiration.ratio..P.R. <- Photosynthesis.respiration.ratio..P.R. %>%
  mutate(standardized_value = case_when(
    standard_unit == "x:y" ~ value,              # Already in x:y
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers 
Photosynthesis.respiration.ratio..P.R. <- Photosynthesis.respiration.ratio..P.R. %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Photosynthesis.respiration.ratio..P.R. <- Photosynthesis.respiration.ratio..P.R. %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Photosynthesis.respiration.ratio..P.R. <- Photosynthesis.respiration.ratio..P.R. %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Photosynthesis.respiration.ratio..P.R. <- Photosynthesis.respiration.ratio..P.R. %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Photosynthesis_respiration_ratio_P_R <- Photosynthesis.respiration.ratio..P.R. %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

#There are not outliers 

# Count the number of outlier categories
num_outliers1 <- nrow(Photosynthesis.respiration.ratio..P.R. %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Photosynthesis.respiration.ratio..P.R. %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Photosynthesis.respiration.ratio..P.R. %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p31 <- ggplot(Photosynthesis.respiration.ratio..P.R., aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Photosynthesis.respiration.ratio..P.R.$trait_name), "- N:", nrow(Photosynthesis.respiration.ratio..P.R.)),
       x = "Trait Name",
       y = "x:y",  # Y label changed to "x:y"
       color = "Outlier Status") +
  theme_minimal()

p31

#32) Planktonic phase duration in larvae

Planktonic.phase.duration.in.larvae.$value=as.numeric(Planktonic.phase.duration.in.larvae.$value)

#We check standard units and standardized in case needed
#in this case, month to years
Planktonic.phase.duration.in.larvae.$standard_unit=as.factor(Planktonic.phase.duration.in.larvae.$standard_unit)
levels(Planktonic.phase.duration.in.larvae.$standard_unit)
Planktonic.phase.duration.in.larvae.$standardized_unit = "days"

# Convert standard units
Planktonic.phase.duration.in.larvae. <- Planktonic.phase.duration.in.larvae. %>%
  mutate(standardized_value = case_when(
    standard_unit == "months" ~ value * 30,      # Convert months to days (assuming an average month is ~30 days)
    standard_unit == "days" ~ value,              # Already in days
    TRUE ~ NA_real_                                # Assign NA for unsupported units
  ))

# Identify outliers 
Planktonic.phase.duration.in.larvae. <- Planktonic.phase.duration.in.larvae. %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Planktonic.phase.duration.in.larvae. <- Planktonic.phase.duration.in.larvae. %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Planktonic.phase.duration.in.larvae. <- Planktonic.phase.duration.in.larvae. %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Planktonic.phase.duration.in.larvae. <- Planktonic.phase.duration.in.larvae. %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Planktonic_phase_duration <- Planktonic.phase.duration.in.larvae. %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Planktonic_phase_duration)

# Count the number of outlier categories
num_outliers1 <- nrow(Planktonic.phase.duration.in.larvae. %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Planktonic.phase.duration.in.larvae. %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Planktonic.phase.duration.in.larvae. %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p32 <- ggplot(Planktonic.phase.duration.in.larvae., aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Planktonic.phase.duration.in.larvae.$trait_name), "- N:", nrow(Planktonic.phase.duration.in.larvae.)),
       x = "Trait Name",
       y = "days",  # Y label set to "days"
       color = "Outlier Status") +
  theme_minimal()

p32

#33) Polyp density

Polyp.density$value=as.numeric(Polyp.density$value)

#We check standard units and standardized in case needed
#in this case, month to years
Polyp.density$standard_unit=as.factor(Polyp.density$standard_unit)
levels(Polyp.density$standard_unit)
Polyp.density$standardized_unit = "units cm^-1"

# We remove density/area units 
Polyp.density_b <- Polyp.density[Polyp.density$standard_unit!="units cm^-2",]
# Convert standard units
Polyp.density_b <- Polyp.density_b %>%
  mutate(standardized_value = case_when(
    standard_unit == "units 10 cm^-1" ~ value / 10,   # Convert to units cm^-1
    standard_unit == "units cm^-1" ~ value,            # Already in units cm^-1
     standard_unit == "units mm^-1" ~ value * 10,      # Convert mm^-1 to cm^-1
    TRUE ~ NA_real_                                    # Assign NA for unsupported units
  ))

# Identify outliers 
Polyp.density_b <- Polyp.density_b %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Polyp.density_b <- Polyp.density_b %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Polyp.density_b <- Polyp.density_b %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Polyp.density_b <- Polyp.density_b %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Polyp_density_b <- Polyp.density_b %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Polyp_density_b)

# Count the number of outlier categories
num_outliers1 <- nrow(Polyp.density_b %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Polyp.density_b %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Polyp.density_b %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p33 <- ggplot(Polyp.density_b, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Polyp.density_b$trait_name), "- N:", nrow(Polyp.density_b)),
       x = "Trait Name",
       y = "units cm^-1",  # Y label set to "units cm^-1"
       color = "Outlier Status") +
  theme_minimal()

# Display the plot
p33


#34) Polyp diameter

Polyp.diameter$value=as.numeric(Polyp.diameter$value)

#We check standard units and standardized in case needed
#in this case, month to years
Polyp.diameter$standard_unit=as.factor(Polyp.diameter$standard_unit)
levels(Polyp.diameter$standard_unit)
Polyp.diameter$standardized_unit = "mm"

# Convert standard units
Polyp.diameter <- Polyp.diameter %>%
  mutate(standardized_value = case_when(
    standard_unit == "cm" ~ value * 10,            # Convert cm to mm
    standard_unit == "mm" ~ value,                  # Already in mm
    standard_unit == "µm" ~ value / 1000,           # Convert micras (µm) to mm
    TRUE ~ NA_real_                                  # Assign NA for unsupported units
  ))

# Identify outliers 
Polyp.diameter <- Polyp.diameter %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Polyp.diameter <- Polyp.diameter %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both outlier types
Polyp.diameter <- Polyp.diameter %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Polyp.diameter <- Polyp.diameter %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Polyp_diameter <- Polyp.diameter %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Polyp_diameter)

# Count the number of outlier categories
num_outliers1 <- nrow(Polyp.diameter %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Polyp.diameter %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Polyp.diameter %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p34 <- ggplot(Polyp.diameter, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Polyp.diameter$trait_name), "- N:", nrow(Polyp.diameter)),
       x = "Trait Name",
       y = "mm",  # Y label set to "mm"
       color = "Outlier Status") +
  theme_minimal()

p34

#35) Polyp fecundity

Polyp.fecundity$value=as.numeric(Polyp.fecundity$value)

#We check standard units and standardized in case needed
#in this case, month to years
Polyp.fecundity$standard_unit=as.factor(Polyp.fecundity$standard_unit)
levels(Polyp.fecundity$standard_unit)
Polyp.fecundity$standardized_unit = "units"

##We will focus in the most numerous "units" group to search for outliers
Polyp.fecundity_b <- Polyp.fecundity[Polyp.fecundity$standard_unit=="units",]

# Transform units (Assuming units are already in the desired form)
Polyp.fecundity_b <- Polyp.fecundity_b %>%
  mutate(standardized_value = case_when(
    standard_unit == "units" ~ value,              # Already in units
    TRUE ~ NA_real_                                 # Assign NA for unsupported units
  ))

# Identify outliers 
Polyp.fecundity_b <- Polyp.fecundity_b %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Polyp.fecundity_b <- Polyp.fecundity_b %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both outlier types
Polyp.fecundity_b <- Polyp.fecundity_b %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Polyp.fecundity_b <- Polyp.fecundity_b %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Polyp_fecundity_b <- Polyp.fecundity_b %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Polyp_fecundity_b)

# Count the number of outlier categories
num_outliers1 <- nrow(Polyp.fecundity_b %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Polyp.fecundity_b %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Polyp.fecundity_b %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p35 <- ggplot(Polyp.fecundity_b, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Polyp.fecundity_b$trait_name), "- N:", nrow(Polyp.fecundity_b)),
       x = "Trait Name",
       y = "units",  # Y label set to "units"
       color = "Outlier Status") +
  theme_minimal()

p35

#36) Polyp height

Polyp.height$value=as.numeric(Polyp.height$value)

#We check standard units and standardized in case needed
#in this case, month to years
Polyp.height$standard_unit=as.factor(Polyp.height$standard_unit)
levels(Polyp.height$standard_unit)
Polyp.height$standardized_unit = "mm"

# Convert standard units
Polyp.height <- Polyp.height %>%
  mutate(standardized_value = case_when(
    standard_unit == "cm" ~ value * 10,            # Convert cm to mm
    standard_unit == "mm" ~ value,                  # Already in mm
     TRUE ~ NA_real_                                  # Assign NA for unsupported units
  ))

# Identify outliers 
Polyp.height <- Polyp.height %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Polyp.height <- Polyp.height %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both outlier types
Polyp.height <- Polyp.height %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Polyp.height <- Polyp.height %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Polyp_height <- Polyp.height %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Polyp_height)

# Count the number of outlier categories
num_outliers1 <- nrow(Polyp.height %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Polyp.height %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Polyp.height %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p36 <- ggplot(Polyp.height, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Polyp.height$trait_name), "- N:", nrow(Polyp.height)),
       x = "Trait Name",
       y = "mm",  # Y label set to "mm"
       color = "Outlier Status") +
  theme_minimal()

p36

#37) Polyp whorl diameter

Polyp.whorl.diameter$value=as.numeric(Polyp.whorl.diameter$value)

#We check standard units and standardized in case needed
#in this case, month to years
Polyp.whorl.diameter$standard_unit=as.factor(Polyp.whorl.diameter$standard_unit)
levels(Polyp.whorl.diameter$standard_unit)
Polyp.whorl.diameter$standardized_unit = "mm"

# Convert standard units
Polyp.whorl.diameter <- Polyp.whorl.diameter %>%
  mutate(standardized_value = case_when(
     standard_unit == "mm" ~ value,                  # Already in mm
    TRUE ~ NA_real_                                  # Assign NA for unsupported units
  ))

# Identify outliers 
Polyp.whorl.diameter <- Polyp.whorl.diameter %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Polyp.whorl.diameter <- Polyp.whorl.diameter %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both outlier types
Polyp.whorl.diameter <- Polyp.whorl.diameter %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Polyp.whorl.diameter <- Polyp.whorl.diameter %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Polyp_whorl_diameter <- Polyp.whorl.diameter %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Polyp_whorl_diameter)

# Count the number of outlier categories
num_outliers1 <- nrow(Polyp.whorl.diameter %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Polyp.whorl.diameter %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Polyp.whorl.diameter %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p37 <- ggplot(Polyp.whorl.diameter, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Polyp.whorl.diameter$trait_name), "- N:", nrow(Polyp.whorl.diameter)),
       x = "Trait Name",
       y = "mm",  # Y label set to "mm"
       color = "Outlier Status") +
  theme_minimal()

p37

#38) Polyp whorl/pairs per length unit

Polyp.whorls.pairs.per.lenght.unit$value=as.numeric(Polyp.whorls.pairs.per.lenght.unit$value)

#We check standard units and standardized in case needed
#in this case, month to years
Polyp.whorls.pairs.per.lenght.unit$standard_unit=as.factor(Polyp.whorls.pairs.per.lenght.unit$standard_unit)
levels(Polyp.whorls.pairs.per.lenght.unit$standard_unit)
Polyp.whorls.pairs.per.lenght.unit$standardized_unit = "units cm^-1"

# Convert standard units
# Transforming units to "units cm^-1"
Polyp.whorls.pairs.per.lenght.unit <- Polyp.whorls.pairs.per.lenght.unit %>%
  mutate(standardized_value = case_when(
    standard_unit == "units 2cm^-1" ~ value / 2,       # Convert to units cm^-1
    standard_unit == "units 3cm^-1" ~ value / 3,       # Convert to units cm^-1
    standard_unit == "units cm^-1" ~ value,             # Already in units cm^-1
    TRUE ~ NA_real_                                     # Assign NA for unsupported units
  ))

# Identify outliers
Polyp.whorls.pairs.per.lenght.unit <- Polyp.whorls.pairs.per.lenght.unit %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers
Polyp.whorls.pairs.per.lenght.unit <- Polyp.whorls.pairs.per.lenght.unit %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Polyp.whorls.pairs.per.lenght.unit <- Polyp.whorls.pairs.per.lenght.unit %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Polyp.whorls.pairs.per.lenght.unit <- Polyp.whorls.pairs.per.lenght.unit %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers
outlier_Polyp_whorls_pairs <- Polyp.whorls.pairs.per.lenght.unit %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Polyp_whorls_pairs)

# Count the number of outlier categories
num_outliers1 <- nrow(Polyp.whorls.pairs.per.lenght.unit %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Polyp.whorls.pairs.per.lenght.unit %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Polyp.whorls.pairs.per.lenght.unit %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p38 <- ggplot(Polyp.whorls.pairs.per.lenght.unit, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Polyp.whorls.pairs.per.lenght.unit$trait_name), "- N:", nrow(Polyp.whorls.pairs.per.lenght.unit)),
       x = "Trait Name",
       y = "units cm^-1",  # Y label set to "units cm^-1"
       color = "Outlier Status") +
  theme_minimal()

p38

#39) Polyps per whorl

Polyps.per.whorl$value=as.numeric(Polyps.per.whorl$value)

#We check standard units and standardized in case needed
#in this case, month to years
Polyps.per.whorl$standard_unit=as.factor(Polyps.per.whorl$standard_unit)
levels(Polyps.per.whorl$standard_unit)
Polyps.per.whorl$standardized_unit = "units"

# Convert standard units
# Transforming units to "units"
Polyps.per.whorl <- Polyps.per.whorl %>%
  mutate(standardized_value = case_when(
     standard_unit == "units" ~ value,             # Already in units cm^-1
    TRUE ~ NA_real_                                     # Assign NA for unsupported units
  ))

# Identify outliers
Polyps.per.whorl <- Polyps.per.whorl %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers
Polyps.per.whorl <- Polyps.per.whorl %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Polyps.per.whorl <- Polyps.per.whorl %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Polyps.per.whorl <- Polyps.per.whorl %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers
outlier_Polyps_per_whorl <- Polyps.per.whorl %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Polyps_per_whorl)

# Count the number of outlier categories
num_outliers1 <- nrow(Polyps.per.whorl %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Polyps.per.whorl %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Polyps.per.whorl %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p39 <- ggplot(Polyps.per.whorl, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Polyps.per.whorl$trait_name), "- N:", nrow(Polyps.per.whorl)),
       x = "Trait Name",
       y = "units",  # Y label set to "units"
       color = "Outlier Status") +
  theme_minimal()

p39

#40) Projected side area

Projected.side.area$value=as.numeric(Projected.side.area$value)

#We check standard units and standardized in case needed
#in this case, month to years
Projected.side.area$standard_unit=as.factor(Projected.side.area$standard_unit)
levels(Projected.side.area$standard_unit)
Projected.side.area$standardized_unit = "cm^2"

# Convert standard units
# Transforming units to "units cm^2"
Projected.side.area <- Projected.side.area %>%
  mutate(standardized_value = case_when(
    standard_unit == "cm^2" ~ value,   
    standard_unit == "mm^2" ~ value / 100,  # Convert mm² to cm²# Already in units cm^2
    TRUE ~ NA_real_                                     # Assign NA for unsupported units
  ))

# Identify outliers 
Projected.side.area <- Projected.side.area %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Projected.side.area <- Projected.side.area %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Projected.side.area <- Projected.side.area %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Projected.side.area <- Projected.side.area %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Project.side.area <- Projected.side.area %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Project.side.area)

# Count the number of outlier categories
num_outliers1 <- nrow(Projected.side.area %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Projected.side.area %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Projected.side.area %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p40 <- ggplot(Projected.side.area, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Projected.side.area$trait_name), "- N:", nrow(Projected.side.area)),
       x = "Trait Name",
       y = "cm²",  # Y label set to "cm²"
       color = "Outlier Status") +
  theme_minimal()

# Print the plot
p40


#41) Propagule size on realease

Propagule.size.on.release$value=as.numeric(Propagule.size.on.release$value)

#We check standard units and standardized in case needed
#in this case, month to years
Propagule.size.on.release$standard_unit=as.factor(Propagule.size.on.release$standard_unit)
levels(Propagule.size.on.release$standard_unit)
Propagule.size.on.release$standardized_unit = "µm"

# Convert standard units
# Transforming units to standard
Propagule.size.on.release <- Propagule.size.on.release %>%
  mutate(standardized_value = case_when(
    standard_unit == "mm" ~ value * 1000,      # Convert mm to µm
    standard_unit == "µm" ~ value,          # Already in µm
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers 
Propagule.size.on.release <- Propagule.size.on.release %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Propagule.size.on.release <- Propagule.size.on.release %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Propagule.size.on.release <- Propagule.size.on.release %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Propagule.size.on.release <- Propagule.size.on.release %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Propagule_size_on_release <- Propagule.size.on.release %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Propagule_size_on_release)

# Count the number of outlier categories
num_outliers1 <- nrow(Propagule.size.on.release %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Propagule.size.on.release %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Propagule.size.on.release %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p41 <- ggplot(Propagule.size.on.release, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Propagule.size.on.release$trait_name), "- N:", nrow(Propagule.size.on.release)),
       x = "Trait Name",
       y = "µm",  # Y label set to "µm"
       color = "Outlier Status") +
  theme_minimal()

p41

#42) Rachis length

Rachis.length$value=as.numeric(Rachis.length$value)

#We check standard units and standardized in case needed
#in this case, month to years
Rachis.length$standard_unit=as.factor(Rachis.length$standard_unit)
levels(Rachis.length$standard_unit)
Rachis.length$standardized_unit = "cm"

# Convert standard units
# Transforming units to standard
Rachis.length <- Rachis.length %>%
  mutate(standardized_value = case_when(
    standard_unit == "mm" ~ value /10,      # Convert mm to cm
    standard_unit == "cm" ~ value,          # Already in cm
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers 
Rachis.length <- Rachis.length %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Rachis.length <- Rachis.length %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Rachis.length <- Rachis.length %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Rachis.length <- Rachis.length %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Rachis_length <- Rachis.length %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Rachis_length)

# Count the number of outlier categories
num_outliers1 <- nrow(Rachis.length %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Rachis.length %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Rachis.length %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p42 <- ggplot(Rachis.length, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Rachis.length$trait_name), "- N:", nrow(Rachis.length)),
       x = "Trait Name",
       y = "cm",  # Y label set to "cm"
       color = "Outlier Status") +
  theme_minimal()

# Print the plot
p42


#43) Size at maturity

Size.at.maturity$value=as.numeric(Size.at.maturity$value)

#We check standard units and standardized in case needed
#in this case, month to years
Size.at.maturity$standard_unit=as.factor(Size.at.maturity$standard_unit)
levels(Size.at.maturity$standard_unit)
Size.at.maturity$standardized_unit = "cm"

#We focuse on size on length only
Size.at.maturity_b<-Size.at.maturity[Size.at.maturity$standard_unit!="cm^3",]

# Convert standard units
# Transforming units to standard
Size.at.maturity_b <- Size.at.maturity_b %>%
  mutate(standardized_value = case_when(
    standard_unit == "mm" ~ value /10,      # Convert mm to cm
    standard_unit == "cm" ~ value,          # Already in cm
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))


# Identify outliers 
Size.at.maturity_b <- Size.at.maturity_b %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers 
Size.at.maturity_b <- Size.at.maturity_b %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Size.at.maturity_b <- Size.at.maturity_b %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Size.at.maturity_b <- Size.at.maturity_b %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Size_at_maturity_b <- Size.at.maturity_b %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Size_at_maturity_b)

# Count the number of outlier categories
num_outliers1 <- nrow(Size.at.maturity_b %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Size.at.maturity_b %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Size.at.maturity_b %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p43 <- ggplot(Size.at.maturity_b, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Size.at.maturity_b$trait_name), "- N:", nrow(Size.at.maturity_b)),
       x = "Trait Name",
       y = "cm",  # Y label set to "cm"
       color = "Outlier Status") +
  theme_minimal()

# Print the plot
p43


#44) Size of prey items

Size.of.prey.items$value=as.numeric(Size.of.prey.items$value)

#We check standard units and standardized in case needed
#in this case, month to years
Size.of.prey.items$standard_unit=as.factor(Size.of.prey.items$standard_unit)
levels(Size.of.prey.items$standard_unit)
Size.of.prey.items$standardized_unit = "µm"

# Convert standard units
# Transforming units to standard
Size.of.prey.items <- Size.of.prey.items %>%
  mutate(standardized_value = case_when(
    standard_unit == "mm" ~ value *1000,      # Convert cm to mm
    standard_unit == "µm" ~ value,          # Already in µm
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers 
Size.of.prey.items <- Size.of.prey.items %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))  # Use original value for outlier detection

# Identify extreme outliers 
Size.of.prey.items <- Size.of.prey.items %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))  # Use original value for extreme outlier detection

# Merge both
Size.of.prey.items <- Size.of.prey.items %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Size.of.prey.items <- Size.of.prey.items %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Size_of_prey_items <- Size.of.prey.items %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Print the outliers table
print(outlier_Size_of_prey_items)

# Count the number of outlier categories
num_outliers1 <- nrow(Size.of.prey.items %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Size.of.prey.items %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Size.of.prey.items %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p44 <- ggplot(Size.of.prey.items, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Size.of.prey.items$trait_name), "- N:", nrow(Size.of.prey.items)),
       x = "Trait Name",
       y = "µm",  # Y label set to "mm"
       color = "Outlier Status") +
  theme_minimal()

# Print the plot
p44

#45) Skeletal density. Only one observation was compiled for that trait

#46) Spermatogenic cycle duration

Spermatogenic.cycle.duration$value=as.numeric(Spermatogenic.cycle.duration$value)

#We check standard units and standardized in case needed
#in this case, month to years
Spermatogenic.cycle.duration$standard_unit=as.factor(Spermatogenic.cycle.duration$standard_unit)
levels(Spermatogenic.cycle.duration$standard_unit)
Spermatogenic.cycle.duration$standardized_unit = "months"

# Convert standard units
# Transforming units to standard
Spermatogenic.cycle.duration <- Spermatogenic.cycle.duration %>%
  mutate(standardized_value = case_when(
    standard_unit == "months" ~ value,          # Already in mm
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers 
Spermatogenic.cycle.duration <- Spermatogenic.cycle.duration %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))  # Use original value for outlier detection

# Identify extreme outliers 
Spermatogenic.cycle.duration <- Spermatogenic.cycle.duration %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))  # Use original value for extreme outlier detection

# Merge both
Spermatogenic.cycle.duration <- Spermatogenic.cycle.duration %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Spermatogenic.cycle.duration <- Spermatogenic.cycle.duration %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Spermatogenic_cycle_duration <- Spermatogenic.cycle.duration %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Count the number of outlier categories
num_outliers1 <- nrow(Spermatogenic.cycle.duration %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Spermatogenic.cycle.duration %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Spermatogenic.cycle.duration %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p46 <- ggplot(Spermatogenic.cycle.duration, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Spermatogenic.cycle.duration$trait_name), "- N:", nrow(Spermatogenic.cycle.duration)),
       x = "Trait Name",
       y = "months",  # Y label set to "months"
       color = "Outlier Status") +
  theme_minimal()

# Display the plot object
p46

#47) Stiffness Young Modulus 

Stiffness_Young.s_modulus$value=as.numeric(Stiffness_Young.s_modulus$value)

#We check standard units and standardized in case needed
#in this case, month to years
Stiffness_Young.s_modulus$standard_unit=as.factor(Stiffness_Young.s_modulus$standard_unit)
levels(Stiffness_Young.s_modulus$standard_unit)


Stiffness_Young.s_modulus$standardized_unit = "G dynes cm^-2"

##We focus on the standard with more observations
Stiffness_Young.s_modulus_b <- Stiffness_Young.s_modulus[Stiffness_Young.s_modulus$standard_id=="st89",]

# Convert standard units
# Transforming units to standard
Stiffness_Young.s_modulus_b <- Stiffness_Young.s_modulus_b %>%
  mutate(standardized_value = case_when(
    standard_unit == "G dynes cm^-2" ~ value,          # Already in G dynes cm^-2
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers 
Stiffness_Young.s_modulus_b <- Stiffness_Young.s_modulus_b %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))  # Use original value for outlier detection

# Identify extreme outliers 
Stiffness_Young.s_modulus_b <- Stiffness_Young.s_modulus_b %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))  # Use original value for extreme outlier detection

# Merge both
Stiffness_Young.s_modulus_b <- Stiffness_Young.s_modulus_b %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Stiffness_Young.s_modulus_b <- Stiffness_Young.s_modulus_b %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Stiffness_Young_s_modulus_b <- Stiffness_Young.s_modulus_b %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Count the number of outlier categories
num_outliers1 <- nrow(Stiffness_Young.s_modulus_b %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Stiffness_Young.s_modulus_b %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Stiffness_Young.s_modulus_b %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p47 <- ggplot(Stiffness_Young.s_modulus_b, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Stiffness_Young.s_modulus_b$trait_name), "- N:", nrow(Stiffness_Young.s_modulus_b)),
       x = "Trait Name",
       y = "G dynes cm^-2",  # Y label set to "G dynes cm^-2"
       color = "Outlier Status") +
  theme_minimal()

# Display the plot object
p47

#48) Symbiont density 

Symbiont.density$value=as.numeric(Symbiont.density$value)

#We check standard units and standardized in case needed
#in this case, month to years
Symbiont.density$standard_unit=as.factor(Symbiont.density$standard_unit)
levels(Symbiont.density$standard_unit)

#There are few observation and with multiple valid standards so, we manually check that all values
#were correct in scientific sources 

#49) Tentacle length 

Tentacle.length$value=as.numeric(Tentacle.length$value)

#We check standard units and standardized in case needed
#in this case, month to years
Tentacle.length$standard_unit=as.factor(Tentacle.length$standard_unit)
levels(Tentacle.length$standard_unit)
Tentacle.length$standardized_unit = "mm"

# Convert standard units
# Transforming units to standard
Tentacle.length <- Tentacle.length %>%
  mutate(standardized_value = case_when(
    standard_unit == "cm" ~ value *10,      # Convert cm to mm
    standard_unit == "mm" ~ value,          # Already in mm
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers 
Tentacle.length <- Tentacle.length %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))  # Use original value for outlier detection

# Identify extreme outliers 
Tentacle.length <- Tentacle.length %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))  # Use original value for extreme outlier detection

# Merge both
Tentacle.length <- Tentacle.length %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Tentacle.length <- Tentacle.length %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Tentacle_length <- Tentacle.length %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Count the number of outlier categories
num_outliers1 <- nrow(Tentacle.length %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Tentacle.length %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Tentacle.length %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p49 <- ggplot(Tentacle.length, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Tentacle.length$trait_name), "- N:", nrow(Tentacle.length)),
       x = "Trait Name",
       y = "mm",  # Y label set to "mm"
       color = "Outlier Status") +
  theme_minimal()

# Display the plot object
p49

#50) Tissue thickness. There are few (4) observations for this trait, so we checked that the introduced values are correct in the primary sources

#51) Total biomass. There are few (11) observations for this trait and with multiple valid standards so we checked that the introduced values are correct in the primary sources

#52) Trunk diameter

Trunk.diameter$value=as.numeric(Trunk.diameter$value)

#We check standard units and standardized in case needed
#in this case, month to years
Trunk.diameter$standard_unit=as.factor(Trunk.diameter$standard_unit)
levels(Trunk.diameter$standard_unit)
Trunk.diameter$standardized_unit = "mm"

# Convert standard units
# Transforming units to standard
Trunk.diameter <- Trunk.diameter %>%
  mutate(standardized_value = case_when(
    standard_unit == "cm" ~ value *10,      # Convert cm to mm
    standard_unit == "mm" ~ value,          # Already in mm
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers 
Trunk.diameter <- Trunk.diameter %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))  # Use original value for outlier detection

# Identify extreme outliers 
Trunk.diameter <- Trunk.diameter %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))  # Use original value for extreme outlier detection

# Merge both
Trunk.diameter <- Trunk.diameter %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Trunk.diameter <- Trunk.diameter %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers 
outlier_Trunk_diameter <- Trunk.diameter %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Count the number of outlier categories
num_outliers1 <- nrow(Trunk.diameter %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Trunk.diameter %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Trunk.diameter %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p52 <- ggplot(Trunk.diameter, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Trunk.diameter$trait_name), "- N:", nrow(Trunk.diameter)),
       x = "Trait Name",
       y = "mm",  # Y label set to "mm"
       color = "Outlier Status") +
  theme_minimal()

# Display the plot object
p52

#53) Type locality depth

Type.locality.depth$value=as.numeric(Type.locality.depth$value)

#We check standard units and standardized in case needed
#in this case, month to years
Type.locality.depth$standard_unit=as.factor(Type.locality.depth$standard_unit)
levels(Type.locality.depth$standard_unit)
Type.locality.depth$standardized_unit = "m"

# Convert standard units
# Transforming units to standard
Type.locality.depth <- Type.locality.depth %>% 
  mutate(standardized_value = case_when(
    standard_unit == "cm" ~ value / 100,        # Convert cm to m
    standard_unit == "ft" ~ value * 0.3048,     # Convert ft to m
    standard_unit == "fm" ~ value * 0.3048,     # Convert fm to m
    standard_unit == "m" ~ value,                # Already in m
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers
Type.locality.depth <- Type.locality.depth %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers
Type.locality.depth <- Type.locality.depth %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
Type.locality.depth <- Type.locality.depth %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
Type.locality.depth <- Type.locality.depth %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers
outlier_Type_locality_depth <- Type.locality.depth %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Count the number of outlier categories
num_outliers1 <- nrow(Type.locality.depth %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(Type.locality.depth %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(Type.locality.depth %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p53 <- ggplot(Type.locality.depth, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(Type.locality.depth$trait_name), "- N:", nrow(Type.locality.depth)),
       x = "Trait Name",
       y = "m",  # Y label set to "m"
       color = "Outlier Status") +
  theme_minimal()

# Display the plot
p53

#54) CaCO3 in sclerites has only 5 observations, so we manually checked that the compiled values were correct using the original sources

#55) CaCO3 in the axis

X.CaCO3..in.the.axis$value=as.numeric(X.CaCO3..in.the.axis$value)

#We check standard units and standardized in case needed
#in this case, month to years
X.CaCO3..in.the.axis$standard_unit=as.factor(X.CaCO3..in.the.axis$standard_unit)
levels(X.CaCO3..in.the.axis$standard_unit)
X.CaCO3..in.the.axis$standardized_unit = "mmol Ca2+ g^-1 dry wt"
#Most observations are from the "mmol Ca2+ g^-1 dry wt" standard so we will focus on that one
X.CaCO3..in.the.axis_b <-X.CaCO3..in.the.axis[X.CaCO3..in.the.axis$standard_unit=="mmol Ca2+ g^-1 dry wt",]

# Convert standard units
# Transforming units to standard
X.CaCO3..in.the.axis_b <- X.CaCO3..in.the.axis_b %>% 
  mutate(standardized_value = case_when(
     standard_unit == "mmol Ca2+ g^-1 dry wt" ~ value,                # Already in mmol Ca2+ g^-1 dry wt
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers
X.CaCO3..in.the.axis_b <- X.CaCO3..in.the.axis_b %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers
X.CaCO3..in.the.axis_b <- X.CaCO3..in.the.axis_b %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
X.CaCO3..in.the.axis_b <- X.CaCO3..in.the.axis_b %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
X.CaCO3..in.the.axis_b <- X.CaCO3..in.the.axis_b %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers
outlier_X_CaCO3 <- X.CaCO3..in.the.axis_b %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Count the number of outlier categories
num_outliers1 <- nrow(X.CaCO3..in.the.axis_b %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(X.CaCO3..in.the.axis_b %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(X.CaCO3..in.the.axis_b %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p55 <- ggplot(X.CaCO3..in.the.axis_b, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(X.CaCO3..in.the.axis_b$trait_name), "- N:", nrow(X.CaCO3..in.the.axis_b)),
       x = "Trait Name",
       y = "mmol Ca2+ g^-1 dry wt",  # Y label set to "mmol Ca2+ g^-1 dry wt"
       color = "Outlier Status") +
  theme_minimal()

# Display the plot
p55

#56) Mg in the axis

X.Mg..in.axis$value=as.numeric(X.Mg..in.axis$value)

#We check standard units and standardized in case needed
#in this case, month to years
X.Mg..in.axis$standard_unit=as.factor(X.Mg..in.axis$standard_unit)
levels(X.Mg..in.axis$standard_unit)
X.Mg..in.axis$standardized_unit = "mol %"
#Most observations are from the "mol %" standard so we will focus on that one
X.Mg..in.axis_b <-X.Mg..in.axis[X.Mg..in.axis$standard_unit=="mol %",]

# Convert standard units
# Transforming units to standard
X.Mg..in.axis_b <- X.Mg..in.axis_b %>% 
  mutate(standardized_value = case_when(
    standard_unit == "mol %" ~ value,                # Already in mol %
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers
X.Mg..in.axis_b <- X.Mg..in.axis_b %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers
X.Mg..in.axis_b <- X.Mg..in.axis_b %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
X.Mg..in.axis_b <- X.Mg..in.axis_b %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
X.Mg..in.axis_b <- X.Mg..in.axis_b %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers
outlier_X_Mg <- X.Mg..in.axis_b %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Count the number of outlier categories
num_outliers1 <- nrow(X.Mg..in.axis_b %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(X.Mg..in.axis_b %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(X.Mg..in.axis_b %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p56 <- ggplot(X.Mg..in.axis_b, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(X.Mg..in.axis_b$trait_name), "- N:", nrow(X.Mg..in.axis_b)),
       x = "Trait Name",
       y = "mol%",  # Y label set to "mol%"
       color = "Outlier Status") +
  theme_minimal()

# Display the plot
p56

#57) MgCO3 in axis

X.MgCO3..in.axis$value=as.numeric(X.MgCO3..in.axis$value)

#We check standard units and standardized in case needed
#in this case, month to years
X.MgCO3..in.axis$standard_unit=as.factor(X.MgCO3..in.axis$standard_unit)
levels(X.MgCO3..in.axis$standard_unit)
X.MgCO3..in.axis$standardized_unit = "mol % in MgCa(CO3)2"

# Convert standard units
# Transforming units to standard
X.MgCO3..in.axis <- X.MgCO3..in.axis %>% 
  mutate(standardized_value = case_when(
    standard_unit == "mol % in MgCa(CO3)2" ~ value,                # Already in mol % in MgCa(CO3)2
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers
X.MgCO3..in.axis <- X.MgCO3..in.axis %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers
X.MgCO3..in.axis <- X.MgCO3..in.axis %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
X.MgCO3..in.axis <- X.MgCO3..in.axis %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
X.MgCO3..in.axis <- X.MgCO3..in.axis %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers
outlier_X_MgCO3 <- X.MgCO3..in.axis %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Count the number of outlier categories
num_outliers1 <- nrow(X.MgCO3..in.axis %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(X.MgCO3..in.axis %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(X.MgCO3..in.axis %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p57 <- ggplot(X.MgCO3..in.axis, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(X.MgCO3..in.axis$trait_name), "- N:", nrow(X.MgCO3..in.axis)),
       x = "Trait Name",
       y = "mol % in MgCa(CO3)2",  # Y label set to "mol % in MgCa(CO3)2"
       color = "Outlier Status") +
  theme_minimal()

# Display the plot
p57

#58) MgCO3 in sclerites

X.MgCO3..in.sclerites$value=as.numeric(X.MgCO3..in.sclerites$value)

#We check standard units and standardized in case needed
#in this case, month to years
X.MgCO3..in.sclerites$standard_unit=as.factor(X.MgCO3..in.sclerites$standard_unit)
levels(X.MgCO3..in.sclerites$standard_unit)
#We focus on the standard with more observations
X.MgCO3..in.sclerites$standardized_unit = "mol %"

X.MgCO3..in.sclerites_b <-X.MgCO3..in.sclerites[X.MgCO3..in.sclerites$standard_unit=="mol %",]
# Convert standard units
# Transforming units to standard
X.MgCO3..in.sclerites_b <- X.MgCO3..in.sclerites_b %>% 
  mutate(standardized_value = case_when(
    standard_unit == "mol %" ~ value,                # Already in mol %
    TRUE ~ NA_real_                              # Assign NA for unsupported units
  ))

# Identify outliers
X.MgCO3..in.sclerites_b <- X.MgCO3..in.sclerites_b %>%
  mutate(outlier = identify_outliers(as.numeric(standardized_value)))

# Identify extreme outliers
X.MgCO3..in.sclerites_b <- X.MgCO3..in.sclerites_b %>%
  mutate(extreme_outlier = identify_extreme_outliers(as.numeric(standardized_value)))

# Merge both
X.MgCO3..in.sclerites_b <- X.MgCO3..in.sclerites_b %>%
  mutate(outliers = case_when(
    outlier == FALSE & extreme_outlier == FALSE ~ "FALSE",   # Neither outlier nor extreme outlier
    outlier == TRUE & extreme_outlier == FALSE ~ "OUTLIER",  # Outlier but not extreme
    extreme_outlier == TRUE ~ "EXTREME"                       # Extreme outlier
  ))

# Ensure the outliers column has correct categories: "Not Outlier", "Outlier", and "Extreme"
X.MgCO3..in.sclerites_b <- X.MgCO3..in.sclerites_b %>%
  mutate(outliers = case_when(
    outliers == "FALSE" ~ "Not Outlier",
    outliers == "OUTLIER" ~ "Outlier",
    outliers == "EXTREME" ~ "Extreme"
  ))

# Create a table of outliers
outlier_X_MgCO3_sclerites <- X.MgCO3..in.sclerites_b %>%
  filter(outliers %in% c("Outlier", "Extreme")) %>%  # Filter for outliers and extreme outliers
  select(resource_id, standardized_value, standardized_unit, trait_name, outliers)  # Select relevant columns

# Count the number of outlier categories
num_outliers1 <- nrow(X.MgCO3..in.sclerites_b %>% filter(outliers == "Outlier"))
num_extremes1 <- nrow(X.MgCO3..in.sclerites_b %>% filter(outliers == "Extreme"))
num_normal1 <- nrow(X.MgCO3..in.sclerites_b %>% filter(outliers == "Not Outlier"))

# Create the scatter plot with outliers
p58 <- ggplot(X.MgCO3..in.sclerites_b, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
  geom_jitter(width = 0.1, size = 2, alpha = 0.6) +  # Jitter points
  scale_color_manual(values = c("black", "orange", "red"),  # Correct color mapping
                     breaks = c("Not Outlier", "Outlier", "Extreme"),  # Ensure breaks match unique values
                     labels = c(paste("Not Outlier (N =", num_normal1, ")"),
                                paste("Outlier (N =", num_outliers1, ")"),
                                paste("Extreme (N =", num_extremes1, ")"))) +
  labs(title = paste("Plot for", unique(X.MgCO3..in.sclerites_b$trait_name), "- N:", nrow(X.MgCO3..in.sclerites_b)),
       x = "Trait Name",
       y = "mol %",  # Y label set to "mol %"
       color = "Outlier Status") +
  theme_minimal()

# Display the plot
p58

#59) Larval Swimming Speed. Very few observations and in different standards, so we simply checked that the compiled value were correct in original publications
#60) Dark respiration: Given the few number of values and different standards used, we directly checked that the compiled values were correct in the original publication


# Now we can generate a summary file

# Combine all trait datasets into a single data frame
all_traits <- list(
  Age.at.maturity,
  Axis.radial.growth.rate,
  Branch.diameter,
  Ca.Mg.ratio.in.sclerites,
  Calcification.rate_a,
  Calyx.height,
  Calyx.width,
  Colony.age,
  Colony.area,
  Colony.fecundity,
  Colony.height,
  Colony.width,
  Depth.lower,
  Depth.upper,
  Diameter.of.axis.skeleton,
  Distance.between.polyps,
  Fractal.dimension,
  Growth.rate_b,
  Larval.length,
  Larval.width,
  Longevity,
  Mg.Ca.ratio.in.sclerites,
  Number.of.pinnules.per.tentacle,
  Number.of.tentacles.per.polyp,
  Oocyte.size.at.maturity,
  Oogenic.cycle.duration,
  Peduncle.length,
  Photosynthesis.respiration.ratio..P.R.,
  Planktonic.phase.duration.in.larvae.,
  Polyp.density_b,
  Polyp.diameter,
  Polyp.fecundity_b,
  Polyp.height,
  Polyp.whorl.diameter,
  Polyp.whorls.pairs.per.lenght.unit,
  Polyps.per.whorl,
  Projected.side.area,
  Propagule.size.on.release,
  Rachis.length,
  Size.at.maturity_b,
  Size.of.prey.items,
  Spermatogenic.cycle.duration,
  Stiffness_Young.s_modulus_b,
  Tentacle.length,
  Trunk.diameter,
  Type.locality.depth
) %>%
  bind_rows()

# Create a PDF for plots
pdf("Trait_plots.pdf", width = 6, height = 5)  # PDF size for plots

# Initialize an empty data frame for all outlier data to save to CSV later
all_outlier_data <- data.frame()

# Loop through each unique trait_name for plots
for (trait in unique(all_traits$trait_name)) {
  
  # Filter the data for the current trait
  trait_data <- all_traits %>% filter(trait_name == trait)
  
  # Identify the outliers for the current trait
  outlier_trait_data <- trait_data %>%
    filter(outliers %in% c("Outlier", "Extreme")) %>%
    select(resource_id, value, standard_unit,standardized_value, standardized_unit, trait_name, outliers)
  
  # Append outlier data to the all_outlier_data data frame
  all_outlier_data <- bind_rows(all_outlier_data, outlier_trait_data)
  
  # Count the number of outliers, extremes, and normal points for this trait
  num_outliers <- nrow(trait_data %>% filter(outliers == "Outlier"))
  num_extremes <- nrow(trait_data %>% filter(outliers == "Extreme"))
  num_normal <- nrow(trait_data %>% filter(outliers == "Not Outlier"))
  
  # Create the scatter plot for the current trait
  p <- ggplot(trait_data, aes(x = trait_name, y = as.numeric(standardized_value), color = outliers)) +
    geom_jitter(width = 0.1, size = 1.5, alpha = 0.6) +  # Adjusted point size
    scale_color_manual(values = c("black", "orange", "red"),
                       breaks = c("Not Outlier", "Outlier", "Extreme"),
                       labels = c(paste("Not Outlier (N =", num_normal, ")"),
                                  paste("Outlier (N =", num_outliers, ")"),
                                  paste("Extreme (N =", num_extremes, ")"))) +
    labs(title = paste("Plot for", unique(trait_data$trait_name), "- N:", nrow(trait_data)),
         x = "Trait Name", y = unique(trait_data$standardized_unit), color = "Outlier Status") +
    theme_minimal(base_size = 10) +  # Adjusted base font size
    theme(axis.text.x = element_text(size = 8),  # X-axis text size
          axis.text.y = element_text(size = 8),  # Y-axis text size
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20))  # Set plot margins
  
  # Print the plot to the PDF
  print(p)
  
  # No need to add a new page, ggplot will handle that automatically for each plot
}

dev.off()  # Close the PDF device

## Add metadata for resources

resource_id<-read.csv("/Users/dani/Desktop/new versions OTD/V2_final_03_10_2024/DATASET ZENODO updated/Table_ids/resource_id.csv")


combined_data <- all_outlier_data %>%
  left_join(resource_id, by = "resource_id")

final_data <- combined_data %>%
  select(value,standard_unit,
    standardized_value,
    standardized_unit,
    trait_name,
    outliers,
    resource_id,
    primary_secondary,  # Include this if it exists in all_outlier_data
    author,
    year,
    title,
    resource_type,
    doi_ISBN,
    resource_journal,
    resource_volume_pages
  )
# Write the combined outlier data to a CSV file
write.csv(final_data, "Trait_outlier_data.csv", row.names = FALSE)


### b) Detect wrong category names for Categorical traits

##For categorical traits, we will check whether the values in the database are among the accepted ones

#1) Asexual reproduction mode
Asexual.reproduction.mode$value <- as.factor(Asexual.reproduction.mode$value)
#We manually introduce the expected values (accepted values)
Accepted_values1 <- c("autotomy", 
                      "budding", 
                      "colony fission", 
                      "colony fragmentation", 
                      "fragmentation", 
                      "parthenogenetic", 
                      "polyp bail-out", 
                      "production of runners", 
                      "translocation")

#1st) We look whether there is any level that is not in the reproductive_methods
missing_levels <- setdiff(levels(Asexual.reproduction.mode$value), Accepted_values1)
if (length(missing_levels) > 0) {
  cat("The following levels are not in reproductive_methods:", paste(missing_levels, collapse = ", "), "\n")
} else {
  cat("All levels match the reproductive_methods.\n")
}
missing_levels

#2nd) We check that all levels in the database match the accepted values
levels_match <- all(levels(Asexual.reproduction.mode$value) %in% Accepted_values1)
levels_match

#2) Branching architecture (unrestricted descriptive categories)
#3) Climate zone 

Climate.zone$value <- as.factor(Climate.zone$value)
#We manually introduce the expected values (accepted values)
Accepted_values3 <- c("polar","temperate","tropical")

#1st) We look whether there is any level that is not in the reproductive_methods
missing_levels3 <- setdiff(levels(Climate.zone$value), Accepted_values3)
if (length(missing_levels3) > 0) {
  cat("The following levels are not in reproductive_methods:", paste(missing_levels3, collapse = ", "), "\n")
} else {
  cat("All levels match the reproductive_methods.\n")
}
missing_levels3

#2nd) We check that all levels in the database match the accepted values
levels_match3 <- all(levels(Climate.zone$value) %in% Accepted_values3)
levels_match3

#4) Coordination of polyps (unrestricted categories)

#5) Depth zone

Depth.zone$value <- as.factor(Depth.zone$value)
#We manually introduce the expected values (accepted values)
Accepted_values5 <- c("deep","moderately deep","shallow")

#1st) We look whether there is any level that is not in the reproductive_methods
missing_levels5 <- setdiff(levels(Depth.zone$value), Accepted_values5)
if (length(missing_levels5) > 0) {
  cat("The following levels are not in reproductive_methods:", paste(missing_levels5, collapse = ", "), "\n")
} else {
  cat("All levels match the reproductive_methods.\n")
}
missing_levels5

#2nd) We check that all levels in the database match the accepted values
levels_match5 <- all(levels(Depth.zone$value) %in% Accepted_values5)
levels_match5

#6) Food sources (unrestricted categories)

#7) Growth form (unrestricted categories)

#8) Growth_form basic
Growth.form_.basic$value <- as.factor(Growth.form_.basic$value)
#We manually introduce the expected values (accepted values)
Accepted_values8 <- c(
  "arborescent",
  "encrusting",
  "encrusting (stoloniferous)",
  "erect branched (tangled)",
  "erect unbranched (capitate)",
  "erect unbranched (simple)",
  "feather-like",
  "horizontal branched (with axial skeletal support)",
  "massive",
  "membranous (sheet-like)",
  "solitary/pseudosolitary",
  "sub-massive (branched)",
  "sub-massive (columnar)",
  "sub-massive (foliose)",
  "sub-massive (lobate/digitate)",
  "sub-massive (plate-like)"
)

#1st) We look whether there is any level that is not in the reproductive_methods
missing_levels8 <- setdiff(levels(Growth.form_.basic$value), Accepted_values8)
if (length(missing_levels8) > 0) {
  cat("The following levels are not in reproductive_methods:", paste(missing_levels8, collapse = ", "), "\n")
} else {
  cat("All levels match the reproductive_methods.\n")
}
missing_levels8

#2nd) We check that all levels in the database match the accepted values
levels_match8 <- all(levels(Growth.form_.basic$value) %in% Accepted_values8)
levels_match8

#9) Growth form_detailed

Growth.form_detailed$value <- as.factor(Growth.form_detailed$value)
#We manually introduce the expected values (accepted values)

Accepted_values9<-c(
  "arborescent (with domed capitula)",
  "arborescent (with terminal branching)",
  "arborescent lobate",
  "bottlebrush",
  "branched planar (multi-plane)",
  "branched planar",
  "bushy/tree-like (TYPE rigid/semirigid)",
  "bushy/tree-like (TYPE soft)",
  "digitiform",
  "encrusting (digitate/lobate/folded)",
  "encrusting (flat)",
  "encrusting (stoloniferous)",
  "encrusting (with verrucae)",
  "erect branched tangled (stoloniferous)",
  "erect branched tangled (with unconsolidated axis)",
  "erect unbranched (pen-like)",
  "erect unbranched (with disc-shaped capitulum)",
  "erect unbranched (with domed capitulum)",
  "erect unbranched (with umbellate cluster of polyps)",
  "feather-like",
  "flagelliform/filiform",
  "horizontal branched (with axial skeletal support)",
  "massive (hemispherical)",
  "membranous (sheet-like)",
  "solitary/pseudosolitary",
  "sparsely branched",
  "sub-massive (branched)",
  "sub-massive (columnar)",
  "sub-massive (foliose)",
  "sub-massive (lobate/digitate)",
  "sub-massive (plate-like)"
)

#1st) We look whether there is any level that is not in the reproductive_methods
missing_levels9 <- setdiff(levels(Growth.form_detailed$value), Accepted_values9)
if (length(missing_levels9) > 0) {
  cat("The following levels are not in reproductive_methods:", 
      paste(missing_levels9, collapse = ", "), "\n")
} else {
  cat("All levels match the reproductive_methods.\n")
}
missing_levels9

#2nd) We check that all levels in the database match the accepted values
levels_match9 <- all(levels(Growth.form_detailed$value) %in% Accepted_values9)
levels_match9

#10) IUCN Red list category
IUCN.Red.List.category$value <- as.factor(IUCN.Red.List.category$value)
#We manually introduce the expected values (accepted values)
Accepted_values10<-c("CR", "DD", "EN", "LC", "NT", "VU")

#1st) We look whether there is any level that is not in the reproductive_methods
missing_levels10 <- setdiff(levels(IUCN.Red.List.category$value), Accepted_values10)
if (length(missing_levels10) > 0) {
  cat("The following levels are not in reproductive_methods:", 
      paste(missing_levels10, collapse = ", "), "\n")
} else {
  cat("All levels match the reproductive_methods.\n")
}
missing_levels10

#2nd) We check that all levels in the database match the accepted values
levels_match10 <- all(levels(IUCN.Red.List.category$value) %in% Accepted_values10)
levels_match10

#11) Light zone
Light.zone$value <- as.factor(Light.zone$value) 

# We manually introduce the expected values (accepted values)
Accepted_values11 <- c("aphotic", "mesophotic", "photic")

# 1st) We look whether there is any level that is not in the accepted values
missing_levels11 <- setdiff(levels(Light.zone$value), Accepted_values11)

if (length(missing_levels11) > 0) {
  cat("The following levels are not in the accepted values:", paste(missing_levels11, collapse = ", "), "\n")
} else {
  cat("All levels match the accepted values.\n")
}

missing_levels11

#2nd) We check that all levels in the database match the accepted values
levels_match11 <- all(levels(Light.zone$value) %in% Accepted_values11)
levels_match11


#12) Main non-scleritic skeletal component (unrestricted values)

#13) Marine province

Marine.province$value <- as.factor(Marine.province$value)

# We manually introduce the expected values (accepted values)
Accepted_values13 <- c(
  "Agulhas",
  "Andaman",
  "Arctic",
  "Bay of Bengal",
  "Benguela",
  "Central Indian Ocean Islands",
  "Central Polynesia",
  "Cold Temperate Northeast Pacific",
  "Cold Temperate Northwest Atlantic",
  "Cold Temperate Northwest Pacific",
  "Continental High Antarctic",
  "East Central Australian Shelf",
  "Eastern Coral Triangle",
  "Galapagos",
  "Gulf of Guinea",
  "Hawaii",
  "Java Transitional",
  "Juan Fernandez and Desventuradas",
  "Lord Howe and Norfolk Islands",
  "Lusitanian",
  "Magellanic",
  "Marquesas",
  "Marshall, Gilbert, and Ellis Islands",
  "Mediterranean Sea",
  "North Brazil Shelf",
  "Northeast Australian Shelf",
  "Northern European Seas",
  "Northern New Zealand",
  "Northwest Australian Shelf",
  "Red Sea and Gulf of Aden",
  "Sahul Shelf",
  "Scotia Sea",
  "Somali/Arabian",
  "South China Sea",
  "South Kuroshio",
  "Southeast Australian Shelf",
  "Southeast Polynesia",
  "Southern New Zealand",
  "Southwest Australian Shelf",
  "St. Helena and Ascension Islands",
  "Subantarctic Islands",
  "Subantarctic New Zealand",
  "Sunda Shelf",
  "Temperate Australasia",
  "Temperate Northern Pacific",
  "Tristan Gough",
  "Tropical East Pacific",
  "Tropical Northwestern Atlantic",
  "Tropical Northwestern Pacific",
  "Tropical Southwestern Atlantic",
  "Tropical Southwestern Pacific",
  "Warm Temperate Northeast Pacific",
  "Warm Temperate Northern Atlantic",
  "Warm Temperate Northwest Atlantic",
  "Warm Temperate Northwest Pacific",
  "Warm Temperate Southeastern Pacific",
  "Warm Temperate Southwestern Atlantic",
  "West African Transition",
  "West and South Indian Shelf",
  "West Central Australian Shelf",
  "Western Coral Triangle",
  "Western Indian Ocean"
)


# 1st) We look whether there is any level that is not in the accepted values
missing_levels13 <- setdiff(levels(Marine.province$value), Accepted_values13)

if (length(missing_levels13) > 0) {
  cat("The following levels are not in the accepted values:", paste(missing_levels13, collapse = ", "), "\n")
} else {
  cat("All levels match the accepted values.\n")
}
missing_levels13

#2nd) We check that all levels in the database match the accepted values
levels_match13 <- all(levels(Marine.province$value) %in% Accepted_values13)
levels_match13

#14) Marine realm

Marine.realm$value <- as.factor(Marine.realm$value)

# We manually introduce the expected values (accepted values)
Accepted_values14 <- c(
  "Arctic", 
  "Central Indo-Pacific", 
  "Cosmopolitan", 
  "Eastern Indo-Pacific", 
  "Southern Ocean", 
  "Temperate Australasia", 
  "Temperate Northern Atlantic", 
  "Temperate Northern Pacific", 
  "Temperate South America", 
  "Temperate Southern Africa", 
  "Tropical Atlantic", 
  "Tropical Eastern Pacific", 
  "Western Indo-Pacific"
)

# 1st) We look whether there is any level that is not in the accepted values
missing_levels14 <- setdiff(levels(Marine.realm$value), Accepted_values14)

if (length(missing_levels14) > 0) {
  cat("The following levels are not in the accepted values:", paste(missing_levels14, collapse = ", "), "\n")
} else {
  cat("All levels match the accepted values.\n")
}
missing_levels14

#2nd) We check that all levels in the database match the accepted values
levels_match14 <- all(levels(Marine.realm$value) %in% Accepted_values14)
levels_match14

#15) Mode of larval development

Mode.of.larval.development$value <- as.factor(Mode.of.larval.development$value)

# We manually introduce the expected values (accepted values)
Accepted_values15 <- c("brooder","external brooder", "spawner")


# 1st) We look whether there is any level that is not in the accepted values
missing_levels15 <- setdiff(levels(Mode.of.larval.development$value), Accepted_values15)

if (length(missing_levels15) > 0) {
  cat("The following levels are not in the accepted values:", paste(missing_levels15, collapse = ", "), "\n")
} else {
  cat("All levels match the accepted values.\n")
}
missing_levels15

#2nd) We check that all levels in the database match the accepted values
levels_match15 <- all(levels(Mode.of.larval.development$value) %in% Accepted_values15)
levels_match15

#16) Ocean basin

Ocean.basin$value <- as.factor(Ocean.basin$value)

# We manually introduce the expected values (accepted values)
Accepted_values16 <- c("Arctic","Atlantic", "Pacific","Southern Ocean")


# 1st) We look whether there is any level that is not in the accepted values
missing_levels16 <- setdiff(levels(Ocean.basin$value), Accepted_values16)

if (length(missing_levels16) > 0) {
  cat("The following levels are not in the accepted values:", paste(missing_levels16, collapse = ", "), "\n")
} else {
  cat("All levels match the accepted values.\n")
}
missing_levels16

#2nd) We check that all levels in the database match the accepted values
levels_match16 <- all(levels(Ocean.basin$value) %in% Accepted_values16)
levels_match16

##17) Reproductive period (unrestricted categories)

#18) Sexual system

Sexual.system$value <- as.factor(Sexual.system$value)

# We manually introduce the expected values (accepted values)
Accepted_values18 <- c("gonochore","gonochoric with rare hermaphroditic colonies ", "hermaphroditic","mixed sexuality")


# 1st) We look whether there is any level that is not in the accepted values
missing_levels18 <- setdiff(levels(Sexual.system$value), Accepted_values18)

if (length(missing_levels18) > 0) {
  cat("The following levels are not in the accepted values:", paste(missing_levels18, collapse = ", "), "\n")
} else {
  cat("All levels match the accepted values.\n")
}
missing_levels18

#2nd) We check that all levels in the database match the accepted values
levels_match18 <- all(levels(Sexual.system$value) %in% Accepted_values18)
levels_match18

##19) Skeletal rigidity
Skeletal.rigidity$value <- as.factor(Skeletal.rigidity$value)

# We manually introduce the expected values (accepted values)
Accepted_values19 <- c("rigid","semi-rigid", "soft")


# 1st) We look whether there is any level that is not in the accepted values
missing_levels19 <- setdiff(levels(Skeletal.rigidity$value), Accepted_values19)

if (length(missing_levels19) > 0) {
  cat("The following levels are not in the accepted values:", paste(missing_levels19, collapse = ", "), "\n")
} else {
  cat("All levels match the accepted values.\n")
}
missing_levels19

#2nd) We check that all levels in the database match the accepted values
levels_match19 <- all(levels(Skeletal.rigidity$value) %in% Accepted_values19)
levels_match19

#20) Spawning date (unrestricted categories)

#21) Symbiodinium clade (unrestricted categories)

#22) Type of growth

Type.of.growth$value <- as.factor(Type.of.growth$value)

# We manually introduce the expected values (accepted values)
Accepted_values22 <- c("encrusting", 
                       "erect branched", 
                       "erect unbranched", 
                       "horizontal branched", 
                       "horizontal unbranched", 
                       "massive", 
                       "solitary/pseudosolitary")


# 1st) We look whether there is any level that is not in the accepted values
missing_levels22 <- setdiff(levels(Type.of.growth$value), Accepted_values22)

if (length(missing_levels22) > 0) {
  cat("The following levels are not in the accepted values:", paste(missing_levels22, collapse = ", "), "\n")
} else {
  cat("All levels match the accepted values.\n")
}
missing_levels22

#2nd) We check that all levels in the database match the accepted values
levels_match22 <- all(levels(Type.of.growth$value) %in% Accepted_values22)
levels_match22

#23) Type of skeleton

Type.of.skeleton$value <- as.factor(Type.of.skeleton$value)

# We manually introduce the expected values (accepted values)
Accepted_values23 <- c(
  "Axial-like layer (i.e., unconsolidated axis) formed by sclerites, with coenenchyme subdivided into medulla and cortex",
  "Axis of consolidated but unfused sclerites",
  "Axis of fused sclerites cemented together within a solid calcareous matrix",
  "Axis of partially fused sclerites embedded in protein",
  "Axis of unfused sclerites embedded in protein, with hollow cross-chambered central core",
  "Hard aragonite skeleton",
  "Jointed axis with proteinaceous nodes that lack sclerites alternating with solid, non-scleritic calcareous internodes",
  "Jointed axis with solid calcareous inter-nodes (fused or consolidated sclerites within a calcite matrix) alternating with proteinaceous nodes with sclerites",
  "Lack skeletal structure",
  "No axis, calcified structures (e.g., polyp tubes, branches and/or membranes) derived of fused sclerites",
  "No axis, longitudinally disposed sclerites surround canals to form rigid branches",
  "No axis, no sclerites",
  "No axis, no sclerites, tissue covered by thin periderm or cuticular envelope",
  "No axis, unconsolidated sclerites embedded in tissue",
  "Proteinaceous axis with hollow, cross chambered central core and non-scleritic calcareous loculi",
  "Proteinaceous axis with hollow, cross-chambered central core (unmineralized)",
  "Proteinaceous axis with hollow, cross-chambered central core containing non-scleritic carbonate hydroxylapatite (CHAp)",
  "Solid (scleroproteinaceous) axis with concentrically-layered non-scleritic aragonite",
  "Solid (scleroproteinaceous) axis with concentrically-layered non-scleritic calcareous material",
  "Solid (scleroproteinaceous) axis with radially-layered non-scleritic calcareous material",
  "Solid proteinaceous axis (unmineralized), colony lacking sclerites or other calcareous material"
)


# 1st) We look whether there is any level that is not in the accepted values
missing_levels23 <- setdiff(levels(Type.of.skeleton$value), Accepted_values23)

if (length(missing_levels23) > 0) {
  cat("The following levels are not in the accepted values:", paste(missing_levels23, collapse = ", "), "\n")
} else {
  cat("All levels match the accepted values.\n")
}
missing_levels23

#2nd) We check that all levels in the database match the accepted values
levels_match23 <- all(levels(Type.of.skeleton$value) %in% Accepted_values23)
levels_match23

#24) Water clarity preference

Water.clarity.preference$value <- as.factor(Water.clarity.preference$value)

# We manually introduce the expected values (accepted values)
Accepted_values24 <- c("both","clear","turbid")


# 1st) We look whether there is any level that is not in the accepted values
missing_levels24 <- setdiff(levels(Water.clarity.preference$value), Accepted_values24)

if (length(missing_levels24) > 0) {
  cat("The following levels are not in the accepted values:", paste(missing_levels24, collapse = ", "), "\n")
} else {
  cat("All levels match the accepted values.\n")
}
missing_levels24

#2nd) We check that all levels in the database match the accepted values
levels_match24 <- all(levels(Water.clarity.preference$value) %in% Accepted_values24)
levels_match24

##25) Wave exposure preference


Wave.exposure.preference$value <- as.factor(Wave.exposure.preference$value)

# We manually introduce the expected values (accepted values)
Accepted_values25 <- c("both","exposed","protected")


# 1st) We look whether there is any level that is not in the accepted values
missing_levels25 <- setdiff(levels(Wave.exposure.preference$value), Accepted_values25)

if (length(missing_levels25) > 0) {
  cat("The following levels are not in the accepted values:", paste(missing_levels25, collapse = ", "), "\n")
} else {
  cat("All levels match the accepted values.\n")
}
missing_levels25

#2nd) We check that all levels in the database match the accepted values
levels_match25 <- all(levels(Wave.exposure.preference$value) %in% Accepted_values25)
levels_match25


###c) Duplicated rows detection 

# Identify all duplicates (including the first occurrence)
all_duplicates <- df[duplicated(df) | duplicated(df, fromLast = TRUE), ]

# Display all duplicate rows
print(all_duplicates)

#There are zero row duplicates in the dataset.


### d) Ensure that each observation is associated to only one location, resource and species

test1 <- df %>%
  group_by(observation_id) %>%
  summarise(
    UniqueLocations = n_distinct(location_id),
    UniqueResources = n_distinct(resource_id),
    UniqueSpecies = n_distinct(specie_id)
  ) %>%
  filter(UniqueLocations > 1 | UniqueResources > 1 | UniqueSpecies > 1)

print(test1)

# There are zero errors


