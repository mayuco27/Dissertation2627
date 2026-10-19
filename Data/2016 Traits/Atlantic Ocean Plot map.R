library(leaflet)
library(sf)

# 1. Clean the dataset by explicitly dropping the spatial geometry metadata
# This converts it into a pure, clean data frame table
plot_ready_data <- result_df %>% st_drop_geometry()

# 2. Build the map using the correct spelling: "longitude"
Atlantic_Coordinates_map <- leaflet(data = plot_ready_data) %>%
  addTiles() %>% 
  addCircleMarkers(
    lng = plot_ready_data$longitude,  # Spelled exactly as shown in your summary
    lat = plot_ready_data$latitude,   
    radius = 5,
    color = "blue",
    stroke = FALSE, 
    fillOpacity = 0.7,
    popup = paste("Lat:", plot_ready_data$latitude, "<br>Lon:", plot_ready_data$longitude)
  )
# 3. Save the widget
install.packages("htmlwidgets")
library(htmlwidgets)
saveWidget(Atlantic_Coordinates_map, file = "atlantic_ocean_map.html", selfcontained = TRUE)

