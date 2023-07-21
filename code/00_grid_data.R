library(tidyverse)
library(sf)

color_scale <- c("Underreported"  = "#FF9000",
                 "Nonsignificant" = "#858585",
                 "Overreported"   = "#3888ff")
region_colors <- c("east" = "#0ead7f",
                   "west" = "#5a55a1")
region_shapes <- c("east" = "square",
                   "west" = "circle")

# Define taxonomic families that are butterflies, to exclude moths in inat
butterfly_families <- c(
  "Papilionidae", "Pieridae", "Lycaenidae", "Hesperiidae",
  "Riodinidae", "Libytheidae", "Nymphalidae"
)
east <- c("Alabama", "Arkansas", "Connecticut", "Delaware", 
          "District of Columbia", "Florida", "Georgia", "Illinois", "Indiana", 
          "Iowa", "Kansas", "Kentucky", "Louisiana", "Maine", "Maryland", 
          "Massachusetts", "Michigan", "Minnesota", "Mississippi", "Missouri", 
          "Nebraska", "New Hampshire", "New Jersey", "New York", 
          "North Carolina", "North Dakota", "Ohio", "Oklahoma",
          "Pennsylvania", "Rhode Island", "South Carolina", "South Dakota", 
          "Tennessee", "Texas", "Vermont", "Virginia", "West Virginia", 
          "Wisconsin", "Manitoba", "New Brunswick", 
          "Newfoundland and Labrador", "Nova Scotia", "Ontario", 
          "Prince Edward Island", "Québec", "Saskatchewan")

west <- c("Arizona", "California", "Colorado", "Idaho", "Montana", "Nevada",
          "New Mexico", "Oregon", "Utah", "Washington", "Wyoming",
          "Alberta", "British Columbia")


# Make a table to match inat and eButterfly taxonomy
species_tbl_raw <- readxl::read_xlsx("data/SpeciesList.xlsx")
species_tbl_list <- list()
ct <- 0
for (i in 1:nrow(species_tbl_raw)) {
  if (!is.na(species_tbl_raw$GBIF_species[i])) {
    ct <- ct + 1
    species_tbl_list[[ct]] <- data.frame(
      eButterfly_species = species_tbl_raw$eButterfly_species[i],
      iNaturalist_species = str_split(species_tbl_raw$GBIF_species[i], pattern = ", ")[[1]],
      type = "species"
    )
  }
  if (!is.na(species_tbl_raw$GBIF_verbatimScientificName[i])) {
    ct <- ct + 1
    species_tbl_list[[ct]] <- data.frame(
      eButterfly_species = species_tbl_raw$eButterfly_species[i],
      iNaturalist_species = str_split(species_tbl_raw$GBIF_verbatimScientificName[i], pattern = ", ")[[1]],
      type = "verbatim"
    )
  }
}
species_tbl <- bind_rows(species_tbl_list)

unzip("data/eButterfly_alldat.zip", 
      exdir = "data")

# Read in eButterfly data
ebud_dat <- read_csv("data/eButterfly_alldat.csv") %>%
  filter(`Checklist Completed` == "Yes") %>% 
  filter(!is.na(Latitude), !is.na(Longitude)) %>%
  mutate(species = paste(Genus, Species),
         year = lubridate::year(`Date Observed`)) %>% 
  select(OccuranceID, decimalLongitude = Longitude, 
         decimalLatitude = Latitude, stateProvince = `Province/State`,
         species, eventDate = `Date Observed`, year) %>% 
  filter(year %in% 2000:2021) %>%
  mutate(stateProvince = recode(stateProvince, "Quebec" = "Québec")) %>% 
  filter(decimalLatitude < 65, decimalLatitude > 22,
         decimalLongitude > -140, decimalLongitude < -50) %>% 
  filter(species %in% species_tbl$eButterfly_species)

inat_dat <- read_tsv(stop(
  "Please download the iNaturalist data from GBIF here:
      https://doi.org/10.15468/dl.rhmxtn
  Then modify this filepath to point to the file"
)) %>%
# inat_dat <- read_tsv("data/0436135-210914110416597.csv") %>%
  filter(occurrenceStatus == "PRESENT",
         family %in% butterfly_families) %>%
  mutate(year = lubridate::year(eventDate)) %>% 
  filter(year %in% 2000:2021) %>% 
  mutate(stateProvince = recode(stateProvince, "Quebec" = "Québec")) %>% 
  mutate(eButterfly_species = NA)

# Correct the iNaturalist names to match eButterfly taxonomy
for (i in 1:nrow(species_tbl)) {
  if (species_tbl$type[i] == "verbatim") {
    
    inat_dat$eButterfly_species[
      inat_dat$verbatimScientificName == species_tbl$iNaturalist_species[i]
    ] <- species_tbl$eButterfly_species[i]
    
  } else if (species_tbl$type[i] == "species") {
    
    inat_dat$eButterfly_species[
      inat_dat$species == species_tbl$iNaturalist_species[i]
    ] <- species_tbl$eButterfly_species[i]
    
  }
}

inat_dat <- inat_dat %>% 
  filter(!is.na(eButterfly_species)) %>% 
  mutate(raw_species = species) %>% 
  mutate(species = eButterfly_species,
         inatVerbatim = verbatimScientificName)

inat_dat %>% 
  distinct(species, family, inatVerbatim) %>% 
  mutate(genus = unlist(lapply(species, function(x) {
    str_split(x, " ")[[1]][1]
  }))) %>%
  mutate(inat_genus = unlist(lapply(inatVerbatim, function(x) {
    str_split(x, " ")[[1]][1]
  }))) %>%
  write_csv("intermediate/species_taxonomy.csv")

for (region in c("east", "west")) {
  
  if (region == "east") {
    ebud_reg <- ebud_dat %>% 
      filter(stateProvince %in% east) %>% 
      filter(decimalLongitude > -111)
    inat_reg <- inat_dat %>% 
      filter(stateProvince %in% east) %>% 
      filter(decimalLongitude > -111)
  } else if (region == "west") {
    ebud_reg <- ebud_dat %>% 
      filter(stateProvince %in% west) %>% 
      filter(decimalLongitude < -100)
    inat_reg <- inat_dat %>% 
      filter(stateProvince %in% west) %>% 
      filter(decimalLongitude < -100)
  } else {
    stop("Invalid region")
  }
  
  
  ebud_locs <- count(ebud_reg, decimalLongitude, decimalLatitude)
  inat_locs <- count(inat_reg, decimalLongitude, decimalLatitude)
  
  ebud_pts <- 
    st_as_sf(ebud_locs, 
             crs = st_crs("+proj=longlat"),
             coords = c("decimalLongitude", "decimalLatitude")) %>% 
    st_transform(st_crs("EPSG:3573"))
  inat_pts <- 
    st_as_sf(inat_locs, 
             crs = st_crs("+proj=longlat"),
             coords = c("decimalLongitude", "decimalLatitude")) %>% 
    st_transform(st_crs("EPSG:3573"))
  
  # Regular hexagons with SHORT DIAMETER 40km
  grid_poly <- st_make_grid(st_bbox(bind_rows(ebud_pts, inat_pts)) +
                              c(-80000, -80000, 80000, 80000), 
                            cellsize = c(40000), 
                            square = FALSE) %>% 
    sf::st_as_sf() %>% 
    mutate(gridcell = row_number())
  
  
  ebud_gridcell <- unlist(lapply(st_within(ebud_pts, grid_poly), function(x) {
    if (length(x) == 0) return(NA)
    else return(x)
  }))
  ebud_locs$gridcell <- grid_poly$gridcell[ebud_gridcell]
  
  inat_gridcell <- unlist(lapply(st_within(inat_pts, grid_poly), function(x) {
    if (length(x) == 0) return(NA)
    else return(x)
  }))
  inat_locs$gridcell <- grid_poly$gridcell[unlist(inat_gridcell)]
  
  
  ebud_dat_wcell <- left_join(ebud_reg, ebud_locs) %>% 
    mutate(dataset = "eButterfly")
  inat_dat_wcell <- left_join(inat_reg, inat_locs) %>% 
    mutate(dataset = "iNat")
  
  ## Sara trying something
  
  #inat_dat_wcell %>% count(recordedBy, species, year, gridcell) %>% write.csv(., "intermediate/inatEastByUser.csv", row.names = F)
  #inat_dat_wcell %>% count(recordedBy, species, year, gridcell) %>% write.csv(., "intermediate/inatWestByUser.csv", row.names = F)
  
  # Now we have gridcell IDs associated w each 
  # Target output is the following:
  # Each row is a species / cell / year / dataset
  # One column gives counts of that species; 2nd column gives total # obs across species
  
  gridcts <- bind_rows(inat_dat_wcell[, c("species", "year", "dataset", "gridcell")], 
                       ebud_dat_wcell[, c("species", "year", "dataset", "gridcell")]) %>% 
    count(species, year, dataset, gridcell) %>% 
    ungroup()
  
  # Pivot wider to get a column for each dataset count, zero-fill, then pivot back
  gridcts_wide <- gridcts %>% 
    pivot_wider(names_from = "dataset", values_from = "n") %>% 
    mutate(iNat = ifelse(is.na(iNat), 0, iNat),
           eButterfly = ifelse(is.na(eButterfly), 0, eButterfly))
  
  gridcts <- gridcts_wide %>% 
    pivot_longer(cols = c("iNat", "eButterfly"), 
                 names_to = "dataset", values_to = "n") %>%
    group_by(year, dataset, gridcell) %>%
    mutate(ncelltotal = sum(n)) %>% 
    ungroup()
  
  cellyears_to_drop <- gridcts %>% 
    filter(ncelltotal == 0) %>% 
    distinct(gridcell, year) %>% 
    mutate(gcy = paste0("C",gridcell, "_Y", year))
  
  gridcts <- gridcts %>% 
    filter(!paste0("C",gridcell, "_Y", year) %in% cellyears_to_drop$gcy)
  
  
  # we also need lat/long centroids for each gridcell ID. 
  # (Let's save both projection + long/lat coords)
  grid_summary <- st_centroid(grid_poly)
  grid_summary <- bind_cols(grid_summary, st_coordinates(grid_summary)[, 1:2])
  grid_summary <- st_transform(grid_summary, st_crs("+proj=longlat"))
  grid_summary <- st_drop_geometry(bind_cols(grid_summary, st_coordinates(grid_summary)[, 1:2]))
  colnames(grid_summary) <- c("gridcell", "x", "y", "longitude", "latitude")
  
  write_csv(grid_summary, paste0("intermediate/grid_summary_", region, ".csv"))
  write_csv(left_join(gridcts, grid_summary, by = "gridcell"), 
            paste0("intermediate/gridcts_", region, ".csv"))
  
  
  # Produce some summaries of species detections in this region,
  # for use in selecting target species
  spec_summary <- gridcts %>% 
    select(-ncelltotal) %>% 
    pivot_wider(names_from = 'dataset', values_from = 'n') %>% 
    group_by(species) %>% 
    summarize(ncellyear_inat = sum(iNat > 0),
              nobs_inat = sum(iNat),
              ncellyear_ebud = sum(eButterfly > 0),
              nobs_ebud = sum(eButterfly)
    )
  write_csv(spec_summary, paste0("intermediate/spec_summary_", region, ".csv"))
}


# choose species to run
specs_west <- read_csv("intermediate/spec_summary_west.csv") %>% 
  mutate(region = "west")
specs_east <- read_csv("intermediate/spec_summary_east.csv") %>% 
  mutate(region = "east")

target_specs <- bind_rows(specs_west, specs_east) %>% 
  filter(ncellyear_inat > 25 & ncellyear_ebud > 25)
write_csv(target_specs, "intermediate/target_species.csv")
