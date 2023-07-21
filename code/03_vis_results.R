
library(tidyverse)

source("code/maketraits.R")

color_scale <- c(
  "Underreported"  = "#FF9000",
  "Nonsignificant" = "#858585",
  "Overreported"   = "#3888ff"
)

family_shapes <- c(
  "Nymphalidae" = 21,
  "Papilionidae" = 22,
  "Hesperiidae" = 23,
  "Lycaenidae" = 24,
  "Pieridae" = 25
)

region_colors <- c("East" = "#0ead7f",
                   "West" = "#5a55a1")
region_shapes <- c("East" = "square",
                   "West" = "circle")



#### Visualize results of meta-analyses ####

summary_df_traits <- read_csv("output/summary_df_traits.csv")
summary_df_color <- read_csv("output/summary_df_color.csv")
summary_df_feature <- read_csv("output/summary_df_feature.csv")

traitplot <- summary_df_traits %>%
  filter(param != "Intercept") %>%
  filter(region == "joint") %>%
  mutate(param = recode(param, "wingspanscaled" = "Wingspan",
                        "genusIDratescaled" = "Ease of ID",
                        "colorDiversityscaled" = "Color diversity",
                        "eButterflylogcountscaled" = "Prevalence",
                        "featureDiversityscaled" = "Feature diversity")) %>%
  ggplot(aes(paste0("\t", param), Estimate, ymin = `l-95% CI`, ymax = `u-95% CI`,
             shape = sign(`l-95% CI`) == sign(`u-95% CI`))) +
  geom_pointrange(position = position_dodge(width = 0.5),
                  show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  scale_shape_manual(values = c(1, 19)) +
  coord_flip() +
  # scale_color_manual("", values = region_colors) +
  # scale_shape_manual("", values = region_shapes) +
  xlab("") + ylab("Estimated effect (95%CI)") +
  ggtitle("A. Traits")

colorplot <- summary_df_color %>%
  filter(param != "Intercept") %>%
  ggplot(aes(str_to_title(param), Estimate,
             ymin = `l-95% CI`, ymax = `u-95% CI`,
             #color = str_to_title(param),
             shape = sign(`l-95% CI`) == sign(`u-95% CI`))) +
  geom_pointrange(position = position_dodge(width = 0.5),
                  show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  coord_flip() +
  scale_shape_manual(values = c(1, 19)) +
  # scale_color_manual("", values = region_colors) +
  # scale_shape_manual("", values = region_shapes) +
  xlab("") + ylab("Estimated effect (95%CI)") +
  ggtitle("B. Wing colors") #+
  # scale_color_manual(values = c(
  #   "Yellow" = "#ffd000",
  #   "White" = "#bbbbbb",
  #   "Gray" = "#444444",
  #   "Black" = "black",
  #   "Red" = "darkred",
  #   "Purple" = "purple",
  #   "Pink" = "#ff7aff",
  #   "Orange" = "orange",
  #   "Green" = "darkgreen",
  #   "Brown" = "brown",
  #   "Blue" = "blue"
  # ))

featureplot <- summary_df_feature %>% 
  filter(param != "Intercept") %>% 
  mutate(param = recode(param, "band_stripe" = "Band/Stripe")) %>% 
  ggplot(aes(str_to_title(param), Estimate, 
             ymin = `l-95% CI`, ymax = `u-95% CI`,
             shape = sign(`l-95% CI`) == sign(`u-95% CI`))) +
  geom_pointrange(position = position_dodge(width = 0.5),
                  show.legend = FALSE) +
  geom_hline(yintercept = 0) +
  theme_minimal() +
  coord_flip() +
  scale_shape_manual(values = c(1, 19)) +
  # scale_color_manual("", values = region_colors) +
  # scale_shape_manual("", values = region_shapes) +
  xlab("") + ylab("Estimated effect (95%CI)") +
  ggtitle("C. Wing features")


effects_arranged <- gridExtra::grid.arrange(
  traitplot, colorplot, featureplot,
  heights = c(1,1.5,1)
)


ggsave(effects_arranged,
       filename = "plots/fig3_effects_plot.jpg", width = 4, height = 8)


#### Make a phylogenetic tree ####
# install.packages("BiocManager")
# BiocManager::install("ggtree")
library(ggtree)
library(treeio)

ct <- 0 

this_inds <- inds_wtraits %>% 
  group_by(species) %>% 
  mutate(specs_in_both = n() > 1) %>% 
  ungroup() %>% 
  filter(!specs_in_both | region == "East")

tree<- read.tree("data/ultra_bi.tre") ## use this one
source("code/helper_fns.R")
group_info <- data.frame(
  id = get.tree(tree)$tip.label,
  species = lapply(gsub("_", " ", tree$tip.label), dropFamily) %>% 
    unlist()
) %>% 
  # left_join(read_csv("intermediate/species_taxonomy.csv")) %>% 
  left_join(this_inds)


groups_list <- list(
  Nymphalidae =  group_info$id[group_info$family == "Nymphalidae"],
  Papilionidae = group_info$id[group_info$family == "Papilionidae"],
  Hesperiidae =  group_info$id[group_info$family == "Hesperiidae"],
  Lycaenidae =   group_info$id[group_info$family == "Lycaenidae"],
  Pieridae =     group_info$id[group_info$family == "Pieridae"]
)

# tree <- tidytree::groupOTU(tree, 
#                            groups_list,
#                            group_name = "family")

p <- ggtree(tree, layout = "circular")
p <- p %<+% group_info + 
  geom_tippoint(aes(shape = family, color = sig_type, fill = sig_type)) +
  scale_color_manual("Significance", values = color_scale) +
  scale_fill_manual("Significance", values = color_scale) +
  scale_shape_manual("Family", values = family_shapes) +
  theme(
    panel.background = element_rect(fill='transparent'), #transparent panel bg
    plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    panel.grid.major = element_blank(), #remove major gridlines
    panel.grid.minor = element_blank(), #remove minor gridlines
    legend.background = element_rect(fill='transparent'), #transparent legend bg
    legend.box.background = element_rect(fill='transparent') #transparent legend panel
  )

# get_legend <- function(myggplot) {
#   tmp <- ggplot_gtable(ggplot_build(myggplot))
#   leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
#   legend <- tmp$grobs[[leg]]
#   return(legend)
# }

ggsave(p, filename = "plots/species_phyloplot.png", 
       width = 6, height = 8, bg='transparent')



#### Visualize the regions we used & data ####
library(sf)

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
          "Alberta","British Columbia")

grid_cts <- bind_rows(
  read_csv('intermediate/gridcts_east.csv') %>% mutate(region = "East"),
  read_csv('intermediate/gridcts_west.csv') %>% mutate(region = "Wast")
)

grid_locs <- bind_rows(
  read_csv('intermediate/grid_summary_east.csv') %>% mutate(region = "East"),
  read_csv('intermediate/grid_summary_west.csv') %>% mutate(region = "Wast")
)

grid_n <- grid_cts %>% 
  distinct(region, dataset, gridcell, year, ncelltotal) %>% 
  mutate(dataset = recode(dataset, "iNat" = "iNaturalist")) %>% 
  group_by(region, dataset, gridcell) %>% 
  summarize(log_n = log(sum(ncelltotal))) %>% 
  left_join(grid_locs[, c("gridcell", "region", "longitude", "latitude")])


states_provinces <- bind_rows(
  geodata::gadm(country="USA", level=1, path = "data") %>% st_as_sf(),
  geodata::gadm(country="CAN", level=1, path = "data") %>% st_as_sf()
) %>% 
  filter(NAME_1 %in% c(east, west)) %>% 
  mutate(region = ifelse(NAME_1 %in% east, "East", "West")) %>% 
  st_transform(crs = st_crs("+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +datum=NAD83 +units=m +no_defs")) %>% 
  st_simplify()

east_region <- states_provinces %>% 
  filter(region == "East") %>% 
  st_simplify() %>% 
  st_union()
west_region <- states_provinces %>% 
  filter(region == "West") %>% 
  st_simplify() %>% 
  st_union()
  



datpt <- st_as_sf(
  grid_n, coords = c("longitude", "latitude"),
  crs = "+proj=longlat"
)

region_plot <- ggplot() +
  geom_sf(data = states_provinces, col = "#777777", linetype = 2,
          fill = "white") +
  geom_sf(data = east_region, fill = 'transparent', col = "#000000") +
  geom_sf(data = west_region, fill = 'transparent', col = "#000000") +
  geom_sf(data = datpt, aes(col = log_n), size = 0.2) +
  theme_minimal() +
  theme(legend.position = "bottom") +
  facet_wrap(~dataset) +
  scale_color_viridis_c("Log number of observations", option = "plasma")


ggsave("plots/fig1_regions.jpg", region_plot, 
       width = 6, height = 4)


#### Supplemental: make a plot of overreporting surface for one species ####
warning("this will only work if the species' GAM has already been estimated")
pred <- read_csv("intermediate/ind_results/predicted_indices_Cupido comyntas_east.csv") %>% 
  select(longitude, latitude, year, dataset, index, se) %>% 
  pivot_wider(names_from = "dataset", values_from = c("index", "se")) %>% 
  mutate(difference = index_eButterfly - index_iNaturalist)
predpt <- st_as_sf(
  pred, coords = c("longitude", "latitude"),
  crs = "+proj=longlat"
)


diffplot <- ggplot() +
  geom_sf(data = states_provinces[states_provinces$region == "East", ],
          fill = "white") +
  geom_sf(data = predpt[predpt$year %in% 2013:2017,], aes(color = difference),
          size = 0.5) +
  # scale_color_viridis_c() +
  facet_grid(~year) +
  theme_void() +
  scale_color_gradient2(
    "Difference",
    mid = "#cccccc", midpoint = 0,
    low = color_scale[1], high = "#024dbd") +
  ggtitle("Estimated OI of Cupido comyntas")

ebutterfly_errorplot <- ggplot() +
  geom_sf(data = states_provinces[states_provinces$region == "East", ],
          fill = "white") +
  geom_sf(data = predpt[predpt$year %in% 2013:2017,], aes(color = se_eButterfly),
          size = 0.5) +
  # scale_color_viridis_c() +
  facet_grid(~year) +
  scale_color_viridis_c("SE")+
  theme_void() +
  ggtitle("Estimated SE in eButterfly reporting rate of Cupido comyntas")

iNaturalist_errorplot <- ggplot() +
  geom_sf(data = states_provinces[states_provinces$region == "East", ],
          fill = "white") +
  geom_sf(data = predpt[predpt$year %in% 2013:2017,], aes(color = se_iNaturalist),
          size = 0.5) +
  # scale_color_viridis_c() +
  facet_grid(~year) +
  scale_color_viridis_c("SE") +
  theme_void() +
  ggtitle("Estimated SE in eButterfly reporting rate of Cupido comyntas")


ggsave(cowplot::plot_grid(
  diffplot, ebutterfly_errorplot, iNaturalist_errorplot,
  ncol = 1
), filename = "plots/cupido_vis.jpg", width = 9, height = 7)



### Make supplemental tables

output_files <- list.files("output", pattern = "summary_df_", full.names = TRUE)

for (i in 1:length(output_files)) {
  thistbl <- read_csv(output_files[i])
  
  thistbl %>% 
    mutate(`95% CI` = paste0("(", round(`l-95% CI`, 3), ", ", round(`u-95% CI`, 3), ")")) %>% 
    mutate(Estimate = round(Estimate, 3),
           Est.Error = round(Est.Error, 3)) %>% 
    select(Parameter = param,
           Estimate,
           `Std. error` = Est.Error,
           `95% CI`) %>% 
    write_csv(gsub(pattern = "summary_df_",
                   replacement = "supplemental_tbl_",
                   output_files[i]))
}
