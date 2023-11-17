library(tidyverse)
library(brms)     ## meta-analysis
## instructions for setting up brms
## https://learnb4ss.github.io/learnB4SS/articles/install-brms.html

color_scale <- c(
  "Underreported"  = "#FF9000",
  "Nonsignificant" = "#858585",
  "Overreported"   = "#3888ff"
)

region_colors <- c("east" = "#0ead7f",
                   "west" = "#5a55a1")
region_shapes <- c("east" = "square",
                   "west" = "circle")

#### Read in overreporting indices with trait data ####
source("code/maketraits.R")

trait_cols <- c("colorDiversity_scaled", "featureDiversity_scaled",
                "wingspan_scaled", "genus_IDrate_scaled", "is_migratory")

# Make a summary of which fields we do and don't have for each species we modeled
species_traits <- inds_wtraits %>% 
  select(species, all_of(trait_cols)) %>% 
  distinct()
colMeans(is.na(species_traits))
cor(species_traits[, trait_cols])


inds_east <- inds_wtraits %>% 
  filter(region == "East") %>% 
  select(species, region, difference, diff_SE, 
         featureDiversity_scaled, colorDiversity_scaled, 
         eButterfly_logcount_scaled,
         wingspan_scaled, genus_IDrate_scaled,
         is_migratory) %>% 
  na.omit()
inds_west <- inds_wtraits %>% 
  filter(region == "West") %>% 
  select(species, region, difference, diff_SE, 
         featureDiversity_scaled, colorDiversity_scaled, 
         eButterfly_logcount_scaled,
         wingspan_scaled, genus_IDrate_scaled,
         is_migratory) %>% 
  na.omit()



#### Traits metaanalysis ####

# # Estimate a brms model separately for each region
# brm_mod_east <- brm(
#   difference | se(diff_SE) ~ 
#     wingspan_scaled + featureDiversity_scaled + 
#     colorDiversity_scaled + genus_IDrate_scaled +
#     eButterfly_logcount_scaled +
#     (1|species), 
#   prior = prior_string("normal(0, 5)", class = "b") +
#     prior_string("cauchy(0, 5)", class = "sd"),
#   data = inds_east, cores = 1) ## could add cores = ...
# 
# brm_mod_west <- brm(
#   difference | se(diff_SE) ~ 
#     wingspan_scaled + featureDiversity_scaled + 
#     colorDiversity_scaled + genus_IDrate_scaled +
#     eButterfly_logcount_scaled +
#     (1|species),
#   prior = prior_string("normal(0, 5)", class = "b") +
#     prior_string("cauchy(0, 5)", class = "sd"),
#   data = inds_west,
#   cores = 1)


inds_both <- bind_rows(inds_east, inds_west) %>% 
  mutate(species_region = as.factor(paste0(region, "_", species)))


brm_mod_both <- brm(
  difference | se(diff_SE) ~ 
    wingspan_scaled + featureDiversity_scaled + 
    colorDiversity_scaled + genus_IDrate_scaled +
    eButterfly_logcount_scaled + 
    is_migratory + 
    (1|species/species_region),
  prior = prior_string("normal(0, 5)", class = "b") +
    prior_string("cauchy(0, 5)", class = "sd"),
  data = inds_both,
  cores = 1)


summary_df_traits <- summary(brm_mod_both)$fixed
summary_df_traits$param <- gsub("[^[:alpha:]]", "", rownames(summary_df_traits))
rownames(summary_df_traits) <- NULL

write_csv(summary_df_traits, "output/summary_df_traits.csv")


#### Meta-analysis by taxonomic family ####

inds_wfam <- left_join(inds_both, taxonomy, by = "species") %>% 
  mutate(
    Hesperiidae  = as.numeric(family == "Hesperiidae"),
    Lycaenidae   = as.numeric(family == "Lycaenidae"),
    Nymphalidae  = as.numeric(family == "Nymphalidae"),
    Papilionidae = as.numeric(family == "Papilionidae"),
    Pieridae     = as.numeric(family == "Pieridae")
  )

brm_mod_both_fam <- brm(
  difference | se(diff_SE) ~ 
    0 + Hesperiidae + Lycaenidae + Nymphalidae + 
    Papilionidae + Pieridae + (1|species/species_region),
  data = inds_wfam,
  prior = prior_string("normal(0, 5)", class = "b") +
    prior_string("cauchy(0, 5)", class = "sd"),
  cores = 1)

summary_df_fam <- summary(brm_mod_both_fam)$fixed
summary_df_fam$param <- rownames(summary_df_fam)
rownames(summary_df_fam) <- NULL

write_csv(summary_df_fam, "output/summary_df_fam.csv")

#### Meta-analysis by colors ####

inds_wtraits$species_region <- paste0(inds_wtraits$species, "_", 
                                      inds_wtraits$region)
brm_mod_both_color <- brm(
  difference | se(diff_SE) ~ 
    black + gray + white + brown + red + orange +
    yellow + green + blue + pink + purple +
    (1|species/species_region),
  prior = prior_string("normal(0, 5)", class = "b") +
    prior_string("cauchy(0, 5)", class = "sd"),
  data = inds_wtraits,
  cores = 1
)

summary_df_color <- summary(brm_mod_both_color)$fixed
summary_df_color$param <- rownames(summary_df_color)
rownames(summary_df_color) <- NULL
write_csv(summary_df_color, "output/summary_df_color.csv")


#### Meta-analysis by features ####
inds_wtraits$band_stripe <- as.numeric(inds_wtraits$band | inds_wtraits$stripe)
brm_mod_both_feature <- brm(
  difference | se(diff_SE) ~ 
    eyespot + checker + tail + band_stripe + spot +
    (1|species/species_region),
  prior = prior_string("normal(0, 5)", class = "b") +
    prior_string("cauchy(0, 5)", class = "sd"),
  data = inds_wtraits,
  cores = 1)

summary_df_feature <- summary(brm_mod_both_feature)$fixed
summary_df_feature$param <- rownames(summary_df_feature)
rownames(summary_df_feature) <- NULL
write_csv(summary_df_feature, "output/summary_df_feature.csv")
