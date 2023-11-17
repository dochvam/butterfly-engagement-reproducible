library(tidyverse)


color_scale <- c(
  "Underreported"  = "#FF9000",
  "Nonsignificant" = "#858585",
  "Overreported"   = "#3888ff"
)

region_colors <- c("east" = "#0ead7f",
                   "west" = "#5a55a1")
region_shapes <- c("east" = "square",
                   "west" = "circle")

#### Read in overreporting indices and join to trait data ####

inds <- read_csv("output/estimated_indices_fromGAMs.csv")

traits_combined <- readxl::read_xlsx("data/SpeciesListForTraits.xlsx", na = "NA")


# flight_traits <- read_csv("data/flightStyleConsensusTraits.csv")

taxonomy <- read_csv("data/iNat_taxonomy.csv") %>%
  select(-inatVerbatim) %>%
  distinct() %>%
  filter(inat_genus != "Euremia") %>%  # Manually resolve problem cases
  filter(!(inat_genus == "Hemiargus" & species == "Echinargus isola"))

recognizability <- read_csv("data/butterfly_recognizability_bygenus.csv") %>%
  select(inat_genus = Genus, nAll, nResearchGrade) %>%
  mutate(genus_IDrate = nResearchGrade / nAll)

prevalence <- read_csv("data/all_species_counts_ebutterfly.csv")

inds_wtraits <- 
  left_join(inds, traits_combined, by = c("species" = "eButterfly_species")) %>% 
  left_join(taxonomy, by = "species") %>%
  left_join(recognizability, by = "inat_genus") %>%
  left_join(prevalence, by = "species") %>%
  mutate(colorDiversity_scaled   = as.numeric(scale(colorDiversity)),
         featureDiversity_scaled = as.numeric(scale(featureDiversity)),
         wingspan_scaled         = as.numeric(scale(aveWingspan)),
         genus_IDrate_scaled     = as.numeric(scale(genus_IDrate)),
         eButterfly_logcount_scaled = as.numeric(scale(log(n)))
  ) %>%
  mutate(ymin = difference + qnorm(0.025)*diff_SE,
         ymax = difference + qnorm(0.975)*diff_SE,
         pval_uncorr = 2*ifelse(pnorm(difference / diff_SE)>0.5, 1-pnorm(difference / diff_SE), pnorm(difference / diff_SE)),
  ) %>%
  mutate(
    pval_adj = p.adjust(pval_uncorr, "fdr"),
    sig_type = ifelse(pval_adj < 0.05,
                      ifelse(difference < 0, "Underreported", "Overreported"),
                      "Nonsignificant")
  ) %>%
  mutate(region = recode(region, "east" = "East", "west" = "West")) %>%
  mutate(sig_type = factor(sig_type, levels = names(color_scale))) %>%
  arrange(region, -difference) %>%
  mutate(sp_fac = factor(species, levels = unique(species)))



