# Identification ease, wing pattern diversity, and family explain taxonomic bias in butterfly observations

This subdirectory contains input code for reproducing the analyses described in the paper. It includes a few files harmonizing species taxonomy across different data sources; a taxonomic tree, "ultra_bi.tre", used for visualizations; and the eButterfly observation data.

Due to Github file size limits, you will need to manually retrieve the iNaturalist data from GBIF at [this link](https://doi.org/10.15468/dl.rhmxtn) and add it to this directory.

## Descriptions of the data files

Below, we describe each data file provided.

SpeciesList.xlsx\
This .xlsx file contains a list of butterfly species considered for analysis. It is used for resolving taxonomy between iNaturalist and eButterfly. It contains the following columns:

-    eButterfly_species - the name of the species in eButterfly

-   GBIF_species - the name or names of the species in iNaturalist/GBIF, if those species appear in the "species" field

-   GBIF_verbatimScientificName - the name or names of the species in iNaturalist/GBIF, if those species appear in the "verbatimScientificName" field of GBIF data

SpeciesListForTraits.xlsx

SpeciesList.xlsx\
This .xlsx file contains a list of butterfly species considered for analysis. It is used for resolving taxonomy between iNaturalist, eButterfly, and the trait data. It also contains the raw trait data we aggregated for our species. It contains the following columns:

-    eButterfly_species - the name of the species in eButterfly

-   GBIF_species - the name or names of the species in iNaturalist/GBIF, if those species appear in the "species" field

-   GBIF_verbatimScientificName - the name or names of the species in iNaturalist/GBIF, if those species appear in the "verbatimScientificName" field of GBIF data

-   LepTraits_name - the scientific name of the species according to the LepTraits dataset

-   BAMONA_name - the scientific name of the species according to BAMONA

-   aveWingspan - the species' typical wingtip-to-wingtip length, in cm

-   black - whether or not the species has black on its wing

-   gray - whether or not the species has gray on its wing

-   white - whether or not the species has white on its wing

-   brown - whether or not the species has brown on its wing

-   red - whether or not the species has red on its wing

-   orange - whether or not the species has orange on its wing

-   yellow - whether or not the species has yellow on its wing

-   green - whether or not the species has green on its wing

-   blue - whether or not the species has blue on its wing

-   pink - whether or not the species has pink on its wing

-   purple - whether or not the species has purple on its wing

-   eyespot - whether or not the species has an eyespot feature on its wing

-   stripe - whether or not the species has a stripe feature on its wing

-   checker - whether or not the species has a checker feature on its wing

-   tail - whether or not the species has a wing tail

-   band - whether or not the species has a banding feature on its wing

-   spot - whether or not the species has a spot feature on its wing

-   colorDiversity - the number of colors present on the species' wing

-   featureDiversity - the number of features present on the species' wing

all_species_counts_ebutterfly.csv

This file contains the total counts of all the butterfly species in eButterfly at the time this study was conducted. It contains the following columns:

-   species - the scientific name of the species

-   n - the number of detection events of this species in eButterfly

butterfly_recognizability_bygenus.csv

The "recognizability index" for each iNaturalist genus produced for this study. It contains the following columns:

-   Genus - the iNaturalist genus

-   nAll - The total number of observations of this genus in iNaturalist

-   nResearchGrade - The number of research grade observations of this genus in iNaturalist

eButterfly_alldat.zip

A .zip file containing the raw eButterfly data used in this study. It contains the following columns

-   OccurranceID - the unique ID for this observation

-   Family - taxonomic family of the species ID'ed

-   Genus - taxonomic genus

-   Species - specific epithet

-   Scientific Name - binomial scientific name

-   English Name - english common name

-   Number of individuals - How many individuals of the species were observed

-   Country - the country where the observation took place

-   Province/State - the province or state where the observation took place

-   Latitude

-   Longitude

-   Checklist ID - the unique ID for the checklist sampling event

-   Distance Traveled - the distance traveled by the observer during the sampling event

-   Date Observed - the date of sampling

-   Start time - the time of day of sampling

-   Duration - the length of the sampling event

-   Sampling protocol - the type of sampling conducted

-   Checklist completed - indicates whether the checklist event was "complete", meaning all ID'ed species were reported

-   Action - indicates whether latitude and longitude were obscured (for sensitive species)

iNat_taxonomy.csv

Since we use eButterfly species names, we use this file to track which iNaturalist genus each species belongs to, for association with the recognizability data.

-   species - The species' scientific name in eButterfly

-   family - the species' taxonomic family

-   inatVerbatim - the verbatimScientificName from GBIF/iNaturalist

-   genus - the species' eButterfly genus

-   inat_genus - the species' iNaturalist genus

ultra_bi.tre

This file contains the phylogenetic tree used to produce figure 4.
