# Amphibious fishes follow different assembly rules across Neotropical temporary pools

This repository contains the data and R code used in the manuscript **“Amphibious fishes follow different assembly rules across Neotropical temporary pools.”**

## Repository contents

### `DATABASE.xlsx`

Compiled and taxonomically standardized database used in the analyses. It contains:

- assemblage records;
- sampling-unit and dataset identifiers;
- geographic coordinates and sampling information;
- species taxonomy and amphibious-lifestyle classification;
- pool dimensions and environmental variables used as model inputs.

### `sites_hii_2020.RData`

R object containing values of the **2020 Human Influence Index (HII)** extracted for the unique sampling coordinates. These values are merged with the main database in the analysis code and used as the spatial proxy for anthropogenic pressure.

### `code.R`

Complete R workflow used to:

- merge the 2020 HII values with sampling coordinates;
- calculate assemblage-level richness and abundance metrics;
- calculate amphibious richness and abundance proportions;
- prepare environmental predictors;
- fit the isolation and habitat-availability generalized linear mixed models;
- reproduce the main and supplementary figures.

## Data provenance

The database was assembled through a collaborative effort involving researchers and members of governmental environmental agencies working across the Neotropical region. Detailed dataset provenance, sampling information, collector attribution, and species-level references supporting amphibious classifications are provided in the manuscript and its Data Availability Statement.

## Citation

The manuscript is currently in preparation. Citation information will be added following publication.
