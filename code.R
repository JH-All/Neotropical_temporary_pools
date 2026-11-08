# Packages --------------------------
library(readxl)
library(tidyverse)
library(purrr)
library(tibble)
library(glmmTMB)
library(performance)
library(DHARMa)
library(broom.mixed)
library(cowplot)
library(patchwork)
library(scales)
library(rnaturalearth)
library(ggspatial)
library(sf)

# Data -------------------------------
species = read_excel("DATABASE.xlsx", sheet = "species")
env = read_excel("DATABASE.xlsx", sheet = "env")
assemblages = read_excel("DATABASE.xlsx", sheet = "assemblage")

# Figure 1 ---------------------------------------------------------------

## World map --------------------------------------------------------------

mundo <- rnaturalearth::ne_countries(
  scale = "medium",
  returnclass = "sf"
) |>
  sf::st_make_valid()

## Temporary-pool locations -----------------------------------------------

locais_sf <- env |>
  filter(
    is.finite(Lat),
    is.finite(Long)
  ) |>
  distinct(
    Lat,
    Long
  ) |>
  st_as_sf(
    coords = c(
      "Long",
      "Lat"
    ),
    crs = 4326
  )

## Latitude lines ---------------------------------------------------------

latitudes <- seq(
  -40,
  20,
  by = 20
)

longitudes_linha <- seq(
  -125,
  -25,
  length.out = 500
)

paralelos <- st_sf(
  latitude = latitudes,
  geometry = st_sfc(
    lapply(
      latitudes,
      function(lat) {
        st_linestring(
          cbind(
            longitudes_linha,
            rep(
              lat,
              length(
                longitudes_linha
              )
            )
          )
        )
      }
    ),
    crs = 4326
  )
)

## Latitude labels --------------------------------------------------------

rotulos_latitude <- tibble(
  Long = -117,
  Lat = latitudes,
  rotulo = case_when(
    latitudes < 0 ~
      paste0(
        abs(latitudes),
        "°S"
      ),
    
    latitudes > 0 ~
      paste0(
        latitudes,
        "°N"
      ),
    
    TRUE ~ "0°"
  )
) |>
  st_as_sf(
    coords = c(
      "Long",
      "Lat"
    ),
    crs = 4326
  )

## Equal-area projection --------------------------------------------------

crs_neotropical <- st_crs(
  paste0(
    "+proj=laea ",
    "+lat_0=-12 ",
    "+lon_0=-70 ",
    "+datum=WGS84 ",
    "+units=m ",
    "+no_defs"
  )
)

## Final map --------------------------------------------------------------

fig1 <- ggplot() +
  
  geom_sf(
    data = mundo,
    fill = "grey70",
    color = "white",
    linewidth = 0.25
  ) +
  
  geom_sf(
    data = paralelos,
    color = "grey25",
    linewidth = 0.35,
    linetype = "dashed",
    alpha = 0.7
  ) +
  geom_sf_text(
    data = rotulos_latitude,
    aes(
      label = rotulo
    ),
    family = "serif",
    fontface = "italic",
    size = 3.2,
    color = "grey20",
    hjust = 0
  ) +
  
  geom_sf(
    data = locais_sf,
    shape = 21,
    fill = "black",
    color = "white",
    stroke = 0.22,
    size = 2.8,
    alpha = 0.80
  ) +
  
  annotation_scale(
    location = "bl",
    width_hint = 0.20,
    line_width = 0.5,
    pad_x = grid::unit(
      0.4,
      "cm"
    ),
    pad_y = grid::unit(
      0.4,
      "cm"
    )
  ) +
  
  annotation_north_arrow(
    location = "tr",
    which_north = "true",
    height = grid::unit(
      1.1,
      "cm"
    ),
    width = grid::unit(
      1.1,
      "cm"
    ),
    pad_x = grid::unit(
      0.4,
      "cm"
    ),
    pad_y = grid::unit(
      0.4,
      "cm"
    ),
    style = north_arrow_minimal
  ) +
  
  coord_sf(
    crs = crs_neotropical,
    default_crs = st_crs(
      4326
    ),
    xlim = c(
      -120,
      -30
    ),
    ylim = c(
      -60,
      35
    ),
    expand = FALSE,
    lims_method = "box",
    datum = NA
  ) +
  
  theme_void(
    base_size = 20,
    base_family = "serif"
  ) +
  
  theme(
    panel.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    plot.background = element_rect(
      fill = "white",
      color = NA
    ),
    
    panel.border = element_rect(
      fill = NA,
      color = "grey30",
      linewidth = 0.5
    ),
    
    plot.margin = margin(
      12,
      14,
      10,
      12
    )
  )

fig1

## Save -------------------------------------------------------------------

ggsave(
  filename = "Figure_1.jpg",
  plot = fig1,
  width = 9,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

# Human Influence Index (2020) ------------------------

# Load previously extracted HII values for 2020

load("sites_hii_2020.RData")

# Match HII values to sampling sites using coordinates

env <- env |>
  left_join(
    sites_hii_2020 |>
      select(
        Lat,
        Long,
        HII_2020
      ),
    by = c("Lat", "Long")
  )

# Check HII assignment

env |>
  summarise(
    N_sites = n(),
    Missing_HII = sum(is.na(HII_2020)),
    Min_HII = min(HII_2020, na.rm = TRUE),
    Max_HII = max(HII_2020, na.rm = TRUE)
  )

# Average taxonomic distinctness -----------------------

# Taxonomic hierarchy: species, genus, family, order.
# Taxonomic distances represent the number of branches
# connecting each pair of species (2, 4, 6, or 8).
# Average taxonomic distinctness is undefined for S < 2.

## Prepare taxonomic classification ---------------------

taxonomia <- species |>
  transmute(
    Species,
    Genus = word(Species, 1),
    Family,
    Order
  ) |>
  distinct(Species, .keep_all = TRUE)

## Match species to assemblage columns ------------------

stopifnot(
  all(colnames(assemblages) %in% taxonomia$Species),
  nrow(assemblages) == nrow(env)
)

taxonomia <- taxonomia |>
  slice(match(colnames(assemblages), Species))

stopifnot(
  !anyNA(taxonomia[, c("Genus", "Family", "Order")])
)

## Calculate pairwise taxonomic distances ---------------

dist_tax <- outer(
  seq_len(nrow(taxonomia)),
  seq_len(nrow(taxonomia)),
  Vectorize(function(i, j) {
    if (i == j) return(0)
    if (taxonomia$Genus[i] == taxonomia$Genus[j]) return(2)
    if (taxonomia$Family[i] == taxonomia$Family[j]) return(4)
    if (taxonomia$Order[i] == taxonomia$Order[j]) return(6)
    return(8)
  })
)

## Calculate average taxonomic distinctness per site ----

calcular_delta <- function(x) {
  
  presentes <- which(!is.na(x) & x > 0)
  
  if (length(presentes) < 2) return(NA_real_)
  
  d <- dist_tax[presentes, presentes, drop = FALSE]
  
  mean(d[upper.tri(d)])
}

env$Average_taxonomic_distinctness <- apply(
  as.matrix(assemblages),
  1,
  calcular_delta
)

## Check calculated values ------------------------------

env |>
  summarise(
    N_sites = n(),
    Missing_distinctness = sum(
      is.na(Average_taxonomic_distinctness)
    ),
    Min_distinctness = min(
      Average_taxonomic_distinctness,
      na.rm = TRUE
    ),
    Mean_distinctness = mean(
      Average_taxonomic_distinctness,
      na.rm = TRUE
    ),
    Max_distinctness = max(
      Average_taxonomic_distinctness,
      na.rm = TRUE
    )
  )

# Figure 2 ---------------------------------------------------------------

## Amphibious species -----------------------------------------------------

stopifnot(
  nrow(env) == nrow(assemblages),
  nrow(species) == ncol(assemblages),
  identical(species$Species, colnames(assemblages)),
  all(species$Amphibious_lifestype %in% c("Y", "N"))
)

resumo_amphibious <- species |>
  summarise(
    Total_species = n(),
    Amphibious_species = sum(Amphibious_lifestype == "Y"),
    Amphibious_percentage = 100 * mean(Amphibious_lifestype == "Y")
  )

resumo_amphibious


## Recalculate richness and abundance ------------------------------------

matriz_abundancia <- as.matrix(assemblages)

stopifnot(
  is.numeric(matriz_abundancia),
  !anyNA(matriz_abundancia),
  all(matriz_abundancia >= 0),
  all(matriz_abundancia == round(matriz_abundancia))
)

matriz_PA <- (matriz_abundancia > 0) * 1L

colunas_amphibious <- species$Amphibious_lifestype == "Y"

S_total <- rowSums(matriz_PA)

S_amphibious_count <- rowSums(
  matriz_PA[, colunas_amphibious, drop = FALSE]
)

abundance_total <- rowSums(matriz_abundancia)

abundance_amphibious_count <- rowSums(
  matriz_abundancia[, colunas_amphibious, drop = FALSE]
)

stopifnot(
  all(S_total > 0),
  all(S_amphibious_count <= S_total),
  all(abundance_amphibious_count <= abundance_total)
)

env <- env |>
  mutate(
    richness = as.numeric(S_total),
    
    S_amphibious_count = as.integer(
      S_amphibious_count
    ),
    
    S_amphibious = 100 *
      S_amphibious_count / richness,
    
    abundance = if_else(
      Data_Type == "Abundance",
      as.numeric(abundance_total),
      NA_real_
    ),
    
    abundance_amphibious_count = as.integer(
      abundance_amphibious_count
    ),
    
    abundance_amphibious = if_else(
      Data_Type == "Abundance" & abundance_total > 0,
      100 * abundance_amphibious_count / abundance_total,
      NA_real_
    )
  )

checagem_respostas <- env |>
  summarise(
    N_sites = n(),
    Richness_min = min(richness),
    Richness_max = max(richness),
    N_abundance = sum(is.finite(abundance)),
    Abundance_min = min(abundance, na.rm = TRUE),
    Abundance_max = max(abundance, na.rm = TRUE)
  )

checagem_respostas


## Check updated predictors ----------------------------------------------

checagem_preditores <- env |>
  summarise(
    N_sites = n(),
    HII_missing = sum(is.na(HII_2020)),
    HII_zero = sum(HII_2020 == 0, na.rm = TRUE),
    TaxDistinct_missing = sum(
      is.na(Average_taxonomic_distinctness)
    ),
    TaxDistinct_min = min(
      Average_taxonomic_distinctness,
      na.rm = TRUE
    ),
    TaxDistinct_max = max(
      Average_taxonomic_distinctness,
      na.rm = TRUE
    )
  )

checagem_preditores


# Shared functions --------------------------------------------------------

## Prepare predictors -----------------------------------------------------

preparar_preditores <- function(dados) {
  
  dados |>
    mutate(
      log_distance = log10(Distance),
      log_hii = log10(HII),
      log_taxdist = log10(TaxDistinct),
      log_bio15 = log10(Bio15),
      
      distance_z = as.numeric(scale(log_distance)),
      hii_z = as.numeric(scale(log_hii)),
      taxdist_z = as.numeric(scale(log_taxdist)),
      bio15_z = as.numeric(scale(log_bio15)),
      
      Dataset = droplevels(Dataset),
      Basin_ID = droplevels(Basin_ID)
    )
}


## Extract random-effect variance ----------------------------------------

extrair_variancia <- function(modelo, grupo) {
  
  vc <- VarCorr(modelo)$cond
  
  if (!grupo %in% names(vc)) {
    return(NA_real_)
  }
  
  matriz <- as.matrix(vc[[grupo]])
  
  unname(
    matriz["(Intercept)", "(Intercept)"]
  )
}


## Extract Nakagawa R2 ----------------------------------------------------

extrair_r2 <- function(modelo) {
  
  resultado <- suppressWarnings(
    tryCatch(
      performance::r2_nakagawa(
        modelo,
        tolerance = 1e-8,
        approximation = "delta"
      ),
      error = function(e) NULL
    )
  )
  
  if (is.null(resultado)) {
    return(tibble(
      R2_marginal = NA_real_,
      R2_conditional = NA_real_
    ))
  }
  
  tibble(
    R2_marginal = unname(resultado$R2_marginal),
    R2_conditional = unname(resultado$R2_conditional)
  )
}


## Model summary table ----------------------------------------------------

resumir_modelos <- function(modelos, dados) {
  
  imap_dfr(
    modelos,
    function(modelo, resposta) {
      
      r2 <- extrair_r2(modelo)
      
      grupos <- names(VarCorr(modelo)$cond)
      
      tibble(
        Response = resposta,
        Family = family(modelo)$family,
        
        Random_effect_structure = paste(
          grupos,
          collapse = " + "
        ),
        
        N = nobs(modelo),
        Datasets = n_distinct(dados$Dataset),
        Basins = n_distinct(dados$Basin_ID),
        
        AIC = AIC(modelo),
        
        R2_marginal = r2$R2_marginal,
        R2_conditional = r2$R2_conditional,
        
        Variance_Dataset = extrair_variancia(
          modelo, "Dataset"
        ),
        
        Variance_Basin = extrair_variancia(
          modelo, "Basin_ID"
        )
      )
    }
  ) |>
    mutate(
      across(
        c(
          AIC,
          R2_marginal,
          R2_conditional,
          Variance_Dataset,
          Variance_Basin
        ),
        ~ round(.x, 4)
      )
    )
}


## Model diagnostics ------------------------------------------------------

diagnosticar_modelos <- function(modelos) {
  
  diagnosticos <- imap_dfr(
    modelos,
    function(modelo, resposta) {
      
      tibble(
        Response = resposta,
        Convergence = modelo$fit$convergence,
        Positive_Hessian = modelo$sdr$pdHess,
        
        Singular = isTRUE(
          performance::check_singularity(
            modelo,
            tolerance = 1e-8
          )
        )
      )
    }
  )
  
  set.seed(123)
  
  residuos <- map(
    modelos,
    ~ DHARMa::simulateResiduals(
      fittedModel = .x,
      n = 1000,
      plot = FALSE
    )
  )
  
  testes <- imap_dfr(
    residuos,
    function(residuo, resposta) {
      
      uniformidade <- DHARMa::testUniformity(
        residuo, plot = FALSE
      )
      
      dispersao <- DHARMa::testDispersion(
        residuo, plot = FALSE
      )
      
      outliers <- DHARMa::testOutliers(
        residuo, plot = FALSE
      )
      
      quantis <- tryCatch(
        DHARMa::testQuantiles(
          residuo, plot = FALSE
        ),
        error = function(e) NULL
      )
      
      tibble(
        Response = resposta,
        Uniformity_p = uniformidade$p.value,
        Dispersion_ratio = unname(dispersao$statistic),
        Dispersion_p = dispersao$p.value,
        
        Quantiles_p = if (is.null(quantis)) {
          NA_real_
        } else {
          quantis$p.value
        },
        
        Outliers_p = outliers$p.value
      )
    }
  ) |>
    mutate(
      Residual_problem = if_any(
        c(
          Uniformity_p,
          Dispersion_p,
          Quantiles_p,
          Outliers_p
        ),
        ~ !is.na(.x) & .x < 0.05
      )
    )
  
  diagnosticos |>
    left_join(testes, by = "Response")
}


## Extract coefficients ---------------------------------------------------

extrair_coeficientes <- function(modelos) {
  
  imap_dfr(
    modelos,
    function(modelo, resposta) {
      
      broom.mixed::tidy(
        modelo,
        effects = "fixed",
        component = "cond",
        conf.int = TRUE,
        conf.level = 0.95
      ) |>
        filter(term != "(Intercept)") |>
        mutate(Response = resposta)
    }
  ) |>
    mutate(
      term_label = case_when(
        term == "distance_z" ~
          "River distance",
        
        term == "hii_z" ~
          "Human Influence Index (2020)",
        
        term == "taxdist_z" ~
          "Average taxonomic distinctness",
        
        term == "bio15_z" ~
          "Precipitation seasonality",
        
        term == "distance_z:hii_z" ~
          "River distance × human influence",
        
        term == "distance_z:taxdist_z" ~
          "River distance × taxonomic distinctness",
        
        term == "distance_z:bio15_z" ~
          "River distance × precipitation seasonality",
        
        TRUE ~ term
      ),
      
      significant = !is.na(p.value) & p.value < 0.05,
      
      coefficient_label = if_else(
        significant,
        sprintf("%.2f", estimate),
        NA_character_
      )
    )
}


## Predictor order --------------------------------------------------------

ordem_termos <- c(
  "River distance",
  "Human Influence Index (2020)",
  "Average taxonomic distinctness",
  "Precipitation seasonality",
  "River distance × human influence",
  "River distance × taxonomic distinctness",
  "River distance × precipitation seasonality"
)


## Model control ----------------------------------------------------------

controle_modelo <- glmmTMBControl(
  optCtrl = list(
    iter.max = 10000,
    eval.max = 10000
  )
)


# Figure 2A - Richness ----------------------------------------------------

## Model data -------------------------------------------------------------

dados_riqueza <- env |>
  transmute(
    Dataset = factor(Dataset),
    Basin_ID = factor(Basin_ID),
    
    S = as.numeric(richness),
    
    S_amphibious_n = as.integer(
      S_amphibious_count
    ),
    
    Distance = as.numeric(River_distance_m),
    HII = as.numeric(HII_2020),
    
    TaxDistinct = as.numeric(
      Average_taxonomic_distinctness
    ),
    
    Bio15 = as.numeric(Precipitation_seasonality)
  ) |>
  filter(
    !is.na(Dataset),
    !is.na(Basin_ID),
    is.finite(S),
    S >= 2,
    is.finite(S_amphibious_n),
    is.finite(Distance),
    Distance > 0,
    is.finite(HII),
    HII > 0,
    is.finite(TaxDistinct),
    TaxDistinct > 0,
    is.finite(Bio15),
    Bio15 > 0
  ) |>
  mutate(
    S = as.integer(round(S)),
    S_minus_1 = S - 1L,
    
    S_non_amphibious_n =
      S - S_amphibious_n,
    
    S_amphibious_prop =
      S_amphibious_n / S
  ) |>
  preparar_preditores()

stopifnot(
  nrow(dados_riqueza) > 0,
  all(dados_riqueza$S_minus_1 > 0),
  all(dados_riqueza$S_amphibious_n >= 0),
  all(dados_riqueza$S_non_amphibious_n >= 0),
  all(
    dados_riqueza$S_amphibious_n +
      dados_riqueza$S_non_amphibious_n ==
      dados_riqueza$S
  )
)

resumo_dados <- dados_riqueza |>
  summarise(
    N = n(),
    Datasets = n_distinct(Dataset),
    Basins = n_distinct(Basin_ID),
    
    Zero_amphibious = sum(S_amphibious_n == 0),
    With_amphibious = sum(S_amphibious_n > 0),
    All_amphibious = sum(S_non_amphibious_n == 0)
  )

resumo_dados


## Transformations --------------------------------------------------------

parametros_transformacao <- dados_riqueza |>
  summarise(
    across(
      c(
        log_distance,
        log_hii,
        log_taxdist,
        log_bio15
      ),
      list(mean = mean, sd = sd)
    )
  )

parametros_transformacao


## Collinearity -----------------------------------------------------------

formula_vif <- S_minus_1 ~
  distance_z * (
    hii_z +
      taxdist_z +
      bio15_z
  )

modelo_vif <- lm(
  formula_vif,
  data = dados_riqueza
)

vif_riqueza <- performance::check_collinearity(
  modelo_vif
)

vif_riqueza
plot(vif_riqueza)

correlacoes_preditores <- dados_riqueza |>
  select(
    distance_z,
    hii_z,
    taxdist_z,
    bio15_z
  ) |>
  cor() |>
  round(2)

correlacoes_preditores


## Models -----------------------------------------------------------------

modelo_S_total <- glmmTMB(
  S_minus_1 ~
    distance_z * (
      hii_z +
        taxdist_z +
        bio15_z
    ) +
    (1 | Dataset) +
    (1 | Basin_ID),
  
  data = dados_riqueza,
  
  family = truncated_nbinom2(
    link = "log"
  ),
  
  control = controle_modelo
)

modelo_S_amphibious <- glmmTMB(
  cbind(
    S_amphibious_n,
    S_non_amphibious_n
  ) ~
    distance_z * (
      hii_z +
        taxdist_z +
        bio15_z
    ) +
    (1 | Dataset),
  
  data = dados_riqueza,
  
  family = binomial(
    link = "logit"
  ),
  
  control = controle_modelo
)

modelos_riqueza <- list(
  "Total richness" = modelo_S_total,
  "Amphibious richness (%)" = modelo_S_amphibious
)

summary(modelo_S_total)
summary(modelo_S_amphibious)


## Model summary ----------------------------------------------------------

tabela_modelos <- resumir_modelos(
  modelos_riqueza,
  dados_riqueza
)

tabela_modelos

## Coefficients -----------------------------------------------------------

coeficientes_riqueza <- extrair_coeficientes(
  modelos_riqueza
)

coeficientes_riqueza


# Figure 2B - Abundance ---------------------------------------------------

## Model data -------------------------------------------------------------

dados_abundancia <- env |>
  transmute(
    Dataset = factor(Dataset),
    Basin_ID = factor(Basin_ID),
    Data_Type,
    
    Abundance = as.numeric(abundance),
    
    Abundance_amphibious_n = as.integer(
      abundance_amphibious_count
    ),
    
    Distance = as.numeric(River_distance_m),
    HII = as.numeric(HII_2020),
    
    TaxDistinct = as.numeric(
      Average_taxonomic_distinctness
    ),
    
    Bio15 = as.numeric(Precipitation_seasonality)
  ) |>
  filter(
    Data_Type == "Abundance",
    !is.na(Dataset),
    !is.na(Basin_ID),
    is.finite(Abundance),
    Abundance > 0,
    is.finite(Abundance_amphibious_n),
    is.finite(Distance),
    Distance > 0,
    is.finite(HII),
    HII > 0,
    is.finite(TaxDistinct),
    TaxDistinct > 0,
    is.finite(Bio15),
    Bio15 > 0
  ) |>
  mutate(
    Abundance = as.integer(round(Abundance)),
    
    Abundance_non_amphibious_n =
      Abundance - Abundance_amphibious_n,
    
    Abundance_amphibious_prop =
      Abundance_amphibious_n / Abundance
  ) |>
  preparar_preditores()

stopifnot(
  nrow(dados_abundancia) > 0,
  all(dados_abundancia$Abundance > 0),
  all(dados_abundancia$Abundance_amphibious_n >= 0),
  all(dados_abundancia$Abundance_non_amphibious_n >= 0),
  all(
    dados_abundancia$Abundance_amphibious_n +
      dados_abundancia$Abundance_non_amphibious_n ==
      dados_abundancia$Abundance
  )
)

resumo_dados_abundancia <- dados_abundancia |>
  summarise(
    N = n(),
    Datasets = n_distinct(Dataset),
    Basins = n_distinct(Basin_ID),
    
    Zero_amphibious = sum(
      Abundance_amphibious_n == 0
    ),
    
    With_amphibious = sum(
      Abundance_amphibious_n > 0
    ),
    
    All_amphibious = sum(
      Abundance_non_amphibious_n == 0
    )
  )

resumo_dados_abundancia


## Transformations --------------------------------------------------------

parametros_transformacao_abundancia <- dados_abundancia |>
  summarise(
    across(
      c(
        log_distance,
        log_hii,
        log_taxdist,
        log_bio15
      ),
      list(mean = mean, sd = sd)
    )
  )

parametros_transformacao_abundancia


## Collinearity -----------------------------------------------------------

formula_vif_abundancia <- log10(Abundance) ~
  distance_z * (
    hii_z +
      taxdist_z +
      bio15_z
  )

modelo_vif_abundancia <- lm(
  formula_vif_abundancia,
  data = dados_abundancia
)

vif_abundancia <- performance::check_collinearity(
  modelo_vif_abundancia
)

vif_abundancia
plot(vif_abundancia)

correlacoes_abundancia <- dados_abundancia |>
  select(
    distance_z,
    hii_z,
    taxdist_z,
    bio15_z
  ) |>
  cor() |>
  round(2)

correlacoes_abundancia


## Models -----------------------------------------------------------------

# Total abundance: Dataset + Basin_ID.

modelo_A_total <- glmmTMB(
  Abundance ~
    distance_z * (
      hii_z +
        taxdist_z +
        bio15_z
    ) +
    (1 | Dataset) +
    (1 | Basin_ID),
  
  data = dados_abundancia,
  
  family = truncated_nbinom2(
    link = "log"
  ),
  
  control = controle_modelo
)

# Amphibious abundance: Dataset only.
# Basin_ID was removed because its estimated variance was zero.

modelo_A_amphibious <- glmmTMB(
  cbind(
    Abundance_amphibious_n,
    Abundance_non_amphibious_n
  ) ~
    distance_z * (
      hii_z +
        taxdist_z +
        bio15_z
    ) +
    (1 | Dataset),
  
  data = dados_abundancia,
  
  family = betabinomial(
    link = "logit"
  ),
  
  control = controle_modelo
)

modelos_abundancia <- list(
  "Total abundance" = modelo_A_total,
  "Amphibious abundance (%)" = modelo_A_amphibious
)

summary(modelo_A_total)
summary(modelo_A_amphibious)


## Model summary ----------------------------------------------------------

tabela_modelos_abundancia <- resumir_modelos(
  modelos_abundancia,
  dados_abundancia
)

tabela_modelos_abundancia

## Coefficients -----------------------------------------------------------

coeficientes_abundancia <- extrair_coeficientes(
  modelos_abundancia
)

coeficientes_abundancia


# Figure 2 - Plot preparation --------------------------------------------

## Shared plotting function -----------------------------------------------

preparar_coeficientes_plot <- function(
    coeficientes,
    ordem_respostas) {
  
  deslocamentos <- setNames(
    c(0.13, -0.13),
    ordem_respostas
  )
  
  dados <- coeficientes |>
    mutate(
      Response = factor(
        Response,
        levels = ordem_respostas
      ),
      
      term_y = length(ordem_termos) -
        match(term_label, ordem_termos) + 1,
      
      y_plot = term_y +
        unname(
          deslocamentos[as.character(Response)]
        )
    )
  
  stopifnot(!anyNA(dados$term_y))
  
  amplitude_x <- diff(
    range(
      c(dados$conf.low, dados$conf.high),
      na.rm = TRUE
    )
  )
  
  rotulos <- dados |>
    filter(significant) |>
    mutate(
      label_x = estimate + 0.045 * amplitude_x,
      
      label_y = y_plot +
        if_else(
          Response == ordem_respostas[1],
          0.18,
          -0.18
        ),
      
      line_xend = label_x - 0.01 * amplitude_x
    )
  
  list(
    dados = dados,
    rotulos = rotulos
  )
}


## Prepare panel A --------------------------------------------------------

ordem_respostas <- c(
  "Total richness",
  "Amphibious richness (%)"
)

preparacao_A <- preparar_coeficientes_plot(
  coeficientes_riqueza,
  ordem_respostas
)

coeficientes_plot <- preparacao_A$dados
dados_rotulos <- preparacao_A$rotulos


## Prepare panel B --------------------------------------------------------

ordem_respostas_abundancia <- c(
  "Total abundance",
  "Amphibious abundance (%)"
)

preparacao_B <- preparar_coeficientes_plot(
  coeficientes_abundancia,
  ordem_respostas_abundancia
)

coeficientes_plot_abundancia <- preparacao_B$dados
dados_rotulos_abundancia <- preparacao_B$rotulos


## Colors and symbols -----------------------------------------------------

cores_respostas <- c(
  "Total richness" = "#77A88D",
  "Amphibious richness (%)" = "#B7866B"
)

simbolos_respostas <- c(
  "Total richness" = 21,
  "Amphibious richness (%)" = 22
)

cores_respostas_abundancia <- c(
  "Total abundance" = "#77A88D",
  "Amphibious abundance (%)" = "#B7866B"
)

simbolos_respostas_abundancia <- c(
  "Total abundance" = 21,
  "Amphibious abundance (%)" = 22
)


## Common x-axis limits ---------------------------------------------------

# Reserve additional space on the right for both legends.

valores_x <- c(
  coeficientes_plot$conf.low,
  coeficientes_plot$conf.high,
  dados_rotulos$label_x,
  coeficientes_plot_abundancia$conf.low,
  coeficientes_plot_abundancia$conf.high,
  dados_rotulos_abundancia$label_x,
  0
)

amplitude_total <- diff(
  range(valores_x, na.rm = TRUE)
)

limites_x_fig2 <- c(
  min(valores_x, na.rm = TRUE) -
    0.06 * amplitude_total,
  
  max(valores_x, na.rm = TRUE) +
    0.65 * amplitude_total
)

limites_x_fig2


# Figure 2 - Shared forest plot function ---------------------------------

criar_forest_plot <- function(
    dados,
    rotulos,
    cores,
    simbolos,
    rotulos_legenda,
    tag_painel) {
  
  ggplot(dados) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey55",
      linewidth = 0.7
    ) +
    
    geom_segment(
      aes(
        x = conf.low,
        xend = conf.high,
        y = y_plot,
        yend = y_plot,
        color = Response
      ),
      linewidth = 1.3,
      lineend = "round",
      alpha = 0.70,
      show.legend = FALSE
    ) +
    
    geom_segment(
      data = rotulos,
      aes(
        x = estimate,
        xend = estimate,
        y = y_plot,
        yend = label_y,
        color = Response
      ),
      linewidth = 0.9,
      lineend = "round",
      alpha = 0.9,
      show.legend = FALSE
    ) +
    
    geom_segment(
      data = rotulos,
      aes(
        x = estimate,
        xend = line_xend,
        y = label_y,
        yend = label_y,
        color = Response
      ),
      linewidth = 0.9,
      lineend = "round",
      alpha = 0.9,
      show.legend = FALSE
    ) +
    
    geom_point(
      aes(
        x = estimate,
        y = y_plot,
        fill = Response,
        shape = Response
      ),
      size = 5.6,
      color = "black",
      stroke = 0.65,
      alpha = 0.98,
      show.legend = TRUE
    ) +
    
    geom_text(
      data = rotulos,
      aes(
        x = label_x,
        y = label_y,
        label = coefficient_label,
        color = Response
      ),
      family = "serif",
      fontface = "bold",
      size = 6.2,
      hjust = 0,
      vjust = 0.5,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    
    scale_x_continuous(
      expand = expansion(mult = c(0, 0))
    ) +
    
    scale_y_continuous(
      breaks = seq_along(ordem_termos),
      labels = rev(ordem_termos),
      
      limits = c(
        0.45,
        length(ordem_termos) + 0.55
      ),
      
      expand = expansion(
        mult = c(0.01, 0.01)
      )
    ) +
    
    scale_color_manual(
      values = cores,
      guide = "none"
    ) +
    
    scale_fill_manual(
      values = cores,
      labels = rotulos_legenda,
      name = NULL
    ) +
    
    scale_shape_manual(
      values = simbolos,
      labels = rotulos_legenda,
      name = NULL
    ) +
    
    guides(
      shape = "none",
      
      fill = guide_legend(
        override.aes = list(
          shape = unname(simbolos),
          size = 5.5,
          color = "black",
          alpha = 1
        )
      )
    ) +
    
    coord_cartesian(
      xlim = limites_x_fig2,
      clip = "off"
    ) +
    
    labs(
      x = "Standardized coefficient (95% CI)",
      y = NULL,
      tag = tag_painel
    ) +
    
    theme_classic(
      base_size = 22,
      base_family = "serif"
    ) +
    
    theme(
      axis.title.x = element_text(
        size = 21,
        margin = margin(t = 12)
      ),
      
      axis.text.x = element_text(
        size = 18,
        color = "black"
      ),
      
      axis.text.y = element_text(
        face = "bold",
        size = 16,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.6
      ),
      
      plot.tag = element_text(
        family = "serif",
        face = "bold",
        size = 25
      ),
      
      plot.tag.position = c(
        0.015,
        0.985
      ),
      
      legend.position = "inside",
      
      legend.position.inside = c(
        1.05,
        0.97
      ),
      
      legend.justification = c(
        1,
        1
      ),
      
      legend.text = element_text(
        family = "serif",
        size = 16
      ),
      
      legend.key.height = grid::unit(
        0.85,
        "cm"
      ),
      
      legend.key.width = grid::unit(
        0.85,
        "cm"
      ),
      
      legend.spacing.y = grid::unit(
        0.05,
        "cm"
      ),
      
      legend.background = element_rect(
        fill = scales::alpha("white", 0.95),
        color = "grey75",
        linewidth = 0.3
      ),
      
      legend.margin = margin(
        6,
        8,
        6,
        8
      ),
      
      plot.margin = margin(
        18,
        20,
        16,
        16
      )
    )
}


# Figure 2A - Final plot --------------------------------------------------

rotulos_legenda_A <- c(
  "Total richness" = "Community richness",
  "Amphibious richness (%)" = "Amphibious richness (%)"
)

fig2_A <- criar_forest_plot(
  dados = coeficientes_plot,
  rotulos = dados_rotulos,
  cores = cores_respostas,
  simbolos = simbolos_respostas,
  rotulos_legenda = rotulos_legenda_A,
  tag_painel = "A"
)

fig2_A


# Figure 2B - Final plot --------------------------------------------------

rotulos_legenda_B <- c(
  "Total abundance" = "Community abundance",
  "Amphibious abundance (%)" = "Amphibious abundance (%)"
)

fig2_B <- criar_forest_plot(
  dados = coeficientes_plot_abundancia,
  rotulos = dados_rotulos_abundancia,
  cores = cores_respostas_abundancia,
  simbolos = simbolos_respostas_abundancia,
  rotulos_legenda = rotulos_legenda_B,
  tag_painel = "B"
)

fig2_B


# Figure 2 - Complete -----------------------------------------------------

figura_2 <- (
  fig2_A | fig2_B
) +
  patchwork::plot_layout(
    widths = c(1, 1)
  )

figura_2


## Save -------------------------------------------------------------------

ggsave(
  filename = "Figure_2.jpg",
  plot = figura_2,
  width = 20,
  height = 8,
  units = "in",
  dpi = 600,
  bg = "white"
)

# Figure 3 ---------------------------------------------------------------

## Interaction: river distance x taxonomic distinctness -------------------

# Predict total community abundance using the updated model from Figure 2.
# Other predictors are held at their standardized mean (z = 0).
# Predictions exclude random effects (population-level predictions).


## Prediction function ----------------------------------------------------

predizer_interacao <- function(modelo, dados) {
  
  niveis_z <- c(-1, 0, 1)
  
  grade <- tidyr::crossing(
    distance_z = seq(
      min(dados$distance_z),
      max(dados$distance_z),
      length.out = 300
    ),
    taxdist_z = niveis_z
  ) |>
    mutate(
      hii_z = 0,
      bio15_z = 0,
      
      Dataset = factor(
        levels(dados$Dataset)[1],
        levels = levels(dados$Dataset)
      ),
      
      Basin_ID = factor(
        levels(dados$Basin_ID)[1],
        levels = levels(dados$Basin_ID)
      )
    )
  
  # Back-transform standardized log10 river distance.
  
  media_log_distance <- mean(
    dados$log_distance
  )
  
  sd_log_distance <- sd(
    dados$log_distance
  )
  
  grade <- grade |>
    mutate(
      Distance = 10^(
        media_log_distance +
          distance_z * sd_log_distance
      ),
      
      Moderator_level = factor(
        taxdist_z,
        levels = c(-1, 0, 1),
        labels = c(
          "Low (−1 SD)",
          "Average",
          "High (+1 SD)"
        )
      )
    )
  
  # Population-level predictions on the link scale.
  
  predicao <- predict(
    modelo,
    newdata = grade,
    type = "link",
    se.fit = TRUE,
    re.form = NA,
    allow.new.levels = TRUE
  )
  
  inversa_link <- family(modelo)$linkinv
  
  grade |>
    mutate(
      estimate = inversa_link(
        as.numeric(predicao$fit)
      ),
      
      CI_low = inversa_link(
        as.numeric(
          predicao$fit -
            1.96 * predicao$se.fit
        )
      ),
      
      CI_high = inversa_link(
        as.numeric(
          predicao$fit +
            1.96 * predicao$se.fit
        )
      )
    )
}


## Predictions ------------------------------------------------------------

predicao_fig3 <- predizer_interacao(
  modelo = modelo_A_total,
  dados = dados_abundancia
)

head(predicao_fig3)


## Raw moderator values ---------------------------------------------------

# Convert -1 SD, mean and +1 SD to the original scale
# of average taxonomic distinctness.

niveis_moderadores <- tibble(
  Moderator = "Average taxonomic distinctness",
  
  Level = c(
    "Low (−1 SD)",
    "Average",
    "High (+1 SD)"
  ),
  
  z = c(-1, 0, 1)
) |>
  mutate(
    Raw_value = 10^(
      mean(dados_abundancia$log_taxdist) +
        z * sd(dados_abundancia$log_taxdist)
    )
  )

niveis_moderadores


## Colors and line types --------------------------------------------------

cores_interacao <- c(
  "Low (−1 SD)" = "#BFA9D4",
  "Average" = "#7899B5",
  "High (+1 SD)" = "#6B625E"
)

tipos_linha <- c(
  "Low (−1 SD)" = "dotted",
  "Average" = "dashed",
  "High (+1 SD)" = "solid"
)


## Moderator label position -----------------------------------------------

intervalo_log_x <- range(
  log10(predicao_fig3$Distance),
  na.rm = TRUE
)

posicao_moderador_x <- 10^(
  intervalo_log_x[1] +
    0.72 * diff(intervalo_log_x)
)


## Figure 3 ---------------------------------------------------------------

Figure_3 <- ggplot(
  predicao_fig3,
  aes(
    x = Distance,
    y = estimate,
    color = Moderator_level,
    fill = Moderator_level,
    linetype = Moderator_level,
    group = Moderator_level
  )
) +
  
  # 95% confidence intervals.
  
  geom_ribbon(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    color = NA,
    linetype = "blank",
    alpha = 0.14
  ) +
  
  # Predicted abundance.
  
  geom_line(
    linewidth = 1.7,
    alpha = 0.95
  ) +
  
  # Moderator annotation.
  
  annotate(
    geom = "label",
    x = posicao_moderador_x,
    y = Inf,
    label = "Moderator:\nAverage taxonomic distinctness",
    hjust = 0.5,
    vjust = 1.15,
    family = "serif",
    fontface = "bold",
    size = 5.2,
    lineheight = 1.05,
    color = "grey20",
    fill = scales::alpha("white", 0.92),
    linewidth = 0.35,
    label.padding = grid::unit(
      0.35,
      "lines"
    ),
    label.r = grid::unit(
      0.12,
      "lines"
    )
  ) +
  
  # Log10 distance axis.
  
  scale_x_log10(
    labels = scales::label_number(
      big.mark = ","
    )
  ) +
  
  # Colors.
  
  scale_color_manual(
    values = cores_interacao,
    breaks = names(cores_interacao),
    name = "Moderator level",
    drop = FALSE
  ) +
  
  scale_fill_manual(
    values = cores_interacao,
    guide = "none",
    drop = FALSE
  ) +
  
  scale_linetype_manual(
    values = tipos_linha,
    guide = "none",
    drop = FALSE
  ) +
  
  guides(
    color = guide_legend(
      nrow = 1,
      byrow = TRUE,
      override.aes = list(
        linewidth = 2,
        linetype = unname(tipos_linha),
        alpha = 1
      )
    )
  ) +
  
  labs(
    x = "Distance to nearest river or stream (m)",
    y = "Estimated community abundance"
  ) +
  
  theme_classic(
    base_size = 20,
    base_family = "serif"
  ) +
  
  theme(
    axis.title.x = element_text(
      size = 22,
      margin = margin(t = 14)
    ),
    
    axis.title.y = element_text(
      size = 22,
      margin = margin(r = 14)
    ),
    
    axis.text.x = element_text(
      size = 18,
      color = "black"
    ),
    
    axis.text.y = element_text(
      size = 18,
      color = "black"
    ),
    
    axis.ticks = element_line(
      linewidth = 0.7
    ),
    
    axis.line = element_line(
      linewidth = 0.7
    ),
    
    legend.position = "bottom",
    legend.direction = "horizontal",
    
    legend.title = element_text(
      family = "serif",
      face = "bold",
      size = 18
    ),
    
    legend.text = element_text(
      family = "serif",
      size = 17
    ),
    
    legend.key.width = grid::unit(
      1.8,
      "cm"
    ),
    
    legend.key.height = grid::unit(
      0.65,
      "cm"
    ),
    
    legend.spacing.x = grid::unit(
      0.3,
      "cm"
    ),
    
    legend.margin = margin(
      t = 10
    ),
    
    plot.margin = margin(
      20,
      25,
      20,
      20
    )
  )


## Display ---------------------------------------------------------------

Figure_3


## Save ------------------------------------------------------------------

ggsave(
  filename = "Figure_3.jpg",
  plot = Figure_3,
  width = 10,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# Figure 4 ---------------------------------------------------------------

# Amphibious classification ----------------------------------------------

stopifnot(
  nrow(env) == nrow(assemblages),
  nrow(species) == ncol(assemblages),
  identical(species$Species, colnames(assemblages)),
  all(species$Amphibious_lifestype %in% c("Y", "N"))
)

resumo_amphibious_volume <- species |>
  summarise(
    Total_species = n(),
    Amphibious_species = sum(Amphibious_lifestype == "Y"),
    Amphibious_percentage =
      100 * mean(Amphibious_lifestype == "Y")
  )

resumo_amphibious_volume


# Recalculate response variables -----------------------------------------

matriz_assemblages <- as.matrix(assemblages)

stopifnot(
  is.numeric(matriz_assemblages),
  !anyNA(matriz_assemblages),
  all(matriz_assemblages >= 0),
  all(matriz_assemblages == round(matriz_assemblages))
)

matriz_PA <- (matriz_assemblages > 0) * 1L

colunas_amphibious <-
  species$Amphibious_lifestype == "Y"

S_total <- rowSums(matriz_PA)

S_amphibious_count <- rowSums(
  matriz_PA[, colunas_amphibious, drop = FALSE]
)

abundance_total <- rowSums(matriz_assemblages)

abundance_amphibious_count <- rowSums(
  matriz_assemblages[, colunas_amphibious, drop = FALSE]
)

stopifnot(
  all(S_total > 0),
  all(S_amphibious_count >= 0),
  all(S_amphibious_count <= S_total),
  all(abundance_amphibious_count >= 0),
  all(abundance_amphibious_count <= abundance_total)
)

env <- env |>
  mutate(
    richness = as.numeric(S_total),
    
    S_amphibious_count =
      as.integer(S_amphibious_count),
    
    S_amphibious =
      100 * S_amphibious_count / richness,
    
    abundance = if_else(
      Data_Type == "Abundance",
      as.numeric(abundance_total),
      NA_real_
    ),
    
    abundance_amphibious_count =
      as.integer(abundance_amphibious_count),
    
    abundance_amphibious = if_else(
      Data_Type == "Abundance" & abundance_total > 0,
      100 * abundance_amphibious_count / abundance_total,
      NA_real_
    )
  )


# Descriptive response tables --------------------------------------------

checagem_riqueza_volume <- tibble(
  Metric = c(
    "Total richness",
    "Amphibious richness"
  ),
  
  Minimum = c(
    min(S_total),
    min(S_amphibious_count)
  ),
  
  Maximum = c(
    max(S_total),
    max(S_amphibious_count)
  ),
  
  Mean = c(
    mean(S_total),
    mean(S_amphibious_count)
  )
)

checagem_riqueza_volume

checagem_abundancia_volume <- env |>
  summarise(
    N_abundance = sum(
      Data_Type == "Abundance" &
        is.finite(abundance)
    ),
    
    Abundance_min = min(
      abundance,
      na.rm = TRUE
    ),
    
    Abundance_max = max(
      abundance,
      na.rm = TRUE
    ),
    
    Amphibious_percent_min = min(
      abundance_amphibious,
      na.rm = TRUE
    ),
    
    Amphibious_percent_max = max(
      abundance_amphibious,
      na.rm = TRUE
    )
  )

checagem_abundancia_volume


# Check updated predictors -----------------------------------------------

checagem_preditores_volume <- env |>
  summarise(
    N_sites = n(),
    
    HII_missing = sum(is.na(HII_2020)),
    
    HII_zero = sum(
      HII_2020 == 0,
      na.rm = TRUE
    ),
    
    TaxDistinct_missing = sum(
      is.na(Average_taxonomic_distinctness)
    ),
    
    TaxDistinct_min = min(
      Average_taxonomic_distinctness,
      na.rm = TRUE
    ),
    
    TaxDistinct_max = max(
      Average_taxonomic_distinctness,
      na.rm = TRUE
    ),
    
    Volume_missing = sum(is.na(Volume_m3)),
    
    Volume_zero = sum(
      Volume_m3 == 0,
      na.rm = TRUE
    )
  )

checagem_preditores_volume


# Shared predictor preparation -------------------------------------------

preparar_preditores_volume <- function(dados) {
  
  dados |>
    mutate(
      log_volume = log10(Volume),
      log_distance = log10(Distance),
      log_hii = log10(HII),
      log_taxdist = log10(TaxDistinct),
      log_bio15 = log10(Bio15),
      
      volume_z = as.numeric(scale(log_volume)),
      distance_z = as.numeric(scale(log_distance)),
      hii_z = as.numeric(scale(log_hii)),
      taxdist_z = as.numeric(scale(log_taxdist)),
      bio15_z = as.numeric(scale(log_bio15)),
      
      Dataset = droplevels(Dataset)
    )
}


# Figure 4A - Richness data -----------------------------------------------

dados_riqueza_volume <- env |>
  transmute(
    Dataset = factor(Dataset),
    
    S = as.numeric(richness),
    
    S_amphibious_n =
      as.integer(S_amphibious_count),
    
    Volume = as.numeric(Volume_m3),
    
    Distance = as.numeric(River_distance_m),
    
    HII = as.numeric(HII_2020),
    
    TaxDistinct = as.numeric(
      Average_taxonomic_distinctness
    ),
    
    Bio15 = as.numeric(
      Precipitation_seasonality
    )
  ) |>
  filter(
    !is.na(Dataset),
    is.finite(S),
    S >= 2,
    is.finite(S_amphibious_n),
    is.finite(Volume),
    Volume > 0,
    is.finite(Distance),
    Distance > 0,
    is.finite(HII),
    HII > 0,
    is.finite(TaxDistinct),
    TaxDistinct > 0,
    is.finite(Bio15),
    Bio15 > 0
  ) |>
  mutate(
    S = as.integer(round(S)),
    
    S_minus_1 = S - 1L,
    
    S_non_amphibious_n =
      S - S_amphibious_n,
    
    S_amphibious_prop =
      S_amphibious_n / S
  ) |>
  preparar_preditores_volume()

stopifnot(
  nrow(dados_riqueza_volume) > 0,
  all(dados_riqueza_volume$S_minus_1 > 0),
  all(dados_riqueza_volume$S_amphibious_n >= 0),
  all(dados_riqueza_volume$S_non_amphibious_n >= 0),
  all(
    dados_riqueza_volume$S_amphibious_n +
      dados_riqueza_volume$S_non_amphibious_n ==
      dados_riqueza_volume$S
  )
)

resumo_dados_riqueza_volume <- dados_riqueza_volume |>
  summarise(
    N = n(),
    Datasets = n_distinct(Dataset),
    
    Zero_amphibious = sum(
      S_amphibious_n == 0
    ),
    
    With_amphibious = sum(
      S_amphibious_n > 0
    ),
    
    All_amphibious = sum(
      S_non_amphibious_n == 0
    )
  )

resumo_dados_riqueza_volume


# Figure 4B - Abundance data ----------------------------------------------

dados_abundancia_volume <- env |>
  transmute(
    Dataset = factor(Dataset),
    
    Data_Type = Data_Type,
    
    Abundance = as.numeric(abundance),
    
    Abundance_amphibious_n =
      as.integer(abundance_amphibious_count),
    
    Volume = as.numeric(Volume_m3),
    
    Distance = as.numeric(River_distance_m),
    
    HII = as.numeric(HII_2020),
    
    TaxDistinct = as.numeric(
      Average_taxonomic_distinctness
    ),
    
    Bio15 = as.numeric(
      Precipitation_seasonality
    )
  ) |>
  filter(
    Data_Type == "Abundance",
    !is.na(Dataset),
    is.finite(Abundance),
    Abundance > 0,
    is.finite(Abundance_amphibious_n),
    is.finite(Volume),
    Volume > 0,
    is.finite(Distance),
    Distance > 0,
    is.finite(HII),
    HII > 0,
    is.finite(TaxDistinct),
    TaxDistinct > 0,
    is.finite(Bio15),
    Bio15 > 0
  ) |>
  mutate(
    Abundance = as.integer(round(Abundance)),
    
    Abundance_non_amphibious_n =
      Abundance - Abundance_amphibious_n,
    
    Abundance_amphibious_prop =
      Abundance_amphibious_n / Abundance
  ) |>
  preparar_preditores_volume()

stopifnot(
  nrow(dados_abundancia_volume) > 0,
  all(dados_abundancia_volume$Abundance > 0),
  all(dados_abundancia_volume$Abundance_amphibious_n >= 0),
  all(dados_abundancia_volume$Abundance_non_amphibious_n >= 0),
  all(
    dados_abundancia_volume$Abundance_amphibious_n +
      dados_abundancia_volume$Abundance_non_amphibious_n ==
      dados_abundancia_volume$Abundance
  )
)

resumo_dados_abundancia_volume <- dados_abundancia_volume |>
  summarise(
    N = n(),
    Datasets = n_distinct(Dataset),
    
    Zero_amphibious = sum(
      Abundance_amphibious_n == 0
    ),
    
    With_amphibious = sum(
      Abundance_amphibious_n > 0
    ),
    
    All_amphibious = sum(
      Abundance_non_amphibious_n == 0
    )
  )

resumo_dados_abundancia_volume


# Transformation parameters ----------------------------------------------

parametros_riqueza_volume <- dados_riqueza_volume |>
  summarise(
    across(
      c(
        log_volume,
        log_distance,
        log_hii,
        log_taxdist,
        log_bio15
      ),
      list(
        mean = mean,
        sd = sd
      )
    )
  )

parametros_riqueza_volume

parametros_abundancia_volume <- dados_abundancia_volume |>
  summarise(
    across(
      c(
        log_volume,
        log_distance,
        log_hii,
        log_taxdist,
        log_bio15
      ),
      list(
        mean = mean,
        sd = sd
      )
    )
  )

parametros_abundancia_volume


# Collinearity ------------------------------------------------------------

## Richness ---------------------------------------------------------------

formula_vif_riqueza_volume <-
  S_minus_1 ~
  volume_z * (
    distance_z +
      hii_z +
      taxdist_z +
      bio15_z
  )

modelo_vif_riqueza_volume <- lm(
  formula_vif_riqueza_volume,
  data = dados_riqueza_volume
)

vif_riqueza_volume <- performance::check_collinearity(
  modelo_vif_riqueza_volume
)

vif_riqueza_volume
plot(vif_riqueza_volume)

correlacoes_riqueza_volume <- dados_riqueza_volume |>
  select(
    volume_z,
    distance_z,
    hii_z,
    taxdist_z,
    bio15_z
  ) |>
  cor() |>
  round(2)

correlacoes_riqueza_volume


## Abundance --------------------------------------------------------------

formula_vif_abundancia_volume <-
  log10(Abundance) ~
  volume_z * (
    distance_z +
      hii_z +
      taxdist_z +
      bio15_z
  )

modelo_vif_abundancia_volume <- lm(
  formula_vif_abundancia_volume,
  data = dados_abundancia_volume
)

vif_abundancia_volume <- performance::check_collinearity(
  modelo_vif_abundancia_volume
)

vif_abundancia_volume
plot(vif_abundancia_volume)

correlacoes_abundancia_volume <- dados_abundancia_volume |>
  select(
    volume_z,
    distance_z,
    hii_z,
    taxdist_z,
    bio15_z
  ) |>
  cor() |>
  round(2)

correlacoes_abundancia_volume


# Models ------------------------------------------------------------------

controle_modelo_volume <- glmmTMBControl(
  optCtrl = list(
    iter.max = 10000,
    eval.max = 10000
  )
)


## Total richness ---------------------------------------------------------

modelo_S_total_volume <- glmmTMB(
  S_minus_1 ~
    volume_z * (
      distance_z +
        hii_z +
        taxdist_z +
        bio15_z
    ) +
    (1 | Dataset),
  
  data = dados_riqueza_volume,
  
  family = truncated_nbinom2(
    link = "log"
  ),
  
  control = controle_modelo_volume
)


## Amphibious richness proportion ----------------------------------------

modelo_S_amphibious_volume <- glmmTMB(
  cbind(
    S_amphibious_n,
    S_non_amphibious_n
  ) ~
    volume_z * (
      distance_z +
        hii_z +
        taxdist_z +
        bio15_z
    ) +
    (1 | Dataset),
  
  data = dados_riqueza_volume,
  
  family = binomial(
    link = "logit"
  ),
  
  control = controle_modelo_volume
)


## Total abundance --------------------------------------------------------

modelo_A_total_volume <- glmmTMB(
  Abundance ~
    volume_z * (
      distance_z +
        hii_z +
        taxdist_z +
        bio15_z
    ) +
    (1 | Dataset),
  
  data = dados_abundancia_volume,
  
  family = truncated_nbinom2(
    link = "log"
  ),
  
  control = controle_modelo_volume
)


## Amphibious abundance proportion ---------------------------------------

modelo_A_amphibious_volume <- glmmTMB(
  cbind(
    Abundance_amphibious_n,
    Abundance_non_amphibious_n
  ) ~
    volume_z * (
      distance_z +
        hii_z +
        taxdist_z +
        bio15_z
    ) +
    (1 | Dataset),
  
  data = dados_abundancia_volume,
  
  family = betabinomial(
    link = "logit"
  ),
  
  control = controle_modelo_volume
)


## Model lists ------------------------------------------------------------

modelos_riqueza_volume <- list(
  "Total richness" = modelo_S_total_volume,
  "Amphibious richness (%)" = modelo_S_amphibious_volume
)

modelos_abundancia_volume <- list(
  "Total abundance" = modelo_A_total_volume,
  "Amphibious abundance (%)" = modelo_A_amphibious_volume
)


## Model summaries --------------------------------------------------------

summary(modelo_S_total_volume)
summary(modelo_S_amphibious_volume)
summary(modelo_A_total_volume)
summary(modelo_A_amphibious_volume)


# Model summary functions -------------------------------------------------

## Random-effect variance -------------------------------------------------

extrair_variancia_volume <- function(
    modelo,
    grupo = "Dataset") {
  
  vc <- VarCorr(modelo)$cond
  
  if (!grupo %in% names(vc)) {
    return(NA_real_)
  }
  
  matriz <- as.matrix(vc[[grupo]])
  
  if (!"(Intercept)" %in% rownames(matriz)) {
    return(NA_real_)
  }
  
  unname(
    matriz["(Intercept)", "(Intercept)"]
  )
}


## Nakagawa R2 ------------------------------------------------------------

extrair_r2_volume <- function(modelo) {
  
  resultado <- suppressWarnings(
    tryCatch(
      performance::r2_nakagawa(
        modelo,
        tolerance = 1e-8,
        approximation = "delta"
      ),
      error = function(e) NULL
    )
  )
  
  if (
    is.null(resultado) ||
    is.null(resultado$R2_marginal)
  ) {
    return(
      tibble(
        R2_marginal = NA_real_,
        R2_conditional = NA_real_
      )
    )
  }
  
  tibble(
    R2_marginal = unname(
      resultado$R2_marginal
    ),
    
    R2_conditional = if (
      is.null(resultado$R2_conditional)
    ) {
      NA_real_
    } else {
      unname(resultado$R2_conditional)
    }
  )
}


## Model summary table ----------------------------------------------------

criar_tabela_modelos_volume <- function(
    modelos,
    dados) {
  
  purrr::imap_dfr(
    modelos,
    function(modelo, resposta) {
      
      r2 <- extrair_r2_volume(modelo)
      
      tibble(
        Response = resposta,
        Family = family(modelo)$family,
        
        Random_effect_structure =
          "Random intercept: Dataset",
        
        N = nobs(modelo),
        
        Datasets = n_distinct(dados$Dataset),
        
        AIC = AIC(modelo),
        
        R2_marginal = r2$R2_marginal,
        
        R2_conditional = r2$R2_conditional,
        
        Variance_Dataset = extrair_variancia_volume(
          modelo,
          "Dataset"
        )
      )
    }
  ) |>
    mutate(
      across(
        c(
          AIC,
          R2_marginal,
          R2_conditional,
          Variance_Dataset
        ),
        ~ round(.x, 4)
      )
    )
}


# Model summary tables ----------------------------------------------------

tabela_modelos_riqueza_volume <-
  criar_tabela_modelos_volume(
    modelos_riqueza_volume,
    dados_riqueza_volume
  )

tabela_modelos_abundancia_volume <-
  criar_tabela_modelos_volume(
    modelos_abundancia_volume,
    dados_abundancia_volume
  )

tabela_modelos_riqueza_volume
tabela_modelos_abundancia_volume

tabela_modelos_figure4 <- bind_rows(
  tabela_modelos_riqueza_volume,
  tabela_modelos_abundancia_volume
)

tabela_modelos_figure4


# Model diagnostics -------------------------------------------------------

criar_diagnosticos_volume <- function(modelos) {
  
  diagnosticos <- purrr::imap_dfr(
    modelos,
    function(modelo, resposta) {
      
      tibble(
        Response = resposta,
        
        Convergence = modelo$fit$convergence,
        
        Positive_Hessian = modelo$sdr$pdHess,
        
        Singular = isTRUE(
          performance::check_singularity(
            modelo,
            tolerance = 1e-8
          )
        )
      )
    }
  )
  
  set.seed(123)
  
  residuos <- purrr::map(
    modelos,
    ~ DHARMa::simulateResiduals(
      fittedModel = .x,
      n = 1000,
      plot = FALSE
    )
  )
  
  testes <- purrr::imap_dfr(
    residuos,
    function(residuo, resposta) {
      
      uniformidade <- DHARMa::testUniformity(
        residuo,
        plot = FALSE
      )
      
      dispersao <- DHARMa::testDispersion(
        residuo,
        plot = FALSE
      )
      
      outliers <- DHARMa::testOutliers(
        residuo,
        plot = FALSE
      )
      
      quantis <- tryCatch(
        DHARMa::testQuantiles(
          residuo,
          plot = FALSE
        ),
        error = function(e) NULL
      )
      
      tibble(
        Response = resposta,
        
        Uniformity_p = uniformidade$p.value,
        
        Dispersion_ratio = unname(
          dispersao$statistic
        ),
        
        Dispersion_p = dispersao$p.value,
        
        Quantiles_p = if (is.null(quantis)) {
          NA_real_
        } else {
          quantis$p.value
        },
        
        Outliers_p = outliers$p.value
      )
    }
  ) |>
    mutate(
      Residual_problem = if_any(
        c(
          Uniformity_p,
          Dispersion_p,
          Quantiles_p,
          Outliers_p
        ),
        ~ !is.na(.x) & .x < 0.05
      )
    )
  
  list(
    table = left_join(
      diagnosticos,
      testes,
      by = "Response"
    ),
    
    residuals = residuos
  )
}


# Coefficient functions ---------------------------------------------------

## Predictor order --------------------------------------------------------

ordem_termos_volume <- c(
  "Pool volume",
  "River distance",
  "Human Influence Index (2020)",
  "Average taxonomic distinctness",
  "Precipitation seasonality",
  "Pool volume × river distance",
  "Pool volume × human influence",
  "Pool volume × taxonomic distinctness",
  "Pool volume × precipitation seasonality"
)


## Extract coefficients ---------------------------------------------------

extrair_coeficientes_volume <- function(modelos) {
  
  purrr::imap_dfr(
    modelos,
    function(modelo, resposta) {
      
      broom.mixed::tidy(
        modelo,
        effects = "fixed",
        component = "cond",
        conf.int = TRUE,
        conf.level = 0.95
      ) |>
        filter(term != "(Intercept)") |>
        mutate(Response = resposta)
    }
  ) |>
    mutate(
      term_label = case_when(
        term == "volume_z" ~
          "Pool volume",
        
        term == "distance_z" ~
          "River distance",
        
        term == "hii_z" ~
          "Human Influence Index (2020)",
        
        term == "taxdist_z" ~
          "Average taxonomic distinctness",
        
        term == "bio15_z" ~
          "Precipitation seasonality",
        
        term %in% c(
          "volume_z:distance_z",
          "distance_z:volume_z"
        ) ~
          "Pool volume × river distance",
        
        term %in% c(
          "volume_z:hii_z",
          "hii_z:volume_z"
        ) ~
          "Pool volume × human influence",
        
        term %in% c(
          "volume_z:taxdist_z",
          "taxdist_z:volume_z"
        ) ~
          "Pool volume × taxonomic distinctness",
        
        term %in% c(
          "volume_z:bio15_z",
          "bio15_z:volume_z"
        ) ~
          "Pool volume × precipitation seasonality",
        
        TRUE ~ term
      ),
      
      significant =
        !is.na(p.value) & p.value < 0.05,
      
      coefficient_label = if_else(
        significant,
        sprintf("%.2f", estimate),
        NA_character_
      )
    )
}


## Prepare coefficient positions -----------------------------------------

preparar_coeficientes_volume <- function(
    coeficientes,
    ordem_respostas) {
  
  deslocamentos <- setNames(
    c(0.13, -0.13),
    ordem_respostas
  )
  
  coeficientes_plot <- coeficientes |>
    mutate(
      Response = factor(
        Response,
        levels = ordem_respostas
      ),
      
      term_y = length(ordem_termos_volume) -
        match(term_label, ordem_termos_volume) + 1,
      
      y_plot = term_y +
        unname(
          deslocamentos[as.character(Response)]
        )
    )
  
  stopifnot(!anyNA(coeficientes_plot$term_y))
  
  amplitude_x <- diff(
    range(
      c(
        coeficientes_plot$conf.low,
        coeficientes_plot$conf.high
      ),
      na.rm = TRUE
    )
  )
  
  dados_rotulos <- coeficientes_plot |>
    filter(significant) |>
    mutate(
      label_x = estimate + 0.045 * amplitude_x,
      
      label_y = y_plot +
        if_else(
          Response == ordem_respostas[1],
          0.18,
          -0.18
        ),
      
      line_xend =
        label_x - 0.01 * amplitude_x
    )
  
  list(
    coefficients = coeficientes_plot,
    labels = dados_rotulos
  )
}


# Prepare coefficients ----------------------------------------------------

ordem_respostas_riqueza_volume <- c(
  "Total richness",
  "Amphibious richness (%)"
)

ordem_respostas_abundancia_volume <- c(
  "Total abundance",
  "Amphibious abundance (%)"
)

coeficientes_riqueza_volume <-
  extrair_coeficientes_volume(
    modelos_riqueza_volume
  )

coeficientes_abundancia_volume <-
  extrair_coeficientes_volume(
    modelos_abundancia_volume
  )

coeficientes_riqueza_volume
coeficientes_abundancia_volume

preparo_riqueza_volume <-
  preparar_coeficientes_volume(
    coeficientes_riqueza_volume,
    ordem_respostas_riqueza_volume
  )

preparo_abundancia_volume <-
  preparar_coeficientes_volume(
    coeficientes_abundancia_volume,
    ordem_respostas_abundancia_volume
  )

coeficientes_plot_riqueza_volume <-
  preparo_riqueza_volume$coefficients

rotulos_riqueza_volume <-
  preparo_riqueza_volume$labels

coeficientes_plot_abundancia_volume <-
  preparo_abundancia_volume$coefficients

rotulos_abundancia_volume <-
  preparo_abundancia_volume$labels


# Common x-axis limits ----------------------------------------------------

# Use the same limits for both panels.
# Reserve additional space on the right for legends.

valores_x_fig4 <- c(
  coeficientes_plot_riqueza_volume$conf.low,
  coeficientes_plot_riqueza_volume$conf.high,
  rotulos_riqueza_volume$label_x,
  
  coeficientes_plot_abundancia_volume$conf.low,
  coeficientes_plot_abundancia_volume$conf.high,
  rotulos_abundancia_volume$label_x,
  
  0
)

amplitude_fig4 <- diff(
  range(valores_x_fig4, na.rm = TRUE)
)

limites_x_fig4 <- c(
  min(valores_x_fig4, na.rm = TRUE) -
    0.06 * amplitude_fig4,
  
  max(valores_x_fig4, na.rm = TRUE) +
    0.65 * amplitude_fig4
)

limites_x_fig4


# Colors, symbols and legend labels --------------------------------------

cores_riqueza_volume <- c(
  "Total richness" = "#77A88D",
  "Amphibious richness (%)" = "#B7866B"
)

simbolos_riqueza_volume <- c(
  "Total richness" = 21,
  "Amphibious richness (%)" = 22
)

cores_abundancia_volume <- c(
  "Total abundance" = "#77A88D",
  "Amphibious abundance (%)" = "#B7866B"
)

simbolos_abundancia_volume <- c(
  "Total abundance" = 21,
  "Amphibious abundance (%)" = 22
)

rotulos_legenda_riqueza_volume <- c(
  "Total richness" = "Community richness",
  "Amphibious richness (%)" = "Amphibious richness (%)"
)

rotulos_legenda_abundancia_volume <- c(
  "Total abundance" = "Community abundance",
  "Amphibious abundance (%)" = "Amphibious abundance (%)"
)


# Shared coefficient plot function ---------------------------------------

criar_plot_coeficientes_volume <- function(
    coeficientes,
    rotulos,
    ordem_respostas,
    cores,
    simbolos,
    rotulos_legenda,
    tag) {
  
  ggplot(coeficientes) +
    
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      color = "grey55",
      linewidth = 0.7
    ) +
    
    geom_segment(
      aes(
        x = conf.low,
        xend = conf.high,
        y = y_plot,
        yend = y_plot,
        color = Response
      ),
      linewidth = 1.3,
      lineend = "round",
      alpha = 0.70,
      show.legend = FALSE
    ) +
    
    geom_segment(
      data = rotulos,
      aes(
        x = estimate,
        xend = estimate,
        y = y_plot,
        yend = label_y,
        color = Response
      ),
      linewidth = 0.9,
      lineend = "round",
      alpha = 0.90,
      show.legend = FALSE
    ) +
    
    geom_segment(
      data = rotulos,
      aes(
        x = estimate,
        xend = line_xend,
        y = label_y,
        yend = label_y,
        color = Response
      ),
      linewidth = 0.9,
      lineend = "round",
      alpha = 0.90,
      show.legend = FALSE
    ) +
    
    geom_point(
      aes(
        x = estimate,
        y = y_plot,
        fill = Response,
        shape = Response
      ),
      size = 5.6,
      color = "black",
      stroke = 0.65,
      alpha = 0.98,
      show.legend = TRUE
    ) +
    
    geom_text(
      data = rotulos,
      aes(
        x = label_x,
        y = label_y,
        label = coefficient_label,
        color = Response
      ),
      family = "serif",
      fontface = "bold",
      size = 6.2,
      hjust = 0,
      vjust = 0.5,
      alpha = 0.95,
      show.legend = FALSE
    ) +
    
    scale_x_continuous(
      expand = expansion(mult = c(0, 0))
    ) +
    
    scale_y_continuous(
      breaks = seq_along(ordem_termos_volume),
      labels = rev(ordem_termos_volume),
      
      limits = c(
        0.45,
        length(ordem_termos_volume) + 0.55
      ),
      
      expand = expansion(
        mult = c(0.01, 0.01)
      )
    ) +
    
    scale_color_manual(
      values = cores,
      guide = "none"
    ) +
    
    scale_fill_manual(
      values = cores,
      breaks = ordem_respostas,
      labels = rotulos_legenda[ordem_respostas],
      name = NULL
    ) +
    
    scale_shape_manual(
      values = simbolos,
      guide = "none"
    ) +
    
    guides(
      fill = guide_legend(
        override.aes = list(
          shape = unname(
            simbolos[ordem_respostas]
          ),
          size = 5.5,
          color = "black",
          alpha = 1
        )
      )
    ) +
    
    coord_cartesian(
      xlim = limites_x_fig4,
      clip = "off"
    ) +
    
    labs(
      tag = tag,
      x = "Standardized coefficient (95% CI)",
      y = NULL
    ) +
    
    theme_classic(
      base_size = 22,
      base_family = "serif"
    ) +
    
    theme(
      axis.title.x = element_text(
        size = 21,
        margin = margin(t = 12)
      ),
      
      axis.text.x = element_text(
        size = 18,
        color = "black"
      ),
      
      axis.text.y = element_text(
        face = "bold",
        size = 16,
        color = "black"
      ),
      
      axis.ticks = element_line(
        linewidth = 0.6
      ),
      
      plot.tag = element_text(
        family = "serif",
        face = "bold",
        size = 25,
        color = "black"
      ),
      
      plot.tag.position = c(
        0.015,
        0.985
      ),
      
      legend.position = "inside",
      
      legend.position.inside = c(
        1.10,
        0.97
      ),
      
      legend.justification = c(
        1,
        1
      ),
      
      legend.text = element_text(
        family = "serif",
        size = 16
      ),
      
      legend.key.height = grid::unit(
        0.85,
        "cm"
      ),
      
      legend.key.width = grid::unit(
        0.85,
        "cm"
      ),
      
      legend.spacing.y = grid::unit(
        0.05,
        "cm"
      ),
      
      legend.background = element_rect(
        fill = scales::alpha("white", 0.95),
        color = "grey75",
        linewidth = 0.3
      ),
      
      legend.margin = margin(
        6,
        8,
        6,
        8
      ),
      
      plot.margin = margin(
        18,
        50,
        16,
        16
      )
    )
}


# Figure 4A - Richness ----------------------------------------------------

fig4_A <- criar_plot_coeficientes_volume(
  coeficientes = coeficientes_plot_riqueza_volume,
  rotulos = rotulos_riqueza_volume,
  ordem_respostas = ordem_respostas_riqueza_volume,
  cores = cores_riqueza_volume,
  simbolos = simbolos_riqueza_volume,
  rotulos_legenda = rotulos_legenda_riqueza_volume,
  tag = "A"
)

fig4_A


# Figure 4B - Abundance ---------------------------------------------------

fig4_B <- criar_plot_coeficientes_volume(
  coeficientes = coeficientes_plot_abundancia_volume,
  rotulos = rotulos_abundancia_volume,
  ordem_respostas = ordem_respostas_abundancia_volume,
  cores = cores_abundancia_volume,
  simbolos = simbolos_abundancia_volume,
  rotulos_legenda = rotulos_legenda_abundancia_volume,
  tag = "B"
)

fig4_B


# Complete Figure 4 -------------------------------------------------------

Figure_4 <- (
  fig4_A | fig4_B
) +
  patchwork::plot_layout(
    widths = c(1, 1)
  )

Figure_4


# Save --------------------------------------------------------------------

ggsave(
  filename = "Figure_4.jpg",
  plot = Figure_4,
  width = 20,
  height = 9,
  units = "in",
  dpi = 600,
  bg = "white"
)


# Figure 5 ---------------------------------------------------------------

# Interaction prediction function ----------------------------------------

predizer_interacao_volume <- function(
    modelo,
    dados,
    moderador = "taxdist_z",
    percentual = TRUE) {
  
  niveis_z <- c(-1, 0, 1)
  
  # Prediction grid: volume varies across its observed range,
  # while taxonomic distinctness is fixed at -1, 0 and +1 SD.
  
  grade <- tidyr::crossing(
    volume_z = seq(
      min(dados$volume_z, na.rm = TRUE),
      max(dados$volume_z, na.rm = TRUE),
      length.out = 300
    ),
    moderator_z = niveis_z
  ) |>
    mutate(
      distance_z = 0,
      hii_z = 0,
      taxdist_z = 0,
      bio15_z = 0,
      
      Dataset = factor(
        levels(dados$Dataset)[1],
        levels = levels(dados$Dataset)
      )
    )
  
  grade[[moderador]] <- grade$moderator_z
  
  # Back-transform log10-standardized volume.
  
  media_log_volume <- mean(
    dados$log_volume,
    na.rm = TRUE
  )
  
  sd_log_volume <- sd(
    dados$log_volume,
    na.rm = TRUE
  )
  
  grade <- grade |>
    mutate(
      Volume = 10^(
        media_log_volume +
          volume_z * sd_log_volume
      ),
      
      Moderator_level = factor(
        moderator_z,
        levels = c(-1, 0, 1),
        labels = c(
          "Low (−1 SD)",
          "Average",
          "High (+1 SD)"
        )
      )
    )
  
  # Population-level predictions, excluding random effects.
  
  predicao <- predict(
    modelo,
    newdata = grade,
    type = "link",
    se.fit = TRUE,
    re.form = NA,
    allow.new.levels = TRUE
  )
  
  inversa_link <- family(modelo)$linkinv
  
  multiplicador <- if (percentual) 100 else 1
  
  grade |>
    mutate(
      estimate = multiplicador *
        inversa_link(
          as.numeric(predicao$fit)
        ),
      
      CI_low = multiplicador *
        inversa_link(
          as.numeric(
            predicao$fit -
              1.96 * predicao$se.fit
          )
        ),
      
      CI_high = multiplicador *
        inversa_link(
          as.numeric(
            predicao$fit +
              1.96 * predicao$se.fit
          )
        )
    )
}


# Predictions -------------------------------------------------------------

## Figure 5A - Amphibious richness proportion ----------------------------

predicao_fig5_A <- predizer_interacao_volume(
  modelo = modelo_S_amphibious_volume,
  dados = dados_riqueza_volume,
  moderador = "taxdist_z",
  percentual = TRUE
)


## Figure 5B - Amphibious abundance proportion ---------------------------

predicao_fig5_B <- predizer_interacao_volume(
  modelo = modelo_A_amphibious_volume,
  dados = dados_abundancia_volume,
  moderador = "taxdist_z",
  percentual = TRUE
)

head(predicao_fig5_A)
head(predicao_fig5_B)


# Raw moderator values ----------------------------------------------------

# Values of average taxonomic distinctness on its original scale,
# calculated separately for richness and abundance datasets.

niveis_moderadores_fig5 <- bind_rows(
  tibble(
    Panel = "A",
    Moderator = "Average taxonomic distinctness",
    
    Level = c(
      "Low (−1 SD)",
      "Average",
      "High (+1 SD)"
    ),
    
    z = c(-1, 0, 1),
    
    Raw_value = 10^(
      mean(
        dados_riqueza_volume$log_taxdist,
        na.rm = TRUE
      ) +
        z * sd(
          dados_riqueza_volume$log_taxdist,
          na.rm = TRUE
        )
    )
  ),
  
  tibble(
    Panel = "B",
    Moderator = "Average taxonomic distinctness",
    
    Level = c(
      "Low (−1 SD)",
      "Average",
      "High (+1 SD)"
    ),
    
    z = c(-1, 0, 1),
    
    Raw_value = 10^(
      mean(
        dados_abundancia_volume$log_taxdist,
        na.rm = TRUE
      ) +
        z * sd(
          dados_abundancia_volume$log_taxdist,
          na.rm = TRUE
        )
    )
  )
)

niveis_moderadores_fig5


# Colors and line types ---------------------------------------------------

cores_interacao_fig5 <- c(
  "Low (−1 SD)" = "#BFA9D4",
  "Average" = "#7899B5",
  "High (+1 SD)" = "#6B625E"
)

tipos_linha_fig5 <- c(
  "Low (−1 SD)" = "dotted",
  "Average" = "dashed",
  "High (+1 SD)" = "solid"
)


# Plot function -----------------------------------------------------------

criar_plot_interacao_volume <- function(
    dados,
    nome_moderador,
    eixo_y,
    percentual = TRUE) {
  
  intervalo_log_x <- range(
    log10(dados$Volume),
    na.rm = TRUE
  )
  
  posicao_moderador_x <- 10^(
    intervalo_log_x[1] +
      0.72 * diff(intervalo_log_x)
  )
  
  grafico <- ggplot(
    dados,
    aes(
      x = Volume,
      y = estimate,
      color = Moderator_level,
      fill = Moderator_level,
      linetype = Moderator_level,
      group = Moderator_level
    )
  ) +
    
    # 95% confidence intervals.
    
    geom_ribbon(
      aes(
        ymin = CI_low,
        ymax = CI_high
      ),
      color = NA,
      linetype = "blank",
      alpha = 0.14
    ) +
    
    # Predicted amphibious proportions.
    
    geom_line(
      linewidth = 1.6,
      alpha = 0.95
    ) +
    
    # Moderator annotation.
    
    annotate(
      geom = "label",
      x = posicao_moderador_x,
      y = Inf,
      label = paste0(
        "Moderator:\n",
        nome_moderador
      ),
      hjust = 0.5,
      vjust = 1.15,
      family = "serif",
      fontface = "bold",
      size = 5,
      lineheight = 1.05,
      color = "grey20",
      fill = scales::alpha("white", 0.92),
      linewidth = 0.35,
      label.padding = grid::unit(
        0.32,
        "lines"
      ),
      label.r = grid::unit(
        0.12,
        "lines"
      )
    ) +
    
    # Log10 pool volume axis.
    
    scale_x_log10(
      labels = scales::label_number(
        big.mark = ","
      )
    ) +
    
    # Color scale.
    
    scale_color_manual(
      values = cores_interacao_fig5,
      breaks = names(cores_interacao_fig5),
      name = "Moderator level",
      drop = FALSE
    ) +
    
    scale_fill_manual(
      values = cores_interacao_fig5,
      guide = "none",
      drop = FALSE
    ) +
    
    scale_linetype_manual(
      values = tipos_linha_fig5,
      guide = "none",
      drop = FALSE
    ) +
    
    guides(
      color = guide_legend(
        nrow = 1,
        byrow = TRUE,
        override.aes = list(
          linewidth = 2,
          linetype = unname(tipos_linha_fig5),
          alpha = 1
        )
      )
    ) +
    
    labs(
      x = expression(
        "Pool volume (" * m^3 * ")"
      ),
      y = eixo_y
    ) +
    
    theme_classic(
      base_size = 19,
      base_family = "serif"
    ) +
    
    theme(
      axis.title.x = element_text(
        size = 20,
        margin = margin(t = 12)
      ),
      
      axis.title.y = element_text(
        size = 20,
        margin = margin(r = 12)
      ),
      
      axis.text.x = element_text(
        size = 17,
        color = "black"
      ),
      
      axis.text.y = element_text(
        size = 17,
        color = "black"
      ),
      
      axis.line = element_line(
        linewidth = 0.7
      ),
      
      axis.ticks = element_line(
        linewidth = 0.7
      ),
      
      legend.title = element_text(
        family = "serif",
        face = "bold",
        size = 17
      ),
      
      legend.text = element_text(
        family = "serif",
        size = 16
      ),
      
      plot.margin = margin(
        20,
        20,
        15,
        15
      )
    )
  
  if (percentual) {
    grafico <- grafico +
      scale_y_continuous(
        limits = c(0, 100),
        breaks = seq(
          0,
          100,
          by = 20
        ),
        labels = scales::label_percent(
          scale = 1
        )
      )
  }
  
  grafico
}


# Common panel tag theme --------------------------------------------------

tema_tags_fig5 <- theme(
  plot.tag = element_text(
    family = "serif",
    face = "bold",
    size = 25,
    color = "black"
  ),
  
  plot.tag.position = c(
    0.015,
    0.985
  ),
  
  plot.margin = margin(
    20,
    20,
    15,
    15
  )
)


# Figure 5A ---------------------------------------------------------------

fig5_A <- criar_plot_interacao_volume(
  dados = predicao_fig5_A,
  nome_moderador = "Average taxonomic\ndistinctness",
  eixo_y = "Estimated amphibious richness (%)",
  percentual = TRUE
) +
  labs(tag = "A") +
  tema_tags_fig5

fig5_A


# Figure 5B ---------------------------------------------------------------

fig5_B <- criar_plot_interacao_volume(
  dados = predicao_fig5_B,
  nome_moderador = "Average taxonomic\ndistinctness",
  eixo_y = "Estimated amphibious abundance (%)",
  percentual = TRUE
) +
  labs(tag = "B") +
  tema_tags_fig5

fig5_B


# Complete Figure 5 -------------------------------------------------------

Figure_5 <- (
  fig5_A | fig5_B
) +
  patchwork::plot_layout(
    ncol = 2,
    widths = c(1, 1),
    guides = "collect"
  ) &
  theme(
    legend.position = "bottom",
    legend.direction = "horizontal",
    
    legend.key.width = grid::unit(
      1.8,
      "cm"
    ),
    
    legend.key.height = grid::unit(
      0.65,
      "cm"
    ),
    
    legend.spacing.x = grid::unit(
      0.3,
      "cm"
    ),
    
    legend.text = element_text(
      family = "serif",
      size = 16
    ),
    
    legend.title = element_text(
      family = "serif",
      face = "bold",
      size = 17
    ),
    
    legend.margin = margin(
      t = 12
    )
  )

Figure_5


# Save --------------------------------------------------------------------

ggsave(
  filename = "Figure_5.jpg",
  plot = Figure_5,
  width = 15,
  height = 7.5,
  units = "in",
  dpi = 600,
  bg = "white"
)

