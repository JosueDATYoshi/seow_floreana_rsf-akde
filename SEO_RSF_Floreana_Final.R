###─────────────────────────────────────────────────────────────────────────###
# Short-eared Owl RSF – Floreana Island, Galápagos ####
#
# 
# 
# 
# Resource Selection Function (RSF) analysis of the Galápagos Short-eared Owl
# (Asio flammeus galapagoensis) using continuous-time movement modelling (ctmm).
# Analysis covers 16 individuals tracked on Floreana Island.
#
# Josue Arteaga-Torres | josue.arteaga.t@gmail.com
###─────────────────────────────────────────────────────────────────────────###
#
# SCRIPT INDEX
# ─────────────────────────────────────────────────────────────────────────────
# SEC. 1  Packages
# SEC. 2  System setup (working directory, seed, parallel workers)
# SEC. 3  Data input (telemetry, rasters, shapefiles, results table)
#
# SEC. 4  Modelling
#   4.1   Movement model selection  (ctmm.select)
#   4.2   Autocorrelated KDE        (AKDE)
#   4.3   Resource Selection Function (rsf.fit)
#   4.4   RSF-informed AKDE         (AKDE+RSF)
#
# SEC. 5  Summary tables
#   5.1   Population-level RSF coefficients
#   5.2   Individual-level RSF coefficients
#   5.3   Combined data frame for figures
#
# SEC. 6  Figures
#   Fig. 1    3×3 RSF forest plot (individuals + population)
#   Fig. 2    Panel A – Home range meta-analysis (ctmm::meta)
#             Panel B – AKDE+RSF spatial map, scale bar, Galápagos inset
#   S.Fig. 1  Population-level RSF forest plot (all variables, sorted)
#   S.Fig. 2  Land use + road + AKDE overlap map
# ─────────────────────────────────────────────────────────────────────────────


###─────────────────────────────────────────────────────────────────────────###
# SEC. 1 · PACKAGES ####
###─────────────────────────────────────────────────────────────────────────###

rm(list = ls())
gc()

# ctmm is under active development; update from GitHub if CRAN is outdated:
# remotes::install_github("ctmm-initiative/ctmm")
library(ctmm)           # continuous-time movement modelling, AKDE, RSF
library(future.apply)   # parallel model fitting via futures
library(parallelly)     # query available CPU cores

library(dplyr)          # data manipulation (filter, mutate, summarise)
library(tidyr)          # data reshaping (pivot, nest)
library(tibble)         # rownames_to_column and modern data frames
library(readxl)         # read Excel files (.xlsx)

library(sf)             # vector spatial data (shapefiles, polygons)
library(terra)          # raster handling (modern replacement for raster pkg)
library(raster)         # legacy raster format (still required by ctmm)
library(tidyterra)      # ggplot2 geoms for terra SpatRaster objects
library(ggnewscale)     # multiple fill/colour scales in one ggplot
library(scales)         # alpha() and other scale helpers
library(ggspatial)      # scale bar and north arrow annotations for ggplot maps
library(rnaturalearth)  # world/country outlines for context inset map
library(rnaturalearthdata)

library(ggplot2)        # general plotting framework
library(patchwork)      # compose multiple ggplot panels into one figure
library(cowplot)        # plot_grid() for aligned mixed base-R / ggplot figures
library(RColorBrewer)   # colour palettes for multi-individual maps
library(here)           # stablishing the working directory directly file path

###─────────────────────────────────────────────────────────────────────────###
# SEC. 2 · SYSTEM SETUP ####
###─────────────────────────────────────────────────────────────────────────###

setwd(here())

# Seed ensures reproducibility of Monte Carlo integrator in rsf.fit()
set.seed(3847)

# Distribute model fitting across available CPU cores (cap at 6)
plan("multisession", workers = min(6, parallelly::availableCores()))


###─────────────────────────────────────────────────────────────────────────###
# SEC. 3 · DATA INPUT ####
###─────────────────────────────────────────────────────────────────────────###

# ── 3.1 Telemetry and movement data ####
###─────────────────────────────────────────

# Pre-filtered list (PS0024 and PS0028 already excluded).
# Run SEO_cleanup_publication.R once to generate this file.
load("./outputs/publication/SEO_owls_pub.RData")  # -> owls_tele_pub, owls_move_pub

x <- owls_tele

# cat("Individuals in analysis:", paste(names(x), collapse = ", "), "\n")
# cat("N individuals:", length(x), "\n")

## ── 3.2 Environmental rasters ####
###────────────────────────────────────────────────

# z-scored continuous variables; categorical land-use layers combined into a
# single factor raster. Standardisation preserves cross-variable comparability.
owls_all_raster <- readRDS("./Rasters/SEO_standar_and_combined.rds")
owls_env_list   <- unstack(owls_all_raster)
names(owls_env_list) <- names(owls_all_raster)

# Assign individual layer objects to the global environment for direct access
for (nm in names(owls_all_raster)) assign(nm, owls_all_raster[[nm]])

# Full raster list passed to RSF; Forest (LandUse class 1) is the reference
Rlist            <- owls_env_list
categorical_list <- list(LandUse = 1)

# ── 3.3 Spatial layers ####
###───────────────────────────────────────────────────────

# All spatial layers are pre-processed and stored in a single geopackage.
# Run SEO_cleanup_publication.R once to generate this file.
# Layers: galapagos (all islands, Floreana labelled), floreana_border, road
gpkg_path <- "./outputs/publication/SEO_spatial_layers.gpkg"

# Island outline in UTM Zone 15S (EPSG:32715); used for masking and plotting
floreana <- sf::st_read(gpkg_path, layer = "floreana_border", quiet = TRUE)

# Main road centreline; already reprojected to UTM 15S
road_sf <- sf::st_read(gpkg_path, layer = "road", quiet = TRUE)

# Galápagos archipelago islands; UTM 15S, "Isla Santa María" already labelled
# as "Floreana" for easy filtering.
galapagos_islands <- sf::st_read(gpkg_path, layer = "galapagos", quiet = TRUE)

galapagos_flor <- galapagos_islands |>
  dplyr::filter(NOMBRE == "Floreana")
galapagos_rest <- galapagos_islands |>
  dplyr::filter(is.na(NOMBRE) | NOMBRE != "Floreana")

###─────────────────────────────────────────────────────────────────────────###
# SEC. 4 · MODELLING ####
###─────────────────────────────────────────────────────────────────────────###

###─────────────────────────────────────────────────────────────────────────###
## 4.1 · Movement model selection (ctmm.select) ####
###─────────────────────────────────────────────────────────────────────────###
#
# AICc-based selection across candidate continuous-time movement models
# (OUF, OU, IID, and anisotropic variants). isotropic = TRUE constrains the
# home range to a circular covariance structure, appropriate for these owls.

# guess <- lapply(x, ctmm.guess,
#                 CTMM        = ctmm(isotropic = TRUE),
#                 interactive = FALSE)
#
# ctmm_fit <- future_mapply(ctmm.select, x, guess,
#                            SIMPLIFY    = FALSE,
#                            future.seed = TRUE)
#
# save(ctmm_fit, file = "./outputs/models/ctmm_select_Final.RData")

load("./outputs/models/ctmm_select_Final.RData")   # loads object named 'ctmm'
ctmm_fit <- ctmm; rm(ctmm)   # rename to avoid masking the ctmm() function

###─────────────────────────────────────────────────────────────────────────###
## 4.2 · Autocorrelated Kernel Density Estimate (AKDE) ####
###─────────────────────────────────────────────────────────────────────────###
#
# A shared grid aligned to Rlist[[1]] ensures all AKDE outputs and
# environmental rasters share the same spatial extent and resolution,
# which is required for RSF estimation.

# akde_fit <- akde(x, ctmm_fit, grid = Rlist[[1]])
# save(akde_fit, file = "./outputs/AKDEs/AKDE_Final.RData")

load("./outputs/AKDEs/AKDE_Final.RData")   # loads object named 'akde'
akde_fit <- akde; rm(akde)

###─────────────────────────────────────────────────────────────────────────###
## 4.3 · Resource Selection Function (rsf.fit) ####
###─────────────────────────────────────────────────────────────────────────###
#
# Monte Carlo integrator is recommended for complex, real-world raster
# landscapes. reference = categorical_list sets Forest as the baseline level
# for the LandUse factor (all other classes are estimated relative to Forest).

# rsf_fit <- future_mapply(rsf.fit, x, akde_fit,
#                           MoreArgs    = list(R          = Rlist,
#                                              integrator = "MonteCarlo",
#                                              reference  = categorical_list),
#                           SIMPLIFY    = FALSE,
#                           future.seed = TRUE)
#
# save(rsf_fit, file = "./outputs/RSFs/RSF_Final_06.24.RData")

load("./outputs/RSFs/RSF_Final_06.24.RData")   # loads object named 'rsf'
rsf_fit <- rsf; rm(rsf)

###─────────────────────────────────────────────────────────────────────────###
## 4.4 · RSF-informed AKDE (AKDE+RSF) ####
###─────────────────────────────────────────────────────────────────────────###
#
# Incorporates estimated habitat selection coefficients into the utilisation
# distribution, producing a home range that reflects both movement autocorrelation
# and resource preferences. This is the primary output used in home range maps.

# akde_rsf_fit <- akde(x, CTMM = rsf_fit, R = Rlist, grid = Rlist[[1]])
# save(akde_rsf_fit, file = "./outputs/AKDE+RSF/AKDERSF_Final_06.24.RData")

load("./outputs/AKDE+RSF/AKDERSF_Final_06.24.RData")   # loads 'akde_rsf'
akde_rsf_fit <- akde_rsf; rm(akde_rsf)


###─────────────────────────────────────────────────────────────────────────###
# SEC. 5 · SUMMARY TABLES ####
###─────────────────────────────────────────────────────────────────────────###

###─────────────────────────────────────────────────────────────────────────###
## 5.1 · Population-level RSF coefficients ####
###─────────────────────────────────────────────────────────────────────────###
#' Calculation of population-level rsf
# mean_rsf <- mean(rsf)

# save(mean_rsf, file = here::here("outputs", "/Mean_Population_RSF_Final.RData"))

load("./outputs/Mean_Population_RSF_Final.RData")   # -> mean_rsf

mean_rsf_summary <- summary(mean_rsf)

# CI matrix: rows = model parameters (RSF + movement); columns = low, est, high
ci_raw <- as.data.frame(mean_rsf_summary$CI)
colnames(ci_raw)[1:3] <- c("low", "estimate", "high")

# Retain only RSF variable rows (skip movement and autocorrelation parameters)
ci_df <- ci_raw[c(1, 3, 4, 6, 8, 10:11, 13, 15, 17), ]

# Assign human-readable labels (reference category Forest excluded from table)
pop_var_labels <- c(
  "Land Type-Ocean", "Land Type-Urban", "Land Type-Low vegetation",
  "Land Type-Lava",  "Distance to Road", "Ruggedness Index",
  "Slope Angle",     "Elevation",        "Tree Cover",
  "Vegetation Index (NDVI)"
)
rownames(ci_df) <- pop_var_labels
ci_df$term      <- pop_var_labels

# write.csv(ci_df, "./outputs/RSF_pop_final.csv", row.names = TRUE)

# ── 3.4 Pre-processed RSF variable name table ####
###────────────────────────────────

# Long-format table of RSF variable labels harmonised across individuals
# (used to replace ctmm's internal row names with readable names)
rsf_long <- read_excel("./outputs/RSF_results_final_long.xlsx", sheet = "Sheet1")


###─────────────────────────────────────────────────────────────────────────###
## 5.2 · Individual-level RSF coefficients ####
###─────────────────────────────────────────────────────────────────────────###

rsf_variables <- c("NDVI", "TreeCover", "Elevation", "Slope",
                   "Ruggedness", "RoadProximity", "LandUse")

result_ind <- do.call(rbind, lapply(rsf_fit, function(fit) {
  df <- as.data.frame(summary(fit)$CI)
  df <- df[grep(paste(rsf_variables, collapse = "|"), rownames(df)), ]
  tibble::rownames_to_column(df, var = "variable")
}))
result_ind$animal_id <- rep(names(x), each = 11)

# Replace ctmm's internal parameter names with readable variable labels
result_ind$variable <- rsf_long$Variable

# write.csv(result_ind, "./outputs/RSF_Results_MC_Final.csv", row.names = FALSE)

###─────────────────────────────────────────────────────────────────────────###
## 5.3 · Combined data frame (individuals + population) ####
###─────────────────────────────────────────────────────────────────────────###

# Mapping from raw labels to harmonised display names used across all figures
name_map <- c(
  "Lava Field"     = "Land Type-Lava",
  "Ocean"          = "Land Type-Ocean",
  "Urban Area"     = "Land Type-Urban",
  "Low Vegetation" = "Land Type-Low vegetation",
  "Slope angle"    = "Slope Angle",
  "Fresh Water"    = "Land Type-Fresh Water",
  "NDVI"           = "Vegetation Index (NDVI)"
)

result_ind <- result_ind |>
  dplyr::mutate(term = dplyr::recode(variable, !!!name_map, .default = variable))

ci_df <- ci_df |>
  dplyr::mutate(term = dplyr::recode(term, !!!name_map, .default = term))

# Population-based CI caps prevent extreme individual estimates from distorting
# axis scales (individual cap = population estimate ± 5 × population CI width)
cap_table <- ci_df |>
  dplyr::mutate(
    ci_width = high - low,
    cap_low  = estimate - 5 * ci_width,
    cap_high = estimate + 5 * ci_width
  ) |>
  dplyr::select(term, cap_low, cap_high)

# Fresh Water has no population-level estimate; borrow caps from Ocean
cap_table <- dplyr::bind_rows(
  cap_table,
  cap_table |>
    dplyr::filter(term == "Land Type-Ocean") |>
    dplyr::mutate(term = "Land Type-Fresh Water")
)

ind_df <- result_ind |>
  dplyr::rename(estimate = est) |>
  dplyr::mutate(source = "Individual")

pop_df <- ci_df |>
  dplyr::mutate(animal_id = "Population Mean", source = "Population")

plot_df <- dplyr::bind_rows(ind_df, pop_df) |>
  dplyr::left_join(cap_table, by = "term") |>
  dplyr::mutate(
    # Flag values that were clipped to the population-derived cap
    truncated = (low < cap_low | high > cap_high |
                   estimate < cap_low | estimate > cap_high),
    low       = pmax(low,      cap_low,  na.rm = TRUE),
    high      = pmin(high,     cap_high, na.rm = TRUE),
    estimate  = pmin(pmax(estimate, cap_low, na.rm = TRUE), cap_high, na.rm = TRUE),
    significance = dplyr::case_when(
      low > 0  ~ "Positive",
      high < 0 ~ "Negative",
      TRUE     ~ "Non-significant"
    )
  )


###─────────────────────────────────────────────────────────────────────────###
# SEC. 6 · FIGURES ####
###─────────────────────────────────────────────────────────────────────────###

# Shared significance colour palette (Okabe-Ito; well-separated luminance
# values ensure the three categories remain distinguishable in greyscale):
#   Positive:        sky blue   #56B4E9  (luminance ≈ 0.42)
#   Negative:        vermillion #D55E00  (luminance ≈ 0.23)
#   Non-significant: gray80              (luminance ≈ 0.60)
sig_colors <- c(
  "Positive"        = "#56B4E9",
  "Negative"        = "#D55E00",
  "Non-significant" = "gray80"
)
src_shapes <- c("Individual" = 16, "Population" = 18)


###─────────────────────────────────────────────────────────────────────────###
## Fig. 1 · RSF forest plot: 3 × 3 individual + population grid ####
###─────────────────────────────────────────────────────────────────────────###
#
# Layout (row × column):
#   Row 1 – Vegetation Index | Tree Cover   | Low Vegetation
#   Row 2 – Distance to Road | Urban        | Lava
#   Row 3 – Elevation        | Ruggedness   | Slope
#
# Within each panel: circles = individual estimates; diamonds = population mean.
# Faded points/bars indicate values clipped by the population-based axis cap.
# Ocean and Fresh Water omitted (not retained in final model).

plot_df_fig1 <- plot_df |>
  dplyr::filter(!term %in% c("Land Type-Fresh Water", "Land Type-Ocean"))

var_order_fig1 <- c(
  "Vegetation Index (NDVI)", "Tree Cover",          "Land Type-Low vegetation",
  "Distance to Road",        "Land Type-Urban",      "Land Type-Lava",
  "Elevation",               "Ruggedness Index",     "Slope Angle"
)

var_labels_fig1 <- c(
  "Vegetation Index (NDVI)"  = "Vegetation Index",
  "Tree Cover"               = "Tree Cover",
  "Land Type-Low vegetation" = "Low Vegetation",
  "Distance to Road"         = "Distance to Road",
  "Land Type-Urban"          = "Urban",
  "Land Type-Lava"           = "Lava",
  "Elevation"                = "Elevation",
  "Ruggedness Index"         = "Ruggedness Index",
  "Slope Angle"              = "Slope Angle"
)

# Helper function: produce one forest panel per environmental variable
make_panel_fig1 <- function(vname, df) {
  sub <- df |> dplyr::filter(term == vname)
  ggplot(sub, aes(x = estimate, y = animal_id)) +
    geom_point(
      data = \(d) dplyr::filter(d, source == "Individual"),
      aes(color = significance, shape = source, alpha = truncated),
      size = 1.8
    ) +
    geom_errorbarh(
      data = \(d) dplyr::filter(d, source == "Individual"),
      aes(xmin = low, xmax = high, color = significance, alpha = truncated),
      height = 0.2
    ) +
    geom_point(
      data = \(d) dplyr::filter(d, source == "Population"),
      aes(color = significance, shape = source),
      size = 4
    ) +
    geom_errorbarh(
      data = \(d) dplyr::filter(d, source == "Population"),
      aes(xmin = low, xmax = high, color = significance),
      linewidth = 1.2, height = 0.3
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    scale_color_manual(
      values = sig_colors,
      limits = c("Positive", "Negative", "Non-significant"),
      drop   = FALSE,
      guide  = "none"
    ) +
    scale_shape_manual(values = src_shapes, guide = "none") +
    scale_alpha_manual(
      values = c(`TRUE` = 0.35, `FALSE` = 1),
      guide  = "none"
    ) +
    labs(x = NULL, y = NULL, title = var_labels_fig1[vname]) +
    theme_minimal(base_size = 9) +
    theme(
      legend.position = "none",
      plot.title      = element_text(size = 8, face = "bold"),
      axis.text.y     = element_text(size = 7)
    )
}

fig1 <- wrap_plots(
  lapply(var_order_fig1, make_panel_fig1, df = plot_df_fig1),
  ncol = 3
) +
  plot_annotation(
    caption = "RSF Coefficient Estimate",
    theme   = theme(plot.caption = element_text(hjust = 0.5, size = 9))
  )

fig1

# ggsave("./outputs/figures/Fig1_RSF_3x3.pdf",  fig1, width = 10, height = 10)
# ggsave("./outputs/figures/Fig1_RSF_3x3.tiff", fig1, width = 10, height = 10,
#         dpi = 600, compression = "lzw")


###─────────────────────────────────────────────────────────────────────────###
## Fig. 2 · Home range meta-analysis (A) + AKDE+RSF spatial map (B) ####
###─────────────────────────────────────────────────────────────────────────###

# ── Individual colour palette (shared between Panels A and B) ────────────────
# Draws from Set1 + Dark2; red removed and the 5th colour muted to avoid
# confusion with the significance palette used in Fig. 1.
pal_ind <- c(
  RColorBrewer::brewer.pal(name = "Set1",  n = 9),
  RColorBrewer::brewer.pal(name = "Dark2", n = 8)
)
pal_ind <- pal_ind[-3]        # remove red
pal_ind[5] <- "#888866"       # mute clashing 5th colour

# Assign palette colours to individual IDs in the same order as ctmm's colouring
pal <- ctmm::color(akde_rsf_fit, by = "individual")
for (i in seq_along(pal)) pal[i] <- pal_ind[i]

ind_colors <- setNames(pal_ind[seq_along(akde_rsf_fit)], names(akde_rsf_fit))

###─────────────────────────────────────────────────────────────────────────###
### Fig. 2A · Home range size meta-analysis ####
###─────────────────────────────────────────────────────────────────────────###
#
# ctmm::meta() ranks individuals by home range size (with uncertainty) and
# overlays the population mean (black diamond). Wrapped with wrap_elements()
# so patchwork can compose it alongside the ggplot map.

fig2a <- wrap_elements(full = ~ ctmm::meta(
  akde_rsf_fit,
  col     = c(pal, "black"),
  verbose = FALSE,
  sort    = TRUE,
  labels  = TRUE
))

###─────────────────────────────────────────────────────────────────────────###
### Fig. 2B · AKDE+RSF spatial map with inset ####
###─────────────────────────────────────────────────────────────────────────###
#
# Layers (bottom to top):
#   1. Floreana island outline (light grey fill)
#   2. Per-individual continuous UD gradient clipped to 95% HR boundary:
#      core (UD ≈ 0, most used) near-opaque; edge (UD ≈ 0.95) transparent.
#      Visual style mirrors ctmm::plot(akde_rsf).
#   3. Main road in fixed grey (described in figure caption; no legend entry).
#   4. Scale bar (bottom-left).
# Inset: all Galápagos islands in uniform grey; Floreana study area marked
#   with an orange bounding-box rectangle (no contrasting island colour).

# ── Continuous UD rasters ────────────────────────────────────────────────────
# ctmm stores UD as p-values: 0 = most intensely used, 1 = outside HR.
akde_ud_list <- lapply(akde_rsf_fit, function(ud) {
  r        <- raster::raster(ud)
  r[r > 0.95] <- NA    # mask cells beyond the 95% contour
  terra::rast(r)
})

# One geom_spatraster block per individual.
# NOTE: ggnewscale layers cannot be reused across ggplot() calls; always
#   build fresh with this mapply() pattern.
ind_overlay_layers <- unlist(
  mapply(function(r, col) {
    list(
      new_scale_fill(),
      geom_spatraster(data = r),
      scale_fill_gradient(
        low      = scales::alpha(col, 0.9),
        high     = scales::alpha(col, 0.05),
        na.value = NA,
        guide    = "none"
      )
    )
  },
  akde_ud_list,
  ind_colors[names(akde_rsf_fit)],
  SIMPLIFY = FALSE),
  recursive = FALSE
)

# ── Galápagos context inset ───────────────────────────────────────────────────
# All islands in uniform grey. An orange bounding-box rectangle locates
# Floreana without singling out its polygon colour.
# Expand Floreana's bounding box by 15 km on each side so the indicator
# rectangle is clearly visible at the full-archipelago inset scale.
flor_bb <- sf::st_bbox(galapagos_flor)
pad     <- 15000   # metres (UTM 15S)
floreana_bbox_sf <- sf::st_as_sfc(
  structure(
    c(xmin = unname(flor_bb["xmin"]) - pad,
      ymin = unname(flor_bb["ymin"]) - pad,
      xmax = unname(flor_bb["xmax"]) + pad,
      ymax = unname(flor_bb["ymax"]) + pad),
    class = "bbox",
    crs   = sf::st_crs(galapagos_flor)
  )
)

fig2b_inset <- ggplot() +
  geom_sf(
    data      = galapagos_islands,
    fill      = "gray75",
    color     = "gray40",
    linewidth = 0.2
  ) +
  geom_sf(
    data      = floreana_bbox_sf,
    fill      = NA,
    color     = "#D55E00",
    linewidth = 0.8
  ) +
  theme_void() +
  theme(
    panel.border = element_rect(color = "gray30", fill = NA, linewidth = 0.5)
  )

# ── Main map (Panel B) ───────────────────────────────────────────────────────
# No tag added here — labels are handled by cowplot::plot_grid() below,
# which converts both panels to grobs and guarantees left-edge alignment.
# NOTE: ind_overlay_layers must be rebuilt fresh for each ggplot() call;
#   ggnewscale layers are stateful and cannot be reused across calls.
ind_overlay_layers <- unlist(
  mapply(function(r, col) {
    list(
      new_scale_fill(),
      geom_spatraster(data = r),
      scale_fill_gradient(
        low      = scales::alpha(col, 0.9),
        high     = scales::alpha(col, 0.05),
        na.value = NA,
        guide    = "none"
      )
    )
  },
  akde_ud_list,
  ind_colors[names(akde_rsf_fit)],
  SIMPLIFY = FALSE),
  recursive = FALSE
)

fig2b_main <- ggplot() +
  geom_sf(data = floreana, fill = "gray92", color = "gray50", linewidth = 0.4) +
  ind_overlay_layers +
  geom_sf(data = road_sf, color = "gray20", linewidth = 0.6, inherit.aes = FALSE) +
  annotation_scale(location = "bl", width_hint = 0.25) +
  coord_sf(crs = 32715, expand = FALSE) +
  theme_minimal(base_size = 10) +
  theme(
    axis.title      = element_blank(),
    axis.text       = element_text(size = 8),
    legend.position = "none"
  )

fig2b <- fig2b_main +
  inset_element(
    fig2b_inset,
    left = 0.65, bottom = 0.65,
    right = 1.00, top   = 1.00
  )

# ── Compose Fig. 2: Panel A over Panel B ────────────────────────────────────
# cowplot::plot_grid() converts both elements to grobs before laying them out,
# so the outer left/right edges of A and B are guaranteed to be flush.
# rel_heights = c(1, 1.25) gives the requested proportional area split.
fig2 <- cowplot::plot_grid(
  ~ ctmm::meta(
      akde_rsf_fit,
      col     = c(pal, "black"),
      verbose = FALSE,
      sort    = TRUE,
      labels  = TRUE
    ),
  fig2b,
  ncol           = 1,
  rel_heights    = c(1, 1.25),
  labels         = c("A", "B"),
  label_size     = 12,
  label_fontface = "bold"
)

fig2

# ggsave("./outputs/figures/Fig2_HomeRange_Map.pdf",  fig2, width = 8, height = 14)
# ggsave("./outputs/figures/Fig2_HomeRange_Map.tiff", fig2, width = 8, height = 14,
#         dpi = 600, compression = "lzw")


###─────────────────────────────────────────────────────────────────────────###
## S.Fig. 1 · Population-level RSF: all variables, sorted by estimate ####
###─────────────────────────────────────────────────────────────────────────###
#
# Full population-level (meta-analysis mean) forest plot for supplementary
# material. Ocean dropped (no meaningful habitat use); remaining variables
# sorted from most negative to most positive selection coefficient.

sfig1_df <- ci_df |>
  dplyr::filter(term != "Land Type-Ocean") |>
  dplyr::mutate(
    significance = dplyr::case_when(
      low > 0  ~ "Positive",
      high < 0 ~ "Negative",
      TRUE     ~ "Non-significant"
    )
  )

sfig1 <- ggplot(
  sfig1_df,
  aes(x = estimate, y = reorder(term, estimate), color = significance)
) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_errorbarh(
    aes(xmin = low, xmax = high),
    height    = 0.25,
    linewidth = 0.8
  ) +
  geom_point(shape = 18, size = 4) +
  scale_color_manual(
    values = sig_colors,
    limits = c("Positive", "Negative", "Non-significant"),
    name   = "Effect",
    drop   = FALSE
  ) +
  labs(
    x = "RSF Coefficient Estimate",
    y = NULL
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10, face = "bold")
  )

sfig1

# ggsave("./outputs/figures/SFig1_PopRSF.pdf",  sfig1, width = 7, height = 5)
# ggsave("./outputs/figures/SFig1_PopRSF.tiff", sfig1, width = 7, height = 5,
#         dpi = 600, compression = "lzw")


###─────────────────────────────────────────────────────────────────────────###
## S.Fig. 2 · Land use + road + AKDE overlap map ####
###─────────────────────────────────────────────────────────────────────────###
#
# Base layer: Floreana land-use raster (Okabe-Ito colours).
# Overlay: per-individual 95% HR rasters stacked with α = 0.25 (dark grey);
#   darker cells indicate higher spatial overlap among individuals.
# Road centreline in white for contrast against the coloured base.
# Galápagos context inset (same as Fig. 2) and scale bar added.

lu_terra <- terra::rast(LandUse)
levels(lu_terra)[[1]]$class[4] <- "Urban"
levels(lu_terra)[[1]]$class[5] <- "Lagoon"

landuse_cols <- c(
  "Forest"         = "#009E73",  # Okabe-Ito bluish green
  "Lava"           = "#999999",  # grey
  "Low vegetation" = "#F0E442",  # yellow
  "Urban"          = "#E69F00",  # amber
  "Lagoon"         = "#56B4E9",  # sky blue
  "Ocean"          = "#0072B2"   # blue
)

layer_alpha <- 0.25   # per-layer opacity; compounding reveals overlap intensity

# Binary 95% HR rasters (1 = inside UD, NA = outside / transparent)
akde_bin_sfig2 <- lapply(akde_rsf_fit, function(ud) {
  r     <- raster::raster(ud)
  r_bin <- r <= 0.95
  r_bin[!raster::values(r_bin)] <- NA
  terra::rast(r_bin)
})

# Stack one dark-grey geom per individual; opacity compounds naturally.
# NOTE: ggnewscale layers cannot be reused across separate ggplot() calls;
# these layers are built fresh here (not shared with Fig. 2B).
overlap_layers_sfig2 <- unlist(lapply(akde_bin_sfig2, function(r) {
  list(
    new_scale_fill(),
    geom_spatraster(data = r),
    scale_fill_gradient(
      low      = scales::alpha("#1A1A1A", 0),
      high     = scales::alpha("#1A1A1A", layer_alpha),
      na.value = NA,
      guide    = "none"
    )
  )
}), recursive = FALSE)

# Main map
# Road plotted with fixed colour (described in caption; no legend entry).
# fig2b_inset reused from Fig. 2 — same bbox-rectangle style, uniform grey islands.
sfig2_main <- ggplot() +
  # Land use base layer
  geom_spatraster(data = lu_terra, aes(fill = class)) +
  scale_fill_manual(values = landuse_cols, name = "Land Use", na.value = NA) +
  # AKDE overlap overlay
  overlap_layers_sfig2 +
  # Road — fixed white colour; described in caption, no legend entry
  geom_sf(data = road_sf, color = "white", linewidth = 0.66, inherit.aes = FALSE) +
  annotation_scale(location = "br", width_hint = 0.25) +
  coord_sf(crs = 32715, expand = FALSE) +
  theme_minimal(base_size = 11) +
  theme(
    axis.title       = element_blank(),
    axis.text        = element_text(size = 8),
    legend.position  = "right",
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    # Ocean cells outside the island are filled with the ocean colour
    panel.background = element_rect(fill = "#0072B2", color = NA)
  )

# Embed the same Galápagos inset used in Fig. 2 (bbox rectangle, uniform grey)
sfig2 <- sfig2_main +
  inset_element(
    fig2b_inset,
    left = 0.65, bottom = 0.65,
    right = 1.00, top   = 1.00
  )

sfig2

# ggsave("./outputs/figures/SFig2_LandUse_AKDE.pdf",  sfig2, width = 8, height = 7)
# ggsave("./outputs/figures/SFig2_LandUse_AKDE.tiff", sfig2, width = 8, height = 7,
#         dpi = 600, compression = "lzw")
