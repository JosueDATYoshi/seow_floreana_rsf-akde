###─────────────────────────────────────────────────────────────────────────###
# Short-eared Owl RSF – Floreana Island, Galápagos ####
#
# Resource Selection Function (RSF) analysis of the Galápagos Short-eared Owl
# (Asio flammeus galapagoensis) using continuous-time movement modelling (ctmm).
# Analysis covers 16 individuals tracked on Floreana Island.
#
# Code Author: Josue Arteaga-Torres 
# Research authors: Josué D. Arteaga-Torres, Shane C. Sumasgutner, 
# Paula A. Castaño, Sonia Kleindorfer, Petra Sumasgutner
# 

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

library(dplyr)          # data manipulation
library(tibble)         # rownames_to_column

library(sf)             # vector spatial data
library(terra)          # raster handling
library(raster)         # legacy raster format (required by ctmm)
library(tidyterra)      # ggplot2 geoms for terra SpatRaster objects
library(ggnewscale)     # multiple fill/colour scales in one ggplot
library(scales)         # alpha() helper
library(ggspatial)      # scale bar annotations for ggplot maps

library(ggplot2)        # plotting
library(patchwork)      # compose ggplot panels
library(cowplot)        # plot_grid() for mixed base-R / ggplot figures
library(RColorBrewer)   # colour palettes
library(here)           # project-relative file paths


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

## 3.1 · Telemetry data ####

load("./inputs/SEO_owls_data.RData")  # -> owls_tele
x <- owls_tele

## 3.2 · Environmental rasters ####
# Continuous raster layers (NDVI, TreeCover, Elevation, Slope, Ruggedness,
# Urban, FarmLand, RoadProximity). Z-score before passing to rsf.fit().

load("./inputs/SEO_rasterlist_NDVI_tree_elev_slope_rug_urb_farm_road.RData")  # -> owls_env_list
Rlist            <- owls_env_list
categorical_list <- list(LandUse = 1)  # Forest = reference category

## 3.3 · Land use raster ####
# Rasterized from ecosystem polygons to the same 30 m UTM 15S grid.
# Five classes: Forest (1, RSF reference), Lava (2), Low vegetation (3),
# Urban (4), Lagoon (5). Ocean cells are NA (rendered as background in S.Fig. 2).

floreana_eco <- sf::st_read("./inputs/floreana_ecosystems.gpkg", quiet = TRUE) |>
  sf::st_transform(32715) |>
  dplyr::mutate(
    lu_int = dplyr::case_when(
      simplified_class %in% c("Deciduous Forest", "Evergreen Seasonal Forest",
                               "Evergreen Shrub and Forest")              ~ 1L,
      simplified_class == "Lava Field"                                    ~ 2L,
      simplified_class %in% c("Deciduous Shrub", "Deciduous Grassland",
                               "Invasive Plants", "Agriculture", 
                              "Agriculture Buffer")                       ~ 3L,
      simplified_class == "Urban"                                         ~ 4L,
      simplified_class == "Water"                                         ~ 5L,
      TRUE                                                                ~ 3L
    )
  )

LandUse <- terra::rasterize(
  terra::vect(floreana_eco),
  terra::rast(owls_env_list[[1]]),
  field = "lu_int"
)
LandUse[is.na(LandUse)] <- 6L
levels(LandUse) <- data.frame(
  value = 1:6,
  class = c("Forest", "Lava", "Low vegetation", "Urban", "Lagoon", "Ocean")
)

## 3.4 · Spatial layers ####
# GeoPackage layers: floreana_border, road, galapagos (all islands)

gpkg_path         <- "./inputs/SEO_spatial_layers.gpkg"
floreana          <- sf::st_read(gpkg_path, layer = "floreana_border", quiet = TRUE)
road_sf           <- sf::st_read(gpkg_path, layer = "road",            quiet = TRUE)
galapagos_islands <- sf::st_read(gpkg_path, layer = "galapagos",       quiet = TRUE)
galapagos_flor    <- galapagos_islands |> dplyr::filter(NOMBRE == "Floreana")


###─────────────────────────────────────────────────────────────────────────###
# SEC. 4 · MODELLING ####
###─────────────────────────────────────────────────────────────────────────###

## 4.1 · Movement model selection (ctmm.select) ####
# AICc-based selection across OUF, OU, IID, and anisotropic variants.
# isotropic = TRUE constrains the range distribution to a circular covariance structure.

# guess    <- lapply(x, ctmm.guess, CTMM = ctmm(isotropic = TRUE), interactive = FALSE)
# ctmm_fit <- future_mapply(ctmm.select, x, guess, SIMPLIFY = FALSE, future.seed = TRUE)
# save(ctmm_fit, file = "./outputs/ctmm_select_Final.RData")

load("./outputs/ctmm_select_Final.RData")
ctmm_fit <- ctmm; rm(ctmm)

## 4.2 · Autocorrelated KDE (AKDE) ####
# Shared grid aligned to Rlist[[1]] ensures AKDE and raster extents match.

# akde_fit <- akde(x, ctmm_fit, grid = Rlist[[1]])
# save(akde_fit, file = "./outputs/AKDE_Final.RData")

load("./outputs/AKDE_Final.RData")
akde_fit <- akde; rm(akde)

## 4.3 · Resource Selection Function (rsf.fit) ####
# Monte Carlo integrator for complex real-world raster landscapes.

# rsf_fit <- future_mapply(rsf.fit, x, akde_fit,
#                          MoreArgs    = list(R          = Rlist,
#                                             integrator = "MonteCarlo",
#                                             reference  = categorical_list),
#                          SIMPLIFY    = FALSE,
#                          future.seed = TRUE)
# save(rsf_fit, file = "./outputs/RSF_Final.RData")

load("./outputs/RSF_Final.RData")
rsf_fit <- rsf; rm(rsf)

## 4.4 · RSF-informed AKDE (AKDE+RSF) ####
# Incorporates habitat selection coefficients into the utilisation distribution.

# akde_rsf_fit <- akde(x, CTMM = rsf_fit, R = Rlist, grid = Rlist[[1]])
# save(akde_rsf_fit, file = "./outputs/AKDERSF_Final.RData")

load("./outputs/AKDERSF_Final.RData")
akde_rsf_fit <- akde_rsf; rm(akde_rsf)

## 4.5 · Population-level RSF (mean) ####
# Meta-analysis mean across individual RSF fits; used to derive population-level coefficients.

# mean_rsf <- mean(rsf_fit)
# save(mean_rsf, file = "./outputs/Mean_Population_RSF_Final.RData")

load("./outputs/Mean_Population_RSF_Final.RData")  # -> mean_rsf


###─────────────────────────────────────────────────────────────────────────###
# SEC. 5 · SUMMARY TABLES ####
###─────────────────────────────────────────────────────────────────────────###

## 5.1 · Population-level RSF coefficients ####

ci_raw <- as.data.frame(summary(mean_rsf)$CI)
colnames(ci_raw)[1:3] <- c("low", "estimate", "high")

# Rows corresponding to RSF variables (skips movement/autocorrelation parameters)
ci_df <- ci_raw[c(1, 3, 4, 6, 8, 10:11, 13, 15, 17), ]

pop_var_labels <- c(
  "Land Type-Ocean", "Land Type-Urban", "Land Type-Low vegetation",
  "Land Type-Lava",  "Distance to Road", "Ruggedness Index",
  "Slope Angle",     "Elevation",        "Tree Cover",
  "Vegetation Index (NDVI)"
)
rownames(ci_df) <- pop_var_labels
ci_df$term      <- pop_var_labels

## 5.2 · Individual-level RSF coefficients ####

# Maps ctmm's internal parameter names to readable labels
var_labels <- c(
  "LandUse.6_1 (1/LandUse.6_1)"    = "Ocean",
  "LandUse.5_1 (1/LandUse.5_1)"    = "Fresh Water",
  "LandUse.4_1 (1/LandUse.4_1)"    = "Urban Area",
  "LandUse.3_1 (1/LandUse.3_1)"    = "Low Vegetation",
  "LandUse.2_1 (1/LandUse.2_1)"    = "Lava Field",
  "RoadProximity (1/RoadProximity)" = "Distance to Road",
  "Ruggedness (1/Ruggedness)"       = "Ruggedness Index",
  "Slope (1/Slope)"                 = "Slope Angle",
  "Elevation (1/Elevation)"         = "Elevation",
  "TreeCover (1/TreeCover)"         = "Tree Cover",
  "NDVI (1/NDVI)"                   = "Vegetation Index (NDVI)"
)

result_ind <- do.call(rbind, lapply(rsf_fit, function(fit) {
  df <- as.data.frame(summary(fit)$CI)
  tibble::rownames_to_column(df[names(var_labels), ], var = "variable")
}))
result_ind$animal_id <- rep(names(x), each = 11)
result_ind$variable  <- unname(var_labels[result_ind$variable])

## 5.3 · Combined data frame (individuals + population) ####

name_map <- c(
  "Ocean"          = "Land Type-Ocean",
  "Fresh Water"    = "Land Type-Fresh Water",
  "Urban Area"     = "Land Type-Urban",
  "Low Vegetation" = "Land Type-Low vegetation",
  "Lava Field"     = "Land Type-Lava"
)

result_ind <- result_ind |>
  dplyr::mutate(term = dplyr::recode(variable, !!!name_map, .default = variable))

# CI caps prevent extreme individual estimates from distorting axis scales
# (cap = population estimate ± 5 × population CI width)
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

# Shared significance colour palette (Okabe-Ito)
sig_colors <- c(
  "Positive"        = "#56B4E9",
  "Negative"        = "#D55E00",
  "Non-significant" = "gray80"
)
src_shapes <- c("Individual" = 16, "Population" = 18)


###─────────────────────────────────────────────────────────────────────────###
## Fig. 1 · RSF forest plot: 3 × 3 individual + population grid ####
###─────────────────────────────────────────────────────────────────────────###

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
      width = 0.2
    ) +
    geom_point(
      data = \(d) dplyr::filter(d, source == "Population"),
      aes(color = significance, shape = source),
      size = 4
    ) +
    geom_errorbarh(
      data = \(d) dplyr::filter(d, source == "Population"),
      aes(xmin = low, xmax = high, color = significance),
      linewidth = 1.2, width = 0.3
    ) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
    scale_color_manual(
      values = sig_colors,
      limits = c("Positive", "Negative", "Non-significant"),
      drop   = FALSE,
      guide  = "none"
    ) +
    scale_shape_manual(values = src_shapes, guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 0.35, `FALSE` = 1), guide = "none") +
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

ggsave("./outputs/figures/Fig1_RSF_3x3.pdf",  fig1, width = 10, height = 10)
ggsave("./outputs/figures/Fig1_RSF_3x3.tiff", fig1, width = 10, height = 10,
        dpi = 600, compression = "lzw")


###─────────────────────────────────────────────────────────────────────────###
## Fig. 2 · Range distribution meta-analysis (A) + AKDE+RSF spatial map (B) ####
###─────────────────────────────────────────────────────────────────────────###

# Individual colour palette (shared between Panels A and B)
# Draws from Set1 + Dark2; red removed and 5th colour muted
pal_ind    <- c(RColorBrewer::brewer.pal(9, "Set1"), RColorBrewer::brewer.pal(8, "Dark2"))
pal_ind    <- pal_ind[-3]
pal_ind[5] <- "#888866"
ind_colors <- setNames(pal_ind[seq_along(akde_rsf_fit)], names(akde_rsf_fit))

# UD rasters clipped to 95% HR boundary (ctmm stores UD as p-values: 0 = core)
akde_ud_list <- lapply(akde_rsf_fit, function(ud) {
  r           <- raster::raster(ud)
  r[r > 0.95] <- NA
  terra::rast(r)
})

# Galápagos context inset: uniform grey islands, orange bbox marks Floreana
flor_bb <- sf::st_bbox(galapagos_flor)
pad     <- 15000  # 15 km padding in UTM 15S metres
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
  geom_sf(data = galapagos_islands, fill = "gray75", color = "gray40", linewidth = 0.2) +
  geom_sf(data = floreana_bbox_sf,  fill = NA, color = "#D55E00", linewidth = 0.8) +
  theme_void() +
  theme(panel.border = element_rect(color = "gray30", fill = NA, linewidth = 0.5))

# ggnewscale layers are stateful and must be rebuilt fresh for each ggplot() call
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
  akde_ud_list, ind_colors[names(akde_rsf_fit)],
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
  inset_element(fig2b_inset, left = 0.65, bottom = 0.65, right = 1.00, top = 1.00)

fig2 <- cowplot::plot_grid(
  ~ ctmm::meta(akde_rsf_fit, col = c(ind_colors, "black"),
               verbose = FALSE, sort = TRUE, labels = TRUE),
  fig2b,
  ncol           = 1,
  rel_heights    = c(1, 1.25),
  labels         = c("A", "B"),
  label_size     = 12,
  label_fontface = "bold"
)

fig2

ggsave("./outputs/figures/Fig2_Range_Map.pdf",  fig2, width = 8, height = 14)
ggsave("./outputs/figures/Fig2_Range_Map.tiff", fig2, width = 8, height = 14,
        dpi = 600, compression = "lzw")


###─────────────────────────────────────────────────────────────────────────###
## S.Fig. 1 · Population-level RSF: all variables, sorted by estimate ####
###─────────────────────────────────────────────────────────────────────────###

sfig1_df <- ci_df |>
  dplyr::filter(term != "Land Type-Ocean") |>
  dplyr::mutate(
    significance = dplyr::case_when(
      low > 0  ~ "Positive",
      high < 0 ~ "Negative",
      TRUE     ~ "Non-significant"
    )
  )

sfig1 <- ggplot(sfig1_df, aes(x = estimate, y = reorder(term, estimate), color = significance)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray60") +
  geom_errorbarh(aes(xmin = low, xmax = high), width = 0.25, linewidth = 0.8) +
  geom_point(shape = 18, size = 4) +
  scale_color_manual(
    values = sig_colors,
    limits = c("Positive", "Negative", "Non-significant"),
    name   = "Effect",
    drop   = FALSE
  ) +
  labs(x = "RSF Coefficient Estimate", y = NULL) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position = "right",
    legend.title    = element_text(size = 10, face = "bold")
  )

sfig1

ggsave("./outputs/figures/SFig1_PopRSF.pdf",  sfig1, width = 7, height = 5)
ggsave("./outputs/figures/SFig1_PopRSF.tiff", sfig1, width = 7, height = 5,
        dpi = 600, compression = "lzw")


###─────────────────────────────────────────────────────────────────────────###
## S.Fig. 2 · Land use + road + AKDE overlap map ####
###─────────────────────────────────────────────────────────────────────────###

lu_terra <- LandUse

landuse_cols <- c(
  "Forest"         = "#009E73",
  "Lava"           = "#999999",
  "Low vegetation" = "#F0E442",
  "Urban"          = "#E69F00",
  "Lagoon"         = "#56B4E9",
  "Ocean"          = "#0072B2"
)

# Binary 95% HR rasters (1 = inside UD, NA = outside)
akde_bin_sfig2 <- lapply(akde_rsf_fit, function(ud) {
  r     <- raster::raster(ud)
  r_bin <- r <= 0.95
  r_bin[!raster::values(r_bin)] <- NA
  terra::rast(r_bin)
})

# Per-individual dark-grey overlay; α = 0.25 compounds to reveal overlap intensity
overlap_layers_sfig2 <- unlist(lapply(akde_bin_sfig2, function(r) {
  list(
    new_scale_fill(),
    geom_spatraster(data = r),
    scale_fill_gradient(
      low      = scales::alpha("#1A1A1A", 0),
      high     = scales::alpha("#1A1A1A", 0.25),
      na.value = NA,
      guide    = "none"
    )
  )
}), recursive = FALSE)

sfig2_main <- ggplot() +
  geom_spatraster(data = lu_terra, aes(fill = class)) +
  scale_fill_manual(values = landuse_cols, name = "Land Use", na.value = NA) +
  overlap_layers_sfig2 +
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
    panel.background = element_rect(fill = "#0072B2", color = NA)
  )

sfig2 <- sfig2_main +
  inset_element(fig2b_inset, left = 0.65, bottom = 0.65, right = 1.00, top = 1.00)

sfig2

ggsave("./outputs/figures/SFig2_LandUse_AKDE.pdf",  sfig2, width = 8, height = 7)
ggsave("./outputs/figures/SFig2_LandUse_AKDE.tiff", sfig2, width = 8, height = 7,
        dpi = 600, compression = "lzw")
