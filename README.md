<<<<<<< HEAD
Short-eared Owl RSF – Floreana Island, Galápagos
================
Josue Arteaga-Torres
2026-05-22

## Overview

This project performs a **Resource Selection Function (RSF) analysis**
of the Galápagos Short-eared Owl (*Asio flammeus galapagoensis*) on
Floreana Island, Galápagos, Ecuador. Movement data from 16 GPS-tracked
individuals are combined with environmental raster layers to model
habitat selection using continuous-time movement modelling (ctmm).

The analysis estimates which environmental features (vegetation, land
use, topography, road proximity) owls preferentially use or avoid
relative to available habitat, accounting for autocorrelation in
movement trajectories.

------------------------------------------------------------------------

## Author

Josue Arteaga-Torres —

**AI disclaimer:** Code annotations and debugging were assisted by the
Posit Assistant tool (v0.3.2) in RStudio 2026.04.0+526 “Globemaster
Allium” (2026-04-18) for Windows.

------------------------------------------------------------------------

<a href="https://github.com/JosueDATYoshi/seow_floreana_rsf-akde/tree/main">Short-eared
Owl RSF – Floreana Island, Galápagos</a> © 2026 by
<a href="https://example.com">Josue D. Arteaga-Torres</a> is licensed
under <a href="https://creativecommons.org/licenses/by/4.0/">CC BY
4.0</a><img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;">

## Analysis Workflow

    Telemetry data (16 individuals)
            │
            ▼
    SEC. 4.1  Movement model selection (ctmm.select, AICc)
            │
            ▼
    SEC. 4.2  Autocorrelated KDE (AKDE) — home range estimation
            │
            ▼
    SEC. 4.3  Resource Selection Function (rsf.fit, Monte Carlo integrator)
            │   Forest (LandUse class 1) as reference category
            ▼
    SEC. 4.4  RSF-informed AKDE (AKDE + RSF) — habitat-weighted utilisation distribution

------------------------------------------------------------------------

## Data Inputs

| Input | Description |
|----|----|
| `inputs/SEO_owls_data.RData` | Pre-filtered telemetry data (`owls_tele`, `owls_move`). Two individuals (PS0024, PS0028) excluded prior to this script. |
| `inputs/SEO_rasterlist_NDVI_tree_elev_slope_rug_urb_farm_road.RData` | List of continuous environmental rasters (`owls_env_list`): NDVI, TreeCover, Elevation, Slope, Ruggedness, Urban, FarmLand, RoadProximity. Z-score before passing to `rsf.fit()`. |
| `inputs/floreana_ecosystems.gpkg` | Ecosystem polygon layer; rasterized at runtime to produce the `LandUse` categorical raster (Forest, Lava, Low vegetation, Urban, Lagoon). |
| `inputs/SEO_spatial_layers.gpkg` | GeoPackage with three layers: `floreana_border` (island outline, UTM 15S), `road` (main road centreline), `galapagos` (all archipelago islands). |

Pre-computed model objects (`.RData`) are loaded directly to avoid
re-running computationally expensive fits:

| Object | File |
|----|----|
| Movement models (`ctmm_fit`) | `outputs/ctmm_select_Final.RData` |
| AKDE home ranges (`akde_fit`) | `outputs/AKDE_Final.RData` |
| RSF fits (`rsf_fit`) | `outputs/RSF_Final.RData` |
| RSF-informed AKDE (`akde_rsf_fit`) | `outputs/AKDERSF_Final.RData` |
| Population-mean RSF (`mean_rsf`) | `outputs/Mean_Population_RSF_Final.RData` |
| `outputs/RSF_results_final_long.xlsx` | Long-format table of harmonised RSF variable labels used in figures. |

------------------------------------------------------------------------

## Environmental Variables

| Variable | Type | Notes |
|----|----|----|
| Land Use | Categorical | Forest (reference), Low Vegetation, Urban, Lava, Lagoon, Ocean |
| Distance to Road | Continuous | Z-scored |
| Ruggedness Index | Continuous | Z-scored |
| Slope Angle | Continuous | Z-scored |
| Elevation | Continuous | Z-scored |
| Tree Cover | Continuous | Z-scored |
| Vegetation Index (NDVI) | Continuous | Z-scored |

------------------------------------------------------------------------

## Figures

| Figure | Description |
|----|----|
| **Fig. 1** | 3×3 RSF forest plot — individual estimates (circles) and population mean (diamonds) for 9 environmental variables. Points are colour-coded by significance direction (blue = positive selection, red = negative, grey = non-significant). |
| **Fig. 2A** | Home range size meta-analysis (`ctmm::meta()`), ranked by 95% AKDE+RSF area with population mean overlaid. |
| **Fig. 2B** | Spatial map of per-individual AKDE+RSF utilisation distributions (95% contour) overlaid on Floreana island outline, with Galápagos archipelago context inset. |
| **S.Fig. 1** | Population-level RSF forest plot for all variables, sorted by selection coefficient (supplementary material). |
| **S.Fig. 2** | Land-use base map overlaid with stacked 95% AKDE+RSF home ranges (opacity compounds to show spatial overlap among individuals). |

Figures are saved to `outputs/figures/` in both PDF and TIFF (600 dpi,
LZW) formats.

------------------------------------------------------------------------

## Key Packages

| Package | Role |
|----|----|
| `ctmm` | Continuous-time movement modelling, AKDE, RSF |
| `future.apply` / `parallelly` | Parallel model fitting |
| `sf`, `terra`, `raster` | Spatial vector and raster handling |
| `tidyterra`, `ggnewscale`, `ggspatial` | ggplot2 extensions for spatial data |

`ggplot2`, `patchwork`, `cowplot` \| Figure composition \|

> `ctmm` is under active development. If the CRAN version is outdated,
> install from GitHub: `remotes::install_github("ctmm-initiative/ctmm")`

------------------------------------------------------------------------

## Reproducibility

A fixed random seed (`set.seed(3847)`) ensures the Monte Carlo
integrator in `rsf.fit()` produces identical results across runs. Model
fitting is parallelised across up to 6 CPU cores via
`future::plan("multisession")`.

------------------------------------------------------------------------

## Session Info

    R version 4.5.3 RC (2026-03-03 r89528 ucrt)
    Platform: x86_64-w64-mingw32/x64
    Running under: Windows 11 x64 (build 26200)

    Matrix products: default
      LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
    [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
    [5] LC_TIME=English_United States.utf8    

    time zone: America/Denver
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] here_1.0.2          RColorBrewer_1.1-3  cowplot_1.2.0       patchwork_1.3.2    
     [5] ggplot2_4.0.2       ggspatial_1.1.10    scales_1.4.0        ggnewscale_0.5.2   
     [9] tidyterra_1.1.0     raster_3.6-32       sp_2.2-1            terra_1.8-93       
    [13] sf_1.1-0            readxl_1.4.5        tibble_3.3.1        dplyr_1.2.0        
    [17] parallelly_1.46.1   future.apply_1.20.2 future_1.69.0       ctmm_1.3.1         

    loaded via a namespace (and not attached):
     [1] gtable_0.3.6        xfun_0.57           lattice_0.22-9      numDeriv_2016.8-1.1
     [5] vctrs_0.7.1         tools_4.5.3         generics_0.1.4      parallel_4.5.3     
     [9] proxy_0.4-29        pkgconfig_2.0.3     Matrix_1.7-4        KernSmooth_2.23-26 
    [13] data.table_1.18.4   S7_0.2.1            lifecycle_1.0.5     compiler_4.5.3     
    [17] farver_2.1.2        textshaping_1.0.5   statmod_1.5.1       codetools_0.2-20   
    [21] Bessel_0.7-0        htmltools_0.5.9     class_7.3-23        gmp_0.7-5.1        
    [25] yaml_2.3.12         pillar_1.11.1       tidyr_1.3.2         classInt_0.4-11    
    [29] tidyselect_1.2.1    digest_0.6.39       purrr_1.2.1         listenv_0.10.1     
    [33] labeling_0.4.3      rprojroot_2.1.1     fastmap_1.2.0       grid_4.5.3         
    [37] expm_1.0-0          cli_3.6.5           magrittr_2.0.4      e1071_1.7-17       
    [41] Rmpfr_1.1-2         withr_3.0.2         rmarkdown_2.31      globals_0.19.1     
    [45] otel_0.2.0          cellranger_1.1.0    ragg_1.5.2          evaluate_1.0.5     
    [49] knitr_1.51          gridGraphics_0.5-1  rlang_1.1.7         Rcpp_1.1.1         
    [53] glue_1.8.0          DBI_1.3.0           rstudioapi_0.18.0   R6_2.6.1           
    [57] systemfonts_1.3.2   units_1.0-1  
=======
Short-eared Owl RSF – Floreana Island, Galápagos
================
Josue Arteaga-Torres
2026-05-22

## Overview

This project performs a **Resource Selection Function (RSF) analysis**
of the Galápagos Short-eared Owl (*Asio flammeus galapagoensis*) on
Floreana Island, Galápagos, Ecuador. Movement data from 16 GPS-tracked
individuals are combined with environmental raster layers to model
habitat selection using continuous-time movement modelling (ctmm).

The analysis estimates which environmental features (vegetation, land
use, topography, road proximity) owls preferentially use or avoid
relative to available habitat, accounting for autocorrelation in
movement trajectories.

------------------------------------------------------------------------

## Author

Josue Arteaga-Torres —

**AI disclaimer:** Code annotations and debugging were assisted by the
Posit Assistant tool (v0.3.2) in RStudio 2026.04.0+526 “Globemaster
Allium” (2026-04-18) for Windows.

------------------------------------------------------------------------

<a href="https://github.com/JosueDATYoshi/seow_floreana_rsf-akde/tree/main">Short-eared
Owl RSF – Floreana Island, Galápagos</a> © 2026 by
<a href="https://example.com">Josue D. Arteaga-Torres</a> is licensed
under <a href="https://creativecommons.org/licenses/by/4.0/">CC BY
4.0</a><img src="https://mirrors.creativecommons.org/presskit/icons/cc.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;"><img src="https://mirrors.creativecommons.org/presskit/icons/by.svg" alt="" style="max-width: 1em;max-height:1em;margin-left: .2em;">

## Analysis Workflow

    Telemetry data (16 individuals)
            │
            ▼
    SEC. 4.1  Movement model selection (ctmm.select, AICc)
            │
            ▼
    SEC. 4.2  Autocorrelated KDE (AKDE) — home range estimation
            │
            ▼
    SEC. 4.3  Resource Selection Function (rsf.fit, Monte Carlo integrator)
            │   Forest (LandUse class 1) as reference category
            ▼
    SEC. 4.4  RSF-informed AKDE (AKDE + RSF) — habitat-weighted utilisation distribution

------------------------------------------------------------------------

## Data Inputs

| Input | Description |
|----|----|
| `inputs/SEO_owls_data.RData` | Pre-filtered telemetry data (`owls_tele`, `owls_move`). Two individuals (PS0024, PS0028) excluded prior to this script. |
| `inputs/SEO_rasterlist_NDVI_tree_elev_slope_rug_urb_farm_road.RData` | List of continuous environmental rasters (`owls_env_list`): NDVI, TreeCover, Elevation, Slope, Ruggedness, Urban, FarmLand, RoadProximity. Z-score before passing to `rsf.fit()`. |
| `inputs/floreana_ecosystems.gpkg` | Ecosystem polygon layer; rasterized at runtime to produce the `LandUse` categorical raster (Forest, Lava, Low vegetation, Urban, Lagoon). |
| `inputs/SEO_spatial_layers.gpkg` | GeoPackage with three layers: `floreana_border` (island outline, UTM 15S), `road` (main road centreline), `galapagos` (all archipelago islands). |

Pre-computed model objects (`.RData`) are loaded directly to avoid
re-running computationally expensive fits:

| Object | File |
|----|----|
| Movement models (`ctmm_fit`) | `outputs/ctmm_select_Final.RData` |
| AKDE home ranges (`akde_fit`) | `outputs/AKDE_Final.RData` |
| RSF fits (`rsf_fit`) | `outputs/RSF_Final.RData` |
| RSF-informed AKDE (`akde_rsf_fit`) | `outputs/AKDERSF_Final.RData` |
| Population-mean RSF (`mean_rsf`) | `outputs/Mean_Population_RSF_Final.RData` |
| `outputs/RSF_results_final_long.xlsx` | Long-format table of harmonised RSF variable labels used in figures. |

------------------------------------------------------------------------

## Environmental Variables

| Variable | Type | Notes |
|----|----|----|
| Land Use | Categorical | Forest (reference), Low Vegetation, Urban, Lava, Lagoon, Ocean |
| Distance to Road | Continuous | Z-scored |
| Ruggedness Index | Continuous | Z-scored |
| Slope Angle | Continuous | Z-scored |
| Elevation | Continuous | Z-scored |
| Tree Cover | Continuous | Z-scored |
| Vegetation Index (NDVI) | Continuous | Z-scored |

------------------------------------------------------------------------

## Figures

| Figure | Description |
|----|----|
| **Fig. 1** | 3×3 RSF forest plot — individual estimates (circles) and population mean (diamonds) for 9 environmental variables. Points are colour-coded by significance direction (blue = positive selection, red = negative, grey = non-significant). |
| **Fig. 2A** | Home range size meta-analysis (`ctmm::meta()`), ranked by 95% AKDE+RSF area with population mean overlaid. |
| **Fig. 2B** | Spatial map of per-individual AKDE+RSF utilisation distributions (95% contour) overlaid on Floreana island outline, with Galápagos archipelago context inset. |
| **S.Fig. 1** | Population-level RSF forest plot for all variables, sorted by selection coefficient (supplementary material). |
| **S.Fig. 2** | Land-use base map overlaid with stacked 95% AKDE+RSF home ranges (opacity compounds to show spatial overlap among individuals). |

Figures are saved to `outputs/figures/` in both PDF and TIFF (600 dpi,
LZW) formats.

------------------------------------------------------------------------

## Key Packages

| Package | Role |
|----|----|
| `ctmm` | Continuous-time movement modelling, AKDE, RSF |
| `future.apply` / `parallelly` | Parallel model fitting |
| `sf`, `terra`, `raster` | Spatial vector and raster handling |
| `tidyterra`, `ggnewscale`, `ggspatial` | ggplot2 extensions for spatial data |

`ggplot2`, `patchwork`, `cowplot` \| Figure composition \|

> `ctmm` is under active development. If the CRAN version is outdated,
> install from GitHub: `remotes::install_github("ctmm-initiative/ctmm")`

------------------------------------------------------------------------

## Reproducibility

A fixed random seed (`set.seed(3847)`) ensures the Monte Carlo
integrator in `rsf.fit()` produces identical results across runs. Model
fitting is parallelised across up to 6 CPU cores via
`future::plan("multisession")`.

------------------------------------------------------------------------

## Session Info

    R version 4.5.3 RC (2026-03-03 r89528 ucrt)
    Platform: x86_64-w64-mingw32/x64
    Running under: Windows 11 x64 (build 26200)

    Matrix products: default
      LAPACK version 3.12.1

    locale:
    [1] LC_COLLATE=English_United States.utf8  LC_CTYPE=English_United States.utf8   
    [3] LC_MONETARY=English_United States.utf8 LC_NUMERIC=C                          
    [5] LC_TIME=English_United States.utf8    

    time zone: America/Denver
    tzcode source: internal

    attached base packages:
    [1] stats     graphics  grDevices utils     datasets  methods   base     

    other attached packages:
     [1] here_1.0.2          RColorBrewer_1.1-3  cowplot_1.2.0       patchwork_1.3.2    
     [5] ggplot2_4.0.2       ggspatial_1.1.10    scales_1.4.0        ggnewscale_0.5.2   
     [9] tidyterra_1.1.0     raster_3.6-32       sp_2.2-1            terra_1.8-93       
    [13] sf_1.1-0            readxl_1.4.5        tibble_3.3.1        dplyr_1.2.0        
    [17] parallelly_1.46.1   future.apply_1.20.2 future_1.69.0       ctmm_1.3.1         

    loaded via a namespace (and not attached):
     [1] gtable_0.3.6        xfun_0.57           lattice_0.22-9      numDeriv_2016.8-1.1
     [5] vctrs_0.7.1         tools_4.5.3         generics_0.1.4      parallel_4.5.3     
     [9] proxy_0.4-29        pkgconfig_2.0.3     Matrix_1.7-4        KernSmooth_2.23-26 
    [13] data.table_1.18.4   S7_0.2.1            lifecycle_1.0.5     compiler_4.5.3     
    [17] farver_2.1.2        textshaping_1.0.5   statmod_1.5.1       codetools_0.2-20   
    [21] Bessel_0.7-0        htmltools_0.5.9     class_7.3-23        gmp_0.7-5.1        
    [25] yaml_2.3.12         pillar_1.11.1       tidyr_1.3.2         classInt_0.4-11    
    [29] tidyselect_1.2.1    digest_0.6.39       purrr_1.2.1         listenv_0.10.1     
    [33] labeling_0.4.3      rprojroot_2.1.1     fastmap_1.2.0       grid_4.5.3         
    [37] expm_1.0-0          cli_3.6.5           magrittr_2.0.4      e1071_1.7-17       
    [41] Rmpfr_1.1-2         withr_3.0.2         rmarkdown_2.31      globals_0.19.1     
    [45] otel_0.2.0          cellranger_1.1.0    ragg_1.5.2          evaluate_1.0.5     
    [49] knitr_1.51          gridGraphics_0.5-1  rlang_1.1.7         Rcpp_1.1.1         
    [53] glue_1.8.0          DBI_1.3.0           rstudioapi_0.18.0   R6_2.6.1           
    [57] systemfonts_1.3.2   units_1.0-1  
>>>>>>> 290c6a5fd58aec4da1c097ed1d56f05862c8a65c
