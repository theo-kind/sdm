# ================================= Notes =================================
#
# 01b_species_panel.R
# Overview figure of the virtual species, sized for an A4 page:
# true response curves (one column per variable in `vars`) and the
# suitability map per species (rows), with prevalence in the row label.
# Reads the outputs of 01_generate_species.R.
#
# ==========================================================================

## ================================= Load set up =================================
# (rootDir / envimaR setup, package loading, and function sourcing
#  are handled centrally by src/00_setup.R — source that first)
library(envimaR)

rootDir <- paste0(getwd(),"/project_folder")

rootDir <- envimaR::alternativeEnvi(
  root_folder = rootDir,
  alt_env_id = "COMPUTERNAME",
  alt_env_value = "PCRZP",
  alt_env_root_folder = "F:/BEN/edu"
)

path <- file.path(rootDir, "src", "00_setup.R")
source(path, echo = TRUE)


# ================================ Parameters ==============================

vars       <- c("bio_3", "bio_5", "bio_12", "bio_15")
var_labels <- c(bio_3  = "BIO3\nIsothermality (%)",
                bio_5  = "BIO5\nMax. temperature (\u00b0C)",
                bio_12 = "BIO12\nAnnual precip. (mm)",
                bio_15 = "BIO15\nPrecip. season. (CV)")
line_col   <- "#3B4A7A"
map_fact   <- 6          # aggregation of the maps for display only
fig_size   <- c(w = 6.3, h = 9.2)   # inches; ~16 x 23 cm (A4 text area)

env      <- terra::rast(file.path(envrmt$path_raster, "env_rou_30s.tif"))
rng      <- terra::minmax(env)
overview <- readRDS(file.path(envrmt$path_virtual_species, "species_overview.RDS"))
sp_lab   <- setNames(sprintf("%s\nPrev. %.2f", overview$id, overview$prevalence),
                     overview$id)


# ================================ Data ====================================

# true responses: Gaussian, min-max scaled over the range in Romania
curves <- purrr::map_dfr(overview$id, \(sp) {
  pars <- readRDS(file.path(envrmt$path_virtual_species,
                            paste0(sp, "_meta.RDS")))$true_pars
  purrr::map_dfr(vars, \(v) {
    p <- pars[pars$variable == v, ]
    x <- seq(rng[1, v], rng[2, v], length.out = 300)
    y <- dnorm(x, p$mean, p$sd)
    tibble::tibble(species = sp, variable = v, x = x,
                   y = (y - min(y)) / (max(y) - min(y)))
  })
}) |>
  dplyr::mutate(species  = factor(species, levels = overview$id),
                variable = factor(variable, levels = vars))

# suitability maps (aggregated for display)
maps <- purrr::map_dfr(overview$id, \(sp) {
  terra::rast(file.path(envrmt$path_virtual_species, paste0(sp, ".tif")))[["suitab"]] |>
    terra::aggregate(fact = map_fact, fun = "mean", na.rm = TRUE) |>
    terra::as.data.frame(xy = TRUE, na.rm = TRUE) |>
    stats::setNames(c("x", "y", "suitab")) |>
    dplyr::mutate(species = sp)
}) |>
  dplyr::mutate(species = factor(species, levels = overview$id),
                panel = "Suitability")


# ================================= Figure =================================

theme_panel <- theme_bw(base_size = 7) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", size = 6.5),
        panel.spacing.y = unit(2, "pt"),
        plot.margin = margin(2, 2, 2, 2))

p_curves <- ggplot(curves, aes(x, y)) +
  geom_line(colour = line_col, linewidth = 0.5) +
  facet_grid(species ~ variable, scales = "free_x",
             labeller = labeller(variable = var_labels)) +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1)) +
  scale_x_continuous(n.breaks = 4) +
  labs(x = NULL, y = "Relative suitability") +
  theme_panel +
  theme(strip.text.y = element_blank(),
        axis.text = element_text(size = 5.5))

p_maps <- ggplot(maps, aes(x, y, fill = suitab)) +
  geom_raster() +
  coord_quickmap(expand = FALSE) +
  facet_grid(species ~ panel, labeller = labeller(species = sp_lab)) +
  scale_fill_viridis_c(limits = c(0, 1), breaks = c(0, 0.5, 1),
                       name = "Suitability") +
  theme_panel +
  theme(axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title = element_blank(), panel.border = element_blank(),
        strip.text.y = element_text(angle = 0, hjust = 0, size = 6),
        legend.position = "bottom",
        legend.key.height = unit(0.2, "cm"), legend.key.width = unit(0.6, "cm"),
        legend.title = element_text(size = 6, vjust = 0.9),
        legend.text = element_text(size = 5.5),
        legend.margin = margin(0, 0, 0, 0))

p <- patchwork::wrap_plots(p_curves, p_maps, nrow = 1, widths = c(4, 1.35))

ggsave(file.path(envrmt$path_figures, "01_species_summary_panel.png"),
       p, width = fig_size[["w"]], height = fig_size[["h"]], dpi = 400)
