# ================================= Notes =================================
#
# 01b_species_panel.R
# Overview figure of the virtual species: true response curves
# (variables in `vars`) + suitability map per species.
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

vars     <- c("bio_3", "bio_5", "bio_12", "bio_15")
line_col <- "#3B4A7A"
map_col  <- hcl.colors(100, "viridis")

env      <- terra::rast(file.path(envrmt$path_raster, "env_rou_30s.tif"))
rng      <- terra::minmax(env)
overview <- readRDS(file.path(envrmt$path_virtual_species, "species_overview.RDS"))
n_sp     <- nrow(overview)


# ================================ Functions ===============================

# true response: Gaussian, min-max scaled over the range in Romania
true_curve <- function(pars, v, n = 500) {
  p <- pars[pars$variable == v, ]
  x <- seq(rng[1, v], rng[2, v], length.out = n)
  y <- dnorm(x, p$mean, p$sd)
  list(x = x, y = (y - min(y)) / (max(y) - min(y)))
}


# ================================= Figure =================================

png(file.path(envrmt$path_figures, "01_species_summary_panel.png"),
    width = 4000, height = 600 * n_sp, res = 250)

par(mfrow = c(n_sp, length(vars) + 1),
    oma = c(0.5, 3, 2, 0.5), mgp = c(2, 0.6, 0))

for (i in seq_len(n_sp)) {
  sp_id <- overview$id[i]
  meta  <- readRDS(file.path(envrmt$path_virtual_species, paste0(sp_id, "_meta.RDS")))
  
  # Response Curves
  for (j in seq_along(vars)) {
    rc <- true_curve(meta$true_pars, vars[j])
    par(mar = c(2.5, 2.8, 1.2, 0.6))
    plot(rc$x, rc$y, type = "l", lwd = 2, col = line_col,
         ylim = c(0, 1), xlab = "", ylab = "", las = 1, cex.axis = 0.9)
    if (i == 1) mtext(vars[j], side = 3, line = 1.2, font = 2)
  }
  
  # Suitability map
  suit <- terra::rast(file.path(envrmt$path_virtual_species,
                                paste0(sp_id, ".tif")))[["suitab"]]
  terra::plot(suit, col = map_col, range = c(0, 1),
              mar = c(2.5, 1.5, 1.2, 3.5),
              pax = list(cex.axis = 0.7), plg = list(cex = 0.7))
  if (i == 1) mtext("Suitability map", side = 3, line = 1.2, font = 2)
  mtext(sprintf("Prevalence = %.2f", overview$prevalence[i]),
        side = 1, line = 2, cex = 0.6)
  
  # Row label
  mtext(sp_id, side = 2, outer = TRUE, las = 1, font = 2, line = 0.5,
        at = 1 - (i - 0.5) / n_sp)
}

dev.off()