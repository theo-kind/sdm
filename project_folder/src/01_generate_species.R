# ================================= Notes =================================
#
# 01_generate_species.R
# Generates n_species random virtual species (Gaussian responses,
# multiplicative) over the bioclim variables in `vars` for Romania
# (WorldClim 2.1, 30").
# Parameters are drawn randomly within controlled bounds:
#   - optimum between the 10th and 90th percentile of the variable
#   - sd = 10–40 % of the variable range
# Output per species: GeoTIFF (suitab + pa) + metadata RDS (true
# parameters for the later response curve comparison).
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

vars      <- c("bio_3", "bio_5", "bio_12", "bio_15")
n_species <- 12
prev_min  <- 0.05          # minimum prevalence
prev_max  <- 0.50          # upper limit so that Buffer (out) keeps enough area
max_tries <- 50            # max. redraws per species
seed_base <- 2026
mean_q    <- c(0.10, 0.90) # optimum between these percentiles
sd_rel    <- c(0.10, 0.40) # sd as a fraction of the variable range

# ======================= Environmental data (masked to Romania) ===========

env_file <- file.path(envrmt$path_raster, "env_rou_30s.tif")

if (!file.exists(env_file)) {
  r <- terra::rast(file.path(envrmt$path_raster,
                             "climate/wc2.1_country/ROU_wc2.1_30s_bio.tif"))
  names(r) <- substr(names(r), 11, 16)
  rou <- geodata::gadm(country = "ROU", level = 0, path = envrmt$path_raster)
  env <- r[[vars]] |> terra::crop(rou) |> terra::mask(rou)
  terra::writeRaster(env, env_file, overwrite = TRUE)
}
env <- terra::rast(env_file)

# Check: cell count after masking (~400,000 expected) and collinearity
print(terra::global(!is.na(env[[1]]), "sum"))
smp <- terra::spatSample(env, 50000, na.rm = TRUE)
print(round(cor(smp), 2))

# Ranges and percentiles for drawing the parameters
rng <- terra::minmax(env)
q   <- apply(smp, 2, quantile, probs = mean_q)


# ================================ Functions ===============================

draw_pars <- function() {
  purrr::map_dfr(vars, \(v) tibble::tibble(
    variable = v,
    mean     = runif(1, q[1, v], q[2, v]),
    sd       = runif(1, sd_rel[1], sd_rel[2]) * (rng[2, v] - rng[1, v])
  ))
}


# ============================== Generate species ==========================

species_meta <- vector("list", n_species)
suit_list    <- vector("list", n_species)

for (i in seq_len(n_species)) {
  sp_id <- sprintf("sp%02d", i)
  set.seed(seed_base + i)
  
  prev <- NA_real_
  for (k in seq_len(max_tries)) {
    true_pars <- draw_pars()
    fun_list  <- purrr::pmap(true_pars, \(variable, mean, sd)
                             c(fun = "dnorm", mean = mean, sd = sd))
    names(fun_list) <- true_pars$variable
    
    vsp <- tryCatch({
      sp <- virtualspecies::generateSpFromFun(
        raster.stack = env,
        parameters   = do.call(virtualspecies::formatFunctions, fun_list),
        species.type = "multiplicative",
        plot         = FALSE
      )
      virtualspecies::convertToPA(sp, plot = FALSE)
    }, error = function(e) NULL)
    
    if (is.null(vsp)) next
    prev <- terra::global(vsp$pa.raster, "mean", na.rm = TRUE)[1, 1]
    if (prev >= prev_min && prev <= prev_max) break
  }
  if (is.na(prev) || prev < prev_min || prev > prev_max) {
    stop(sp_id, ": no species within the prevalence window after ", max_tries, " attempts")
  }
  
  # relative niche breadth per variable (covariate only, not a design factor)
  true_pars$rel_sd <- true_pars$sd /
    (rng[2, true_pars$variable] - rng[1, true_pars$variable])
  
  # save rasters as GeoTIFF (SpatRaster objects do not survive saveRDS)
  out <- c(vsp$suitab.raster, vsp$pa.raster)
  names(out) <- c("suitab", "pa")
  terra::writeRaster(out,
                     file.path(envrmt$path_virtual_species, paste0(sp_id, ".tif")),
                     overwrite = TRUE)
  
  # metadata incl. true response parameters
  meta <- list(id = sp_id, seed = seed_base + i, tries = k,
               prevalence = prev, true_pars = true_pars)
  saveRDS(meta, file.path(envrmt$path_virtual_species, paste0(sp_id, "_meta.RDS")))
  
  species_meta[[i]] <- meta
  suit_list[[i]]    <- vsp$suitab.raster
  message(sp_id, " done | prevalence = ", round(prev, 3), " | attempts = ", k)
}


# ================================= Overview ===============================

overview <- purrr::map_dfr(species_meta, \(m) tibble::tibble(
  id          = m$id,
  prevalence  = m$prevalence,
  mean_rel_sd = mean(m$true_pars$rel_sd),
  tries       = m$tries
))
print(overview)
saveRDS(overview, file.path(envrmt$path_virtual_species, "species_overview.RDS"))


# ================================= Figure =================================

suit_stack <- terra::rast(suit_list)
names(suit_stack) <- overview$id

png(file.path(envrmt$path_figures, "01_virtual_species_suitability.png"),
    width = 2400, height = 1400, res = 200)
terra::plot(suit_stack, nc = 3, range = c(0, 1))
dev.off()
