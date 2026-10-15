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
prev_min  <- 0.1          # minimum prevalence
prev_max  <- 0.4         # upper limit so that Buffer (out) keeps enough area
max_tries <- 50            # max. redraws per species
seed_base <- 2026
mean_q    <- c(0.10, 0.90) # optimum between these percentiles
sd_rel    <- c(0.05, 0.60) # sd as a fraction of the variable range

# ======================= Environmental data (masked to Romania) ===========

set.seed(seed_base)

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

q_sd <- apply(smp, 2, quantile, probs = c(0.025, 0.975))  # core range for sd

# ================================ Functions ===============================

draw_pars <- function() {
  purrr::map_dfr(vars, \(v) tibble::tibble(
    variable = v,
    mean     = runif(1, q[1, v], q[2, v]),
    sd       = runif(1, sd_rel[1], sd_rel[2]) * (q_sd[2, v] - q_sd[1, v])
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
      virtualspecies::convertToPA(sp, PA.method = "probability",
                                  prob.method = "linear", a = 1, b = 0,
                                  plot = FALSE)
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
    (q_sd[2, true_pars$variable] - q_sd[1, true_pars$variable])
  
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



# =================== Variable selection: correlation figure ===============
# Spearman correlations of all 19 bioclim variables within Romania (same
# sample as used for the variable selection), ordered by hierarchical
# clustering (1 - |r|). Bold: |r| >= 0.7. Outlined: pairs of selected vars.

bio_all <- terra::rast(file.path(envrmt$path_raster,
                                 "climate/wc2.1_country/ROU_wc2.1_30s_bio.tif"))
names(bio_all) <- substr(names(bio_all), 11, 16)
rou     <- geodata::gadm(country = "ROU", level = 0, path = envrmt$path_raster)
bio_all <- bio_all |> terra::crop(rou) |> terra::mask(rou)

set.seed(seed_base)
smp_all <- terra::spatSample(bio_all, 20000, na.rm = TRUE)
cm_all  <- cor(smp_all, method = "spearman")

# order variables by clustering so that correlated groups sit together
hc  <- hclust(as.dist(1 - abs(cm_all)), method = "complete")
ord <- hc$labels[hc$order]

# statistics of the selected set
cm_sel  <- cm_all[vars, vars]
max_r   <- max(abs(cm_sel[upper.tri(cm_sel)]))
max_vif <- max(diag(solve(cm_sel)))

bio_names <- c(
  bio_1  = "Annual mean temp.",          bio_2  = "Mean diurnal range",
  bio_3  = "Isothermality",              bio_4  = "Temp. seasonality",
  bio_5  = "Max temp. warmest month",    bio_6  = "Min temp. coldest month",
  bio_7  = "Temp. annual range",         bio_8  = "Mean temp. wettest qtr.",
  bio_9  = "Mean temp. driest qtr.",     bio_10 = "Mean temp. warmest qtr.",
  bio_11 = "Mean temp. coldest qtr.",    bio_12 = "Annual precipitation",
  bio_13 = "Prec. wettest month",        bio_14 = "Prec. driest month",
  bio_15 = "Prec. seasonality",          bio_16 = "Prec. wettest qtr.",
  bio_17 = "Prec. driest qtr.",          bio_18 = "Prec. warmest qtr.",
  bio_19 = "Prec. coldest qtr."
)
axis_lab <- setNames(
  paste0(toupper(sub("_", "", names(bio_names))), " ", bio_names,
         ifelse(names(bio_names) %in% vars, " *", "")),
  names(bio_names))

cm_df <- as.data.frame(as.table(cm_all)) |>
  setNames(c("var1", "var2", "r")) |>
  dplyr::mutate(
    i    = match(var1, ord),
    j    = match(var2, ord),
    var1 = factor(var1, levels = ord),
    var2 = factor(var2, levels = rev(ord)),
    high = abs(r) >= 0.7,
    sel  = var1 %in% vars & var2 %in% vars
  ) |>
  dplyr::filter(j > i)   # lower triangle without diagonal

p_cor <- ggplot(cm_df, aes(var1, var2, fill = r)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = sprintf("%.2f", r),
                fontface = ifelse(high, "bold", "plain")), size = 2.1) +
  geom_tile(data = dplyr::filter(cm_df, sel),
            fill = NA, colour = "black", linewidth = 0.8) +
  scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                       limits = c(-1, 1), name = "Spearman r") +
  scale_x_discrete(labels = axis_lab) +
  scale_y_discrete(labels = axis_lab) +
  coord_fixed() +
  labs(x = NULL, y = NULL,
       caption = sprintf(paste0("Bold: absolute r \u2265 0.7. Outlined / *: selected ",
                                "variables; max. absolute r = %.2f, max. VIF = %.1f."),
                         max_r, max_vif)) +
  theme_minimal(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank(),
        legend.position = "inside",
        legend.position.inside = c(0.85, 0.8))

ggsave(file.path(envrmt$path_figures, "01_bioclim_correlation.png"),
       p_cor, width = 8.5, height = 8, dpi = 300)

# matrix for the supplement
res_dir <- file.path(envrmt$path_data, "results")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(round(cm_all[ord, ord], 3),
          file.path(res_dir, "bioclim_correlation_spearman.csv"))



path <- file.path(rootDir, "src", "01b_species_panel.R")
source(path, echo = TRUE)