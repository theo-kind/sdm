# ================================= Notes =================================
#
# 02_run_models.R
# Draw one presence set per species × n_pres × replicate, then fit
# maxnet for 6 background strategies × 4 ratios and evaluate against
# the truth. Results are saved as one RDS per presence set.
# Resumable: existing result files are skipped.
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

test_mode <- FALSE   # TRUE: 1 species, n_pres = 100 & 500, 1 replicate

vars          <- c("bio_3", "bio_5", "bio_12", "bio_15")
n_pres_levels <- c(50, 100, 300, 500)
ratio_levels  <- c("2:1" = 0.5, "1:1" = 1, "1:10" = 10, "1:100" = 100)  # PA per presence
strategies    <- c("ecopa", "use", "rand_geo", "rand_env", "buf_in", "buf_out")
n_rep         <- 20      # replicates
buffer_width  <- 50000   # m
use_thres     <- 0.75    # KDE quantile (USE default)
grid_res      <- 20      # grid in PC space (USE, Random env)
ecopa_bins    <- 30      # bins per variable (EcoPA)
ecopa_in      <- 0.5     # share of PAs inside the niche (EcoPA)
n_eval        <- 50000   # fixed evaluation cells per species
seed_base     <- 2026

res_dir <- file.path(envrmt$path_data, "results", if (test_mode) "test" else "runs")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)


# ============================ Prepare data ================================

env       <- terra::rast(file.path(envrmt$path_raster, "env_rou_30s.tif"))
tmpl      <- env[[1]]
cell_tab  <- terra::as.data.frame(env[[vars]], cells = TRUE, na.rm = TRUE)
rng       <- sapply(cell_tab[vars], range)
env_means <- colMeans(cell_tab[vars])
q_trim <- sapply(cell_tab[vars], quantile, probs = c(0.025, 0.975))

# PCA once (basis for USE and Random env)
pca   <- prcomp(cell_tab[vars], scale. = TRUE)
print(summary(pca))
pc    <- as.data.frame(pca$x[, 1:2])
pc_sf <- sf::st_as_sf(pc, coords = c("PC1", "PC2"))

overview <- readRDS(file.path(envrmt$path_virtual_species, "species_overview.RDS"))

if (test_mode) {
  overview      <- overview[1, ]
  n_pres_levels <- c(100, 500)
  n_rep         <- 1
}

runs <- tidyr::expand_grid(strategy = strategies, ratio = names(ratio_levels))


# ================================ Main loop ===============================

for (sp_id in overview$id) {
  
  sp_num    <- as.integer(sub("sp", "", sp_id))
  meta      <- readRDS(file.path(envrmt$path_virtual_species, paste0(sp_id, "_meta.RDS")))
  sp_val    <- terra::rast(file.path(envrmt$path_virtual_species,
                                     paste0(sp_id, ".tif")))[cell_tab$cell]
  pres_pool <- which(sp_val$pa == 1)
  
  # fixed evaluation cells per species
  set.seed(seed_base + sp_num)
  eval_rows <- sample(nrow(cell_tab), n_eval)
  eval_env  <- cell_tab[eval_rows, vars]
  eval_pa   <- sp_val$pa[eval_rows]
  eval_suit <- sp_val$suitab[eval_rows]
  
  for (n_pres in n_pres_levels) {
    for (rep_id in seq_len(n_rep)) {
      
      out_file <- file.path(res_dir, sprintf("%s_p%03d_r%02d.RDS", sp_id, n_pres, rep_id))
      if (file.exists(out_file)) next
      message(format(Sys.time(), "%H:%M"), " | ", sp_id,
              " | n_pres = ", n_pres, " | replicate ", rep_id)
      
      # Presences: identical for all strategies × ratios (paired design)
      set.seed(seed_base + sp_num * 1e5 + n_pres * 100 + rep_id)
      pres_rows <- sample_rows(pres_pool, n_pres)
      pres_env  <- cell_tab[pres_rows, vars]
      
      prep <- prep_background(pres_rows, cell_tab, pc, pc_sf, tmpl, rng, vars,
                              buffer_width, use_thres, ecopa_bins)
      
      res <- purrr::pmap(runs, \(strategy, ratio) {
        n_bg <- n_pres * ratio_levels[[ratio]]
        t0   <- Sys.time()
        bg   <- sample_background(strategy, n_bg, prep, cell_tab, pc, vars,
                                  grid_res, ecopa_in)
        out  <- run_model(pres_env, bg, eval_env, eval_pa, eval_suit,
                          meta$true_pars, rng, env_means, q_trim)
        info <- tibble::tibble(
          species = sp_id, n_pres = n_pres, replicate = rep_id,
          strategy = strategy, ratio = ratio,
          n_bg_target = n_bg, n_bg_real = nrow(bg),
          time_s = round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
        )
        list(
          metrics = dplyr::bind_cols(info, out$metrics),
          curves  = if (!is.null(out$curves))
            dplyr::bind_cols(info[rep(1, nrow(out$curves)), ], out$curves)
        )
      })
      
      saveRDS(list(metrics = purrr::map_dfr(res, "metrics"),
                   curves  = purrr::map_dfr(res, "curves")),
              out_file)
    }
  }
}


# =============================== Test output ==============================

if (test_mode) {
  test <- purrr::map_dfr(list.files(res_dir, full.names = TRUE),
                         \(f) readRDS(f)$metrics)
  print(test[, c("n_pres", "strategy", "ratio", "n_bg_target", "n_bg_real",
                 "time_s", "auc", "tss", "boyce", "rc_rmse_mean", "error")],
        n = Inf)
}