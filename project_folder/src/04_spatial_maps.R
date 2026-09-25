# ================================= Notes =================================
#
# 04_spatial_maps.R
# Spatial comparison of modelled vs. true distributions.
# For selected species, the models from 02 are reproduced for one
# n_pres × ratio combination across all replicates (same seeds and
# drawing order), predicted across Romania and summarised per strategy
# as mean and SD.
# Output: GeoTIFF per species (mean + SD per strategy), figure per
# species, table with summary statistics of the averaged maps.
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

map_species <- c("sp01")
map_n_pres  <- 100
map_ratio   <- "1:10"
plot_fact   <- 3        # aggregation for display only (30" -> 1.5')

# identical to 02_run_models.R (required for exact reproduction of the draws)
vars          <- c("bio_3", "bio_5", "bio_12", "bio_15")
ratio_levels  <- c("2:1" = 0.5, "1:1" = 1, "1:10" = 10, "1:100" = 100)
strategies    <- c("ecopa", "use", "rand_geo", "rand_env", "buf_in", "buf_out")
n_rep         <- 10
buffer_width  <- 50000
use_thres     <- 0.75
grid_res      <- 20
ecopa_bins    <- 30
ecopa_in      <- 0.5
seed_base     <- 2026

strategy_labels <- c(ecopa = "EcoPA", use = "USE", rand_geo = "Random (geo)",
                     rand_env = "Random (env)", buf_in = "Buffer (in)",
                     buf_out = "Buffer (out)")
lab_map      <- c(true = "True", strategy_labels)
panel_levels <- unname(lab_map)

agree_levels <- c("correct absence", "correct presence",
                  "overprediction", "underprediction")
agree_cols   <- c("correct absence"  = "grey90",
                  "correct presence" = "#1B7837",
                  "overprediction"   = "#E08214",
                  "underprediction"  = "#762A83")

map_dir <- file.path(envrmt$path_data, "results", "maps")
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)


# ============================ Prepare data ================================

env      <- terra::rast(file.path(envrmt$path_raster, "env_rou_30s.tif"))
tmpl     <- env[[1]]
cell_tab <- terra::as.data.frame(env[[vars]], cells = TRUE, na.rm = TRUE)
rng      <- sapply(cell_tab[vars], range)

pca   <- prcomp(cell_tab[vars], scale. = TRUE)
pc    <- as.data.frame(pca$x[, 1:2])
pc_sf <- sf::st_as_sf(pc, coords = c("PC1", "PC2"))

runs <- tidyr::expand_grid(strategy = strategies, ratio = names(ratio_levels))


# ================================ Functions ===============================

# Predict in chunks (saves memory with ~400,000 cells)
predict_chunked <- function(mod, newdata, chunk = 1e5) {
  idx <- split(seq_len(nrow(newdata)), ceiling(seq_len(nrow(newdata)) / chunk))
  unlist(lapply(idx, \(i) as.vector(
    predict(mod, newdata[i, , drop = FALSE], type = "cloglog"))), use.names = FALSE)
}

# value vector (rows of cell_tab) -> raster
to_rast <- function(v, name) {
  r <- tmpl * NA
  r[cell_tab$cell] <- v
  names(r) <- name
  r
}

# Reproduce one presence set from 02 and predict the models of the chosen
# ratio across Romania. All background draws run in the same order as in
# 02 so that the random number stream stays identical.
replay_set <- function(sp_num, rep_id, pres_pool) {
  set.seed(seed_base + sp_num * 1e5 + map_n_pres * 100 + rep_id)
  pres_rows <- sample_rows(pres_pool, map_n_pres)
  pres_env  <- cell_tab[pres_rows, vars]
  prep <- prep_background(pres_rows, cell_tab, pc, pc_sf, tmpl, rng, vars,
                          buffer_width, use_thres, ecopa_bins)
  preds <- list()
  for (k in seq_len(nrow(runs))) {
    s  <- runs$strategy[k]
    r  <- runs$ratio[k]
    bg <- sample_background(s, map_n_pres * ratio_levels[[r]], prep, cell_tab,
                            pc, vars, grid_res, ecopa_in)
    if (r != map_ratio) next
    mod <- tryCatch(fit_maxnet(pres_env, bg), error = function(e) NULL)
    if (!is.null(mod)) preds[[s]] <- predict_chunked(mod, cell_tab[vars])
  }
  preds
}

# raster -> long data.frame for ggplot (aggregated for display)
to_df <- function(r, fun = "mean") {
  terra::aggregate(r, fact = plot_fact, fun = fun, na.rm = TRUE) |>
    terra::as.data.frame(xy = TRUE, na.rm = FALSE) |>
    tidyr::pivot_longer(-c(x, y), names_to = "layer", values_to = "value") |>
    dplyr::filter(!is.na(value)) |>
    dplyr::mutate(panel = factor(lab_map[layer], levels = panel_levels))
}

# one map row (7 columns: true + 6 strategies, empty panels stay blank)
plot_row <- function(d, fill_scale, title, strips = FALSE) {
  ggplot(d, aes(x, y, fill = value)) +
    geom_raster() +
    coord_quickmap(expand = FALSE) +
    facet_wrap(~ panel, nrow = 1, drop = FALSE) +
    fill_scale +
    labs(title = title, fill = NULL) +
    theme_void(base_size = 9) +
    theme(strip.text = if (strips) element_text(face = "bold") else element_blank(),
          plot.title = element_text(size = 9, face = "bold"))
}


# ================================ Main loop ===============================

summary_list <- list()

for (sp_id in map_species) {
  
  sp_num    <- as.integer(sub("sp", "", sp_id))
  sp_r      <- terra::rast(file.path(envrmt$path_virtual_species, paste0(sp_id, ".tif")))
  sp_val    <- sp_r[cell_tab$cell]
  pres_pool <- which(sp_val$pa == 1)
  tag       <- sprintf("%s_p%03d_%s", sp_id, map_n_pres, gsub(":", "-", map_ratio))
  tif_file  <- file.path(map_dir, paste0(tag, ".tif"))
  
  # ---- mean and SD maps (computed once, then read from file) ----
  if (file.exists(tif_file)) {
    maps <- terra::rast(tif_file)
  } else {
    n  <- setNames(rep(0, length(strategies)), strategies)
    s1 <- s2 <- setNames(lapply(strategies, \(.) numeric(nrow(cell_tab))), strategies)
    
    for (rep_id in seq_len(n_rep)) {
      message(format(Sys.time(), "%H:%M"), " | ", sp_id, " | replicate ", rep_id)
      preds <- replay_set(sp_num, rep_id, pres_pool)
      for (s in names(preds)) {
        s1[[s]] <- s1[[s]] + preds[[s]]
        s2[[s]] <- s2[[s]] + preds[[s]]^2
        n[s]    <- n[s] + 1
      }
    }
    
    maps <- terra::rast(lapply(strategies, \(s) {
      m  <- s1[[s]] / n[s]
      sd <- sqrt(pmax(s2[[s]] / n[s] - m^2, 0) * n[s] / (n[s] - 1))
      c(to_rast(m, paste0(s, "_mean")), to_rast(sd, paste0(s, "_sd")))
    }))
    terra::writeRaster(maps, tif_file, overwrite = TRUE)
  }
  
  mean_r <- maps[[paste0(strategies, "_mean")]]; names(mean_r) <- strategies
  sd_r   <- maps[[paste0(strategies, "_sd")]];   names(sd_r)   <- strategies
  diff_r <- mean_r - sp_r[["suitab"]];           names(diff_r) <- strategies
  true_r <- sp_r[["suitab"]];                    names(true_r) <- "true"
  
  # ---- binary agreement + summary statistics ----
  mv <- mean_r[cell_tab$cell]
  agree_list <- list()
  for (s in strategies) {
    m <- mv[[s]]
    if (all(is.na(m))) next
    at  <- eval_auc_tss(m, sp_val$pa)
    bin <- m >= at$thr
    cat <- ifelse(sp_val$pa == 1 &  bin, 1L,
                  ifelse(sp_val$pa == 0 &  bin, 2L,
                         ifelse(sp_val$pa == 1 & !bin, 3L, 0L)))
    agree_list[[s]] <- to_rast(cat, s)
    
    summary_list[[paste(sp_id, s)]] <- tibble::tibble(
      species = sp_id, strategy = strategy_labels[[s]],
      n_pres = map_n_pres, ratio = map_ratio,
      auc = at$auc, tss = at$tss,
      schoener_d = schoener_d(m, sp_val$suitab),
      overpred   = mean(cat == 2),   # share of the country area
      underpred  = mean(cat == 3)
    )
  }
  
  # ---- figure ----
  df_mean  <- dplyr::bind_rows(to_df(true_r), to_df(mean_r))
  df_sd    <- to_df(sd_r)
  df_diff  <- to_df(diff_r)
  df_agree <- to_df(terra::rast(agree_list), fun = "modal") |>
    dplyr::mutate(value = factor(agree_levels[value + 1], levels = agree_levels))
  
  p <- patchwork::wrap_plots(
    plot_row(df_mean, scale_fill_viridis_c(limits = c(0, 1)),
             "Suitability (mean across replicates)", strips = TRUE),
    plot_row(df_sd, scale_fill_viridis_c(option = "magma", limits = c(0, NA)),
             "Uncertainty (SD across replicates)"),
    plot_row(df_diff, scale_fill_gradient2(low = "#2166AC", mid = "grey95",
                                           high = "#B2182B", limits = c(-1, 1)),
             "Difference from truth (model − true)"),
    plot_row(df_agree, scale_fill_manual(values = agree_cols, drop = FALSE),
             "Agreement (threshold: maxTSS)"),
    ncol = 1
  ) +
    patchwork::plot_annotation(
      title = paste0(sp_id, " | n = ", map_n_pres, " | Ratio ", map_ratio))
  
  ggsave(file.path(envrmt$path_figures, paste0("04_maps_", tag, ".png")),
         p, width = 16, height = 9.5, dpi = 300)
}


# ================================ Summary statistics ======================

map_summary <- dplyr::bind_rows(summary_list)
readr::write_csv(map_summary,
                 file.path(map_dir, sprintf("map_summary_p%03d_%s.csv",
                                            map_n_pres, gsub(":", "-", map_ratio))))
print(map_summary, n = Inf)