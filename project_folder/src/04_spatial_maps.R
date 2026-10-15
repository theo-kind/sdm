# ================================= Notes =================================
#
# 04_spatial_maps.R
# Spatial comparison of modelled vs. true distributions.
# For selected species, the models from 02 are reproduced for the chosen
# ratios across all replicates (same seeds and drawing order), predicted
# across Romania and summarised per strategy x ratio:
#   - mean, SD and difference from truth of the continuous suitability
#     (one figure per ratio)
#   - replicate agreement: every replicate is converted to presence /
#     absence at its own maxTSS threshold; per cell the number of
#     replicates predicting presence is counted (0-n_rep). Per panel:
#     Fleiss' kappa (agreement among replicates) and mean sensitivity /
#     specificity against the true distribution (correctness).
#     Layout after Fig. 7 in the example paper: rows = strategies,
#     columns = ratios, true range as contour line.
# Output: GeoTIFF per species (mean, SD, count per strategy x ratio),
# per-replicate sensitivity/specificity, figures, summary tables.
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
map_n_pres  <- 250
map_ratios  <- c("2:1", "1:10", "1:100")   # columns of the agreement figure
plot_fact   <- 3        # aggregation for display only (30" -> 1.5')

# identical to 02_run_models.R (required for exact reproduction of the draws)
vars          <- c("bio_3", "bio_5", "bio_12", "bio_15")
ratio_levels  <- c("2:1" = 0.5, "1:1" = 1, "1:10" = 10, "1:40" = 40, "1:100" = 100)
strategies    <- c("ecopa", "use", "rand_geo", "rand_env",
                   "buf_in_10", "buf_out_10", "buf_in_50", "buf_out_50")
n_rep         <- 20
buffer_width  <- c(10000, 50000)
use_thres     <- 0.75
grid_res      <- 20
ecopa_bins    <- 30
ecopa_in      <- 0.5
seed_base     <- 2026

strategy_labels <- c(rand_geo = "Random (geo)",
                     buf_in_10 = "Buffer 10 km (in)", buf_out_10 = "Buffer 10 km (out)",
                     buf_in_50 = "Buffer 50 km (in)", buf_out_50 = "Buffer 50 km (out)",
                     rand_env = "Random (env)", use = "USE", ecopa = "EcoPA")
# display order only; `strategies` above keeps the draw order of 02
lab_map      <- c(true = "True", strategy_labels)
panel_levels <- unname(lab_map)

truth_line <- "#e34948"   # contour of the true range in the agreement maps

map_dir <- file.path(envrmt$path_data, "results", "maps")
dir.create(map_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(all(map_ratios %in% names(ratio_levels)))


# ============================ Prepare data ================================

env      <- terra::rast(file.path(envrmt$path_raster, "env_rou_30s.tif"))
tmpl     <- env[[1]]
cell_tab <- terra::as.data.frame(env[[vars]], cells = TRUE, na.rm = TRUE)
rng      <- sapply(cell_tab[vars], range)
sc       <- list(mu = colMeans(cell_tab[vars]),            # z-scaling for maxnet
                 sd = sapply(cell_tab[vars], sd))          # (as in 02_run_models.R)

pca   <- prcomp(cell_tab[vars], scale. = TRUE)
pc    <- as.data.frame(pca$x[, 1:2])
pc_sf <- sf::st_as_sf(pc, coords = c("PC1", "PC2"))

runs <- tidyr::expand_grid(strategy = strategies, ratio = names(ratio_levels))

# one layer key per strategy x mapped ratio, e.g. "rand_geo__1-100"
key  <- function(s, r) paste0(s, "__", gsub(":", "-", r))
keys <- tidyr::expand_grid(strategy = strategies, ratio = map_ratios) |>
  dplyr::mutate(key = key(strategy, ratio))


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

# Reproduce one presence set from 02 and predict the models of the mapped
# ratios across Romania. All background draws run in the same order as in
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
    if (!r %in% map_ratios) next
    mod <- tryCatch(fit_maxnet(pres_env, bg, sc), error = function(e) NULL)
    if (!is.null(mod)) preds[[key(s, r)]] <- predict_chunked(mod, cell_tab[vars])
  }
  preds
}

# Fleiss' kappa for binary ratings: cnt = number of "presence" votes per
# cell, m = number of raters (replicates)
fleiss_binary <- function(cnt, m) {
  P_i <- (cnt * (cnt - 1) + (m - cnt) * (m - cnt - 1)) / (m * (m - 1))
  p1  <- mean(cnt) / m
  Pe  <- p1^2 + (1 - p1)^2
  (mean(P_i) - Pe) / (1 - Pe)
}

# raster -> long data.frame for ggplot (aggregated for display)
to_df <- function(r, fun = "mean") {
  terra::aggregate(r, fact = plot_fact, fun = fun, na.rm = TRUE) |>
    terra::as.data.frame(xy = TRUE, na.rm = FALSE) |>
    tidyr::pivot_longer(-c(x, y), names_to = "layer", values_to = "value") |>
    dplyr::filter(!is.na(value)) |>
    dplyr::mutate(panel = factor(lab_map[layer], levels = panel_levels))
}

# one map row (9 columns: true + 8 strategies, empty panels stay blank)
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
agree_list   <- list()

for (sp_id in map_species) {

  sp_num    <- as.integer(sub("sp", "", sp_id))
  sp_r      <- terra::rast(file.path(envrmt$path_virtual_species, paste0(sp_id, ".tif")))
  sp_val    <- sp_r[cell_tab$cell]
  pres_pool <- which(sp_val$pa == 1)
  tag       <- sprintf("%s_p%03d", sp_id, map_n_pres)
  tif_file  <- file.path(map_dir, paste0(tag, "_maps.tif"))
  rep_file  <- file.path(map_dir, paste0(tag, "_replicate_sens_spec.csv"))

  # ---- mean, SD and presence count (computed once, then read from file) ----
  cache_ok <- file.exists(tif_file) && file.exists(rep_file) &&
    all(paste0(keys$key, "_mean") %in% names(terra::rast(tif_file)))
  if (cache_ok) {
    maps      <- terra::rast(tif_file)
    rep_stats <- utils::read.csv(rep_file)
  } else {
    n   <- setNames(rep(0, nrow(keys)), keys$key)
    s1  <- s2 <- setNames(lapply(keys$key, \(.) numeric(nrow(cell_tab))), keys$key)
    cnt <- setNames(lapply(keys$key, \(.) integer(nrow(cell_tab))), keys$key)
    rep_stats <- list()

    for (rep_id in seq_len(n_rep)) {
      message(format(Sys.time(), "%H:%M"), " | ", sp_id, " | replicate ", rep_id)
      preds <- replay_set(sp_num, rep_id, pres_pool)
      for (k in names(preds)) {
        p <- preds[[k]]
        s1[[k]] <- s1[[k]] + p
        s2[[k]] <- s2[[k]] + p^2
        n[k]    <- n[k] + 1
        # binary prediction of this replicate at its maxTSS threshold
        bin      <- p >= eval_auc_tss(p, sp_val$pa)$thr
        cnt[[k]] <- cnt[[k]] + bin
        rep_stats[[length(rep_stats) + 1]] <- data.frame(
          key = k, replicate = rep_id,
          sens = mean(bin[sp_val$pa == 1]),
          spec = mean(!bin[sp_val$pa == 0]))
      }
    }
    rep_stats <- dplyr::bind_rows(rep_stats)

    maps <- terra::rast(lapply(keys$key, \(k) {
      m  <- s1[[k]] / n[k]
      sd <- sqrt(pmax(s2[[k]] / n[k] - m^2, 0) * n[k] / (n[k] - 1))
      c(to_rast(m, paste0(k, "_mean")), to_rast(sd, paste0(k, "_sd")),
        to_rast(cnt[[k]], paste0(k, "_count")))
    }))
    terra::writeRaster(maps, tif_file, overwrite = TRUE)
    utils::write.csv(rep_stats, rep_file, row.names = FALSE)
  }

  true_r <- sp_r[["suitab"]]; names(true_r) <- "true"

  # common colour range of the SD maps across all mapped ratios (so that the
  # figures of different ratios are comparable)
  sd_max <- max(terra::global(maps[[paste0(keys$key, "_sd")]], "max",
                              na.rm = TRUE)$max, na.rm = TRUE)
  sd_scale   <- scale_fill_viridis_c(option = "magma", limits = c(0, sd_max))
  mean_scale <- scale_fill_viridis_c(limits = c(0, 1))
  diff_scale <- scale_fill_gradient2(low = "#2166AC", mid = "grey95",
                                     high = "#B2182B", limits = c(-1, 1))
  all_mean <- all_sd <- all_diff <- list()   # for the combined figure


  # ---- per ratio: mean / SD / difference maps + summary statistics ----
  for (r in map_ratios) {
    lyr    <- key(strategies, r)
    mean_r <- maps[[paste0(lyr, "_mean")]]; names(mean_r) <- strategies
    sd_r   <- maps[[paste0(lyr, "_sd")]];   names(sd_r)   <- strategies
    diff_r <- mean_r - true_r;              names(diff_r) <- strategies

    mv <- mean_r[cell_tab$cell]
    for (s in strategies) {
      m <- mv[[s]]
      if (all(is.na(m))) next
      at  <- eval_auc_tss(m, sp_val$pa)
      bin <- m >= at$thr
      summary_list[[paste(sp_id, s, r)]] <- tibble::tibble(
        species = sp_id, strategy = strategy_labels[[s]],
        n_pres = map_n_pres, ratio = r,
        auc = at$auc, tss = at$tss,
        schoener_d = schoener_d(m, sp_val$suitab),
        overpred   = mean(bin & sp_val$pa == 0),   # share of the country area
        underpred  = mean(!bin & sp_val$pa == 1)
      )
    }

    d_mean <- dplyr::bind_rows(to_df(true_r), to_df(mean_r))
    d_sd   <- to_df(sd_r)
    d_diff <- to_df(diff_r)
    rlab   <- paste("P:PA", r)
    all_mean[[r]] <- dplyr::mutate(
      if (r == map_ratios[1]) d_mean else dplyr::filter(d_mean, layer != "true"),
      ratio = rlab)
    all_sd[[r]]   <- dplyr::mutate(d_sd, ratio = rlab)
    all_diff[[r]] <- dplyr::mutate(d_diff, ratio = rlab)

    p <- patchwork::wrap_plots(
      plot_row(d_mean, mean_scale,
               "Suitability (mean across replicates)", strips = TRUE),
      plot_row(d_sd, sd_scale, "Uncertainty (SD across replicates)"),
      plot_row(d_diff, diff_scale, "Difference from truth (model \u2212 true)"),
      ncol = 1
    ) +
      patchwork::plot_annotation(
        title = paste0(sp_id, " | n = ", map_n_pres, " | Ratio ", r))

    ggsave(file.path(envrmt$path_figures,
                     sprintf("04_maps_%s_%s.png", tag, gsub(":", "-", r))),
           p, width = 20, height = 7.2, dpi = 300, bg = "white")
  }


  # ---- combined figure: all mapped ratios, common colour ranges ----
  # three blocks (mean, SD, difference); rows = ratios, columns = true +
  # strategies (true suitability only in the first row)
  ratio_fac <- function(d) dplyr::mutate(
    d, ratio = factor(ratio, levels = paste("P:PA", map_ratios)))
  plot_block <- function(d, fill_scale, title, strips = FALSE) {
    ggplot(ratio_fac(d), aes(x, y, fill = value)) +
      geom_raster() +
      coord_quickmap(expand = FALSE) +
      facet_grid(ratio ~ panel, drop = FALSE, switch = "y") +
      fill_scale +
      labs(title = title, fill = NULL) +
      theme_void(base_size = 9) +
      theme(strip.text.x = if (strips) element_text(face = "bold") else element_blank(),
            strip.text.y.left = element_text(angle = 0, hjust = 1),
            panel.spacing.y = unit(4, "pt"),
            plot.title = element_text(size = 9, face = "bold"))
  }
  p_all <- patchwork::wrap_plots(
    plot_block(dplyr::bind_rows(all_mean), mean_scale,
               "Suitability (mean across replicates)", strips = TRUE),
    plot_block(dplyr::bind_rows(all_sd), sd_scale,
               "Uncertainty (SD across replicates)"),
    plot_block(dplyr::bind_rows(all_diff), diff_scale,
               "Difference from truth (model \u2212 true)"),
    ncol = 1
  ) +
    patchwork::plot_annotation(title = paste0(sp_id, " | n = ", map_n_pres))

  ggsave(file.path(envrmt$path_figures, sprintf("04_maps_%s_all_ratios.png", tag)),
         p_all, width = 20, height = 1.5 + 4.6 * length(map_ratios), dpi = 300,
         bg = "white", limitsize = FALSE)


  # ---- replicate agreement (strategies x ratios) ----
  agree <- rep_stats |>
    dplyr::group_by(key) |>
    dplyr::summarise(n_ok = dplyr::n(), sens = mean(sens), spec = mean(spec),
                     .groups = "drop") |>
    dplyr::inner_join(keys, by = "key") |>
    dplyr::mutate(
      fleiss = purrr::map2_dbl(key, n_ok, \(k, m)
        fleiss_binary(maps[[paste0(k, "_count")]][cell_tab$cell][[1]], m)),
      species = sp_id, .before = 1)
  agree_list[[sp_id]] <- agree

  cnt_r <- maps[[paste0(agree$key, "_count")]]
  names(cnt_r) <- agree$key
  df_cnt <- terra::aggregate(cnt_r, fact = plot_fact, fun = "mean", na.rm = TRUE) |>
    terra::as.data.frame(xy = TRUE, na.rm = FALSE) |>
    tidyr::pivot_longer(-c(x, y), names_to = "key", values_to = "value") |>
    dplyr::filter(!is.na(value)) |>
    dplyr::left_join(keys, by = "key") |>
    dplyr::mutate(value = ifelse(value < 0.5, NA, value))   # 0 votes -> grey

  # true range: contour of the true suitability at 0.5 (= probability of
  # presence 0.5 with the linear PA conversion in 01)
  df_true <- terra::aggregate(true_r, fact = plot_fact, fun = "mean",
                              na.rm = TRUE) |>
    terra::as.data.frame(xy = TRUE, na.rm = FALSE) |>
    dplyr::mutate(true = dplyr::coalesce(true, 0))

  lab <- agree |>
    dplyr::mutate(label = sprintf("Fleiss: %.2f\nSens: %.2f\nSpec: %.2f",
                                  fleiss, sens, spec))

  fac <- function(d) dplyr::mutate(d,
    strategy = factor(strategy_labels[strategy], levels = strategy_labels),
    ratio    = factor(paste("P:PA", ratio), levels = paste("P:PA", map_ratios)))

  p_agree <- ggplot(fac(df_cnt), aes(x, y)) +
    geom_raster(aes(fill = value)) +
    geom_contour(data = df_true, aes(x, y, z = true), breaks = 0.5,
                 colour = truth_line, linewidth = 0.25, inherit.aes = FALSE) +
    geom_text(data = fac(lab), aes(x = -Inf, y = -Inf, label = label),
              hjust = -0.05, vjust = -0.15, size = 2, lineheight = 0.9) +
    coord_quickmap(expand = FALSE) +
    facet_grid(strategy ~ ratio, switch = "y") +
    scale_fill_viridis_c(limits = c(0.5, n_rep), breaks = c(1, 5, 10, 15, 20),
                         na.value = "grey80",
                         name = "Replicates predicting presence") +
    guides(fill = guide_colourbar(title.position = "top")) +
    labs(title = paste0(sp_id, " | n = ", map_n_pres),
         caption = paste0("Each replicate thresholded at its maxTSS. Grey: no ",
                          "replicate predicts presence; red line: true suitability = 0.5.\n",
                          "Fleiss' kappa: agreement among replicates; Sens / Spec: ",
                          "mean sensitivity / specificity against the true range.")) +
    theme_void(base_size = 9) +
    theme(strip.text.x = element_text(face = "bold", margin = margin(b = 3)),
          strip.text.y.left = element_text(angle = 0, hjust = 1,
                                           margin = margin(r = 4)),
          legend.position = "bottom",
          legend.key.width = unit(1.2, "cm"), legend.key.height = unit(0.3, "cm"),
          plot.caption = element_text(hjust = 0, size = 7),
          plot.title = element_text(face = "bold"),
          plot.background = element_rect(fill = "white", colour = NA))

  ggsave(file.path(envrmt$path_figures, sprintf("04_replicate_agreement_%s.png", tag)),
         p_agree, width = 8.5, height = 13, dpi = 300, bg = "white")
}


# ================================ Summary statistics ======================

map_summary <- dplyr::bind_rows(summary_list)
readr::write_csv(map_summary,
                 file.path(map_dir, sprintf("map_summary_p%03d.csv", map_n_pres)))
print(map_summary, n = Inf)

agree_summary <- dplyr::bind_rows(agree_list) |>
  dplyr::mutate(strategy = factor(strategy_labels[strategy], levels = strategy_labels),
                ratio    = factor(ratio, levels = map_ratios)) |>
  dplyr::arrange(species, strategy, ratio) |>
  dplyr::select(species, strategy, ratio, n_ok, fleiss, sens, spec)
readr::write_csv(agree_summary,
                 file.path(map_dir, sprintf("replicate_agreement_p%03d.csv", map_n_pres)))
print(agree_summary, n = Inf)
