# ================================= Notes =================================
#
# 03c_rmse_replicates.R
# Accuracy and spread of the response curves, figure after Bazzichetto
# et al. (2023, Fig. with RMSE vs. sampling effort and replicate curves):
#   (a) RMSE between modelled and true response curve vs. ratio, one line
#       per strategy, one panel per variable.
#       RMSE = sqrt(mean over replicates and x of (pred - true)^2) per
#       species (bias and spread across replicates combined), then
#       averaged over all species.
#   (b) Example species: modelled curves of all replicates (grey) and true
#       curve (orange) for the best and the worst strategy (lowest /
#       highest mean RMSE in (a) at the highest ratio), at the lowest and
#       the highest ratio.
# Reads the run files of the current design (data/results/runs).
#
# ==========================================================================

## ================================= Load set up =================================
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

n_pres_sel   <- 250                                        # as in 02 / 03
ratio_order  <- c("2:1", "1:1", "1:10", "1:40", "1:100")  # as in 02 / 03
spag_ratios  <- ratio_order[c(1, length(ratio_order))]     # rows in (b)
spag_species <- "sp05"   # NULL = automatic (species with the median mean RMSE),
                       # or e.g. "sp05"
spag_strats <- c("rand_geo", "buf_out_50")   # NULL = automatic (best / worst strategy in (a)),
                       # or e.g. c("rand_geo", "ecopa")

vars       <- c("bio_3", "bio_5", "bio_12", "bio_15")
var_labels <- c(bio_3 = "BIO3", bio_5 = "BIO5", bio_12 = "BIO12", bio_15 = "BIO15")
strategy_labels <- c(rand_geo = "Random (geo)",
                     buf_in_10 = "Buffer 10 km (in)", buf_out_10 = "Buffer 10 km (out)",
                     buf_in_50 = "Buffer 50 km (in)", buf_out_50 = "Buffer 50 km (out)",
                     rand_env = "Random (env)", use = "USE", ecopa = "EcoPA")
truth_col <- "#D55E00"   # true response (Okabe-Ito vermillion, contrasts with grey)

res_dir <- file.path(envrmt$path_data, "results", "runs")
out_dir <- file.path(envrmt$path_data, "results", "summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ================================ Load data ===============================

files <- list.files(res_dir, pattern = sprintf("_p%03d_r\\d+\\.RDS$", n_pres_sel),
                    full.names = TRUE)
message(length(files), " result files found")

curves <- purrr::map_dfr(files, \(f) readRDS(f)$curves) |>
  dplyr::filter(ratio %in% ratio_order) |>
  dplyr::mutate(strategy = factor(strategy, levels = names(strategy_labels)),
                ratio    = factor(ratio, levels = ratio_order),
                variable = factor(variable, levels = vars))


# ===================== (a) RMSE across replicates =========================

# per species: RMSE over replicates and x (bias + spread)
rmse_sp <- curves |>
  dplyr::group_by(species, strategy, ratio, variable) |>
  dplyr::summarise(rmse = sqrt(mean((pred - true)^2)),
                   n_rep = dplyr::n_distinct(replicate), .groups = "drop")
readr::write_csv(rmse_sp, file.path(out_dir, "rc_rmse_replicates.csv"))

# ---- bias and spread of the whole response curve ----
# Per species x strategy x ratio x variable, across the replicates:
#   bias = distance between the mean modelled curve and the true curve
#          (root mean square over x)
#   SD   = spread of the replicate curves around their mean curve
#          (root mean square over x)
# so that rmse^2 = bias^2 + SD^2 (rmse as in rmse_sp above).
curve_bv <- curves |>
  dplyr::group_by(species, strategy, ratio, variable, x) |>
  dplyr::summarise(dev = mean(pred) - true[1],
                   var = mean((pred - mean(pred))^2), .groups = "drop") |>
  dplyr::group_by(species, strategy, ratio, variable) |>
  dplyr::summarise(bias = sqrt(mean(dev^2)), sd = sqrt(mean(var)),
                   .groups = "drop")

# table: mean over species and variables
tab_curve_bv <- curve_bv |>
  dplyr::group_by(strategy, ratio) |>
  dplyr::summarise(bias = mean(bias), sd = mean(sd), .groups = "drop")
readr::write_csv(tab_curve_bv, file.path(out_dir, "table_curve_bias_sd.csv"))
tab_curve_bv |>
  dplyr::mutate(strategy = strategy_labels[as.character(strategy)],
                cell = sprintf("%.3f / %.3f", bias, sd)) |>
  dplyr::select(strategy, ratio, cell) |>
  tidyr::pivot_wider(names_from = ratio, values_from = cell) |>
  readr::write_csv(file.path(out_dir, "table_curve_bias_sd_wide.csv"))

# mean over species
rmse_all <- rmse_sp |>
  dplyr::group_by(strategy, ratio, variable) |>
  dplyr::summarise(rmse = mean(rmse), .groups = "drop") |>
  dplyr::mutate(family = factor(strategy_family[as.character(strategy)],
                                levels = names(family_lty)),
                strategy_lab = factor(strategy_labels[as.character(strategy)],
                                      levels = strategy_labels))

p_rmse <- ggplot(rmse_all, aes(ratio, rmse, colour = strategy_lab,
                               fill = strategy_lab, shape = strategy_lab,
                               group = strategy_lab)) +
  geom_line(aes(linetype = family), linewidth = 0.45, alpha = 0.8) +
  geom_point(size = 1.8, stroke = 0.6) +
  facet_wrap(~ variable, nrow = 1,
             labeller = labeller(variable = var_labels)) +
  scale_colour_manual(values = by_label(strategy_cols, strategy_labels)) +
  scale_fill_manual(values = by_label(strategy_fills, strategy_labels)) +
  scale_shape_manual(values = by_label(strategy_shapes, strategy_labels)) +
  scale_linetype_manual(values = family_lty) +
  # common y axis for all variables: 0-0.5 (extended if a value is larger)
  scale_y_continuous(limits = c(0, max(0.5, rmse_all$rmse)),
                     breaks = seq(0, 1, 0.1), expand = c(0, 0)) +
  labs(x = "Ratio (P:PA)", y = "RMSE (response curve)", title = "(a)",
       colour = "Strategy", fill = "Strategy", shape = "Strategy",
       linetype = "Space") +
  theme_bw(base_size = 9) +
  theme(strip.background = element_rect(colour = "grey40"),
        strip.text = element_text(size = 9),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(face = "bold"),
        legend.position = "bottom", legend.box = "vertical")


# ======================= (b) Replicate curves =============================

# best and worst strategy: mean RMSE (over species and variables) at the
# highest ratio
if (is.null(spag_strats)) {
  r <- rmse_all |>
    dplyr::filter(ratio == ratio_order[length(ratio_order)]) |>
    dplyr::group_by(strategy) |>
    dplyr::summarise(rmse = mean(rmse)) |>
    dplyr::arrange(rmse)
  spag_strats <- as.character(c(r$strategy[1], r$strategy[nrow(r)]))
}

# example species: mean RMSE closest to the median over species
if (is.null(spag_species)) {
  s <- rmse_sp |>
    dplyr::group_by(species) |>
    dplyr::summarise(rmse = mean(rmse))
  spag_species <- s$species[which.min(abs(s$rmse - median(s$rmse)))]
}

spag_block <- function(sp, strat, show_y = TRUE) {
  d <- dplyr::filter(curves, species == sp, strategy == strat,
                     ratio %in% spag_ratios) |>
    dplyr::mutate(ratio = droplevels(ratio))
  truth <- dplyr::distinct(d, ratio, variable, x, true)
  ggplot(d, aes(x, pred, group = replicate)) +
    geom_line(colour = "grey20", alpha = 0.3, linewidth = 0.3) +
    geom_line(data = truth, aes(x, true), inherit.aes = FALSE,
              colour = truth_col, linewidth = 0.7) +
    facet_grid(ratio ~ variable, scales = "free_x",
               labeller = labeller(variable = var_labels,
                                   ratio = \(r) paste("P:PA", r))) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(x = strategy_labels[[strat]],
         y = if (show_y) "Response (scaled)" else NULL) +
    theme_bw(base_size = 8) +
    theme(strip.background = element_rect(colour = "grey40"),
          axis.title.x = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1, size = 6),
          # ratio strips only on the right-hand block
          strip.text.y = if (show_y) element_blank() else element_text(angle = 0))
}



# ================================ Figure ==================================

title_b <- ggplot() + theme_void() +
  labs(title = sprintf("(b) Example species %s: best vs. worst strategy at %s",
                       spag_species, ratio_order[length(ratio_order)])) +
  theme(plot.title = element_text(face = "bold", size = 9))

fig <- patchwork::wrap_plots(
  p_rmse, title_b,
  spag_block(spag_species, spag_strats[1], TRUE) |
    spag_block(spag_species, spag_strats[2], FALSE),
  ncol = 1, heights = c(3.2, 0.2, 2.6)) +
  patchwork::plot_annotation(
    caption = paste0("(a) RMSE between modelled and true response curve across ",
                     "replicates (bias and spread), mean over all species.\n",
                     "(b) Grey: modelled curves of all replicates; orange: true ",
                     "response. Rows: ratio (P:PA)."),
    theme = theme(plot.caption = element_text(hjust = 0, size = 8)))

ggsave(file.path(envrmt$path_figures, "03c_rmse_replicates.png"),
       fig, width = 12, height = 9, dpi = 300)
