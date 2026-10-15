# ================================= Notes =================================
#
# 06_bias_comparison.R
# Comparison of the two presence sampling scenarios:
#   - unbiased: presences from all of Romania   (src/,        summary/)
#   - biased:   presences from 17 GADM counties (src_blocks/, summary_blocks/)
# Rows: all evaluation metrics (AUC, TSS, Kappa, Boyce, Spearman, RMSE of
# the response curves); columns: scenario. Same y axis within a row, so that the change
# between the scenarios can be read directly.
# Reads metrics_all.RDS written by 03_analysis.R of both pipelines.
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

scenario_files <- c(
  "Unbiased presences" = file.path(envrmt$path_data, "results", "summary",
                                   "metrics_all.RDS"),
  "Biased presences (17 counties)" = file.path(envrmt$path_data, "results",
                                               "summary_blocks", "metrics_all.RDS"))
# all metrics of 03 (performance metrics first, then agreement with the truth)
metric_set <- c(auc          = "AUC",
                tss          = "TSS",
                kappa        = "Kappa",
                boyce        = "Boyce",
                spearman     = "Spearman (suitability)",
                rc_rmse_mean = "RMSE (response curves)")
strategy_labels <- c(rand_geo = "Random (geo)",
                     buf_in_10 = "Buffer 10 km (in)", buf_out_10 = "Buffer 10 km (out)",
                     buf_in_50 = "Buffer 50 km (in)", buf_out_50 = "Buffer 50 km (out)",
                     rand_env = "Random (env)", use = "USE", ecopa = "EcoPA")


# ================================ Load data ===============================

d <- purrr::imap_dfr(scenario_files, \(f, sc) readRDS(f) |>
                       dplyr::mutate(scenario = sc)) |>
  dplyr::filter(is.na(error)) |>
  dplyr::mutate(scenario = factor(scenario, levels = names(scenario_files))) |>
  tidyr::pivot_longer(dplyr::all_of(names(metric_set)),
                      names_to = "metric", values_to = "value") |>
  dplyr::mutate(metric = factor(metric_set[metric], levels = metric_set),
                family = factor(strategy_family[names(strategy_labels)[strategy]],
                                levels = names(family_lty)))

# medians per scenario x strategy x ratio (also as table)
med <- d |>
  dplyr::group_by(scenario, metric, strategy, ratio) |>
  dplyr::summarise(median = median(value), .groups = "drop")
readr::write_csv(med, file.path(envrmt$path_data, "results",
                                "bias_comparison_medians.csv"))


# ================================= Figure =================================

pd <- position_dodge(width = 0.6)
p <- ggplot(d, aes(ratio, value, colour = strategy, fill = strategy,
                   shape = strategy, group = strategy)) +
  stat_summary(aes(linetype = family), fun = median, geom = "line",
               position = pd, linewidth = 0.4) +
  stat_summary(fun = median,
               fun.min = \(x) quantile(x, 0.25),
               fun.max = \(x) quantile(x, 0.75),
               position = pd, size = 0.45, stroke = 0.6) +
  facet_grid(metric ~ scenario, scales = "free_y", switch = "y",
             labeller = labeller(metric = label_wrap_gen(14))) +
  scale_colour_manual(values = by_label(strategy_cols, strategy_labels)) +
  scale_fill_manual(values = by_label(strategy_fills, strategy_labels)) +
  scale_shape_manual(values = by_label(strategy_shapes, strategy_labels)) +
  scale_linetype_manual(values = family_lty) +
  guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
  labs(x = "Ratio (P:PA)", y = NULL, colour = "Strategy", fill = "Strategy",
       shape = "Strategy", linetype = "Space",
       caption = "Median and interquartile range across species \u00d7 replicates") +
  theme_bw() +
  theme(legend.position = "bottom", strip.placement = "outside",
        strip.background.y = element_blank(),
        strip.text.y.left = element_text(angle = 90, size = 10),
        strip.text.x = element_text(face = "bold"))

ggsave(file.path(envrmt$path_figures, "06_bias_comparison.png"),
       p, width = 10, height = 13, dpi = 300)
