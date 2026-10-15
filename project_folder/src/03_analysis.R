# ================================= Notes =================================
#
# 03_analysis.R
# Combines the results from 02_run_models.R, computes failure and
# shortfall rates, creates figures (metrics, response curves) and
# fits linear mixed models (species + presence set as random effects).
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

res_dir <- file.path(envrmt$path_data, "results", "runs")
out_dir <- file.path(envrmt$path_data, "results", "summary")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

strategy_labels <- c(rand_geo = "Random (geo)",
                     buf_in_10 = "Buffer 10 km (in)", buf_out_10 = "Buffer 10 km (out)",
                     buf_in_50 = "Buffer 50 km (in)", buf_out_50 = "Buffer 50 km (out)",
                     rand_env = "Random (env)", use = "USE", ecopa = "EcoPA")
# order: geographic, then environmental strategies (style: functions/plot_style.R)
ratio_order  <- c("2:1", "1:1", "1:10", "1:40", "1:100")   # must match ratio_levels in 02
metrics_perf <- c(auc = "AUC", tss = "TSS", kappa = "Kappa", boyce = "Boyce")
metrics_true <- c(spearman = "Spearman (suitability)",
                  rc_rmse_mean = "RMSE (response curves)")
n_pres_sel   <- 250   # current design (must match n_pres_levels in 02)
rc_n_pres    <- n_pres_sel   # n_pres for the response curve figures

emmeans::emm_options(lmer.df = "asymptotic")  # fast for large numbers of observations


# ============================ Load results ================================

# only result files of the current design (older runs with other n_pres are ignored)
files <- list.files(res_dir, pattern = sprintf("_p%03d_r\\d+\\.RDS$", n_pres_sel),
                    full.names = TRUE)
message(length(files), " result files found")
res <- purrr::map(files, readRDS)

prep_factors <- function(d) {
  d |> dplyr::mutate(
    strategy = factor(strategy_labels[strategy], levels = strategy_labels),
    ratio    = factor(ratio, levels = ratio_order),
    n_pres   = factor(n_pres),
    species  = factor(species),
    set      = interaction(species, n_pres, replicate, drop = TRUE)  # presence set
  )
}

metrics <- purrr::map_dfr(res, "metrics") |> prep_factors()
curves  <- purrr::map_dfr(res, "curves")  |> prep_factors()
saveRDS(metrics, file.path(out_dir, "metrics_all.RDS"))


# ====================== Failure rates & capacity shortfalls ===============

fail <- metrics |>
  dplyr::group_by(strategy, ratio, n_pres) |>
  dplyr::summarise(
    n_runs    = dplyr::n(),
    fail_rate = mean(!is.na(error)),
    shortfall = mean(n_bg_real < n_bg_target),
    .groups   = "drop")
readr::write_csv(fail, file.path(out_dir, "failure_rates.csv"))
print(dplyr::filter(fail, fail_rate > 0 | shortfall > 0), n = Inf)

p_fail <- ggplot(fail, aes(ratio, strategy, fill = fail_rate)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = scales::percent(fail_rate, accuracy = 1)), size = 3) +
  scale_fill_gradient(low = "grey95", high = "firebrick", limits = c(0, 1),
                      labels = scales::percent) +
  facet_wrap(~ n_pres, labeller = labeller(n_pres = \(x) paste("n =", x))) +
  labs(x = "Ratio (P:PA)", y = NULL, fill = "Failure rate") +
  theme_bw()
ggsave(file.path(envrmt$path_figures, "03_failure_rates.png"),
       p_fail, width = 10, height = 4, dpi = 300)


# ============================ Metric figures ==============================

plot_metrics <- function(d, metric_set, file) {
  d_long <- d |>
    dplyr::filter(is.na(error)) |>
    tidyr::pivot_longer(dplyr::all_of(names(metric_set)),
                        names_to = "metric", values_to = "value") |>
    dplyr::mutate(metric = factor(metric_set[metric], levels = metric_set))
  
  d_long$family <- factor(strategy_family[names(strategy_labels)[d_long$strategy]],
                          levels = names(family_lty))
  
  # colour-blind friendly: colour + shape per strategy (filled = geographic,
  # open = environmental) + line type per strategy family
  pd <- position_dodge(width = 0.6)
  p <- ggplot(d_long, aes(ratio, value, colour = strategy, fill = strategy,
                          shape = strategy, group = strategy)) +
    stat_summary(aes(linetype = family), fun = median, geom = "line",
                 position = pd, linewidth = 0.4) +
    stat_summary(fun = median,
                 fun.min = \(x) quantile(x, 0.25),
                 fun.max = \(x) quantile(x, 0.75),
                 position = pd, size = 0.45, stroke = 0.6) +
    facet_grid(metric ~ ., scales = "free_y") +   # n_pres is fixed (250)
    scale_colour_manual(values = by_label(strategy_cols, strategy_labels)) +
    scale_fill_manual(values = by_label(strategy_fills, strategy_labels)) +
    scale_shape_manual(values = by_label(strategy_shapes, strategy_labels)) +
    scale_linetype_manual(values = family_lty) +
    guides(colour = guide_legend(nrow = 2), linetype = guide_legend(nrow = 2)) +
    labs(x = "Ratio (P:PA)", y = NULL, colour = "Strategy", fill = "Strategy",
         shape = "Strategy", linetype = "Space",
         caption = "Median and interquartile range across species × replicates") +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(file, p, width = 10, height = 8, dpi = 300)
  p
}

p_perf  <- plot_metrics(metrics, metrics_perf,
                        file.path(envrmt$path_figures, "03_metrics_performance.png"))
p_truth <- plot_metrics(metrics, metrics_true,
                        file.path(envrmt$path_figures, "03_metrics_truth.png"))

# combined figure: a) performance metrics, b) agreement with the truth
# (the legend fills the free space below b, so that the rows of b match the
# height of a)
no_cap <- labs(caption = NULL)
p_metrics <- patchwork::wrap_plots(
  A = p_perf + no_cap + labs(title = "a) Model performance"),
  B = p_truth + no_cap + labs(title = "b) Agreement with the truth"),
  C = patchwork::guide_area(),
  design = "AB\nAC"   # a over both rows; b top right; legend bottom right
) +
  patchwork::plot_layout(guides = "collect") &
  guides(colour = guide_legend(ncol = 2), linetype = guide_legend(ncol = 1)) &
  theme(legend.position = "right", legend.box = "horizontal",
        plot.title = element_text(face = "bold", size = 11))

ggsave(file.path(envrmt$path_figures, "03_metrics_combined.png"),
       p_metrics, width = 12, height = 8, dpi = 300)


# ========================= Response curve figures =========================

rc_mean <- curves |>
  dplyr::filter(n_pres == rc_n_pres) |>
  dplyr::group_by(species, strategy, ratio, variable, x) |>
  dplyr::summarise(sd = sd(pred), pred = mean(pred), true = dplyr::first(true),
                   .groups = "drop")
# ratios shown with a +/- 1 SD band across replicates (variance); all ratios
# are shown as mean lines
rc_band_ratios <- c("2:1", "1:100")

rc_mean$variable <- factor(rc_mean$variable, levels = c("bio_3", "bio_5", "bio_12", "bio_15"))
rc_var_labels <- c(bio_3  = "BIO3 Isothermality (%)",
                   bio_5  = "BIO5 Max. temperature (\u00b0C)",
                   bio_12 = "BIO12 Annual precipitation (mm)",
                   bio_15 = "BIO15 Precip. seasonality (CV)")

for (sp in levels(rc_mean$species)) {
  d <- dplyr::filter(rc_mean, species == sp)
  
  # colour-blind friendly: one-hue ramp + line type per ratio; truth as a
  # wide grey band behind the model curves
  p <- ggplot(d, aes(x, pred, colour = ratio, linetype = ratio)) +
    geom_line(data = dplyr::distinct(d, variable, x, true),
              aes(x, true), inherit.aes = FALSE,
              colour = truth_style$colour, linewidth = truth_style$linewidth) +
    geom_ribbon(data = dplyr::filter(d, ratio %in% rc_band_ratios),
                aes(x = x, ymin = pmax(pred - sd, 0), ymax = pmin(pred + sd, 1),
                    fill = ratio),
                colour = NA, alpha = 0.2, inherit.aes = FALSE) +
    geom_line(linewidth = 0.6) +
    facet_grid(strategy ~ variable, scales = "free_x",
               labeller = labeller(variable = rc_var_labels)) +
    scale_colour_manual(values = ratio_cols) +
    scale_linetype_manual(values = ratio_lty) +
    scale_fill_manual(values = ratio_cols, guide = "none") +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    labs(x = NULL, y = "Response (scaled)", colour = "Ratio (P:PA)",
         linetype = "Ratio (P:PA)", title = sp) +   # n and truth: see figure caption
    theme_bw() +
    theme(legend.position = "bottom",
          legend.key.width = unit(1.5, "cm"),   # long keys: line types readable
          panel.grid.minor = element_blank(),
          strip.text.y = element_text(angle = 0, hjust = 0),
          strip.background.y = element_blank(),
          plot.title = element_text(face = "bold"))
  
  ggsave(file.path(envrmt$path_figures, paste0("03_response_curves_", sp, ".png")),
         p, width = 10, height = 11, dpi = 300)
}


# ============================== Mixed models ==============================
# value ~ strategy * ratio + (1 | species) + (1 | presence set)
# (n_pres is fixed at 500 and therefore not a model term)
# Presence set as random effect because of the paired design.

fit_lmm <- function(metric) {
  d <- dplyr::filter(metrics, is.na(error), !is.na(.data[[metric]]))
  f <- as.formula(paste(metric,
                        "~ strategy * ratio + (1 | species) + (1 | set)"))
  lmerTest::lmer(f, data = d)
}

model_metrics <- c(names(metrics_perf), names(metrics_true))
models <- purrr::map(purrr::set_names(model_metrics), fit_lmm)
saveRDS(models, file.path(out_dir, "lmm_models.RDS"))

# ANOVA (type III, Satterthwaite)
anova_tab <- purrr::imap_dfr(models, \(m, nm)
                             as.data.frame(anova(m)) |>
                               tibble::rownames_to_column("term") |>
                               dplyr::mutate(metric = nm, .before = 1))
readr::write_csv(anova_tab, file.path(out_dir, "lmm_anova.csv"))
print(anova_tab)

# Estimated marginal means per strategy within ratio
emm_tab <- purrr::imap_dfr(models, \(m, nm)
                           emmeans::emmeans(m, ~ strategy | ratio) |>
                             as.data.frame() |>
                             dplyr::mutate(metric = nm, .before = 1))
readr::write_csv(emm_tab, file.path(out_dir, "lmm_emmeans.csv"))

# Pairwise strategy comparisons per ratio (Tukey)
pairs_tab <- purrr::imap_dfr(models, \(m, nm)
                             emmeans::emmeans(m, pairwise ~ strategy | ratio, adjust = "tukey")$contrasts |>
                               as.data.frame() |>
                               dplyr::mutate(metric = nm, .before = 1))
readr::write_csv(pairs_tab, file.path(out_dir, "lmm_pairwise_strategy.csv"))

# ---- Result tables for the paper (functions/stats_tables.R) ----
# A: strategies compared within each ratio (estimated marginal means +
#    letters; same letter = no significant difference, "a" = best)
# B: change between successive ratios within each strategy
#    (negative = decrease; for the RMSE a decrease is an improvement)
lower_better <- c("rc_rmse_mean")

tab_A <- purrr::imap_dfr(models, \(m, nm)
  table_strategy_within_ratio(m, lower_is_better = nm %in% lower_better) |>
    dplyr::mutate(metric = nm, .before = 1))
tab_B <- purrr::imap_dfr(models, \(m, nm)
  table_ratio_steps(m) |> dplyr::mutate(metric = nm, .before = 1))
readr::write_csv(tab_A, file.path(out_dir, "table_strategies_within_ratio.csv"))
readr::write_csv(tab_B, file.path(out_dir, "table_ratio_steps.csv"))

# wide versions, as in the paper (one row per metric x strategy)
strat_order <- unname(strategy_labels)
tab_A |>
  dplyr::mutate(cell = sprintf("%.3f %s", emmean, letters),
                ratio = factor(ratio, ratio_order),
                strategy = factor(strategy, strat_order)) |>
  dplyr::select(metric, strategy, ratio, cell) |>
  tidyr::pivot_wider(names_from = ratio, values_from = cell) |>
  dplyr::arrange(metric, strategy) |>
  readr::write_csv(file.path(out_dir, "table_strategies_within_ratio_wide.csv"))
tab_B |>
  dplyr::mutate(cell = sprintf("%+.3f %s", change, signif),
                strategy = factor(strategy, strat_order)) |>
  dplyr::select(metric, strategy, step, cell) |>
  tidyr::pivot_wider(names_from = step, values_from = cell) |>
  dplyr::arrange(metric, strategy) |>
  readr::write_csv(file.path(out_dir, "table_ratio_steps_wide.csv"))

# ============ Response-curve parameters (after Bazzichetto et al. 2023) ===
# For every model and variable two parameters are read from the (min-max
# scaled) modelled and true response curves:
#   - optimum: position of the maximum
#   - breadth: width of the range with response >= 0.5 (half maximum)
# Errors are expressed in % of the variable range in Romania. Across the
# replicates of each species x strategy x ratio x variable they are split
# into bias (mean error), SD (spread) and RMSE (sqrt(bias^2 + variance)).
# Flat curves (variable not used) have no optimum (NA) and a breadth of the
# whole range.

curve_par <- curves |>
  dplyr::group_by(species, set, replicate, strategy, ratio, variable) |>
  dplyr::summarise(
    rng   = diff(range(x)),
    opt_m = if (any(flat)) NA_real_ else x[which.max(pred)],
    opt_t = x[which.max(true)],
    wid_m = mean(pred >= 0.5) * rng,
    wid_t = mean(true >= 0.5) * rng,
    .groups = "drop") |>
  dplyr::mutate(err_opt = 100 * (opt_m - opt_t) / rng,
                err_wid = 100 * (wid_m - wid_t) / rng)

# true relative niche breadth per species x variable (from 01): moderator
meta_breadth <- purrr::map_dfr(levels(metrics$species), \(sp) {
  m <- readRDS(file.path(envrmt$path_virtual_species, paste0(sp, "_meta.RDS")))
  tibble::tibble(species = sp, variable = m$true_pars$variable,
                 rel_sd = m$true_pars$rel_sd)
})
curve_par <- curve_par |>
  dplyr::left_join(meta_breadth, by = c("species", "variable"))

# bias / SD / RMSE per species x strategy x ratio x variable
par_summary <- curve_par |>
  tidyr::pivot_longer(c(err_opt, err_wid), names_to = "parameter",
                      values_to = "err") |>
  dplyr::mutate(parameter = factor(parameter, c("err_opt", "err_wid"),
                                   c("Optimum", "Breadth"))) |>
  dplyr::group_by(parameter, species, strategy, ratio, variable, rel_sd) |>
  dplyr::summarise(bias = mean(err, na.rm = TRUE),
                   sd   = sd(err, na.rm = TRUE),
                   rmse = sqrt(mean(err^2, na.rm = TRUE)),
                   n_ok = sum(!is.na(err)),
                   .groups = "drop")
readr::write_csv(par_summary, file.path(out_dir, "curve_parameters_bias_sd.csv"))

# ---- table: bias vs. variance per strategy x ratio ----
# Mean over species and variables of |bias| and SD (both in % of the
# variable range). bias_share = bias^2 / (bias^2 + SD^2): share of the
# squared error that is systematic (bias) rather than random (spread).
tab_bv <- par_summary |>
  dplyr::group_by(parameter, strategy, ratio) |>
  dplyr::summarise(abs_bias   = mean(abs(bias), na.rm = TRUE),
                   # bias_share before sd, which is overwritten below
                   bias_share = mean(bias^2, na.rm = TRUE) /
                     (mean(bias^2, na.rm = TRUE) + mean(sd^2, na.rm = TRUE)),
                   sd         = mean(sd, na.rm = TRUE),
                   .groups = "drop")
readr::write_csv(tab_bv, file.path(out_dir, "table_bias_variance.csv"))

# wide: one row per parameter x strategy; cells "|bias| / SD" per ratio,
# plus the bias share at the lowest and the highest ratio
tab_bv |>
  dplyr::mutate(cell = sprintf("%.1f / %.1f", abs_bias, sd)) |>
  dplyr::select(parameter, strategy, ratio, cell) |>
  tidyr::pivot_wider(names_from = ratio, values_from = cell) |>
  dplyr::left_join(
    tab_bv |>
      dplyr::filter(ratio %in% ratio_order[c(1, length(ratio_order))]) |>
      dplyr::mutate(ratio = paste("Bias share", ratio),
                    bias_share = sprintf("%.0f %%", 100 * bias_share)) |>
      dplyr::select(parameter, strategy, ratio, bias_share) |>
      tidyr::pivot_wider(names_from = ratio, values_from = bias_share),
    by = c("parameter", "strategy")) |>
  readr::write_csv(file.path(out_dir, "table_bias_variance_wide.csv"))

# ---- figure: bias vs. spread (one line per strategy, points = ratios) ----
par_plot <- par_summary |>
  dplyr::group_by(parameter, strategy, ratio) |>
  dplyr::summarise(abs_bias = mean(abs(bias), na.rm = TRUE),
                   sd       = mean(sd, na.rm = TRUE), .groups = "drop") |>
  dplyr::arrange(parameter, strategy, ratio)

p_bv <- ggplot(par_plot, aes(sd, abs_bias, colour = strategy, fill = strategy,
                             shape = strategy, group = strategy)) +
  geom_path(linewidth = 0.5) +
  geom_point(aes(size = ratio), stroke = 0.6) +
  facet_wrap(~ parameter, scales = "free") +
  scale_colour_manual(values = by_label(strategy_cols, strategy_labels)) +
  scale_fill_manual(values = by_label(strategy_fills, strategy_labels)) +
  scale_shape_manual(values = by_label(strategy_shapes, strategy_labels)) +
  scale_size_manual(values = setNames(seq(1.2, 3.6, length.out = length(ratio_order)),
                                      ratio_order)) +
  expand_limits(x = 0, y = 0) +
  labs(x = "Spread across replicates (SD, % of variable range)",
       y = "|Bias| (% of variable range)",
       colour = "Strategy", fill = "Strategy", shape = "Strategy",
       size = "Ratio (P:PA)",
       caption = paste0("Mean over species and variables; point size increases ",
                        "with the ratio (more background). Lower left = accurate ",
                        "and precise.")) +
  theme_bw(base_size = 9) +
  theme(legend.position = "right", panel.grid.minor = element_blank(),
        plot.caption = element_text(hjust = 0))

ggsave(file.path(envrmt$path_figures, "03_curve_parameters_bias_sd.png"),
       p_bv, width = 10, height = 4.5, dpi = 300)

# ---- niche breadth as moderator (cf. Bazzichetto et al. 2023) ----
# absolute parameter error per model and variable; rel_sd centred
d_mod <- curve_par |>
  dplyr::mutate(rel_sd_c = rel_sd - mean(rel_sd))

fit_breadth <- function(resp) {
  f <- as.formula(paste0("abs(", resp, ") ~ strategy * ratio * rel_sd_c + ",
                         "(1 | species) + (1 | set) + (1 | variable)"))
  lmerTest::lmer(f, data = dplyr::filter(d_mod, !is.na(.data[[resp]])))
}
models_breadth <- list(optimum = fit_breadth("err_opt"),
                       breadth = fit_breadth("err_wid"))
saveRDS(models_breadth, file.path(out_dir, "lmm_curve_parameters.RDS"))

anova_breadth <- purrr::imap_dfr(models_breadth, \(m, nm)
                                 as.data.frame(anova(m)) |>
                                   tibble::rownames_to_column("term") |>
                                   dplyr::mutate(parameter = nm, .before = 1))
readr::write_csv(anova_breadth, file.path(out_dir, "lmm_curve_parameters_anova.csv"))
print(anova_breadth)

# slope of the error on niche breadth per strategy (does breadth moderate?)
trend_breadth <- purrr::imap_dfr(models_breadth, \(m, nm)
                                 emmeans::emtrends(m, ~ strategy, var = "rel_sd_c") |>
                                   as.data.frame() |>
                                   dplyr::mutate(parameter = nm, .before = 1))
readr::write_csv(trend_breadth, file.path(out_dir, "lmm_curve_parameters_breadth_trends.csv"))
