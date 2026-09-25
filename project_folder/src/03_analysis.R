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

strategy_labels <- c(ecopa = "EcoPA", use = "USE", rand_geo = "Random (geo)",
                     rand_env = "Random (env)", buf_in = "Buffer (in)",
                     buf_out = "Buffer (out)")
ratio_order  <- c("2:1", "1:1", "1:10", "1:100")
metrics_perf <- c(auc = "AUC", tss = "TSS", boyce = "Boyce")
metrics_true <- c(rmse = "RMSE (suitability)", schoener_d = "Schoener's D",
                  rc_rmse_mean = "RMSE (response curves)",
                  rc_rmsecore_mean = "RMSE (response curves, core)")
rc_n_pres    <- 300   # n_pres for the response curve figures

emmeans::emm_options(lmer.df = "asymptotic")  # fast for large numbers of observations


# ============================ Load results ================================

files <- list.files(res_dir, pattern = "\\.RDS$", full.names = TRUE)
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
  
  pd <- position_dodge(width = 0.6)
  p <- ggplot(d_long, aes(ratio, value, colour = strategy, group = strategy)) +
    stat_summary(fun = median, geom = "line", position = pd, linewidth = 0.4) +
    stat_summary(fun = median,
                 fun.min = \(x) quantile(x, 0.25),
                 fun.max = \(x) quantile(x, 0.75),
                 position = pd, size = 0.3) +
    facet_grid(metric ~ n_pres, scales = "free_y",
               labeller = labeller(n_pres = \(x) paste("n =", x))) +
    scale_colour_brewer(palette = "Dark2") +
    labs(x = "Ratio (P:PA)", y = NULL, colour = "Strategy",
         caption = "Median and interquartile range across species × replicates") +
    theme_bw() +
    theme(legend.position = "bottom")
  ggsave(file, p, width = 10, height = 8, dpi = 300)
  p
}

plot_metrics(metrics, metrics_perf,
             file.path(envrmt$path_figures, "03_metrics_performance.png"))
plot_metrics(metrics, metrics_true,
             file.path(envrmt$path_figures, "03_metrics_truth.png"))


# ========================= Response curve figures =========================

rc_mean <- curves |>
  dplyr::filter(n_pres == rc_n_pres) |>
  dplyr::group_by(species, strategy, ratio, variable, x) |>
  dplyr::summarise(pred = mean(pred), true = dplyr::first(true), .groups = "drop")

rc_mean$variable <- factor(rc_mean$variable, levels = c("bio_3", "bio_5", "bio_12", "bio_15"))

for (sp in levels(rc_mean$species)) {
  d <- dplyr::filter(rc_mean, species == sp)
  
  p <- ggplot(d, aes(x, pred, colour = ratio)) +
    geom_line(data = dplyr::distinct(d, variable, x, true),
              aes(x, true), inherit.aes = FALSE,
              colour = "black", linetype = "dashed", linewidth = 0.5) +
    geom_line(linewidth = 0.6) +
    facet_grid(strategy ~ variable, scales = "free_x") +
    scale_colour_viridis_d(end = 0.9) +
    labs(x = NULL, y = "Response (scaled)", colour = "Ratio (P:PA)",
         title = paste0(sp, " | n = ", rc_n_pres,
                        " | mean across replicates, dashed: true response")) +
    theme_bw() +
    theme(legend.position = "bottom")
  
  ggsave(file.path(envrmt$path_figures, paste0("03_response_curves_", sp, ".png")),
         p, width = 10, height = 11, dpi = 300)
}


# ============================== Mixed models ==============================
# value ~ strategy * ratio * n_pres + (1 | species) + (1 | presence set)
# Presence set as random effect because of the paired design.

fit_lmm <- function(metric) {
  d <- dplyr::filter(metrics, is.na(error), !is.na(.data[[metric]]))
  f <- as.formula(paste(metric,
                        "~ strategy * ratio * n_pres + (1 | species) + (1 | set)"))
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

# Estimated marginal means per strategy within ratio × n_pres
emm_tab <- purrr::imap_dfr(models, \(m, nm)
                           emmeans::emmeans(m, ~ strategy | ratio * n_pres) |>
                             as.data.frame() |>
                             dplyr::mutate(metric = nm, .before = 1))
readr::write_csv(emm_tab, file.path(out_dir, "lmm_emmeans.csv"))

# Pairwise strategy comparisons per ratio (averaged over n_pres, Tukey)
pairs_tab <- purrr::imap_dfr(models, \(m, nm)
                             emmeans::emmeans(m, pairwise ~ strategy | ratio, adjust = "tukey")$contrasts |>
                               as.data.frame() |>
                               dplyr::mutate(metric = nm, .before = 1))
readr::write_csv(pairs_tab, file.path(out_dir, "lmm_pairwise_strategy.csv"))