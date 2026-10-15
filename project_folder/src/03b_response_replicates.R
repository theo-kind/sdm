# ================================= Notes =================================
#
# 03b_response_replicates.R
# Modelled response curves of all replicates (grey) against the true
# response (red) for one species × n_pres × ratio. Rows = strategies,
# columns = variables (layout after Bazzichetto et al. 2023).
# Shaded: outside the 2.5–97.5 % range of each variable in Romania.
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

# ================================ Parameters ===============================

sp_id     <- "sp05"
n_pres    <- 250
ratio_sel <- "1:100"

vars       <- c("bio_3", "bio_5", "bio_12", "bio_15")
var_labels <- c(bio_3 = "BIO3 (%)", bio_5 = "BIO5 (\u00b0C)",
                bio_12 = "BIO12 (mm)", bio_15 = "BIO15 (CV)")
strategies <- c(rand_geo = "Random (geo)",
                buf_in_10 = "Buffer 10 km (in)", buf_out_10 = "Buffer 10 km (out)",
                buf_in_50 = "Buffer 50 km (in)", buf_out_50 = "Buffer 50 km (out)",
                rand_env = "Random (env)", use = "USE", ecopa = "EcoPA")
core_q     <- c(0.025, 0.975)
seed_base  <- 2026

res_dir <- file.path(envrmt$path_data, "results", "runs")


# ================================ Load data ================================

files <- list.files(res_dir, full.names = TRUE,
                    pattern = sprintf("^%s_p%03d_r\\d+\\.RDS$", sp_id, n_pres))
message(length(files), " replicate files found")

curves <- purrr::map_dfr(files, \(f) readRDS(f)$curves) |>
  dplyr::filter(ratio == ratio_sel) |>
  dplyr::mutate(strategy = factor(strategies[strategy], levels = strategies),
                variable = factor(variable, levels = vars))

# true response (identical across strategies and replicates)
truth <- dplyr::distinct(curves, variable, x, true)

# number of successful fits per strategy
n_ok <- curves |>
  dplyr::distinct(strategy, replicate) |>
  dplyr::count(strategy) |>
  dplyr::mutate(variable = factor(vars[1], levels = vars),
                label = sprintf("n = %d", n))


# ======================= Core range of each variable =======================

env <- terra::rast(file.path(envrmt$path_raster, "env_rou_30s.tif"))
set.seed(seed_base)
smp <- terra::spatSample(env[[vars]], 50000, na.rm = TRUE)

core <- curves |>
  dplyr::group_by(variable) |>
  dplyr::summarise(xmin = min(x), xmax = max(x), .groups = "drop") |>
  dplyr::mutate(lo = sapply(smp[as.character(variable)], quantile, core_q[1]),
                hi = sapply(smp[as.character(variable)], quantile, core_q[2]))

tails <- dplyr::bind_rows(
  dplyr::transmute(core, variable, xmin = xmin, xmax = lo),
  dplyr::transmute(core, variable, xmin = hi,   xmax = xmax))


# ================================= Figure ==================================

p <- ggplot() +
  geom_rect(data = tails, aes(xmin = xmin, xmax = xmax,
                              ymin = -Inf, ymax = Inf), fill = "grey93") +
  geom_line(data = curves, aes(x, pred, group = replicate),
            colour = "grey25", alpha = 0.35, linewidth = 0.3) +
  geom_line(data = truth, aes(x, true), colour = "black", linetype = "dashed",
            linewidth = 0.8) +
  geom_text(data = n_ok, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.15, vjust = 1.4, size = 2.3, colour = "grey30") +
  facet_grid(strategy ~ variable, scales = "free_x",
             labeller = labeller(variable = var_labels)) +
  scale_y_continuous(breaks = c(0, 0.5, 1)) +
  labs(x = NULL, y = "Response (scaled)",
       title = sprintf("%s, n = %d presences, ratio %s",
                       sp_id, n_pres, ratio_sel),
       caption = paste0("Grey: modelled response per replicate; black dashed: true ",
                        "response. Shaded: outside the 2.5\u201397.5 % range.")) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        strip.text.y = element_text(angle = 0, hjust = 0))

ggsave(file.path(envrmt$path_figures,
                 sprintf("03b_response_%s_p%03d_%s.png", sp_id, n_pres,
                         gsub(":", "-", ratio_sel))),
       p, width = 7.5, height = 12, dpi = 300)