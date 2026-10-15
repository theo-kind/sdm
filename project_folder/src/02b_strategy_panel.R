# ================================= Notes =================================
#
# 02b_strategy_panel.R
# Illustration of the six background strategies for one species and one
# presence set: background points in geographic space (top row) and in
# environmental space (PC1/PC2, bottom row).
#
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

sp_id  <- "sp01"
n_pres <- 500
ratio  <- 1          # background points per presence (1:1)

# identical to 02_run_models.R
vars         <- c("bio_3", "bio_5", "bio_12", "bio_15")
buffer_width <- c(10000, 50000)
use_thres    <- 0.75
grid_res     <- 20
ecopa_bins   <- 30
ecopa_in     <- 0.5
seed_base    <- 2026

strategies <- c(ecopa = "EcoPA", use = "USE", rand_geo = "Random (geo)",
                rand_env = "Random (env)",
                buf_in_10 = "Buffer 10 km (in)", buf_out_10 = "Buffer 10 km (out)",
                buf_in_50 = "Buffer 50 km (in)", buf_out_50 = "Buffer 50 km (out)")
# point colours and shapes: see functions/plot_style.R (pt_style)


# ============================== Prepare data ===============================

env      <- terra::rast(file.path(envrmt$path_raster, "env_rou_30s.tif"))
tmpl     <- env[[1]]
cell_tab <- terra::as.data.frame(env[[vars]], cells = TRUE, na.rm = TRUE)
rownames(cell_tab) <- NULL   # row names = row index, used to locate cells
rng      <- sapply(cell_tab[vars], range)
xy       <- terra::xyFromCell(tmpl, cell_tab$cell)

# PCA as in 02 (basis for USE and Random env)
pca      <- prcomp(cell_tab[vars], scale. = TRUE)
pc       <- as.data.frame(pca$x[, 1:2])
pc_sf    <- sf::st_as_sf(pc, coords = c("PC1", "PC2"))
var_expl <- round(100 * summary(pca)$importance[2, 1:2])

# presences
sp_val <- terra::rast(file.path(envrmt$path_virtual_species,
                                paste0(sp_id, ".tif")))[cell_tab$cell]
set.seed(seed_base)
pres_rows <- sample_rows(which(sp_val$pa == 1), n_pres)
pres <- tibble::tibble(x = xy[pres_rows, 1], y = xy[pres_rows, 2],
                       PC1 = pc$PC1[pres_rows], PC2 = pc$PC2[pres_rows])


# ============================ Draw background ==============================

prep <- prep_background(pres_rows, cell_tab, pc, pc_sf, tmpl, rng, vars,
                        buffer_width, use_thres, ecopa_bins)

bg <- purrr::imap_dfr(strategies, \(lab, s) {
  b   <- sample_background(s, n_pres * ratio, prep, cell_tab, pc, vars,
                           grid_res, ecopa_in)
  p   <- predict(pca, b[vars])
  # EcoPA returns environmental bins, not raster cells -> no coordinates
  geo <- if (s == "ecopa") matrix(NA_real_, nrow(b), 2) else
    xy[as.integer(rownames(b)), , drop = FALSE]
  tibble::tibble(id = s, strategy = lab, x = geo[, 1], y = geo[, 2],
                 PC1 = p[, 1], PC2 = p[, 2])
})
bg$strategy <- factor(bg$strategy, levels = unname(strategies))

# display limits of the environmental space (Romanian cells + 10 %)
lim <- lapply(pc, \(v) range(v) + c(-1, 1) * 0.1 * diff(range(v)))

# share of background points outside the Romanian environmental space
bg |>
  dplyr::group_by(strategy) |>
  dplyr::summarise(n = dplyr::n(),
                   outside = mean(PC1 < lim$PC1[1] | PC1 > lim$PC1[2] |
                                    PC2 < lim$PC2[1] | PC2 > lim$PC2[2])) |>
  print()


# ================================= Figure ==================================
# Columns: geographic strategies | environmental strategies.
# Top row: geographic space, bottom row: environmental space (PC1/PC2).

groups <- list(
  geo = c(rand_geo = "Random", buf_in_10 = "Buffer 10 km (in)",
          buf_out_10 = "Buffer 10 km (out)", buf_in_50 = "Buffer 50 km (in)",
          buf_out_50 = "Buffer 50 km (out)"),
  env = c(rand_env = "Random", use = "USE", ecopa = "EcoPA")
)
group_titles <- c(geo = "Geographic sampling", env = "Environmental sampling")
pt_cols <- pt_style$cols

# background points with the short panel label per group
bg_group <- function(g) {
  dplyr::filter(bg, id %in% names(groups[[g]])) |>
    dplyr::mutate(label = factor(groups[[g]][id], levels = groups[[g]]))
}

# land mask (aggregated for plotting)
land <- terra::as.data.frame(terra::aggregate(tmpl, fact = 5, na.rm = TRUE),
                             xy = TRUE, na.rm = TRUE)

# shared scales and themes
# colour-blind friendly: colour + shape (presence: triangle, background: circle)
scale_pts <- list(
  scale_colour_manual(values = pt_cols, breaks = names(pt_cols), name = NULL),
  scale_shape_manual(values = pt_style$shapes, breaks = names(pt_cols), name = NULL)
)
lab_pc    <- labs(x = sprintf("PC1 (%d %%)", var_expl[1]),
                  y = sprintf("PC2 (%d %%)", var_expl[2]))
coord_pc  <- coord_cartesian(xlim = lim$PC1, ylim = lim$PC2, expand = FALSE)

theme_base <- theme_bw(base_size = 12) +
  theme(panel.grid = element_blank(), strip.background = element_blank(),
        legend.key = element_blank(), legend.text = element_text(size = 12),
        plot.margin = margin(4, 10, 4, 10))
theme_geo <- theme_base +
  theme(axis.title = element_blank(), axis.text = element_blank(),
        axis.ticks = element_blank(), panel.border = element_blank(),
        strip.text = element_text(size = 11, margin = margin(b = 3)),
        plot.title = element_text(face = "bold", size = 13, hjust = 0.5,
                                  margin = margin(b = 2)))
theme_env <- theme_base +
  theme(strip.text = element_blank(),
        aspect.ratio = 1)   # square PCA panels
no_y <- theme(axis.title.y = element_blank(), axis.text.y = element_blank(),
              axis.ticks.y = element_blank())

# the point legend is drawn by one panel only (avoids duplicate legends)
pt_legend <- guides(colour = guide_legend(override.aes = list(size = 2.5,
                                                              alpha = 1)),
                    shape  = guide_legend())   # merged with the colour legend
no_legend <- guides(colour = "none", shape = "none")

# ---- panels ----

geo_panel <- function(g) {
  d <- bg_group(g)
  ggplot() +
    geom_raster(data = land, aes(x, y), fill = "grey88") +
    geom_point(data = dplyr::filter(d, !is.na(x)),
               aes(x, y, colour = "Background", shape = "Background"),
               size = 0.6, alpha = 0.6) +
    geom_point(data = pres, aes(x, y, colour = "Presence", shape = "Presence"),
               size = 0.6) +
    facet_wrap(~ label, nrow = 1, drop = FALSE) +
    scale_pts + no_legend +
    coord_quickmap(expand = FALSE) +
    ggtitle(group_titles[[g]]) +
    theme_geo
}

env_panel <- function(g, first = FALSE) {
  ggplot() +
    geom_bin_2d(data = pc, aes(PC1, PC2), bins = 80) +
    scale_fill_gradient(low = "grey92", high = "grey50", trans = "log10",
                        guide = "none") +
    geom_point(data = pres, aes(PC1, PC2, colour = "Presence", shape = "Presence"),
               size = 0.6) +
    geom_point(data = bg_group(g),
               aes(PC1, PC2, colour = "Background", shape = "Background"),
               size = 0.6, alpha = 0.6) +
    facet_wrap(~ label, nrow = 1, drop = FALSE) +
    scale_pts + coord_pc + lab_pc +
    (if (first) pt_legend else no_legend) +
    theme_env + (if (first) theme() else no_y)
}

# ---- assemble ----

p <- patchwork::wrap_plots(
  geo_panel("geo"), geo_panel("env"),
  env_panel("geo", first = TRUE), env_panel("env"),
  ncol = 2, widths = c(5, 3)
) +
  patchwork::plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        plot.caption = element_text(hjust = 0, size = 8))

ggsave(file.path(envrmt$path_figures, "02b_strategy_panel.png"),
       p, width = 15, height = 4.8, dpi = 300)
