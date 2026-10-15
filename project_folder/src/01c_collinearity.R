# ================================= Notes =================================
#
# 01c_collinearity.R
# Collinearity diagnostics for the selected predictors following
# Dormann et al. (2013, Ecography 36: 27-46, Table 1):
#   - |r| (threshold 0.7), VIF (10), tolerance (0.1), condition indices /
#     condition number (30; < 10 = weak collinearity), determinant of R,
#     variance-decomposition proportions
#   - cluster dendrogram of all 19 bioclim variables (Spearman, 1 - |r|,
#     average linkage) with the |r| = 0.7 cut and the selected variables
#     highlighted
# Same cell sample as the correlation figure in 01_generate_species.R
# (20,000 cells, set.seed(seed_base)).
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

vars      <- c("bio_3", "bio_5", "bio_12", "bio_15")
seed_base <- 2026
r_thres   <- 0.7
col_sel   <- "#2a78d6"   # selected variables (plus bold + asterisk)
col_other <- "grey40"

res_dir <- file.path(envrmt$path_data, "results")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)


# ============================ Data (as in 01) =============================

bio_all <- terra::rast(file.path(envrmt$path_raster,
                                 "climate/wc2.1_country/ROU_wc2.1_30s_bio.tif"))
names(bio_all) <- substr(names(bio_all), 11, 16)
rou     <- geodata::gadm(country = "ROU", level = 0, path = envrmt$path_raster)
bio_all <- bio_all |> terra::crop(rou) |> terra::mask(rou)

set.seed(seed_base)
smp <- terra::spatSample(bio_all, 20000, na.rm = TRUE)


# ======================= Diagnostics (Dormann et al. 2013) ================

diagnose <- function(x, method) {
  R   <- cor(x, method = method)
  e   <- eigen(R)
  vif <- diag(solve(R))
  ci  <- sqrt(max(e$values) / e$values)
  # variance-decomposition proportions on the smallest eigenvalue (Belsley)
  phi <- sweep(e$vectors^2, 2, e$values, "/")
  vd  <- (phi / rowSums(phi))[, which.max(ci)]
  max_r <- sapply(seq_along(vars), \(i) max(abs(R[i, -i])))

  per_var <- tibble::tibble(
    method = method, variable = vars,
    max_abs_r = round(max_r, 2), VIF = round(vif, 2),
    tolerance = round(1 / vif, 2), VD_smallest_eigen = round(vd, 2))
  overall <- tibble::tibble(
    method = method,
    max_abs_r = round(max(abs(R[upper.tri(R)])), 2),
    max_VIF = round(max(vif), 2),
    condition_number = round(max(ci), 1),
    determinant = round(det(R), 3))
  list(per_var = per_var, overall = overall)
}

diag_s <- diagnose(smp[vars], "spearman")
diag_p <- diagnose(smp[vars], "pearson")   # maxnet uses linear/quadratic terms

per_var <- dplyr::bind_rows(diag_s$per_var, diag_p$per_var)
overall <- dplyr::bind_rows(diag_s$overall, diag_p$overall) |>
  dplyr::mutate(threshold = "|r| < 0.7; VIF < 10; CN < 30 (weak: < 10)")

print(per_var); print(overall)
readr::write_csv(per_var, file.path(res_dir, "collinearity_per_variable.csv"))
readr::write_csv(overall, file.path(res_dir, "collinearity_overall.csv"))


# ===================== Cluster dendrogram of all variables =================

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

cm <- cor(smp, method = "spearman")
hc <- hclust(as.dist(1 - abs(cm)), method = "average")
dd <- as.dendrogram(hc)

# label text and style: selected variables blue, bold, with asterisk
dd <- dendrapply(dd, \(n) {
  if (is.leaf(n)) {
    v   <- attr(n, "label")
    sel <- v %in% vars
    attr(n, "label") <- paste0(toupper(sub("_", "", v)), " ", bio_names[[v]],
                               if (sel) " *" else "")
    attr(n, "nodePar") <- list(pch = NA, lab.col = if (sel) col_sel else col_other,
                               lab.font = if (sel) 2 else 1, lab.cex = 0.8)
  }
  n
})

png(file.path(envrmt$path_figures, "01c_collinearity_dendrogram.png"),
    width = 2200, height = 1500, res = 250)
par(mar = c(4, 1, 1, 13))
plot(dd, horiz = TRUE, xlab = "1 \u2212 |Spearman r|", axes = TRUE)
abline(v = 1 - r_thres, lty = 2, col = "grey20")
mtext(sprintf("|r| = %.1f", r_thres), side = 3, at = 1 - r_thres,
      line = -1, cex = 0.75, adj = -0.05)
dev.off()


# ================== Pairwise correlations: figure + table ==================
# A: |r| of the 6 pairs among the selected variables (Spearman and Pearson)
#    against the threshold of 0.7 (shaded: suggested range 0.5-0.7,
#    Dormann et al. 2013).
# B: |Spearman r| of every non-selected variable with each selected one;
#    the highest value per row (bold, framed) shows by which selected
#    variable it is represented.

lab <- function(v) paste0(toupper(sub("_", "", v)), " ", bio_names[v])

r_s <- cor(smp, method = "spearman")
r_p <- cor(smp, method = "pearson")

# ---- table of all pairs among the selected variables ----
pairs <- t(combn(vars, 2))
pair_tab <- tibble::tibble(
  var1 = pairs[, 1], var2 = pairs[, 2],
  spearman = round(r_s[pairs], 2),
  pearson  = round(r_p[pairs], 2)
) |>
  dplyr::mutate(pair = paste(toupper(sub("_", "", var1)), "–",
                             toupper(sub("_", "", var2))),
                max_abs = pmax(abs(spearman), abs(pearson))) |>
  dplyr::arrange(dplyr::desc(max_abs))
print(pair_tab)
readr::write_csv(dplyr::select(pair_tab, var1, var2, spearman, pearson),
                 file.path(res_dir, "collinearity_pairs_selected.csv"))

# ---- panel A ----
d_a <- pair_tab |>
  tidyr::pivot_longer(c(spearman, pearson), names_to = "method",
                      values_to = "r") |>
  dplyr::mutate(method = factor(method, c("spearman", "pearson"),
                                c("Spearman", "Pearson")),
                pair = factor(pair, rev(pair_tab$pair)))

p_a <- ggplot(d_a, aes(abs(r), pair)) +
  annotate("rect", xmin = 0.5, xmax = 0.7, ymin = -Inf, ymax = Inf,
           fill = "grey92") +
  geom_vline(xintercept = r_thres, linetype = "dashed", colour = "grey20") +
  geom_line(aes(group = pair), colour = "grey60") +
  geom_point(aes(shape = method), size = 2.6, colour = col_sel,
             fill = col_sel, stroke = 0.8) +
  geom_text(data = dplyr::mutate(pair_tab, pair = factor(pair, levels(d_a$pair))),
            aes(x = 1, y = pair, label = sprintf("%.2f / %.2f", spearman, pearson)),
            hjust = 1, size = 2.7, colour = "grey20", inherit.aes = FALSE) +
  annotate("text", x = 1, y = Inf, label = "r: Spearman / Pearson", hjust = 1,
           vjust = -0.5, size = 2.6, colour = "grey30") +
  coord_cartesian(clip = "off") +
  scale_shape_manual(values = c(Spearman = 16, Pearson = 2), name = NULL) +
  scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
  labs(x = "|r|", y = NULL, title = "A  Selected variables") +
  theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        legend.position = "bottom", plot.title = element_text(face = "bold"))

# ---- panel B ----
others <- setdiff(colnames(r_s), vars)
d_b <- tidyr::expand_grid(other = others, sel = vars) |>
  dplyr::mutate(r = abs(r_s[cbind(other, sel)])) |>
  dplyr::group_by(other) |>
  dplyr::mutate(best = r == max(r), best_r = max(r)) |>
  dplyr::ungroup() |>
  dplyr::mutate(sel   = factor(sel, vars, toupper(sub("_", "", vars))),
                other = reorder(factor(lab(other)), best_r))

p_b <- ggplot(d_b, aes(sel, other, fill = r)) +
  geom_tile(colour = "white", linewidth = 0.6) +
  geom_tile(data = dplyr::filter(d_b, best), fill = NA, colour = "grey10",
            linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", r), fontface = ifelse(best, "bold", "plain"),
                colour = r > 0.6), size = 2.6) +
  scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"),
                      guide = "none") +
  scale_fill_gradient(low = "#eef4fc", high = "#0d3570", limits = c(0, 1),
                      name = "|Spearman r|") +
  labs(x = "Selected variable", y = NULL,
       title = "B  Representation of non-selected variables") +
  theme_minimal(base_size = 9) +
  theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"),
        legend.position = "bottom")

p_pairs <- patchwork::wrap_plots(p_a, p_b, ncol = 2, widths = c(1, 1.1)) +
  patchwork::plot_annotation(
    caption = paste0("A: dashed line = threshold |r| = 0.7, shaded = suggested ",
                     "range 0.5–0.7 (Dormann et al. 2013).\n",
                     "B: framed/bold = highest |r| per row, i.e. the selected ",
                     "variable that best represents the non-selected one."),
    theme = theme(plot.caption = element_text(hjust = 0, size = 8)))

ggsave(file.path(envrmt$path_figures, "01c_pairwise_correlation.png"),
       p_pairs, width = 10, height = 5.2, dpi = 300)
