# ================================= Notes =================================
#
# plot_style.R
# Shared, colour-blind-friendly plot style for all figure scripts
# (sourced automatically by 00_setup.R).
# Colours were checked for colour-vision deficiency (deutan/protan/tritan);
# identity never relies on colour alone:
#   - strategies: colour + point shape (filled = geographic, open =
#     environmental) + line type per strategy family
#   - ratios:     one-hue blue ramp (light -> dark = more background)
#                 + line type
#   - presence / background: orange triangles vs. blue open circles
#
# ==========================================================================

# ---- Background strategies (order: geographic, then environmental) ----
strategy_order  <- c("rand_geo", "buf_in_10", "buf_out_10", "buf_in_50",
                     "buf_out_50", "rand_env", "use", "ecopa")
strategy_family <- c(rand_geo = "geographic", buf_in_10 = "geographic",
                     buf_out_10 = "geographic", buf_in_50 = "geographic",
                     buf_out_50 = "geographic", rand_env = "environmental",
                     use = "environmental", ecopa = "environmental")
strategy_cols   <- c(rand_geo = "#2a78d6", buf_in_10 = "#eb6834",
                     buf_out_10 = "#1baf7a", buf_in_50 = "#eda100",
                     buf_out_50 = "#e87ba4", rand_env = "#008300",
                     use = "#4a3aa7", ecopa = "#e34948")
strategy_shapes <- c(rand_geo = 21, buf_in_10 = 24, buf_out_10 = 25,
                     buf_in_50 = 22, buf_out_50 = 23,   # filled (geographic)
                     rand_env = 21, use = 24, ecopa = 22) # open (environmental)
strategy_fills  <- ifelse(strategy_family == "geographic", strategy_cols, "white")
family_lty      <- c(geographic = "solid", environmental = "dashed")

# re-key a style vector from strategy ids to the display labels of a script
by_label <- function(style, labels) setNames(unname(style[names(labels)]), labels)

# ---- Ratios (P:PA; light -> dark = more background points) ----
ratio_cols <- c("2:1" = "#74a9cf", "1:1" = "#3690c0", "1:10" = "#0570b0",
                "1:40" = "#045a8d", "1:100" = "#023858")
ratio_lty  <- c("2:1" = "dotted", "1:1" = "dotdash", "1:10" = "dashed",
                "1:40" = "longdash", "1:100" = "solid")

# ---- Presence vs. background points ----
pt_style <- list(
  cols   = c(Presence = "#D55E00", Background = "#0072B2"),
  shapes = c(Presence = 17, Background = 1)
)

# ---- True response line (drawn as a wide grey band behind the models) ----
truth_style <- list(colour = "grey70", linewidth = 1.8)
