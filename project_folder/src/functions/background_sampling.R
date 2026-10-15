# ---- Helper functions ----

# random selection from a pool (without replacement, at most the pool size)
sample_rows <- function(pool, n) {
  pool[sample.int(length(pool), min(n, length(pool)))]
}

# Uniform sampling in PC1-PC2 space (USE logic, Da Re et al. 2023):
# grid over the candidates, random rank per grid cell, then round by round
# (1st point per cell, 2nd point per cell, ...) until n is reached
uniform_env_sample <- function(pool, n, pc, grid_res) {
  gx <- cut(pc$PC1[pool], breaks = grid_res, labels = FALSE)
  gy <- cut(pc$PC2[pool], breaks = grid_res, labels = FALSE)
  g  <- gx + (gy - 1) * grid_res
  rk <- ave(seq_along(pool), g, FUN = function(i) sample.int(length(i)))
  pool[order(rk, runif(length(pool)))[seq_len(min(n, length(pool)))]]
}


# EcoPA fallback: if all occupied niche bins hold the same number of presences
# (typical for small n), EcoPA's in-niche weights (max count - count) are all
# zero and PAGeneration() fails. In that case, draw in-niche PAs with equal
# weights (the limiting case of EcoPA's weighting); out-of-niche PAs as in EcoPA.
ecopa_fallback <- function(niche, n, ratio_in) {
  counts <- niche[[1]]
  counts[is.na(counts)] <- 0
  df <- as.data.frame.table(counts, stringsAsFactors = FALSE)
  names(df) <- c(niche[[2]], "Freq")
  df[niche[[2]]] <- lapply(df[niche[[2]]], as.numeric)
  n_in    <- ceiling(round(n * ratio_in))        # same rounding as PAGeneration()
  n_out   <- ceiling(round(n * (1 - ratio_in)))
  in_bins  <- df[df$Freq > 0, niche[[2]]]
  out_bins <- df[df$Freq == 0, niche[[2]]]
  rbind(in_bins[sample.int(nrow(in_bins), n_in, replace = TRUE), , drop = FALSE],
        out_bins[sample.int(nrow(out_bins), n_out, replace = TRUE), , drop = FALSE])
}

# ---- Preparation per presence set (once per species × n_pres × replicate) ----

prep_background <- function(pres_rows, cell_tab, pc, pc_sf, tmpl, rng, vars,
                            buffer_width, use_thres, ecopa_bins) {
  
  avail <- setdiff(seq_len(nrow(cell_tab)), pres_rows)
  
  # Buffer: cells within each buffer_width (m) of the presences
  # -> candidate pools buf_in_<km> / buf_out_<km> (e.g. buf_in_10, buf_out_50)
  r <- tmpl * NA
  r[cell_tab$cell[pres_rows]] <- 1
  buf <- list()
  for (w in buffer_width) {
    inside <- terra::values(terra::buffer(r, width = w), mat = FALSE)
    inside <- which(inside[cell_tab$cell] == 1)
    buf[[paste0("buf_in_",  w / 1000)]] <- intersect(inside, avail)
    buf[[paste0("buf_out_", w / 1000)]] <- setdiff(avail, inside)
  }

  # USE: KDE of the presences in PC space, exclude the convex hull of the core
  pres_pc <- as.matrix(pc[pres_rows, c("PC1", "PC2")])
  kde     <- ks::kde(x = pres_pc, eval.points = pres_pc)$estimate
  core    <- as.data.frame(pres_pc[kde > quantile(kde, use_thres), , drop = FALSE])
  hull    <- sf::st_convex_hull(sf::st_union(
    sf::st_as_sf(core, coords = c("PC1", "PC2"))))
  in_hull <- sf::st_intersects(hull, pc_sf)[[1]]
  
  # EcoPA: binned niche space over the full range in Romania
  niche <- EcoPA::SpeciesNiche(
    data         = cell_tab[pres_rows, vars],
    niche_border = as.vector(rng[, vars]),
    bins_sizes   = unname(signif((rng[2, vars] - rng[1, vars]) / ecopa_bins, 2))
  )
  
  c(list(
    avail    = avail,
    use_cand = setdiff(avail, in_hull),
    niche    = niche
  ), buf)
}


# ---- Draw background ----

sample_background <- function(strategy, n, prep, cell_tab, pc, vars,
                              grid_res, ecopa_in) {
  
  if (strategy == "ecopa") {
    occupied <- prep$niche[[1]][!is.na(prep$niche[[1]]) & prep$niche[[1]] > 0]
    if (length(unique(occupied)) == 1) {
      message("EcoPA fallback: equal in-niche weights (all occupied bins have equal counts)")
      bg <- ecopa_fallback(prep$niche, n, ecopa_in)
    } else {
      bg <- EcoPA::PAGeneration(data = prep$niche, nb_pa = n,
                                ratio_pa_in = ecopa_in)[[1]]
    }
    return(as.data.frame(bg[, vars]))
  }
  
  # buffer strategies (buf_in_<km>, buf_out_<km>): random draw from the pool
  if (startsWith(strategy, "buf_")) {
    if (is.null(prep[[strategy]])) stop("Unknown buffer strategy: ", strategy,
                                        " (check buffer_width)")
    return(cell_tab[sample_rows(prep[[strategy]], n), vars])
  }

  rows <- switch(strategy,
                 rand_geo = sample_rows(prep$avail, n),
                 rand_env = uniform_env_sample(prep$avail, n, pc, grid_res),
                 use      = uniform_env_sample(prep$use_cand, n, pc, grid_res),
                 stop("Unknown strategy: ", strategy)
  )
  cell_tab[rows, vars]
}