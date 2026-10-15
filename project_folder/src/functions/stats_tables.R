# ---- Result tables from the mixed models (used in 03_analysis.R) ----

# Compact letter display from pairwise p-values (insert-and-absorb).
# levels: groups ordered best first; a, b: members of each pair; p: adjusted
# p-values. Groups sharing a letter do not differ significantly.
cld_letters <- function(levels, a, b, p, alpha = 0.05) {
  M <- matrix(TRUE, length(levels), 1, dimnames = list(levels, NULL))
  for (k in which(p < alpha)) {
    i <- a[k]; j <- b[k]
    for (c in which(M[i, ] & M[j, ])) {
      new <- M[, c]; M[i, c] <- FALSE; new[j] <- FALSE; M <- cbind(M, new)
    }
    keep <- rep(TRUE, ncol(M))            # absorb columns contained in others
    for (x in seq_len(ncol(M))) for (y in seq_len(ncol(M)))
      if (x != y && keep[y] && all(M[, x] <= M[, y]) &&
          (any(M[, x] != M[, y]) || x > y)) { keep[x] <- FALSE; break }
    M <- M[, keep, drop = FALSE]
  }
  M <- M[, order(apply(M, 2, \(v) which(v)[1])), drop = FALSE]
  setNames(apply(M, 1, \(v) paste(letters[which(v)], collapse = "")), levels)
}

# Table A: estimated marginal means per strategy x ratio with letters that
# compare the strategies within each ratio (Tukey); "a" = best group.
table_strategy_within_ratio <- function(model, lower_is_better = FALSE) {
  em <- as.data.frame(emmeans::emmeans(model, ~ strategy * ratio))
  ct <- as.data.frame(emmeans::emmeans(model, pairwise ~ strategy | ratio,
                                       adjust = "tukey")$contrasts)
  ct$a <- sub(" - .*", "", ct$contrast)
  ct$b <- sub(".* - ", "", ct$contrast)
  purrr::map_dfr(levels(em$ratio), \(r) {
    e <- em[em$ratio == r, ]
    e <- e[order(if (lower_is_better) e$emmean else -e$emmean), ]
    q <- ct[ct$ratio == r, ]
    l <- cld_letters(as.character(e$strategy), q$a, q$b, q$p.value)
    data.frame(strategy = as.character(e$strategy), ratio = r,
               emmean = e$emmean, SE = e$SE, letters = unname(l[as.character(e$strategy)]))
  })
}

# Table B: change between successive ratios within each strategy
# (later minus earlier ratio; Tukey-adjusted over all ratio pairs).
table_ratio_steps <- function(model) {
  ct <- as.data.frame(emmeans::emmeans(model, pairwise ~ ratio | strategy,
                                       adjust = "tukey")$contrasts)
  lv    <- levels(model@frame$ratio)
  steps <- paste(head(lv, -1), "-", tail(lv, -1))
  ct |>
    dplyr::filter(contrast %in% steps) |>
    dplyr::transmute(strategy = as.character(strategy),
                     step = factor(sub(" - ", " -> ", contrast),
                                   levels = sub(" - ", " -> ", steps)),
                     change = -estimate, SE, p.value,
                     signif = dplyr::case_when(p.value < 0.001 ~ "***",
                                               p.value < 0.01  ~ "**",
                                               p.value < 0.05  ~ "*",
                                               TRUE ~ "n.s."))
}
