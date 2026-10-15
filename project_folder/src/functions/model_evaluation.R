# ---- Standardisation of the predictors ----
# Predictors are z-transformed with the Romania-wide mean and sd (sc) before
# fitting. Without centring, x and x^2 are almost collinear (e.g. BIO3 ranges
# 17-36), the lasso keeps only one of them and a unimodal (Gaussian) response
# cannot be fitted. The scaling is stored in the model; predict() applies it
# automatically, so all predictions can use the original units.
zscale <- function(d, sc) {
  as.data.frame(scale(d[names(sc$mu)], center = sc$mu, scale = sc$sd))
}

predict.maxnet_z <- function(object, newdata, ...) {
  class(object) <- setdiff(class(object), "maxnet_z")
  predict(object, zscale(newdata, object$zscale), ...)
}

# ---- Fit maxnet ----
fit_maxnet <- function(pres_env, bg_env, sc) {
  dat <- zscale(rbind(pres_env, bg_env), sc)
  p   <- c(rep(1, nrow(pres_env)), rep(0, nrow(bg_env)))
  mod <- maxnet::maxnet(p, dat,
                        f = maxnet::maxnet.formula(p, dat, classes = "lq"),
                        addsamplestobackground = TRUE)
  mod$zscale <- sc
  class(mod) <- c("maxnet_z", class(mod))
  mod
}

# ---- AUC (Mann-Whitney), maxTSS and Cohen's kappa against true PA ----
eval_auc_tss <- function(pred, obs) {
  n_p <- as.numeric(sum(obs == 1))
  n_a <- as.numeric(sum(obs == 0))
  auc <- (sum(rank(pred)[obs == 1]) - n_p * (n_p + 1) / 2) / (n_p * n_a)
  o   <- order(pred, decreasing = TRUE)
  tpr <- cumsum(obs[o] == 1) / n_p
  fpr <- cumsum(obs[o] == 0) / n_a
  k   <- which.max(tpr - fpr)
  # Cohen's kappa at the maxTSS threshold
  tp <- tpr[k] * n_p; fp <- fpr[k] * n_a
  fn <- n_p - tp;     tn <- n_a - fp
  n  <- n_p + n_a
  pe <- ((tp + fp) * (tp + fn) + (fn + tn) * (fp + tn)) / n^2
  kappa <- ((tp + tn) / n - pe) / (1 - pe)
  list(auc = auc, tss = tpr[k] - fpr[k], kappa = kappa, thr = pred[o][k])
}

# ---- Schoener's D between prediction and true suitability ----
schoener_d <- function(a, b) 1 - 0.5 * sum(abs(a / sum(a) - b / sum(b)))


# ---- Modelled vs. true response curves ----
# Other variables at the landscape mean; both curves min-max scaled.
# The modelled curve uses the exponential output (proportional to the
# relative occurrence rate; unlike cloglog it does not saturate).
# Near-constant curves (variable effectively unused) are set to a constant
# response of 1 (no effect of the variable) and flagged with flat = TRUE, so
# that they enter the RMSE with a real error instead of being dropped.
eval_response <- function(mod, true_pars, rng, env_means, n = 100) {
  purrr::map_dfr(names(env_means), \(v) {
    nd <- as.data.frame(as.list(env_means))[rep(1, n), ]
    nd[[v]] <- seq(rng[1, v], rng[2, v], length.out = n)
    p    <- as.vector(predict(mod, nd, type = "exponential"))
    flat <- diff(range(p)) <= 1e-3 * max(p)
    ps   <- if (flat) rep(1, n) else (p - min(p)) / diff(range(p))
    tibble::tibble(variable = v, x = nd[[v]], pred = ps, flat = flat,
                   true = true_response(true_pars, v, rng, n)$y)
  })
}


# ---- Fit one model and evaluate it completely ----
run_model <- function(pres_env, bg_env, eval_env, eval_pa, eval_suit,
                      true_pars, rng, env_means, q_trim, sc) {
  
  mod <- tryCatch(fit_maxnet(pres_env, bg_env, sc), error = function(e) e)
  if (inherits(mod, "error")) {
    return(list(metrics = tibble::tibble(error = conditionMessage(mod)),
                curves  = NULL))
  }
  
  pred  <- as.vector(predict(mod, eval_env, type = "cloglog"))
  at    <- eval_auc_tss(pred, eval_pa)
  boyce <- tryCatch(
    ecospat::ecospat.boyce(fit = pred, obs = pred[eval_pa == 1],
                           PEplot = FALSE)$cor,
    error = function(e) NA_real_)
  
  curves <- eval_response(mod, true_pars, rng, env_means)
  rc <- curves |>
    dplyr::mutate(core = x >= q_trim[1, variable] & x <= q_trim[2, variable]) |>
    dplyr::group_by(variable) |>
    dplyr::summarise(
      rmse      = sqrt(mean((pred - true)^2)),
      rmsecore  = sqrt(mean((pred[core] - true[core])^2)),
      opt_shift = if (all(flat)) NA_real_ else
        (x[which.max(pred)] - x[which.max(true)]) / (max(x) - min(x)),
      .groups   = "drop") |>
    tidyr::pivot_wider(names_from = variable,
                       values_from = c(rmse, rmsecore, opt_shift),
                       names_glue = "rc_{.value}_{variable}")
  
  metrics <- tibble::tibble(
    auc          = at$auc,
    tss          = at$tss,
    kappa        = at$kappa,
    boyce        = boyce,
    rmse         = sqrt(mean((pred - eval_suit)^2)),
    spearman     = cor(pred, eval_suit, method = "spearman"),
    schoener_d   = schoener_d(pred, eval_suit),
    prev_bias    = mean(pred >= at$thr) - mean(eval_pa),
    rc_rmse_mean     = mean(unlist(rc[grepl("^rc_rmse_bio", names(rc))])),
    rc_rmsecore_mean = mean(unlist(rc[grepl("^rc_rmsecore_bio", names(rc))])),
    error        = NA_character_
  )
  
  list(metrics = dplyr::bind_cols(metrics, rc), curves = curves)
}