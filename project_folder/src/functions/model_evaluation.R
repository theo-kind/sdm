# ---- Fit maxnet ----
fit_maxnet <- function(pres_env, bg_env) {
  dat <- rbind(pres_env, bg_env)
  p   <- c(rep(1, nrow(pres_env)), rep(0, nrow(bg_env)))
  maxnet::maxnet(p, dat,
                 f = maxnet::maxnet.formula(p, dat, classes = "lqh"),
                 addsamplestobackground = TRUE)
}

# ---- AUC (Mann-Whitney) and maxTSS against true PA ----
eval_auc_tss <- function(pred, obs) {
  n_p <- as.numeric(sum(obs == 1))
  n_a <- as.numeric(sum(obs == 0))
  auc <- (sum(rank(pred)[obs == 1]) - n_p * (n_p + 1) / 2) / (n_p * n_a)
  o   <- order(pred, decreasing = TRUE)
  tpr <- cumsum(obs[o] == 1) / n_p
  fpr <- cumsum(obs[o] == 0) / n_a
  k   <- which.max(tpr - fpr)
  list(auc = auc, tss = tpr[k] - fpr[k], thr = pred[o][k])
}

# ---- Schoener's D between prediction and true suitability ----
schoener_d <- function(a, b) 1 - 0.5 * sum(abs(a / sum(a) - b / sum(b)))

# ---- Modelled vs. true response curves ----
# other variables at the landscape mean; both curves min-max scaled
eval_response <- function(mod, true_pars, rng, env_means, n = 100) {
  purrr::map_dfr(names(env_means), \(v) {
    nd <- as.data.frame(as.list(env_means))[rep(1, n), ]
    nd[[v]] <- seq(rng[1, v], rng[2, v], length.out = n)
    p  <- as.vector(predict(mod, nd, type = "cloglog"))
    ps <- if (diff(range(p)) > 0) (p - min(p)) / diff(range(p)) else rep(0, n)
    tibble::tibble(variable = v, x = nd[[v]], pred = ps,
                   true = true_response(true_pars, v, rng, n)$y)
  })
}

# ---- Fit one model and evaluate it completely ----
run_model <- function(pres_env, bg_env, eval_env, eval_pa, eval_suit,
                      true_pars, rng, env_means, q_trim) {
  
  mod <- tryCatch(fit_maxnet(pres_env, bg_env), error = function(e) e)
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
      opt_shift = (x[which.max(pred)] - x[which.max(true)]) / (max(x) - min(x)),
      .groups   = "drop") |>
    tidyr::pivot_wider(names_from = variable,
                       values_from = c(rmse, rmsecore, opt_shift),
                       names_glue = "rc_{.value}_{variable}")
  
  metrics <- tibble::tibble(
    auc          = at$auc,
    tss          = at$tss,
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