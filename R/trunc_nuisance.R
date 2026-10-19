################################################################################
# Nuisance estimation for left-truncated semi-competing risks data (AIPW-trunc).
#
# Data: all functions work on cleaned data with columns
#   Q, X1, X2, delta1, delta2, A, Z1, ..., Zp, Delta, status
# where Delta = delta1 + (1 - delta1)*delta2 indicates that T = min(T1, T2) is
# observed, and status = 0 (censored), 1 (non-terminal event) or 2 (terminal event
# without the non-terminal event) is the event type at X1.
#
# Nuisance parameters (all for the population before truncation, i.e., the starred
# variables in the paper):
#   Lambda_j(t | a, Z)   cumulative transition rates, j = 1, 2, 3
#   G(t | A, Z)          P(Q <= t | A, Z), distribution function of the truncation time
#   S_D(t | A)           P(D > t | A), survival function of the residual censoring time
#                        D = C - Q; evaluated as the left limit P(D >= t | A)
#   pi(Z)                P(A = 1 | Z), propensity score
#
# Naming conventions:
#   dat_train    data used to fit a nuisance model
#   dat_predict  data at which a fitted model is evaluated, with n_pred rows
#   weights      observation (e.g., bootstrap) weights of the rows of dat_train
#   *_fit        a fitted nuisance model, e.g., SD_fit for S_D and G_fit for G
#   fit_*        fits a nuisance model and returns *_fit
#   predict_*    returns a time-free quantity, one value per row of dat_predict,
#                e.g., the linear predictor or pi(Z)
#   eval_*       evaluates an estimated function of time, e.g., S_D or G, at given
#                time points
#   grid         sorted unique values of (Q, X1, X2) in the full data, of length m
#   Arm-specific quantities are lists over a = 0, 1, where element a + 1 is for arm a.
#   Time-varying quantities are n_pred x m matrices, with rows for the subjects in
#   dat_predict and columns for the time points in grid.
################################################################################


#' Evaluate a cumulative step function
#'
#' Evaluates the right-continuous step function that starts at 0 and jumps by
#' `jump_size[j]` at `jump_time[j]`.
#'
#' @param jump_time Numeric vector of jump times, sorted increasingly.
#' @param jump_size Numeric vector of jump sizes, of the same length as `jump_time`.
#' @param t_eval Numeric vector of time points at which to evaluate the function.
#' @return A numeric vector of the same length as `t_eval`, with
#'   `sum_{j: jump_time[j] <= t} jump_size[j]` for each `t` in `t_eval`.
#' @noRd
cum_step = function(jump_time, jump_size, t_eval){
  c(0, cumsum(jump_size))[findInterval(t_eval, jump_time) + 1]
}


#' Cumulative sums along the rows of a matrix
#'
#' @param M Numeric matrix.
#' @return A numeric matrix of the same dimension as `M`, whose `(i, j)` entry is
#'   `sum(M[i, 1:j])`.
#' @noRd
row_cumsum = function(M){
  if (ncol(M) == 1) return(M)
  t(apply(M, 1, cumsum))
}


#' Jumps along the rows of a matrix of cumulative quantities
#'
#' The inverse of `row_cumsum()`; the first column is the jump from 0.
#'
#' @param M Numeric matrix whose rows are cumulative quantities on a grid.
#' @return A numeric matrix of the same dimension as `M`, whose `(i, j)` entry is
#'   `M[i, j] - M[i, j - 1]`, with `M[i, 0] = 0`.
#' @noRd
row_jump = function(M){
  M - cbind(0, M[, -ncol(M), drop = F])
}


#' Weighted Breslow estimator for (start, stop] data
#'
#' Estimates the baseline cumulative hazard of a Cox model, where the risk set at
#' time `s` is `{start < s <= stop}`. The jump at each event time `t` is the
#' weighted number of events at `t` divided by the sum of `weights*exp(lp)` over
#' the risk set at `t`.
#'
#' `glmnet` returns only the coefficients of the penalized Cox model, so the
#' baseline cumulative hazard is computed here by plugging in its linear
#' predictors. The same function is used for the unpenalized Cox model, so both
#' models give the baseline hazard in the same form, uncentered (at `lp = 0`) and
#' consistent with `predict_lp()`.
#'
#' @param start Numeric vector of entry times; a subject is not at risk at or before
#'   its entry time.
#' @param stop Numeric vector of exit times, i.e., event or censoring times, of the
#'   same length as `start`, with `start < stop`.
#' @param event Numeric (0/1) vector, the event indicator at `stop`, of the same
#'   length as `start`.
#' @param lp Numeric vector of uncentered linear predictors of the Cox model, of the
#'   same length as `start`.
#' @param weights Nonnegative numeric vector of observation weights, of the same
#'   length as `start`; events with zero weight are ignored.
#' @return A list with two numeric vectors of the same length: the sorted distinct
#'   event times (`time`) and the jumps of the baseline cumulative hazard at these
#'   times (`jump`).
#' @noRd
breslow = function(start, stop, event, lp, weights){
  risk_score = weights*exp(lp)
  is_event = event == 1 & weights > 0
  time = sort(unique(stop[is_event]))
  n_event = as.vector(tapply(weights[is_event], factor(stop[is_event], levels = time), sum))
  n_event[is.na(n_event)] = 0
  # sum of risk_score over {v >= s} for each s, with v = start or stop
  sum_at_or_after = function(v, s){
    o = order(v)
    suffix = c(rev(cumsum(rev(risk_score[o]))), 0)
    suffix[findInterval(s, v[o], left.open = T) + 1]
  }
  # risk set at s: start < s <= stop, i.e. {stop >= s} minus {start >= s}
  at_risk = sum_at_or_after(stop, time) - sum_at_or_after(start, time)
  list(time = time, jump = ifelse(at_risk > 0, n_event/at_risk, 0))
}


#' Design matrix for the hazard models
#'
#' For `model = 'pCox'`, every variable other than `A` with more than 5 unique values
#' is expanded by natural splines with 3 degrees of freedom; otherwise each variable
#' enters linearly. The spline knots are returned, so that the same basis can be
#' applied to new data at prediction time.
#'
#' @param dat Data frame containing the variables in `vars`.
#' @param vars Character vector of covariate names.
#' @param model Character scalar, `'Cox'` or `'pCox'`.
#' @param basis Named list over `vars` from a previous call on the training data,
#'   used for prediction; each element is `NULL` (linear term) or a list with the
#'   numeric vectors `knots` and `Boundary.knots`. If `NULL`, the basis is computed
#'   from `dat`.
#' @return A list with
#'   * `X`: numeric design matrix with one row per row of `dat`, whose attribute
#'     `'var'` (character vector) gives the variable of each column;
#'   * `basis`: the named list of spline bases, in the form of the argument `basis`.
#' @noRd
make_design = function(dat, vars, model, basis = NULL){
  if (is.null(basis)) {
    basis = map(vars, function(v){
      x = dat[[v]]
      if (model == 'pCox' && v != 'A' && length(unique(x)) > 5) {
        b = splines::ns(x, df = 3)
        list(knots = attr(b, 'knots'), Boundary.knots = attr(b, 'Boundary.knots'))
      } else NULL
    })
    names(basis) = vars
  }
  X_ls = map(vars, function(v){
    x = dat[[v]]
    if (is.null(basis[[v]])) {
      out = matrix(x, ncol = 1)
    } else {
      out = splines::ns(x, knots = basis[[v]]$knots,
                        Boundary.knots = basis[[v]]$Boundary.knots)
      out = matrix(out, nrow = length(x))
    }
    colnames(out) = if (ncol(out) == 1) v else paste0(v, '_ns', 1:ncol(out))
    out
  })
  X = do.call(cbind, X_ls)
  attr(X, 'var') = unlist(map2(vars, X_ls, ~rep(.x, ncol(.y))))
  list(X = X, basis = basis)
}


#' Fit a (penalized) Cox model to (start, stop] data
#'
#' Fits a Cox model (`model = 'Cox'`) or an L2-penalized Cox model with spline terms
#' (`model = 'pCox'`, with the treatment unpenalized and the penalty chosen by 5-fold
#' cross-validation), and obtains the Breslow baseline cumulative hazard. Rows with
#' zero weight are dropped.
#'
#' @param dat_train Data frame containing the variables in `vars`.
#' @param start Numeric vector of entry times, one per row of `dat_train`.
#' @param stop Numeric vector of exit times, i.e., event or censoring times, one per
#'   row of `dat_train`, with `start < stop`.
#' @param event Numeric (0/1) vector, the event indicator at `stop`, one per row of
#'   `dat_train`.
#' @param weights Nonnegative numeric vector of observation weights, one per row of
#'   `dat_train`.
#' @param vars Character vector of covariate names.
#' @param model Character scalar, `'Cox'` or `'pCox'`.
#' @return A list with
#'   * `beta`: named numeric vector of coefficients, one per column of the design
#'     matrix (`NA` coefficients from `coxph()` are set to 0);
#'   * `var`: character vector, the variable of each coefficient;
#'   * `basis`: the spline bases, see `make_design()`;
#'   * `vars`, `model`: the arguments of the same names;
#'   * `bh`: the Breslow baseline cumulative hazard, a list with numeric vectors
#'     `time` and `jump`, see `breslow()`.
#' @noRd
fit_hazard = function(dat_train, start, stop, event, weights, vars, model = 'Cox'){
  keep = weights > 0
  dat_train = dat_train[keep, , drop = F]
  start = start[keep]; stop = stop[keep]; event = event[keep]; weights = weights[keep]
  design = make_design(dat_train, vars, model)
  X = design$X

  if (model == 'Cox') {
    fit = coxph(Surv(start, stop, event) ~ X, weights = weights, ties = 'breslow', timefix = F)
    beta = as.vector(fit$coefficients)
    beta[is.na(beta)] = 0
  }
  if (model == 'pCox') {
    # L2 penalty, the treatment is not penalized
    penalty_factor = ifelse(attr(X, 'var') == 'A', 0, 1)
    # a short lambda path keeps (start, stop] fits fast: nlambda = 20 values
    # (default 100), ending at lambda.min.ratio = 1e-3 times the largest (default
    # 1e-4 when n > p), which also stays away from the nearly unpenalized end.
    # Near the small-lambda end of the path, with spline terms the
    # partial likelihood can be flat in some directions, so glmnet may hit its
    # iteration limit and warn 'did not converge'. These lambdas are not selected
    # by CV, so the warning does not affect beta at lambda.min; since fit_hazard()
    # is called for every nuisance model, fold and bootstrap replicate, the
    # warning is muffled to avoid flooding the output. Only this warning is
    # muffled (withCallingHandlers, not suppressWarnings); all others are shown.
    cv_fit = withCallingHandlers(
      glmnet::cv.glmnet(X, Surv(start, stop, event), family = 'cox',
                        alpha = 0, weights = weights, nfolds = 5,
                        nlambda = 20, lambda.min.ratio = 1e-3,
                        penalty.factor = penalty_factor,
                        # Breslow ties, as in coxph() above and breslow()
                        cox.ties = 'breslow'),
      warning = function(cond){
        if (grepl('did not converge', conditionMessage(cond))) invokeRestart('muffleWarning')
      })
    # the muffling above relies on CV not selecting the small-lambda end of the path
    if (cv_fit$lambda.min == min(cv_fit$lambda)) {
      warning('pCox: CV selected the smallest lambda on the path, where glmnet may ',
              'not have converged; the fit is nearly unpenalized, consider ',
              'model = "Cox".', call. = F)
    }
    beta = as.vector(stats::coef(cv_fit, s = 'lambda.min'))
  }
  names(beta) = colnames(X)
  lp = as.vector(X %*% beta)
  list(beta = beta, var = attr(X, 'var'), basis = design$basis, vars = vars,
       model = model, bh = breslow(start, stop, event, lp, weights))
}


#' Linear predictor of a fitted hazard model
#'
#' Computes `X beta` without centering, so that it matches the baseline hazard
#' returned by `breslow()`. Restricting `vars` gives the contribution of a subset
#' of the covariates, e.g., the term of `T1` in the model of `Lambda_3`.
#'
#' @param fit List, a model fitted by `fit_hazard()`.
#' @param dat_predict Data frame containing all covariates in `fit$vars`.
#' @param vars Character vector, a subset of `fit$vars` whose terms are included;
#'   by default all covariates of `fit`.
#' @return A numeric vector with one linear predictor per row of `dat_predict`.
#' @noRd
predict_lp = function(fit, dat_predict, vars = fit$vars){
  X = make_design(dat_predict, fit$vars, fit$model, fit$basis)$X
  keep = attr(X, 'var') %in% vars
  as.vector(X[, keep, drop = F] %*% fit$beta[keep])
}


#' Kaplan-Meier estimator of the residual censoring distribution
#'
#' Estimates `S_D(t | A)` within each treatment group. The residual censoring time
#' `D = C - Q` is right censored by `T2 - Q`, so the observed time is `X2 - Q` with
#' event indicator `1 - delta2`.
#'
#' @param dat_train Data frame, cleaned data with at least the columns `A`, `Q`, `X2`
#'   and `delta2`.
#' @param weights Nonnegative numeric vector of observation weights, one per row of
#'   `dat_train`.
#' @return A list of length 2 over `A = 0, 1`, each a list with two numeric vectors of
#'   the same length: the jump times (`time`) and the Kaplan-Meier estimates at these
#'   times (`surv`). If an arm has no censoring, `time = 0` and `surv = 1`.
#' @noRd
fit_S_D = function(dat_train, weights){
  map(0:1, function(a){
    idx = dat_train$A == a & weights > 0
    if (sum(1 - dat_train$delta2[idx]) == 0) return(list(time = 0, surv = 1))
    km = survival::survfit(Surv(X2 - Q, 1 - delta2) ~ 1, data = dat_train[idx, ],
                           weights = weights[idx])
    list(time = km$time, surv = km$surv)
  })
}


#' Evaluate the residual censoring survival function
#'
#' Evaluates the left limit `S_D(t- | A) = P(D >= t | A)`, which is the probability
#' that a subject with residual time `t` is still uncensored at `t`.
#'
#' @param SD_fit List, a model fitted by `fit_S_D()`.
#' @param t_eval Numeric vector of subject-specific residual times, e.g., `X1 - Q`.
#' @param A Numeric (0/1) vector of treatments, of the same length as `t_eval`.
#' @return A numeric vector of the same length as `t_eval`, not bounded below.
#' @noRd
eval_S_D = function(SD_fit, t_eval, A){
  out = numeric(length(t_eval))
  for (a in 0:1){
    idx = A == a
    pos = findInterval(t_eval[idx], SD_fit[[a + 1]]$time, left.open = T)
    out[idx] = c(1, SD_fit[[a + 1]]$surv)[pos + 1]
  }
  out
}


#' Estimate the truncation distribution
#'
#' Estimates `G(t | A, Z)` by a (penalized) Cox model on the reversed time scale
#' `s = tau0 - t`: `Q` is right truncated by `T`, so `tau0 - Q` is left truncated by
#' `tau0 - T`, and the model is fitted to the data `(tau0 - X1, tau0 - Q]` among the
#' subjects with `T` observed (`Delta = 1`), with IPCW weights `1/S_D(X1 - Q | A)`.
#'
#' @param dat_train Data frame, cleaned data with at least the columns `Q`, `X1`, `A`,
#'   `Delta` and the covariates in `Zs`.
#' @param weights Nonnegative numeric vector of observation weights, one per row of
#'   `dat_train`.
#' @param SD_fit List, a model fitted by `fit_S_D()`.
#' @param Zs Character vector of covariate names.
#' @param model Character scalar, `'Cox'` or `'pCox'`.
#' @param tau0 Numeric scalar, the origin of the reversed time scale, larger than all
#'   observed times.
#' @param S_min Numeric scalar in `(0, 1)`, lower bound of `S_D` in the weights.
#' @return A list, the model fitted by `fit_hazard()` on the reversed time scale, with
#'   the additional element `tau0`.
#' @noRd
fit_G = function(dat_train, weights, SD_fit, Zs, model, tau0, S_min){
  idx = dat_train$Delta == 1
  dat = dat_train[idx, , drop = F]
  weights_G = weights[idx]/pmax(eval_S_D(SD_fit, dat$X1 - dat$Q, dat$A), S_min)
  G_fit = fit_hazard(dat, start = tau0 - dat$X1, stop = tau0 - dat$Q,
                     event = rep(1, nrow(dat)), weights = weights_G,
                     vars = c('A', Zs), model = model)
  G_fit$tau0 = tau0
  G_fit
}


#' Evaluate the truncation distribution
#'
#' Evaluates `G(t | A, Z) = exp{-exp(lp(A, Z)) sum_{q_j > t} dH_j}`, where `dH_j` are
#' the jumps of the baseline reverse time cumulative hazard at the truncation times
#' `q_j = tau0 - s_j`.
#'
#' @param G_fit List, a model fitted by `fit_G()`.
#' @param dat_predict Data frame containing `A` and the covariates of `G_fit`.
#' @param t_eval Numeric vector. If `on_grid = FALSE`, one time point per row of
#'   `dat_predict`; if `on_grid = TRUE`, a grid of time points shared by all subjects.
#' @param on_grid Logical scalar; see `t_eval`.
#' @return If `on_grid = FALSE`, a numeric vector with `G(t_eval[i] | A_i, Z_i)`, one
#'   per row of `dat_predict`. If `on_grid = TRUE`, a list with
#'   * `G`: `n_pred x length(t_eval)` numeric matrix of `G(t | A, Z)`;
#'   * `G_first_minus`: numeric vector of length `n_pred`, the value of `G` before
#'     its first jump, i.e., the estimated `P(Q < q_1 | A, Z)` with `q_1` the
#'     smallest truncation time at which the baseline hazard jumps.
#'
#'   The values are not bounded below.
#' @noRd
eval_G = function(G_fit, dat_predict, t_eval, on_grid = F){
  # truncation times q_j = tau0 - s_j, in increasing order, and the jumps dH_j
  q = rev(G_fit$tau0 - G_fit$bh$time)
  dH = rev(G_fit$bh$jump)
  total = sum(dH)
  elp = exp(predict_lp(G_fit, dat_predict))
  # G(t | A, Z) = exp{-exp(lp) sum_{q_j > t} dH_j}
  H_upper = total - cum_step(q, dH, t_eval)
  if (on_grid) {
    list(G = exp(-outer(elp, H_upper)), G_first_minus = exp(-elp*total))
  } else {
    exp(-elp*H_upper)
  }
}


#' Estimate the propensity score
#'
#' Estimates `pi(Z)` by weighted logistic regression (`model = 'logit'`) or weighted
#' gradient boosting (`model = 'gbm'`) of `A` on `Z` among the subjects with `T`
#' observed (`Delta = 1`), with the inverse probability of truncation and censoring
#' weights `1/{G(X1 | A, Z) S_D(X1 - Q | A)}`, normalized to mean 1.
#'
#' The weights are needed because the observed subjects are a biased sample of the
#' population before truncation, and `pi(Z)` is defined in that population.
#'
#' @param dat_train Data frame, cleaned data with at least the columns `Q`, `X1`, `A`,
#'   `Delta` and the covariates in `Zs`.
#' @param weights Nonnegative numeric vector of observation weights, one per row of
#'   `dat_train`.
#' @param SD_fit List, a model fitted by `fit_S_D()`.
#' @param G_fit List, a model fitted by `fit_G()`.
#' @param Zs Character vector of covariate names.
#' @param model Character scalar, `'logit'` or `'gbm'`.
#' @param S_min Numeric scalar in `(0, 1)`, lower bound of `S_D` in the weights.
#' @param G_min Numeric scalar in `(0, 1)`, lower bound of `G` in the weights.
#' @return The fitted `glm` or `gbm` object; for `gbm`, with the additional element
#'   `best_tree`, an integer scalar, the number of trees selected by
#'   cross-validation.
#' @noRd
fit_PS = function(dat_train, weights, SD_fit, G_fit, Zs, model, S_min, G_min){
  idx = dat_train$Delta == 1
  dat = dat_train[idx, , drop = F]
  weights_PS = weights[idx]/pmax(eval_G(G_fit, dat, dat$X1), G_min)/
    pmax(eval_S_D(SD_fit, dat$X1 - dat$Q, dat$A), S_min)
  weights_PS = weights_PS/mean(weights_PS)
  dat = dat[, c('A', Zs), drop = F]
  if (model == 'logit') {
    PS_fit = glm(A ~ ., data = dat, family = quasibinomial, weights = weights_PS)
  }
  if (model == 'gbm') {
    # gbm prints the cross-validation progress even when verbose = FALSE
    utils::capture.output(suppressMessages(
      PS_fit <- gbm::gbm(A ~ ., data = dat, distribution = 'bernoulli',
                         n.trees = 3000, interaction.depth = 2, shrinkage = 0.01,
                         weights = weights_PS, cv.folds = 5, n.cores = 1, verbose = F)))
    PS_fit$best_tree = gbm::gbm.perf(PS_fit, method = 'cv', plot.it = F)
  }
  PS_fit
}


#' Evaluate the propensity score
#'
#' @param PS_fit A `glm` or `gbm` object fitted by `fit_PS()`.
#' @param dat_predict Data frame containing the covariates in `Zs`.
#' @param Zs Character vector of covariate names.
#' @param model Character scalar, `'logit'` or `'gbm'`, as in `fit_PS()`.
#' @return A numeric vector with one propensity score per row of `dat_predict`, not
#'   clipped.
#' @noRd
predict_PS_trunc = function(PS_fit, dat_predict, Zs, model){
  newdata = dat_predict[, Zs, drop = F]
  if (model == 'logit') return(as.vector(predict(PS_fit, newdata = newdata, type = 'response')))
  as.vector(predict(PS_fit, newdata, n.trees = PS_fit$best_tree, type = 'response'))
}


#' Fit and evaluate all nuisance models for one fold
#'
#' Fits all nuisance models on `dat_train` and evaluates them for the subjects in
#' `dat_predict` on `grid`. The models are fitted separately, so that the transition
#' rates do not depend on `(pi, G, S_D)`: `S_D`, then `G` (with IPCW), then `pi`
#' (with IPTW and IPCW); `Lambda_1`, `Lambda_2` and `Lambda_3` by (penalized) Cox
#' models with left-truncated risk sets.
#'
#' @param dat_train Data frame, cleaned data used to fit the nuisance models.
#' @param dat_predict Data frame, cleaned data at which the nuisance models are
#'   evaluated, with `n_pred` rows.
#' @param grid Numeric vector, the sorted unique values of `(Q, X1, X2)` in the full
#'   data, of length `m`.
#' @param weights_train Nonnegative numeric vector of observation weights, one per row
#'   of `dat_train`.
#' @param Zs Character vector of covariate names.
#' @param PS_model Character scalar, the model for `pi`: `'logit'` or `'gbm'`.
#' @param trunc_model Character scalar, the model for `G`: `'Cox'` or `'pCox'`.
#' @param time_model Character scalar, the model for the transition rates: `'Cox'`
#'   or `'pCox'`.
#' @param S_min Numeric scalar in `(0, 1)`, lower bound of `S_D` (and `S`) in
#'   denominators.
#' @param G_min Numeric scalar in `(0, 1)`, lower bound of `G` in denominators.
#' @return A list with elements
#'   * `PS`: numeric vector of length `n_pred`, `pi(Z)`, not clipped;
#'   * `G`: `n_pred x m` numeric matrix, `G(t | A, Z)` at the observed treatment;
#'   * `G_first_minus`: numeric vector of length `n_pred`, `G(t | A, Z)` before its
#'     first jump, so that the jumps of `G` on the grid are
#'     `G - cbind(G_first_minus, G[, -m])`;
#'   * `L1`, `L2`: lists over `a = 0, 1` of `n_pred x m` numeric matrices,
#'     `Lambda_j(t | a, Z)`;
#'   * `lp3`: list over `a = 0, 1` of `n_pred x m` numeric matrices, the linear
#'     predictor `lp3(u, a, Z)` of `Lambda_3` given `T1 = u` for `u` in `grid`, so
#'     that `Lambda_3(t | u, a, Z) = {B3(t) - B3(u)} exp{lp3(u, a, Z)}`;
#'   * `B3`, `dB3`: `n_pred x m` numeric matrices with identical rows, the baseline
#'     cumulative hazard of `Lambda_3` on the grid and its jumps;
#'   * `SD_X1`: numeric vector of length `n_pred`, `S_D(X1 - Q | A)`;
#'   * `SD_t`: `n_pred x m` numeric matrix, `S_D(t - Q | A)`;
#'   * `SD_u`: `n_pred x m` numeric matrix, `S_D(X_tilde(t) - Q | A)`, with
#'     `X_tilde(t) = max{X1, min(X2, t)}`;
#'
#'   where `S_D` is evaluated as `P(D >= t | A)` (see `eval_S_D()`) and bounded below
#'   by `S_min`.
#' @noRd
nuisance_trunc_fold = function(dat_train, dat_predict, grid, weights_train, Zs,
                               PS_model, trunc_model, time_model,
                               S_min, G_min){
  n_pred = nrow(dat_predict)
  m = length(grid)
  tau0 = max(grid) + 1

  # residual censoring S_D(t | A): Kaplan-Meier within each treatment group
  SD_fit = fit_S_D(dat_train, weights_train)

  # truncation distribution G(t | A, Z)
  G_fit = fit_G(dat_train, weights_train, SD_fit, Zs, trunc_model, tau0, S_min)
  G_grid = eval_G(G_fit, dat_predict, grid, on_grid = T)

  # propensity score
  PS_fit = fit_PS(dat_train, weights_train, SD_fit, G_fit, Zs, PS_model, S_min, G_min)
  PS = predict_PS_trunc(PS_fit, dat_predict, Zs, PS_model)

  # Lambda1 and Lambda2: cumulative cause-specific hazards of T, left truncated by Q
  L12_fit = map(1:2, ~fit_hazard(dat_train, start = dat_train$Q, stop = dat_train$X1,
                                 event = 1*(dat_train$status == .x), weights = weights_train,
                                 vars = c('A', Zs), model = time_model))
  # L_ls[[j]][[a + 1]]: n_pred x m matrix of Lambda_j(grid | a, Z)
  L_ls = map(L12_fit, function(fit){
    base = cum_step(fit$bh$time, fit$bh$jump, grid)
    map(0:1, ~exp(predict_lp(fit, mutate(dat_predict, A = .x))) %o% base)
  })

  # Lambda3(t | t1, A, Z): cumulative hazard of T2 among subjects with delta1 = 1,
  # left truncated by X1 = T1, with T1 as a covariate (no conditional Markov assumption)
  idx3 = dat_train$delta1 == 1 & dat_train$X2 > dat_train$X1
  dat3 = mutate(dat_train[idx3, , drop = F], T1 = X1)
  L3_fit = fit_hazard(dat3, start = dat3$X1, stop = dat3$X2, event = dat3$delta2,
                      weights = weights_train[idx3], vars = c('A', Zs, 'T1'),
                      model = time_model)
  B3 = cum_step(L3_fit$bh$time, L3_fit$bh$jump, grid)
  dB3 = B3 - c(0, B3[-m])
  # the linear predictor is additive in T1: lp(u, a, z) = lp(0, a, z) + f(u)
  f_T1 = predict_lp(L3_fit, mutate(dat_predict[rep(1, m), , drop = F], T1 = grid), vars = 'T1')
  lp3 = map(0:1, ~outer(predict_lp(L3_fit, mutate(dat_predict, A = .x, T1 = 0),
                                   vars = c('A', Zs)), f_T1, '+'))

  # S_D at the time points needed by the estimators
  SD_X1 = eval_S_D(SD_fit, dat_predict$X1 - dat_predict$Q, dat_predict$A)
  t_mat = matrix(grid, n_pred, m, byrow = T)
  A_mat = rep(dat_predict$A, m)
  SD_t = matrix(eval_S_D(SD_fit, as.vector(t_mat - dat_predict$Q), A_mat), n_pred, m)
  X_u = pmax(dat_predict$X1, pmin(dat_predict$X2, t_mat))
  SD_u = matrix(eval_S_D(SD_fit, as.vector(X_u - dat_predict$Q), A_mat), n_pred, m)

  list(PS = PS,
       G = G_grid$G, G_first_minus = G_grid$G_first_minus,
       L1 = L_ls[[1]], L2 = L_ls[[2]],
       lp3 = lp3,
       B3 = matrix(B3, n_pred, m, byrow = T), dB3 = matrix(dB3, n_pred, m, byrow = T),
       SD_X1 = pmax(SD_X1, S_min), SD_t = pmax(SD_t, S_min), SD_u = pmax(SD_u, S_min))
}


#' Estimate the nuisance parameters, with or without cross-fitting
#'
#' Estimates the nuisance parameters for all subjects in `data`, either on the whole
#' sample (`crossfit = FALSE`) or by `k`-fold cross-fitting, where the nuisances of
#' the subjects in each fold are estimated from the other folds.
#'
#' @param data Data frame, cleaned data with `n` rows.
#' @param grid Numeric vector, the sorted unique values of `(Q, X1, X2)` in `data`.
#' @param weights Nonnegative numeric vector of observation weights, one per row of
#'   `data`.
#' @param Zs Character vector of covariate names.
#' @param PS_model Character scalar, the model for `pi`: `'logit'` or `'gbm'`.
#' @param trunc_model Character scalar, the model for `G`: `'Cox'` or `'pCox'`.
#' @param time_model Character scalar, the model for the transition rates: `'Cox'`
#'   or `'pCox'`.
#' @param crossfit Logical scalar; whether to use cross-fitting.
#' @param k Integer scalar, the number of folds, between 2 and `n`.
#' @param shuffle Logical scalar; if `TRUE`, the rows are shuffled before forming the
#'   folds.
#' @param PS_min Numeric scalar in `[0, 0.5)`, lower bound of the propensity score;
#'   `pi(Z)` is clipped to `[PS_min, 1 - PS_min]`.
#' @param S_min Numeric scalar in `(0, 1)`, lower bound of `S_D` (and `S`) in
#'   denominators.
#' @param G_min Numeric scalar in `(0, 1)`, lower bound of `G` in denominators.
#' @param silent Logical scalar; if `TRUE`, suppresses progress messages.
#' @return The list returned by `nuisance_trunc_fold()` for all `n` subjects, in the
#'   order of the rows of `data`, with `PS` clipped, and the additional element
#'   `runtime`, a `difftime` object.
#' @noRd
nuisance_trunc = function(data, grid, weights, Zs,
                          PS_model = 'logit', trunc_model = 'Cox', time_model = 'Cox',
                          crossfit = F, k = 5, shuffle = T,
                          PS_min = 0.1, S_min = 0.05, G_min = 0.01, silent = T){
  start_time = Sys.time()
  n = nrow(data)

  if (!crossfit) {
    if (!silent) message('Fitting nuisance models on the whole sample')
    nuis = nuisance_trunc_fold(data, data, grid, weights, Zs,
                               PS_model, trunc_model, time_model, S_min, G_min)
  } else {
    folds = cut(1:n, breaks = k, labels = F)
    if (shuffle) folds = sample(folds, length(folds))
    fold_ls = map(1:k, function(fold){
      if (!silent) message('Cross-fitting fold: ', fold)
      infold = folds == fold
      nuisance_trunc_fold(dat_train = data[!infold, , drop = F],
                          dat_predict = data[infold, , drop = F], grid,
                          weights_train = weights[!infold], Zs,
                          PS_model, trunc_model, time_model, S_min, G_min)
    })
    # put the out-of-fold predictions back in the original order
    stack_folds = function(get_piece){
      pieces = map(fold_ls, get_piece)
      if (is.matrix(pieces[[1]])) {
        out = matrix(NA_real_, n, ncol(pieces[[1]]))
        for (fold in 1:k) out[folds == fold, ] = pieces[[fold]]
      } else {
        out = numeric(n)
        for (fold in 1:k) out[folds == fold] = pieces[[fold]]
      }
      out
    }
    nuis = fold_ls[[1]]
    for (name in names(nuis)){
      if (is.list(nuis[[name]])) {
        # arm-specific nuisances, a list over a = 0, 1
        nuis[[name]] = map(1:2, function(j) stack_folds(function(f) f[[name]][[j]]))
      } else {
        nuis[[name]] = stack_folds(function(f) f[[name]])
      }
    }
  }

  nuis$PS = pmin(pmax(nuis$PS, PS_min), 1 - PS_min)
  has_NA = map_lgl(nuis, ~if (is.list(.x)) any(map_lgl(.x, anyNA)) else anyNA(.x))
  if (any(has_NA)) stop('Nuisance estimates contain NA. Please check!')
  nuis$runtime = round(Sys.time() - start_time, 2)
  nuis
}
