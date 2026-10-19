################################################################################
# AIPW-trunc estimators of Lambda_j^a(t), j = 1, 2, 3, a = 0, 1, and of the CIFs,
# for semi-competing risks data under covariate-dependent left truncation and right
# censoring, together with the IPW, naive and full data estimators.
#
# Notation follows the Appendix of the paper: omega, W (for Lambda1, Lambda2) and
# omega3, W3, V3 (for Lambda3), with the truncation augmentation terms m1-m4. The data,
# the nuisance parameters and the naming conventions are as in trunc_nuisance.R; in
# particular, nuis denotes the list returned by nuisance_trunc.
################################################################################


#' Jumps of the estimated cumulative transition rates for one arm
#'
#' Computes the jumps of `Lambda_j^a`, `j = 1, 2, 3`, on `grid[1:kmax]` by solving the
#' AIPW-trunc estimating equations at each grid point, together with the IPW and
#' naive estimators. Each estimating equation is linear in the jump, so the jump is
#' a ratio of a weighted sum over subjects (numerator) to another (denominator). The
#' AIPW-trunc numerator and denominator are the IPW ones plus augmentation terms.
#' For `Lambda_3`, the augmentation terms are computed by `trunc_aug_cpp()`.
#'
#' @param data Data frame, cleaned data with `n` rows.
#' @param grid Numeric vector, the sorted unique values of `(Q, X1, X2)` in `data`,
#'   of length `m`.
#' @param nuis List returned by `nuisance_trunc()` for the subjects in `data`.
#' @param a Integer scalar, the treatment arm, 0 or 1.
#' @param kmax Integer scalar in `1, ..., m`, the number of grid points at which the
#'   jumps are computed.
#' @param weights Nonnegative numeric vector of observation weights, one per row of
#'   `data`.
#' @param S_min Numeric scalar in `(0, 1)`, lower bound of `S(t | A, Z)` in
#'   denominators.
#' @param G_min Numeric scalar in `(0, 1)`, lower bound of `G(t | A, Z)` in
#'   denominators.
#' @return A list with elements `AIPW`, `IPW` and `naive`, each a `kmax x 3` numeric
#'   matrix of jumps (columns: `Lambda1`, `Lambda2`, `Lambda3`). A jump may be `NaN`
#'   where its denominator is 0, e.g., beyond the last event.
#' @noRd
dLambda_trunc = function(data, grid, nuis, a, kmax, weights, S_min, G_min){
  n = nrow(data)
  m = length(grid)
  k_seq = 1:kmax
  # grid indices of the observed times
  iQ = match(data$Q, grid)
  iX1 = match(data$X1, grid)
  iX2 = match(data$X2, grid)
  # value of an n x m matrix M at a subject-specific grid index j_i for each row i
  at_index = function(M, j) M[cbind(1:n, j)]
  # n x kmax matrix whose (i, k) entry is the grid index k
  k_mat = matrix(k_seq, n, kmax, byrow = T)

  # S(t | a, Z) for the arm of interest, and S(t | A, Z) at the observed treatment
  S_arm = map(1:2, ~exp(-nuis$L1[[.x]] - nuis$L2[[.x]]))
  S_obs = S_arm[[1]]
  S_obs[data$A == 1, ] = S_arm[[2]][data$A == 1, ]
  S_a = S_arm[[a + 1]]
  S_obs_f = pmax(S_obs, S_min)

  # truncation distribution at the observed treatment
  G = pmax(nuis$G, G_min)
  dG = nuis$G - cbind(nuis$G_first_minus, nuis$G[, -m, drop = F])

  # C(x) = sum_{v < x} dG(v)/{S(v) G(v)^2}, and the same with the extra factor F(v)
  integrand = dG/(S_obs_f*G^2)
  C = cbind(0, row_cumsum(integrand)[, -m, drop = F])
  CF = cbind(0, row_cumsum((1 - S_obs)*integrand)[, -m, drop = F])
  SG_Q = at_index(S_obs_f, iQ)*at_index(G, iQ)

  # m1(Q, X1, A, Z, t) on columns 1..kmax, and m2(Q, X1, A, Z)
  C_max = matrix(C[cbind(rep(1:n, kmax), as.vector(pmax(k_mat, iQ)))], n, kmax)
  m1 = (iQ >= k_mat)/SG_Q - (iX1 > k_mat)*(at_index(C, iX1) - C_max)
  m2 = (1 - at_index(S_obs, iQ))/SG_Q - (at_index(CF, iX1) - at_index(CF, iQ))

  # treatment weights 1(A = a)/pi_a(Z)
  pi_a = if (a == 1) nuis$PS else 1 - nuis$PS
  coef_A = (data$A == a)/pi_a
  Delta = data$Delta
  G_X1 = at_index(G, iX1)
  base2 = 1/G_X1 + m2

  ## Lambda1 and Lambda2
  Y1 = (iQ < k_mat) & (iX1 >= k_mat)
  at_X1 = iX1 == k_mat
  dN1 = data$delta1*at_X1
  dN2 = (1 - data$delta1)*data$delta2*at_X1
  omega = coef_A*Delta/(G[, k_seq, drop = F]*nuis$SD_X1)
  W = Delta/nuis$SD_X1*(coef_A*m1 + (1 - coef_A)*base2)
  S_minus = cbind(1, S_a[, -m, drop = F])[, k_seq, drop = F]
  dL1 = row_jump(nuis$L1[[a + 1]])
  dL2 = row_jump(nuis$L2[[a + 1]])

  num_IPW = map(list(dN1, dN2), ~colSums(weights*omega*.x))
  den_IPW = colSums(weights*omega*Y1)
  num_AIPW = map2(num_IPW, list(dL1, dL2),
                  ~.x + colSums(weights*W*S_minus*.y[, k_seq, drop = F]))
  den_AIPW = den_IPW + colSums(weights*W*S_minus)

  ## Lambda3: only time points where an event or a baseline jump occurs contribute
  is_event3 = data$delta1*data$delta2 == 1
  k3 = sort(unique(c(iX2[is_event3], which(colSums(nuis$dB3 > 0) > 0))))
  k3 = k3[k3 <= kmax]
  dL3_IPW = dL3_AIPW = numeric(kmax)

  if (length(k3) > 0) {
    # the C++ routine uses 0-based grid indices
    aug = trunc_aug_cpp(S = S_a, dL1 = dL1, lp3 = nuis$lp3[[a + 1]],
                        B3 = nuis$B3, dB3 = nuis$dB3, G = nuis$G, dG = dG,
                        iQ = iQ - 1L, iX1 = iX1 - 1L,
                        need_m = data$A == a & Delta == 1, k_target = k3 - 1L,
                        S_min = S_min, G_min = G_min)
    k3_mat = matrix(k3, n, length(k3), byrow = T)
    Y3 = data$delta1*((iX1 < k3_mat) & (iX2 >= k3_mat))
    dN3 = data$delta1*data$delta2*(iX2 == k3_mat)
    Delta3 = Delta*(1 - (1 - data$delta2)*(iX2 <= k3_mat))
    omega3 = coef_A/(G_X1*nuis$SD_t[, k3, drop = F])
    SD_u = nuis$SD_u[, k3, drop = F]
    V3 = coef_A*Delta3/SD_u
    W3 = Delta3/SD_u*(coef_A*m1[, k3, drop = F] + (1 - coef_A)*base2)

    num3_IPW = colSums(weights*omega3*dN3)
    den3_IPW = colSums(weights*omega3*Y3)
    num3 = num3_IPW + colSums(weights*(W3*aug$rhod + V3*aug$m4))
    den3 = den3_IPW + colSums(weights*(W3*aug$Rd + V3*aug$m3))
    dL3_IPW[k3] = num3_IPW/den3_IPW
    dL3_AIPW[k3] = num3/den3
  }

  ## naive Nelson-Aalen estimators with left-truncated risk sets
  weights_a = weights*(data$A == a)
  Y3_all = data$delta1*((iX1 < k_mat) & (iX2 >= k_mat))
  dN3_all = data$delta1*data$delta2*(iX2 == k_mat)
  naive = cbind(colSums(weights_a*dN1)/colSums(weights_a*Y1),
                colSums(weights_a*dN2)/colSums(weights_a*Y1),
                colSums(weights_a*dN3_all)/colSums(weights_a*Y3_all))

  list(AIPW = cbind(num_AIPW[[1]]/den_AIPW, num_AIPW[[2]]/den_AIPW, dL3_AIPW),
       IPW = cbind(num_IPW[[1]]/den_IPW, num_IPW[[2]]/den_IPW, dL3_IPW),
       naive = naive)
}


#' CIFs from the jumps of the cumulative transition rates
#'
#' Computes the CIFs of both arms, with `S(u) = exp{-Lambda_1(u) - Lambda_2(u)}`,
#' `CIF_j(t) = sum_{u <= t} S(u-) dLambda_j(u)` for `j = 1, 2`, and
#' `CIF_3(t | t1) = 1 - exp{Lambda_3(t1) - Lambda_3(t)}` for `t > t1` and 0
#' otherwise. NA or negative jumps are set to 0.
#'
#' @param dL_ls List of length 2 over `a = 0, 1` of numeric matrices of jumps, each
#'   with one row per time point in `grid_k` and 3 columns (`Lambda1`, `Lambda2`,
#'   `Lambda3`), e.g., one method of the output of `dLambda_trunc()` for each arm.
#' @param grid_k Numeric vector, sorted, the time points of the rows of the matrices
#'   in `dL_ls`.
#' @param time_interest Numeric vector of time points, sorted, at which to evaluate
#'   the CIFs.
#' @param time1_interest Numeric vector of non-terminal event times to condition on
#'   for `CIF3`, or `NULL` if `CIF3` is not requested.
#' @param estimand Character vector, one or more of `'all'`, `'CIF1'`, `'CIF2'`,
#'   `'CIF3'`.
#' @param add_contrast Logical scalar; if `TRUE`, also computes risk differences and
#'   ratios.
#' @return A data frame as described in `assemble_CIF()`.
#' @noRd
CIF_from_dLambda = function(dL_ls, grid_k, time_interest, time1_interest,
                            estimand, add_contrast){
  if (any(estimand == 'all')) estimand = c('CIF1', 'CIF2', 'CIF3')
  pos = findInterval(time_interest, grid_k)
  at = function(x, p) c(0, x)[p + 1]

  arm = map(dL_ls, function(dL){
    # negative or undefined jumps (e.g., 0/0 beyond the last event) are set to 0
    dL[is.na(dL) | dL < 0] = 0
    L = apply(dL, 2, cumsum)
    if (!is.matrix(L)) L = matrix(L, nrow = 1)
    S = exp(-L[, 1] - L[, 2])
    S_minus = c(1, S[-length(S)])
    list(L = L,
         CIF1 = cumsum(S_minus*dL[, 1]),
         CIF2 = cumsum(S_minus*dL[, 2]))
  })

  CIF12 = map(c(CIF1 = 'CIF1', CIF2 = 'CIF2'), function(j) map(arm, ~at(.x[[j]], pos)))
  CIF3 = map(time1_interest, function(t1){
    pos1 = findInterval(t1, grid_k)
    map(arm, function(x){
      L3 = x$L[, 3]
      ifelse(time_interest <= t1, 0, 1 - exp(at(L3, pos1) - at(L3, pos)))
    })
  })
  assemble_CIF(time_interest, estimand, add_contrast, CIF12, CIF3)
}


#' Arrange the CIF estimates of both arms into a data frame
#'
#' Shared by `CIF_from_dLambda()` and `full_data_est_trunc()`, so that the estimators
#' and the full data benchmark have the same columns.
#'
#' @param time_interest Numeric vector of time points at which the CIFs are
#'   evaluated.
#' @param estimand Character vector, one or more of `'all'`, `'CIF1'`, `'CIF2'`,
#'   `'CIF3'`.
#' @param add_contrast Logical scalar; if `TRUE`, also reports risk differences and
#'   ratios.
#' @param CIF12 Named list with elements `CIF1` and `CIF2`, each a list of length 2
#'   over `a = 0, 1` of numeric vectors with one value per time point in
#'   `time_interest`.
#' @param CIF3 List with one element per time point in `time1_interest`, each a list
#'   of length 2 over `a = 0, 1` as in `CIF12`; may be empty if `CIF3` is not
#'   requested.
#' @return A data frame with one row per time point in `time_interest` and columns
#'   `time` and `est_<name>_0`, `est_<name>_1` (and `est_<name>_diff`,
#'   `est_<name>_ratio` if `add_contrast = TRUE`) for each `name` in `estimand`,
#'   where the `CIF3` columns have the suffix `.T1.l` for the `l`-th time point in
#'   `time1_interest` if there are several.
#' @noRd
assemble_CIF = function(time_interest, estimand, add_contrast, CIF12, CIF3){
  if (any(estimand == 'all')) estimand = c('CIF1', 'CIF2', 'CIF3')
  out = data.frame(time = time_interest)
  add_est = function(out, name, v, suffix = ''){
    col = function(type) paste0('est_', name, '_', type, suffix)
    out[[col(0)]] = v[[1]]
    out[[col(1)]] = v[[2]]
    if (add_contrast) {
      out[[col('diff')]] = v[[2]] - v[[1]]
      out[[col('ratio')]] = v[[2]]/v[[1]]
    }
    out
  }
  for (j in intersect(c('CIF1', 'CIF2'), estimand)){
    out = add_est(out, j, CIF12[[j]])
  }
  if ('CIF3' %in% estimand) {
    for (l in seq_along(CIF3)){
      suffix = if (length(CIF3) == 1) '' else paste0('.T1.', l)
      out = add_est(out, 'CIF3', CIF3[[l]], suffix)
    }
  }
  out
}


#' Full data estimator for left-truncated semi-competing risks
#'
#' \code{full_data_est_trunc()} computes the full data estimators of the treatment-specific
#' risks from the potential outcomes \eqn{(T_1^a, T_2^a)}, \eqn{a = 0, 1}, of all subjects
#' in the population before truncation and censoring. They are unattainable in practice
#' and serve as a benchmark in simulations, e.g., with the full data \code{dat_full}
#' returned by \code{\link{sim_data_trunc}}.
#'
#' With \eqn{N} subjects in the full data, the estimators are the empirical probabilities
#' \deqn{\widehat{\mathrm{CIF}}_1^a(t) = N^{-1}\sum_i 1(T_{1i}^a \le t), \quad
#'       \widehat{\mathrm{CIF}}_2^a(t) = N^{-1}\sum_i 1(T_{1i}^a = \infty, T_{2i}^a \le t),}
#' and, under the marginal Markov assumption,
#' \eqn{\widehat{\mathrm{CIF}}_3^a(t \mid t_1) = 1 - \exp\{\hat\Lambda_3^a(t_1) - \hat\Lambda_3^a(t)\}}
#' for \eqn{t > t_1} (and 0 otherwise), where \eqn{\hat\Lambda_3^a} is the Nelson-Aalen
#' estimator with at-risk indicator \eqn{1(T_1^a < u \le T_2^a)}, i.e., the solution of the
#' full data estimating equation for \eqn{\Lambda_3^a}. This is the same functional of
#' \eqn{\Lambda_3^a} as in \code{\link{causalAIPWtrunc}}.
#'
#' @param data A data.frame with the potential outcomes, one row per subject in the full
#'   data. The non-terminal event time is \code{Inf} if the terminal event occurs first.
#' @param T1_0 Character scalar giving the column name of the potential non-terminal
#'   event time under \eqn{a = 0}.
#' @param T2_0 Character scalar giving the column name of the potential terminal event
#'   time under \eqn{a = 0}.
#' @param T1_1 Character scalar giving the column name of the potential non-terminal
#'   event time under \eqn{a = 1}.
#' @param T2_1 Character scalar giving the column name of the potential terminal event
#'   time under \eqn{a = 1}.
#' @param estimand Character vector, one or more of \code{"all"}, \code{"CIF1"},
#'   \code{"CIF2"}, \code{"CIF3"}.
#' @param time_interest Numeric vector of time points at which to evaluate the
#'   estimands; nonpositive values are dropped.
#' @param time1_interest Numeric vector of non-terminal event times to condition on.
#'   Required for \code{"CIF3"}.
#' @param add_contrast Logical scalar; if \code{TRUE}, also reports risk differences
#'   and ratios.
#' @param name Optional character scalar, the method label; if given, it is added as
#'   the first column \code{method}, so that the output can be combined with the
#'   \code{est} component of \code{\link{causalAIPWtrunc}}.
#'
#' @return A data frame with the same columns as the \code{est} component of
#'   \code{\link{causalAIPWtrunc}}: \code{time} and \code{est_CIF*}, and \code{method} if
#'   \code{name} is given.
#'
#' @examples
#' sim = sim_data_trunc(n = 300)
#' full_data_est_trunc(sim$dat_full, time_interest = c(2, 3, 4),
#'                     time1_interest = 1.5)
#' @export
full_data_est_trunc = function(data, T1_0 = 'T1_0', T2_0 = 'T2_0',
                               T1_1 = 'T1_1', T2_1 = 'T2_1',
                               estimand = 'all', time_interest,
                               time1_interest = NULL,
                               add_contrast = F, name = NULL){
  if (!is.data.frame(data)) stop('Input dataset should be a dataframe!')
  vars = c(T1_0, T2_0, T1_1, T2_1)
  if (any(!vars %in% names(data))) {
    stop('Columns not found in data: ', paste(setdiff(vars, names(data)), collapse = ', '))
  }
  if (any(!estimand %in% c('all', 'CIF1', 'CIF2', 'CIF3'))) {
    stop('The estimand should be one or some of the following: all, CIF1, CIF2, CIF3')
  }
  if (any(estimand == 'all')) estimand = c('CIF1', 'CIF2', 'CIF3')
  if (is.null(time1_interest) && 'CIF3' %in% estimand) {
    stop('Please specify time to non-terminal event to condition on!')
  }
  time_interest = sort(time_interest[time_interest > 0])
  if (length(time_interest) == 0) stop('time_interest should contain positive time points!')
  time1_interest = if ('CIF3' %in% estimand) sort(time1_interest) else NULL

  arm = map(list(vars[1:2], vars[3:4]), function(v){
    T1 = data[[v[1]]]
    T2 = data[[v[2]]]
    ill = is.finite(T1)
    if (anyNA(c(T1, T2)) || any(T1[ill] > T2[ill])) {
      stop('Potential outcomes should be non-missing with T1 <= T2 (T1 = Inf without the non-terminal event)!')
    }
    # Nelson-Aalen for Lambda3 with risk set T1 < u <= T2
    is3 = ill & T2 > T1
    ev = sort(T2[is3])
    u = unique(ev)
    d = tabulate(match(ev, u), nbins = length(u))
    at_risk = findInterval(u, sort(T1[is3]), left.open = T) -
      findInterval(u, ev, left.open = T)
    L3 = c(0, cumsum(d/at_risk))
    list(CIF1 = map_dbl(time_interest, ~mean(ill & T1 <= .x)),
         CIF2 = map_dbl(time_interest, ~mean(!ill & T2 <= .x)),
         L3 = function(x) L3[findInterval(x, u) + 1])
  })

  CIF12 = map(c(CIF1 = 'CIF1', CIF2 = 'CIF2'), function(j) map(arm, ~.x[[j]]))
  CIF3 = map(time1_interest, function(t1){
    map(arm, ~ifelse(time_interest <= t1, 0, 1 - exp(.x$L3(t1) - .x$L3(time_interest))))
  })
  out = assemble_CIF(time_interest, estimand, add_contrast, CIF12, CIF3)
  if (!is.null(name)) out = cbind(method = name, out)
  out
}


#' Sanity checks of the cleaned data before estimation
#'
#' Checks that there are no missing values, that `delta1`, `delta2` and `A` are 0/1,
#' that `0 <= Q < X1 <= X2` with `X1 = X2` when `delta1 = 0`, and that both
#' treatment arms are observed.
#'
#' @param data Data frame, cleaned data with columns `Q`, `X1`, `X2`, `delta1`,
#'   `delta2`, `A` and the covariates.
#' @return `NULL` invisibly; stops with an informative message if the data are
#'   invalid.
#' @noRd
check_trunc_data = function(data){
  if (anyNA(data)) stop('data contains NA in the variables used for estimation!')
  for (v in c('delta1', 'delta2', 'A')){
    if (any(!data[[v]] %in% c(0, 1))) stop(v, ' should be coded as 0/1!')
  }
  if (any(data$Q < 0)) stop('The truncation time Q should be nonnegative!')
  if (any(data$Q >= data$X1)) stop('Every subject must have Q < X1!')
  if (any(data$X1 > data$X2)) stop('Every subject must have X1 <= X2!')
  if (any(data$delta1 == 0 & data$X1 != data$X2)) {
    stop('X1 should equal X2 for subjects without the non-terminal event (delta1 = 0)!')
  }
  if (any(table(factor(data$A, levels = 0:1)) == 0)) stop('Both treatment arms should be observed!')
  invisible(NULL)
}


#' Point estimates from cleaned data
#'
#' Estimates the nuisance parameters, the jumps of `Lambda_j^a` on the grid and the
#' CIFs, for each requested method. Called by `causalAIPWtrunc()` on the data and on
#' each bootstrap sample.
#'
#' All arguments are as in `causalAIPWtrunc()`, after its argument checks and data
#' cleaning.
#'
#' @param data Data frame, cleaned data, with the columns `Delta` and `status` added.
#' @param Zs Character vector of covariate names in `data`, i.e., `Z1, ..., Zp`.
#' @param estimand Character vector, one or more of `'all'`, `'CIF1'`, `'CIF2'`,
#'   `'CIF3'`.
#' @param time_interest Numeric vector of positive time points, sorted, or `NULL` for
#'   all observed time points.
#' @param time1_interest Numeric vector of non-terminal event times, sorted, or `NULL`
#'   if `CIF3` is not requested.
#' @param PS_model Character scalar, `'logit'` or `'gbm'`.
#' @param trunc_model Character scalar, `'Cox'` or `'pCox'`.
#' @param time_model Character scalar, `'Cox'` or `'pCox'`.
#' @param crossfit Logical scalar (`'auto'` already resolved).
#' @param k Integer scalar, the number of folds.
#' @param weights Nonnegative numeric vector of observation weights, one per row of
#'   `data`.
#' @param shuffle Logical scalar.
#' @param PS_min Numeric scalar in `[0, 0.5)`.
#' @param S_min Numeric scalar in `(0, 1)`.
#' @param G_min Numeric scalar in `(0, 1)`.
#' @param add_IPW Logical scalar.
#' @param add_naive_est Logical scalar.
#' @param name Character vector of length 2, the labels of the AIPW-trunc and IPW
#'   methods.
#' @param add_trans Logical scalar; if `TRUE`, also returns the transformed estimates
#'   used for the bootstrap confidence intervals (see `transform_est()`).
#' @param silent Logical scalar.
#' @param add_contrast Logical scalar.
#' @return A list with
#'   * `est`: data frame as described in `assemble_CIF()`, with the additional first
#'     column `method` and one block of rows per method;
#'   * `dLambda`: named list over methods of data frames with columns `time` and
#'     `dLambda<j>_<a>`, the jumps of the cumulative transition rates on
#'     `grid[1:kmax]`;
#'   * `nuisance`: the nuisance estimates, see `nuisance_trunc()`;
#'   * `runtime`: a `difftime` object.
#' @noRd
causalAIPWtrunc_est = function(data, Zs, estimand, time_interest, time1_interest,
                               PS_model, trunc_model, time_model,
                               crossfit, k, weights, shuffle,
                               PS_min, S_min, G_min,
                               add_IPW, add_naive_est, name,
                               add_trans, silent, add_contrast){
  start_time = Sys.time()
  grid = sort(unique(c(data$Q, data$X1, data$X2)))
  if (is.null(time_interest)) time_interest = grid[grid > 0]
  # only time points up to the largest time point of interest are needed
  kmax = max(1, findInterval(max(c(time_interest, time1_interest)), grid))

  nuis = nuisance_trunc(data, grid, weights, Zs, PS_model, trunc_model, time_model,
                        crossfit, k, shuffle, PS_min, S_min, G_min, silent)
  if (!silent) message('Estimate nuisance parameters: completed in ',
                       format(nuis$runtime))

  dL = map(0:1, ~dLambda_trunc(data, grid, nuis, .x, kmax, weights, S_min, G_min))
  methods = c('AIPW', if (add_IPW) 'IPW', if (add_naive_est) 'naive')
  labels = c(name[1], if (add_IPW) name[2], if (add_naive_est) 'naive')

  est = map2(methods, labels, function(method, label){
    dL_ls = map(dL, ~.x[[method]])
    if (method == 'AIPW' && anyNA(unlist(dL_ls))) {
      warning('Method name: ', label, ', NA estimates of dLambda have occurred! Please check!')
    }
    out = CIF_from_dLambda(dL_ls, grid[1:kmax], time_interest, time1_interest,
                           estimand, add_contrast)
    if (add_trans) out = transform_est(out)
    cbind(method = label, out)
  }) %>% bind_rows()

  dLambda = map(setNames(methods, labels), function(method){
    out = data.frame(time = grid[1:kmax])
    for (a in 0:1){
      d = dL[[a + 1]][[method]]
      for (j in 1:3) out[[paste0('dLambda', j, '_', a)]] = d[, j]
    }
    out
  })

  list(est = est, dLambda = dLambda, nuisance = nuis,
       runtime = round(Sys.time() - start_time, 2))
}


#' Causal estimation for left-truncated semi-competing risks via AIPW-trunc
#'
#' \code{causalAIPWtrunc()} estimates treatment-specific risks in semi-competing risks
#' data subject to covariate-dependent left truncation and right censoring. Subjects are
#' observed only if the truncation time \code{Q} is smaller than
#' \eqn{T = \min(T_1, T_2)}. The censoring time is \eqn{C = Q + D}, where the residual
#' censoring time \eqn{D} may depend on the treatment. The estimators augment inverse
#' probability of treatment, truncation and censoring weights, and are doubly robust with
#' respect to the two sets of nuisance models: (i) the cumulative transition rates of
#' \eqn{(T_1, T_2)} given \eqn{(A, Z)}, and (ii) the propensity score together with
#' the truncation distribution \eqn{G(t \mid A, Z) = P(Q \le t \mid A, Z)} and the
#' residual censoring distribution \eqn{S_D(t \mid A)}.
#'
#' The nuisance models are estimated separately so that one set does not rely on the
#' other, where "Cox" below refers to the Cox or penalized Cox model chosen by
#' \code{time_model} or \code{trunc_model}: \eqn{\Lambda_1} and \eqn{\Lambda_2} by Cox
#' models with left-truncated risk sets \eqn{(Q, X_1]}; \eqn{\Lambda_3(t \mid t_1, A, Z)}
#' by a Cox model among subjects with \code{delta1 = 1}, with entry time \code{X1}, exit
#' time \code{X2} and \eqn{T_1} as an additional covariate; \eqn{S_D} by Kaplan-Meier
#' within each treatment group; \eqn{G} by a Cox model on the reversed time scale with
#' weights \eqn{\Delta / S_D(X_1 - Q \mid A)}; and the propensity score by regressing
#' \eqn{A} on \eqn{Z} (model \code{PS_model}) with weights
#' \eqn{\Delta / \{G(X_1 \mid A, Z) S_D(X_1 - Q \mid A)\}}, where \eqn{\Delta} indicates that
#' \eqn{T = \min(T_1, T_2)} is observed. With cross-fitting, all nuisance models of the
#' subjects in each fold are estimated from the other folds.
#'
#' Supported estimands are
#' \itemize{
#'   \item \code{"CIF1"}: risk of the non-terminal event,
#'   \item \code{"CIF2"}: risk of the terminal event without the non-terminal event,
#'   \item \code{"CIF3"}: risk of the terminal event given the non-terminal event at
#'         \code{time1_interest}, under the marginal Markov assumption,
#'   \item \code{"all"}: all of the above.
#' }
#'
#' @param data A data.frame containing the variables specified by \code{Q}, \code{X1},
#'   \code{X2}, \code{delta1}, \code{delta2}, \code{treatment} and \code{covars}.
#' @param Q Character scalar giving the column name of the left truncation time. Every
#'   subject must have \code{0 <= Q < X1}.
#' @param X1 Character scalar giving the column name of the observed time to the first
#'   event, \eqn{\min(T_1, T_2, C)}.
#' @param X2 Character scalar giving the column name of the observed time to the
#'   terminal event, \eqn{\min(T_2, C)}.
#' @param delta1 Character scalar giving the column name of the 0/1 indicator that the
#'   non-terminal event is observed at \code{X1}.
#' @param delta2 Character scalar giving the column name of the 0/1 indicator that the
#'   terminal event is observed at \code{X2}.
#' @param treatment Character scalar giving the column name of the 0/1 treatment.
#' @param covars Character vector of baseline covariate column names used in the
#'   nuisance models.
#' @param estimand Character vector, one or more of \code{"all"}, \code{"CIF1"},
#'   \code{"CIF2"}, \code{"CIF3"}.
#' @param time_interest Optional numeric vector of time points at which to evaluate the
#'   estimands; nonpositive values are dropped. If \code{NULL}, all unique observed
#'   time points are used.
#' @param time1_interest Numeric vector of non-terminal event times to condition on.
#'   Required for \code{"CIF3"}.
#' @param PS_model Character scalar, the propensity score model: \code{"logit"}
#'   (weighted logistic regression) or \code{"gbm"} (weighted gradient boosting).
#' @param trunc_model Character scalar, the model for the truncation distribution
#'   \eqn{G}: \code{"Cox"} or \code{"pCox"} (L2-penalized Cox regression with natural
#'   spline basis functions of the continuous covariates, tuned by cross-validation).
#' @param time_model Character scalar, the model for the cumulative transition rates:
#'   \code{"Cox"} or \code{"pCox"}.
#' @param k Integer scalar, the number of folds for cross-fitting, between 2 and
#'   \code{nrow(data)}.
#' @param crossfit Logical scalar or \code{"auto"} (default). \code{"auto"} cross-fits
#'   all nuisance models when any of them is \code{"gbm"} or \code{"pCox"}.
#' @param weights Optional nonnegative numeric vector of observation weights, one per
#'   row of \code{data}; if \code{NULL}, all weights are 1.
#' @param shuffle Logical scalar; if \code{TRUE}, rows are shuffled before forming
#'   folds.
#' @param PS_min Numeric scalar in \eqn{[0, 0.5)}, lower bound of the propensity score;
#'   scores are clipped to \code{[PS_min, 1 - PS_min]}.
#' @param S_min Numeric scalar in \eqn{(0, 1)}, lower bound for \eqn{S_D} and for
#'   \eqn{S(t \mid A, Z)} in denominators.
#' @param G_min Numeric scalar in \eqn{(0, 1)}, lower bound for \eqn{G(t \mid A, Z)} in
#'   denominators.
#' @param add_IPW Logical scalar; if \code{TRUE}, also reports the IPW estimators.
#' @param add_naive_est Logical scalar; if \code{TRUE}, also reports the naive
#'   Nelson-Aalen estimators with left-truncated risk sets.
#' @param n_boot Integer scalar, the number of bootstrap replicates; \code{0} skips
#'   inference.
#' @param boot_raw Logical scalar; if \code{TRUE}, returns the raw bootstrap replicates.
#' @param alpha Numeric scalar in \eqn{(0, 1)}, the significance level of the
#'   confidence intervals.
#' @param seed Optional integer scalar, the random seed.
#' @param name Optional character vector of length 1 or 2, the method label(s),
#'   defaulting to \code{"AIPWtrunc"} and \code{"IPW"}.
#' @param CI_raw Logical scalar; if \code{TRUE}, also returns Wald intervals on the
#'   original scale.
#' @param silent Logical scalar; if \code{TRUE}, suppresses progress messages.
#' @param add_contrast Logical scalar; if \code{TRUE}, also reports risk differences
#'   and ratios.
#' @param parallel Logical scalar; if \code{TRUE}, runs the bootstrap in parallel via
#'   \pkg{future.apply}.
#' @param n_cores Optional integer scalar, the number of workers when
#'   \code{parallel = TRUE}; if \code{NULL}, all detected cores are used.
#'
#' @return A list with components
#' \describe{
#'   \item{\code{est}}{Point estimates with columns \code{time} and \code{est_CIF*}
#'     (and \code{method} when several methods are requested).}
#'   \item{\code{est_se}}{Bootstrap standard errors and confidence intervals in long format
#'     (or the raw replicates when \code{boot_raw = TRUE}); \code{NULL} if \code{n_boot = 0}.}
#'   \item{\code{dLambda}}{Jumps of the estimated cumulative transition rates on the
#'     observed time grid, one data frame per method.}
#'   \item{\code{runtime}}{Runtime of the point estimation.}
#'   \item{\code{call}}{The matched call.}
#' }
#'
#' @examples
#' \dontrun{
#' dat = sim_data_trunc(n = 500)$dat
#' fit = causalAIPWtrunc(dat, covars = c('Z1', 'Z2'),
#'                       time_interest = c(2, 3, 4), time1_interest = 1.5,
#'                       add_IPW = TRUE, add_naive_est = TRUE, n_boot = 0)
#' fit$est
#' }
#' @export
causalAIPWtrunc = function(data, Q = 'Q', X1 = 'X1', X2 = 'X2',
                           delta1 = 'delta1', delta2 = 'delta2',
                           treatment = 'A', covars = 'Z',
                           estimand = 'all',
                           time_interest = NULL,
                           time1_interest = NULL,
                           PS_model = 'logit', trunc_model = 'Cox', time_model = 'Cox',
                           k = 5, crossfit = 'auto', weights = NULL, shuffle = T,
                           PS_min = 0.1, S_min = 0.05, G_min = 0.01,
                           add_IPW = F, add_naive_est = F,
                           n_boot = 100, boot_raw = F, alpha = 0.05,
                           seed = NULL, name = NULL, CI_raw = F, silent = T,
                           add_contrast = F, parallel = F, n_cores = NULL){
  if (!is.null(seed)) set.seed(seed)
  if (!is.data.frame(data)) stop('Input dataset should be a dataframe!')
  check_option = function(value, arg, choices){
    if (length(value) != 1 || !value %in% choices) {
      stop(arg, ' should be one of the following: ', paste(choices, collapse = ', '))
    }
  }
  check_option(PS_model, 'PS_model', c('logit', 'gbm'))
  check_option(trunc_model, 'trunc_model', c('Cox', 'pCox'))
  check_option(time_model, 'time_model', c('Cox', 'pCox'))
  if (any(!estimand %in% c('all', 'CIF1', 'CIF2', 'CIF3'))) {
    stop('The estimand should be one or some of the following: all, CIF1, CIF2, CIF3')
  }
  if (is.null(time1_interest) && any(estimand %in% c('all', 'CIF3'))) {
    stop('Please specify time to non-terminal event to condition on!')
  }
  vars = c(Q, X1, X2, delta1, delta2, treatment, covars)
  if (any(!vars %in% names(data))) {
    stop('Columns not found in data: ', paste(setdiff(vars, names(data)), collapse = ', '))
  }
  n = nrow(data)
  if (is.null(weights)) weights = rep(1, n)
  if (length(weights) != n || anyNA(weights) || any(weights < 0)) {
    stop('weights should be a nonnegative vector with one entry per row of data!')
  }
  if (crossfit == 'auto') crossfit = any(c(PS_model, trunc_model, time_model) %in% c('gbm', 'pCox'))
  if (crossfit && (length(k) != 1 || k < 2 || k > n)) stop('k should be an integer between 2 and nrow(data)!')

  specify_name = !is.null(name)
  # fill in default labels, e.g., name = 'myAIPW' with add_IPW = TRUE gives c('myAIPW', 'IPW')
  if (length(name) < 2) name = c(name, c('AIPWtrunc', 'IPW')[(length(name) + 1):2])
  add_trans = !boot_raw && n_boot > 0

  # clean data
  Zs = paste0('Z', seq_along(covars))
  data = data %>%
    dplyr::select(all_of(vars)) %>%
    setNames(c('Q', 'X1', 'X2', 'delta1', 'delta2', 'A', Zs)) %>%
    as.data.frame()
  check_trunc_data(data)
  data = data %>%
    mutate(Delta = delta1 + (1 - delta1)*delta2,
           # status: 0=censored, 1=non-terminal event, 2=terminal event without non-terminal event
           status = ifelse(delta1 == 1, 1, 2*delta2))
  if (!is.null(time_interest)) {
    time_interest = sort(time_interest[time_interest > 0])
    if (length(time_interest) == 0) stop('time_interest should contain positive time points!')
  }
  if (!is.null(time1_interest)) time1_interest = sort(time1_interest)

  est_fun = function(dat, dat_weights){
    causalAIPWtrunc_est(dat, Zs, estimand, time_interest, time1_interest,
                        PS_model, trunc_model, time_model,
                        crossfit, k, dat_weights, shuffle, PS_min, S_min, G_min,
                        add_IPW, add_naive_est, name, add_trans, silent, add_contrast)
  }
  fit = est_fun(data, weights)
  est = fit$est

  est_se = NULL
  if (n_boot > 0) {
    boot_one = function(i){
      boot_index = sample(1:n, n, replace = T)
      cbind(boot_id = i, est_fun(data[boot_index, ], weights[boot_index])$est)
    }
    if (parallel) {
      if (is.null(n_cores)) n_cores = parallel::detectCores()
      old_plan = future::plan()
      on.exit(future::plan(old_plan), add = TRUE)
      future::plan('multisession', workers = n_cores)
      boot_results = future.apply::future_lapply(seq_len(n_boot), boot_one,
                                                 future.seed = T,
                                                 future.packages = c('causalDR')) %>%
        bind_rows()
    } else {
      boot_results = map(seq_len(n_boot), boot_one) %>% bind_rows()
    }

    if (boot_raw) {
      est_se = boot_results
    } else {
      est_se = est %>%
        pivot_longer(-c(method, time), names_to = 'estimand', values_to = 'est',
                     names_prefix = 'est_') %>%
        get_CI(get_se(boot_results), alpha = alpha, CI_raw = CI_raw)
    }
    if (!add_naive_est & !add_IPW & !specify_name) est_se = select(est_se, -method)
  }
  if (!add_naive_est & !add_IPW & !specify_name) est = select(est, -method)

  est = est %>% select(any_of(c('method', 'time')), contains('CIF')) %>%
    select(-contains('_trans'))
  if (!is.null(est_se) && !boot_raw) est_se = filter(est_se, str_detect(estimand, 'CIF'))

  list(est = est,
       est_se = est_se,
       dLambda = fit$dLambda,
       runtime = fit$runtime,
       call = match.call())
}
