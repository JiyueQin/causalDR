library(survival)

test_that('cum_step evaluates a right-continuous step function', {
  expect_equal(cum_step(c(1, 2, 4), c(0.1, 0.2, 0.3), c(0, 1, 1.5, 2, 5)),
               c(0, 0.1, 0.1, 0.3, 0.6))
})

test_that('breslow and fit_hazard match survival for (start, stop] data', {
  set.seed(1)
  n = 80
  start = runif(n)
  stop = start + rexp(n)
  event = rbinom(n, 1, 0.7)
  weights = runif(n, 0.5, 2)
  dat = data.frame(A = rbinom(n, 1, 0.5), Z1 = rnorm(n))

  # with lp = 0 the Breslow estimator is the weighted Nelson-Aalen estimator
  base_cumhaz = breslow(start, stop, event, rep(0, n), weights)
  km = survfit(Surv(start, stop, event) ~ 1, weights = weights)
  expect_equal(cum_step(base_cumhaz$time, base_cumhaz$jump, km$time), km$cumhaz)

  hz = fit_hazard(dat, start, stop, event, weights, vars = c('A', 'Z1'))
  cox = coxph(Surv(start, stop, event) ~ A + Z1, data = dat, weights = weights,
              ties = 'breslow', timefix = F)
  base = basehaz(cox, centered = F)
  expect_equal(unname(hz$beta), unname(coef(cox)))
  expect_equal(cum_step(hz$base_cumhaz$time, hz$base_cumhaz$jump, base$time), base$hazard)
})

test_that('eval_S_D gives the left limit of the Kaplan-Meier curve', {
  dat = data.frame(Q = 0, X2 = c(1, 2, 3, 4, 1, 2), delta2 = c(0, 1, 0, 1, 0, 0),
                   A = c(0, 0, 0, 0, 1, 1))
  S_D_fit = fit_S_D(dat, rep(1, 6))
  # arm 0: censorings at 1 and 3 among risk sets of size 4 and 2
  expect_equal(eval_S_D(S_D_fit, c(1, 1.5, 3, 3.5), rep(0, 4)), c(1, 3/4, 3/4, 3/8))
  expect_equal(eval_S_D(S_D_fit, c(1, 2, 2.5), rep(1, 3)), c(1, 1/2, 0))

  # no censoring in an arm: S_D = 1 without a warning
  dat_uncensored = mutate(dat, delta2 = replace(delta2, A == 1, 1))
  expect_no_warning(S_D_fit <- fit_S_D(dat_uncensored, rep(1, 6)))
  expect_equal(eval_S_D(S_D_fit, c(0.5, 3), rep(1, 2)), c(1, 1))

  # an empty arm: S_D = 1 with a warning
  expect_warning(S_D_fit <- fit_S_D(dat, as.numeric(dat$A == 0)), 'No subjects with A = 1')
  expect_equal(eval_S_D(S_D_fit, c(0.5, 3), rep(1, 2)), c(1, 1))
})

test_that('eval_G on a grid agrees with subject-specific evaluation', {
  set.seed(2)
  dat = sim_data_trunc(n = 150)$dat %>%
    mutate(Delta = delta1 + (1 - delta1)*delta2)
  grid = sort(unique(c(dat$Q, dat$X1, dat$X2)))
  S_D_fit = fit_S_D(dat, rep(1, nrow(dat)))
  G_fit = fit_G(dat, rep(1, nrow(dat)), S_D_fit, c('Z1', 'Z2'), 'Cox',
                tau0 = max(grid) + 1, S_min = 0.05)
  G_grid = eval_G(G_fit, dat, grid, on_grid = T)
  index_X1 = match(dat$X1, grid)
  expect_equal(G_grid[cbind(seq_len(nrow(dat)), index_X1)], eval_G(G_fit, dat, dat$X1))
  expect_true(all(G_grid >= 0 & G_grid <= 1))
  expect_true(all(apply(G_grid, 1, diff) >= -1e-12))

  # G is 0 below the smallest truncation time and 1 at the end of the grid, so its
  # jumps on the grid, starting from 0, carry the full mass
  q1 = min(dat$Q[dat$Delta == 1])
  expect_true(all(G_grid[, grid < q1] == 0))
  expect_true(all(G_grid[, grid >= q1] > 0))
  dG = G_grid - cbind(0, G_grid[, -length(grid)])
  expect_equal(rowSums(dG), rep(1, nrow(dat)))
})

test_that('eval_G puts the remaining mass at a point mass of Q at 0', {
  set.seed(3)
  dat = sim_data_trunc(n = 200, Q_model = 'reverse_PH')$dat %>%
    mutate(Delta = delta1 + (1 - delta1)*delta2)
  grid = sort(unique(c(dat$Q, dat$X1, dat$X2)))
  S_D_fit = fit_S_D(dat, rep(1, nrow(dat)))
  G_fit = fit_G(dat, rep(1, nrow(dat)), S_D_fit, c('Z1', 'Z2'), 'Cox',
                tau0 = max(grid) + 1, S_min = 0.05)
  G_grid = eval_G(G_fit, dat, grid, on_grid = T)
  # grid[1] = 0 is the point mass, and the reverse hazard dG(0)/G(0) is 1
  expect_equal(grid[1], 0)
  expect_true(all(G_grid[, 1] > 0))
  expect_equal(eval_G(G_fit, dat, rep(-1e-8, nrow(dat))), rep(0, nrow(dat)))
})

test_that('trunc_aug_cpp matches a brute-force evaluation of its definitions', {
  set.seed(3)
  n = 4; m = 7
  S = t(apply(matrix(runif(n*m, 0.9, 1), n), 1, cumprod))
  dLambda1 = matrix(rbinom(n*m, 1, 0.6)*runif(n*m, 0, 0.2), n)
  lp3 = matrix(rnorm(n*m, 0, 0.3), n)
  dB3 = matrix(runif(m, 0, 0.2), n, m, byrow = T)
  B3 = row_cumsum(dB3)
  G = t(apply(matrix(runif(n*m, 0.5, 1), n), 1, function(x) sort(x)))
  dG = G - cbind(0.3, G[, -m])
  index_Q = c(0L, 1L, 2L, 0L)
  index_X1 = c(3L, 5L, 6L, 2L)
  need_m = c(T, T, F, T)
  k_target = c(1L, 3L, 4L, 6L)
  S_min = 0.05; G_min = 0.01

  out = trunc_aug_cpp(S, dLambda1, lp3, B3, dB3, G, dG, index_Q, index_X1, need_m, k_target, S_min, G_min)

  # 0-based grid index j is column j + 1
  for (i in 1:n){
    b = function(u, k){
      S_minus = if (u > 0) S[i, u] else 1
      B3_minus = if (k > 0) B3[i, k] else 0
      S_minus*dLambda1[i, u + 1]*exp(-(B3_minus - B3[i, u + 1])*exp(lp3[i, u + 1]))
    }
    R = function(v, k) sum(map_dbl(0:(m - 1), ~if (.x <= v && .x < k) b(.x, k) else 0))
    rho = function(v, k) sum(map_dbl(0:(m - 1), ~if (.x <= v && .x < k)
      b(.x, k)*exp(lp3[i, .x + 1])*dB3[i, k + 1] else 0))
    g = function(v) dG[i, v + 1]/(max(S[i, v + 1], S_min)*max(G[i, v + 1], G_min)^2)
    SG_Q = max(S[i, index_Q[i] + 1], S_min)*max(G[i, index_Q[i] + 1], G_min)
    for (kk in seq_along(k_target)){
      k = k_target[kk]
      expect_equal(out$R_diagonal[i, kk], R(k, k))
      expect_equal(out$rho_diagonal[i, kk], rho(k, k))
      if (need_m[i] && index_Q[i] < k) {
        v_seq = index_Q[i]:(min(index_X1[i], k) - 1)
        m3 = R(index_Q[i], k)/SG_Q - sum(map_dbl(v_seq, ~R(.x, k)*g(.x)))
        m4 = rho(index_Q[i], k)/SG_Q - sum(map_dbl(v_seq, ~rho(.x, k)*g(.x)))
      } else {
        m3 = m4 = 0
      }
      expect_equal(out$m3[i, kk], m3)
      expect_equal(out$m4[i, kk], m4)
    }
  }
})

test_that('CIF_from_dLambda computes the CIFs from the hazard jumps', {
  dL = cbind(c(0.1, 0.2), c(0.05, 0.1), c(0.3, 0.4))
  out = CIF_from_dLambda(list(dL, dL), grid_k = c(1, 2), time_interest = c(1, 2),
                         time1_interest = 1, estimand = 'all', add_contrast = T)
  S1 = exp(-0.15)
  expect_equal(out$est_CIF1_0, c(0.1, 0.1 + S1*0.2))
  expect_equal(out$est_CIF2_1, c(0.05, 0.05 + S1*0.1))
  expect_equal(out$est_CIF3_0, c(0, 1 - exp(-0.4)))
  expect_equal(out$est_CIF1_diff, c(0, 0))
})

test_that('causalAIPWtrunc returns valid estimates', {
  set.seed(4)
  dat = sim_data_trunc(n = 200)$dat
  fit = causalAIPWtrunc(dat, covars = c('Z1', 'Z2'), time_interest = c(2, 3),
                        time1_interest = 1.5, add_IPW = T, add_naive_est = T,
                        name = 'myAIPW', n_boot = 0)
  expect_equal(unique(fit$est$method), c('myAIPW', 'IPW', 'naive'))
  expect_named(fit$est, c('method', 'time', paste0('est_CIF', rep(1:3, each = 2), '_', 0:1)))
  CIF = as.matrix(fit$est[, -(1:2)])
  expect_true(all(CIF >= 0 & CIF <= 1))
  expect_named(fit$dLambda, c('myAIPW', 'IPW', 'naive'))
})

test_that('causalAIPWtrunc rejects invalid input', {
  set.seed(5)
  dat = sim_data_trunc(n = 50)$dat
  run = function(d = dat, covars = c('Z1', 'Z2'), ...){
    causalAIPWtrunc(d, covars = covars, estimand = 'CIF1', n_boot = 0, ...)
  }
  expect_error(run(covars = 'W'), 'Columns not found')
  expect_error(run(weights = rep(1, 3)), 'weights')
  expect_error(run(PS_model = 'RF'), 'PS_model')
  expect_error(run(mutate(dat, Q = X1)), 'Q < X1')
  expect_error(run(mutate(dat, A = A + 1)), '0/1')
  expect_error(run(mutate(dat, Z1 = replace(Z1, 1, NA))), 'NA')
  expect_error(run(time_interest = -1), 'positive')
  expect_error(causalAIPWtrunc(dat, covars = c('Z1', 'Z2'), n_boot = 0), 'non-terminal')
})

test_that('a tie X1 = X2 with delta1 = delta2 = 1 is not a Lambda3 event', {
  set.seed(7)
  dat = sim_data_trunc(n = 150)$dat
  i = which(dat$delta1 == 1 & dat$delta2 == 1)[1]
  # subject i with both events at X1 (a tie), or censored at X1 (delta2 = 0)
  tie = mutate(dat, X2 = replace(X2, i, X1[i]))
  censored = mutate(tie, delta2 = replace(delta2, i, 0))
  run = function(d){
    causalAIPWtrunc(d, covars = c('Z1', 'Z2'), time1_interest = 1.5,
                    add_naive_est = TRUE, n_boot = 0)
  }
  fit_tie = expect_no_error(run(tie))
  fit_censored = run(censored)
  # neither has risk time after X1, so the naive Lambda3 jumps agree; the naive
  # estimator does not use S_D, which differs as delta2 also enters its Kaplan-Meier
  expect_equal(fit_tie$dLambda$naive[, c('dLambda3_0', 'dLambda3_1')],
               fit_censored$dLambda$naive[, c('dLambda3_0', 'dLambda3_1')])
})

test_that('dLambda_trunc rejects nuisance estimates on a different grid or data', {
  set.seed(6)
  dat = sim_data_trunc(n = 80)$dat %>%
    select(Q, X1, X2, delta1, delta2, A, Z1, Z2) %>%
    as.data.frame() %>%
    mutate(Delta = delta1 + (1 - delta1)*delta2,
           status = ifelse(delta1 == 1, 1, 2*delta2))
  grid = sort(unique(c(dat$Q, dat$X1, dat$X2)))
  nuisance = nuisance_trunc(dat, grid, rep(1, nrow(dat)), c('Z1', 'Z2'))
  run = function(d = dat, g = grid){
    dLambda_trunc(d, g, nuisance, a = 1, kmax = 5, weights = rep(1, nrow(d)),
                  S_min = 0.05, G_min = 0.01)
  }
  expect_no_error(run())
  expect_error(run(g = grid[-1]), 'grid must contain')
  expect_error(run(d = dat[-1, ]), 'do not match')
})

test_that('the misspecified models use the squared covariates', {
  dat = data.frame(A = c(0, 1, 1, 0), Z1 = c(-1, 2, 0.5, 3), T1 = c(1, 2, 3, 4))
  X = make_design(dat, c('A', 'Z1', 'T1'), 'Cox.quad')$X
  # Z1 is squared, A and T1 enter linearly; 'var' keeps the variable names
  expect_equal(unname(X[, 2]), dat$Z1^2)
  expect_equal(unname(X[, 1]), dat$A)
  expect_equal(unname(X[, 3]), dat$T1)
  expect_equal(attr(X, 'var'), c('A', 'Z1', 'T1'))

  set.seed(8)
  sim = sim_data_trunc(n = 200)$dat
  for (model in c('logit.quad', 'Cox.quad')) {
    args = list(covars = c('Z1', 'Z2'), time_interest = c(2, 3), time1_interest = 1.5,
                add_IPW = TRUE, n_boot = 0)
    if (model == 'logit.quad') args$PS_model = model
    else args[c('trunc_model', 'time_model')] = model
    fit = expect_no_error(do.call(causalAIPWtrunc, c(list(sim), args)))
    CIF = as.matrix(fit$est[, -(1:2)])
    expect_true(all(CIF >= 0 & CIF <= 1))
  }
  # other model names are still rejected
  expect_error(causalAIPWtrunc(sim, covars = c('Z1', 'Z2'), PS_model = 'logit.cube',
                               time1_interest = 1.5, n_boot = 0), 'PS_model')
})

test_that('breslow handles event times that differ only beyond 15 digits', {
  # as.character() maps both stop times to '1', which once gave duplicated factor
  # levels; such times arise on the reversed time scale tau0 - Q of fit_G()
  stop = c(1, 1 + 4e-16, 2, 3)
  expect_true(length(unique(stop)) == 4 && as.character(stop[1]) == as.character(stop[2]))
  bh = breslow(start = rep(0, 4), stop = stop, event = c(1, 1, 1, 0), lp = rep(0, 4),
               weights = rep(1, 4))
  # risk sets 4, 3, 2 at the three distinct event times
  expect_equal(bh$time, stop[1:3])
  expect_equal(bh$jump, c(1/4, 1/3, 1/2))
})
