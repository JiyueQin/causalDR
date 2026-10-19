test_that("full_data_est_trunc matches hand calculations", {
  # arm 0: two subjects with the non-terminal event, two without
  # arm 1: no non-terminal events
  full = data.frame(T1_0 = c(1, 2, Inf, Inf), T2_0 = c(3, 2.5, 1.5, 4),
                    T1_1 = Inf, T2_1 = c(1, 2, 3, 4))
  est = full_data_est_trunc(full, time_interest = c(1, 2, 3), time1_interest = 1.5,
                            add_contrast = TRUE)

  expect_equal(est$time, c(1, 2, 3))
  expect_equal(est$est_CIF1_0, c(0.25, 0.5, 0.5))
  expect_equal(est$est_CIF2_0, c(0, 0.25, 0.25))
  expect_equal(est$est_CIF1_1, c(0, 0, 0))
  expect_equal(est$est_CIF2_1, c(0.25, 0.5, 0.75))
  # Lambda3 jumps 1/2 at u = 2.5 (both at risk) and 1/1 at u = 3
  expect_equal(est$est_CIF3_0, c(0, 0, 1 - exp(-1.5)))
  expect_equal(est$est_CIF3_1, c(0, 0, 0))
  expect_equal(est$est_CIF2_diff, est$est_CIF2_1 - est$est_CIF2_0)
})

test_that("full_data_est_trunc has the same columns as causalAIPWtrunc estimates", {
  full = data.frame(T1_0 = c(1, 2, Inf, Inf), T2_0 = c(3, 2.5, 1.5, 4),
                    T1_1 = c(0.5, Inf, 2, Inf), T2_1 = c(1, 2, 3, 4))
  time_interest = c(1, 2, 3)
  time1_interest = c(0.5, 1.5)
  dL = matrix(0.1, 4, 3)
  for (contrast in c(FALSE, TRUE)){
    ref = CIF_from_dLambda(list(dL, dL), c(0.5, 1, 2, 3), time_interest, time1_interest,
                           'all', contrast)
    est = full_data_est_trunc(full, time_interest = time_interest,
                              time1_interest = time1_interest, add_contrast = contrast)
    expect_equal(names(est), names(ref))
  }
  est = full_data_est_trunc(full, estimand = 'CIF2', time_interest = time_interest,
                            name = 'full data')
  expect_equal(names(est), c('method', 'time', 'est_CIF2_0', 'est_CIF2_1'))
  expect_error(full_data_est_trunc(full, time_interest = 1), 'condition on')
  expect_error(full_data_est_trunc(full[, -1], estimand = 'CIF1', time_interest = 1),
               'not found')
})

test_that("sim_data_trunc returns the full data that generated the observed data", {
  set.seed(1)
  sim = sim_data_trunc(n = 200)
  dat = sim$dat
  full = sim$dat_full

  expect_named(sim, c('dat', 'dat_full', 'prop_truncated'))
  expect_equal(nrow(dat), 200)
  expect_true(sim$prop_truncated > 0 && sim$prop_truncated < 1)
  expect_equal(sum(full$observed), 200)
  expect_true(full$observed[nrow(full)])
  obs = full[full$observed, ]
  expect_equal(obs$id, dat$id)
  expect_equal(obs$Q, dat$Q)
  expect_equal(obs$A, dat$A)
  expect_equal(obs$Z1, dat$Z1)
  # factual outcomes of the observed subjects agree with the observed times
  T1 = ifelse(obs$A == 1, obs$T1_1, obs$T1_0)
  T2 = ifelse(obs$A == 1, obs$T2_1, obs$T2_0)
  expect_true(all(dat$X1 <= pmin(T1, T2) + 1e-12))
  expect_true(all(obs$Q < pmin(T1, T2)))
})

test_that("full data estimates are close to the true values of sim_data_trunc", {
  set.seed(2026)
  sim = sim_data_trunc(n = 500, Q_model = 'reverse_PH')
  est = full_data_est_trunc(sim$dat_full, time_interest = c(0.5, 3),
                            time1_interest = 0.5)
  # Lambda_0^a(t) = 0.2 t and 0.3 t; CIF1 = CIF2 = {1 - exp(-2 Lambda_0)}/2
  L0 = list(function(t) 0.2*t, function(t) 0.3*t)
  for (a in 0:1){
    truth12 = (1 - exp(-2*L0[[a + 1]](c(0.5, 3))))/2
    truth3 = c(0, 1 - exp(L0[[a + 1]](0.5) - L0[[a + 1]](3)))
    expect_equal(est[[paste0('est_CIF1_', a)]], truth12, tolerance = 0.1)
    expect_equal(est[[paste0('est_CIF2_', a)]], truth12, tolerance = 0.1)
    expect_equal(est[[paste0('est_CIF3_', a)]][2], truth3[2], tolerance = 0.1)
    expect_equal(est[[paste0('est_CIF3_', a)]][1], 0)
  }
})
