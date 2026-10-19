test_that("Q_model = 'uniform' generates Q from the reverse-time PH model of Wang et al.", {
  eta = c(-0.6, 0.4, 0.2, 0.6)
  set.seed(3)
  sim = sim_data_trunc(n = 500, Q_model = 'uniform', eta = eta, q_max = 5, T_min = 1)
  full = sim$dat_full

  # continuous Q on (0, q_max): no ties at 0
  expect_true(all(full$Q > 0 & full$Q < 5))
  expect_false(any(duplicated(full$Q)))
  # G(Q | A, Z) = (Q/q_max)^exp(lp) is uniform in the population before truncation
  lp = eta[1] + eta[2]*full$Z1 + eta[3]*full$Z2 + eta[4]*full$A
  expect_gt(stats::ks.test((full$Q/5)^exp(lp), 'punif')$p.value, 0.001)
  # observed subjects satisfy Q < T
  expect_true(all(sim$dat$Q < sim$dat$X1))
})

test_that("T_min shifts the event times", {
  set.seed(4)
  sim = sim_data_trunc(n = 300, Q_model = 'uniform', eta = c(-0.6, 0.4, 0.2, 0.6),
                       q_max = 5, T_min = 1)
  full = sim$dat_full
  for (v in c('T1_0', 'T2_0', 'T1_1', 'T2_1')) expect_true(all(full[[v]] >= 1))

  # truth: CIF1 = CIF2 = {1 - exp(-2 Lambda_0(t - T_min))}/2 for t >= T_min
  set.seed(5)
  sim = sim_data_trunc(n = 500, Q_model = 'uniform', eta = c(-0.6, 0.4, 0.2, 0.6),
                       q_max = 5, T_min = 1)
  est = full_data_est_trunc(sim$dat_full, time_interest = c(0.5, 3), time1_interest = 1.5)
  L0 = list(function(t) 0.2*t, function(t) 0.3*t)
  for (a in 0:1){
    expect_equal(est[[paste0('est_CIF1_', a)]][1], 0)
    truth12 = (1 - exp(-2*L0[[a + 1]](3 - 1)))/2
    truth3 = 1 - exp(L0[[a + 1]](1.5 - 1) - L0[[a + 1]](3 - 1))
    expect_equal(est[[paste0('est_CIF1_', a)]][2], truth12, tolerance = 0.1)
    expect_equal(est[[paste0('est_CIF2_', a)]][2], truth12, tolerance = 0.1)
    expect_equal(est[[paste0('est_CIF3_', a)]][2], truth3, tolerance = 0.1)
  }
})

test_that("the default design is the uniform design of Wang et al.", {
  set.seed(6)
  sim_default = sim_data_trunc(n = 100)
  set.seed(6)
  sim_wang = sim_data_trunc(n = 100, Q_model = 'uniform', eta = c(-0.6, 0.4, 0.2, 0.6),
                            q_max = 5, T_min = 1)
  expect_equal(sim_default, sim_wang)

  # the defaults of reverse_PH are those of earlier versions: point mass at 0, no shift
  set.seed(6)
  sim_rph = sim_data_trunc(n = 100, Q_model = 'reverse_PH')
  set.seed(6)
  sim_rph2 = sim_data_trunc(n = 100, Q_model = 'reverse_PH', eta = c(0.3, -0.3, 0.3),
                            p0 = 0.2, q_max = 2, T_min = 0)
  expect_equal(sim_rph, sim_rph2)
  expect_true(any(sim_rph$dat$Q == 0))
})

test_that("sim_data_trunc checks its arguments", {
  expect_error(sim_data_trunc(n = 10, Q_model = 'weibull'), 'Q_model')
  expect_error(sim_data_trunc(n = 10, eta = c(1, 2, 3)), 'length 4')
  expect_error(sim_data_trunc(n = 10, Q_model = 'reverse_PH', eta = c(1, 2, 3, 4)), 'length 3')
  expect_error(sim_data_trunc(n = 10, Q_model = 'exp'), 'Q_model')
  expect_error(sim_data_trunc(n = 10, T_min = -1), 'T_min')
  expect_error(sim_data_trunc(n = 10, Q_model = 'reverse_PH', p0 = 1), 'p0')
})
