#' Simulate semi-competing risks data with treatment, confounding, and censoring
#'
#' Generates a simulated dataset with a binary treatment \eqn{A}, baseline covariates
#' \eqn{Z_1, Z_2}, semi-competing risks event times \eqn{T_1} (nonterminal event, e.g.,
#' illness) and \eqn{T_2} (terminal event, e.g., death), and right censoring \eqn{C}.
#' then the observed data are constructed as:
#' \deqn{X_1 = \min(T_1, T_2, C), \quad X_2 = \min(T_2, C),}
#' with event indicators \eqn{\delta_1 = I(T_1 = X_1)} and \eqn{\delta_2 = I(T_2 = X_2)}.
#'
#' Treatment is assigned by a logistic propensity score model
#' \eqn{\text{logit}\{P(A=1 \mid Z)\} = (1, Z_1, Z_2)^\top \alpha}. Censoring time follows an exponential distribution
#'  with rate \eqn{\exp\{(1, Z_1, Z_2, A)^\top \gamma\}}.
#'
#' The semi-competing risks structure allows \eqn{T_1=\infty} with positive probability.
#' Parameters \code{k2} and \code{k3} control multiplicative scaling of hazards relative
#' to a baseline hazard \eqn{\lambda_0} (specified separately for \eqn{A=0} and \eqn{A=1}).
#'
#' @param n Integer sample size. Default is \code{500}.
#'
#' @param k2 Numeric vector of length 2. Scaling ratios \eqn{(\,k_{2,0}, k_{2,1}\,)}
#' for \eqn{a=0,1}, representing the ratio \eqn{\lambda_{2a}/\lambda_{0a}} used in the
#' event-time generation.
#'
#' @param k3 Numeric vector of length 2. Scaling ratios \eqn{(\,k_{3,0}, k_{3,1}\,)}
#' for \eqn{a=0,1}, representing the ratio \eqn{\lambda_{3a}/\lambda_{0a}} used in the
#' event-time generation.
#'
#' @param alpha Numeric vector of coefficients \eqn{(\alpha_0,\alpha_1,\alpha_2)}
#' for the propensity score logistic model using design matrix \code{cbind(1, Z1, Z2)}.
#'
#' @param gamma Numeric vector of coefficients \eqn{(\gamma_0,\gamma_1,\gamma_2,\gamma_A)}
#' for the censoring model. Censoring times are generated as \code{rexp(n, rate = exp(eta))}
#' with \eqn{\eta=(1, Z_1, Z_2, A)^\top \gamma}.
#'
#' @param Lambda0_0 Function giving the baseline cumulative hazard function \eqn{\Lambda_{0,0}(t)}
#' for \eqn{A=0}. Must accept a numeric vector \code{t} and return a numeric vector.
#'
#' @param Lambda0_inv_0 Function giving the inverse \eqn{\Lambda_{0,0}^{-1}(u)} used for
#' inverse-transform sampling under \eqn{A=0}.
#'
#' @param Lambda0_1 Function giving the baseline cumulative hazard function \eqn{\Lambda_{0,1}(t)}
#' for \eqn{A=1}. If \code{NULL}, it is set equal to \code{Lambda0_0}.
#'
#' @param Lambda0_inv_1 Function giving the inverse \eqn{\Lambda_{0,1}^{-1}(u)} used for
#' inverse-transform sampling under \eqn{A=1}. If \code{NULL}, it is set equal to
#' \code{Lambda0_inv_0}.
#'
#' @details
#' The returned tibble contains the observed time variables \code{X1} and \code{X2} and
#' their corresponding indicators \code{delta1} and \code{delta2}. The column \code{event_type}
#' is a human-readable label for the observed event type:
#' \itemize{
#'   \item \code{"0.censored before illness or death"}: censored before \eqn{T_1} and \eqn{T_2}
#'   \item \code{"1.illness then censor"}: \eqn{T_1} observed, then censored before \eqn{T_2}
#'   \item \code{"2.death without illness"}: \eqn{T_2} observed before \eqn{T_1}
#'   \item \code{"3.illness then death"}: both \eqn{T_1} and \eqn{T_2} observed with \eqn{T_1 < T_2}
#' }
#'
#' @return A tibble with \code{n} rows and the following columns:
#' \describe{
#'   \item{\code{id}}{Subject identifier \code{1:n}.}
#'   \item{\code{A}}{Treatment indicator (0/1).}
#'   \item{\code{X1}}{Observed time \eqn{X_1=\min(T_1,T_2,C)}.}
#'   \item{\code{X2}}{Observed time \eqn{X_2=\min(T_2,C)}.}
#'   \item{\code{delta1}}{Indicator \eqn{I(T_1=X_1)} for observing the nonterminal event.}
#'   \item{\code{delta2}}{Indicator \eqn{I(T_2=X_2)} for observing the terminal event.}
#'   \item{\code{Z1, Z2}}{Baseline covariates.}
#'   \item{\code{event_type}}{Event type label derived from \code{delta1} and \code{delta2}.}
#' }
#'
#' @examples
#' # Simulate data with defaults
#' dat <- sim_data()
#' head(dat)
#'
#' # Change hazard scaling parameters
#' dat2 <- sim_data(n = 200, k2 = c(1, 2), k3 = c(1, 1.5))
#'
#' # Use same baseline hazard under A=1 as under A=0, i.e., no treatment effect
#' dat3 <- sim_data(n = 200, Lambda0_1 = NULL, Lambda0_inv_1 = NULL)
#' @export
sim_data = function(n = 500,
                    # k2 = (k2_0, k2_1), the ratio of lambda2_a over lambda0_a
                    k2 = c(1,1),
                    # k3 = (k3_0, k3_1), the ratio of lambda3_a over lambda0_a
                    k3 = c(1,1),
                    # coefficients of PS logit model:(intercept, Z1,Z2)
                    alpha = c(0, 1, 1),
                    # coefficients of censoring model(hazard of C given Z and A): (intercept, Z1, Z2, Z3, A)
                    gamma = c(-2.5, 1, 1, 0.5),
                    Lambda0_0 = function(t){0.2*t},
                    Lambda0_inv_0 = function(a){5*a},
                    Lambda0_1 = function(t){0.1*t^2} ,
                    Lambda0_inv_1 = function(a){sqrt(10*a)})
{
  U1 <- runif(n, 0, 1)
  U2 <- runif(n, 0, 1)

  if (is.null(Lambda0_1)){
    Lambda0_1 = Lambda0_0
    Lambda0_inv_1 = Lambda0_inv_0
  }



  Z1 = U1 -0.5
  Z2 = U2 - 0.5

  Z_matrix = cbind(1, Z1, Z2)
  pA = as.numeric(inv.logit(Z_matrix %*% alpha))
  A = rbinom(n, 1, prob = pA)

  ZA_matrix = cbind(Z_matrix, A)
  Z_matrix_1 = cbind(Z_matrix, 1)
  Z_matrix_0 = cbind(Z_matrix, 0)
  # get the list of quantities for denominators

  d_ls = map(list(1,k2[1], k3[1], 1, k2[2], k3[2]), ~rep(.x, n))
  names(d_ls) = paste0('d', rep(1:3, 2), '_', rep(0:1, each = 3))


  p_T1_infty_0 = d_ls$d2_0 / (d_ls$d1_0 + d_ls$d2_0)
  T1_inf_ind_0 =  rbinom(n, 1, p_T1_infty_0)
  T1_0 <- NULL
  T2_0 <- NULL
  for (i in 1:n) {
    if (T1_inf_ind_0[i] == 1) {
      T1_0[i] = Inf
      T2_0[i] = Lambda0_inv_0(-log(U2[i]) / (d_ls$d1_0[i] + d_ls$d2_0[i]))
    }
    if (T1_inf_ind_0[i] == 0) {
      T1_0[i] = Lambda0_inv_0(-log(U1[i]) / (d_ls$d1_0[i] + d_ls$d2_0[i]))
      T2_0[i] = Lambda0_inv_0(-(log(U2[i]) / d_ls$d3_0[i]) + Lambda0_0(T1_0[i]))
    }
  }



  p_T1_infty_1 = d_ls$d2_1 / (d_ls$d1_1 + d_ls$d2_1)
  T1_inf_ind_1 =  rbinom(n, 1, p_T1_infty_1)
  T1_1 <- NULL
  T2_1 <- NULL
  for (i in 1:n) {
    if (T1_inf_ind_1[i] == 1) {
      T1_1[i] = Inf
      T2_1[i] = Lambda0_inv_1(-log(U2[i]) / (d_ls$d1_1[i] + d_ls$d2_1[i]))
    }
    if (T1_inf_ind_1[i] == 0) {
      T1_1[i] <- Lambda0_inv_1(-log(U1[i]) / (d_ls$d1_1[i] + d_ls$d2_1[i]))
      T2_1[i] <- Lambda0_inv_1(-(log(U2[i]) / d_ls$d3_1[i]) + Lambda0_1(T1_1[i]))
    }
  }

  C_0 = rexp(n, rate = exp(Z_matrix_0 %*% gamma))
  C_1 = rexp(n, rate = exp(Z_matrix_1 %*% gamma))



  dat = tibble(
    id = 1:n,
    A,
    Z1,
    Z2,
    T1_0,
    T2_0,
    T1_inf_ind_0,
    T1_1,
    T2_1,
    T1_inf_ind_1,
    C_0,
    C_1,
    pA,
    U1,
    U2
  ) %>%
    mutate(
      T1 = ifelse(A, T1_1, T1_0),
      T2 = ifelse(A, T2_1, T2_0),
      C = ifelse(A, C_1, C_0),
      T1_inf_ind = ifelse(A, T1_inf_ind_1, T1_inf_ind_0)
    ) %>%
    mutate(
      X1 = pmin(T1, T2, C),
      X2 = pmin(T2, C),
      delta1 = 1 * (T1 == X1),
      delta2 = 1 * (T2 == X2)) %>%
    mutate(event_type = case_when(
      (1-delta1) & (1-delta2) ~ '0.censored before illness or death',
      delta1 & (1-delta2) ~ '1.illness then censor',
      (1-delta1) & delta2 ~ '2.death without illness',
      delta1 & delta2 ~ '3.illness then death'
    )) %>%
    select(id, A, X1, X2, delta1, delta2, starts_with('Z'), event_type)

  dat

}
