.plot_estimates = function(fit, estimand = 'all_CIF', CI = T,
                           xlab='Time',
                           ylab = 'Value',
                           time_cut=NULL, title=NULL, truth=NULL,
                           point_size = 1, pointrange_size=0.1){

  if (estimand == 'all_Lambda') {
    estimand_vec = paste0('Lambda', c(1,2,3,1,2,3), '_', c(0,0,0,1,1,1))}
  if (estimand == 'all_CIF') {
    estimand_vec = paste0('CIF', rep(1:3, each = 4), '_', rep(c(0,1, 'diff', 'ratio'),2))}

  if(CI){
    dat_plot = fit$est_se
  }else{
    dat_plot = fit$est
  }

  if (!'method' %in% colnames(dat_plot)) {dat_plot = dat_plot %>%
    mutate(method = 'causalAIPCW')}

  if(!CI){
    dat_plot = dat_plot %>%
      pivot_longer(-c(time, method),
                   names_to = 'estimand',
                   values_to = 'est',
                   names_prefix = 'est_')
  }


  if(is.null(ylab)) ylab = param

  if (!is.null(truth)){
    if(is.vector(truth)){
      dat_plot = dat_plot %>%
        mutate(truth = truth)
    }else if (is.data.frame(truth)){
      # use bind_rows instead of rbind since truth does not have variables about CI
      dat_plot = bind_rows(dat_plot,truth)
    }
    else{
      dat_plot = dat_plot %>%
        mutate(truth = map_dbl(time, truth))

    }
  }
  dat_plot = dat_plot %>%
    filter(estimand %in% estimand_vec)

  if (!is.null(time_cut)) {
    dat_plot = dat_plot %>%
      filter(time<= time_cut)
  }

  p = dat_plot %>%
    ggplot(aes(x = time, y = est, color = method)) +
    geom_point(size =point_size) +
    geom_line() +
    labs(x = xlab, y = ylab, title = title) +
    theme_bw() +
    theme(legend.position = 'bottom') +
    ggsci::scale_color_d3(name = 'Method') +
    facet_wrap(~estimand, scales = 'free')
  if(CI) {
    p = p+geom_pointrange(aes(ymin = lower, ymax = upper),
                          size = pointrange_size)
  }
  p
}
#' Plot estimated causal estimands over time
#'
#' Produces publication-ready plots of time-varying causal estimands (e.g., cumulative
#' incidence functions or cumulative hazards) from a fitted object returned by
#' \code{causalAIPCW()}.
#'
#' The function supports plotting point estimates alone or point estimates with
#' pointwise confidence intervals, and automatically facets plots by estimand.
#'
#' @param fit A fitted object returned by \code{causalAIPCW()} (or a compatible object)
#'   containing \code{fit$est} (point estimates) and optionally \code{fit$est_se}
#'   (estimates with confidence intervals).
#'
#' @param estimand Character string specifying which family of estimands to plot.
#'   Supported values include:
#'   \itemize{
#'     \item \code{"all_CIF"} (default): cumulative incidence functions for all event types
#'       (and available contrasts).
#'     \item \code{"all_Lambda"}: cumulative hazard functions for all transitions.
#'   }
#'
#' @param CI Logical; if \code{TRUE} (default), plots pointwise confidence intervals using
#'   \code{fit$est_se}. If \code{FALSE}, plots point estimates using \code{fit$est} only.
#'
#' @param xlab Character x-axis label. Default is \code{"Time"}.
#'
#' @param ylab Character y-axis label. Default is \code{"Value"}.
#'
#' @param time_cut Optional numeric value. If provided, restricts the plot to time points
#'   \code{time <= time_cut}.
#'
#' @param title Optional character string giving the plot title.
#'
#' @param point_size Numeric point size passed to \code{ggplot2::geom_point()}.
#'
#' @param pointrange_size Numeric line width for confidence interval bars produced by
#'   \code{ggplot2::geom_pointrange()} when \code{CI = TRUE}.
#'
#'
#' @return A \pkg{ggplot2} object.
#'
#' @seealso \code{\link{causalAIPCW}}, \code{\link{sim_data}}
#'
#' @examples
#' \dontrun{
#' ## Simulate data and fit the model
#' dat <- sim_data(n = 500)
#' fit <- causalAIPCW(
#'   data = dat,
#'   covars = c("Z1", "Z2"),
#'   estimand = "all",
#'   time1_interest = c(1, 2),
#'   n_boot = 100,
#'   seed = 1
#' )
#'
#' ## Plot CIFs with confidence intervals
#' plot_estimates(fit, estimand = "all_CIF", CI = TRUE)}
#' @export

plot_estimates = function(fit, estimand = 'all_CIF', CI = T,
                          xlab='Time',
                          ylab = 'Value',
                          time_cut=NULL, title=NULL,
                          point_size = 1, pointrange_size=0.1){
  .plot_estimates(fit=fit, estimand = estimand, CI = CI,
                  xlab=xlab,
                  ylab = ylab,
                  time_cut=time_cut, title=title,
                  point_size = point_size, pointrange_size=pointrange_size)
}

plot_estimates_ls = function(fit, params = 'all_Lambda', xlab='Time',ylabs = NULL,
                             time_cut=NULL, titles=NULL, truth_ls=NULL,
                             ncol=3, nrow = 2, point_size = 1){
  plot_ls = list()
  if (params == 'all_Lambda') {params = paste0('Lambda', c(1,2,3,1,2,3), '_', c(0,0,0,1,1,1))}
  else if (params == 'all_CIF') {params = paste0('CIF', c(1,2,1,2), '_', c(0,0,1,1))}
  n_plot = length(params)
  if (length(titles)==1 & n_plot>1) {titles = rep(titles, n_plot)}

  for (i in 1:n_plot){
    plot_ls[[i]] = plot_estimates(fit, params[i], xlab, ylabs[i], time_cut, titles[i], truth_ls[[i]], point_size)
  }
  ggpubr::ggarrange(plotlist = plot_ls, ncol= ncol, nrow = nrow)

}
