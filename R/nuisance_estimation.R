################################################################################
# Function predict_S_c - Predict the survival function of C given A and Z
# return: a list of two matrices of dimension:
#         number of subjects for prediction * number of time points of interest
#
################################################################################

predict_S_c = function(dat_predict, dat_train,
                       delta_c, X_c, time_interest,
                       model='Cox',
                       weights=NULL, nodesize = 15, ntree= 500) {
  dat_train = dat_train %>%
    select(all_of(c(X_c, delta_c)), A, contains('Z')) %>%
    rename(all_of(c(X_c = X_c, delta_c = delta_c)))

  dat_predict0 = cbind(A=0, select(dat_predict, contains('Z')))
  dat_predict1 = cbind(A=1, select(dat_predict, contains('Z')))

  #X_c_sort = sort(unique(dat_train %>% pull(X_c)))


  if (model %in% c('Cox', 'Cox.quad', 'Cox.small')) {

    if (model == 'Cox.quad') {
      dat_train = dat_train %>%
        mutate_at(vars(contains('Z')), ~.x^2)

      dat_predict0 = dat_predict0 %>% mutate_at(vars(contains('Z')), ~.x^2)
      dat_predict1 = dat_predict1 %>% mutate_at(vars(contains('Z')), ~.x^2)}
    if (model == 'Cox.small') {
      dat_train = dat_train %>% select(X_c, delta_c, A, Z1)
      dat_predict0 = dat_predict0 %>% select(A, Z1)
      dat_predict1 = dat_predict1 %>% select(A, Z1)
    }

    fit = coxph(
      Surv(X_c, delta_c) ~ .,
      timefix = F,
      data = dat_train,
      weights = weights
    )
    # cumulative hazard of C0:
    #  test data size * number of unique time points in train data
    basehaz_df = basehaz(fit, centered = F)
    Lambda_c0 = exp(as.matrix(dat_predict0) %*% fit$coefficients) %*%
      t(basehaz_df$hazard)
    S_c_0 = exp(-Lambda_c0)
    # cumulative hazard of C1:
    # test data size * number of unique time points in train data
    Lambda_c1 = exp(as.matrix(dat_predict1) %*% fit$coefficients) %*%
      t(basehaz_df$hazard)
    S_c_1 = exp(-Lambda_c1)
    # unique time points, including censoring
    # in_time = X_c_sort
    in_time = basehaz_df$time
  }

  if (model == 'Spline') {
    fit = hare(dat_train$X_c, dat_train$delta_c,
               as.matrix(dat_train %>% select(A, contains('Z'))))
    # same: do.call(rbind, map(time_interest, ~phare(q = .x,cov = dat_predict0,fit = fit))) %>% t()
    S_c_0 = 1 - sapply(time_interest, function(x)
      phare(
        q = x,
        cov = dat_predict0,
        fit = fit
      ))
    #S_c_0[is.nan(S_c_0)] = NA
    S_c_1 = 1 - sapply(time_interest, function(x)
      phare(
        q = x,
        cov = dat_predict1,
        fit = fit
      ))
    #S_c_1[is.nan(S_c_1)] = NA
    in_time = time_interest
  }
  if (model == 'RSF') {
    fit = rfsrc(
      Surv(X_c, delta_c) ~ .,
      data = dat_train,
      ntime = 0,
      #splitrule = 'logrank',
      splitrule = 'bs.gradient',
      ntree = ntree,
      nodesize = nodesize,
      mtry = 2,
      case.wt = weights
    )

    S_c_0 = predict(fit, newdata = dat_predict0)$survival
    S_c_1 = predict(fit, newdata = dat_predict1)$survival
    # unique event time points
    in_time = fit$time.interest
  }

  list(S_c_0 = S_c_0, S_c_1 = S_c_1) %>%
    map(~fill_time_matrix(.x, in_time, time_interest))

}



####################################################################################
# Function predict_Lambda - Predict the three cumulative transition rates
# given A and Z
# return: a list of 6(when model is cox) or 4(when model is non-parametric) matrices of dimension:
#         number of subjects for prediction* number of time points of interest
#
####################################################################################


predict_Lambda = function(dat_predict,
                          dat_train,
                          time,
                          status,
                          time_interest,
                          model = 'Cox',
                          weights = NULL,
                          nodesize = 15,
                          ntree = 500) {
  dat_train1 = dat_train %>%
    select(all_of(c(time, status)), A, contains('Z')) %>%
    rename(all_of(c(time = time, status = status)))

  with_illness = (dat_train$delta1 == 1)
  dat_train3 = dat_train %>%
    mutate(event3 = delta1*delta2) %>%
    filter(delta1 ==1) %>%
    select(X1, X2, event3, A, contains('Z'))

  dat_predict0 = cbind(A = 0, select(dat_predict, contains('Z')))
  dat_predict1 = cbind(A = 1, select(dat_predict, contains('Z')))


  if (model == 'Cox') {
    # Lambda1, Lambda2
    fit_ls = map(c(1, 2),
                 ~ coxph(
                   Surv(time, status == .x) ~ .,
                   timefix = F,
                   data = dat_train1,
                   weights = weights
                 ))
    # Lambda3 (assume Markov)
    fit_ls[[3]] = coxph(
      Surv(time=X1, time2=X2, event= event3) ~ .,
      timefix = F,
      data = dat_train3,
      weights = weights[with_illness]
    )

    #basehaz_df_ls = map(fit_ls, ~basehaz(.x, centered = F))
    Lambda_ls = c(map(
      fit_ls,
      ~ exp(as.matrix(dat_predict0) %*% .x$coefficients) %*%
        t(basehaz(.x, centered = F)$hazard)
    ), map(
      fit_ls,
      ~ exp(as.matrix(dat_predict1) %*% .x$coefficients) %*%
        t(basehaz(.x, centered = F)$hazard)
    ))
    # unique time points, including censored time points
    in_time1 = basehaz(fit_ls[[1]], centered = F)$time
    in_time3 = basehaz(fit_ls[[3]], centered = F)$time

  }
  if (model == 'RSF') {
    # a non-parametric way to estimate lambda 3 under Markov:
    # truncation forest LTRCforests: still truncating X1, off CRAN and no weights option
    # or LTRCtrees, which has weights option
    # but here, with a non-parametric approach, we do not want to assume Markov

    fit = rfsrc(
      Surv(time, status) ~ .,
      data = dat_train1,
      ntime = 0,
      splitrule = 'logrankCR',
      ntree = ntree,
      nodesize = nodesize,
      mtry = 2,
      case.wt = weights
    )

    Lambda0 = predict(fit, newdata = dat_predict0)$chf
    Lambda1 = predict(fit, newdata = dat_predict1)$chf
    Lambda_ls = list(Lambda0[, , 1], Lambda0[, , 2],
                     Lambda1[, , 1], Lambda1[, , 2])
    # unique time points of event 1 or event 2
    in_time1 = fit$time.interest
  }
  if(model == 'Cox'){
    # Lambda1_0: cumulative CS-hazard of illness for those without treatment
    # Lambda2_0: cumulative CS-hazard of death without illness without treatment
    names(Lambda_ls) = paste0('Lambda', c('1_0', '2_0', '3_0', '1_1', '2_1', '3_1'))
    map2(Lambda_ls,rep(list(in_time1, in_time1, in_time3),2),
         ~ fill_time_matrix(.x, .y, time_interest, 'Lambda'))

  }else{
    # Lambda1_0: cumulative CS-hazard of illness for those without treatment
    # Lambda2_0: cumulative CS-hazard of death without illness without treatment
    names(Lambda_ls) = paste0('Lambda', c('1_0', '2_0', '1_1', '2_1'))
    map(Lambda_ls,
        ~ fill_time_matrix(.x, in_time1, time_interest, 'Lambda'))
  }

}


################################################################################
# Function predict_S_T2_nonpara - Predict the survival function of T2 given A and Z non-parametrically
#                     if assume Markov, it can be estimated from three lambda
# return: a list of two matrices of dimension:
#         number of subjects for prediction * number of time points of interest
#
################################################################################

predict_S_T2_nonpara = function(dat_predict, dat_train,
                                X2, delta2,
                                time_interest,
                                model='RSF',
                                weights=NULL, nodesize = 15, ntree= 500) {
  dat_train = dat_train %>%
    select(all_of(c(X2, delta2)), A, contains('Z'))


  dat_predict0 = cbind(A=0, select(dat_predict, contains('Z')))
  dat_predict1 = cbind(A=1, select(dat_predict, contains('Z')))

  #X_c_sort = sort(unique(dat_train %>% pull(X_c)))

  if (model == 'RSF') {
    fit = rfsrc(
      Surv(X2, delta2) ~ .,
      data = dat_train,
      ntime = 0,
      # ?why not the default split rule of logrank
      # splitrule = 'logrank',
      splitrule = 'bs.gradient',
      ntree = ntree,
      nodesize = nodesize,
      mtry = 2,
      case.wt = weights
    )

    S_T2_0 = predict(fit, newdata = dat_predict0)$survival
    S_T2_1 = predict(fit, newdata = dat_predict1)$survival
    # unique event time points
    in_time = fit$time.interest
  }

  list(S_T2_0 = S_T2_0, S_T2_1 = S_T2_1) %>%
    map(~fill_time_matrix(.x, in_time, time_interest))

}



####################################################################################
# Function predict_S_t - predict the survival function of T = min(T1, T2) from a Lambda list
# return: a list of two matrices of dimension:
#         number of subjects * number of time points of interest
#
####################################################################################

#todo: can also estimate S_t non-parametrically directly instead of through Lambda.

predict_S_t = function(Lambda_ls) {
  out_ls = map(list(Lambda_ls$Lambda1_0 + Lambda_ls$Lambda2_0,
                    Lambda_ls$Lambda1_1 + Lambda_ls$Lambda2_1),
               ~ exp(-.x))
  names(out_ls) = paste0('S_t_', 0:1)
  out_ls
}



####################################################################################
# Function predict_S_T12_S_T2_markov - predict the S_T2 and S_T12 (P(T1<r<=T2))from S_t list and a Lambda list (assuming markov)
# return: a list of four matrices of dimension:
#         number of subjects * number of time points of interest
#
####################################################################################


predict_S_T12_S_T2_markov = function(S_t, Lambda) {
  dLambda1_0 = get_jump_size(Lambda$Lambda1_0)
  dLambda1_1 = get_jump_size(Lambda$Lambda1_1)
  #m = ncol(dLambda1_0)
  # exp(x) will be infinity when x is large such as 1000
  # could also use the trick to prevent underflow and overflow, here simply set to 0
  product0 = S_t[[1]]*exp(Lambda$Lambda3_0)*dLambda1_0
  product0 = ifelse(is.nan(product0)|!is.finite(product0), 0,product0)
  S_T12_0 = exp(-Lambda$Lambda3_0)*t(apply(product0, 1, cumsum))
  S_T2_0 = S_t[[1]] + S_T12_0
  # all(S_T2_0>=0&S_T2_0<=1)
  product1 = S_t[[2]]*exp(Lambda$Lambda3_1)*dLambda1_1
  product1 = ifelse(is.nan(product1)|!is.finite(product1), 0,product1)
  # S_T12_1 may be exactly 0
  S_T12_1 = exp(-Lambda$Lambda3_1)*t(apply(product1, 1, cumsum))
  S_T2_1 = S_t[[2]] + S_T12_1
  # all(S_T2_1>=0&S_T2_1<=1)

  out_ls = list(S_T12_0 = S_T12_0, S_T2_0 = S_T2_0, S_T12_1 = S_T12_1, S_T2_1 = S_T2_1)
  out_ls
}



####################################################################################
# Function predict_CIF1 - predict CIF1 from S_t list and Lambda1 list
# return: a list of two matrices of dimension:
#         number of subjects * number of time points of interest
#
####################################################################################


predict_CIF1 = function(S_t, Lambda) {
  dLambda1_0 = get_jump_size(Lambda[[1]])
  dLambda1_1 = get_jump_size(Lambda[[2]])
  product0 = S_t[[1]]*dLambda1_0
  CIF1_0 = t(apply(product0, 1, cumsum))
  product1 = S_t[[2]]*dLambda1_1
  CIF1_1 = t(apply(product1, 1, cumsum))

  out_ls = list(CIF1_0 = CIF1_0, CIF1_1 = CIF1_1)
  out_ls
}

####################################################################################
# Function predict_S_T12 - predict the S_T12 (P(T1<r<T2)) from S_t and S_t2
# return: a list of two matrices of dimension:
#         number of subjects * number of time points of interest
#
####################################################################################

predict_S_T12 = function(S_t, S_T2, epsilon=0.01) {
  # the estimates may be negative
  S_T12_0 = S_T2$S_T2_0-S_t$S_t_0
  S_T12_0[S_T12_0<epsilon] = epsilon
  S_T12_1 = S_T2$S_T2_1-S_t$S_t_1
  S_T12_1[S_T12_1<epsilon] = epsilon
  list(S_T12_0 =S_T12_0,
       S_T12_1 =S_T12_1)
}


####################################################################################
#  Function predict_PS - Predict propensity score
#  return: a vector of length: size of data for prediction
#
####################################################################################

predict_PS = function(dat_predict,
                      dat_train,
                      model = 'logit',
                      weights = NULL) {
  dat_train = dat_train %>% select(A, contains('Z'))
  dat_predict = dat_predict %>% select(A, contains('Z'))
  if (model %in% c('logit', 'logit.quad', 'logit.small')) {
    if (model == 'logit.quad'){
      dat_train = dat_train %>% select(A, contains('Z')) %>%
        mutate_at(vars(-A), ~.x^2)
      dat_predict = dat_predict %>% select(A, contains('Z')) %>%
        mutate_at(vars(-A), ~.x^2)
    }

    if (model == 'logit.small'){
      dat_train = dat_train %>% select(A, Z1)
      dat_predict = dat_predict %>% select(A, Z1)
    }
    fit = glm(A ~ .,
              data = dat_train,
              # same prediction when using binomial or quasibinomial
              family = binomial,
              weights = weights)
    PS = predict(fit, newdata = dat_predict, type = 'response')
  }

  if (model == 'RF') {
    fit = ranger::ranger(
      A ~ .,
      data = dat_train %>% mutate(A = factor(A)),
      probability = TRUE,
      num.trees = 2000,
      case.weights = weights,
    )
    PS = predict(fit, dat_predict)$predictions[, 2]
  }
  if (model == 'gbm') {
    # todo: the consideration of over-fitting.
    # Using the last tree imposes a high risk of over-fitting and bad test performance
    fit = gbm::gbm(
      A ~ .,
      data = dat_train,
      distribution = 'bernoulli',
      n.trees = 10000,
      interaction.depth = 1,
      weights = weights,
      cv.folds = 5,
      shrinkage = 0.01,
      # using only one core to avoid issues when running on OSG server
      n.cores = 1
    )
    best_tree = gbm::gbm.perf(fit, method = 'cv', plot.it = F)
    PS = predict(fit, dat_predict, n.trees = best_tree, type = 'response')
  }
  PS
}


nuisance_est = function(data,
                        n,
                        time_all,
                        PS_model = 'logit', censor_model = 'Cox',
                        time_model = 'Cox', c_nodesize = 15,
                        c_ntree = 500,  t_nodesize = 15, t_ntree = 500,
                        k = 5, crossfit = 'auto', weights = NULL,
                        shuffle = T,
                        PS_min = 0.1, S_min = 0.05
){
  start_time = Sys.time()

  m_all = length(time_all)
  if (is.null(weights)) weights = rep(1, n)

  # initialize nuisance estimates with NA
  # nuisance are estimated at all available time points so later they can be re-used without re-estimating

  if(time_model == 'Cox'){
    # estimate Lambda3 with cox under Markov assumption
    nuisuance_names_core = c('S_c', 'Lambda1', 'Lambda2', 'Lambda3', 'S_t', 'S_T12', 'S_T2')
  }else{
    nuisuance_names_core = c('S_c', 'Lambda1', 'Lambda2', 'S_T2', 'S_t', 'S_T12')
  }

  nuisuance_names = paste0(rep(nuisuance_names_core, 2), '_',
                           rep(0:1, each = length(nuisuance_names_core)))

  nuisance_ls = map(1:length(nuisuance_names), ~matrix(NA, nrow = n, ncol = m_all))

  names(nuisance_ls) = nuisuance_names

  # initialize a numeric vector of NAs to store propensity scores
  PS = rep(NA_real_, n)

  # get estimates of specified nuisance models by cross-fitting
  models = c(PS_model, censor_model, time_model)
  # can specify what models to crossfit
  crossfit_models = crossfit
  if (crossfit == 'all') {crossfit_models = c('PS','censor','time')}
  if (crossfit == 'none') {crossfit_models = NULL}
  if (crossfit == 'non-param') {
    #non_param_indic = (models != c('logit', 'Cox', 'Cox'))
    non_param_indic = !str_detect(models, c('logit', 'Cox', 'Cox'))
    crossfit_models = c('PS','censor','time')[non_param_indic]
  }
  if (crossfit == 'auto') {
    non_param_indic = !str_detect(models, c('logit', 'Cox', 'Cox'))
    if (any(non_param_indic)){crossfit_models = c('PS','censor','time')}
    else {crossfit_models = NULL}
    # Rotnitzky, A., E. Smucler, and J. M. Robins (2021). Characterization of parameters with a mixed bias property. Biometrika 108(1), 231–238.
  }

  if (length(crossfit_models)>0){
    folds = cut(1:n, breaks = k, labels = F)
    # recommend to shuffle the data to create folds to prevent any inherent pattern
    # set.seed(1)
    if (shuffle) folds = sample(folds,length(folds))
    for (fold in 1:k) {
      infold = (folds == fold)
      dat_infold = data[infold, ]
      dat_outfold = data[!infold, ]
      fold_size = nrow(dat_infold)
      #infold_matrix = matrix(rep(infold, m_all), nrow = n)
      weights_outfold = weights[!infold]

      if ('censor' %in% crossfit_models){
        cat('Crossfitting fold:', fold, 'Censor model:', models[2],'\n')
        S_c_infold = predict_S_c(dat_infold, dat_outfold, 'delta_c', 'X2',
                                 time_all, model = censor_model,
                                 weights = weights_outfold,
                                 nodesize = c_nodesize,
                                 ntree = c_ntree)

        nuisance_select = c('S_c_0',  'S_c_1')

        # xx = matrix(c(T,T,T,F,F,F,T,T,T), nrow= 3, byrow = T)
        # yy = matrix(1:12, nrow= 3)
        # replace(yy, xx, matrix(100:106, nrow = 2))
        # nuisance_ls[nuisance_select] = map2(nuisance_ls[nuisance_select],
        #                                     S_c_infold,
        #                                     ~replace(.x, infold_matrix, .y))
        nuisance_ls[nuisance_select] = map2(nuisance_ls[nuisance_select],
                                            S_c_infold,
                                            ~`[<-`(.x, infold,,.y))

      }

      if ('time' %in% crossfit_models){
        cat('Crossfitting fold:', fold, 'Time model:', models[3],'\n')
        Lambda_infold = predict_Lambda(dat_infold, dat_outfold, 'X1', 'status',
                                       time_all, model = time_model,
                                       weights = weights_outfold,
                                       nodesize = t_nodesize,
                                       ntree = t_ntree)
        #print('Lambda has been estimated')

        S_t_infold = predict_S_t(Lambda_infold)

        if (time_model == 'Cox') {
          S_T12_S_T2_infold = predict_S_T12_S_T2_markov(S_t_infold, Lambda_infold)
          nuisance_select = c('Lambda1_0', 'Lambda2_0','Lambda3_0',
                              'S_t_0', 'S_T12_0', 'S_T2_0',
                              'Lambda1_1', 'Lambda2_1', 'Lambda3_1',
                              'S_t_1', 'S_T12_1', 'S_T2_1')

          nuisance_ls[nuisance_select] = map2(nuisance_ls[nuisance_select],
                                              c( Lambda_infold[1:3],
                                                 S_t_infold[1],
                                                 S_T12_S_T2_infold[1:2],
                                                 Lambda_infold[4:6],
                                                 S_t_infold[2],
                                                 S_T12_S_T2_infold[3:4]),
                                              ~`[<-`(.x, infold,,.y))
        }else{
          S_T2_infold = predict_S_T2_nonpara(dat_infold, dat_outfold,
                                             'X2', 'delta2',
                                             time_interest = time_all,
                                             model = time_model,
                                             weights = weights_outfold,
                                             nodesize = t_nodesize,
                                             ntree = t_ntree)
          S_T12_infold = predict_S_T12(S_t_infold, S_T2_infold)
          nuisance_select = c('Lambda1_0', 'Lambda2_0',
                              'S_t_0', 'S_T2_0', 'S_T12_0',
                              'Lambda1_1', 'Lambda2_1',
                              'S_t_1', 'S_T2_1', 'S_T12_1')

          nuisance_ls[nuisance_select] = map2(nuisance_ls[nuisance_select],
                                              c( Lambda_infold[1:2],
                                                 S_t_infold[1],
                                                 S_T2_infold[1],
                                                 S_T12_infold[1],
                                                 Lambda_infold[3:4],
                                                 S_t_infold[2],
                                                 S_T2_infold[2],
                                                 S_T12_infold[2]),
                                              ~`[<-`(.x, infold,,.y))

        }

      }
      if('PS' %in% crossfit_models){
        cat('Crossfitting fold:', fold, 'PS Model:', models[1],'\n')
        PS_infold = predict_PS(dat_infold, dat_outfold,
                               model = PS_model,
                               weights = weights_outfold)
        PS[infold] = PS_infold
      }
    }

  }

  # get nuisance estimates without cross-fitting (i.e., use all the data)
  if (!'censor' %in% crossfit_models){
    cat('Fitting nuisance model on the whole sample;', 'Censor model:', models[2],'\n')
    S_c = predict_S_c(dat_predict = data, dat_train = data,
                      delta_c = 'delta_c',
                      X_c = 'X2',
                      time_interest = time_all, model = censor_model,
                      weights = weights, nodesize = c_nodesize,
                      ntree = c_ntree)

    nuisance_select = c('S_c_0',  'S_c_1')

    nuisance_ls[nuisance_select] = S_c

  }
  if (!'time' %in% crossfit_models) {
    cat('Fitting nuisance model on the whole sample;', 'Time model:', models[3],'\n')
    Lambda = predict_Lambda(
      dat_predict = data,
      dat_train = data,
      time = 'X1',
      status = 'status',
      time_interest = time_all,
      model = time_model,
      weights = weights,
      nodesize = t_nodesize,
      ntree = t_ntree
    )

    S_t = predict_S_t(Lambda)
    if (time_model == 'Cox') {
      S_T12_S_T2 = predict_S_T12_S_T2_markov(S_t, Lambda)
      nuisance_select = c(
        'Lambda1_0',
        'Lambda2_0',
        'Lambda3_0',
        'S_t_0',
        'S_T12_0',
        'S_T2_0',
        'Lambda1_1',
        'Lambda2_1',
        'Lambda3_1',
        'S_t_1',
        'S_T12_1',
        'S_T2_1'
      )
      nuisance_ls[nuisance_select] = c(Lambda[1:3], S_t[1], S_T12_S_T2[1:2], Lambda[4:6], S_t[2], S_T12_S_T2[3:4])
    }
    else{
      S_T2 = predict_S_T2_nonpara(
        data,
        data,
        'X2',
        'delta2',
        time_interest = time_all,
        model = time_model,
        weights = weights,
        nodesize = t_nodesize,
        ntree = t_ntree
      )
      S_T12 = predict_S_T12(S_t, S_T2)
      nuisance_select = c(
        'Lambda1_0',
        'Lambda2_0',
        'S_t_0',
        'S_T2_0',
        'S_T12_0',
        'Lambda1_1',
        'Lambda2_1',
        'S_t_1',
        'S_T2_1',
        'S_T12_1'
      )

      nuisance_ls[nuisance_select] = c(Lambda[1:2], S_t[1], S_T2[1], S_T12[1], Lambda[3:4], S_t[2], S_T2[2], S_T12[2])

    }
  }

  if (!'PS' %in% crossfit_models){
    cat('Fitting nuisance model on the whole sample;', 'PS model:', models[1],'\n')
    PS = predict_PS(data, data,
                    model = PS_model,
                    weights = weights)}

  # truncate the IPW weights
  # in addition to lower bound PS, also need to lower bound upper bound PS since
  # we need to calculate (1-A)/(1-PS) for those who did not get the treatment
  PS[PS < PS_min] = PS_min
  PS[PS > 1-PS_min] = 1-PS_min

  #nuisance_select = str_detect(nuisuance_names, '^S')
  nuisance_select = c('S_c_0','S_c_1')
  nuisance_ls[nuisance_select] = map(nuisance_ls[nuisance_select], ~replace(.x, .x<S_min, S_min))

  # add a small number to survival function to avoid 0
  nuisance_select = c('S_t_0','S_t_1', 'S_T2_0', 'S_T2_1', 'S_T12_0', 'S_T12_1')
  nuisance_ls[nuisance_select] = map(nuisance_ls[nuisance_select], ~replace(.x, .x<0.01, 0.01))

  nuisance_ls$PS = PS
  if (any(is.na(unlist(nuisance_ls)))) stop('Nuisance estimates contain NA. Please check!')
  # if (any(is.na(unlist(nuisance_ls)))) {
  #   return(weights)
  #   stop('Nuisance estimates contain NA. Please check!')}
  #

  end_time = Sys.time()
  runtime = round(end_time - start_time,2)


  list(nuisance_ls = nuisance_ls,
       runtime = runtime
  )
}
