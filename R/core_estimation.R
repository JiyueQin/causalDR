
sourceCpp("get_integral.cpp")



####################################################################################
#  Function get_CIF - get CIF from S_t and Lambda
#  return: a vector with length: number of time points of interest
#
####################################################################################

get_CIF = function(S_t, input, type = 'dLambda'){
  if (type == 'Lambda'){
    dLambda  = input - lag(input)
    dLambda[1] = 0
  }else{
    dLambda = input
  }
  cumsum(S_t*dLambda)
  # todo: need to make sure the output is less than 1
}




####################################################################################
#  Function get_terminal_risk - get terminal risk from Lambda 3 in estimate df and T1 to condition on
#  return: a vector with length: number of time points of interest
#
####################################################################################

get_terminal_risk = function(est, time1_interest, time_all,
                             time1_interest_single_idx, simple_name = T,
                             add_contrast =F){
  time1_interest_single = time1_interest[time1_interest_single_idx]
  # find the index of the largest time point prior to the target time point
  time1_interest_match_idx = length(time_all[time_all<=time1_interest_single])
  if (time1_interest_match_idx == 0){
    stop('time1 specified is too small! Please check!')
  }else{
    est_Lambda3_given_T1 = c(est$est_Lambda3_0[time1_interest_match_idx],
                             est$est_Lambda3_1[time1_interest_match_idx])
  }

  if(simple_name) {suffix = NULL}
  else{suffix = paste0('.T1.',time1_interest_single_idx)}

  n_cols = ncol(est)
  out = est %>%
    mutate(risk0 = ifelse(time<= time1_interest_single, 0,
                          1-exp(est_Lambda3_given_T1[1]-est_Lambda3_0))) %>%
    mutate(risk1 = ifelse(time<= time1_interest_single, 0,
                          1-exp(est_Lambda3_given_T1[2]-est_Lambda3_1)))
  if(add_contrast){
    out = out %>%
      mutate(diff = risk1-risk0,
             # note, NaN can occur with 0/0
             ratio = risk1/risk0)
    colnames(out)[(n_cols+1):(n_cols+4)] = paste0('est_risk_', c('0','1', 'diff', 'ratio'), suffix)
  }else{
    colnames(out)[(n_cols+1):(n_cols+2)] = paste0('est_risk_', c('0','1'), suffix)
  }

  out
}


####################################################################################
#  Function get_naive_est - get naive estimator from cleaned data
#  return: a dataframe
#
####################################################################################

get_naive_est = function(data, time_all,estimand, time_interest,
                         time1_interest,weights, add_trans, add_contrast){
  start_time = Sys.time()

  if(is.null(time_interest)){
    time_interest = time_all
    flag_specify_time_interest = F
  }else{
    # extrapolate with the closest time point in the sample
    # need to compute dLambda at all the time points before the largest time point of interest
    flag_specify_time_interest = T
    time_interest_user = sort(time_interest)
    time_interest_match_idx =  map_dbl(time_interest_user, ~length(time_all[time_all<=.x]))
  }

  m_time1 = length(time1_interest)

  time_ls1 = map(0:1, ~ data %>% filter(A==.x) %>% pull(X1) %>% unique() %>% sort())
  time_ls3 = map(0:1, ~ data %>% filter(A==.x, delta1==1) %>% pull(X2) %>% unique() %>% sort())
  time_ls = c(rep(time_ls1, 2), time_ls3)

  # get naive estimates which are invalid under the assumptions
  est_Lambda_naive_df_ls = c(map(0:1, ~coxph(Surv(X1, status ==1)~1, data %>% filter(A==.x),
                                             weights = weights[data$A==.x]) %>%
                                   basehaz(centered = F)),
                             map(0:1, ~coxph(Surv(X1, status ==2)~1, data %>% filter(A==.x),
                                             weights = weights[data$A==.x]) %>%
                                   basehaz(centered = F)),
                             map(0:1, ~coxph(Surv(X1, X2, delta1*delta2)~1,
                                             data %>% filter(A==.x, delta1==1),
                                             weights = weights[data$A==.x & data$delta1]) %>%
                                   basehaz(centered = F)))
  names(est_Lambda_naive_df_ls) = paste0('est_Lambda', c(1,1,2,2,3,3),'_',c(0,1,0,1,0,1))
  time_ls = map(est_Lambda_naive_df_ls, ~.x %>% pull(time))

  est_Lambda_naive_ls = map(est_Lambda_naive_df_ls, ~.x %>% pull(hazard))

  est_Lambda_naive_ls = map2(est_Lambda_naive_ls, time_ls,
                             ~fill_time_matrix(matrix(.x, nrow=1),
                                               .y, time_all, 'Lambda') %>% as.vector)

  est_S_t_0_naive = exp(-est_Lambda_naive_ls$est_Lambda1_0-est_Lambda_naive_ls$est_Lambda2_0)
  est_S_t_1_naive = exp(-est_Lambda_naive_ls$est_Lambda1_1-est_Lambda_naive_ls$est_Lambda2_1)
  est_naive = cbind(data.frame(time = time_all), est_Lambda_naive_ls, est_S_t_0 = est_S_t_0_naive, est_S_t_1 = est_S_t_1_naive)
  if (estimand == 'all'|'CIF1' %in% estimand){
    est_CIF1_0_naive= get_CIF(est_S_t_0_naive, est_Lambda_naive_ls$est_Lambda1_0, 'Lambda')
    est_CIF1_1_naive = get_CIF(est_S_t_1_naive, est_Lambda_naive_ls$est_Lambda1_1, 'Lambda')

    est_naive = cbind(est_naive, est_CIF1_0 = est_CIF1_0_naive,
                      est_CIF1_1 = est_CIF1_1_naive)
    if(add_contrast){
      est_naive = mutate(est_naive,
                         est_CIF1_diff = est_CIF1_1 -est_CIF1_0,
                         est_CIF1_ratio = est_CIF1_1/est_CIF1_0)
    }

  }

  if (estimand == 'all' |'CIF2' %in% estimand){
    est_CIF2_0_naive = get_CIF(est_S_t_0_naive, est_Lambda_naive_ls$est_Lambda2_0, 'Lambda')
    est_CIF2_1_naive = get_CIF(est_S_t_1_naive, est_Lambda_naive_ls$est_Lambda2_1, 'Lambda')
    est_CIF2_diff_naive = est_CIF2_1_naive - est_CIF2_0_naive
    est_naive = cbind(est_naive,
                      est_CIF2_0 = est_CIF2_0_naive,
                      est_CIF2_1 = est_CIF2_1_naive)
    if(add_contrast){
      est_naive = mutate(est_naive,
                         est_CIF2_diff = est_CIF2_1 - est_CIF2_0,
                         est_CIF2_ratio = est_CIF2_1/est_CIF2_0)
    }
  }

  if (estimand == 'all' |'terminal_risk' %in% estimand){
    if (m_time1==1) est_naive = get_terminal_risk(est_naive, time1_interest,
                                                  time_all,1, simple_name = T,
                                                  add_contrast = add_contrast)
    else{
      for (j in 1:m_time1){
        est_naive = get_terminal_risk(est_naive, time1_interest, time_all,
                                      j, simple_name = F,
                                      add_contrast = add_contrast)
      }
    }
  }
  if (flag_specify_time_interest){
    est_naive = est_naive[time_interest_match_idx, ] %>%
      mutate(time = time_interest_user)

  }

  rownames(est_naive) = NULL
  if(add_trans) est_naive = transform_est(est_naive)

  end_time = Sys.time()
  runtime = round(end_time - start_time,2)

  list(est = est_naive,
       runtime = runtime)

}

# add or minus a small epsilon so the probability is not exactly 0 or 1
transform_est = function(est, epsilon = 1e-6){
  log_log = function(x){
    y = case_when(x ==1 ~ x-epsilon,
                  x==0 ~ x+epsilon,
                  .default = x)
    log(-log(y))
  }
  # transform estimates so that their values are not restricted for CI
  log_log_vars = est %>%
    select(matches('CIF[12]_[01]'), matches('risk_[01]')) %>% names

  log_vars = est %>%
    select(contains('ratio')) %>% names()

  est %>%
    mutate_at(log_log_vars, list(trans = ~log_log(.))) %>%
    mutate_at(log_vars, list(trans = log))

}




# some arguments are redundant because they are already defined in the main function
causalAIPCW_est = function(data,
                           n,
                           time_all,
                           estimand,
                           time_interest,
                           time1_interest,
                           PS_model = 'logit', censor_model = 'Cox',
                           time_model = 'Cox', c_nodesize = 15,
                           c_ntree = 500,  t_nodesize = 15, t_ntree = 500,
                           k = 5, crossfit = 'auto', weights = NULL,
                           shuffle = T,
                           PS_min = 0.1, S_min = 0.05, add_IPW = F,
                           name = NULL, add_trans = F, silent = T,
                           add_contrast = F
){
  start_time = Sys.time()

  if(is.null(time_interest)){
    time_interest = time_all
    flag_specify_time_interest = F
  }else{
    # extrapolate with the closest time point in the sample
    # need to compute dLambda at all the time points before the largest time point of interest
    flag_specify_time_interest = T
    time_interest_user = sort(time_interest[time_interest>0])
    time_interest_match_idx =  map_dbl(time_interest_user, ~length(time_all[time_all<=.x]))
    time_interest = time_all[1:max(time_interest_match_idx)]
  }


  m = length(time_interest)
  m_all = length(time_all)
  m_time1 = length(time1_interest)

  if(is.null(name)){
    if(add_IPW){name = c('causalAIPCW','IPW')}
    else{name = 'causalAIPCW'}
  }

  if (is.null(weights)) weights = rep(1, n)
  # need nuisance estimates for all time points due to the involvement of nuisance at X2
  nuisance_out = nuisance_est(data,
                              n,
                              time_all,
                              PS_model , censor_model ,
                              time_model, c_nodesize,
                              c_ntree,  t_nodesize, t_ntree,
                              k, crossfit, weights,
                              shuffle,
                              PS_min, S_min)
  nuisance_ls = nuisance_out$nuisance_ls
  if (!silent) {cat('Estimate nuisance parameters: Completed in', as.numeric(nuisance_out$runtime, units = "secs"), "seconds!\n")}
  #nuisance_runtime = nuisance_out$runtime
  #print(nuisance_runtime)

  ## note dS_c are negative
  dS_c_0 = get_jump_size(nuisance_ls$S_c_0, 'S')
  dS_c_1 = get_jump_size(nuisance_ls$S_c_1, 'S')

  # each matrix is of n*m
  P_at_risk_ls_Lambda1 = list(nuisance_ls$S_t_0[,1:m], nuisance_ls$S_t_1[,1:m])

  P_at_risk_ls_Lambda3 = list(nuisance_ls$S_T12_0[,1:m], nuisance_ls$S_T12_1[,1:m])

  # for Cox model, use lambda3
  # if(time_model == 'Cox'){
  #   Lambda_ls = nuisance_ls[c('Lambda1_0', 'Lambda2_0', 'Lambda3_0',
  #                             'Lambda1_1', 'Lambda2_1', 'Lambda3_1')]
  #   dLambda_ls = map(Lambda_ls, get_jump_size)
  #   names(dLambda_ls) = c('dLambda1_0', 'dLambda2_0', 'dLambda3_0',
  #                         'dLambda1_1', 'dLambda2_1', 'dLambda3_1')
  #   dP_event_ls_Lambda3 =list(nuisance_ls$S_T12_0*dLambda_ls$dLambda3_0,
  #                             nuisance_ls$S_T12_1*dLambda_ls$dLambda3_1)
  # }else{
  #   Lambda_ls = nuisance_ls[c('Lambda1_0', 'Lambda2_0',
  #                             'Lambda1_1', 'Lambda2_1')]
  #   dLambda_ls = map(Lambda_ls, get_jump_size)
  #   names(dLambda_ls) = c('dLambda1_0', 'dLambda2_0',
  #                         'dLambda1_1', 'dLambda2_1')
  #
  #   CIF1_ls = predict_CIF1(list(nuisance_ls$S_t_0, nuisance_ls$S_t_1),
  #                          list(Lambda_ls$Lambda1_0, Lambda_ls$Lambda1_1))
  #   P_event_ls_Lambda3 =list(CIF1_ls[[1]]-nuisance_ls$S_T12_0,
  #                            CIF1_ls[[2]]-nuisance_ls$S_T12_1)
  #   dP_event_ls_Lambda3 = map(P_event_ls_Lambda3, get_jump_size)
  #   # force the jumps to be positive
  #   dP_event_ls_Lambda3 = map(dP_event_ls_Lambda3, ~replace(.x, .x<0, 0))
  # }


  Lambda_ls = nuisance_ls[c('Lambda1_0', 'Lambda2_0',
                            'Lambda1_1', 'Lambda2_1')]
  dLambda_ls = map(Lambda_ls, get_jump_size)
  names(dLambda_ls) = c('dLambda1_0', 'dLambda2_0',
                        'dLambda1_1', 'dLambda2_1')
  dP_event_ls_Lambda1 =list(nuisance_ls$S_t_0[,1:m]*dLambda_ls$dLambda1_0[,1:m],
                            nuisance_ls$S_t_1[,1:m]*dLambda_ls$dLambda1_1[,1:m])
  dP_event_ls_Lambda2 =list(nuisance_ls$S_t_0[,1:m]*dLambda_ls$dLambda2_0[,1:m],
                            nuisance_ls$S_t_1[,1:m]*dLambda_ls$dLambda2_1[,1:m])


  #CIF1_ls = predict_CIF1(list(nuisance_ls$S_t_0, nuisance_ls$S_t_1),
  #                         list(Lambda_ls$Lambda1_0, Lambda_ls$Lambda1_1))
  #P_event_ls_Lambda3 =list(CIF1_ls[[1]]-nuisance_ls$S_T12_0,
  #                           CIF1_ls[[2]]-nuisance_ls$S_T12_1)
  #dP_event_ls_Lambda3 = map(P_event_ls_Lambda3, get_jump_size)

  ## note dS_T2 are negative

  dS_T2_0 = get_jump_size(nuisance_ls$S_T2_0, 'S')
  dS_T2_1 = get_jump_size(nuisance_ls$S_T2_1, 'S')

  dP_event_ls_Lambda3 = list(-dS_T2_0[,1:m]-nuisance_ls$S_t_0[,1:m]*dLambda_ls$dLambda2_0[,1:m],
                             -dS_T2_1[,1:m]-nuisance_ls$S_t_1[,1:m]*dLambda_ls$dLambda2_1[,1:m])
  # force the jumps to be positive
  dP_event_ls_Lambda3 = map(dP_event_ls_Lambda3, ~replace(.x, .x<0, 0))


  # compute estimates of parameters of interest
  # est_dLambda1_1: cause-specific hazard of illness for potential outcomes with treatment
  # est_dLambda2_1: cause-specific hazard of death without illness for potential outcomes with treatment
  # est_dLambda3_1: conditional hazard of death given illness for potential outcomes with treatment
  # est_dLambda1_0 = est_dLambda2_0 =
  #   est_dLambda1_1 = est_dLambda2_1 =
  #   est_dLambda3_0 = est_dLambda3_1 = numeric(m)
  #
  # if (add_IPW){
  #   est_dLambda1_0_IPW = est_dLambda2_0_IPW =
  #     est_dLambda1_1_IPW = est_dLambda2_1_IPW =
  #     est_dLambda3_0_IPW = est_dLambda3_1_IPW = numeric(m)
  #
  # }
  #
  X2_index = match(data$X2, time_all)

  #S_c1_x_star = get_vec_from_mat(nuisance_ls$S_c1, X_star, time_interest)
  nuisance_X2_select = c('S_c_0', 'S_c_1', 'S_t_0', 'S_t_1', 'S_T2_0', 'S_T2_1')

  nuisance_X2_ls = map(nuisance_ls[nuisance_X2_select], ~get_vec_from_mat(.x, X2_index))
  names(nuisance_X2_ls) = paste0(nuisance_X2_select, '_X2')

  # get_integral = function(X_star, nuisance_name){
  #   nuisance_name0 = paste0(nuisance_name, '_0')
  #   nuisance_name1 = paste0(nuisance_name, '_1')
  #   nuisance_name0_X2 = paste0(nuisance_name, '_0_X2')
  #   nuisance_name1_X2 = paste0(nuisance_name, '_1_X2')
  #   indicator_X_star = do.call(cbind, map(time_all, ~1*(.x<=X_star)))
  #   integral_0 = (1-data$delta2)*(data$X2<=X_star)/nuisance_X2_ls[[nuisance_name0_X2]]/nuisance_X2_ls[['S_c_0_X2']]+
  #     rowSums(indicator_X_star/nuisance_ls[[nuisance_name0]]/(nuisance_ls$S_c_0^2)*dS_c_0)
  #   integral_1 = (1-data$delta2)*(data$X2<=X_star)/nuisance_X2_ls[[nuisance_name1_X2]]/nuisance_X2_ls[['S_c_1_X2']]+
  #     rowSums(indicator_X_star/nuisance_ls[[nuisance_name1]]/(nuisance_ls$S_c_1^2)*dS_c_1)
  #   list(integral_0 = integral_0 , integral_1 = integral_1)
  # }
  #

  integrand_core = list(1/(nuisance_ls$S_c_0[,1:m]^2)*dS_c_0[,1:m], 1/(nuisance_ls$S_c_1[,1:m]^2)*dS_c_1[,1:m])
  # matrix: n*m
  X_star_mat_Lambda1 = outer(data$X1, time_interest, pmin)
  X_star_mat_Lambda3 = outer(data$X2, time_interest, pmin)


  integral_mat_0_Lambda1 = get_integral_cpp(time_interest, X_star_mat_Lambda1, 1/nuisance_ls[['S_t_0']][,1:m]*integrand_core[[1]])
  integral_mat_1_Lambda1 = get_integral_cpp(time_interest, X_star_mat_Lambda1, 1/nuisance_ls[['S_t_1']][,1:m]*integrand_core[[2]])
  integral_mat_0_Lambda3 = get_integral_cpp(time_interest, X_star_mat_Lambda3, 1/nuisance_ls[['S_T2_0']][,1:m]*integrand_core[[1]])
  integral_mat_1_Lambda3 = get_integral_cpp(time_interest, X_star_mat_Lambda3, 1/nuisance_ls[['S_T2_1']][,1:m]*integrand_core[[2]])
  # # for Lambda1 and Lambda2
  # integral_mat_0_Lambda1 =  integral_mat_1_Lambda1 = integral_mat_0_Lambda3 = integral_mat_1_Lambda3 = matrix(NA_real_, nrow = n, ncol = m)
  #
  # for (j in 1:m){
  #   time = time_interest[j]
  #   if(j ==1 & time != min(data$X1)) stop('The first time point is not the smallest time point! Please check!')
  #   X_star1 = pmin(data$X1, time)
  #   X_star3 = pmin(data$X2, time)
  #   integral_ls_Lambda1 = get_integral(X_star1, 'S_t')
  #   integral_ls_Lambda3 = get_integral(X_star3, 'S_T2')
  #   integral_mat_0_Lambda1[,j] = integral_ls_Lambda1[[1]]
  #   integral_mat_1_Lambda1[,j] = integral_ls_Lambda1[[2]]
  #   integral_mat_0_Lambda3[,j] = integral_ls_Lambda3[[1]]
  #   integral_mat_1_Lambda3[,j] = integral_ls_Lambda3[[2]]
  # }

  integral_mat_0_Lambda1 = integral_mat_0_Lambda1 + (1-data$delta2)*(data$X2<=X_star_mat_Lambda1)/nuisance_X2_ls[['S_t_0_X2']]/nuisance_X2_ls[['S_c_0_X2']]
  integral_mat_1_Lambda1 = integral_mat_1_Lambda1 + (1-data$delta2)*(data$X2<=X_star_mat_Lambda1)/nuisance_X2_ls[['S_t_1_X2']]/nuisance_X2_ls[['S_c_1_X2']]
  integral_mat_0_Lambda3 = integral_mat_0_Lambda3 + (1-data$delta2)*(data$X2<=X_star_mat_Lambda3)/nuisance_X2_ls[['S_T2_0_X2']]/nuisance_X2_ls[['S_c_0_X2']]
  integral_mat_1_Lambda3 = integral_mat_1_Lambda3 + (1-data$delta2)*(data$X2<=X_star_mat_Lambda3)/nuisance_X2_ls[['S_T2_1_X2']]/nuisance_X2_ls[['S_c_1_X2']]

  #integral_all = list(integral_mat_0_Lambda1, integral_mat_1_Lambda1, integral_mat_0_Lambda3, integral_mat_1_Lambda3)
  # compare
  # map(1:6, ~all.equal(est_dLambda_ls[[.x]], est_dLambda_ls_raw[[.x]])
  get_weights2 = function(integral_ls){
    #weights2_0 = (1+ (1-data$A)/(1-nuisance_ls$PS)*(integral_ls[[1]]-1))*multiplier_ls[[1]][,time_index]
    weights2_0 = (1+ (1-data$A)/(1-nuisance_ls$PS)*(integral_ls[[1]]-1))
    weights2_1 = (1+ data$A/nuisance_ls$PS*(integral_ls[[2]]-1))
    list(weights2_0 = weights2_0, weights2_1 = weights2_1)
  }

  weights2_mat_ls_Lambda1 = get_weights2(list(integral_mat_0_Lambda1,
                                              integral_mat_1_Lambda1))
  weights2_mat_ls_Lambda3 = get_weights2(list(integral_mat_0_Lambda3,
                                              integral_mat_1_Lambda3))

  # each matrix is of n*m
  weights1_mat_ls = list(weights1_0 = (1-data$A)/(1-nuisance_ls$PS)/nuisance_ls$S_c_0[,1:m],
                         weights1_1 = data$A/nuisance_ls$PS/nuisance_ls$S_c_1[,1:m])


  illness = data$delta1*do.call(cbind, map(time_interest, ~(data$X1==.x)))
  death = do.call(cbind, map(time_interest, ~data$X2==.x))

  death_without_illness = (1-data$delta1)*data$delta2*death
  death_after_illness = data$delta1*data$delta2*death
  at_risk0 = do.call(cbind, map(time_interest, ~1*(.x<=data$X1)))
  at_risk1 = data$delta1*do.call(cbind, map(time_interest, ~data$X1<.x & .x<= data$X2))
  weights1_illness = list(weights1_mat_ls[[1]]*illness, weights1_mat_ls[[2]]*illness)
  weights1_death_without_illness = list(weights1_mat_ls[[1]]*death_without_illness, weights1_mat_ls[[2]]*death_without_illness)
  weights1_death_after_illness = list(weights1_mat_ls[[1]]*death_after_illness, weights1_mat_ls[[2]]*death_after_illness)
  weights1_at_risk0 = list(weights1_mat_ls[[1]]*at_risk0, weights1_mat_ls[[2]]*at_risk0)
  weights1_at_risk1 = list(weights1_mat_ls[[1]]*at_risk1, weights1_mat_ls[[2]]*at_risk1)

  est_dLambda1_0 = weights %*%(weights1_illness[[1]] + weights2_mat_ls_Lambda1[[1]]*dP_event_ls_Lambda1[[1]])/weights %*%(weights1_at_risk0[[1]] + weights2_mat_ls_Lambda1[[1]]*P_at_risk_ls_Lambda1[[1]])
  est_dLambda1_1 = weights %*%(weights1_illness[[2]] + weights2_mat_ls_Lambda1[[2]]*dP_event_ls_Lambda1[[2]])/weights %*%(weights1_at_risk0[[2]] + weights2_mat_ls_Lambda1[[2]]*P_at_risk_ls_Lambda1[[2]])
  est_dLambda2_0 = weights %*%(weights1_death_without_illness[[1]] + weights2_mat_ls_Lambda1[[1]]*dP_event_ls_Lambda2[[1]])/weights %*%(weights1_at_risk0[[1]] + weights2_mat_ls_Lambda1[[1]]*P_at_risk_ls_Lambda1[[1]])
  est_dLambda2_1 = weights %*%(weights1_death_without_illness[[2]] + weights2_mat_ls_Lambda1[[2]]*dP_event_ls_Lambda2[[2]])/weights %*%(weights1_at_risk0[[2]] + weights2_mat_ls_Lambda1[[2]]*P_at_risk_ls_Lambda1[[2]])
  est_dLambda3_0 = weights %*%(weights1_death_after_illness[[1]] + weights2_mat_ls_Lambda3[[1]]*dP_event_ls_Lambda3[[1]])/weights %*%(weights1_at_risk1[[1]] + weights2_mat_ls_Lambda3[[1]]*P_at_risk_ls_Lambda3[[1]])
  est_dLambda3_1 = weights %*%(weights1_death_after_illness[[2]] + weights2_mat_ls_Lambda3[[2]]*dP_event_ls_Lambda3[[2]])/weights %*%(weights1_at_risk1[[2]] + weights2_mat_ls_Lambda3[[2]]*P_at_risk_ls_Lambda3[[2]])

  est_dLambda_ls = list(est_dLambda1_0, est_dLambda2_0, est_dLambda3_0,
                        est_dLambda1_1, est_dLambda2_1,
                        est_dLambda3_1)

  if(anyNA(unlist(est_dLambda_ls))){
    warning(paste0('Method name:', name[1], ', NA estimates of dLambda have occurred! Please Check!'))}
  else if (!silent){print("dLambda has been estimated sucessfully!")}

  # map(est_dLambda_ls, anyNA)

  if (add_IPW){
    est_dLambda1_0_IPW =  weights %*%(weights1_illness[[1]])/ weights %*%(weights1_at_risk0[[1]])
    est_dLambda1_1_IPW =  weights %*%(weights1_illness[[2]])/ weights %*%(weights1_at_risk0[[2]])
    est_dLambda2_0_IPW =  weights %*%(weights1_death_without_illness[[1]])/ weights %*%(weights1_at_risk0[[1]])
    est_dLambda2_1_IPW =  weights %*%(weights1_death_without_illness[[2]])/ weights %*%(weights1_at_risk0[[2]])
    est_dLambda3_0_IPW =  weights %*%(weights1_death_after_illness[[1]])/ weights %*%(weights1_at_risk1[[1]])
    est_dLambda3_1_IPW =  weights %*%(weights1_death_after_illness[[2]])/ weights %*%(weights1_at_risk1[[2]])
    est_dLambda_ls_IPW = list(est_dLambda1_0_IPW, est_dLambda2_0_IPW, est_dLambda3_0_IPW,
                              est_dLambda1_1_IPW, est_dLambda2_1_IPW,
                              est_dLambda3_1_IPW)
  }

  risk_contrasts = function(est_dLambda_ls, name, add_trans=F,
                            add_contrast = F){
    neg_time_index_ls = map(est_dLambda_ls, ~which(.x<0))

    # in AIPW: force est_dLambda to be non-negative
    # in IPW:  for large time points, 0/0 = NaN can occur
    est_dLambda_ls = map(est_dLambda_ls, ~ifelse(.x<0|is.na(.x), 0, .x))
    names(est_dLambda_ls) = paste0('est_dLambda', c('1_0', '2_0', '3_0', '1_1', '2_1', '3_1'))
    est_Lambda_ls = map(est_dLambda_ls, cumsum)
    names(est_Lambda_ls) = paste0('est_Lambda', c('1_0', '2_0', '3_0', '1_1', '2_1', '3_1'))

    est_S_t_0 = exp(-est_Lambda_ls$est_Lambda1_0-est_Lambda_ls$est_Lambda2_0)
    est_S_t_1 = exp(-est_Lambda_ls$est_Lambda1_1-est_Lambda_ls$est_Lambda2_1)

    est = cbind(data.frame(time = time_interest),est_Lambda_ls, est_S_t_0, est_S_t_1)

    if (any(estimand == 'all')) estimand_vec = c('CIF1', 'CIF2', 'CIF3')
    else estimand_vec = estimand
    if ('CIF1' %in% estimand_vec){
      est_CIF1_1 = get_CIF(est_S_t_1, est_dLambda_ls$est_dLambda1_1)
      est_CIF1_0 = get_CIF(est_S_t_0, est_dLambda_ls$est_dLambda1_0)
      est = cbind(est, est_CIF1_0, est_CIF1_1)
      if(add_contrast){
        est_CIF1_diff = est_CIF1_1 - est_CIF1_0
        est_CIF1_ratio = est_CIF1_1/est_CIF1_0
        est = cbind(est, est_CIF1_diff, est_CIF1_ratio)
      }

    }

    if ('CIF2' %in% estimand_vec){
      est_CIF2_1 = get_CIF(est_S_t_1, est_dLambda_ls$est_dLambda2_1)
      est_CIF2_0 = get_CIF(est_S_t_0, est_dLambda_ls$est_dLambda2_0)
      est = cbind(est, est_CIF2_0, est_CIF2_1)
      if(add_contrast){
        est_CIF2_diff = est_CIF2_1 - est_CIF2_0
        est_CIF2_ratio = est_CIF2_1 / est_CIF2_0
        est = cbind(est, est_CIF2_diff, est_CIF2_ratio)
      }
    }
    #return(list(est, nuisance_ls, crossfit_models))
    # if(max(est_CIF1_1+est_CIF2_1, est_CIF1_0+est_CIF2_0)>1) {
    #   print('CIF sum >1 occurred')}



    if ('CIF3' %in% estimand_vec){
      m_time1 = length(time1_interest)
      if (m_time1==1) est = get_terminal_risk(est, time1_interest, time_all,
                                              1, simple_name = T,
                                              add_contrast = add_contrast)
      else{
        for (j in 1:m_time1){
          est = get_terminal_risk(est, time1_interest, time_all,
                                  j, simple_name = F,
                                  add_contrast= add_contrast)
        }
      }
    }

    # only keep results of time points of interest
    if (flag_specify_time_interest) {
      est = est[time_interest_match_idx, ] %>%
        mutate(time = time_interest_user)
    }

    rownames(est) = NULL
    if (add_trans) est = transform_est(est)
    out = list(est = cbind(method = name, est),
               neg_time_index_ls = neg_time_index_ls)
  }


  risk_out = risk_contrasts(est_dLambda_ls, name[1], add_trans = add_trans,
                            add_contrast = add_contrast)
  est = risk_out$est
  if(add_IPW){
    risk_out_IPW = risk_contrasts(est_dLambda_ls_IPW, name[2],
                                  add_trans = add_trans,
                                  add_contrast = add_contrast)
    est = rbind(est, risk_out_IPW$est)
  }


  end_time = Sys.time()
  runtime = round(end_time - start_time,2)


  list(est = est,
       nuisance = nuisance_ls,
       runtime = runtime,
       neg_time_index_ls = risk_out$neg_time_index_ls,
       call = match.call()
  )
}






####################################################################################
#  Function get_se - summarize the results from bootstrap to get SE
#  return: a data frame with SE for each estimate
#
###################################################################################

get_se = function(boot_results) {
  boot_results %>%
    pivot_longer(starts_with('est_'),
                 names_to = 'estimand',
                 values_to = 'est',
                 names_prefix = 'est_') %>%
    group_by(method, time, estimand) %>%
    summarise(se = sd(est, na.rm = T), .groups = 'drop')

}


####################################################################################
#  Function get_CI - get confidence intervals by bootstrap and transforming
#  return: a data frame with CI for each estimate
#
###################################################################################

get_CI = function(est, se_results, CI_raw = T, alpha = 0.05){
  q = qnorm(1 - alpha / 2)
  est_se = suppressMessages(left_join(est, se_results)) %>%
    mutate(lower = est - q * se, upper = est + q * se)


  est_se_trans = est_se %>% filter(str_detect(estimand, '_trans')) %>%
    mutate(estimand = str_remove(estimand, '_trans')) %>%
    select(-c(est, se)) %>%
    rename_at(vars(lower, upper), ~paste0(., '_trans'))

  out = suppressMessages(right_join(est_se_trans,
                                    est_se %>%
                                      filter(!str_detect(estimand, '_trans'))))

  if(CI_raw){out = out %>%
    mutate(lower_raw = lower, upper_raw =upper)}

  out %>%
    mutate(lower = case_when(
      str_detect(estimand, 'CIF[12]_[01]')|str_detect(estimand, 'risk_[01]') ~
        exp(-exp(upper_trans)),
      str_detect(estimand, 'ratio') ~ exp(lower_trans),
      .default = lower),
      upper = case_when(
        str_detect(estimand, 'CIF[12]_[01]')|str_detect(estimand, 'risk_[01]') ~
          exp(-exp(lower_trans)),
        str_detect(estimand, 'ratio') ~ exp(upper_trans),
        .default = upper)) %>%
    select(-c(lower_trans, upper_trans)) %>%
    arrange(method, time, estimand)
}


#' Causal estimation for semi-competing risks via augmented IPCW (AIPCW)
#'
#' \code{causalAIPCW()} estimates causal estimands in semi-competing risks settings
#' with covariate-dependent censoring. The method involves
#' (i) a treatment model (propensity score), (ii) a censoring model, and
#' (iii) an event-time model. The function supports parametric/semi-parametric models
#' (e.g., logistic regression, Cox PH) and machine-learning models (e.g., random forests, GBM,
#' random survival forests, splines), with optional cross-fitting and bootstrap inference.
#'
#' The input data are assumed to contain two observed time variables:
#' \itemize{
#'   \item \code{X1}: observed time for the non-terminal event (e.g., illness),
#'   \item \code{X2}: observed time for the terminal event (e.g., death),
#' }
#' with indicators \code{delta1} and \code{delta2}. Treatment assignment is binary (\code{A}),
#' and baseline covariates are supplied via \code{covars}.
#'
#' Supported estimands include cumulative incidence functions (CIFs) for key event types:
#' \itemize{
#'   \item \code{"CIF1"}: risk for the non-terminal event,
#'   \item \code{"CIF2"}: risk for the terminal event without prior non-terminal event,
#'   \item \code{"CIF3"}: risk of the terminal event
#'         conditional on surviving past a specified non-terminal time (requires \code{time1_interest}),
#'   \item \code{"all"}: compute all estimands above.
#' }
#'
#' @param data A data.frame containing the variables specified by \code{X1}, \code{X2},
#'   \code{delta1}, \code{delta2}, \code{treatment}, and \code{covars}.
#'
#' @param X1,X2 Character scalars giving the column names for the observed time variables.
#'
#' @param delta1,delta2 Character scalars giving the column names for event indicators.
#'   \code{delta1 = 1} indicates the non-terminal event is observed at \code{X1}. \code{delta2 = 1}
#'   indicates the terminal event is observed at \code{X2}. Right censoring is inferred from \code{delta2}.
#'
#' @param treatment Character scalar giving the column name of the binary treatment indicator.
#'
#' @param covars Character vector of baseline covariate column names used in nuisance models.
#'
#' @param estimand Character vector specifying which estimand(s) to compute. Must be one or more of
#'   \code{"all"}, \code{"CIF1"}, \code{"CIF2"}, \code{"CIF3"}.
#'
#' @param tau Optional numeric truncation time. If supplied and \code{tau < max(X2)}, times larger than
#'   \code{tau} are truncated to \code{tau} and indicators are set to 0 for events beyond \code{tau}.
#'
#' @param time_interest Optional numeric vector of time points at which to evaluate estimands. If \code{NULL},
#'   uses all unique observed times from \code{X1} and \code{X2}.
#'
#' @param time1_interest Optional numeric vector of time points to condition on the non-terminal process.
#'   Required when \code{estimand} includes \code{"terminal_risk"} (or \code{"all"}).
#'
#' @param PS_model Model for the propensity score \eqn{P(A=1 \mid Z)}. Must be one of
#'   \code{"logit"}, \code{"logit.quad"}, \code{"logit.small"}, \code{"RF"}, \code{"gbm"}.
#'
#' @param censor_model Model for the censoring distribution used to construct IPC weights. Must be one of
#'   \code{"Cox"}, \code{"Cox.quad"}, \code{"Cox.small"}, \code{"RSF"}, \code{"Spline"}.
#'
#' @param time_model Model for event-time nuisance components used for augmentation. Must be one of
#'   \code{"Cox"}, \code{"RSF"}, \code{"Spline"}.
#'
#' @param c_nodesize Integer node size for censoring random survival forest (\code{censor_model = "RSF"}).
#'
#' @param c_ntree Integer number of trees for censoring random survival forest (\code{censor_model = "RSF"}).
#'
#' @param t_nodesize Integer node size for event-time random survival forest (\code{time_model = "RSF"}).
#'
#' @param t_ntree Integer number of trees for event-time random survival forest (\code{time_model = "RSF"}).
#'
#' @param k Integer number of folds for cross-fitting. Default is \code{5}.
#'
#' @param crossfit Cross-fitting control. Default \code{"auto"} enables cross-fitting when at least one nuisance
#'   model is nonparametric / ML (e.g., \code{"RF"}, \code{"gbm"}, \code{"RSF"}, \code{"Spline"}) and may disable it
#'   when all nuisance models are parametric.
#'
#' @param weights Optional numeric vector of observation weights (length \code{nrow(data)}). If \code{NULL}, uses
#'   equal weights.
#'
#' @param shuffle Logical; if \code{TRUE}, shuffles rows before constructing folds for cross-fitting.
#'
#' @param PS_min Numeric lower bound for propensity scores used in weighting (clipping) to stabilize weights.
#'
#' @param S_min Numeric lower bound for censoring survival probabilities used in weighting to avoid extreme IPC weights.
#'
#' @param add_naive_est Logical; if \code{TRUE}, also computes a naive (unadjusted) estimator via \code{get_naive_est()}
#'   for comparison.
#'
#' @param n_boot Integer number of bootstrap replicates for uncertainty quantification. Use \code{0} to skip inference.
#'
#' @param boot_raw Logical; if \code{TRUE}, returns raw bootstrap replicates instead of summarized SE/CI.
#'
#' @param alpha Significance level for confidence intervals (e.g., \code{0.05} yields 95\% CIs).
#'
#' @param seed Optional integer random seed for reproducibility.
#'
#' @param add_IPW Logical; if \code{TRUE}, also computes an IPW-only estimator in addition to AIPCW.
#'
#' @param name Optional character name(s) for method labels. If \code{NULL}, defaults to
#'   \code{"causalAIPCW"} (and \code{"IPW"} if \code{add_IPW = TRUE}).
#'
#' @param CI_raw Logical; if \code{TRUE}, returns confidence intervals on the original (non-transformed) scale when
#'   transformations are used internally.
#'
#' @param silent Logical; if \code{TRUE}, suppresses progress messages.
#'
#' @param add_contrast Logical; if \code{TRUE}, additionally reports contrasts between treatment arms as implemented
#'   internally.
#'
#' @param parallel Logical; if \code{TRUE}, uses parallel bootstrap via \pkg{future} / \pkg{future.apply}.
#'
#' @param n_cores Optional integer number of cores to use when \code{parallel = TRUE}. If \code{NULL},
#'   uses \code{parallel::detectCores()}.
#'
#' @details
#'  When \code{boot_raw = FALSE}, bootstrap replicates are summarized
#' via internal helper functions (\code{get_se()} and \code{get_CI()}).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{est}}{A data frame of point estimates indexed by \code{time}. Contains one or more \code{est_*} columns
#'     (e.g., \code{est_CIF1}, \code{est_CIF2}, \code{est_CIF3}) depending on \code{estimand}. A \code{method} column may
#'     appear when multiple methods are requested.}
#'   \item{\code{est_se}}{If \code{n_boot > 0}: either summarized SE/CI results (default) or raw bootstrap replicates when
#'     \code{boot_raw = TRUE}. If \code{n_boot = 0}, \code{NULL}.}
#'   \item{\code{nuisance}}{A list of fitted nuisance components (e.g., fitted PS,
#'     censoring, and time models).}
#'   \item{\code{runtime}}{Runtime for obtaining the point estimates.}
#'   \item{\code{call}}{The matched call with defaults filled in.}
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate some data
#' dat = sim_data(n=200)
#' # Estimate CIF1, CIF2 and CIF3 at a user-specified grid of times with
#' default nuisance models without SEs
#' fit <- causalAIPCW(
#'   data = dat,
#'   covars = c("Z1", "Z2"),
#'   time_interest = seq(0.5, 5, by = 0.5),
#'   time1_interest = c(0.5, 1),
#'   n_boot = 0,
#'   seed = 1
#' )
#' fit$est
#'
#' # Include terminal risk (reported as CIF3) and add IPW for comparison
#'fit2 <- causalAIPCW(
#'  data = dat,
#'  covars = c("Z1", "Z2"),
#'  estimand = c('CIF1', 'CIF2'),
#'  time_interest = seq(0.5, 5, by = 0.5),
#'  add_IPW = TRUE,
#'  n_boot = 10
#')
#' }
#'
#' @export
causalAIPCW = function(data, X1='X1',
                       X2='X2', delta1='delta1', delta2='delta2',
                       treatment='A', covars='Z',
                       estimand = 'all',
                       tau = NULL,
                       time_interest = NULL,
                       time1_interest = NULL,
                       PS_model = 'logit', censor_model = 'Cox',
                       time_model = 'Cox', c_nodesize = 15,
                       c_ntree = 500,  t_nodesize = 15, t_ntree = 500,
                       k = 5, crossfit = 'auto', weights = NULL,
                       shuffle = T,
                       PS_min = 0.1, S_min = 0.05, add_naive_est = F,
                       n_boot = 100, boot_raw = F, alpha = 0.05,
                       # name: name of the estimating method, eg: i use DR:Cox/CoX_logit to denote the
                       # models for T/C_A.
                       seed = NULL, add_IPW = F, name = NULL,
                       CI_raw = F, silent = T,
                       add_contrast = F, parallel = F, n_cores = NULL
){
  if (!is.null(seed)) set.seed(seed)
  # check input values
  if (!is.data.frame(data)) stop('Input dataset should be a dataframe!')
  if (!PS_model %in% c('logit','logit.quad', 'logit.small', 'RF', 'gbm')) stop('The model for Propensity Score should be one of the following: logit, RF, gbm')
  if (!censor_model %in% c('Cox','Cox.quad', 'Cox.small', 'RSF', 'Spline')) stop('The model for censoring should be one of the following: Cox, RSF, Spline')
  if (!time_model %in% c('Cox','RSF', 'Spline')) stop('The model for event time should be one of the following: Cox, RSF, Spline')
  if (is.null(time1_interest) & ('CIF3' %in% estimand|any(estimand == 'all'))) stop("Please specify time to non-terminal event to condition on!")
  if (any(!estimand %in% c('all', 'CIF1','CIF2','CIF3'))) stop('The estimand should be one or some of the following: all, CIF1, CIF2, CIF3')

  if (is.null(name)) {
    specify_name = F
    if (add_IPW) name = c('causalAIPCW','IPW')
    else name = 'causalAIPCW'}
  else{specify_name = T}

  if(!boot_raw & n_boot>0) add_trans = T
  else add_trans = F

  n = nrow(data)

  if (is.null(weights)) weights = rep(1, n)

  # clean data
  data = data %>%
    dplyr::select(all_of(c(X1, X2, delta1, delta2, treatment, covars))) %>%
    setNames(c('X1', 'X2', 'delta1', 'delta2', 'A', paste0('Z', 1:length(covars))))

  # the max of time_interest needs to be <tau<max(X2)
  if (!is.null(tau)){
    if(tau<max(data$X2)) {
      data = data %>%
        mutate_at(vars(X1, X2), ~ifelse(. >tau, tau, .)) %>%
        mutate(delta1 = ifelse(X1>tau, 0, delta1),
               delta2 = ifelse(X2>tau, 0, delta2))

    }
  }


  data = data %>%
    mutate(delta_c = 1-delta2) %>%
    # status: 0=censored, 1=illness, 2=death without illness
    mutate(status = ifelse(!delta1 & delta2, 2, delta1))

  time_all =  unique(c(data$X1, data$X2)) %>% sort()

  if(!is.null(time1_interest) ){
    time1_interest = sort(time1_interest)}

  fit = causalAIPCW_est(data,
                        n,
                        time_all,
                        estimand,
                        time_interest,
                        time1_interest,
                        PS_model, censor_model,
                        time_model, c_nodesize,
                        c_ntree,  t_nodesize , t_ntree ,
                        k, crossfit, weights ,
                        shuffle,
                        PS_min, S_min, add_IPW, name, add_trans, silent, add_contrast)
  est = fit$est
  if(add_naive_est){
    est_naive  = get_naive_est(data, time_all, estimand, time_interest,
                               time1_interest,weights, add_trans, add_contrast)$est
    est = rbind(est, cbind(method = 'naive', est_naive))
  }

  # n_boot = 0 means no SE will be computed
  est_se = NULL

  if (n_boot > 0) {
    if (parallel) {
      if (is.null(n_cores)) {n_cores = parallel::detectCores()}
      future::plan("multisession", workers = n_cores)
      boot_results = future.apply::future_lapply(seq_len(n_boot), function(i) {
        sourceCpp("get_integral.cpp")
        boot_index = sample(1:n, n, replace = T)
        boot_data = data[boot_index, ]
        boot_weights = NULL
        fit_boot = causalAIPCW_est(
          boot_data,
          n,
          time_all,
          estimand,
          time_interest,
          time1_interest,
          PS_model,
          censor_model,
          time_model,
          c_nodesize,
          c_ntree,
          t_nodesize ,
          t_ntree ,
          k,
          crossfit,
          boot_weights,
          shuffle,
          PS_min,
          S_min,
          add_IPW,
          name,
          add_trans,
          silent,
          add_contrast
        )
        cbind(boot_id = i, fit_boot$est)
      }, future.seed=T) %>% bind_rows()
    }
    else{
      boot_results = vector('list', n_boot)
      for (i in 1:n_boot) {
        boot_index = sample(1:n, n, replace = T)
        boot_data = data[boot_index, ]
        boot_weights = NULL


        fit_boot = causalAIPCW_est(
          data = boot_data,
          n,
          time_all,
          estimand,
          time_interest,
          time1_interest,
          PS_model,
          censor_model,
          time_model,
          c_nodesize,
          c_ntree,
          t_nodesize ,
          t_ntree ,
          k,
          crossfit,
          boot_weights,
          shuffle,
          PS_min,
          S_min,
          add_IPW,
          name,
          add_trans,
          silent,
          add_contrast
        )
        # if(length(fit_boot)==n) {
        #   return(fit_boot)}
        boot_results[[i]] = cbind(boot_id = i, fit_boot$est)
      }
      boot_results = bind_rows(boot_results)
    }

    if(boot_raw){
      est_se = boot_results
    }else{
      se_results = get_se(boot_results)
    }

    #print(se_results)
    if (add_naive_est) {
      boot_results_naive = vector('list', n_boot)
      for (i in 1:n_boot) {
        boot_index = sample(1:n, n, replace = T)
        boot_data = data[boot_index, ]
        boot_weights = NULL
        est_boot_naive = get_naive_est(data = boot_data,
                                       time_all = time_all,
                                       estimand = estimand,
                                       time_interest = time_interest,
                                       time1_interest = time1_interest,
                                       weights = boot_weights,
                                       add_trans,
                                       add_contrast)$est

        boot_results_naive[[i]] = cbind(boot_id = i, est_boot_naive)
      }
      boot_results_naive = bind_rows(boot_results_naive)
      if(boot_raw){
        est_se = rbind(est_se,
                       cbind(method = 'naive', boot_results_naive))
      }else{
        se_results_naive = get_se(boot_results_naive)
        se_results = rbind(
          se_results,
          cbind(method = 'naive', se_results_naive)
        )
      }

    }

    if(!boot_raw){
      est_se = est %>%
        pivot_longer(-c(method, time), names_to = 'estimand', values_to = 'est',
                     names_prefix = 'est_') %>%
        get_CI(se_results, alpha = alpha, CI_raw=CI_raw)
    }

    if(!add_naive_est & !add_IPW & !specify_name) {est_se = select(est_se, -method)}
  }

  if(!add_naive_est & !add_IPW & !specify_name) {est = select(est, -method)}


  full_call = as.list(match.call())
  defaults  = formals(causalAIPCW)
  for (arg in names(defaults)) {
    if (!arg %in% names(full_call)) {
      if(is.null(defaults[[arg]])) {
        # specify a list element to be null
        full_call[arg] = list(NULL)
      }else{
        full_call[[arg]] <- defaults[[arg]]
      }
    }
  }

  # rename estimand names to be consistent with the paper
  est = est %>% select(-contains('_trans')) %>% rename_at(vars(contains('risk')), ~str_replace(., 'risk', 'CIF3'))
  if (!is.null(est_se)) est_se = est_se %>%
    mutate(estimand = str_replace(estimand, 'risk', 'CIF3'))

  list(est = est,
       est_se = est_se,
       nuisance = fit$nuisance,
       runtime = fit$runtime,
       call = as.call(full_call)
  )

}

