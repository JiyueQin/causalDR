
# remove the parent environment which could be large in memory
strip_env <- function(obj) {
  if (is.function(obj)) {
    environment(obj) <- baseenv()
    return(obj)
  } else if (is.list(obj)) {
    return(lapply(obj, strip_env))
  } else {
    return(obj)
  }
}

inv.logit = function(x){
  exp(x)/(1 + exp(x))
}


list_obj_size = function(env = .GlobalEnv, cut = 1024*100){
  obj_names = ls(envir = env)

  out_df = map_df(obj_names, ~tibble(
    name = .x,
    size_raw = object.size(get(.x)),
    size = format(size_raw, units = 'auto'),
  )) %>%
    filter(size_raw>cut) %>%
    arrange(desc(size_raw)) %>%
    select(-size_raw)

  out_df
}



get_error = function(vec1, vec2, type = 'mse', return = F){
  if (type== 'mse') out = mean((vec1-vec2)^2)
  else if (type == 'misclass') out = mean(vec1!=vec2)
  cat(type, ':', out, '\n')
  if (return) out
}




get_par_from_str = function(str){
  par_vec = str_split(str_remove_all(str, '[ \n]'), ',')[[1]]
  cat(paste0(par_vec, collapse = '\n'))
}

get_event_type = function(data){
  data %>% mutate(
    event_type = case_when(
      (1-delta1) & (1-delta2) ~ '0.censored before illness or death',
      delta1 & (1-delta2) ~ '1.illness then censor',
      (1-delta1) & delta2 ~ '2.death without illness',
      delta1 & delta2 ~ '3.illness then death'
    ))
}

get_vec_from_mat = function(in_matrix, new_time_index){
  in_matrix[cbind(1:nrow(in_matrix),new_time_index)]
}


# in_matrix: n subjects * m time points
get_jump_size = function(in_matrix, type='Lambda'){
  out =  t(t(in_matrix) - lag(t(in_matrix)))
  if (type == 'S') {
    # for estimates of survival functions
    out[,1] = in_matrix[,1] -1
  }else{
    out[,1] = in_matrix[,1]
  }
  out
}

####################################################################################
# Function fill_time_matrix - fill the values for out_time by repeating the previous value
# return: a matrix with length(out_time) columns
#
####################################################################################


fill_time_matrix = function(in_matrix, in_time, out_time, type = 'S') {
  if (length(setdiff(in_time, out_time))) {
    stop('in_time must be a subset of out_time')
  }
  boolean = out_time %in% in_time
  out_matrix = matrix(NA, nrow(in_matrix), length(out_time))
  out_matrix[, boolean] = in_matrix
  if (type == 'S') {
    out_matrix[, 1] = ifelse(out_matrix[, 1] == 0 |
                               is.na(out_matrix[, 1]), 1, out_matrix[, 1])
  }
  else if (type == 'Lambda') {
    out_matrix[, 1] = ifelse(is.na(out_matrix[, 1]), 0, out_matrix[, 1])
  }
  t(fill(data.frame(t(out_matrix)), everything()))
}
