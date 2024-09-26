#' Check for mistakes in modelling function arguments
#'
#' @import magrittr
#' @import dplyr
#' @importFrom rlang abort quo_name caller_env
#'
#' @param df data frame
#' @param vars_in_df list of expressions
#' @param vars_na list of expressions
#' @param vars_numeric list of expressions
#' @param vars_duplicate list of expressions
#' @param var_tmp quosure
#' @param timestep numeric
#' @param period list
#' @param vars_pas list of expressions
#' @param call
#'
#' @return df
#' @keywords internal
#' @noRd
int_checkfunction <- function(df, vars_in_df, vars_na, vars_numeric, vars_duplicate, var_tmp, timestep, period, vars_pas=NULL, call = caller_env()) {
  if (FALSE %in% is.element(unlist(vars_in_df[!vars_in_df %in% "NULL"]),colnames(df))) {
    abort("must use existing variables", call = call)
  }
  df <- df[,which(colnames(df) %in% unlist(vars_in_df))] %>% distinct()
  #-----------------------------------------------------------------------------
  vars <- df[,which(colnames(df) %in% unlist(vars_na[!vars_na %in% NULL]))]
  if (TRUE %in% is.na(vars)) {
    vnames <- paste0(names(which(colSums(is.na(vars)) > 0)), collapse=", ")
    abort(paste0("NA or NaN are not allowed in '", vnames, "'", collapse=" "), call = call)
  }
  #-----------------------------------------------------------------------------
  vars <- lapply(df[,which(colnames(df) %in% unlist(vars_numeric[!vars_numeric %in% NULL]))], is.numeric)
  if (FALSE %in% vars) {
    vnames <- paste0(names(vars)[vars %in% TRUE], collapse=", ")
    abort(paste0("'", vnames, "' must be numeric", collapse = " "), call = call)
  }
  #-----------------------------------------------------------------------------
  if (TRUE %in% duplicated(df[,which(colnames(df) %in% unlist(vars_duplicate))])) {
    abort("must not have duplicate entries with multiple matches", call = call)
  }
  #-----------------------------------------------------------------------------
  if (max(diff(sort(unique(pull(df, !!var_tmp))))) > timestep) {
    abort("missing values are not allowed in 'var_tmp'", call = call)
  }
  #-----------------------------------------------------------------------------
  if (!is.null(period)) {
    if (FALSE %in% is.list(period)) {
      abort("'period' must be a list")
    }
    if (FALSE %in% (lapply(period, length) == 2) | FALSE %in% is.numeric(unlist(period))) {
      abort("'period' must be a numeric vector list of length 2", call = call)
    }
    if (FALSE %in% is.element(unlist(period), pull(df, !!var_tmp))) {
      abort("period values out of range", call = call)
    }
  }
  #-----------------------------------------------------------------------------
  if (!is.null(vars_pas)) {
    pasvars <- mapply(function(i) as.vector(vars_pas[[i]], mode = "character"), 2:length(vars_pas), SIMPLIFY = "vector")
    dpas <- reframe(df, across(unlist(vars_pas[[1]]), ~length(.x) == max(.x), .names = "pas"), .by = !!pasvars)
    if (FALSE %in% dpas$pas) {
      abort("missing passes in 'var_pas'", call = call)
    }
  }
  #-----------------------------------------------------------------------------
  return(df)
}


