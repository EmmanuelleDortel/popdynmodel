#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @import magrittr
#' @import dplyr
#' @import nimble
#' @importFrom rlang abort caller_env check_required enquo quo_is_null quo_name syms set_names
#' @importFrom stats na.omit setNames
#' @importFrom stringr str_locate str_split_fixed
#' @importFrom tidyr unite
#' @importFrom tidyselect all_of
#' @importFrom MCMCvis MCMCsummary
globalVariables(c("ntax", "nreg", "nperiod", "reg", "start", "end", "step", "ntime", "nidreg","z","idreg","S","nid","idtax","y","C","ngvar","alpha_N","beta_N","gvar","ngui","nregui","regui","ngtax","gtax","W","w","alpha_B","beta_B","ndate","p.per","epsilon_per","p.col","epsilon_col","nvar","alpha_per","beta_per","var","alpha_col","beta_col","alpha","pro","muN","muB"))
globalVariables(c("p2","p3","p1","parameter","covariate","detection_covariate","paratax","paraid","paratime","para","pos","."))
globalVariables(c("Function","Para","tmp"))
## usethis namespace: end
NULL
