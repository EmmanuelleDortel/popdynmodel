#' Writes model codes and creates input date files
#' 
#' @importFrom rlang check_required enquo syms abort
#'
#' @param df data frame
#' @param model character
#' @param var_id character
#' @param var_tmp character
#' @param var_tax character
#' @param var_cnt character
#' @param var_wei character
#' @param var_surf character
#' @param var_reg character or character vector
#' @param var_guild character or character vector
#' @param var_env character or character vector
#' @param var_envO character or character vector
#' @param var_envP character or character vector
#' @param var_envC character or character vector
#' @param var_pro character
#' @param var_pas character
#' @param var_det character or character vector
#' @param period list
#' @param timestep numeric
#'
#' @return output list
#' @export
#'
#' @examples
#' df <- riverfish[riverfish$pass == 1,] 
#' modpopdyn <- wri_popmodel(df,
#' model = "mod_popdyn",
#' var_id = "pop_id",
#' var_tmp = "year",
#' var_tax = "taxa",
#' var_cnt = "headcount")
#' modpopdynAlt <- wri_popmodel(df,
#' model = "mod_popdynAlt",
#' var_id = "pop_id",
#' var_tmp = "year",
#' var_tax = "taxa",
#' var_cnt = "headcount",
#' var_pro = "fishing")
#' modenvpopgrow <- wri_popmodel(df[df$taxa %in% "eel",],
#' model = "modenv_popgrow",
#' var_id = "pop_id",
#' var_tmp = "year",
#' var_env = "temperature",
#' var_cnt = "headcount")
wri_popmodel <- function(df, model, var_id, var_tmp, var_tax=NULL, var_cnt=NULL, var_wei=NULL, var_surf=NULL, var_reg=NULL, var_guild=NULL, var_env=NULL, var_envO=NULL, var_envP=NULL, var_envC=NULL, var_pro=NULL, var_pas=NULL, var_det=NULL, period=NULL, timestep=1) {
  #-----------------------------------------------------------------------------
  # Check for missing required arguments
  check_required(model)
  check_required(var_id)
  check_required(var_tmp)
  if (quo_is_null(enquo(var_cnt)) & quo_is_null(enquo(var_wei))) {
    abort("'var_cnt' or 'var_wei' must be supplied")
  }
  #-----------------------------------------------------------------------------
  # Select model and check for mistakes
  occup <- grow <- modenv <- modenvG <- alt <- RS <- FALSE
  if (model %in% c("mod_popoccup","modenv_popoccup","mod_popoccupRS","modenv_popoccupRS")) { 
    occup <- TRUE
  } else if (model %in% c("mod_popgrow","modenv_popgrow","mod_popgrowAlt")) {
    grow <- TRUE
  } else { occup <- grow <- TRUE }
  if (model %in% c("mod_popoccupRS","modenv_popoccupRS")) { 
    check_required(var_pas)
    RS <- TRUE
    vars_pas <- syms(c(var_pas, var_id, var_tmp, var_tax))
  } else { vars_pas <- NULL }
  if (model %in% c("mod_popgrowAlt","mod_popdynAlt")) { 
    alt <- TRUE
  }  
  if (model %in% c("modenv_popoccup","modenv_popoccupRS","modenv_popgrow","modenv_popdyn")) {
    df <- do.call(int_checkvarenv, list(df, var_id=enquo(var_id),
                                        var_tmp=enquo(var_tmp), 
                                        vars=syms(c(var_env,var_envO,var_envP,var_envC,var_det))))
  }
  df <- do.call(int_checkfunction, list(df,
                                        vars_in_df=syms(c(var_id, var_tmp, var_tax, var_cnt, var_wei, var_env, var_envO, var_envP, var_envC, var_surf, var_reg, var_guild, var_pas, var_pro, var_det)),
                                        vars_na=syms(c(var_id, var_tmp, var_tax, var_env, var_envO, var_envP, var_envC, var_surf, var_pas, var_pro, var_det)),
                                        vars_numeric=syms(c(var_tmp, var_cnt, var_wei, var_env, var_envO, var_envP, var_envC, var_surf, var_pas, var_det)),
                                        vars_duplicate=syms(c(var_id, var_tmp, var_tax, var_pas)),
                                        var_tmp=enquo(var_tmp), timestep, period,
                                        vars_pas=vars_pas))
  #-----------------------------------------------------------------------------
  var_id <- enquo(var_id)
  var_tmp <- enquo(var_tmp)
  var_tax <- enquo(var_tax)
  var_cnt <- enquo(var_cnt)
  var_wei <- enquo(var_wei)
  var_grow <- enquo(var_env)
  var_envO <- enquo(var_envO)
  var_envP <- enquo(var_envP)
  var_envC <- enquo(var_envC)
  var_surf <- enquo(var_surf)
  var_reg <- enquo(var_reg)
  var_guild <- enquo(var_guild)
  var_pro <- enquo(var_pro)
  var_pas <- enquo(var_pas)
  var_det <- enquo(var_det)
  if (quo_is_null(var_cnt)) { var_pres <- var_wei } else { var_pres <- var_cnt }
  #-----------------------------------------------------------------------------
  # Write model data
  if (!quo_is_null(var_grow)) { modenvG <- TRUE }
  if (!quo_is_null(var_envO) & !quo_is_null(var_envP) &  !quo_is_null(var_envC)) { modenv <- TRUE }
  datamodel <- do.call(int_datamodel, list(df, occup, grow, modenv, modenvG, alt, RS, timestep, period, var_id, var_tmp, var_tax, var_pres, var_reg, var_guild, var_cnt, var_wei, var_surf, var_pro, var_pas, var_envO, var_envP, var_envC, var_grow, var_det))
  pop_data <- datamodel$popdyn_data
  pop_constants <- datamodel$popdyn_const
  pop_inits <- datamodel$popdyn_inits
  pop_parameters <- datamodel$popdyn_parameters
  pop_internal <- list(datamodel$popdyn_int,var_id=var_id, var_tmp=var_tmp, var_tax=var_tax, period=period)
  #-----------------------------------------------------------------------------
  # Write model code
  code <- do.call(int_popoccup, list(occup,modenv,RS,var_envO,var_envP,var_envC,var_guild,var_det))
  if (isFALSE(grow)) {
    pop_code <- as.call(c(as.symbol("{"), code))
  } else if (isFALSE(alt)) {
    pop_code <- do.call(int_popgrow, list(code,modenvG,var_cnt,var_wei,var_guild))
  } else {
    pop_code <- do.call(int_popgrowAlt, list(code,var_cnt,var_wei,var_guild))
  }
  #-----------------------------------------------------------------------------
  output <- list(pop_code=pop_code,pop_data=pop_data, pop_constants=pop_constants,pop_inits=pop_inits,pop_parameters=pop_parameters,pop_internal=pop_internal)
  return(output)
}