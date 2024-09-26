#' Runs multi-taxa environmental occupancy models in a Bayesian hierarchical framework
#'
#' @import magrittr
#' @import nimble
#' @importFrom rlang check_required enquo syms abort
#' @importFrom MCMCvis MCMCsummary
#'
#' @param df data frame
#' @param var_id character
#' @param var_tmp character
#' @param var_tax character
#' @param var_cnt character
#' @param var_pas character
#' @param var_pro character
#' @param var_envO character or character vector
#' @param var_envP character or character vector
#' @param var_envC character or character vector
#' @param var_det character or character vector
#' @param var_reg character or character vector
#' @param var_guild character or character vector
#' @param period list
#' @param timestep numeric
#' @param save_parameters character or character vector
#' @param n_chain numeric
#' @param n_iter numeric
#' @param n_thin numeric
#' @param n_burnin numeric
#'
#' @return output list
#' @export
#'
#' @examples
#' \dontrun{
#' mcmc.out <- modenv_popoccupRS(riverfish,
#' var_id = "pop_id",
#' var_tmp = "year",
#' var_tax = "taxa",
#' var_pas = "pass",
#' var_cnt = "headcount",
#' var_envO = "temperature")
#' mcmc.out <- modenv_popoccupRS(riverfish,
#' var_id = "pop_id",
#' var_tmp = "year",
#' var_tax = "taxa",
#' var_pas = "pass",
#' var_cnt = "headcount",
#' var_envO = "temperature",
#' var_det = "depth")
#' df <- mef_convertquali(riverfish, var_quali="fishing")
#' mcmc.out <- modenv_popoccupRS(df, 
#' var_id = "pop_id",
#' var_tmp = "year",
#' var_tax = "taxa",
#' var_pas = "pass",
#' var_cnt = "headcount",
#' var_envO = "temperature",
#' var_det = c("depth","fishing_complete"))
#' mcmc.out <- modenv_popoccupRS(riverfish[riverfish$taxa %in% "pike",],
#' var_id = "pop_id",
#' var_tmp = "year",
#' var_pas = "pass",
#' var_cnt = "headcount",
#' var_det = c("temperature","depth"))
#' }
modenv_popoccupRS <- function(df, var_id, var_tmp, var_cnt, var_pas, var_pro=NULL, var_tax=NULL, var_envO=NULL, var_envP=NULL, var_envC=NULL, var_det=NULL, var_reg=NULL, var_guild=NULL, period=NULL, timestep=1, save_parameters=NULL, n_chain=3, n_iter=102000, n_thin=1000, n_burnin=2000) {
  #-----------------------------------------------------------------------------
  # Check for missing required arguments
  check_required(var_id)
  check_required(var_tmp)
  check_required(var_cnt)
  check_required(var_pas)
  if (quo_is_null(enquo(var_envO)) & quo_is_null(enquo(var_envP)) & quo_is_null(enquo(var_envC)) & quo_is_null(enquo(var_det))) {
    abort("'var_envO', 'var_envP', 'var_envC' or 'var_det' must be supplied")
  }
  #-----------------------------------------------------------------------------
  # Check for mistakes, if failure return an error message and stop
  df <- do.call(int_checkvarenv, list(df, var_id=enquo(var_id), var_tmp=enquo(var_tmp), vars=syms(c(var_envO,var_envP,var_envC,var_det))))
  df <- do.call(int_checkfunction, list(df,
                                        vars_in_df=syms(c(var_id, var_tmp, var_tax, var_cnt, var_envO, var_envP, var_envC, var_reg, var_guild, var_pas, var_pro, var_det)),
                                        vars_na=syms(c(var_id, var_tmp, var_tax, var_envO, var_envP, var_envC, var_pas, var_pro, var_det)),
                                        vars_numeric=syms(c(var_tmp, var_cnt, var_envO, var_envP, var_envC, var_pas, var_det)),
                                        vars_duplicate=syms(c(var_id, var_tmp, var_tax, var_pas)),
                                        var_tmp=enquo(var_tmp), timestep, period,
                                        vars_pas=syms(c(var_pas, var_id, var_tmp, var_tax))))
  #-----------------------------------------------------------------------------
  var_id <- enquo(var_id)
  var_tmp <- enquo(var_tmp)
  var_tax <- enquo(var_tax)
  var_pres <- enquo(var_cnt)
  var_pas <- enquo(var_pas)
  var_pro <- enquo(var_pro)
  var_envO <- enquo(var_envO)
  var_envP <- enquo(var_envP)
  var_envC <- enquo(var_envC)
  var_det <- enquo(var_det)
  var_reg <- enquo(var_reg)
  var_guild <- enquo(var_guild)
  #-----------------------------------------------------------------------------
  # Write model and model data
  if (quo_is_null(enquo(var_envO)) & quo_is_null(enquo(var_envP)) & quo_is_null(enquo(var_envC))) { modenv <- FALSE } else { modenv <- TRUE }
  datamodel <- do.call(int_datamodel, list(df, occup=TRUE, grow=FALSE, modenv, modenvG=FALSE, alt=FALSE, RS=TRUE, timestep, period, var_id, var_tmp, var_tax, var_pres, var_reg, var_guild, var_cnt=NULL, var_wei=NULL, var_surf=NULL, var_pro, var_pas, var_envO, var_envP, var_envC, var_grow=NULL, var_det))
  code <- do.call(int_popoccup, list(occup=TRUE,modenv,RS=TRUE,var_envO,var_envP,var_envC,var_guild,var_det))
  popdyn_code <- as.call(c(as.symbol("{"), code))
  #-----------------------------------------------------------------------------
  # Define requested parameters
  if (is.null(save_parameters)) { save_parameters <- datamodel$popdyn_parameters } else {
    if (FALSE %in% is.element(save_parameters, datamodel$popdyn_parameters)) {
      para_name <- save_parameters[which(is.element(save_parameters, datamodel$popdyn_parameters) %in% FALSE)]
      abort(paste0("Some parameters are not in model: '", para_name, "'", collapse = " "))
    }
  }
  #-----------------------------------------------------------------------------
  # Fit model
  set.seed(123)
  popdyn <- nimbleModel(code = popdyn_code,
                        constants = datamodel$popdyn_const,
                        data = datamodel$popdyn_data,
                        inits = datamodel$popdyn_inits,
                        name = "popdyn", calculate = FALSE)
  popdynConf <- configureMCMC(popdyn, monitors = save_parameters)
  popdynMCMC <- buildMCMC(popdynConf)
  popdynComp <- compileNimble(popdyn)
  popdynModel <- compileNimble(popdynMCMC, project = popdyn, resetFunctions = TRUE)
  mcmc_chain <- runMCMC(popdynModel, nchains = n_chain, niter = n_iter, thin = n_thin, nburnin = n_burnin, setSeed = 123, samplesAsCodaMCMC = TRUE)
  if (n_chain == 1 & !is.null(dim(mcmc_chain))) {
    mcmc_na <- which(is.na(mcmc_chain[1,]))
    if (length(mcmc_na) > 0) { mcmc_chain <- mcmc_chain[,-mcmc_na] }
  }
  if (n_chain > 1 & !is.null(dim(mcmc_chain[[1]]))) { for (i in 1:n_chain) {
    mcmc_chain[[i]] <- mcmc_chain[[i]][,!colnames(mcmc_chain[[i]]) %in% names(which(is.na(mcmc_chain[[i]][1,])))]
  }}
  #-----------------------------------------------------------------------------
  # set summary data frame of taxa and guilds and list of subscripts
  mcmc_summary <- MCMCsummary(mcmc_chain, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), Rhat = TRUE, n.eff = TRUE)
  mcmc_summary <- do.call(int_mpaestimate, list(mcmc_summary, mcmc_chain))
  list_summary <- do.call(int_transformsummary, list(mcmc_summary, datamodel, var_id, var_tmp, var_tax, period))
  output <- list(mcmc_summary = list_summary$mcmc_summary, mcmc_chain = mcmc_chain, subscript = list_summary$subscript)
  #-----------------------------------------------------------------------------
  return(output)
}