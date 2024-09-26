#' Writes summary information from MCMC chains
#' 
#' @importFrom coda is.mcmc.list
#' @importFrom MCMCvis MCMCsummary
#'
#' @param mcmc_chain mcmc or mcmc.list
#' @param pop_internal list
#'
#' @return output list
#' @export
#'
#' @examples
#' \dontrun{
#' df <- riverfish[riverfish$pass == 1,]
#' modpopoccup <- wri_popmodel(df, 
#' model = "mod_popoccup",
#' var_id = "pop_id",
#' var_tmp = "year",
#' var_tax = "taxa",
#' var_cnt = "headcount")
#' nimModel <- nimbleModel(code = modpopoccup$pop_code,
#' constants = modpopoccup$pop_constants,
#' data = modpopoccup$pop_data,
#' inits = modpopoccup$pop_inits,
#' name = "nimModel", calculate = FALSE)
#' MCMCConf <- configureMCMC(nimModel, monitors = modpopoccup$pop_parameters)
#' MCMCbuild <- buildMCMC(MCMCConf)
#' MCMCComp <- compileNimble(nimModel)
#' nimModelComp <- compileNimble(MCMCbuild, project = nimModel, resetFunctions = TRUE)
#' mcmc_chain <- runMCMC(nimModelComp, nchains = 3, niter = 100, thin = 1, nburnin = 10, samplesAsCodaMCMC = TRUE)
#' MCMC_summary <- wri_MCMCsummary(mcmc_chain, modpopoccup$pop_internal)  
#' }
wri_MCMCsummary <- function(mcmc_chain, pop_internal) {
  #-----------------------------------------------------------------------------
  # Remove missing values from MCMC chain
  if (isTRUE(is.mcmc.list(mcmc_chain))) { 
    n_chain <- length(mcmc_chain)
  } else { n_chain <- 1 }
  if (n_chain == 1 & !is.null(dim(mcmc_chain))) {
    mcmc_na <- which(is.na(mcmc_chain[1,]))
    if (length(mcmc_na) > 0) { mcmc_chain <- mcmc_chain[,-mcmc_na] }
  }
  if (n_chain > 1 & !is.null(dim(mcmc_chain[[1]]))) { for (i in 1:n_chain) {
    mcmc_chain[[i]] <- mcmc_chain[[i]][,!colnames(mcmc_chain[[i]]) %in% names(which(is.na(mcmc_chain[[i]][1,])))]
  }}
  #-----------------------------------------------------------------------------
  # Extract summary information from MCMC chain
  mcmc_summary <- MCMCsummary(mcmc_chain, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), Rhat = TRUE, n.eff = TRUE)
  mcmc_summary <- do.call(int_mpaestimate, list(mcmc_summary, mcmc_chain))
  #-----------------------------------------------------------------------------
  # Set summary data frame and list of subscripts
  var_id <- pop_internal$var_id
  var_tmp <- pop_internal$var_tmp
  var_tax <- pop_internal$var_tax
  period <- pop_internal$period
  datamodel$popdyn_int <- pop_internal[[1]]
  list_summary <- do.call(int_transformsummary, list(mcmc_summary, datamodel, var_id, var_tmp, var_tax, period))
  #-----------------------------------------------------------------------------
  # set summary data frame and list of subscripts
  output <- list(mcmc_summary = list_summary$mcmc_summary, mcmc_chain = mcmc_chain, subscript = list_summary$subscript)
  return(output)
}