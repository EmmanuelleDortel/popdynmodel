#' Estimates length-weight relationships in a Bayesian framework
#' 
#' @import magrittr
#' @import dplyr
#' @import nimble
#' @importFrom rlang check_required enquo quo_name quo_is_null :=
#' @importFrom MCMCvis MCMCsummary
#'
#' @param df data frame
#' @param var_len character
#' @param var_wei character
#' @param var_tax character
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
#' out <- mod_lenwei(indriverfish,
#' var_len = "length",
#' var_wei = "weight",
#' var_tax = "taxa")
#' df <- indriverfish[indriverfish$taxa %in% "eel",]
#' out <- mod_lenwei(df,
#' var_len = "length",
#' var_wei = "weight")
#' }
mod_lenwei <- function(df, var_len, var_wei, var_tax=NULL, n_chain=3, n_iter=102000, n_thin=100, n_burnin=2000) {
  #-----------------------------------------------------------------------------
  # Check for missing required arguments
  check_required(var_len)
  check_required(var_wei)
  #-----------------------------------------------------------------------------
  var_len <- enquo(var_len)
  var_wei <- enquo(var_wei)
  var_tax <- enquo(var_tax)
  #-----------------------------------------------------------------------------
  # Check for mistakes, if failure return an error message and stop
  if (FALSE %in% lapply(df[,c(quo_name(var_wei),quo_name(var_len))], is.numeric)) {
    abort(paste0("'",quo_name(var_len),"' and '",quo_name(var_wei),"' must be numeric", collapse = " "))
  }
  if (TRUE %in% is.na(df[,c(quo_name(var_wei),quo_name(var_len))])) {
    abort(paste0("NA or NaN are not allowed in '", quo_name(var_len),"' and '",quo_name(var_wei), "'", collapse=" "))
  }
  if (!quo_is_null(var_tax)) {
    if (TRUE %in% is.na(pull(df, !!var_tax))) {
      abort(paste0("NA or NaN are not allowed in '", quo_name(var_tax),"'", collapse=" "))
    }
  } 
  #-----------------------------------------------------------------------------
  # Write model data
  if (quo_is_null(var_tax)) {
    df$taxa <- 1
    var_tax <- "taxa"
    var_tax <- enquo(var_tax)
  }
  tax <- sort(unique(pull(df, !!var_tax)))
  ntax <- length(tax)
  m <- mapply(function(i) nrow(filter_at(df, quo_name(var_tax), all_vars(. %in% i))), tax, SIMPLIFY = "vector")
  L <- mapply(function(i) { 
    mat <- vector("numeric",length = max(m))
    dt <- filter_at(df, quo_name(var_tax), all_vars(. %in% i)) %>% pull(!!var_len) 
    mat[1:length(dt)] <- dt
    return(mat) }, tax, SIMPLIFY = "array")
  W <- mapply(function(i) { 
    mat <- vector("numeric",length = max(m))
    dt <- filter_at(df, quo_name(var_tax), all_vars(. %in% i)) %>% pull(!!var_wei) 
    mat[1:length(dt)] <- dt
    return(mat) }, tax, SIMPLIFY = "array")
  mod_constants <- list(ntax=ntax, m=m)
  mod_data <- list(L=L, W=W)
  mod_inits <- list(a = rep(0.05, ntax), b = rep(2, ntax))
  #-----------------------------------------------------------------------------
  # Write model
  mod_code <- nimbleCode({
    for (i in 1:ntax) {
      for (j in 1:m[i]) {
        Log_W[j,i] <- log(a[i]) + b[i] * log(L[j,i])
        W[j,i] ~ dlnorm(Log_W[j,i], tau[i])
      }
      a[i] ~ dlnorm(0, 1)
      b[i] ~ dlnorm(0, 1)
      tau[i] ~ dgamma(0.01, 0.01)
    }
  })
  #-----------------------------------------------------------------------------
  # Fit model
  set.seed(123)
  modlenwei <- nimbleModel(code = mod_code,
                           constants = mod_constants,
                           data = mod_data,
                           inits = mod_inits,
                           name = "modlenwei", calculate = FALSE)
  modConf <- configureMCMC(modlenwei, monitors = c("a","b"))
  modMCMC <- buildMCMC(modConf)
  modComp <- compileNimble(modlenwei)
  modModel <- compileNimble(modMCMC, project = modlenwei, resetFunctions = TRUE)
  mcmc_chain <- runMCMC(modModel, nchains = n_chain, niter = n_iter, thin = n_thin, nburnin = n_burnin, setSeed = 123, samplesAsCodaMCMC = TRUE)
  #-----------------------------------------------------------------------------
  # Set summary data frame
  mcmc_summary <- MCMCsummary(mcmc_chain, probs = c(0.025, 0.25, 0.5, 0.75, 0.975), Rhat = TRUE, n.eff = TRUE)
  mcmc_summary <- do.call(int_mpaestimate, list(mcmc_summary, mcmc_chain))
  mcmc_summary$parameter <- substr(row.names(mcmc_summary), 1, 1)
  if (ntax > 1) {
    mcmc_summary$pos <- as.numeric(substr(row.names(mcmc_summary), 3, 3))
    mcmc_summary$paratax <- tax[mcmc_summary$pos]
    mcmc_summary <- mutate(mcmc_summary, !!quo_name(var_tax) := paratax) %>%
      mutate(!!paste("subscript",quo_name(var_tax),sep="_") := pos) %>%
      select(-pos, -paratax)
  }
  #-----------------------------------------------------------------------------
  output <- list(mcmc_summary = mcmc_summary, mcmc_chain = mcmc_chain)
  return(output)
}