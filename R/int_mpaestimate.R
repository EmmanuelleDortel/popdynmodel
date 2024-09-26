#' Estimates the maximum a posteriori probabilities
#' 
#' @import magrittr
#' @import dplyr
#' @import coda
#' @importFrom stats density
#'
#' @param mcmc_summary data.frame
#' @param mcmc_chain mcmc.list
#'
#' @return data.frame
#' @keywords internal
#' @noRd
int_mpaestimate <- function(mcmc_summary, mcmc_chain) {
  mcmc <- as.matrix(mcmc_chain)
  mpa <- do.call("rbind", lapply(row.names(mcmc_summary), function(i) {
    mcmc.val <- unique(round(mcmc[,i],10))
    if (length(mcmc.val) > 1) {
      mod <- density(mcmc[,i])$x[which(density(mcmc[,i])$y == max(density(mcmc[,i])$y))]
      if (length(mod) > 1) { mod <- NA }
    } else {
      mod <- mcmc.val 
    }
    mat.mod <- matrix(mod, dimnames =  list(i,"mode"))
    return(mat.mod)
  }))
  #-----------------------------------------------------------------------------
  mcmc_summary <- merge(mcmc_summary,mpa,by="row.names",all.x=TRUE) %>%
    relocate(mode, .before = mean) %>%
    set_rownames(.$Row.names) %>% select(-Row.names)
  #-----------------------------------------------------------------------------
  return(mcmc_summary)
}