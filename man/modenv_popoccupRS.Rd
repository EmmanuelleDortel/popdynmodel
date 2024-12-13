% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modenv_popoccupRS.R
\name{modenv_popoccupRS}
\alias{modenv_popoccupRS}
\title{Population occupancy modelling accounting for detection errors in a Bayesian hierarchical framework with environmental stochasticity}
\usage{
modenv_popoccupRS(
  df,
  var_id,
  var_tmp,
  var_cnt,
  var_pas,
  var_pro = NULL,
  var_tax = NULL,
  var_envO = NULL,
  var_envP = NULL,
  var_envC = NULL,
  var_det = NULL,
  var_reg = NULL,
  var_guild = NULL,
  period = NULL,
  timestep = 1,
  save_parameters = NULL,
  n_chain = 3,
  n_iter = 102000,
  n_thin = 1000,
  n_burnin = 2000
)
}
\arguments{
\item{df}{data frame with variables used in the model}

\item{var_id}{[character] variable name with sampling location identifiers (must not contain missing values)}

\item{var_tmp}{[character] variable name with continuous sampling times (must be numeric without missing values)}

\item{var_cnt}{[character] variable name used to define presence/absence (must be numeric)}

\item{var_pas}{[character] variable name with pass numbers in a multiple-pass removal sampling (must be numeric without missing values)}

\item{var_pro}{[optional, character] variable name with sampling protocol (must not contain missing values)}

\item{var_tax}{[optional, character] variable name with taxa identifiers to be specified for the multi-taxa models (must not contain missing values)}

\item{var_envO}{[optional, character] environmental covariate(s) name(s) related to occupancy probabilities whose values must be supplied by location and/or by time (must be numeric without missing values)}

\item{var_envP}{[optional, character] environmental covariate(s) name(s) related to persistence probabilities whose values must be supplied by location and/or by time (must be numeric without missing values)}

\item{var_envC}{[optional, character] environmental covariate(s) name(s) related to colonisation probabilities whose values must be supplied by location and/or by time (must be numeric without missing values)}

\item{var_det}{[optional, character] environmental covariate(s) name(s) related to taxa detection probabilities whose values must be supplied by location and/or by time (must be numeric without missing values)}

\item{var_reg}{[optional, character] variable name(s) with spatial level identifiers. If supplied, estimates of occupancy changes rates at intermediate spatial levels are added to the model. Otherwise, occupancy change rates are only estimated across all locations}

\item{var_guild}{[optional, character] variable name(s) with taxa guilds identifiers. If supplied, estimates of guild occupancy changes rates are added to the model}

\item{period}{[optional, numeric] list of vectors with start and end times for estimating intermediate average growth rates}

\item{timestep}{[optional, numeric] interval between two consecutive times to account for regular gaps in time series in the rate estimates}

\item{save_parameters}{[optional, character] name(s) of the parameter(s) to save (by default, all model parametes are saved)}

\item{n_chain}{[optional, numeric] number of MCMC chains}

\item{n_iter}{[optional, numeric] number of iterations per MCMC chain}

\item{n_thin}{[optional, numeric] backup interval for MCMC samples}

\item{n_burnin}{[optional, numeric] number of initial iterations to discard per Markov chain}
}
\details{
\code{modenv_popoccupRS} relies on Bayesian hierarchical state-space models that describe the occupancy dynamics of one or more taxa from presence-absence data. The successive occupancy states of locations defined by \code{var_id} depend on taxon-specific responses to environmental covariates. At the global level (i.e. across all locations) and at intermediate spatial levels defined by \code{var_reg}, the taxa occupancy changes rates are estimated from the ratio of the number of occuped location from one time to the next. Their product and their geometric mean are respectively used as estimates of occupancy change rate and average occupancy change rate over the whole time modelling period and over intermediate time periods defined by \code{period}. If \code{var_guild} is supplied, the guild occupancy change rates are also assessed from the geometric mean of taxa occupancy change rates. 

If \code{var_envO} is supplied, the occupancy states of the locations are modelled from the linear effect of environmental covariates on occupancy probabilities.

If \code{var_envP} and \code{var_envC} are supplied, the occupancy states of the locations are modelled from the linear effect of environmental covariates on persistence and colonisation probabilities.

If only \code{var_envP} is supplied, the occupancy states of the locations are modelled from the linear effect of environmental covariates on persistence probabilities and from taxon- and location-specific colonisation probabilities.

If only \code{var_envC} is supplied, the occupancy states of the locations are modelled from the linear effect of environmental covariates on colonisation probabilities and from taxon- andlocation-specific persistence probabilities.

\code{modenv_popoccupRS} accounts for errors in taxa detection. If \code{var_det} is supplied, the taxon-specific detection probabilities are modelled from the linear effect of environmental covariates. Otherwise, the taxon-specific detection probabilities depend on locations and sampling protocols \code{var_pro} is supplied.

The models are run using the \code{runMCMC} function of \code{nimble} package.
}
\value{
\code{modenv_popoccupRS} returns a list with the summaries of main statistics associated with the estimation of the requested parameters (\code{mcmc_summary}), the MCMC chains (\code{mcmc_chain}) and the list of parameters subscripts (\code{subscript}).
}
\seealso{
\code{\link{modenv_popoccup}} and \code{\link{get_modparameters}}
}
\description{
Fits Bayesian hierarchical occupancy models to presence-absence time series for estimating the occurrence of one or more taxa from relationships with environmental covariates. The modelling framework accounts for taxa detection biases using multiple-pass removal sampling method.
}
\examples{
\dontrun{
## Occupancy dynamics with linear environmental effects on occupancy probability
mcmc.out <- modenv_popoccupRS(riverfish,
                              var_id = "pop_id",
                              var_tmp = "year",
                              var_tax = "taxa",
                              var_pas = "pass",
                              var_cnt = "headcount",
                              var_envO = "temperature")

# Add linear environmental effects on detection probability
mcmc.out <- modenv_popoccupRS(riverfish,
                              var_id = "pop_id",
                              var_tmp = "year",
                              var_tax = "taxa",
                              var_pas = "pass",
                              var_cnt = "headcount",
                              var_envO = "temperature",
                              var_det = "depth")

# Add sampling protocol effect on detection probability
df <- mef_convertquali(riverfish, var_quali="fishing")
mcmc.out <- modenv_popoccupRS(df, 
                              var_id = "pop_id",
                              var_tmp = "year",
                              var_tax = "taxa",
                              var_pas = "pass",
                              var_cnt = "headcount",
                              var_envO = "temperature",
                              var_det = c("depth","fishing_complete"))

## Occupancy dynamics of pike with linear environmental effects on detection probability
mcmc.out <- modenv_popoccupRS(riverfish[riverfish$taxa \%in\% "pike",],
                              var_id = "pop_id",
                              var_tmp = "year",
                              var_pas = "pass",
                              var_cnt = "headcount",
                              var_det = c("temperature","depth"))
}
}
