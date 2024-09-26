#' Writes population occupancy models
#'
#' @importFrom rlang quo_is_null
#' @import nimble
#'
#' @param occup logical
#' @param modenv logical
#' @param RS logical
#' @param var_envO quosure
#' @param var_envP quosure
#' @param var_envC quosure
#' @param var_guild quosure
#' @param var_det quosure
#'
#' @return code
#' @keywords internal
#' @noRd
int_popoccup <- function(occup,modenv,RS,var_envO,var_envP,var_envC,var_guild,var_det) {
  #-----------------------------------------------------------------------------
  # Occupancy model without environmental effects
  if (isFALSE(modenv)) {
    occupancy <- nimbleCode(
      for(s in 1:ntax) {
        for(i in 1:nid[s]) {
          for (t in 2:ntime) {
            z[s,idtax[s,i],t] ~ dbern(p[s,idtax[s,i],t])
            p[s,idtax[s,i],t] <- z[s,idtax[s,i],t-1] * p.per_id[s,idtax[s,i]] + (1 - z[s,idtax[s,i],t-1]) * p.col_id[s,idtax[s,i]]
          }
          z[s,idtax[s,i],1] ~ dbern(p.per_id[s,idtax[s,i]])
          tauper[s,idtax[s,i]] ~ dgamma(0.1, 0.1)
          taucol[s,idtax[s,i]] ~ dgamma(0.1, 0.1)
          epsilon_per[s,idtax[s,i]] ~ dnorm(0, sd = tauper[s,idtax[s,i]])
          epsilon_col[s,idtax[s,i]] ~ dnorm(0, sd = taucol[s,idtax[s,i]])
          logit(p.per_id[s,idtax[s,i]]) <- p.per[s] + epsilon_per[s,idtax[s,i]]
          logit(p.col_id[s,idtax[s,i]]) <- p.col[s] + epsilon_col[s,idtax[s,i]]
          p.ext_id[s,idtax[s,i]] <- 1 - p.per_id[s,idtax[s,i]]
        }
        p.per[s] ~ dunif(0,1)
        p.col[s] ~ dunif(0,1)
        p.ext[s] <- 1 - p.per[s]
      }
    )
  } else if (FALSE %in% quo_is_null(var_envP) & FALSE %in% quo_is_null(var_envC)) {
    #---------------------------------------------------------------------------
    # Occupancy model with linear environmental effects on persistence and colonisation probabilities
    occupancy <- nimbleCode(
      for(s in 1:ntax) {
        for(i in 1:nid[s]) {
          for (t in 2:ntime) {
            z[s,idtax[s,i],t] ~ dbern(p[s,idtax[s,i],t])
            p[s,idtax[s,i],t] <- z[s,idtax[s,i],t-1] * p.per_id[s,idtax[s,i],t] + (1 - z[s,idtax[s,i],t-1]) * p.col_id[s,idtax[s,i],t]
            if (nvar > 1) {
              logit(p.per_id[s,idtax[s,i],t]) <- alpha_per[s] + sum(beta_per[s,1:nvar] * var[idtax[s,i],t,1:nvar])
            } else {
              logit(p.per_id[s,idtax[s,i],t]) <- alpha_per[s] + beta_per[s,1] * var[idtax[s,i],t,1]
            }
            if (ncol > 1) {
              logit(p.col_id[s,idtax[s,i],t]) <- alpha_col[s] + sum(beta_col[s,1:ncol] * col[idtax[s,i],t,1:ncol])
            } else {
              logit(p.col_id[s,idtax[s,i],t]) <- alpha_col[s] + beta_col[s,1] * col[idtax[s,i],t,1]
            }
          }
          z[s,idtax[s,i],1] ~ dbern(0.5)
        }
        taualpha_col[s] ~ dgamma(0.1, 0.1)
        alpha_col[s] ~ dnorm(0, sd = taualpha_col[s])
        taualpha_per[s] ~ dgamma(0.1, 0.1)
        alpha_per[s] ~ dnorm(0, sd = taualpha_per[s])
        for (k in 1:nvar) {
          taubeta_per[s,k] ~ dgamma(0.1, 0.1)
          beta_per[s,k] ~ dnorm(0, sd = taubeta_per[s,k])
        }
        for (k in 1:ncol) {
          taubeta_col[s,k] ~ dgamma(0.1, 0.1)
          beta_col[s,k] ~ dnorm(0, sd = taubeta_col[s,k])
        }
      }
    )
  } else if (FALSE %in% quo_is_null(var_envP)) {
    #---------------------------------------------------------------------------
    # Occupancy model with linear environmental effects on persistence probability
    occupancy <- nimbleCode(
      for(s in 1:ntax) {
        for(i in 1:nid[s]) {
          for (t in 2:ntime) {
            z[s,idtax[s,i],t] ~ dbern(p[s,idtax[s,i],t])
            p[s,idtax[s,i],t] <- z[s,idtax[s,i],t-1] * p.per_id[s,idtax[s,i],t] + (1 - z[s,idtax[s,i],t-1]) * p.col_id[s,idtax[s,i]]
            if (nvar > 1) {
              logit(p.per_id[s,idtax[s,i],t]) <- alpha_per[s] + sum(beta_per[s,1:nvar] * var[idtax[s,i],t,1:nvar])
            } else {
              logit(p.per_id[s,idtax[s,i],t]) <- alpha_per[s] + beta_per[s,1] * var[idtax[s,i],t,1]
            }
          }
          z[s,idtax[s,i],1] ~ dbern(0.5)
          taucol[s,idtax[s,i]] ~ dgamma(0.1, 0.1)
          epsilon_col[s,idtax[s,i]] ~ dnorm(0, sd = taucol[s,idtax[s,i]])
          logit(p.col_id[s,idtax[s,i]]) <- p.col[s] + epsilon_col[s,idtax[s,i]]
        }
        p.col[s] ~ dunif(0,1)
        taualpha[s] ~ dgamma(0.1, 0.1)
        alpha_per[s] ~ dnorm(0, sd = taualpha[s])
        for (k in 1:nvar) {
          taubeta[s,k] ~ dgamma(0.1, 0.1)
          beta_per[s,k] ~ dnorm(0, sd = taubeta[s,k])
        }
      }
    )
  } else if (FALSE %in% quo_is_null(var_envC)) {
    #---------------------------------------------------------------------------
    # Occupancy model with linear environmental effects on colonisation probability
    occupancy <- nimbleCode(
      for(s in 1:ntax) {
        for(i in 1:nid[s]) {
          for (t in 2:ntime) {
            z[s,idtax[s,i],t] ~ dbern(p[s,idtax[s,i],t])
            p[s,idtax[s,i],t] <- z[s,idtax[s,i],t-1] * p.per_id[s,idtax[s,i]] + (1 - z[s,idtax[s,i],t-1]) * p.col_id[s,idtax[s,i],t]
            if (nvar > 1) {
              logit(p.col_id[s,idtax[s,i],t]) <- alpha_col[s] + sum(beta_col[s,1:nvar] * var[idtax[s,i],t,1:nvar])
            } else {
              logit(p.col_id[s,idtax[s,i],t]) <- alpha_col[s] + beta_col[s,1] * var[idtax[s,i],t,1]
            }
          }
          z[s,idtax[s,i],1] ~ dbern(0.5)
          tauper[s,idtax[s,i]] ~ dgamma(0.1, 0.1)
          epsilon_per[s,idtax[s,i]] ~ dnorm(0, sd = tauper[s,idtax[s,i]])
          logit(p.per_id[s,idtax[s,i]]) <- p.per[s] + epsilon_per[s,idtax[s,i]]
        }
        p.per[s] ~ dunif(0,1)
        taualpha[s] ~ dgamma(0.1, 0.1)
        alpha_col[s] ~ dnorm(0, sd = taualpha[s])
        for (k in 1:nvar) {
          taubeta[s,k] ~ dgamma(0.1, 0.1)
          beta_col[s,k] ~ dnorm(0, sd = taubeta[s,k])
        }
      }
    )
  } else {
    #---------------------------------------------------------------------------
    # Occupancy model with linear environmental effects on occupancy probability
    occupancy <- nimbleCode(
      for(s in 1:ntax) {
        for(i in 1:nid[s]) {
          for (t in 1:ntime) {
            z[s,idtax[s,i],t] ~ dbern(p[s,idtax[s,i],t])
            if (nvar > 1) {
              logit(p[s,idtax[s,i],t]) <- alpha[s] + sum(beta[s,1:nvar] * var[idtax[s,i],t,1:nvar])
            } else {
              logit(p[s,idtax[s,i],t]) <- alpha[s] + beta[s,1] * var[idtax[s,i],t,1]
            }
          }
        }
        taualpha[s] ~ dgamma(0.1, 0.1)
        alpha[s] ~ dnorm(0, sd = taualpha[s])
        for (k in 1:nvar) {
          taubeta[s,k] ~ dgamma(0.1, 0.1)
          beta[s,k] ~ dnorm(0, sd = taubeta[s,k])
        }
      }
    )
  }
  #-----------------------------------------------------------------------------
  code <- list(occupancy)
  #-----------------------------------------------------------------------------
  # Occupancy rates at regional levels
  if (isTRUE(occup)) {
    #---------------------------------------------------------------------------
    # Occupancy rates at regional levels (taxa)
    occupancy_regional <- nimbleCode(
      for(s in 1:ntax) {
        for (j in 1:nreg[s]) {
          for (t in 1:nperiod) {
            nz_reg[s,reg[s,j],t] <- max(z_reg[s,reg[s,j],start[t]:end[t]])
            z_lambda[s,reg[s,j],t] <- nz_reg[s,reg[s,j],t] * prod(cal.z_lambda[s,reg[s,j],start[t]:end[t]])
            OR[s,reg[s,j],t] <- nz_reg[s,reg[s,j],t] * 100 * (z_lambda[s,reg[s,j],t] - 1)
            z_mulambda[s,reg[s,j],t] <- z_lambda[s,reg[s,j],t]^(1 / (step * ndate[t]))
            muOR[s,reg[s,j],t] <- nz_reg[s,reg[s,j],t] * 100 * (z_mulambda[s,reg[s,j],t] - 1)
          }
          for (t in 2:ntime) {
            z_reg[s,reg[s,j],t] <- max(z_cal[s,reg[s,j],1:nidreg[s,j],t])
            nb_z[s,reg[s,j],t] <- sum(z_cal[s,reg[s,j],1:nidreg[s,j],t])
            z_intlambda[s,reg[s,j],t] <- nb_z[s,reg[s,j],t] / max(1,nb_z[s,reg[s,j],t-1])
            cal.z_lambda[s,reg[s,j],t] <- z_reg[s,reg[s,j],t] * z_intlambda[s,reg[s,j],t] + (1 - z_reg[s,reg[s,j],t])
            turnover[s,reg[s,j],t] <- sum(turn_cal[s,reg[s,j],1:nidreg[s,j],t]) / max(1,nb_z[s,reg[s,j],t])
            for (i in 1:nidreg[s,j]) {
              turn_cal[s,reg[s,j],i,t] <- z[s,idreg[j,i,s],t] * (1 - z[s,idreg[j,i,s],t-1])
              z_cal[s,reg[s,j],i,t] <- z[s,idreg[j,i,s],t]
            }
          }
          nb_z[s,reg[s,j],1] <- sum(z_cal[s,reg[s,j],1:nidreg[s,j],1])
          for (i in 1:nidreg[s,j]) {
            z_cal[s,reg[s,j],i,1] <- z[s,idreg[j,i,s],1]
          }
        }
      }
    )
    code <- c(code,list(occupancy_regional))
    #---------------------------------------------------------------------------
    # Occupancy rates for guilds
    if (FALSE %in% quo_is_null(var_guild)) {
      occupancy_guild <- nimbleCode(
        for(g in 1:ngui) {
          for (j in 1:nregui[g]) {
            for (t in 1:nperiod) {
              z_guireg[g,regui[g,j],t] <- max(z_gui[g,regui[g,j],start[t]:end[t]])
              z_lambda_gui[g,regui[g,j],t] <- z_guireg[g,regui[g,j],t] * prod(cal.z_lambda_gui[g,regui[g,j],start[t]:end[t]])
              GOR[g,regui[g,j],t] <- z_guireg[g,regui[g,j],t] * 100 * (z_lambda_gui[g,regui[g,j],t] - 1)
              z_mulambda_gui[g,regui[g,j],t] <- z_lambda_gui[g,regui[g,j],t]^(1 / (step * ndate[t]))
              muGOR[g,regui[g,j],t] <- z_guireg[g,regui[g,j],t] * 100 * (z_mulambda_gui[g,regui[g,j],t] - 1)
            }
            for (t in 2:ntime) {
              z_intlambda_gui[g,regui[g,j],t] <- prod(z_lambda_tax[g,1:ngtax[g,j],regui[g,j],t])^(1/ngtax[g,j])
              z_gui[g,regui[g,j],t] <- max(z_reg_tax[g,1:ngtax[g,j],regui[g,j],t])
              cal.z_lambda_gui[g,regui[g,j],t] <- z_gui[g,regui[g,j],t] * z_intlambda_gui[g,regui[g,j],t] + (1 - z_gui[g,regui[g,j],t])

              for(s in 1:ngtax[g,j]) {
                z_lambda_tax[g,s,regui[g,j],t] <- cal.z_lambda[gtax[j,s,g],regui[g,j],t]
                z_reg_tax[g,s,regui[g,j],t] <- z_reg[gtax[j,s,g],regui[g,j],t]
              }
            }
          }
        }
      )
      code <- c(code,list(occupancy_guild))
    }
  }
  #-----------------------------------------------------------------------------
  if (isTRUE(RS)) {
    #---------------------------------------------------------------------------
    # Occupancy rates for guilds
    if (FALSE %in% quo_is_null(var_det)) {
      errorRS <- nimbleCode(
        for(s in 1:ntax) {
          for(i in 1:nid[s]) {
            for (t in 1:ntime) {
              if (ndet > 1) {
                logit(p.det[s,idtax[s,i],t]) <- alpha_det[s] + sum(beta_det[s,1:ndet] * vardet[idtax[s,i],t,1:ndet])
              } else {
                logit(p.det[s,idtax[s,i],t]) <- alpha_det[s] + beta_det[s,1] * vardet[idtax[s,i],t,1]
              }
              for (n in 1:npas[s,idtax[s,i],t]) {
                x[s,idtax[s,i],t,n] ~ dbern(z[s,idtax[s,i],t] * p.det[s,idtax[s,i],t])
              }
            }
          }
          taualpha_det[s] ~ dgamma(0.1, 0.1)
          alpha_det[s] ~ dnorm(0, sd = taualpha_det[s])
          for (k in 1:ndet) {
            taubeta_det[s,k] ~ dgamma(0.1, 0.1)
            beta_det[s,k] ~ dnorm(0, sd = taubeta_det[s,k])
          }
        }
      )
      code <- c(code, list(errorRS))
      #-------------------------------------------------------------------------
    } else {
      errorRS1 <- nimbleCode(
        for(s in 1:ntax) {
          for(i in 1:nid[s]) {
            for (k in 1:npro) {
              logit(p.det_id[s,k,idtax[s,i]]) <- p.det[s,k] + epsilon_det[idtax[s,i]]
            }
            for (t in 1:ntime) {
              for (n in 1:npas[s,idtax[s,i],t]) {
                x[s,idtax[s,i],t,n] ~ dbern(z[s,idtax[s,i],t] * p.det_id[s,pro[idtax[s,i],t],idtax[s,i]])
              }
            }
          }
          for (k in 1:npro) {
            p.det[s,k] ~ dunif(0,1)
          }
        }
      )
      errorRS2 <- nimbleCode(
        for (i in 1:nidtot) {
          taudet[i] ~ dgamma(0.1, 0.1)
          epsilon_det[i] ~ dnorm(0, sd = taudet[i])
        }
      )  
      code <- c(code, list(errorRS1, errorRS2))
    }
  }
  #-----------------------------------------------------------------------------
  return(code)
}
