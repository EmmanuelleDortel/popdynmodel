#' Creates datasets in model format for multi-species models
#'
#' @import magrittr
#' @import dplyr
#' @importFrom rlang quo_name quo_is_null
#'
#' @param df data frame
#' @param occup logical
#' @param grow logical
#' @param modenv logical
#' @param modenvG logical
#' @param alt logical
#' @param RS logical
#' @param timestep numeric
#' @param period list
#' @param var_id quosure
#' @param var_tmp quosure
#' @param var_tax quosure
#' @param var_pres quosure
#' @param var_reg quosure
#' @param var_guild quosure
#' @param var_cnt quosure
#' @param var_wei quosure
#' @param var_surf quosure
#' @param var_pro quosure
#' @param var_pas quosure
#' @param var_envO quosure
#' @param var_envP quosure
#' @param var_envC quosure
#' @param var_grow quosure
#' @param var_det quosure
#'
#' @return list
#' @keywords internal
#' @noRd
int_datamodel <- function(df, occup, grow, modenv, modenvG, alt, RS, timestep, period, var_id, var_tmp, var_tax, var_pres, var_reg, var_guild, var_cnt, var_wei, var_surf, var_pro, var_pas, var_envO, var_envP, var_envC, var_grow, var_det) {
  guild <- varenv <- varcol <- vargrow <- protocol <- vardet <- NULL
  #-----------------------------------------------------------------------------
  if (quo_is_null(var_tax)) {
    df$taxa <- 1
    var_tax <- "taxa"
    var_tax <- enquo(var_tax)
  }
  tax <- sort(unique(pull(df, !!var_tax)))
  ntax <- length(tax)
  #-----------------------------------------------------------------------------
  time <- seq(min(pull(df, !!var_tmp)),max(pull(df, !!var_tmp)),timestep)
  ntime <- length(time)
  #-----------------------------------------------------------------------------
  id <- sort(unique(pull(df, !!var_id)))
  nid <-  mapply(function(i) length(unique(pull(df[df[,quo_name(var_tax)] %in% i,], !!var_id))), tax, SIMPLIFY = "vector")
  idtax <- do.call("rbind",lapply(1:ntax, function(i) c(which(id %in% sort(unique(pull(df[df[,quo_name(var_tax)] %in% tax[i],], !!var_id)))), rep(NA, max(nid) - nid[i]))))
  #-----------------------------------------------------------------------------
  df_reg <- select(df, !!var_tax, !!var_id) %>% distinct() %>% mutate(reg = "global")
  if (FALSE %in% quo_is_null(var_reg)) {
    nr <- which(colnames(df) %in% colnames(select(df, !!var_reg)))
    df_reg <- df_reg %>%
      rbind(do.call("rbind",lapply(nr, function(i)  data.frame(filter(df, !is.na(df[,i])) %>%
                                                                 select(!!var_tax, !!var_id, all_of(i)) %>% distinct() %>%
                                                                 rename(reg = colnames(df[i])) %>%
                                                                 rowwise() %>% mutate(reg = paste(colnames(df[i]),reg,sep="_"))))))
  }
  region <- sort(unique(df_reg$reg))
  # regions per taxa
  nreg <- tapply(df_reg$reg, list(df_reg[,quo_name(var_tax)]), function(i) length(unique(i)))
  reg <- do.call("rbind",lapply(1:ntax, function(i) c(which(region %in% sort(unique(pull(df_reg[df_reg[,quo_name(var_tax)] %in% tax[i],], reg)))), rep(NA, max(nreg) - nreg[i]))))
  # id per taxa and regions
  nidreg <- do.call("rbind", lapply(1:ntax, function(i) c(do.call("cbind",lapply(1:nreg[i], function(j)
    length(unique(pull(df_reg[(df_reg[,quo_name(var_tax)] %in% tax[i] & df_reg$reg %in% region[reg[i,j]]),], !!var_id)))
  )), rep(NA, max(nreg) - nreg[i]))))
  idreg <- sapply(1:ntax, function(i) rbind(do.call("rbind",lapply(1:nreg[i], function(j)
    c(which(id %in% sort(unique(pull(df_reg[(df_reg[,quo_name(var_tax)] %in% tax[i] & df_reg$reg %in% region[reg[i,j]]),], !!var_id)))), rep(NA, max(nidreg, na.rm = T) - nidreg[i,j]))
  )), array(NA, dim=c(max(nreg) - nreg[i], max(nidreg, na.rm = T)))), simplify = "array")
  #-----------------------------------------------------------------------------
  popdyn_const <- list(step = timestep, ntax = ntax, nid = nid, ntime = ntime, idtax = idtax, nreg = nreg, reg = reg, nidreg = nidreg, idreg = idreg)
  #-----------------------------------------------------------------------------
  if (is.null(period)) {
    start <- c(2,NA)
    end <- c(ntime,NA)
    ndate <- c(ntime - 1,NA)
    nperiod <- 1
  } else {
    addperiod <- c(list(c(min(time),max(time))),period)
    nperiod <- length(addperiod)
    start <- mapply(function(i) which(time %in% addperiod[[i]][1]) + timestep, 1:nperiod, SIMPLIFY = "vector")
    end <- mapply(function(i) which(time %in% addperiod[[i]][2]), 1:nperiod, SIMPLIFY = "vector")
    ndate <- mapply(function(i) end[i] - start[i], 1:nperiod,  SIMPLIFY = "vector")
  }
  popdyn_const <- c(popdyn_const, list(nperiod = nperiod, start = start, end = end))
  #-----------------------------------------------------------------------------
  if (FALSE %in% quo_is_null(var_guild)) {
    ng <- which(colnames(df) %in% colnames(select(df, !!var_guild)))
    df_guild <- rbind(do.call("rbind",lapply(ng, function(i) data.frame(filter(df, !is.na(df[,i])) %>%
                                                                          select(!!var_tax, !!var_id, all_of(i)) %>% distinct() %>%
                                                                          rename(guild = colnames(df[i])) %>%
                                                                          rowwise() %>% mutate(guild = paste(colnames(df[i]),guild,sep="_")))))) %>% merge(df_reg)
    # number of guilds
    guild <- sort(unique(df_guild$guild))
    ngui <- length(guild)
    # regions per guilds
    nregui <- tapply(df_guild$reg, list(df_guild$guild), function(i) length(unique(i)))
    regui <- do.call("rbind",lapply(1:ngui, function(i) c(which(region %in% sort(unique(df_guild$reg[df_guild$guild %in% guild[i]]))), rep(NA,max(nregui)-nregui[i]))))
    # taxa per guilds and region
    ngtax <- do.call("rbind",lapply(1:ngui, function(i) c(do.call("cbind",lapply(1:nregui[i], function(j)
      length(unique(pull(df_guild[(df_guild$guild %in% guild[i] & df_guild$reg %in% region[regui[i,j]]),], !!var_tax)))
    )), rep(NA, max(nregui) - nregui[i]))))
    gtax <- sapply(1:ngui, function(i) rbind(do.call("rbind",lapply(1:nregui[i], function(j)
      c(which(tax %in% sort(unique(pull(df_guild[(df_guild$guild %in% guild[i] & df_guild$reg %in% region[regui[i,j]]),], !!var_tax)))), rep(NA, max(ngtax, na.rm = T) - ngtax[i,j]))
    )), array(NA, dim=c(max(nregui) - nregui[i], max(ngtax, na.rm = T)))), simplify = "array")
    if (is.null(dim(gtax))) {
      gtax <- sapply(1:ngui, function(i) rbind(do.call("rbind",lapply(1:nregui[i], function(j)
        c(which(tax %in% sort(unique(pull(df_guild[(df_guild$guild %in% guild[i] & df_guild$reg %in% region[regui[i,j]]),], !!var_tax)))), rep(NA, max(ngtax, na.rm = T) - ngtax[i,j]))
      )), array(NA, dim=c(1, max(ngtax, na.rm = T)))), simplify = "array")
    }
    # add constants and parameters for guild
    popdyn_const <- c(popdyn_const, list(ngui = ngui, nregui = nregui, regui = regui, ngtax = ngtax, gtax = gtax))
  }
  #-----------------------------------------------------------------------------
  if (isTRUE(RS)) {
    npas <- tapply(df[,quo_name(var_pas)], list(df[,quo_name(var_tax)],df[,quo_name(var_id)],df[,quo_name(var_tmp)]), length)
    popdyn_const <- c(popdyn_const, list(npas = npas))
    if (TRUE %in% quo_is_null(var_det)) {
      pro <- array(1, dim = c(length(id),ntime))
      npro <- 1
      if (FALSE %in% quo_is_null(var_pro)) {
        protocol <-  sort(unique(pull(df, !!var_pro)))
        npro <- length(protocol)
        pro <- tapply(df[,quo_name(var_pro)], list(df[,quo_name(var_id)],df[,quo_name(var_tmp)]), function(i) unique(which(protocol %in% i)))
      }
      popdyn_const <- c(popdyn_const, list(npro = npro, pro = pro, nidtot = length(id)))
    }
  }
  #-----------------------------------------------------------------------------
  popdyn_parameters <- c("z","p.per","p.col","p.per_id","p.col_id")
  popdyn_inits <- list(p.per = rep(0.5, ntax), p.col = rep(0.5, ntax), epsilon_per = matrix(0, ntax, length(id)), epsilon_col = matrix(0, ntax, length(id)))
  if (isFALSE(RS)) {
    z <- tapply(df[,quo_name(var_pres)], list(df[,quo_name(var_tax)],df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique)
    z <- replace(z, which(z > 0), 1)
    popdyn_data <- list(z = z)
  } else {
    x <- tapply(df[,quo_name(var_pres)], list(df[,quo_name(var_tax)],df[,quo_name(var_id)],df[,quo_name(var_tmp)],df[,quo_name(var_pas)]), unique)
    x <- replace(x, which(x > 0), 1)
    z <- tapply(replace(df[,quo_name(var_pres)], which(is.na(df[,quo_name(var_pres)])),0), list(df[,quo_name(var_tax)],df[,quo_name(var_id)],df[,quo_name(var_tmp)]), max)
    z <- replace(z, which(z > 0), 1)
    z <- replace(z, which(z == 0), NA)
    popdyn_data <- list(z = z, x = x)
    if (FALSE %in% quo_is_null(var_det)) {
      vardet <- colnames(select(df, !!var_det))
      ndet <- length(vardet)
      nrdet <- which(colnames(df) %in% vardet)
      detvar <- sapply(nrdet, function(i) tapply(df[,i], list(df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique),simplify = "array")
      popdyn_parameters <- c(popdyn_parameters, "alpha_det", "beta_det")
      popdyn_data <- c(popdyn_data,list(vardet = detvar))
      popdyn_const <- c(popdyn_const, list(ndet = ndet))
      popdyn_inits <- c(popdyn_inits, list(alpha_det = rep(0, ntax), beta_det = matrix(0, ntax, ndet)))
    } else {
      popdyn_parameters <- c(popdyn_parameters, "p.det", "p.det_id")
      popdyn_inits <- c(popdyn_inits, list(p.det = matrix(1, ntax, npro), epsilon_det = rep(0, length(id)))) 
    }
  }
  if (isTRUE(occup)) {
    popdyn_parameters <- c(popdyn_parameters, "p.ext_id","p.ext","z_mulambda","muOR","z_lambda","OR","z_intlambda","turnover")
    popdyn_const <- c(popdyn_const, list(ndate = ndate))
    if (FALSE %in% quo_is_null(var_guild)) {
      popdyn_parameters <- c(popdyn_parameters,"z_intlambda_gui","z_lambda_gui","z_mulambda_gui","GOR","muGOR")
    }
    if (isTRUE(modenv)) {
      popdyn_parameters <- popdyn_parameters[!popdyn_parameters %in% c("p.col","p.per","p.ext","p.col_id","p.per_id","p.ext_id")]
      if (FALSE %in% quo_is_null(var_envP) & FALSE %in% quo_is_null(var_envC)) {
        varenv <- colnames(select(df, !!var_envP))
        nvar <- length(varenv)
        varcol <- c(varenv[varenv %in% colnames(select(df, !!var_envC))],colnames(select(df, !!var_envC))[!colnames(select(df, !!var_envC)) %in% varenv])
        nr <- which(colnames(df) %in% varcol)
        ncol <- length(nr)
        col <- sapply(nr, function(i) tapply(df[,i], list(df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique),simplify = "array")
        popdyn_parameters <- c(popdyn_parameters, "alpha_per", "beta_per","alpha_col", "beta_col")
        popdyn_const <- c(popdyn_const, list(ncol = ncol))
        popdyn_data <- c(popdyn_data, list(col = col))
        popdyn_inits <- popdyn_inits[names(popdyn_inits) %in% c("p.per","p.col","epsilon_per","epsilon_col") == FALSE]
        popdyn_inits <- c(popdyn_inits, list(alpha_per = rep(0,ntax), beta_per = matrix(0, ntax, nvar), alpha_col = rep(0,ntax), beta_per = matrix(0, ntax, ncol)))
      } else if (FALSE %in% quo_is_null(var_envP)) {
        varenv <- colnames(select(df, !!var_envP))
        nvar <- length(varenv)
        popdyn_parameters <- c(popdyn_parameters, "alpha_per", "beta_per","p.col","p.col_id")
        popdyn_inits <- popdyn_inits[names(popdyn_inits) %in% c("p.col","epsilon_col") == FALSE]
        popdyn_inits <- c(popdyn_inits, list(alpha_per = rep(0,ntax), beta_per = matrix(0, ntax, nvar)))
      } else if (FALSE %in% quo_is_null(var_envC)) {
        varenv <- colnames(select(df, !!var_envC))
        nvar <- length(varenv)
        popdyn_parameters <- c(popdyn_parameters, "alpha_col", "beta_col","p.per","p.per_id")
        popdyn_inits <- popdyn_inits[names(popdyn_inits) %in% c("p.per","epsilon_per") == FALSE]
        popdyn_inits <- c(popdyn_inits, list(alpha_col = rep(0,ntax), beta_col = matrix(0, ntax, ncol)))
      } else {
        varenv <- colnames(select(df, !!var_envO))
        nvar <- length(varenv)
        popdyn_parameters <- c(popdyn_parameters, "alpha", "beta")
        popdyn_inits <- popdyn_inits[names(popdyn_inits) %in% c("p.per","p.col","epsilon_per","epsilon_col") == FALSE]
        popdyn_inits <- c(popdyn_inits, list(alpha = rep(0,ntax), beta = matrix(0, ntax, nvar)))
      }
      nr <- which(colnames(df) %in% varenv)
      var <- sapply(nr, function(i) tapply(df[,i], list(df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique),simplify = "array")
      popdyn_data <- c(popdyn_data,list(var = var))
      popdyn_const <- c(popdyn_const, list(nvar = nvar))
    } 
  }
  #-----------------------------------------------------------------------------
  if (isTRUE(grow)) {
    S <- array(1, dim = c(ntax, length(id), ntime))
    if (FALSE %in% quo_is_null(var_surf)) {
      S <- tapply(df[,quo_name(var_surf)], list(df[,quo_name(var_tax)],df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique)
    }
    popdyn_data <- c(popdyn_data, list(S = S))
    #---------------------------------------------------------------------------
    if (isTRUE(alt)) {
      pro <- array(1, dim = c(ntax,length(id),ntime))
      if (FALSE %in% quo_is_null(var_pro)) {
       npro <- summarise(df, across(quo_name(var_pro), ~length(unique(.x)), .names = "n"), .by=c(!!var_id,!!var_tax)) %>%
         filter(n > 1)
       for (i in 1:nrow(npro)) {
         for (t in 2:ntime) {
           d <- filter(df, df[,quo_name(var_tax)] %in% npro[i,quo_name(var_tax)], df[,quo_name(var_id)] %in% npro[i,quo_name(var_id)], df[,quo_name(var_tmp)] %in% time[c(t-1,t)]) %>%
             pull(!!var_pro) %>% unique()
           if (length(d) > 1) {
             pro[which(tax %in% npro[i,quo_name(var_tax)]),which(id %in% npro[i,quo_name(var_id)]),t] <- 0
           }
         }
       }
      }
      popdyn_data <- c(popdyn_data, list(pro = pro))
    }
    #---------------------------------------------------------------------------
    if (isTRUE(modenvG)) {
      varoccup <- c(varenv,varcol)
      vargrow <- c(varoccup[varoccup %in% colnames(select(df, !!var_grow))],colnames(select(df, !!var_grow))[!colnames(select(df, !!var_grow)) %in% varoccup])
      nr <- which(colnames(df) %in% vargrow)
      ngvar <- length(nr)
      gvar <- sapply(nr, function(i) tapply(df[,i], list(df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique),simplify = "array")
      popdyn_data <- c(popdyn_data,list(gvar = gvar))
      popdyn_const <- c(popdyn_const, list(ngvar = ngvar))
    }
    #---------------------------------------------------------------------------
    if (FALSE %in% quo_is_null(var_cnt)) {
      y <- tapply(df[,quo_name(var_cnt)], list(df[,quo_name(var_tax)],df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique)
      y <- replace(y, which(y == 0), 1)
      C <- log(replace(y, which(is.na(y)), 1))
      popdyn_data <- c(popdyn_data, list(y = y))
      popdyn_parameters <- c(popdyn_parameters, "N_PGR","N_muPGR","N_mulambda","N_lambda","N_intlambda","N_mulambda_id","N_lambda_id","N_intlambda_id","N")
      if (isTRUE(modenvG)) {
         popdyn_inits <- c(popdyn_inits, list(C = C, alpha_N = rep(0, ntax), beta_N = nimMatrix(1, ntax, ngvar), tauN = nimMatrix(1, ntax, length(id))))
         popdyn_parameters <- c(popdyn_parameters, "alpha_N","beta_N")
      } else {
        popdyn_inits <- c(popdyn_inits, list(C = C, muN = matrix(1, ntax, length(id)), tauN = matrix(1, ntax, length(id))))
      }
      if (FALSE %in% quo_is_null(var_guild)) {
        popdyn_parameters <- c(popdyn_parameters, "N_lambda_gui","N_mulambda_gui","N_intlambda_gui","N_GGR","N_muGGR")
      }
    }
    #-----------------------------------------------------------------------------
    if (FALSE %in% quo_is_null(var_wei)) {
      w <- tapply(df[,quo_name(var_wei)], list(df[,quo_name(var_tax)],df[,quo_name(var_id)],df[,quo_name(var_tmp)]), unique)
      w <- replace(w, which(w == 0), 1)
      W <- log(replace(w, is.na(w), 1))
      popdyn_data <- c(popdyn_data, list(w = w))
      popdyn_parameters <- c(popdyn_parameters, "B_PGR","B_muPGR","B_mulambda","B_lambda","B_intlambda","B_mulambda_id","B_lambda_id","B_intlambda_id","B")
      if (isTRUE(modenvG)) {
        popdyn_inits <- c(popdyn_inits, list(W = W, alpha_B = rep(0, ntax), beta_B = nimMatrix(1, ntax, ngvar), tauB = nimMatrix(1, ntax, length(id))))
        popdyn_parameters <- c(popdyn_parameters, "alpha_B","beta_B")
      } else {
        popdyn_inits <- c(popdyn_inits, list(W = W, muB = nimMatrix(1, ntax, length(id)), tauB = nimMatrix(1, ntax, length(id))))
      }
      if (FALSE %in% quo_is_null(var_guild)) {
        popdyn_parameters <- c(popdyn_parameters, "B_lambda_gui","B_mulambda_gui","B_intlambda_gui","B_GGR","B_muGGR")
      }
    }
  }
  #-----------------------------------------------------------------------------
  popdyn_int <- list(taxa = tax, id = id, time = time, region = region, guild = guild, varenv = varenv, varcol = varcol, vargrow = vargrow, start = start, end = end, protocol = protocol, vardet = vardet)
  return(list(popdyn_data = popdyn_data, popdyn_const = popdyn_const, popdyn_parameters = popdyn_parameters, popdyn_inits = popdyn_inits, popdyn_int = popdyn_int))
}
