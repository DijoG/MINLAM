# usethis::use_pipe()
#' Check packages
#' 
#' @return installs missing packages
#' @export
check_PACKS <- function() {
  required = c("tidyverse", "furrr", "multimode")
  installable = required[!(required %in% installed.packages()[,"Package"])]
  if(length(installable) > 0) {
    install.packages(installable)
    cat("All packages installed.")
  }
  else {
    cat("No need, packages already installed.\n")
  }
}

#' Extract mode statistics from a mode forest
#'
#' @param y vector, input data of a distribution
#' @param nmod numeric, number of subpopulations/subgroups  
#' @return data frame 
#' @export
get_MODES <- function(y, nmod) {
  mm = multimode::locmodes(y, mod0 = nmod)
  modes = mm$locations[seq(1, length(mm$locations), by = 2)]
  mode_df = data.frame(Est_Mode = modes,
                       Group = 1:length(modes))
  return(mode_df)
}

#' Group/merge modes if they are within a given range
#'
#' @param df data frame containing samples of a distribution in its 'Est_Mode' variable
#' @param within numeric, range
#' @return data frame 
#' @export
group_MODES <- function(df, within = 0.03) {
  df = df %>%
    dplyr::arrange(Est_Mode) %>%
    dplyr::mutate(Group = cumsum(c(1, diff(Est_Mode) > within)))
  
  df %>%
    dplyr::group_by(Group) %>%
    dplyr::summarise(Est_Mode = mean(Est_Mode), .groups = "drop")
}

#' Obtain the number of groups based on bandwidth selection
#'
#' @param y vector, input data of a distribution
#' @return data frame 
#' @export
get_NGRP <- function(y) {
  bwRT = stats::bw.nrd(y)                 # Scott's 1.06
  bwSJ = stats::bw.SJ(y, method = "dpi")  # Sheather & Jones' direct plug-in method
  bwBCV = stats::bw.bcv(y)                # Biased cross-validation by Jonas and Kappenman (1991)
  
  nrd_nmod = multimode::nmodes(y, bw = bwRT)
  sj_nmod = multimode::nmodes(y, bw = bwSJ)
  bcv_nmod = multimode::nmodes(y, bw = bwBCV)
  
  return(data.frame(Method = c("nrd", "bcv", "dpi"),
                    n_grp = c(nrd_nmod, bcv_nmod, sj_nmod)))
}

#' Compute mode based on density estimation
#'
#' @param x vector, samples of a distribution
#' @return numeric 
#' @export
get_mode <- function(x) {
  x = as.vector(x)
  
  if (all(x == round(x))) {
    return(as.numeric(names(which.max(table(x)))))
  } else {
    dx = density(x)
    return(dx$x[which.max(dx$y)])
  }
}

#' Fit an INLA model given response variable y and group assignments
#' Created by Virgilio Gómez-Rubio, https://becarioprecario.bitbucket.io/inla-gitbook/ch-mixture.html
#' 
#' @param y numeric, dependent variable to predict
#' @param grp numeric, estimated subpopulation/subgroup number
#' @param prior_means numeric, internally defined prior means of each subgroup
#' @return data frame 
#' @export
fit_inla_internal <- function(y, grp, prior_means) {
  
  unique_groups = unique(grp)
  n_grp = length(unique_groups)
  
  yy = matrix(NA, ncol = n_grp, nrow = length(y))
  for (i in seq_along(unique_groups)) {
    idx = which(grp == unique_groups[i])
    yy[idx, i] = y[idx]
  }
  
  x = ifelse(!is.na(yy), 1, NA)
  
  d = list(y = yy, x = x)
  
  mod = INLA::inla(y ~ -1 + x, data = d,
                   family = rep("gaussian", n_grp),
                   control.fixed = list(mean = prior_means, prec = 1))
  
  return(list(model = mod, mlik = mod$mlik[1, 1]))
}

#' Compute probabilities for group membership
#' Created by Virgilio Gómez-Rubio, https://becarioprecario.bitbucket.io/inla-gitbook/ch-mixture.html
#' 
#' @param z vector of subpopulation/subgroup membership
#' @return numeric vector
#' @export
get_probs <- function(z) {
  tab = table(z)
  probs = rep(0, n_grp)
  probs[as.integer(names(tab))] = tab / sum(tab)
  return(probs)
}

#' Random group assignment generator
#' Created by Virgilio Gómez-Rubio, https://becarioprecario.bitbucket.io/inla-gitbook/ch-mixture.html
#' 
#' @param z initial vector containing subpopulation/subgroup membership (proposed distribution)
#' @return list 
#' @export
rq_Z <- function(z) {
  m.aux = z$m$model
  means = m.aux$summary.fixed[, "mean"]
  precs = m.aux$summary.hyperpar[, "mean"]
  
  ws = get_probs(z$z)
  
  z.sim = sapply(seq_along(z$z), function(X) {
    aux = ws * dnorm(y[X], means, scale_sigma * sqrt(1 / precs))
    sample(1:n_grp, 1, prob = aux / sum(aux))
  })
  
  z.model = fit_inla_internal(y, z.sim, prior_means = prior_means)
  
  return(list(z = z.sim, m = z.model))
}

#' Compute density transition probability 
#' Created by Virgilio Gómez-Rubio, https://becarioprecario.bitbucket.io/inla-gitbook/ch-mixture.html
#' 
#' @param z.old initial vector containing subpopulation/subgroup membership (z)
#' @param z.new transitioned vector of subpopulation/subgroup membership (by inla)
#' @param log TRUE triggers INLA to conduct logging
#' @return vector 
#' @export
dq_Z <- function(z.old, z.new, log = TRUE) {
  m.aux = z.old$m$model
  means = m.aux$summary.fixed[, "mean"]
  precs = m.aux$summary.hyperpar[, "mean"]
  
  ww = get_probs(z.old$z)
  
  z.probs = sapply(seq_along(y), function(X) {
    aux = ww * dnorm(y[X], means, scale_sigma * sqrt(1 / precs))
    return((aux / sum(aux))[z.new$z[X]])
  })
  
  return(if (log) sum(log(z.probs)) else prod(z.probs))
}

#' Prior probability of z
#' Created by Virgilio Gómez-Rubio, https://becarioprecario.bitbucket.io/inla-gitbook/ch-mixture.html
#' 
#' @param z initial vector containing subpopulation/subgroup membership 
#' @param log logical, by default TRUE
#' @return vector 
#' @export
prior_Z <- function(z, log = TRUE) {
  res = log(1 / n_grp) * length(z$z)
  return(if (log) res else exp(res))
}

#' Wrapper for INLA model fitting
#' Created by Virgilio Gómez-Rubio, https://becarioprecario.bitbucket.io/inla-gitbook/ch-mixture.html
#' 
#' @param y  numeric, dependent variable to predict
#' @param grp numeric, number of subpopulations/subgroups
#' @return vector 
#' @export
fit_inla <- function(y, grp) {
  return(grp$m)
}

#' INLABMA::INLAMH (Metropolitan-Hastings sampling algorithm)
#' Created by Virgilio Gómez-Rubio, https://becarioprecario.bitbucket.io/inla-gitbook/ch-mixture.html
#' Source: https://rdrr.io/rforge/INLABMA/src/R/INLAMH.R
#' 
#' @param d data to be passed to inla() in a suitable format
#' @param fit.inla see above
#' @param b.init see above, nitial values of the model parameters
#' @param rq output of rq_Z(), sample from the proposal distribution 
#' @param dq output of dq_Z(), density of the proposal distribution
#' @param prior prior on the model parameters
#' @param n.sim number of simulations
#' @param n.burin number simulations for burn in
#' @param n.thin thinning after burn in
#' @param n.errors total number of allowed INLA errors 
#' @param verbose by default FALSE, if TRUE information messages are printed
#' @return list of simulations
#' @export
inla_MH <- function(d, fit.inla, b.init, rq, dq, prior, n.sim = 200, n.burnin = 100, 
                    n.thin = 1, n.errors = 20, verbose = FALSE) {
  b.sim = vector("list", n.sim)
  model.sim = vector("list", n.sim)
  acc.sim = rep(NA, n.sim)
  
  b.cur = b.init
  model.cur = fit.inla(d, b.cur)
  save.idx = 0
  total.iter = n.burnin + n.thin * n.sim
  error.count = 0
  
  for (i in 1:total.iter) {
    b.new = rq(b.cur)
    model.new = try(fit.inla(d, b.new), silent = TRUE)
    
    if (inherits(model.new, "try-error")) {
      error.count = error.count + 1
      if (verbose) cat("INLA error", error.count, "at iteration", i, "\n")
      if (error.count > n.errors) stop("Too many INLA failures.")
      if (i > n.burnin && save.idx > 100) {
        idx = sample(save.idx - 1, 1)
        b.cur = b.sim[[idx]]
        model.cur = model.sim[[idx]]
      } else {
        stop("INLA failed early.")
      }
      next
    }
    
    alpha = exp(min(0, model.new$mlik + prior(b.new) + dq(b.new, b.cur) -
                      (model.cur$mlik + prior(b.cur) + dq(b.cur, b.new))))
    acc.sim[i] = runif(1) < alpha
    
    if (acc.sim[i]) {
      b.cur = b.new
      model.cur = model.new
    }
    
    if (i > n.burnin && (i - n.burnin) %% n.thin == 0) {
      save.idx = save.idx + 1
      b.sim[[save.idx]] = b.cur
      model.sim[[save.idx]] = model.cur
    }
    
    if (verbose && i %% 100 == 0) {
      cat("Iteration", i, "completed at", as.character(Sys.time()), "\n")
    }
  }
  
  return(list(acc.sim = acc.sim, model.sim = model.sim, b.sim = b.sim))
}

#' Obtain sampled values of z
#' 
#' @param inla_MH_object  inla_MH() result
#' @return data frame 
#' @export
get_Z <- function(inla_MH_object) {
  zz_ = purrr::map(inla_MH_object$b.sim, ~ .x$z) %>% 
    reduce(rbind)  
  zz_probs = apply(zz_, 2, get_probs) 
  return(zz_probs)
}

#' Get the weights based on probabilities of class-belonging
#' 
#' @param get_z_object  get_Z() result
#' @param y numeric, dependent variable to predict
#' @param grp numeric, number of subpopulations/subgroups
#' @param main_class numeric, subpopulation/subgroup id 
#' @param df_prob by default FALSE (only summary), whether the estimated samples' probabilities to output
#' @return data frame or list of data farmes
#' @export
get_weighted_probs <- function(get_z_object, y, grp, main_class, df_prob = FALSE) {
  df_probs = 
    data.frame(y = y,
               Group = grp, 
               as.data.frame(t(get_z_object)) %>%
                 set_names(stringr::str_c("Group_", 1:length(unique(grp))))) %>%
    dplyr::mutate(Assigned_Group = max.col(dplyr::across(tidyselect::starts_with("Group_")), ties.method = "first"),
                  Min_Assigned = purrr::map_dbl(Assigned_Group, ~ min(y[grp == .])),
                  Max_Assigned = purrr::map_dbl(Assigned_Group, ~ max(y[grp == .])),
                  Mean_Assigned = purrr::map_dbl(Assigned_Group, ~ mean(y[grp == .])),
                  Mode_Assigned = purrr::map_dbl(Assigned_Group, ~ get_mode(y[grp == .])))
  
  wdf_probs =
    df_probs %>%
    dplyr::group_by(Assigned_Group) %>%
    dplyr::summarise(N = dplyr::n(),
                     Mean = mean(Mean_Assigned),
                     Mode = mean(Mode_Assigned)) %>%
    dplyr::mutate(Weight_ratio = N / sum(N),
                  Main_Class = main_class) %>%
    as.data.frame()
  
  if (df_prob) {
    df_probs = 
      df_probs %>%
      dplyr::mutate(Main_Class = main_class)
    return(list(df = df_probs,
                wdf = wdf_probs))
  } else {
    return(wdf_probs)
  }
}

#' MAIN function 
#' 
#' "dpi" ~ Direct plug-in by Sheather & Jones, 
#' "bcv" ~ Biased cross validation by Jonas and Kappenman (1991)
#' "nrd" ~ Scott's 1.06 (modified Silvermann's rule of thumb 1.34)
#' 
#' @param data  input data
#' @param varCLASS character, variable/column (categorical, character or numeric) name containing the subpopulations/subgroups
#' @param varY character, variable/column (numeric) containing the values assigned to the population/subpopulations
#' @param method character, which density estimator to apply: "nrd", "bcv", "dpi", by default "dpi", see above
#' @param within numeric, range as in group_MODES()
#' @param maxNGROUP numeric, after data inspection -> assumed maximum number of subpopulations/subgroups
#' @param df_prob logical, as in get_weighted_probs()
#' @param out_dir by default NULL, if character, path of the output directory to write csv tables 
#' @return data frame or list of data frames (stored as an R object or written as csvs)
#' @export
get_PROBCLASS_MH <- function(data, varCLASS, varY, method = "dpi", within = 0.03, maxNGROUP = 5, df_prob = FALSE, out_dir = NULL) {   
  
  if (!method %in% c("nrd", "bcv", "dpi")) {
    stop("Stopped, method can be: 'nrd', 'bcv' or 'dpi'.")}
  
  # scale_sigma 
  scale_sigma = 1
  
  # main classes
  mclass = unique(as.character(data[[varCLASS]]))
  
  ##> start loop
  OUT = list()
  for (i in seq_along(mclass)) {
    
    # Data preparation 
    #y = data[data[[varCLASS]] == mclass[i], varY]
    y = data %>% dplyr::filter(.data[[varCLASS]] == mclass[i]) %>% dplyr::pull(.data[[varY]])
    
    ##> Number of groups (subgroups/-populations in a multimodal distribution)
    # Initial groups 
    n_grp = get_NGRP(y) %>% dplyr::filter(Method == method) %>% dplyr::pull(n_grp)
    n_grp = min(max(n_grp, 3), maxNGROUP)
    
    # Get modes
    formodf = get_MODES(y, n_grp) %>% tidyr::drop_na()
    formodf$Group = seq_len(nrow(formodf))
    n_grp = nrow(formodf)
    
    # Group nearby modes
    formonearby = 
      formodf %>%
      dplyr::select(Group, Est_Mode) %>%
      dplyr::inner_join(dplyr::select(formodf, Group2 = Group, Est_Mode2 = Est_Mode), by = character()) %>%
      dplyr::filter(Group != Group2, abs(Est_Mode - Est_Mode2) <= within)
    if (nrow(formonearby) > 0) {
      formodf = group_MODES(formodf, within)
      n_grp = nrow(formodf)
    }
    
    # Obtain grp while updating groups
    #if (n_grp == 1) {grp = rep(1, length(y))}
    
    #if (n_grp > 1) {
    #  breaks = sort((formodf$Est_Mode[-nrow(formodf)] + formodf$Est_Mode[-1]) / 2)
    #  grp = as.numeric(cut(y, n_grp), breaks)
      
    #  ## 1)
    #  expected = 1:n_grp
    #  if (!all(expected %in% unique(grp))) {
    #    missing = setdiff(expected, unique(grp))
    #    min_missing = length(missing)
    #    if (min_missing == 1) {
    #      grp = grp - min_missing
    #      grp[grp == 0] = 1
    #      formodf = formodf[-missing,]
    #      formodf$Group = sort(unique(grp))
    #      n_grp = nrow(formodf)
    #    }
    #    if (min_missing > 1) {
    #      if (min_missing == 2) {
    #        grp = grp - min_missing
    #        if (min(grp) == 0) {
    #          grp[grp == 0] = 1
    #        }
    #        if (min(grp) == -1 & !0 %in% grp) {
    #          grp[grp == -1] = 1
    #        }
    #        if (min(grp) == -1 & 0 %in% grp) {
    #          grp[grp == 0] = 2
    #          grp[grp == -1] = 1
    #        }
    #        n_grp = length(unique(grp))
    #      } else {
    #        n_grp = length(unique(grp))
    #        formodf = get_MODES(y, n_grp) %>%
    #          dplyr::arrange(Est_Mode)
    #        breaks = sort((formodf$Est_Mode[-nrow(formodf)] + formodf$Est_Mode[-1]) / 2)
    #        grp = as.numeric(cut(y, n_grp), breaks)
    #      }
    #    }
    #  }
    #  ## 2)
    #  tab = as.data.frame(table(grp))
    #  if (any(tab$Freq == 1)) {
    #    wone = as.numeric(tab$grp)[tab$Freq == 1]
    #    yval = y[which(grp == wone)]
    #    formodf = formodf[-wone, ]
    #    grp[y == yval] = formodf$Group[which.min(abs(formodf$Est_Mode - yval))]
     #   n_grp = nrow(formodf)
    #  }
    #  ## 3)
    #  ugrp = sort(unique(grp))  
    #  if (!all(diff(ugrp) == 1)) {
    #    newugrp = setNames(seq_along(ugrp), ugrp)
    #    grp = as.integer(newugrp[as.character(grp)])
    #  }
    #}
    
    # Obtain grp while updating groups
    if (n_grp == 1) {
      grp = rep(1, length(y))
    } else {
      # Compute breakpoints for mode separation
      breaks = sort((formodf$Est_Mode[-nrow(formodf)] + formodf$Est_Mode[-1]) / 2)
      grp = as.numeric(cut(y, breaks = c(-Inf, breaks, Inf)))  
      
      # Ensure expected groups exist
      expected = seq_len(n_grp)
      present = unique(grp)
      missing = setdiff(expected, present)
      
      if (length(missing) > 0) {
        grp = grp - length(missing)
        grp[grp < 1] = 1
        formodf = formodf[-missing, , drop = FALSE]
        formodf$Group <- seq_len(nrow(formodf))  
        n_grp = nrow(formodf)
      }
      
      # Adjust groups with single observations
      tab = table(grp)
      singletons = as.numeric(names(tab)[tab == 1])
      
      if (length(singletons) > 0) {
        yvals = y[grp %in% singletons]
        closest_mode = vapply(yvals, function(v) which.min(abs(formodf$Est_Mode - v)), integer(1))
        grp[grp %in% singletons] = formodf$Group[closest_mode]
        n_grp = length(unique(grp))
      }
      
      # Ensure group labels are contiguous
      ugrp = sort(unique(grp))
      if (!all(diff(ugrp) == 1)) {
        grp = match(grp, ugrp)  # Faster remapping of group labels
      }
    }
    
    # Parameters
    prior_means = setNames(as.list(formodf$Est_Mode), stringr::str_c("x", seq_len(nrow(formodf))))
    grp_init = list(z = grp, m = fit_inla_internal(y, grp, prior_means))
    
    ##> Some crazy stuff ->
    assign("y", y, envir = .GlobalEnv)
    assign("n_grp", n_grp, envir = .GlobalEnv)
    assign("prior_means", prior_means, envir = .GlobalEnv)
    assign("scale_sigma", scale_sigma, envir = .GlobalEnv)
    Sys.sleep(1)
    
    ##> Metropolis-Hastings sampling
    cat("Processing", mclass[i], "class ----->\n")
    inlamh_res = inla_MH(y, fit_inla, grp_init, rq_Z, dq_Z, 
                         prior_Z, n.sim = 100, n.burnin = 20, n.thin = 5, verbose = TRUE)
    
    ##> Obtain sampled values of z
    Z = get_Z(inlamh_res)
    
    ##> Weight and get the final Mean and Mode
    OUT[[i]] = get_weighted_probs(Z, y, grp, main_class = mclass[i], df_prob)
    
    if (!is.null(out_dir)) {
      csvname = stringr::str_replace_all(mclass[i], "[^a-zA-Z0-9_-]", "_")
      if (df_prob) {
        readr::write_csv2(OUT[[i]]$df, stringr::str_c(out_dir, "/df_", csvname, ".csv"))
        readr::write_csv2(OUT[[i]]$wdf, stringr::str_c(out_dir, "/wdf_", csvname, ".csv"))
      } else {
        readr::write_csv2(OUT[[i]], stringr::str_c(out_dir, "/wdf_", csvname, ".csv"))
      }
    } else {
      next
    }
    
  } ##> end loop
  
  if (is.null(out_dir)) {
    if (df_prob) {
      df_probs = dplyr::bind_rows(lapply(OUT, `[[`, "df"))
      wdf_probs = dplyr::bind_rows(lapply(OUT, `[[`, "wdf"))
      return(list(df = df_probs,
                  wdf = wdf_probs))
    } else {
      return(dplyr::bind_rows(OUT))
    }
  }
}

#' Parallelizaton FRAME for the MAIN function
#' 
#' @param data  input data
#' @param varCLASS character, variable/column (categorical, character or numeric) name containing the subpopulations/subgroups
#' @param varY character, variable/column (numeric) containing the values assigned to the population/subpopulations
#' @param method character, which density estimator to apply: "nrd", "bcv", "dpi", by default "dpi", see above
#' @param within numeric, range as in group_MODES()
#' @param maxNGROUP numeric, after data inspection -> assumed maximum number of subpopulations/subgroups
#' @param df_prob logical, as in get_weighted_probs()
#' @param out_dir by default NULL, if character, path of the output directory to write csv tables 
#' @param n_workers numeric, number of cores to run simulations on, by default 4
#' @return data frame or list of data frames (stored as an R object or written as csvs)
#' @export
fuss_PARALLEL <- function(data, varCLASS, varY, method = "dpi", within = 0.03, maxNGROUP = 5, df_prob = FALSE, out_dir = NULL, n_workers = 4) {
  
  future::plan(multisession, workers = n_workers)
  
  result_list =
    furrr::future_map(data, ~ get_PROBCLASS_MH(.x, varCLASS, 
                                               varY, 
                                               method,
                                               within, 
                                               maxNGROUP,
                                               df_prob,
                                               out_dir))
  if (is.null(out_dir)) {
    if (df_prob) {
      df = purrr::map_dfr(result_list, "df")
      wdf = purrr::map_dfr(result_list, "wdf")
      return(list(df = df,
                  wdf = wdf))
    } else {
      return(dplyr::bind_rows(result_list))
    }
  }
}
