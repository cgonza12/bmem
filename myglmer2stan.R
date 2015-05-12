myglmer2stan <- function (formula, data, family = "gaussian", varpriors = "flat", 
          sample = TRUE, warmup = 5000, iter = 10000, chains = 1,mymodel=NULL, Ranef=T, 
          initmethod = "zero", extract = FALSE, calcDIC = TRUE, calcWAIC = FALSE, 
          verbose = TRUE, fixed_prefix = "beta_", vary_prefix = "vary_", 
          ...) 
{
  
  if(Ranef){
  parse_formula<- function (formula, data) 
  {
    nobars <- function(term) {
      if (!("|" %in% all.names(term))) 
        return(term)
      if (is.call(term) && term[[1]] == as.name("|")) 
        return(NULL)
      if (length(term) == 2) {
        nb <- nobars(term[[2]])
        if (is.null(nb)) 
          return(NULL)
        term[[2]] <- nb
        return(term)
      }
      nb2 <- nobars(term[[2]])
      nb3 <- nobars(term[[3]])
      if (is.null(nb2)) 
        return(nb3)
      if (is.null(nb3)) 
        return(nb2)
      term[[2]] <- nb2
      term[[3]] <- nb3
      term
    }
    findbars <- function(term) {
      if (is.name(term) || !is.language(term)) 
        return(NULL)
      if (term[[1]] == as.name("(")) 
        return(findbars(term[[2]]))
      if (!is.call(term)) 
        stop("term must be of class call")
      if (term[[1]] == as.name("|")) 
        return(term)
      if (length(term) == 2) 
        return(findbars(term[[2]]))
      c(findbars(term[[2]]), findbars(term[[3]]))
    }
    subbars <- function(term) {
      if (is.name(term) || !is.language(term)) 
        return(term)
      if (length(term) == 2) {
        term[[2]] <- subbars(term[[2]])
        return(term)
      }
      stopifnot(length(term) >= 3)
      if (is.call(term) && term[[1]] == as.name("|")) 
        term[[1]] <- as.name("+")
      for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
      term
    }
    hasintercept <- function(term) {
      attr(terms(term), "intercept") == 1
    }
    f_nobars <- nobars(formula)
    if (class(f_nobars) == "name" & length(f_nobars) == 1) {
      f_nobars <- nobars(as.formula(paste(deparse(formula), 
                                          "+ 1")))
    }
    fixef <- colnames(model.matrix(f_nobars, data))
    mdat <- model.matrix(subbars(formula), data)
    outcome_name <- deparse(f_nobars[[2]])
    outcome <- model.frame(f_nobars, data)[, 1]
    if (class(outcome) == "matrix") {
      outcome_name <- colnames(outcome)[1]
    }
    
    var <- findbars(formula)
    ranef <- list()
    for (i in 1:length(var)) {
      name <- all.vars(var[[i]][[3]])
      if (FALSE) {
        if (TRUE) {
          if (class(data[[name]]) != "integer") {
            stop(paste("Grouping variables must be integer type. '", 
                       name, "' is instead of type: ", class(data[[name]]), 
                       ".", sep = ""))
          }
          if (min(data[[name]]) != 1) 
            stop(paste("Group variable '", name, "' doesn't start at index 1.", 
                       sep = ""))
          ulist <- unique(data[[name]])
          diffs <- ulist[2:length(ulist)] - ulist[1:(length(ulist) - 
                                                       1)]
          if (any(diffs != 1)) 
            stop(paste("Group variable '", name, "' is not contiguous.", 
                       sep = ""))
        }
        else {
          mdat[, name] <- as.integer(as.factor(mdat[, 
                                                    name]))
        }
      }
      v <- var[[i]][[2]]
      if (class(v) == "numeric") {
        ranef[[name]] <- "(Intercept)"
      }
      else {
        f <- as.formula(paste("~", deparse(v), sep = ""))
        ranef[[name]] <- colnames(model.matrix(f, data))
      }
    }
    
    list(y = outcome, yname = outcome_name, fixef = fixef, ranef = ranef, 
         dat = as.data.frame(mdat))
  }
  }else{
    parse_formula<- function (formula, data) 
    {
      nobars <- function(term) {
        if (!("|" %in% all.names(term))) 
          return(term)
        if (is.call(term) && term[[1]] == as.name("|")) 
          return(NULL)
        if (length(term) == 2) {
          nb <- nobars(term[[2]])
          if (is.null(nb)) 
            return(NULL)
          term[[2]] <- nb
          return(term)
        }
        nb2 <- nobars(term[[2]])
        nb3 <- nobars(term[[3]])
        if (is.null(nb2)) 
          return(nb3)
        if (is.null(nb3)) 
          return(nb2)
        term[[2]] <- nb2
        term[[3]] <- nb3
        term
      }
      findbars <- function(term) {
        if (is.name(term) || !is.language(term)) 
          return(NULL)
        if (term[[1]] == as.name("(")) 
          return(findbars(term[[2]]))
        if (!is.call(term)) 
          stop("term must be of class call")
        if (term[[1]] == as.name("|")) 
          return(term)
        if (length(term) == 2) 
          return(findbars(term[[2]]))
        c(findbars(term[[2]]), findbars(term[[3]]))
      }
      subbars <- function(term) {
        if (is.name(term) || !is.language(term)) 
          return(term)
        if (length(term) == 2) {
          term[[2]] <- subbars(term[[2]])
          return(term)
        }
        stopifnot(length(term) >= 3)
        if (is.call(term) && term[[1]] == as.name("|")) 
          term[[1]] <- as.name("+")
        for (j in 2:length(term)) term[[j]] <- subbars(term[[j]])
        term
      }
      hasintercept <- function(term) {
        attr(terms(term), "intercept") == 1
      }
      f_nobars <- nobars(formula)
      if (class(f_nobars) == "name" & length(f_nobars) == 1) {
        f_nobars <- nobars(as.formula(paste(deparse(formula), 
                                            "+ 1")))
      }
      fixef <- colnames(model.matrix(f_nobars, data))
      mdat <- model.matrix(subbars(formula), data)
      outcome_name <- deparse(f_nobars[[2]])
      outcome <- model.frame(f_nobars, data)[, 1]
      if (class(outcome) == "matrix") {
        outcome_name <- colnames(outcome)[1]
      }
      if (formula == nobars(formula)) {
        ranef <- list()
      }
      else {
        var <- findbars(formula)
        ranef <- list()
        for (i in 1:length(var)) {
          name <- all.vars(var[[i]][[3]])
          if (FALSE) {
            if (TRUE) {
              if (class(data[[name]]) != "integer") {
                stop(paste("Grouping variables must be integer type. '", 
                           name, "' is instead of type: ", class(data[[name]]), 
                           ".", sep = ""))
              }
              if (min(data[[name]]) != 1) 
                stop(paste("Group variable '", name, "' doesn't start at index 1.", 
                           sep = ""))
              ulist <- unique(data[[name]])
              diffs <- ulist[2:length(ulist)] - ulist[1:(length(ulist) - 
                                                           1)]
              if (any(diffs != 1)) 
                stop(paste("Group variable '", name, "' is not contiguous.", 
                           sep = ""))
            }
            else {
              mdat[, name] <- as.integer(as.factor(mdat[, 
                                                        name]))
            }
          }
          v <- var[[i]][[2]]
          if (class(v) == "numeric") {
            ranef[[name]] <- "(Intercept)"
          }
          else {
            f <- as.formula(paste("~", deparse(v), sep = ""))
            ranef[[name]] <- colnames(model.matrix(f, data))
          }
        }
      }
      list(y = outcome, yname = outcome_name, fixef = fixef, ranef = ranef, 
           dat = as.data.frame(mdat))
    }
    
  }
  
  indent <- "    "
  REML <- TRUE
  vectorized <- c(1, 1, 1, 1, 0)
  names(vectorized) <- c("normal", "binomial", "poisson", 
                         "gamma", "ordered")
  log_sum_exp <- function(x) {
    xmax <- max(x)
    xsum <- sum(exp(x - xmax))
    xmax + log(xsum)
  }
  class2type <- function(col) {
    classes <- c("integer", "numeric")
    types <- c("int", "real")
    types[which(classes == class(col))]
  }
  undot <- function(astring) {
    astring <- gsub(".", "_", astring, fixed = TRUE)
    astring <- gsub(":", "_X_", astring, fixed = TRUE)
    astring <- gsub("(", "", astring, fixed = TRUE)
    astring <- gsub(")", "", astring, fixed = TRUE)
    astring
  }
  make_vary_text <- function(fp, vary_prefix = "vary_", suffix = "", 
                             f_num) {
    vterms <- c()
    for (i in 1:length(fp$ranef)) {
      name <- undot(names(fp$ranef)[i])
      factors <- fp$ranef[[i]]
      jstart <- 0
      if (f_num > 1) 
        jstart <- sum(cluster_vars[[name]][1:(f_num - 
                                                1)])
      total_factors <- sum(cluster_vars[[name]])
      for (j in 1:length(factors)) {
        if (length(factors) == 1 & total_factors == 
              1) {
          aterm <- paste(vary_prefix, name, "[", name, 
                         suffix, "[i]]", sep = "")
        }
        else {
          aterm <- paste(vary_prefix, name, "[", name, 
                         suffix, "[i],", j + jstart, "]", sep = "")
        }
        if (j > 1) {
          aterm <- paste(aterm, " * ", undot(factors[j]), 
                         suffix, "[i]", sep = "")
        }
        vterms <- c(vterms, aterm)
      }
    }
    vterms
  }
  make_fixed_text <- function(fp, fixed_prefix = "beta_", 
                              suffix = "", drop_intercept = FALSE) {
    fterms <- c()
    for (i in 1:length(fp$fixef)) {
      name <- undot(fp$fixef[i])
      if (name == "Intercept") {
        if (!drop_intercept) {
          aterm <- paste("Intercept", suffix, sep = "")
          fterms <- c(fterms, aterm)
        }
        else {
          aterm <- paste("0", sep = "")
          fterms <- c(fterms, aterm)
        }
      }
      else {
        aterm <- paste(fixed_prefix, name, suffix, " * ", 
                       name, suffix, "[i]", sep = "")
        fterms <- c(fterms, aterm)
      }
    }
    fterms
  }
  logistic <- function(x) {
    p <- 1/(1 + exp(-x))
    p <- ifelse(x == Inf, 1, p)
    p
  }
  logit <- function(x) {
    log(x/(1 - x))
  }
  dordlogit <- function(x, a, phi, log = FALSE) {
    a <- c(as.numeric(a), Inf)
    p <- logistic(a[x] - phi)
    na <- c(-Inf, a)
    np <- logistic(na[x] - phi)
    p <- p - np
    if (log == TRUE) 
      p <- log(p)
    p
  }
  dgamma2 <- function(x, mu, scale = NULL, log = FALSE) {
    dgamma(x, shape = mu/scale, scale = scale, log = log)
  }
  if (warmup >= iter) {
    stop("Number of iterations (iter) must exceed number of warmup steps (warmup).")
  }
  if (class(family) != "list") {
    family <- list(family)
  }
  for (f in 1:length(family)) {
    legal_families <- c("gaussian", "binomial", "poisson", 
                        "ordered", "gamma", "zigamma", "zipoisson")
    if (!(family[[f]] %in% legal_families)) {
      stop(paste("Family '", family, "' not recognized.", 
                 sep = ""))
    }
  }
  if (!(initmethod %in% c("lme4", "random", "zero"))) {
    stop(paste("Init method '", initmethod, "' not recognized.", 
               sep = ""))
  }
  if (!(varpriors %in% c("weak", "flat"))) {
    stop(paste("Variance priors '", varpriors, "' not recognized.", 
               sep = ""))
  }
  if (missing(formula) || missing(data)) {
    stop("Must provide both formula and data objects.")
  }
  fixed_prefix <- undot(fixed_prefix)
  vary_prefix <- undot(vary_prefix)
  if (sample == TRUE) 
    require(rstan)
  if (class(formula) != "list") {
    formula <- list(formula)
  }
  if (family[[1]] == "zigamma" | family[[1]] == "zipoisson") {
    if (length(formula) == 1) {
      outvar <- deparse(formula[[1]][[2]])
      y <- data[[outvar]]
      y1 <- ifelse(y == 0, 1, 0)
      if (family[[1]] == "zigamma") {
        y2 <- ifelse(y == 0, NA, y)
      }
      if (family[[1]] == "zipoisson") {
        y2 <- as.integer(y)
      }
      zname <- paste(outvar, "_zeros", sep = "")
      nzname <- paste(outvar, "_nonzero", sep = "")
      if (family[[1]] == "zipoisson") {
        nzname <- paste(outvar, "_count", sep = "")
      }
      data[[zname]] <- as.integer(y1)
      data[[nzname]] <- y2
      f1 <- as.formula(paste(zname, " ~ ", deparse(formula[[1]][[3]]), 
                             sep = ""))
      f2 <- as.formula(paste(nzname, " ~ ", deparse(formula[[1]][[3]]), 
                             sep = ""))
      formula <- list(f1, f2)
      if (family[[1]] == "zigamma") {
        family <- list("binomial", "gamma")
      }
      if (family[[1]] == "zipoisson") {
        family <- list("binomial", "poisson")
      }
    }
  }
  if (family[[1]] == "zigamma" & length(formula) == 2) {
    family <- list("binomial", "gamma")
  }
  if (family[[1]] == "zipoisson" & length(formula) == 2) {
    family <- list("binomial", "poisson")
  }
  fp <- list()
  for (f in 1:length(formula)) {
    fp[[f]] <- parse_formula(formula[[f]], data)
  }
  num_formulas <- length(formula)
  var_suffix <- ""
  if (num_formulas > 1) {
    var_suffix <- paste("_", 1:num_formulas, sep = "")
  }
  if (length(family) == 2) {
    if (family[[1]] == "binomial" & family[[2]] == "gamma") {
      var_suffix <- c("_pi", "_mu")
    }
    if (family[[1]] == "binomial" & family[[2]] == "poisson") {
      var_suffix <- c("_z", "_n")
    }
  }
  fit <- list()
  if (initmethod == "lme4") {
    require(lme4)
    if (verbose == TRUE) {
      message("Fitting lme4 model")
      flush.console()
    }
    for (f in 1:num_formulas) {
      use.family <- family[[f]]
      if (family[[f]] == "ordered") {
        use.family <- "gaussian"
      }
      if (family[[f]] == "gamma") {
        use.family <- gaussian(link = "log")
      }
      fit[[f]] <- glmer(formula[[f]], data = data, family = use.family)
    }
  }
  m_data <- "data{\n"
  for (i in 1:num_formulas) {
    m_data <- paste(m_data, indent, "int N", var_suffix[i], 
                    ";\n", sep = "")
  }
  m_transdata1 <- "transformed data{\n"
  m_transdata2 <- ""
  flag_transdata <- FALSE
  m_pars <- "parameters{\n"
  m_transpars1 <- "transformed parameters{\n"
  m_transpars2 <- ""
  flag_transpars <- FALSE
  m_model <- "model{\n"
  for (i in 1:num_formulas) {
    m_model <- paste(m_model, indent, "real vary", var_suffix[i], 
                     "[N", var_suffix[i], "];\n", sep = "")
    m_model <- paste(m_model, indent, "real glm", var_suffix[i], 
                     "[N", var_suffix[i], "];\n", sep = "")
  }
  m_gq <- ""
  flag_intercept <- FALSE
  for (f in 1:num_formulas) {
    if (family[[f]] == "binomial") {
      if (class(fp[[f]]$y) == "matrix") {
        vname <- fp[[f]]$yname
        bin_tot <- as.integer(fp[[f]]$y[, 1] + fp[[f]]$y[, 
                                                         2])
        fp[[f]]$y <- as.integer(fp[[f]]$y[, 1])
        fp[[f]]$dat[, "bin_total"] <- bin_tot
      }
      else {
        fp[[f]]$y <- as.integer(fp[[f]]$y)
        fp[[f]]$dat[, "bin_total"] <- as.integer(rep(1, 
                                                     nrow(fp[[f]]$dat)))
      }
    }
    if (family[[f]] == "poisson" || family[[f]] == "ordered") {
      fp[[f]]$y <- as.integer(fp[[f]]$y)
    }
    if (family[[f]] == "ordered") {
      Kname <- paste("K", var_suffix[f], sep = "")
      m_data <- paste(m_data, indent, "int<lower=2> ", 
                      Kname, ";\n", sep = "")
    }
    if (length(fp[[f]]$ranef) > 0) {
      for (i in 1:length(fp[[f]]$ranef)) {
        fname <- names(fp[[f]]$ranef)[i]
        fp[[f]]$dat[, fname] <- as.integer(fp[[f]]$dat[, 
                                                       fname])
      }
    }
  }
  for (i in 1:num_formulas) {
    outvar <- undot(fp[[i]]$yname)
    vtype <- class2type(fp[[i]]$y)
    nname <- paste("N", var_suffix[i], sep = "")
    m_data <- paste(m_data, indent, vtype, " ", outvar, 
                    var_suffix[i], "[", nname, "];\n", sep = "")
    for (j in 1:ncol(fp[[i]]$dat)) {
      vtype <- class2type(fp[[i]]$dat[, j])
      vname <- undot(colnames(fp[[i]]$dat)[j])
      if (vname != "Intercept") {
        nname <- paste("N", var_suffix[i], sep = "")
        m_data <- paste(m_data, indent, vtype, " ", 
                        vname, var_suffix[i], "[", nname, "];\n", 
                        sep = "")
      }
    }
  }
  cluster_vars <- list()
  cluster_size <- list()
  for (i in 1:num_formulas) {
    if (length(fp[[i]]$ranef) > 0) {
      for (j in 1:length(fp[[i]]$ranef)) {
        name <- names(fp[[i]]$ranef)[j]
        cluster_name <- undot(name)
        num_factors <- length(fp[[i]]$ranef[[j]])
        num_groups <- length(unique(fp[[i]]$dat[, name]))
        cluster_vars[[cluster_name]] <- c(cluster_vars[[cluster_name]], 
                                          num_factors)
        cluster_size[[cluster_name]] <- max(cluster_size[[cluster_name]], 
                                            num_groups)
      }
    }
  }
  if (length(cluster_vars) > 0) {
    for (i in 1:length(cluster_vars)) {
      var_name <- names(cluster_vars)[i]
      m_data <- paste(m_data, indent, "int N_", var_name, 
                      ";\n", sep = "")
    }
  }
  pars_list <- ""
  for (f in 1:num_formulas) {
    for (i in 1:length(fp[[f]]$fixef)) {
      var_name <- undot(fp[[f]]$fixef[i])
      if (var_name != "Intercept") {
        var_name <- paste(fixed_prefix, var_name, var_suffix[f], 
                          sep = "")
        m_pars <- paste(m_pars, indent, "real ", var_name, 
                        ";\n", sep = "")
        pars_list <- c(pars_list, var_name)
      }
      else {
        if (family[[f]] != "ordered") {
          flag_intercept <- TRUE
          var_name <- paste(var_name, var_suffix[f], 
                            sep = "")
          m_pars <- paste(m_pars, indent, "real ", var_name, 
                          ";\n", sep = "")
          pars_list <- c(pars_list, var_name)
        }
      }
    }
    if (family[[f]] == "gaussian") {
      signame <- paste("sigma", var_suffix[f], sep = "")
      m_pars <- paste(m_pars, indent, "real<lower=0> ", 
                      signame, ";\n", sep = "")
      pars_list <- c(pars_list, signame)
    }
    if (family[[f]] == "gamma") {
      thetaname <- paste("theta", var_suffix[f], sep = "")
      m_pars <- paste(m_pars, indent, "real<lower=0.001> ", 
                      thetaname, ";\n", sep = "")
      pars_list <- c(pars_list, thetaname)
    }
    if (family[[f]] == "ordered") {
      cutsname <- paste("cutpoints", var_suffix[f], sep = "")
      Kname <- paste("K", var_suffix[f], sep = "")
      m_pars <- paste(m_pars, indent, "ordered[", Kname, 
                      "-1] ", cutsname, ";\n", sep = "")
      pars_list <- c(pars_list, cutsname)
    }
  }
  if (length(cluster_vars) > 0) {
    for (i in 1:length(cluster_vars)) {
      cluster_name <- names(cluster_vars)[i]
      num_effects <- sum(cluster_vars[[i]])
      num_groups <- cluster_size[[i]]
      var_name <- paste(vary_prefix, cluster_name, sep = "")
      var_name <- paste(var_name, "[N_", cluster_name, 
                        "]", sep = "")
      ktype <- "real "
      if (num_effects > 1) 
        ktype <- paste("vector[", num_effects, "] ", 
                       sep = "")
      m_pars <- paste(m_pars, indent, ktype, var_name, 
                      ";\n", sep = "")
      pars_list <- c(pars_list, paste(vary_prefix, cluster_name, 
                                      sep = ""))
      if (num_effects == 1) {
        m_pars <- paste(m_pars, indent, "real<lower=0> sigma_", 
                        cluster_name, ";\n", sep = "")
        pars_list <- c(pars_list, paste("sigma_", cluster_name, 
                                        sep = ""))
      }
      else {
        if (varpriors != "weak") {
          m_pars <- paste(m_pars, indent, "cov_matrix[", 
                          num_effects, "] Sigma_", cluster_name, ";\n", 
                          sep = "")
          pars_list <- c(pars_list, paste("Sigma_", 
                                          cluster_name, sep = ""))
        }
        else {
          m_pars <- paste(m_pars, indent, "vector<lower=0>[", 
                          num_effects, "] sigma_", cluster_name, ";\n", 
                          sep = "")
          m_pars <- paste(m_pars, indent, "corr_matrix[", 
                          num_effects, "] Rho_", cluster_name, ";\n", 
                          sep = "")
          m_transpars1 <- paste(m_transpars1, indent, 
                                "cov_matrix[", num_effects, "] Sigma_", 
                                cluster_name, ";\n", sep = "")
          m_transpars2 <- paste(m_transpars2, indent, 
                                "Sigma_", cluster_name, " <- diag_matrix(sigma_", 
                                cluster_name, ") * Rho_", cluster_name, 
                                " * diag_matrix(sigma_", cluster_name, ");\n", 
                                sep = "")
          flag_transpars <- TRUE
          pars_list <- c(pars_list, paste("Sigma_", 
                                          cluster_name, sep = ""))
        }
      }
    }
  }
  m_pars <- paste(m_pars, "}\n", sep = "")
  m_model <- paste(m_model, indent, "// Priors\n", sep = "")
  for (f in 1:num_formulas) {
    for (i in 1:length(fp[[f]]$fixef)) {
      name <- undot(fp[[f]]$fixef[i])
      prefix <- fixed_prefix
      if (name == "Intercept") {
        if (family[[f]] != "ordered") {
          m_model <- paste(m_model, indent, name, var_suffix[f], 
                           " ~ normal( 0 , 100 );\n", sep = "")
        }
      }
      else {
        m_model <- paste(m_model, indent, prefix, name, 
                         var_suffix[f], " ~ normal( 0 , 100 );\n", 
                         sep = "")
      }
    }
  }
  if (length(cluster_vars) > 0) {
    for (i in 1:length(cluster_vars)) {
      num_effects <- sum(cluster_vars[[i]])
      name <- undot(names(cluster_vars)[i])
      if (num_effects == 1) {
        if (varpriors == "weak") 
          m_model <- paste(m_model, indent, "sigma_", 
                           name, " ~ gamma( 2 , 1e-4 );\n", sep = "")
        if (varpriors == "flat") 
          m_model <- paste(m_model, indent, "sigma_", 
                           name, " ~ uniform( 0 , 100 );\n", sep = "")
      }
      else {
        if (varpriors == "weak") {
          m_model <- paste(m_model, indent, "sigma_", 
                           name, " ~ gamma( 2 , 1e-4 );\n", sep = "")
          m_model <- paste(m_model, indent, "Rho_", 
                           name, " ~ lkj_corr( 1.5 );\n", sep = "")
        }
        if (varpriors == "conjugate") {
          m_model <- paste(m_model, indent, "Sigma_", 
                           name, " ~ inv_wishart( ", num_effects + 
                             1, " , Omega_", name, " );\n", sep = "")
          m_data <- paste(m_data, indent, "cov_matrix[", 
                          num_effects, "] Omega_", name, ";\n", sep = "")
        }
      }
    }
  }
  for (f in 1:num_formulas) {
    if (family[[f]] == "gaussian") {
      signame <- paste("sigma", var_suffix[f], sep = "")
      if (varpriors == "weak") 
        m_model <- paste(m_model, indent, signame, " ~ gamma( 2 , 1e-4 );\n", 
                         sep = "")
      if (varpriors == "flat") 
        m_model <- paste(m_model, indent, signame, " ~ uniform( 0 , 100 );\n", 
                         sep = "")
    }
    if (family[[f]] == "gamma") {
      thetaname <- paste("theta", var_suffix[f], sep = "")
      m_model <- paste(m_model, indent, thetaname, " ~ uniform( 0.001 , 20 );\n", 
                       sep = "")
    }
  }
  if (length(cluster_vars) > 0) {
    m_model <- paste(m_model, indent, "// Varying effects\n", 
                     sep = "")
    for (i in 1:length(cluster_vars)) {
      cluster_name <- undot(names(cluster_vars)[i])
      num_effects <- sum(cluster_vars[[i]])
      if (num_effects == 1) {
        m_model <- paste(m_model, indent, "for ( j in 1:N_", 
                         cluster_name, " ) ", vary_prefix, cluster_name, 
                         "[j] ~ normal( 0 , sigma_", cluster_name, 
                         " );\n", sep = "")
      }
      else {
        m_model <- paste(m_model, indent, "for ( j in 1:N_", 
                         cluster_name, " ) ", vary_prefix, cluster_name, 
                         "[j] ~ multi_normal( zeros_", cluster_name, 
                         " , Sigma_", cluster_name, " );\n", sep = "")
        m_transdata1 <- paste(m_transdata1, indent, 
                              "vector[", num_effects, "] zeros_", cluster_name, 
                              ";\n", sep = "")
        m_transdata2 <- paste(m_transdata2, indent, 
                              "for ( i in 1:", num_effects, " ) zeros_", 
                              cluster_name, "[i] <- 0;\n", sep = "")
        flag_transdata <- TRUE
      }
    }
  }
  m_model <- paste(m_model, indent, "// Fixed effects\n", 
                   sep = "")
  vary_text <- list()
  fixed_text <- list()
  for (f in 1:num_formulas) {
    formula_collapse_text <- "\n                + "
    if (length(fp[[f]]$ranef) > 0) {
      vary_terms <- make_vary_text(fp[[f]], vary_prefix, 
                                   var_suffix[f], f)
      vary_text[[f]] <- paste(vary_terms, collapse = formula_collapse_text)
      vary_text[[f]] <- paste(indent, indent, "vary", 
                              var_suffix[f], "[i] <- ", vary_text[[f]], ";\n", 
                              sep = "")
    }
    else {
      vary_text[[f]] <- ""
    }
    fixed_text[[f]] <- make_fixed_text(fp[[f]], fixed_prefix, 
                                       var_suffix[f], drop_intercept = family[[f]] == "ordered")
    fixed_text[[f]] <- paste(fixed_text[[f]], collapse = formula_collapse_text)
    if (length(fp[[f]]$ranef) > 0) {
      fixed_text[[f]] <- paste(indent, indent, "glm", 
                               var_suffix[f], "[i] <- vary", var_suffix[f], 
                               "[i] + ", fixed_text[[f]], ";\n", sep = "")
    }
    else {
      fixed_text[[f]] <- paste(indent, indent, "glm", 
                               var_suffix[f], "[i] <- ", fixed_text[[f]], ";\n", 
                               sep = "")
    }
    m_model <- paste(m_model, indent, "for ( i in 1:N", 
                     var_suffix[f], " ) {\n", vary_text[[f]], fixed_text[[f]], 
                     sep = "")
    out_var <- undot(fp[[f]]$yname)
    out_var <- paste(out_var, var_suffix[f], sep = "")
    if (family[[f]] == "binomial") {
      bintotname <- paste("bin_total", var_suffix[f], 
                          sep = "")
      if (FALSE) {
        m_model <- paste(m_model, "}\n", sep = "")
        m_model <- paste(m_model, indent, indent, out_var, 
                         " ~ bernoulli_logit( glm );", sep = "")
      }
      else {
        m_model <- paste(m_model, indent, indent, "glm", 
                         var_suffix[f], "[i] <- inv_logit( glm", var_suffix[f], 
                         "[i] );\n", sep = "")
        m_model <- paste(m_model, indent, "}\n", sep = "")
        m_model <- paste(m_model, indent, out_var, " ~ binomial( ", 
                         bintotname, " , glm", var_suffix[f], " );\n", 
                         sep = "")
      }
    }
    if (family[[f]] == "gaussian") {
      signame <- paste("sigma", var_suffix[f], sep = "")
      m_model <- paste(m_model, indent, "}\n", sep = "")
      m_model <- paste(m_model, indent, out_var, " ~ normal( glm", 
                       var_suffix[f], " , ", signame, " );\n", sep = "")
    }
    if (family[[f]] == "poisson") {
      m_model <- paste(m_model, indent, indent, "glm[i] <- exp( glm", 
                       var_suffix[f], "[i] );\n", sep = "")
      m_model <- paste(m_model, indent, "}\n", sep = "")
      m_model <- paste(m_model, indent, out_var, " ~ poisson( glm", 
                       var_suffix[f], " );\n", sep = "")
    }
    if (family[[f]] == "ordered") {
      cutsname <- paste("cutpoints", var_suffix[f], sep = "")
      m_model <- paste(m_model, indent, indent, out_var, 
                       "[i] ~ ordered_logistic( glm", var_suffix[f], 
                       "[i] , ", cutsname, " );\n    }\n", sep = "")
    }
    if (family[[f]] == "gamma") {
      thetaname <- paste("theta", var_suffix[f], sep = "")
      m_model <- paste(m_model, indent, indent, "glm", 
                       var_suffix[f], "[i] <- exp( glm", var_suffix[f], 
                       "[i] )*", thetaname, ";\n", sep = "")
      m_model <- paste(m_model, indent, "}\n", sep = "")
      m_model <- paste(m_model, indent, out_var, " ~ gamma( glm", 
                       var_suffix[f], " , ", thetaname, " );\n", sep = "")
    }
  }
  m_model <- paste(m_model, "}\n", sep = "")
  if (calcDIC == TRUE) {
    m_gq <- "generated quantities{\n    real dev;\n"
    for (i in 1:num_formulas) {
      m_gq <- paste(m_gq, indent, "real vary", var_suffix[i], 
                    "[N", var_suffix[i], "];\n", sep = "")
      m_gq <- paste(m_gq, indent, "real glm", var_suffix[i], 
                    "[N", var_suffix[i], "];\n", sep = "")
    }
    m_gq <- paste(m_gq, indent, "dev <- 0;\n", sep = "")
    for (f in 1:num_formulas) {
      m_gq <- paste(m_gq, indent, "for ( i in 1:N", var_suffix[f], 
                    " ) {\n", sep = "")
      m_gq <- paste(m_gq, vary_text[[f]], fixed_text[[f]], 
                    sep = "")
      out_var <- undot(fp[[f]]$yname)
      out_var <- paste(out_var, var_suffix[f], sep = "")
      if (family[[f]] == "binomial") {
        bintotname <- paste("bin_total", var_suffix[f], 
                            sep = "")
        if (FALSE) {
          m_gq <- paste(m_gq, indent, indent, "dev <- dev + (-2) * bernoulli_log( ", 
                        out_var, "[i] , inv_logit(glm[i]) );\n", 
                        sep = "")
        }
        else {
          m_gq <- paste(m_gq, indent, indent, "dev <- dev + (-2) * binomial_log( ", 
                        out_var, "[i] , ", bintotname, "[i] , inv_logit(glm", 
                        var_suffix[f], "[i]) );\n", sep = "")
        }
      }
      if (family[[f]] == "gaussian") {
        signame <- paste("sigma", var_suffix[f], sep = "")
        m_gq <- paste(m_gq, indent, indent, "dev <- dev + (-2) * normal_log( ", 
                      out_var, "[i] , glm", var_suffix[f], "[i] , ", 
                      signame, " );\n", sep = "")
      }
      if (family[[f]] == "poisson") {
        m_gq <- paste(m_gq, indent, indent, "dev <- dev + (-2) * poisson_log( ", 
                      out_var, "[i] , exp(glm", var_suffix[f], "[i]) );\n", 
                      sep = "")
      }
      if (family[[f]] == "ordered") {
        cutsname <- paste("cutpoints", var_suffix[f], 
                          sep = "")
        m_gq <- paste(m_gq, indent, indent, "dev <- dev + (-2) * ordered_logistic_log( ", 
                      out_var, "[i] , glm", var_suffix[f], "[i] , ", 
                      cutsname, " );\n", sep = "")
      }
      if (family[[f]] == "gamma") {
        thetaname <- paste("theta", var_suffix[f], sep = "")
        m_gq <- paste(m_gq, indent, indent, "dev <- dev + (-2) * gamma_log( ", 
                      out_var, "[i] , exp(glm", var_suffix[f], "[i])*", 
                      thetaname, " , ", thetaname, " );\n", sep = "")
      }
      m_gq <- paste(m_gq, indent, "}\n", sep = "")
    }
    m_gq <- paste(m_gq, "}\n ", sep = "")
    pars_list <- c(pars_list, "dev")
  }
  data_list <- list()
  for (f in 1:num_formulas) {
    index_name <- paste("N", var_suffix[f], sep = "")
    data_list[[index_name]] <- nrow(fp[[f]]$dat)
    outvar <- undot(fp[[f]]$yname)
    outvar <- paste(outvar, var_suffix[f], sep = "")
    data_list[[outvar]] <- fp[[f]]$y
    for (i in 1:ncol(fp[[f]]$dat)) {
      name <- undot(colnames(fp[[f]]$dat)[i])
      if (name != "Intercept") {
        col_name <- paste(name, var_suffix[f], sep = "")
        data_list[[col_name]] <- fp[[f]]$dat[, i]
      }
    }
    if (family[[f]] == "ordered") {
      Kname <- paste("K", var_suffix[f], sep = "")
      data_list[[Kname]] <- as.integer(max(fp[[f]]$y))
    }
  }
  if (length(cluster_vars) > 0) {
    for (i in 1:length(cluster_size)) {
      vname <- paste("N_", undot(names(cluster_size)[i]), 
                     sep = "")
      data_list[[vname]] <- cluster_size[[i]]
    }
    for (i in 1:length(cluster_vars)) {
      nterms <- sum(cluster_vars[[i]])
      if (nterms > 1) {
        name <- undot(names(cluster_vars)[i])
        Oname <- paste("Omega_", name, sep = "")
        data_list[[Oname]] <- diag(nterms)
      }
    }
  }
  init_list <- list()
  for (f in 1:num_formulas) {
    for (i in 1:length(fp[[f]]$fixef)) {
      varname <- undot(fp[[f]]$fixef[i])
      prefix <- fixed_prefix
      init_value <- 0
      if (varname == "Intercept") {
        prefix <- ""
        if (family[[f]] == "binomial") {
          init_value <- 0
        }
        if (family[[f]] == "gaussian") {
          init_value <- mean(fp[[f]]$y)
        }
        if (family[[f]] == "poisson") {
          init_value <- log(mean(fp[[f]]$y))
        }
        if (family[[f]] == "gamma") {
          init_value <- log(mean(fp[[f]]$y))
        }
      }
      parname <- paste(prefix, varname, var_suffix[f], 
                       sep = "")
      if (initmethod == "lme4") {
        init_value <- as.numeric(fixef(fit[[f]])[i])
      }
      init_list[[parname]] <- init_value
    }
  }
  if (length(cluster_vars) > 0) {
    for (i in 1:length(cluster_vars)) {
      name <- undot(names(cluster_vars)[i])
      vname <- paste(vary_prefix, name, sep = "")
      num_effects <- sum(cluster_vars[[i]])
      num_groups <- cluster_size[[i]]
      init_values <- matrix(0, nrow = num_groups, ncol = num_effects)
      if (num_effects == 1) 
        init_values <- as.numeric(init_values)
      name_varcov <- "sigma_"
      init_varcov <- 1
      if (num_effects > 1) {
        if (varpriors == "weak") {
          init_varcov <- rep(1, num_effects)
          name_rho <- paste("Rho_", name, sep = "")
          init_list[[name_rho]] <- diag(num_effects)
        }
        else {
          init_varcov <- diag(num_effects)
          name_varcov <- "Sigma_"
        }
      }
      name_varcov <- paste(name_varcov, name, sep = "")
      if (initmethod == "lme4") {
        if (num_effects == 1) {
        }
        else {
        }
      }
      init_list[[vname]] <- init_values
      init_list[[name_varcov]] <- init_varcov
    }
  }
  for (f in 1:num_formulas) {
    if (family[[f]] == "gaussian") {
      signame <- paste("sigma", var_suffix[f], sep = "")
      sig <- sd(fp[[f]]$y)
      if (initmethod == "lme4") 
        sig <- sigma(fit[[f]])
      if (sig < 0.001) 
        sig <- 1
      init_list[[signame]] <- sig
    }
    if (family[[f]] == "ordered") {
      cutsname <- paste("cutpoints", var_suffix[f], sep = "")
      K <- as.integer(max(fp[[f]]$y))
      init_list[[cutsname]] <- seq(from = -2, to = 1.5, 
                                   length.out = K - 1)
    }
    if (family[[f]] == "gamma") {
      thetaname <- paste("theta", var_suffix[f], sep = "")
      ftheta <- function(pars) {
        parmu <- mean(fp[[f]]$y)
        -sum(dgamma2(x = fp[[f]]$y, mu = parmu, scale = pars[1], 
                     log = TRUE))
      }
      o <- optim(c(1), ftheta, method = "SANN")
      init_list[[thetaname]] <- 1/o$par[1]
    }
  }
  m_data <- paste(m_data, "}\n", sep = "")
  if (flag_transdata == FALSE) {
    m_transdata <- ""
  }
  else {
    m_transdata <- paste(m_transdata1, m_transdata2, "}\n\n", 
                         sep = "")
  }
  if (flag_transpars == FALSE) {
    m_transpars <- ""
  }
  else {
    m_transpars <- paste(m_transpars1, m_transpars2, "}\n\n", 
                         sep = "")
  }
  model <- paste(m_data, "\n", m_transdata, m_pars, "\n", 
                 m_transpars, m_model, "\n", m_gq, sep = "")
  if (sample == FALSE) {
    result <- list(data = data_list, model = model, init = list(init_list), 
                   pars = pars_list[2:length(pars_list)], clusters = cluster_vars, 
                   clusters_size = cluster_size)
  }
  else {
    passinit <- initmethod
    if (initmethod == "zero" | initmethod == "lme4") {
      passinit <- list()
      for (i in 1:chains) {
        passinit[[i]] <- init_list
      }
    }
    pars_list <- pars_list[2:length(pars_list)]
    if (verbose == TRUE) {
      message("Starting Stan model")
      flush.console()
      start.time <- Sys.time()
    }
    if (num_formulas == 1) {
      modelname <- paste(deparse(formula[[1]]), collapse = "")
    }
    else {
      modelname <- paste(deparse(formula[[1]]), collapse = "")
    }
    modelname <- paste(modelname, " [", family[[1]], "]", 
                       sep = "")
    # my shit ####
    if(!is.null(mymodel)){model <- mymodel}
    result <- stan(model_name = modelname, model_code = model, 
                   data = data_list, init = passinit, iter = iter, 
                   warmup = warmup, chains = chains, pars = pars_list, 
                   save_dso = FALSE, ...)
    if (verbose == TRUE) {
      message(show(Sys.time() - start.time))
      flush.console()
    }
  }
  if ((calcDIC == TRUE | calcWAIC == TRUE) & sample == TRUE) {
    gc()
    post <- extract(result, permuted = TRUE)
    postbar <- list()
    for (i in 1:length(post)) {
      dims <- dim(post[[i]])
      if (length(dims) == 1) 
        postbar[[names(post)[i]]] <- mean(post[[i]])
      else postbar[[names(post)[i]]] <- apply(post[[i]], 
                                              2:length(dims), mean)
    }
    gc()
    if (calcDIC) {
      Dbar <- postbar$dev
      Dhat <- 0
    }
    if (verbose == TRUE) {
      if (calcWAIC) {
        message("Computing WAIC")
        message("Warning: WAIC only works for single-formula models right now.")
      }
      flush.console()
    }
    for (f in 1:num_formulas) {
      if (calcWAIC == TRUE) {
        if (family[[f]] == "binomial") {
          bintotvar <- paste("bin_total", var_suffix[f], 
                             sep = "")
          if (any(data_list[[bintotvar]] > 1)) {
            message("Warning: WAIC calculation not correct for aggregated binomial models. Recode data to logistic regression form (0/1 outcomes) instead.")
          }
        }
        pd <- 0
        lppd <- 0
        N_all <- nrow(fp[[f]]$dat)
        update_inc <- floor(N_all/10)
        for (i in 1:N_all) {
          if (floor(i/update_inc) == i/update_inc) {
            cat("\r")
            cat(paste("Progress: ", i, " / ", N_all, 
                      sep = ""))
          }
          N_samples <- length(post[[1]])
          vary <- rep(0, N_samples)
          if (length(fp[[f]]$ranef) > 0) {
            for (r in 1:length(fp[[f]]$ranef)) {
              cluster_name <- undot(names(fp[[f]]$ranef)[r])
              parname <- paste(vary_prefix, cluster_name, 
                               sep = "")
              nterms <- length(fp[[f]]$ranef[[r]])
              total_terms <- sum(cluster_vars[[cluster_name]])
              fnames <- fp[[f]]$ranef[[r]]
              for (j in 1:nterms) {
                jstart <- 0
                if (f > 1) {
                  jstart <- sum(cluster_vars[[cluster_name]][1:(f - 
                                                                  1)])
                }
                dox <- 1
                if (undot(fnames[j]) != "Intercept") {
                  xname <- paste(undot(fnames[j]), var_suffix[f], 
                                 sep = "")
                  dox <- data_list[[xname]][i]
                }
                cname <- paste(cluster_name, var_suffix[f], 
                               sep = "")
                if (total_terms > 1) {
                  vary <- vary + post[[parname]][, data_list[[cname]][i], 
                                                 j + jstart] * dox
                }
                else {
                  vary <- vary + post[[parname]][, data_list[[cname]][i]] * 
                    dox
                }
              }
            }
          }
          nterms <- length(fp[[f]]$fixef)
          fnames <- fp[[f]]$fixef
          glm <- rep(0, N_samples)
          for (j in 1:nterms) {
            prefix <- fixed_prefix
            if (undot(fnames[j]) == "Intercept") {
              if (family[[f]] != "ordered") {
                parname <- paste(undot(fnames[j]), var_suffix[f], 
                                 sep = "")
                glm <- glm + post[[parname]]
              }
            }
            else {
              parname <- paste(prefix, undot(fnames[j]), 
                               var_suffix[f], sep = "")
              xname <- paste(undot(fnames[j]), var_suffix[f], 
                             sep = "")
              dox <- data_list[[xname]][i]
              glm <- glm + post[[parname]] * dox
            }
          }
          outvar <- paste(undot(fp[[f]]$yname), var_suffix[f], 
                          sep = "")
          if (family[[f]] == "gaussian") {
            sigmavar <- paste("sigma", var_suffix[f], 
                              sep = "")
            ll <- dnorm(x = data_list[[outvar]][i], 
                        mean = glm + vary, sd = post[[sigmavar]], 
                        log = TRUE)
          }
          if (family[[f]] == "binomial") {
            bintotvar <- paste("bin_total", var_suffix[f], 
                               sep = "")
            ll <- dbinom(x = data_list[[outvar]][i], 
                         prob = logistic(glm + vary), size = data_list[[bintotvar]][i], 
                         log = TRUE)
          }
          if (family[[f]] == "poisson") {
            ll <- dpois(x = data_list[[outvar]][i], 
                        lambda = exp(glm + vary), log = TRUE)
          }
          if (family[[f]] == "ordered") {
            cutsname <- paste("cutpoints", var_suffix[f], 
                              sep = "")
            ll <- sapply(1:N_samples, function(s) dordlogit(x = data_list[[outvar]][i], 
                                                            phi = glm[s] + vary[s], a = post[[cutsname]][s, 
                                                                                                         ], log = TRUE))
          }
          if (family[[f]] == "gamma") {
            thetavar <- paste("theta", var_suffix[f], 
                              sep = "")
            ll <- dgamma2(x = data_list[[outvar]][i], 
                          mu = exp(glm + vary), scale = 1/post[[thetavar]], 
                          log = TRUE)
          }
          pd <- pd + var(ll)
          lppd <- lppd + log_sum_exp(ll) - log(length(ll))
        }
        cat("\r                             ")
        cat(paste("\rlppd =", round(lppd, 2)))
        cat(paste("\npWAIC =", round(pd, 2)))
        cat(paste("\nWAIC =", round(-2 * (lppd - pd), 
                                    2)))
        cat("\n")
        if (FALSE) {
          if (length(cluster_vars) > 0) {
            for (v in 1:length(cluster_vars)) {
              message(paste("Computing WAIC for cluster '", 
                            names(cluster_vars)[v], "'", sep = ""))
              pd_j <- 0
              lppd_j <- 0
              N_groups <- cluster_size[[v]]
              update_inc <- floor(N_groups/10)
              for (g in 1:N_groups) {
                if (floor(g/update_inc) == g/update_inc) {
                  cat("\r")
                  cat(paste("Progress: ", g, " / ", 
                            N_groups, sep = ""))
                }
                i_in_g <- which(data_list[[names(cluster_vars)[v]]] == 
                                  g)
                N_in_g <- length(i_in_g)
                N_samples <- length(post[[1]])
                vary <- matrix(0, nrow = N_in_g, ncol = N_samples)
                if (length(fp[[f]]$ranef) > 0) {
                  for (r in 1:length(fp[[f]]$ranef)) {
                    cluster_name <- undot(names(fp[[f]]$ranef)[r])
                    cname <- paste(cluster_name, var_suffix[f], 
                                   sep = "")
                    parname <- paste(vary_prefix, cluster_name, 
                                     sep = "")
                    nterms <- length(fp[[f]]$ranef[[r]])
                    total_terms <- sum(cluster_vars[[cluster_name]])
                    fnames <- fp[[f]]$ranef[[r]]
                    for (j in 1:nterms) {
                      jstart <- 0
                      if (f > 1) {
                        jstart <- sum(cluster_vars[[cluster_name]][1:(f - 
                                                                        1)])
                      }
                      dox <- rep(1, N_in_g)
                      if (undot(fnames[j]) != "Intercept") {
                        xname <- paste(undot(fnames[j]), 
                                       var_suffix[f], sep = "")
                        dox <- data_list[[xname]][i_in_g]
                      }
                      if (total_terms > 1) {
                        for (i in 1:N_in_g) {
                          vary[i, ] <- vary[i, ] + post[[parname]][, 
                                                                   g, j + jstart] * dox[i]
                        }
                      }
                      else {
                        for (i in 1:N_in_g) {
                          vary[i, ] <- vary[i, ] + post[[parname]][, 
                                                                   g] * dox[i]
                        }
                      }
                    }
                  }
                }
                glm <- matrix(0, nrow = N_in_g, ncol = N_samples)
                nterms <- length(fp[[f]]$fixef)
                fnames <- fp[[f]]$fixef
                for (j in 1:nterms) {
                  prefix <- fixed_prefix
                  if (undot(fnames[j]) == "Intercept") {
                    if (family[[f]] != "ordered") {
                      parname <- paste(undot(fnames[j]), 
                                       var_suffix[f], sep = "")
                      for (i in 1:N_in_g) glm[i, ] <- glm[i, 
                                                          ] + post[[parname]]
                    }
                  }
                  else {
                    parname <- paste(prefix, undot(fnames[j]), 
                                     var_suffix[f], sep = "")
                    xname <- paste(undot(fnames[j]), 
                                   var_suffix[f], sep = "")
                    for (i in 1:N_in_g) {
                      dox <- data_list[[xname]][i_in_g[i]]
                      glm[i, ] <- glm[i, ] + post[[parname]] * 
                        dox
                    }
                  }
                }
                outvar <- paste(undot(fp[[f]]$yname), 
                                var_suffix[f], sep = "")
                if (family[[f]] == "gaussian") {
                  sigmavar <- paste("sigma", var_suffix[f], 
                                    sep = "")
                  ll <- sapply(1:N_samples, function(s) sum(dnorm(x = data_list[[outvar]][i_in_g], 
                                                                  mean = glm[, s] + vary[, s], sd = post[[sigmavar]][s], 
                                                                  log = TRUE)))
                }
                if (family[[f]] == "binomial") {
                  bintotvar <- paste("bin_total", var_suffix[f], 
                                     sep = "")
                  ll <- sapply(1:N_samples, function(s) sum(dbinom(x = data_list[[outvar]][i_in_g], 
                                                                   prob = logistic(glm[, s] + vary[, 
                                                                                                   s]), size = data_list[[bintotvar]][i_in_g], 
                                                                   log = TRUE)))
                }
                if (family[[f]] == "poisson") {
                  ll <- sapply(1:N_samples, function(s) sum(dpois(x = data_list[[outvar]][i_in_g], 
                                                                  lambda = exp(glm[, s] + vary[, s]), 
                                                                  log = TRUE)))
                }
                if (family[[f]] == "ordered") {
                  cutsname <- paste("cutpoints", var_suffix[f], 
                                    sep = "")
                  ll <- sapply(1:N_samples, function(s) dordlogit(x = data_list[[outvar]][i_in_g], 
                                                                  phi = glm[, s] + vary[, s], a = post[[cutsname]][s, 
                                                                                                                   ], log = TRUE))
                }
                if (family[[f]] == "gamma") {
                  thetavar <- paste("theta", var_suffix[f], 
                                    sep = "")
                  ll <- sapply(1:N_samples, function(s) sum(dgamma2(x = data_list[[outvar]][i_in_g], 
                                                                    mu = exp(glm[, s] + vary[, s]), 
                                                                    scale = 1/post[[thetavar]][s], log = TRUE)))
                }
                pd_j <- pd_j + var(ll)
                lppd_j <- lppd_j + log_sum_exp(ll) - 
                  log(length(ll))
              }
              cat("\r                             ")
              message(paste("\rWAIC[", names(cluster_vars)[v], 
                            "]", sep = ""))
              cat(paste("lppd =", round(lppd_j, 2)))
              cat(paste("\npWAIC =", round(pd_j, 2)))
              cat(paste("\nWAIC[", names(cluster_vars)[v], 
                        "] =", round(-2 * (lppd_j - pd_j), 2)))
              cat("\n")
            }
          }
        }
      }
    }
    if (verbose == TRUE) {
      if (calcDIC) {
        if (calcWAIC) 
          cat("\n")
        message("Computing DIC")
      }
      flush.console()
    }
    for (f in 1:num_formulas) {
      if (calcDIC) {
        N_all <- nrow(fp[[f]]$dat)
        vary <- rep(0, N_all)
        if (length(fp[[f]]$ranef) > 0) {
          for (i in 1:length(fp[[f]]$ranef)) {
            cluster_name <- undot(names(fp[[f]]$ranef)[i])
            parname <- paste(vary_prefix, cluster_name, 
                             sep = "")
            nterms <- length(fp[[f]]$ranef[[i]])
            total_terms <- sum(cluster_vars[[cluster_name]])
            fnames <- fp[[f]]$ranef[[i]]
            for (j in 1:nterms) {
              jstart <- 0
              if (f > 1) {
                jstart <- sum(cluster_vars[[cluster_name]][1:(f - 
                                                                1)])
              }
              dox <- 1
              if (undot(fnames[j]) != "Intercept") {
                xname <- paste(undot(fnames[j]), var_suffix[f], 
                               sep = "")
                dox <- data_list[[xname]]
              }
              cname <- paste(cluster_name, var_suffix[f], 
                             sep = "")
              if (total_terms > 1) {
                vary <- vary + postbar[[parname]][data_list[[cname]], 
                                                  j + jstart] * dox
              }
              else {
                vary <- vary + postbar[[parname]][data_list[[cname]]] * 
                  dox
              }
            }
          }
        }
        nterms <- length(fp[[f]]$fixef)
        fnames <- fp[[f]]$fixef
        glm <- rep(0, N_all)
        for (j in 1:nterms) {
          prefix <- fixed_prefix
          if (undot(fnames[j]) == "Intercept") {
            if (family[[f]] != "ordered") {
              parname <- paste(undot(fnames[j]), var_suffix[f], 
                               sep = "")
              glm <- glm + postbar[[parname]]
            }
          }
          else {
            parname <- paste(prefix, undot(fnames[j]), 
                             var_suffix[f], sep = "")
            xname <- paste(undot(fnames[j]), var_suffix[f], 
                           sep = "")
            dox <- data_list[[xname]]
            glm <- glm + postbar[[parname]] * dox
          }
        }
        outvar <- paste(undot(fp[[f]]$yname), var_suffix[f], 
                        sep = "")
        if (family[[f]] == "gaussian") {
          sigmavar <- paste("sigma", var_suffix[f], 
                            sep = "")
          Dhat <- Dhat + (-2) * sum(dnorm(x = data_list[[outvar]], 
                                          mean = glm + vary, sd = postbar[[sigmavar]], 
                                          log = TRUE))
        }
        if (family[[f]] == "binomial") {
          bintotvar <- paste("bin_total", var_suffix[f], 
                             sep = "")
          Dhat <- Dhat + (-2) * sum(dbinom(x = data_list[[outvar]], 
                                           prob = logistic(glm + vary), size = data_list[[bintotvar]], 
                                           log = TRUE))
        }
        if (family[[f]] == "poisson") {
          Dhat <- Dhat + (-2) * sum(dpois(x = data_list[[outvar]], 
                                          lambda = exp(glm + vary), log = TRUE))
        }
        if (family[[f]] == "ordered") {
          Dhat <- Dhat + (-2) * sum(dordlogit(x = data_list[[outvar]], 
                                              phi = glm + vary, a = postbar[["cutpoints"]], 
                                              log = TRUE))
        }
        if (family[[f]] == "gamma") {
          thetavar <- paste("theta", var_suffix[f], 
                            sep = "")
          Dhat <- Dhat + (-2) * sum(dgamma2(x = data_list[[outvar]], 
                                            mu = exp(glm + vary), scale = 1/postbar[[thetavar]], 
                                            log = TRUE))
        }
      }
    }
    pD <- Dbar - Dhat
    DIC <- Dbar + pD
    cat(paste("Expected deviance =", round(Dbar, 2)))
    cat(paste("\nDeviance of expectation =", round(Dhat, 
                                                   2)))
    cat(paste("\npDIC =", round(pD, 2)))
    cat(paste("\nDIC =", round(DIC, 2)))
    cat("\n")
  }
  if (extract == TRUE & sample == TRUE) {
    gc()
    result <- extract(result, permuted = TRUE)
  }
  if (calcDIC == TRUE & sample == TRUE) {
    attr(result, "DIC") <- list(DIC = DIC, pD = pD, Dhat = Dhat, 
                                Dbar = Dbar)
  }
  if (calcWAIC == TRUE & sample == TRUE) {
    attr(result, "WAIC") <- list(pD = pd, lppd = lppd, WAIC = -2 * 
                                   (lppd - pd))
  }
  calcgini <- function(x, weights = rep(1, length = length(x))) {
    ox <- order(x)
    x <- x[ox]
    weights <- weights[ox]/sum(weights)
    p <- cumsum(weights)
    nu <- cumsum(weights * x)
    n <- length(nu)
    nu <- nu/nu[n]
    sum(nu[-1] * p[-n]) - sum(nu[-n] * p[-1])
  }
  ranefattr <- list()
  if (length(cluster_vars) > 0) {
    for (i in 1:length(cluster_vars)) {
      name <- names(cluster_vars)[i]
      n <- cluster_size[[name]]
      data2 <- as.data.frame(data)
      udCols <- sapply(colnames(data2), undot)
      iOrig <- which(udCols == name)
      n_per_j <- sapply(1:n, function(j) sum(data2[, iOrig] == 
                                               j))
      gini <- calcgini(n_per_j)
      factors <- c()
      for (f in 1:num_formulas) {
        if (!is.null(fp[[f]]$ranef[[name]])) 
          factors <- c(factors, fp[[f]]$ranef[[name]])
      }
      ranefattr[[i]] <- list(factors = factors, n = n, 
                             gini = gini)
      names(ranefattr)[i] <- name
    }
  }
  attr(result, "ranef") <- ranefattr
  attr(result, "formulas") <- list(fp = fp, cluster_vars = cluster_vars, 
                                   cluster_size = cluster_size, var_suffix = var_suffix, 
                                   family = family, formula = formula)
  result
}
