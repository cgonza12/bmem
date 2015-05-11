
parse_formula2<- function (formula, data) 
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
