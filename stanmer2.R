
stanmer2 = function (fit, digits = 4, probs = c(0.025, 0.975), fixed_prefix = "beta_", 
          vary_prefix = "vary_") 
{
  
digits = digits
fixed_prefix = "beta_"
vary_prefix = "vary_"

  spot <- function(pattern, target) length(grep(pattern, target, 
                                                fixed = TRUE)) > 0
  unprefix <- function(astring) {
    astring <- gsub(fixed_prefix, "", astring, fixed = TRUE)
    astring <- gsub(vary_prefix, "", astring, fixed = TRUE)
    astring <- gsub("sigma_", "", astring, fixed = TRUE)
    astring <- gsub("Sigma_", "", astring, fixed = TRUE)
    astring <- gsub("Intercept", "(Intercept)", astring, 
                    fixed = TRUE)
    astring
  }
  unprefix2 <- function(astring) {
    astring <- gsub(fixed_prefix, "", astring, fixed = TRUE)
    astring <- gsub("Intercept", "(Intercept)", astring, 
                    fixed = TRUE)
    astring
  }
  
  post <- extract(fit, permuted = TRUE)
  post.mu <- list()
  post.se <- list()
  for (i in 1:length(post)) {
    dims <- length(dim(post[[i]]))
    name <- names(post)[i]
    if (name != "lp__") {
      if (dims == 1) {
        post.mu[[name]] <- mean(post[[i]])
        post.se[[name]] <- sd(post[[i]])
      }
      else {
        post.mu[[name]] <- apply(post[[i]], 2:dims, 
                                 mean)
        post.se[[name]] <- apply(post[[i]], 2:dims, 
                                 sd)
      }
    }
  }
  fixlist <- c()
  ranlist <- c()
  for (i in 1:length(post.mu)) {
    name <- names(post.mu)[i]
    if (spot(vary_prefix, name)) {
      ranlist <- c(ranlist, name)
    }
    else {
      vname <- paste(vary_prefix, unprefix(name), sep = "")
      if (is.null(post.mu[[vname]]) & !spot("cutpoints", 
                                            name) & name != "dev") {
        fixlist <- c(fixlist, name)
      }
    }
  }
  fixlist = fixlist[!fixlist %in% c('vary','glm')]
  if (length(fixlist) > 0) {
    fix <- data.frame(Expectation = as.numeric(post.mu[fixlist]), 
                      StdDev = as.numeric(post.se[fixlist]))
    rownames(fix) <- sapply(fixlist, unprefix2)
    if (!is.null(probs)) {
      q <- matrix(0, nrow = nrow(fix), ncol = length(probs))
      for (k in 1:length(fixlist)) {
        q[k, ] <- quantile(post[[fixlist[k]]], probs)
      }
      colnames(q) <- paste(round(probs * 100, 1), "%", 
                           sep = "")
      fix <- cbind(fix, q)
    }
    cat(paste("glmer2stan model: ", fit@stanmodel@model_name, 
              "\n\n", sep = ""))
    cat("Level 1 estimates:\n")
    print(round(fix, digits))
  
    return(round(fix, digits))
  }

if (!is.null(attr(fit, "DIC"))) {
  w <- attr(fit, "DIC")
  cat(paste("\nDIC: ", round(w$DIC, 0), "   pDIC: ", round(w$pD, 
                                                           1), "   Deviance: ", round(w$Dhat, 1), "\n", sep = ""))

if (!is.null(attr(fit, "WAIC"))) {
  w <- attr(fit, "WAIC")
  cat(paste("\nWAIC: ", round(w$WAIC, 0), "   pWAIC: ", 
            round(w$pD, 1), "   -2*lppd: ", round(-2 * w$lppd, 
                                                  1), "\n", sep = ""))
}
print(w)
}
}
