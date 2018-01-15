## From Yanhao Shi's thesis -- R/Functions/CLA_analysis.R

findSig <- function(Mu0, result, covar, equal.tol){ # equal.tol > 1e-16
  ms.w <- result$MS_weight
  nd <- order(ms.w[, "Mu"])
  n <- length(nd)

  ## revert order for all (w, mu, sig):
  ms.w <- ms.w[nd, ]
  w <- result$weights_set[, nd]
  ind.uni <- c(TRUE, sapply(1:(n-1), function(i) !
                              isTRUE(all.equal(ms.w[i, ], ms.w[i+1, ], tol = equal.tol))))

  p <- rep(1, n)  # pointer, different weights of same ms have same pointer
  k <- 1
  for(i in 1:n){
    if(ind.uni[i]){
      p[i] <- k
      k <- k + 1
    }else{
      p[i] <- p[i-1]
    }
  }

  mu.w <- ms.w[ind.uni, "Mu"]
  sig.w <- ms.w[ind.uni, "Sig"]
  if(Mu0 < min(mu.w) || Mu0 > max(mu.w)) stop(sprintf("Mu0 must be in [%g, %g]",
                                                     max(mu.w), min(mu.w)))
  i <- findInterval(Mu0, mu.w) # [..)...[..]
  if(i == n) { # Mu0 is max(mu.w)
    list(Sig = sig.w[i], weight = w[, i])
  } else{
    m1 <- which(p == i)
    m2 <- which(p == i+1)
    a <- (Mu0 - mu.w[i+1])/(mu.w[i] - mu.w[i+1])
    # solve for a : Mu0 = a* Mu1 + (1-a)* Mu2 :
    wm <- matrix(0, n, m1*m2)
    for(s in m1){
        wm[,(s-1)*m2 + 1:m2] <- a*w[, s] + (1-a)*w[, m2]
    }
    list(Sig = sqrt(colSums(wm *(covar %*% wm))), #return Sig0
         weight = wm) # might have replicate columns
  }
}


##' A version of findSig() only used inside  findMu():
findSig0 <- function(Mu0, result, covar, equal.tol = 1e-8){ # equal.tol > 1e-16
  ms.w <- result$MS_weight
  nd <- order(ms.w[, "Mu"])
  n <- length(nd)

  ## revert order for all (w, mu, sig):
  ms.w <- ms.w[nd, ]
  w <- result$weights_set[, nd]

   mu.w <- ms.w[ , "Mu"]
  sig.w <- ms.w[ , "Sig"]
  i <- findInterval(Mu0, mu.w) # [..)...[..]
  if(i == n){ # Mu0 is max(mu.w)
    sig.w[i]
  } else if(isTRUE(all.equal(mu.w[i], mu.w[i+1], tol = equal.tol))){
    sig.w[i]
  }else{
    a <- (Mu0 - mu.w[i+1])/(mu.w[i] - mu.w[i+1])
    # solve for a : Mu0 = a* Mu1 + (1-a)* Mu2 :
    w0 <- a*w[, i] + (1-a)*w[, i+1]
    sqrt(colSums(w0 *(covar %*% w0) )) #return Sig0
  }
}

findMu <- function(Sig0, result, covar, tol.unir = 1e-6, equal.tol = 1e-6){
  ms.w <- result$MS_weight
  nd <- order(ms.w[, "Sig"])
  n <- length(nd)
   mu.w <- ms.w[nd, "Mu"]
  sig.w <- ms.w[nd, "Sig"]
  w <- result$weights_set[, nd]
  if(Sig0 < min(sig.w) || Sig0 > max(sig.w)) stop(sprintf("Sig0 must be in [%g, %g]",
                                                          max(sig.w), min(sig.w)))
  i <- findInterval(Sig0, sig.w)
  if(i == n){ # Mu0 is max(mu.w)
    list(Mu = mu.w[i], weight = w[, i])
  }else if(isTRUE(all.equal(mu.w[i], mu.w[i+1], tol = equal.tol))){
    list(Mu = mu.w[i], weight = w[, i])
  } else{
    r <- uniroot(function(mu) findSig0(mu, result, covar) - Sig0,
                 interval = mu.w[c(i, i+1)], tol=tol.unir)
    Mu0 <- r$root
    # solve for a : mu = a* Mu1 + (1-a)* Mu2 :
    a <- (Mu0 - mu.w[i+1])/(mu.w[i] - mu.w[i+1])
    w0 <- a* w[, i] + (1-a)*w[, i+1]
    list(Mu = Mu0, weight = w0)
  }
}










