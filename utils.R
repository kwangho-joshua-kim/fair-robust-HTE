# The following function implements the fair and robust CATE estimator 
# proposed in Kim and Zubizarreta (2023) (https://arxiv.org/abs/2306.03625)
# using the power series expansion and sample splitting with K=2 splits

fr_cate <- function(dat, 
                    nm.s.features,
                    nm.l.factors, 
                    solver="mosek",
                    nuisance.est='rf',
                    ps.trim = "Sturmer.1",
                    delta=0.005, sqr=TRUE, interactions=TRUE) {
  
  #################################################
  # dat: n x (p_x+p_s+2) data frame for (X,S,A,Y)
  #   - A: binary treatment
  #   - Y: continuous outcome 
  #   - S: p_s sensitive features (S1, S2,...,S_{p_s}) 
  #   - X: p_x pre-treatment covariates (X1, X2,...,X_{p_x})
  # nm.s.features: a vector of characters, 
  #                containing the names of the sensitive features:
  #                e.g., c("S1", "S2")  
  # nm.l.factors: a vector of characters, 
  #               containing the names of the legitimate factors
  #               to be used for conditional parity (NA if not)
  #               e.g., c(NA, "L2") 
  # delta: tolerance level (a single value)
  # sqr: add the 2nd moments if TRUE
  # interactions: add the interaction terms if TRUE
  # solver: optimization solver \in {mosek, quadprog}
  # nuisance.est: modeling technique for the nuisance estimation
  #   - only Random Forests ("rf") are available for now
  # ps.trim: removes the individuals with extreme weights.
  #   - "Sturmer.1": the common range method by Stürmer (lower cutpoint = lowest PS in the exposed; upper cutpoint = highest PS in the unexposed)
  #   - "Sturmer.2": the original Stürmer method (lower cutpoint, 5th PS percentile in the exposed; upper cutpoint, 95th PS percentile in the unexposed)
  #   - "no.trim": proceed without PS trimming
  #################################################

  
  # construct basis functions
  b <- dat[,!(colnames(dat) %in% c("A", "Y"))]
  if (sqr) { # add the 2nd moments
    b.sqr <- b^2
    colnames(b.sqr) <- paste(colnames(b),".sqr",sep = "")
    if (interactions) { # add the 1nd-order interactions
      indx <- combn(colnames(b),2)
      int <- as.data.frame(do.call(cbind,
                                   lapply(split(indx, col(indx)), 
                                          function(x) b[,x[1]]*b[,x[2]])))
      colnames(int) <- apply(indx, 2, function(x) paste(x[1],".",x[2],sep=""))
      b.sqr <- cbind(b.sqr, int)
    }
    b.mat <- as.matrix(cbind(b, b.sqr))
  } else {
    b.mat <- as.matrix(b)
  }

  # counterfactual component estimation based on two-fold sample splitting 
  n <- nrow(dat)
  fold_ind <- sample.int(n = n, size = floor(.5*n), replace = F)
  phi.b.hat <- matrix(NA,nrow=2, ncol=ncol(b.mat))
  for (k in 1:2) {
    if (k==1) {
      dat.train <- dat[fold_ind, ]; dat.test <- dat[-fold_ind, ];
      b.test <- b.mat[-fold_ind,]
    } else {
      dat.train <- dat[-fold_ind, ]; dat.test <- dat[fold_ind, ];
      b.test <- b.mat[fold_ind,]
    }
    # nuisance function estimation
    if (nuisance.est=="rf") {
      f.pi <- paste("A ~", paste(colnames(dat.train)[!(colnames(dat) %in% c("A", "Y"))], collapse=" + "))
      f.mu <- paste("Y ~", paste(colnames(dat.train)[!(colnames(dat) %in% c("Y"))], collapse=" + "))
      pi.rf <- ranger(f.pi, dat=dat.train, write.forest = TRUE, verbose = FALSE) 
      mu.rf <- ranger(f.mu, dat=dat.train, write.forest = TRUE, verbose = FALSE) 
      pi.hat <- predict(pi.rf, data=dat.test)$predictions
      mu1.hat <- predict(mu.rf, data=data.frame(dat.test[,!(colnames(dat) %in% c("A", "Y"))],A=1))$predictions
      mu0.hat <- predict(mu.rf, data=data.frame(dat.test[,!(colnames(dat) %in% c("A", "Y"))],A=0))$predictions
      muA.hat <- dat.test$A*mu1.hat + (1-dat.test$A)*mu0.hat
      # PS trimming - removing the individuals with extreme PS values
      idx.exclude <- NULL
      if (ps.trim == "Sturmer.1") {
        # the common range method ver.1 by Stürmer et al. Am J Epidemiol 2021;190:1659–1670.
        idx.exclude <- c(which((pi.hat < min(pi.hat[dat.test$A==1])) == TRUE),
                         which((max(pi.hat[dat.test$A==0]) < pi.hat) == TRUE)
        )
      }
      if (ps.trim == "Sturmer.2") {
        # the common range method ver.2 by Stürmer et al. Am J Epidemiol 2010;172:843–854.
        idx.exclude <- c(which((pi.hat < quantile(pi.hat[dat.test$A==1], probs = 0.05)) == TRUE),
                         which((quantile(pi.hat[dat.test$A==0], probs = 0.95) < pi.hat) == TRUE)
        )
      }
    } else {
      stop("other methods TBD")
    }
    if (is.null(idx.exclude)) {
      phi.b.hat[k,] <- colMeans((dat.test$A/pi.hat * (dat.test$Y - muA.hat) + mu1.hat 
                                 - (1-dat.test$A)/(1-pi.hat) * (dat.test$Y - muA.hat) - mu0.hat) * b.test)
    } else {
      phi.diff <- (dat.test$A/pi.hat * (dat.test$Y - muA.hat) + mu1.hat 
                   - (1-dat.test$A)/(1-pi.hat) * (dat.test$Y - muA.hat) - mu0.hat)
      phi.b.hat[k,] <- colMeans((phi.diff * b.test)[-idx.exclude, ])
    }
  }
  phi.b.hat <- colMeans(phi.b.hat)

  # optimization
  Q.mat <- (t(b.mat) %*% b.mat)/nrow(b.mat)
  C.mat <- NULL
  if (solver=="quadprog") {
    for (s.nm in nm.s.features) {
      temp.C <- rbind((1-dat[,s.nm]) %*% b.mat/(n*mean(1 - dat[,s.nm])) - dat[,s.nm] %*% b.mat/(n*mean(dat[,s.nm])),
                      -(1-dat[,s.nm]) %*% b.mat/(n*mean(1 - dat[,s.nm])) + dat[,s.nm] %*% b.mat/(n*mean(dat[,s.nm])))
      C.mat <- rbind(C.mat, temp.C)
    }
    C.mat <- Matrix(C.mat)
    d.vec <- -rep(delta,2*length(nm.s.features))
    model <- solve.QP(Dmat=Q.mat,dvec=phi.b.hat,Amat=t(C.mat),bvec=d.vec)
    beta.hat <- model$solution
  }

  if (solver=="mosek") {
    k <- ncol(b.mat)
    n.stat.parity <- sum(is.na(nm.l.factors))
    n.cond.stat.parity <- 2 * (length(nm.l.factors) - n.stat.parity)
    n.constraints <- n.stat.parity + n.cond.stat.parity
    prob <- list(sense="min")
    prob$c <- c(-phi.b.hat, 1)
    prob$bx <- rbind(blx=c(rep(-Inf,k),0), bux=c(rep(Inf,k+1)))
    for (j in 1:length(nm.s.features)) {
      s.nm <- nm.s.features[j]
      if (is.na(nm.l.factors[j])) {
        temp.C <- cbind((1-dat[,s.nm]) %*% b.mat/(n*mean(1 - dat[,s.nm])) - dat[,s.nm] %*% b.mat/(n*mean(dat[,s.nm])),0)
      } else {
        l.nm <- nm.l.factors[j]
        temp.C <- rbind(cbind((dat[,l.nm]*(1-dat[,s.nm])) %*% b.mat/(n*mean(dat[,l.nm]*(1-dat[,s.nm]))) 
                                - (dat[,l.nm]*dat[,s.nm]) %*% b.mat/(n*mean(dat[,l.nm]*dat[,s.nm])), 0),
                        cbind(((1-dat[,l.nm])*(1-dat[,s.nm])) %*% b.mat/(n*mean((1-dat[,l.nm])*(1-dat[,s.nm]))) 
                              - ((1-dat[,l.nm])*dat[,s.nm]) %*% b.mat/(n*mean((1-dat[,l.nm])*dat[,s.nm])), 0)
        )
      }
      C.mat <- rbind(C.mat, temp.C)
    }
    prob$A <- Matrix(C.mat)
    prob$bc <- rbind(blc=rep(-delta, n.constraints),
                     buc=rep(delta, n.constraints))
    # Cholesky factorization
    Q <- suppressWarnings(chol(Q.mat, TRUE))
    r <- attr(Q, 'rank')
    # if (r < nrow(x)) Q[(r+1):nrow(x), (r+1):nrow(x)] <- 0
    Q[-(1:r), -(1:r)] <- 0
    oo <- order(attr(Q, 'pivot'))
    unpivQ.mat <- Q[, oo]
    # all.equal(crossprod(unpivQ.mat), Q.mat)
    prob$F <- Matrix(rbind(c(rep(0,k), 1),
                           rep(0,k+1),
                           cbind(unpivQ.mat, as(matrix(0,k,1), "dgCMatrix"))
    ), sparse = TRUE)
    prob$g <- c(0, 1, rep(0,k))
    prob$cones <- matrix(list("RQUAD", k+2, NULL), nrow=3, ncol=1)
    rownames(prob$cones) <- c("type","dim","conepar")
    r <- mosek(prob, list(verbose=0))
    beta.hat <- r$sol$itr$xx[1:k]
  }
  
  return(list(b.mat=b.mat, beta.hat=beta.hat))
}
