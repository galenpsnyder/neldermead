## optimization
# some test functions:
# Rosenbrock banana function--global minimum at (1, 1)
rosenbrock <- function(x){
  y <- x[1]
  x <- x[2]
  (1 - x)^2 + 100 * (y - x^2)^2
}

# Beale function--global minimum at (3, 0.5)
beale <- function(x){
  y <- x[2]
  x <- x[1]
  (1.5 - x + x*y)^2 + (2.25 - x + x*y^2)^2 + (2.625 - x + x*y^3)^2
}

nelder.mead <- function(parm, objective, lower = NULL, upper = NULL, tol = sqrt(.Machine$double.eps), hessian = FALSE, ...){
  # initializing important values
  parm.vec <- as.vector(parm)                        # vector of user-provided parameters
  n        <- length(parm.vec)                       # number of parameters--defines dimension of objective function
  g        <- 1 + (2 / n )                           # expansion parameter
  b        <- 0.75 - (1 / (2 * n))                   # contraction parameter
  d        <- 1 - (1 / n)                            # shrink parameter
  s        <- Inf                                    # initial parameter for convergense checking
  tol      <- tol                                    # point at which procedure converges
  iter     <- 0                                      # number of iterations for output
  parm     <- matrix(0, nrow = n + 1, ncol = n + 1)  # parameter matrix--to be filled with objective function evaluations
  
  ifelse(is.null(lower), lower <- rep(-Inf, n), lower <- lower)
  ifelse(is.null(upper), upper <- rep(Inf,  n), upper <- upper)
  
  for(i in 1:n) if(parm.vec[i] < lower[i]) stop("Start values outside contraints!")
  
  for(i in 1:n) if(parm.vec[i] > upper[i]) stop("Start values outside contraints!")
  
  parm[, 1:n] <- matrix(rep(parm.vec, n + 1),        # all rows of parameter matrix gets initial parameters
                        nrow =  n + 1, byrow = T)    # note that the last column is left empty for objective function evaluations
  
  for(i in 2:(n + 1)){
    if(parm[1, (i - 1)] != 0){
      parm[i, (i - 1)] <- parm[1, (i - 1)] + (parm[1, (i - 1)] * 0.05)
    } else {
      parm[i, (i - 1)] <- 0.0075
    }
  }
  
  for(i in 1:(n + 1)){
    parm[i, (n + 1)] <- objective(parm[i, 1:n], ...)
  }
  parm <- parm[order(parm[ ,(n + 1)]), ]
  
  while(s > tol){
    # estimating centroid of vertices less the worst fitting vertex
    x.m <- vector(mode = "numeric", length = n)
    for(j in 1:n){
      x.m.init <- c()
      for(i in 1:n){
        x.m.init <- c(x.m.init, parm[i, j])
      }
      x.m[[j]] <- sum(x.m.init) / n
    }
    
    # calculating reflection of the centroid
    x.r <- x.m + (x.m - parm[(n + 1), 1:n])
    for(i in 1:n) ifelse(x.r[i] < lower[i], x.r[i] <- lower[i], x.r[i] <- x.r[i])
    for(i in 1:n) ifelse(x.r[i] > upper[i], x.r[i] <- upper[i], x.r[i] <- x.r[i])
    x.r.eval <- objective(x.r, ...)
    
    if(x.r.eval < parm[n, (n + 1)]){
      if(x.r.eval < parm[1, (n + 1)]){
        x.e <- x.r + g * (x.r - x.m)
        for(i in 1:n) ifelse(x.e[i] < lower[i], x.e[i] <- lower[i], x.e[i] <- x.e[i])
        for(i in 1:n) ifelse(x.e[i] > upper[i], x.e[i] <- upper[i], x.e[i] <- x.e[i])
        x.e.eval <- objective(x.e, ...)
        if(x.e.eval < x.r.eval){
          parm[(n + 1), 1:n] <- x.e
          parm[(n + 1), (n + 1)] <- x.e.eval
        } else {
          parm[(n + 1), 1:n] <- x.r
          parm[(n + 1), (n + 1)] <- x.r.eval
        }
      } else {
        parm[(n + 1), 1:n] <- x.r
        parm[(n + 1), (n + 1)] <- x.r.eval
      }
    } else {
      if(x.r.eval >= parm[(n + 1), (n + 1)]){
        x.c <- x.m + b * (parm[(n + 1), 1:n] - x.m)
        for(i in 1:n) ifelse(x.c[i] < lower[i], x.c[i] <- lower[i], x.c[i] <- x.c[i])
        for(i in 1:n) ifelse(x.c[i] > upper[i], x.c[i] <- upper[i], x.c[i] <- x.c[i])
        x.c.eval <- objective(x.c, ...)
        if(x.c.eval < parm[(n + 1), (n + 1)]){
          parm[(n + 1), 1:n] <- x.c
          parm[(n + 1), (n + 1)] <- x.c.eval
        } else {
          for(i in 2:(n + 1)){
            parm[i, 1:n] <- parm[1, 1:n] + d *(parm[i, 1:n] - parm[1, 1:n])
            parm[i, (n + 1)] <- objective(parm[i, 1:n], ...)
          }
        }
      } else {
        x.c <- x.m + b * (x.r - x.m)
        for(i in 1:n) ifelse(x.c[i] < lower[i], x.c[i] <- lower[i], x.c[i] <- x.c[i])
        for(i in 1:n) ifelse(x.c[i] > upper[i], x.c[i] <- upper[i], x.c[i] <- x.c[i])
        x.c.eval <- objective(x.c, ...)
        if(x.c.eval < x.r.eval){
          parm[(n + 1), 1:n] <- x.c
          parm[(n + 1), (n + 1)] <- x.c.eval
        } else {
          for(i in 2:(n + 1)){
            parm[i, 1:n] <- parm[1, 1:n] + d *(parm[i, 1:n] - parm[1, 1:n])
            parm[i, (n + 1)] <- objective(parm[i, 1:n], ...)
          }
        }
      }
    }
    parm <- parm[order(parm[ ,(n + 1)]), ]
    iter <- iter + 1
    s <- 2 * ((parm[(n + 1), (n + 1)] - parm[1, (n + 1)]) / (parm[(n + 1), (n + 1)] + parm[1, (n + 1)] + tol))
  }
  if(hessian){
    ## estimating hessian via quadratic surface approximation (e.g., Nelder & Mead, 1965)
    # create empty matrix 'H' and fill columns with values from best fitting simplex node
    # then coerce jitter between points
    H <- matrix(0, nrow = n+1, ncol = n+1)
    H[, 1:n] <- matrix(rep(parm[1, 1:n], n+1), nrow =  n+1, byrow = T)
    for(i in 2:(n+1)){
      if(H[1, (i - 1)] != 0){
        H[i, (i - 1)] <- H[1, (i - 1)] + (H[1, (i - 1)] * 0.01)
      } else {
        H[i, (i - 1)] <- 0.0075
      }
    }
    
    # evaluate objective function at newly estimated points and order from best to worst
    for(i in 1:(n+1)){
      H[i, n+1] <- objective(H[i, 1:n], ...)
    }
    H <- H[order(H[ ,(n+1)]), ]
    
    # place function evaluation points in vector 'y_i'
    # create row and column indices and order them by row values
    y_i <- H[, ncol(H)]
    k <- length(y_i)
    p <- sequence(1:(k-1))
    q <- rep(2:k, 1:(k-1))
    lind <- matrix(c(p, q), nrow = length(p), ncol = 2)
    if(length(p) > 1) lind <- lind[order(lind[, 1]), ]
    p <- lind[, 1]
    q <- lind[, 2]
    h <- length(p)
    
    # create empty matrix and fill each row with the midpoints of the simplex nodes generated above
    # then evaluate objective function at midpoints
    # finally, create full matrix whose i,jth value is the value of the objective function
    #    evaluated at the midpoint between simplex node i, j
    mids <- matrix(0, nrow = h, ncol = k)
    for(i in 1:h){
      mids[i, ] <- (H[p[i], ] + H[q[i], ])/2
    }
    mids[, k] <- sapply(1:h, function(i) objective(mids[i, -k], ...))
    mids_new <- matrix(0, nrow = k, ncol = k)
    for(i in 1:h){
      mids_new[p[i], q[i]] <- mids[i, k]
    }
    mids <- mids_new
    t.mids <- t(mids)
    diag(t.mids) <- 0
    mids <- mids + t.mids
    
    
    b_i <- matrix(0, k, k)
    for(i in 2:k){
      for(j in 2:k){
        if(i == j){
          b_i[i, j] <- 2*(y_i[i] + y_i[1] - 2*mids[1, i])
        } else {
          b_i[i, j] <- 2*(mids[i, j] + y_i[1] - mids[1, i] - mids[1, j])
        }
      }
    }
    b_i <- b_i[-1, -1]
    
    Q <- matrix(0, k-1, k-1)
    for(i in 2:k){
      Q[, i-1] <- H[i, -k] - H[1, -k]
    }
    
    hessian <- 2*t(solve(Q)) %*% b_i %*% solve(Q)
    
    out <- list(parm = parm[1, 1:n],
                value = parm[1, (n + 1)],
                iterations = iter,
                hessian = hessian)
    
  } else {
    out <- list(parm = parm[1, 1:n],
                value = parm[1, (n + 1)],
                iterations = iter)
  }
  out
}

(nelder.mead(par = c(1, 1), rosenbrock, upper = c(2, Inf)))
(nelder.mead(par = c(0, 0), beale))
