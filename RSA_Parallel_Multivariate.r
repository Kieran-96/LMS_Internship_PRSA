euclid.dist <- function(X1,X2){
  #
  #  Purpose:
  #     Calculate the euclidean distance between two vectors
  #
  #  Input:
  #     X1 - The first vector
  #     X2 - A second vector
  #
  #  Output:
  #     The euclidean distance between the two vectors
  #
  output <- .Call("distMatr", x = t(X1), y = t(X2))
  return(output)
}

calcKR <- function(S1, S2, phi, type = c("Exp", "MatThreeTwo")){
  #
  #  Purpose:
  #     Call some loaded C code to calculate the covariance matrix between two sets of coordinates
  #
  #  Input:
  #     S1 - The first set of coordinates
  #     S2 - A second set of coordinates
  #     phi - A parameter for the covariance function
  #     type - a string indicating which of the available covariance functions should be used
  #
  #  Output:
  #     The covariance matrix
  #
  type <- match.arg(type)
  out <- .Call(paste("calcK", type, "r", sep = ""), x = t(S1), y = t(S2), phi = phi)
  return(out)
}

calcKDerivR <- function(S1, S2, phi, type = c("Exp", "MatThreeTwo")){
  #
  #  Purpose:
  #     Call some loaded C code to calculate the derivative of the covariance matrix
  #     between two sets of coordinates with respect to phi
  #
  #  Input:
  #     S1 - The first set of coordinates
  #     S2 - A second set of coordinates
  #     phi - A parameter for the covariance function
  #     type - a string indicating which of the available covariance functions should be used
  #
  #  Output:
  #     The derivative of the covariance matrix
  #
  type <- match.arg(type)
  out <- .Call(paste("calcK", type, "Derivr", sep = ""), x = t(S1), y = t(S2), phi = phi)
  return(out)
}

H.B <- function(C, Z, M, A.inv, Sigma.inv){
  #
  #  Purpose:
  #     Calculate the value of the updating equation for the mean function parameters
  #
  #  Input:
  #     C - A matrix containing covariate values
  #     Z - A matrix containing response values
  #     M - The current values from the mean function
  #     A.inv - The inverse of the current covariance matrix
  #     Sigma.inv - The inverse of the current covariance matrix for the multiple responses
  #
  #  Output:
  #     The value of the updating equation
  #
  output <- t(C) %*% A.inv %*% (Z - M) %*% Sigma.inv
  return(output)
}

H.Sigma <- function(Sigma.inv, mat1, Z, M, m){
  #
  #  Purpose:
  #     Calculate the value of the updating equation for the covariance matrix(response)
  #
  #  Input:
  #     Sigma.inv - The inverse of the current covariance matrix(response)
  #     mat1 - The current value of the expression t((z - M)) %*% A.inv
  #     Z - A matrix containing response values
  #     M - The current values from the mean function
  #     m - The size of the subsample being used
  #
  #  Output:
  #     The value of the updating equation
  #
  partA <- 0.5 * Sigma.inv %*% mat1 %*% (Z - M) %*% Sigma.inv
  partB <- (-m/2) * Sigma.inv
  output <- partA + partB
  return(output)
}

H.U <- function(U, Sigma.inv, mat1, Z, M, m){
  #
  #  Purpose:
  #     Calculate the value of the updating equation for the covariance matrix(response)
  #
  #  Input:
  #     U - The factor of the decomposition of Sigma
  #     Sigma.inv - The inverse of the current covariance matrix(response)
  #     mat1 - The current value of the expression t((z - M)) %*% A.inv
  #     Z - A matrix containing response values
  #     M - The current values from the mean function
  #     m - The size of the subsample being used
  #
  #  Output:
  #     The value of the updating equation
  #
  partA <- U %*% Sigma.inv %*% mat1 %*% (Z - M) %*% Sigma.inv
  partB <- -m * U %*% Sigma.inv
  output <- partA + partB
  return(output)
}

H.phi <- function(A.inv, dAdphi, Sigma.inv, mat1, mat2, q){
  #
  #  Purpose:
  #     Calculate the value of the updating equation for the covariance function parameter phi
  #
  #  Input:
  #     A.inv - The inverse of the current covariance matrix
  #     dAdphi - The derivative of the covariance matrix
  #     Sigma.inv - The inverse of the current covariance matrix(response)
  #     mat1 - The current value of the expression t((z - M)) %*% A.inv
  #     mat2 - The current value of the expression A.inv %*% (z - M)
  #     q - The number of response variables
  #
  #  Output:
  #     The value of the updating equation
  #
  partA <- 0.5 * sum(colSums(t(Sigma.inv %*% mat1) * (dAdphi %*% mat2)))
  partB <- (-q/2) * sum(colSums(t(A.inv) * dAdphi))
  output <- partA + partB
  return(output)
}


H.tau2 <- function(Sigma.inv, mat1, mat2, q, A.inv){
  #
  #  Purpose:
  #     Calculate the value of the updating equation for the covariance function parameter tau2
  #
  #  Input:
  #     Sigma.inv - The inverse of the current covariance matrix(response)
  #     mat1 - The current value of the expression t((z - M)) %*% A.inv
  #     mat2 - The current value of the expression A.inv %*% (z - M)
  #     q - The number of response variables
  #     A.inv - The inverse of the current covariance matrix
  #
  #  Output:
  #     The value of the updating equation
  #
  partA <- 0.5 * sum(colSums(Sigma.inv * (mat1 %*% mat2)))
  partB <- (-q/2) * sum(diag(A.inv))
  output <- partA + partB
  return(output)
}

K_pi <- function(pi_t){
  #
  #  Purpose:
  #     Define a compact set for the truncation step
  #
  #  Input:
  #     pi_t - A number used to define the bounds of the compact set
  #
  #  Output:
  #     A list containing the bounds of the compact set
  #
  K.B <- 3 + pi_t
  K.Sigma <- 3 + pi_t
  K.phi.min <- -pi_t
  K.phi.max <- 5 + pi_t
  K.tau2 <- 3 + pi_t
  return(list(K.B = K.B, K.Sigma = K.Sigma, K.phi.min = K.phi.min, K.phi.max = K.phi.max, K.tau2 = K.tau2))
}

reinitialise <- function(K.0, p, q) {
  #
  #  Purpose:
  #     Randomly generate a list of parameters from a compact set
  #
  #  Input:
  #     K.0 - The bounds of a compact set
  #     p - The number of covariates + 1
  #     q - The number of responses
  #
  #  Output:
  #     A list containing the generated parameters
  #
  init.B <- matrix(replicate(p*q, runif(1, -K.0$K.B, K.0$K.B)), nrow = p, ncol = q)
  init.U <- matrix(0, q, q)
  init.U[upper.tri(init.U, TRUE)] <- replicate((q*(q+1)/2), runif(1, -K.0$K.Sigma, K.0$K.Sigma))
  init.Sigma <- t(init.U) %*% init.U
  #init.U <- chol(init.Sigma)
  init.phi <- runif(1, K.0$K.phi.min, K.0$K.phi.max)
  init.tau2 <- runif(1, -K.0$K.tau2, K.0$K.tau2)
  return(list(init.B = init.B, init.Sigma = init.Sigma, init.phi = init.phi, init.tau2 = init.tau2, init.U = init.U))
}

in.Kspace <- function(B, Sigma, phi, tau2, K.t) {
  #
  #  Purpose:
  #     Check whether a set of parameters is within a compact set
  #
  #  Input:
  #     B - The parameters for the mean function
  #     U - The upper triangular factor of the Choleski decomposition of Sigma
  #     phi - The covariance function parameter phi
  #     tau2 - The covariance function parameter tau2
  #     K.t - The current compact set
  #
  #  Output:
  #     A logical value to indicate whether the set of parameters is in the compact set or not
  #
  if (any(abs(Sigma) > K.t$K.Sigma)) {
    return(FALSE)
  } else if (phi < K.t$K.phi.min || phi > K.t$K.phi.max) {
    return(FALSE)
  } else if (abs(tau2) > K.t$K.tau2) {
    return(FALSE)
  } else if (any(abs(B) > K.t$K.B)) {
    return(FALSE)
  } else {
    return(TRUE)
  }
}

mu <- function(C, B){
  #
  #  Purpose:
  #     Calculate the mean function
  #
  #  Input:
  #     C - A matrix containing the covariate values
  #     B - A matrix containing the parameters for the mean function
  #
  #  Output:
  #     The mean function values
  #
  mu.out <- C %*% B
  return(mu.out)
}

a_t <- function(a.0, t.0, t, v){
  #
  #  Purpose:
  #     Calculate the gradient parameter for the stochastic approximation
  #
  #  Input:
  #     a.0 - The base gradient
  #     t.0 - The first step where the gradient begins to decrease
  #     t - The current step number
  #     v - A parameter to adjust the rate at which the gradient decreases
  #
  #  Output:
  #     The gradient parameter for the stochastic approximation
  #
  output <- a.0 * (t.0/max(t**v, t.0))
  return(output)
}

b_t <- function(b.0, t.0, t, eta) {
  #
  #  Purpose:
  #     Calculate the difference criterion for the truncation step
  #
  #  Input:
  #     b.0 - The base difference criterion
  #     t.0 - The first step where the difference criterion begins to decrease
  #     t - The current step number
  #     eta - A parameter to adjust the rate at which the difference criterion decreases
  #
  #  Output:
  #     The difference criterion for the truncation step
  #
  output <- b.0 * (t.0/max(t, t.0))**eta
  return(output)
}

rsa.iteration <- function(samp.i, Y, S, covars, m, theta, p, q, a, type) {
  #
  #  Purpose:
  #     Calculate the updated parameter values for a single iteration of the stochastic approximation
  #
  #  Input:
  #     samp.i - The indicies of the subsample in the full sample
  #     Y - The response values from the full sample
  #     S - The coordinate values from the full sample
  #     covars - The covariate values from the full sample
  #     m - The size of the subsample
  #     theta - A list containing the current parameter values
  #     p - The number of covariates + 1
  #     q - The number of response variables
  #     a - The current gradient parameter for the stochastic approximation
  #     type - A string indicating which of the available covariance functions should be used
  #
  #  Output:
  #     A list containing the updated parameter values
  #
  z <- Y[samp.i,]
  s <- S[samp.i,]
  c <- covars[samp.i,]
  if (!is.matrix(z)) z <- matrix(z, ncol = 1)
  if (!is.matrix(c)) c <- matrix(c, ncol = 1)
  B <- theta[[1]]
  Sigma <- theta[[2]]
  phi <- exp(theta[[3]])
  tau2 <- exp(theta[[4]])
  U <- theta[[5]]
  M  <- mu(c, B)
  R <- calcKR(s, s, phi = phi, type = type)
  dAdphi  <- calcKDerivR(s, s, phi = phi, type = type)
  A <- R + tau2 * diag(1, m, m)
  Sigma.inv <- solve(Sigma)
  A.inv <- solve(A)
  mat1 <- t((z - M)) %*% A.inv
  mat2 <- A.inv %*% (z - M)
  H.theta <- vector(mode = "list", length = 4)
  H.theta[[1]] <- H.B(c, z, M, A.inv, Sigma.inv)
  H.theta[[2]] <- H.U(U, Sigma.inv, mat1, z, M, m)
  H.theta[[3]] <- H.phi(A.inv, dAdphi, Sigma.inv, mat1, mat2, q)
  H.theta[[4]] <- H.tau2(Sigma.inv, mat1, mat2, q, A.inv)
  theta.half <- theta
  theta.half[[1]] <- B + a * H.theta[[1]]
  theta.half[[5]] <- U + a * H.theta[[2]]
  theta.half[[2]] <- t(theta.half[[5]]) %*% theta.half[[5]] 
  theta.half[[3]] <- theta[[3]] * exp(a * theta[[3]] * H.theta[[3]])
  theta.half[[4]] <- theta[[4]] * exp(a * theta[[4]] * H.theta[[4]])
  theta.half[[5]] <- chol(theta.half[[2]])
  return(theta.half)
}


RSA.parallel.multi <- function(Y, S, covars = NULL, m , K_pi, t.0, a_t, a.0, v, b_t, b.0, eta, J, parallel = FALSE, cores = NULL, max.iter = 2500, type = c("Exp", "MatThreeTwo")){
  #
  #  Purpose:
  #     Perform population resampling-based stochastic approximation to estimate parameter values for
  #     a gaussian process regression
  #
  #  Input:
  #     Y - The response values from the full sample
  #     S - The coordinate values from the full sample
  #     covars - The covariate values from the full sample
  #     m - The size of the subsamples to be used
  #     K_pi - The function to be used when defining the compact set
  #     t.0 - The step number when the gradient/difference criterion parameters start decreasing
  #     a_t - The function to be used to calculate the gradient parameter for the stochastic approximation
  #     a.0 - The initial gradient parameter for the stochastic approximation
  #     v -  A parameter to adjust the rate at which the gradient decreases
  #     b_t - The function to be used to calculate the difference criterion for the truncation step
  #     b.0 - The initial difference criterion for the truncation step
  #     eta -  A parameter to adjust the rate at which the difference criterion decreases
  #     J - The number of parallel chains which should be evaluated
  #     parallel - A logical value indicating whether to use R's parallel backend.  
  #                Defaults to FALSE however the function is likely to be much faster with the parallel when
  #                m is large(ie. >250) and J>1
  #     cores - The number of parallel processors to be used if using the parallel backend
  #     max.iter - The total number of iterations of the stochastic approximation which should be performed
  #     type - A string indicating which of the available covariance functions should be used
  #
  #  Output:
  #     A list containing:
  #                       trace - A matrix whose ith row contains the parameter estimates at the ith iteration
  #                       B - The parameters for the mean function
  #                       Sigma - The the current covariance matrix(response)
  #                       phi - The covariance function parameter phi
  #                       tau2 - The covariance function parameter tau2 
  #
  type <- match.arg(type)
  if (parallel == TRUE) {
    require(parallel)
    require(foreach)
    require(doParallel)
    cl <- makeCluster(cores)
    if (.Platform$OS.type == "unix") {
      clusterCall(cl, function() dyn.load(paste("C/Unix/calcK", type, "2r.so", sep = "")))
      clusterCall(cl, function() dyn.load(paste("C/Unix/calcK", type, "2Derivr.so", sep = "")))
      clusterCall(cl, function() dyn.load("C/Unix/distMat2R.so"))
    } else {
      clusterCall(cl, function() dyn.load(paste("C/Windows/calcK", type, "2r.dll", sep = "")))
      clusterCall(cl, function() dyn.load(paste("C/Windows/calcK", type, "2Derivr.dll", sep = "")))
      clusterCall(cl, function() dyn.load("C/Windows/distMat2R.dll"))
    }
    registerDoParallel(cl)
  } else {
    if (.Platform$OS.type == "unix") {
      dyn.load(paste("C/Unix/calcK", type, "2r.so", sep = ""))
      dyn.load(paste("C/Unix/calcK", type, "2Derivr.so", sep = ""))
      dyn.load("C/Unix/distMat2R.so")
    } else {
      dyn.load(paste("C/Windows/calcK", type, "2r.dll", sep = ""))
      dyn.load(paste("C/Windows/calcK", type, "2Derivr.dll", sep = ""))
      dyn.load("C/Windows/distMat2R.dll")
    }
  }
  t <- 0
  n <- nrow(Y)
  covars <- cbind(rep(1, n), covars)
  p <- ncol(covars)
  q <- ncol(Y)
  pi_t <- 0
  a_t <- get(a_t, mode = "function")
  b_t <- get(b_t, mode = "function")
  K_pi <- get(K_pi, mode = "function")
  K.0 <- do.call(K_pi, list(pi_t = pi_t))
  K.t <- K.0
  theta <- reinitialise(K.0, p, q)
  thetas.cols <- (2*p*q + q**2 + q + 4)/2
  thetas <- matrix(ncol = thetas.cols, nrow = max.iter + 1)
  thetas[1,] <- c(theta[[1]], theta[[2]][upper.tri(theta[[2]], TRUE)], exp(theta[[3]]), exp(theta[[4]]))
  pb = txtProgressBar(min = 0, max = max.iter, initial = 0, style = 3)
  repeat {
    samples <- matrix(NA, ncol = J, nrow = m)
    for (j in 1:J) {
      samples[,j] <- sample.int(n = n, size = m, replace = FALSE)
    }
    a <- do.call(a_t, list(a.0 = a.0, t.0 = t.0, t = t, v = v))
    if (parallel == TRUE) {
      theta.half <- foreach(j = 1:J, .combine = function(u, v) {return(mapply(`+`, u, v))}, .export = c("rsa.iteration", "calcKR", "calcKDerivR", "H.B", "H.phi", "H.tau2", "H.Sigma", "mu", "euclid.dist")) %dopar% {
        rsa.iteration(samples[,j], Y, S, covars, m, theta, p, q, a, type)
      }
      theta.half <- lapply(theta.half, function(X) {return(X/J)})
    } else {
      theta.half = 0
      for (j in 1:J) {
        theta.half <- mapply(`+`, theta.half, rsa.iteration(samples[,j], Y, S, covars, m, theta, p, q, a, type))
      }
      theta.half <- lapply(theta.half, function(X) {return(X/J)})
    }
    b <- do.call(b_t, list(b.0 = b.0, t.0 = t.0, t = t, eta = eta))
    in.K.t <- in.Kspace(theta.half[[1]], theta.half[[2]], theta.half[[3]], theta.half[[4]], K.t)
    if (in.K.t && (euclid.dist(matrix(unlist(theta[1:4]), nrow = 1), matrix(unlist(theta.half[1:4]), nrow = 1))[1,1] <= b)) {
      theta <- theta.half
    } else {
      theta <- reinitialise(K.0, p, q)
      pi_t <- pi_t + 1
      K.t <- do.call(K_pi, list(pi_t = pi_t))
    }
    t <- t + 1
    setTxtProgressBar(pb,t)
    thetas[t + 1,] <- c(theta[[1]], theta[[2]][upper.tri(theta[[2]], TRUE)], exp(theta[[3]]), exp(theta[[4]]))
    if (t >= max.iter) {
      if (parallel == TRUE) {
        if (.Platform$OS.type == "unix") {
          clusterCall(cl, function() dyn.unload(paste("C/Unix/calcK", type, "2r.so", sep = "")))
          clusterCall(cl, function() dyn.unload(paste("C/Unix/calcK", type, "2Derivr.so", sep = "")))
          clusterCall(cl, function() dyn.unload("C/Unix/distMat2R.so"))
        } else {
          clusterCall(cl, function() dyn.unload(paste("C/Windows/calcK", type, "2r.dll", sep = "")))
          clusterCall(cl, function() dyn.unload(paste("C/Windows/calcK", type, "2Derivr.dll", sep = "")))
          clusterCall(cl, function() dyn.unload("C/Windows/distMat2R.dll"))
        }
        stopCluster(cl)
      } else {
        if (.Platform$OS.type == "unix") {
          dyn.unload(paste("C/Unix/calcK", type, "2r.so", sep = ""))
          dyn.unload(paste("C/Unix/calcK", type, "2Derivr.so", sep = ""))
          dyn.unload("C/Unix/distMat2R.so")
        } else {
          dyn.unload(paste("C/Windows/calcK", type, "2r.dll", sep = ""))
          dyn.unload(paste("C/Windows/calcK", type, "2Derivr.dll", sep = ""))
          dyn.unload("C/Windows/distMat2R.dll")
        }
      }
      return(list(trace = thetas, B = theta[[1]], Sigma = theta[[2]], phi = exp(theta[[3]]), tau2 = exp(theta[[4]])))
    }
  }
}

predict.rsa.multi.local <- function(s.train, x.train, y.train, S = s.train, X = x.train, B, Sigma, phi, tau2, delta, type){
  #
  #  Purpose:
  #     Calculate predicted mean covariance functions of a Gaussian process
  #
  #  Input:
  #     s.train - The coordinate values to be used for kriging
  #     x.train - The covariate values to be used for kriging
  #     s.train - The response values to be used for kriging
  #     S - The new coordinate values where the functions are to be evaluated
  #     X - The covariate values at the new coordinates
  #     B - The parameters for the mean function
  #     Sigma - The the current covariance matrix(response)
  #     phi - The covariance function parameter phi
  #     tau2 - The covariance function parameter tau2
  #     type - A string indicating which of the available covariance functions should be used
  #
  #  Output:
  #     A list containing the predicted mean function values  predicted covariance function values
  #
  if (.Platform$OS.type == "unix") {
    dyn.load(paste("C/Unix/calcK", type, "2r.so", sep = ""))
    dyn.load(paste("C/Unix/calcK", type, "2Derivr.so", sep = ""))
    dyn.load("C/Unix/distMat2R.so")
  } else {
    dyn.load(paste("C/Windows/calcK", type, "2r.dll", sep = ""))
    dyn.load(paste("C/Windows/calcK", type, "2Derivr.dll", sep = ""))
    dyn.load("C/Windows/distMat2R.dll")
  }
  D <- euclid.dist(X1 = s.train, X2 = S)
  n2 <- nrow(S)
  pred.mean <- matrix(NA, nrow = n2, ncol = ncol(y.train))
  pred.var <- matrix(NA, nrow = n2, ncol = ncol(y.train))
  for (i in 1:n2) {
    loc.train <- which(D[,i] < delta)
    loc.s.train <- s.train[loc.train,]
    loc.x.train <- x.train[loc.train,]
    loc.y.train <- y.train[loc.train,]
    m <- length(loc.train)
    if (m == 1) {
      loc.s.train <- matrix(loc.s.train, nrow = 1)
      loc.x.train <- matrix(loc.x.train, nrow = 1)
      loc.y.train <- matrix(loc.y.train, nrow = 1)
    }
    DesMat.train <- cbind(rep(1, m), loc.x.train)
    DesMat.X <- cbind(1, matrix(X[i,], nrow = 1))
    mean.func.train <- mu(DesMat.train, B)
    mean.func.X <- mu(DesMat.X, B)
    A.st_t <- calcKR(matrix(S[i,], nrow = 1), loc.s.train, phi, type)
    A.t_st <- calcKR(loc.s.train, matrix(S[i,], nrow = 1), phi, type)
    A.st_st <- calcKR(matrix(S[i,], nrow = 1), matrix(S[i,], nrow = 1), phi, type)
    A.t_t <- calcKR(loc.s.train, loc.s.train, phi, type)
    A <- A.t_t + tau2 * diag(1, m, m)
    resids <- loc.y.train - mean.func.train
    K.st_t <- kronecker(Sigma, A.st_t)
    K.t_st <- kronecker(Sigma, A.t_st)
    K.st_st <- kronecker(Sigma, A.st_st)
    K.chol <- t(kronecker(t(chol(Sigma)), t(chol(A))))
    K.inv <- chol2inv(K.chol)
    mean.func.X <- as.vector(mean.func.X)
    resids <- as.vector(resids)
    pred.mean[i,] <- mean.func.X + K.st_t %*% K.inv %*% resids
    pred.var[i,] <- diag((K.st_st - K.st_t %*% K.inv %*% K.t_st))
  }
  if (.Platform$OS.type == "unix") {
    dyn.unload(paste("C/Unix/calcK", type, "2r.so", sep = ""))
    dyn.unload(paste("C/Unix/calcK", type, "2Derivr.so", sep = ""))
    dyn.unload("C/Unix/distMat2R.so")
  } else {
    dyn.unload(paste("C/Windows/calcK", type, "2r.dll", sep = ""))
    dyn.unload(paste("C/Windows/calcK", type, "2Derivr.dll", sep = ""))
    dyn.unload("C/Windows/distMat2R.dll")
  }
  return(list(pred.mean = pred.mean, pred.var = pred.var))
}

