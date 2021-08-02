
imdh <- function(X, sol = NULL, depth = 8, hmult = 1, epoch = 2, C = 10, alpha = .1, tmax = 100000000, scale = 1, t_init0 = NULL, refine = TRUE, k = NULL){
  n <- nrow(X)
  d <- ncol(X)
  if(is.null(sol)){
    nterms <- 2^depth-1
    b <- numeric(nterms)
    if(is.null(t_init0)) t_init0 <- d
    mu0 <- colMeans(X[1:t_init0,])
    scl <- colMeans(X[1:t_init0,]^2)
    V <- diag(d)[,rep(1:d, ceiling(nterms/d))[1:nterms]] + 1e-50
    #V <- matrix(rnorm(d*nterms), d, nterms)/sqrt(d)
    ss <- numeric(nterms)
    ssss <- ss^4
    mn <- V*0
    sol <- list(clusters = numeric(n), V = V, b = b, ss = ss, ssss = ssss, mean = mn, t_init = rep(0, nterms), hmult = hmult, C = C, alpha = alpha, tmax = tmax, scl = scl, mu0 = mu0, scale = scale, t_init0 = t_init0)
  }
  for(j in 1:epoch){
    print('going')
    updated <- smdh_list_div_cpp(X, sol$V, sol$b, sol$mean, sol$ss, sol$ssss*0, sol$t_init, sol$hmult, sol$C, sol$tmax, sol$alpha, sol$clusters, sol$scl, sol$mu0, sol$scale, sol$t_init0)
    sol$V <- updated[[1]]
    sol$b <- updated[[2]]
    sol$mean <- updated[[3]]
    sol$ss <- updated[[4]]
    sol$ssss <- updated[[5]]
    sol$t_init <- updated[[6]]
    sol$clusters <- updated[[7]] + 1
    sol$scl <- updated[[8]]
    sol$mu0 <- updated[[9]]
    sol$t_init0 <- updated[[10]]
  }
  sol$k4 <- sol$ssss
  if(refine) sol$clusters <- refine(sol, k = k)
  sol
}


imdh_assign <- function(X, sol, refine = TRUE, k = NULL){
  n <- nrow(X)
  d <- ncol(X)
  updated <- smdh_list_div_pass_cpp(X, sol$V, sol$b, sol$mean, sol$ssss*0, sol$scl, sol$mu0, sol$scale)
  sol$clusters <- updated[[1]] + 1
  sol$ssss <- updated[[2]]
  sol$k4 <- sol$ssss
  if(refine) sol$clusters <- refine(sol, k = k)
  sol
}



refine <- function(sol, k = NULL, kmax = Inf){
  dp <- log(length(sol$k4)+1)/log(2)
  leaves <- (2^(dp-1)):length(sol$k4)
  C <- matrix(sol$clusters, length(leaves), length(sol$clusters), byrow = TRUE)
  K <- numeric(length(leaves))
  nk4 <- sol$k4
  k0 <- sum(nk4[leaves])
  preleaves <- (2^(dp-2)):(2^(dp-1)-1)
  K[1] <- k0
  for(j in 1:(length(leaves)-1)){
    kopt <- -Inf
    zopt <- 0
    kids <- c(0,0)
    for(z in preleaves){
      ids <- 2*z+(0:1)
      kdiff <- sum(nk4[ids]) - nk4[z]
      if(kdiff > kopt){
        kopt <- kdiff
        zopt <- z
        kids <- ids
      }
    }
    K[j+1] <- K[j] - sum(nk4[kids]) + nk4[zopt]
    C[j+1,] <- C[j,]
    C[j+1,which(C[j+1,]%in%kids)] <- min(kids)/2
    if(zopt%%2==0){
      if((zopt+1)%in%leaves) preleaves <- c(preleaves[preleaves!=zopt], zopt/2)
      else preleaves <- preleaves[preleaves!=zopt]
    }
    else{
      if((zopt-1)%in%leaves) preleaves <- c(preleaves[preleaves!=zopt], (zopt-1)/2)
      else preleaves <- preleaves[preleaves!=zopt]
    }
    leaves <- c(leaves, zopt)
  }
  if(is.null(k)){
    kmax <- min(kmax, length(K))
    Ko <- K[length(K):(length(K)-kmax+1)]
    if(kmax > 10){
      kopts <- numeric(30)
      kseq <- ceiling(seq(10, kmax, length = 30))
      for(i in 1:30) kopts[i] <- get_elbow(Ko[1:kseq[i]])
      counts <- sapply(1:kmax, function(k) sum(kopts==k))
      kopt <- length(K) - which.max(counts) + 1
    }
    else kopt <- length(K) - get_elbow(Ko) + 1
  }
  else kopt <- length(K) - k + 1
  as.numeric(as.factor(C[kopt,]))
}

get_elbow <- function(yvals){
  d <- length(yvals)
  angles <- sapply(1:(d-1), function(i){
    a1 <- atan((i-1)/(d-1)*(max(yvals)-min(yvals))/(max(yvals)-yvals[i]))
    a2 <- atan((yvals[i]-min(yvals))/(max(yvals)-min(yvals))*(d-1)/(d-i))
    return(a1+a2)
  })
  which.min(angles)
}
