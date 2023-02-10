#' @name simdata
#' @title Simulation Data
#' @description
#' Simulation Data with a sparse structural \code{B} and a sparse \code{A} with orthonormal columns.
#'
#' @param typeA One of \code{sparse0.1}, \code{sparse0.2}, \code{sparse0.3}, \code{sparse0.4}, \code{stiefel} and \code{perfect-simple}.
#'
#' @param typeB One of \code{null}, \code{row}, \code{block}, \code{individual}, \code{joint}, \code{both}, \code{partial} and \code{all}.
#'
#' @param n The number of samples. Default is 50.
#' @param d The number of multi-source datasets. Default is 3.
#' @param q The number of responses. Default is 10.
#' @param pvec A vector representing the number of variables for each dataset. Default is \code{rep(10, 3)}.
#' @param rvec A vector representing the rank of each dataset. Default is \code{rep(2, 3)}.
#' @param es Effect size: choose one of "1", "2", ..., "5" and "norm". Default is "5".
#' @param simplify A logical value for simplifying the matrix of \code{A} using the quartimax rotation. Default is TRUE.
#' @param snr A numeric value for SNR (signal to noise). Default is 1.
#' @param sigma A numeric value of a multiplying constant on the variance of an error term. Default is NULL. If NULL, sigma = |Y0|_F^2 / |E|_F^2) / sqrt(snr), where Y0 = XBA^T.
#' @param rho_X A numeric value for controlling the correlation between predictors.
#' @param rho_E A numeric value for controlling the correlation between errors.
#'
#' @importFrom mnormt rmnorm
#' @author Kipoong Kim <statpng@snu.ac.kr>
#' @references
#' No reference.
#' @return A list with Y, X, C, A, B, sigma, Xsigma, n, d, q, p, rvec, sim.params
#'
#' @examples
#'
#'
#' data <- rrr.sim(typeA="sparse0.4", typeB="all", n=20, d=3, q=10, p=5, rvec=rep(2,3), es="1", simplify=TRUE, snr=1)
#'
#' data$B |> round(3)
#' data$A |> round(3)
#'
#'

#' @export simdata2
simdata2 <- function(typeA="sparse",
                     typeB="all",
                     n = 50,
                     p = 10,
                     q = 10,
                     d = 3,
                     rvec = rep(2,3),
                     nuA=0.2,
                     nuB=0.5,
                     d0 = 3,
                     es = "1",
                     es.B = 1,
                     snr = 1,
                     simplify = TRUE,
                     sigma = NULL,
                     rho_X = 0.5,
                     rho_E = 0){
  
  
  pvec = rep(p, d)
  
  
  es <- as.character(es)

  if(FALSE){
    {
      typeA="sparse0"
      typeB="all"
      n <- 50;
      d <- 3;
      q <- 10;
      q0 <- 10
      pvec <- rep(10,3)
      rvec = rep(1, 10);
      es <- "5" # "norm"
      snr <- 1
      sigma = NULL
      rho_X = 0.5
      rho_E = 0
      nuA = 0.2
      nuB = 0.5
    }



    {
      set.seed(1)
      data <- simdata(typeA="sparse0.2", typeB="random", n=20, d=5, q=20, pvec=rep(50,5), rvec=rep(1,5), es="2", simplify=TRUE, snr=1)

      data$sigma
    }
  }

  CorrAR <- function(p, rho){
    rho^outer(1:p, 1:p, function(x, y) abs(x-y))
  }

  CorrCS <- function(p, rho)
  {
    Sigma <- matrix(nrow = p, ncol = p, rho)
    diag(Sigma) <- 1
    Sigma
  }



  Dist <- switch(es,
                 `1` = function(n) runif(n, 0.1, 0.3),
                 `2` = function(n) runif(n, 0.5, 1.0),
                 `3` = function(n) runif(n, 1.0, 2.0),
                 `4` = function(n) runif(n, 5.0, 10.0),
                 `5` = function(n) runif(n, 10.0, 20.0),
                 `norm` = function(n) rnorm(n))


  pvec = rep(1:length(pvec), pvec)
  p = length(pvec)
  r = sum(rvec)
  r0 = r
  q0 = q


  rvec2 = rep( seq_len(length(rvec)), rvec )


  Sigma = CorrCS
  A1 <- matrix(ncol = r, nrow = q0, rnorm(q0 * r))
  A0 <- matrix(ncol = r, nrow = q - q0, 0)
  A <- rbind(A1, A0)
  A <- svd(A)$u %*% t(svd(A)$v)


  if( typeA == "sparse0.0" ){

    A <- ToOrthogonal(q, r, 0.0, simplify=simplify)

  } else if( typeA == "sparse0.1" ){

    A <- ToOrthogonal(q, r, 0.1, simplify=simplify)

  } else if( typeA == "sparse0.2" ){

    A <- ToOrthogonal(q, r, 0.2, simplify=simplify)

  } else if( typeA == "sparse0.3" ){

    A <- ToOrthogonal(q, r, 0.3, simplify=simplify)

  } else if( typeA == "sparse0.4" ){

    A <- ToOrthogonal(q, r, 0.4, simplify=simplify)

  } else if( typeA == "perfect-simple" ){

    A <- NULL
    gg <- ceiling(q/r)
    for( h in 1:gg ){
      if( h < gg ){
        A <- rbind(A, diag(1, r, r))
      } else {
        A <- rbind(A, diag(1, q%%r, r))
      }
    }

  } else if( typeA == "stiefel" ){

    rustiefel <- function(q, r){
      X <- matrix(rnorm(q*r), q, r)
      tmp <- eigen(t(X) %*% X)
      X %*% (tmp$vec %*% sqrt(diag(1/tmp$val, nrow = R)) %*% t(tmp$vec))
    }

    A <- rstiefel::rustiefel(q, r)
  } else if( typeA == "quartimax" ){
    
    A <- ToOrthogonal(q, r, 0, simplify=simplify)
    A[apply(A, 1, norm, "2")<nuA,] <- 0
    A <- GPArotation::quartimax( with( svd(A), tcrossprod(u,v) ) )$loadings
    
  } else if( typeA == "BlockDiagonal" ){
    
    q_r <- floor(q/r)
    
    V <- matrix(0,ncol=r ,nrow=q)
    
    vec_append <- c(rep(1,q_r),rep(0,q-q_r))
    for( kk in 1:r ){
      
      V[,kk] <- vec_append*sample(c(1,-1),q,replace=TRUE)*runif(q,0.3,1)
      
      vec_append <- c(rep(0,q_r), vec_append)[1:q]
    }
    
    A <- apply(V,2,function(x) x/sqrt(sum(x^2)))
    A[apply(A, 1, norm, "2")<nuA,] <- 0
    A <- with(svd(A), tcrossprod(u,v))
    
  }




  B <- matrix(0, nrow = p, ncol = r)

  # typeB <- cbind(diag(1,d,d), matrix(0,d,r-d) )

  if( typeB == "null" ){
    typeB <- matrix(0, d, r)
  } else if( typeB == "row" ){
    typeB <- rbind(matrix(1, d0, r), matrix(0, d-d0, r))
  } else if( typeB == "block" ){
    typeB <- rbind(
      cbind(matrix(1, d0, r0), matrix(0, d0, r-r0)),
      cbind(matrix(0, d-d0, r0), matrix(0, d-d0, r-r0))
    )
  } else if( typeB == "individual" ){
    typeB <- matrix(0, d, r)
    for( i.d0 in 1:d0 ){
      typeB[i.d0, 1:2+2*(i.d0-1)] <- 1
    }

    # typeB <- rbind(
    #   cbind(diag(1, d0, d0), matrix(0, d0, r-d0)),
    #   cbind(matrix(0, d-d0, d0), matrix(0, d-d0, r-d0))
    # )
  } else if( typeB == "joint" ){
    typeB <- cbind(matrix(1, d, r0), matrix(0, d, r-r0))
  } else if( typeB == "both" ){
    typeB <- cbind( matrix(1,d,1),
                    rbind(
                      cbind(diag(1, d0, d0), matrix(0, d0, r-d0-1)),
                      cbind(matrix(0, d-d0, d0), matrix(0, d-d0, r-d0-1))
                    ) )
  } else if( typeB == "partial" ){
    typeB <- matrix(0,d,r)
    typeB[c(1,2),1:2] <- 1
    typeB[c(1,3),3:4] <- 1
    typeB[c(2,3),5:6] <- 1
  } else if( typeB == "all" ){
    typeB <- matrix(0,d,r)
    typeB[1:3,1] <- 1
    typeB[1,2] <- 1
    typeB[2,3] <- 1
    typeB[c(1,2),4] <- 1
    typeB[c(1,3),5] <- 1
    typeB[c(2,3),6] <- 1
  } else if( typeB == "random" ){
    typeB <- matrix(0,d,r)

    GRID.B <- expand.grid( 1:d, 1:r )
    GRID.B.sel <- GRID.B[ sample(nrow(GRID.B), floor(nrow(GRID.B)*nuB)), ]
    for( jj in 1:nrow(GRID.B.sel) ){
      typeB[GRID.B.sel[jj,1], GRID.B.sel[jj,2]] <- 1
    }

  }

  B.str2 <- kronecker(typeB, rep(1,p/d))
  B <- ifelse( B.str2 == 1, Dist(sum(B.str2==1)), 0 ) * es.B

  attr(B, "pvec") <- as.numeric(table(pvec))




  C <- B %*% t(A)
  Xsigma <- lapply( table(pvec), function(pi) Sigma(pi, rho_X) )
  Xsigma <- Matrix::bdiag(Xsigma)
  # Xsigma <- Sigma(p, rho_X)
  X <- mnormt::rmnorm(n, rep(0, p), Xsigma)
  UU <- mnormt::rmnorm(n, rep(0, q), Sigma(q, rho_E))
  UU.0 <- UU
  svdC <- svd(C)
  C3 <- svdC$u[, r] %*% t(svdC$v[, r]) * svdC$d[r]
  Y3 <- X %*% C3
  if (is.null(sigma)) {
    sigma <- sqrt(sum(as.numeric(Y3)^2)/sum(as.numeric(UU)^2)/snr)
  }
  # snr
  SNR <- sum(as.numeric(Y3)^2)/sum(as.numeric(UU)^2)/sigma^2
  cat( "SNR=", SNR, "\n" )

  UU <- UU * sigma
  Y <- matrix(nrow = n, ncol = q, NA)
  Y <- X %*% C + UU


  sim.params <- list(
    typeA="sparse",
    typeB="all",
    n = n,
    d = d,
    q = q,
    pvec = pvec,
    rvec = rvec,
    es = es,
    snr = snr,
    sigma = sigma,
    rho_X = rho_X,
    rho_E = rho_E
  )



  list(Y = Y, X = X, C = C, A = A, B = B, UU=UU, UU.0 = UU.0, sigma = sigma, snr.out = SNR,
       Xsigma = Xsigma, n=n, d=d, q=q, nrank = r, pvec=pvec, rvec=rvec,
       params = sim.params)



}






#' @export simdata3
simdata3 <- function(typeA="sparse",
                     typeB="all",
                     n = 50,
                     p = 10,
                     q = 10,
                     d = 3,
                     rvec = NULL,
                     nuA=0.2,
                     nuB=0.5,
                     d0 = 3,
                     es = "1",
                     es.B = 1,
                     snr = 1,
                     simplify = TRUE,
                     sigma = NULL,
                     rho_X = 0.5,
                     rho_E = 0){
  
  
  pvec = rep(p, d)
  
  
  es <- as.character(es)
  
  if(FALSE){
    {
      typeA="sparse0"
      typeB="all"
      n <- 50;
      d <- 3;
      q <- 10;
      q0 <- 10
      pvec <- rep(10,3)
      rvec = rep(1, 10);
      es <- "5" # "norm"
      snr <- 1
      sigma = NULL
      rho_X = 0.5
      rho_E = 0
      nuA = 0.2
      nuB = 0.5
    }
    
    
    
    {
      set.seed(1)
      data <- simdata(typeA="sparse0.2", typeB="random", n=20, d=5, q=20, pvec=rep(50,5), rvec=rep(1,5), es="2", simplify=TRUE, snr=1)
      
      data$sigma
    }
  }
  
  CorrAR <- function(p, rho){
    rho^outer(1:p, 1:p, function(x, y) abs(x-y))
  }
  
  CorrCS <- function(p, rho)
  {
    Sigma <- matrix(nrow = p, ncol = p, rho)
    diag(Sigma) <- 1
    Sigma
  }
  
  
  
  Dist <- switch(es,
                 `1` = function(n) runif(n, 0.1, 0.3),
                 `2` = function(n) runif(n, 0.5, 1.0),
                 `3` = function(n) runif(n, 1.0, 2.0),
                 `4` = function(n) runif(n, 5.0, 10.0),
                 `5` = function(n) runif(n, 10.0, 20.0),
                 `norm` = function(n) rnorm(n))
  
  
  pvec = rep(1:length(pvec), pvec)
  p = length(pvec)
  r = sum(rvec)
  r0 = r
  q0 = q
  
  
  rvec2 = rep( seq_len(length(rvec)), rvec )
  
  
  Sigma = CorrCS
  A1 <- matrix(ncol = r, nrow = q0, rnorm(q0 * r))
  A0 <- matrix(ncol = r, nrow = q - q0, 0)
  A <- rbind(A1, A0)
  A <- svd(A)$u %*% t(svd(A)$v)
  
  
  if( typeA == "sparse0.0" ){
    
    A <- ToOrthogonal(q, r, 0.0, simplify=simplify)
    
  } else if( typeA == "sparse0.1" ){
    
    A <- ToOrthogonal(q, r, 0.1, simplify=simplify)
    
  } else if( typeA == "sparse0.2" ){
    
    A <- ToOrthogonal(q, r, 0.2, simplify=simplify)
    
  } else if( typeA == "sparse0.3" ){
    
    A <- ToOrthogonal(q, r, 0.3, simplify=simplify)
    
  } else if( typeA == "sparse0.4" ){
    
    A <- ToOrthogonal(q, r, 0.4, simplify=simplify)
    
  } else if( typeA == "perfect-simple" ){
    
    A <- NULL
    gg <- ceiling(q/r)
    for( h in 1:gg ){
      if( h < gg ){
        A <- rbind(A, diag(1, r, r))
      } else {
        A <- rbind(A, diag(1, q%%r, r))
      }
    }
    
  } else if( typeA == "stiefel" ){
    
    rustiefel <- function(q, r){
      X <- matrix(rnorm(q*r), q, r)
      tmp <- eigen(t(X) %*% X)
      X %*% (tmp$vec %*% sqrt(diag(1/tmp$val, nrow = R)) %*% t(tmp$vec))
    }
    
    A <- rstiefel::rustiefel(q, r)
  } else if( typeA == "quartimax" ){
    
    A <- ToOrthogonal(q, r, 0, simplify=simplify)
    A[apply(A, 1, norm, "2")<nuA,] <- 0
    A <- GPArotation::quartimax( with( svd(A), tcrossprod(u,v) ) )$loadings
    
  } else if( typeA == "BlockDiagonal" ){
    
    q_r <- floor(q/r)
    
    V <- matrix(0,ncol=r ,nrow=q)
    
    vec_append <- c(rep(1,q_r),rep(0,q-q_r))
    for( kk in 1:r ){
      
      V[,kk] <- vec_append*sample(c(1,-1),q,replace=TRUE)*runif(q,0.3,1)
      
      vec_append <- c(rep(0,q_r), vec_append)[1:q]
    }
    
    A <- apply(V,2,function(x) x/sqrt(sum(x^2)))
    A[apply(A, 1, norm, "2")<nuA,] <- 0
    A <- with(svd(A), tcrossprod(u,v))
    
  }
  
  
  
  
  B <- matrix(0, nrow = p, ncol = r)
  
  # typeB <- cbind(diag(1,d,d), matrix(0,d,r-d) )
  
  if( typeB == "null" ){
    typeB <- matrix(0, d, r)
  } else if( typeB == "row" ){
    typeB <- rbind(matrix(1, d0, r), matrix(0, d-d0, r))
  } else if( typeB == "block" ){
    typeB <- rbind(
      cbind(matrix(1, d0, r0), matrix(0, d0, r-r0)),
      cbind(matrix(0, d-d0, r0), matrix(0, d-d0, r-r0))
    )
  } else if( typeB == "individual" ){
    typeB <- matrix(0, d, r)
    for( i.d0 in 1:d0 ){
      typeB[i.d0, 1:1+1*(i.d0-1)] <- 1
    }
    
    # typeB <- rbind(
    #   cbind(diag(1, d0, d0), matrix(0, d0, r-d0)),
    #   cbind(matrix(0, d-d0, d0), matrix(0, d-d0, r-d0))
    # )
  } else if( typeB == "joint" ){
    typeB <- cbind(matrix(1, d, r0), matrix(0, d, r-r0))
  } else if( typeB == "both" ){
    typeB <- cbind( matrix(1,d,1),
                    rbind(
                      cbind(diag(1, d0, d0), matrix(0, d0, r-d0-1)),
                      cbind(matrix(0, d-d0, d0), matrix(0, d-d0, r-d0-1))
                    ) )
  } else if( typeB == "partial" ){
    typeB <- matrix(0,d,r)
    typeB[c(1,2),1] <- 1
    typeB[c(1,3),2] <- 1
    typeB[c(2,3),3] <- 1
  } else if( typeB == "all" ){
    typeB <- matrix(0,d,r)
    typeB[1:3,1] <- 1
    typeB[1,2] <- 1
    typeB[2,3] <- 1
    typeB[3,4] <- 1
    typeB[c(1,2),5] <- 1
    typeB[c(1,3),6] <- 1
    typeB[c(2,3),7] <- 1
  } else if( typeB == "random" ){
    typeB <- matrix(0,d,r)
    
    GRID.B <- expand.grid( 1:d, 1:r )
    GRID.B.sel <- GRID.B[ sample(nrow(GRID.B), floor(nrow(GRID.B)*nuB)), ]
    for( jj in 1:nrow(GRID.B.sel) ){
      typeB[GRID.B.sel[jj,1], GRID.B.sel[jj,2]] <- 1
    }
    
  }
  
  B.str2 <- kronecker(typeB, rep(1,p/d))
  B <- ifelse( B.str2 == 1, Dist(sum(B.str2==1)), 0 ) * es.B
  
  attr(B, "pvec") <- as.numeric(table(pvec))
  
  
  
  
  C <- B %*% t(A)
  Xsigma <- lapply( table(pvec), function(pi) Sigma(pi, rho_X) )
  Xsigma <- Matrix::bdiag(Xsigma)
  # Xsigma <- Sigma(p, rho_X)
  X <- mnormt::rmnorm(n, rep(0, p), Xsigma)
  UU <- mnormt::rmnorm(n, rep(0, q), Sigma(q, rho_E))
  UU.0 <- UU
  svdC <- svd(C)
  C3 <- svdC$u[, r] %*% t(svdC$v[, r]) * svdC$d[r]
  Y3 <- X %*% C3
  if (is.null(sigma)) {
    sigma <- sqrt(sum(as.numeric(Y3)^2)/sum(as.numeric(UU)^2)/snr)
  }
  # snr
  SNR <- sum(as.numeric(Y3)^2)/sum(as.numeric(UU)^2)/sigma^2
  cat( "SNR=", SNR, "\n" )
  
  UU <- UU * sigma
  Y <- matrix(nrow = n, ncol = q, NA)
  Y <- X %*% C + UU
  
  
  sim.params <- list(
    typeA="sparse",
    typeB="all",
    n = n,
    d = d,
    q = q,
    pvec = pvec,
    rvec = rvec,
    es = es,
    snr = snr,
    sigma = sigma,
    rho_X = rho_X,
    rho_E = rho_E
  )
  
  
  
  list(Y = Y, X = X, C = C, A = A, B = B, UU=UU, UU.0 = UU.0, sigma = sigma, snr.out = SNR,
       Xsigma = Xsigma, n=n, d=d, q=q, nrank = r, pvec=pvec, rvec=rvec,
       params = sim.params)
  
  
  
}







#' @export simdata
simdata <- function (typeA="sparse",
                     typeB="all",
                     n = 50,
                     d = 3,
                     q = 10,
                     pvec = rep(10,3),
                     rvec = rep(2,3),
                     # Dist=function(n) rnorm(n),
                     es = "5",
                     simplify = TRUE,
                     snr = 1,
                     sigma = NULL,
                     rho_X = 0.5,
                     rho_E = 0, nuB=0.5){

  if(FALSE){
    {
      typeA="sparse0.2"
      typeB="all"
      n <- 50;
      d <- 3;
      q <- 10;
      q0 <- 10
      pvec <- rep(10,3)
      rvec = rep(1, 10);
      es <- "5" # "norm"
      snr <- 1
      sigma = NULL
      rho_X = 0.5
      rho_E = 0
      nuB = 0.5
    }



    {
      set.seed(1)
      data <- simdata(typeA="sparse0.2", typeB="random", n=20, d=5, q=20, pvec=rep(50,5), rvec=rep(1,5), es="2", simplify=TRUE, snr=1)

      data$sigma
    }
  }

  CorrAR <- function(p, rho){
    rho^outer(1:p, 1:p, function(x, y) abs(x-y))
  }

  CorrCS <- function(p, rho)
  {
    Sigma <- matrix(nrow = p, ncol = p, rho)
    diag(Sigma) <- 1
    Sigma
  }



  Dist <- switch(es,
                 `1` = function(n) runif(n, 0.1, 0.3),
                 `2` = function(n) runif(n, 0.5, 1.0),
                 `3` = function(n) runif(n, 1.0, 2.0),
                 `4` = function(n) runif(n, 5.0, 10.0),
                 `5` = function(n) runif(n, 10.0, 20.0),
                 `norm` = function(n) rnorm(n))


  pvec = rep(1:length(pvec), pvec)
  p = length(pvec)
  r = sum(rvec)
  q0 = q


  rvec2 = rep( seq_len(length(rvec)), rvec )


  Sigma = CorrCS
  A1 <- matrix(ncol = r, nrow = q0, rnorm(q0 * r))
  A0 <- matrix(ncol = r, nrow = q - q0, 0)
  A <- rbind(A1, A0)
  A <- svd(A)$u %*% t(svd(A)$v)


  if( typeA == "sparse0.0" ){

    A <- ToOrthogonal(q, r, 0.0, simplify=simplify)

  } else if( typeA == "sparse0.1" ){

    A <- ToOrthogonal(q, r, 0.1, simplify=simplify)

  } else if( typeA == "sparse0.2" ){

    A <- ToOrthogonal(q, r, 0.2, simplify=simplify)

  } else if( typeA == "sparse0.3" ){

    A <- ToOrthogonal(q, r, 0.3, simplify=simplify)

  } else if( typeA == "sparse0.4" ){

    A <- ToOrthogonal(q, r, 0.4, simplify=simplify)

  } else if( typeA == "perfect-simple" ){

    A <- NULL
    gg <- ceiling(q/r)
    for( h in 1:gg ){
      if( h < gg ){
        A <- rbind(A, diag(1, r, r))
      } else {
        A <- rbind(A, diag(1, q%%r, r))
      }
    }

  } else if( typeA == "stiefel" ){

    rustiefel <- function(q, r){
      X <- matrix(rnorm(q*r), q, r)
      tmp <- eigen(t(X) %*% X)
      X %*% (tmp$vec %*% sqrt(diag(1/tmp$val, nrow = R)) %*% t(tmp$vec))
    }

    A <- rstiefel::rustiefel(q, r)

  }



  B <- matrix(ncol = r, nrow = p, 0)

  if( typeB == "null" ){

    B[1:p,1:r] <- 0

  }

  if( typeB == "row" ){

    B[which(pvec==1), 1:r] <- Dist(sum(pvec==1)*r)

  }

  if( typeB == "block" ){
    B[which(pvec==1), which(rvec2==1)] <- Dist(sum(pvec==1)*sum(rvec2==1))
  }

  if( typeB == "individual" ){
    for( ii in 1:d ){
      B[which(pvec==ii), which(rvec2==ii)] <- Dist(sum(pvec==ii)*sum(rvec2==ii))
    }
  }

  if( typeB == "joint" ){
    B[1:p, which(rvec2==1)] <- Dist(p*sum(rvec2==1))
    B[1:p, which(rvec2==1)] <- Dist(p*sum(rvec2==1))
  }

  if( typeB == "joint2" ){
    B[which(pvec==1), which(rvec2==1)[1]] <- Dist(sum(pvec==1)*1)*1
    B[which(pvec==2), which(rvec2==1)[1]] <- Dist(sum(pvec==2)*1)*2
    B[which(pvec==3), which(rvec2==1)[1]] <- Dist(sum(pvec==3)*1)*3

    B[which(pvec==1), which(rvec2==1)[2]] <- Dist(sum(pvec==1)*1)*4
    B[which(pvec==2), which(rvec2==1)[2]] <- Dist(sum(pvec==2)*1)*5
    B[which(pvec==3), which(rvec2==1)[2]] <- Dist(sum(pvec==3)*1)*6
  }

  if( typeB == "both" ){

    for( i in 1:d ){
      B[1:p, which(rvec2 == 1)] <- Dist(p*sum(rvec2==1))
      B[which(pvec==i),which(rvec2 == (i+1))] <- Dist(sum(pvec==i)*sum(rvec2 == (i+1)))
    }

  }

  if( typeB == "partial" ){
    for( ii in 1:d ){

    }

    B[which(pvec %in% c(1,2)), which(rvec2==1)] <- Dist( sum(pvec %in% c(1,2))*sum(rvec2==1) )
    B[which(pvec %in% c(1,3)), which(rvec2==2)] <- Dist( sum((pvec %in% c(1,3)))*sum(rvec2==2) )
    B[which(pvec %in% c(2,3)), which(rvec2==3)] <- Dist( sum((pvec %in% c(2,3)))*sum(rvec2==3) )
  }

  if( typeB == "all" ){
    B[1:p, 1] <- Dist(p*1)
    B[which(pvec%in%c(1)), 2] <- Dist( sum(pvec%in%c(1)) )
    B[which(pvec%in%c(2)), 3] <- Dist( sum(pvec%in%c(2)) )

    B[which(pvec%in%c(1,2)), 4] <- Dist( sum(pvec%in%c(1,2)) )
    B[which(pvec%in%c(1,3)), 5] <- Dist( sum(pvec%in%c(1,3)) )
    B[which(pvec%in%c(2,3)), 6] <- Dist( sum(pvec%in%c(2,3)) )
  }

  if( typeB == "random" ){

    GRID.B <- expand.grid(  1:d, 1:r  )
    GRID.B.sel <- GRID.B[ sample(nrow(GRID.B), floor(nrow(GRID.B)*nuB)), ]
    for( jj in 1:nrow(GRID.B.sel) ){
      B[which(pvec==GRID.B.sel[jj,1]), GRID.B.sel[jj,2]] <- Dist( sum(pvec==1) )
    }

  }



  C <- B %*% t(A)
  Xsigma <- lapply( table(pvec), function(pi) Sigma(pi, rho_X) )
  Xsigma <- Matrix::bdiag(Xsigma)
  # Xsigma <- Sigma(p, rho_X)
  X <- rmnorm(n, rep(0, p), Xsigma)
  UU <- rmnorm(n, rep(0, q), Sigma(q, rho_E))
  UU.out <- UU
  svdC <- svd(C)
  C3 <- svdC$u[, r] %*% t(svdC$v[, r]) * svdC$d[r]
  Y3 <- X %*% C3
  if (is.null(sigma)) {
    sigma <- sqrt(sum(as.numeric(Y3)^2)/sum(as.numeric(UU)^2)/snr)
  }
  # snr
  SNR <- sum(as.numeric(Y3)^2)/sum(as.numeric(UU)^2)/sigma^2
  cat( "SNR=", SNR, "\n" )

  UU <- UU * sigma
  Y <- matrix(nrow = n, ncol = q, NA)
  Y <- X %*% C + UU


  sim.params <- list(
    typeA="sparse",
    typeB="all",
    n = n,
    d = d,
    q = q,
    pvec = pvec,
    rvec = rvec,
    es = es,
    snr = snr,
    sigma = sigma,
    rho_X = rho_X,
    rho_E = rho_E
  )


  attr(B, "pvec") <- as.numeric(table(pvec))

  list(Y = Y, X = X, C = C, A = A, B = B, UU = UU.out, sigma = sigma, snr.out = SNR,
       Xsigma = Xsigma, n=n, d=d, q=q, pvec=pvec, rvec=rvec,
       params = sim.params)



}






#' @export ToOrthogonal
ToOrthogonal <- function(q, r, threshold=0.2, simplify=TRUE){

  # q <- 10; r <- 6
  L <- matrix(rnorm(q*r), q, r)


  orthogonalize <- function(X, threshold){
    fit.svd <- svd(X)
    orth <- fit.svd$u %*% t(fit.svd$v)
    X <- ifelse(abs(orth) > threshold, orth, 0)
  }



  it <- 0
  while(TRUE){
    it <- it+1

    L <- orthogonalize(L, threshold)



    if( norm(crossprod(L) - diag(1, r), "F") < 1e-15 ){
      # print(it)
      break
    }

    if(it > 1000) break
  }
  # varimax(L, normalize = FALSE)
  # GPArotation::Varimax(L)
  # crossprod(L)

  if( simplify ){
    # L <- as( GPArotation::cfQ(L, kappa = 0.5)$loadings, "matrix" )
    L <- as( GPArotation::quartimax(L)$loadings, "matrix" )
  }

  L
}




# Example of A
function(){
  q=8; r=2; nuA=0.2; simplify=F
  for(nuA in c(0, 0.2, 0.4, 1.0)){
    set.seed(123)
    A <- ToOrthogonal(q, r, 0, simplify=simplify)
    A[apply(A, 1, norm, "2")<nuA,] <- 0
    A <- GPArotation::quartimax( with( svd(A), tcrossprod(u,v) ) )$loadings
    print( round(A,3) )
  }
  
}





