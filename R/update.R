#' @import GPArotation
#' 
#' @export png.svd
png.svd <- function (x, nu = min(n, p), nv = min(n, p), LINPACK = TRUE){
  ## Store the svd function, but with LINPACK = T as default:
  # print("LINPACK:"); print(LINPACK)  ## added so you can see it's changed
  x <- as.matrix(x)
  if (any(!is.finite(x)))
    stop("infinite or missing values in 'x'")
  dx <- dim(x)
  n <- dx[1L]
  p <- dx[2L]
  if (!n || !p)
    stop("a dimension is zero")
  La.res <- La.svd(x, nu, nv)   ## your problem line
  res <- list(d = La.res$d)
  if (nu)
    res$u <- La.res$u
  if (nv) {
    if (is.complex(x))
      res$v <- Conj(t(La.res$vt))
    else res$v <- t(La.res$vt)
  }
  res
}




Update.Rotation <- function(A, rot.method, group=NULL){
  
  bapply <- function(X, block, f, ...){
    
    nrow <- length( unique( block[,1] ) )
    ncol <- length( unique( block[1,] ) )
    
    bl <- unique(as.vector(block))
    mb <- length(bl)
    out <- NULL
    for( i in 1:mb ){
      xij <- bl[i]
      row <- which(block==xij, arr.ind=TRUE)[,1]
      col <- which(block==xij, arr.ind=TRUE)[,2]
      
      out[i] <- f(X[row,col], ...)
    }
    
    matrix(out, nrow, ncol)
  }
  
  
  if( identical(A, diag(1, nrow(A), ncol(A))) ){
    
    return(  list( A=A, rotation=diag(ncol(A)) )  )
    
  }
  
  
  if( identical(A, matrix(0, nrow(A), ncol(A))) ){
    
    return(  list( A=A, rotation=diag(ncol(A)) )  )
    
  }
  
  
  if(!is.null(group)) A <- bapply(A, sapply(1:6, function(x) group+10*x), function(x) norm(x, "2") )
  
  fit.rot.A <- switch(rot.method,
                      "varimax" = png_varimax(A, maxit=100),
                      "quartimax" = png_quartimax(A, maxit=100),
                      "geomin" = geominQ(A, delta=0.5),
                      "oblimin" = oblimin(A, gam=0.5),
                      "oblimax" = oblimax(A),
                      "infomax" = infomaxT(A),
                      "CF" = cfT(A, kappa=0.5),
                      "McCammon" = mccammon(A),
                      "simplimax" = simplimax(A, k = 10),
                      "tandem" = tandemI(A),
                      "entropy" = entropy(A) )
  
  
  A <- fit.rot.A$loadings
  
  # check
  ## off-diagonal elements
  # sum( crossprod(A) %>% {.[abs(row(.)-col(.))>0]} ) < 1e-10
  ## diagonal elements
  # (sum( crossprod(A) %>% {.[abs(row(.)-col(.))==0]} ) - ncol(A)) < 1e-10
  
  list( A=fit.rot.A$loadings, rotation=fit.rot.A$Th )
}






png.quartimax <- function(X){
  
  Q <- diag(ncol(X))
  p <- nrow(X)
  
  Z <- X %*% Q
  dQ <- - Z^3
  G <- crossprod(X, dQ)
  while(TRUE){
    Qold <- Q
    
    Q <- with( svd( Q + 0.1*G ), tcrossprod(u, v) )
    
    Z <- X %*% Q
    dQ <- - Z^3
    # dQ <- 1/p * (Z^3 - Z %*% diag(drop(rep(1, p) %*% Z^2))/p)
    G <- crossprod(X, dQ)
    
    norm(Q - Qold, "F") |> print()
    
    if( norm(Q - Qold, "F") < 1e-10 ) break
  }
  list(loadings = X %*% Q, rotation = Q)
}


png.varimax <- function(X){
  Q <- diag(ncol(X))
  p <- nrow(X)
  
  Z <- X %*% Q
  dQ <- 1/p * (Z^3 - Z %*% diag(drop(rep(1, p) %*% Z^2))/p)
  G <- crossprod(X, dQ)
  while(TRUE){
    Qold <- Q
    
    Q <- with( svd( Q + 0.1*G ), tcrossprod(u, v) )
    
    Z <- X %*% Q
    dQ <- 1/p * (Z^3 - Z %*% diag(drop(rep(1, p) %*% Z^2))/p)
    G <- crossprod(X, dQ)
    
    if( norm(Q - Qold, "F") < 1e-10 ) break
  }
  list(loadings = X %*% Q, rotation = Q)
}











# UpdateA.Grad <- function(X, Y, A, B, pvec, nrank, maxit.B, eps.B, Lambda_ik, verbose=FALSE){
# 
# 
#   omega <- 0
#   n <- nrow(Y)
#   q <- ncol(Y)
# 
# 
#   error <- 1000
#   it <- 1
#   while((it < maxit.B) & (error > eps.B)){
#     it <- it + 1
#     Aold <- A
# 
#     for( k in 1:nrank ){
# 
#       for( h in 1:q ){
# 
#         Ahk <- A[h, k]
#         Yh <- Y[,h]
#         Ymh <- Y[,-h]
#         Amhk <- A[-h,k]
#         Bk <- B[,k]
# 
#         A[h, k] <- (omega*Ahk^2 - crossprod(Yh) )^(-1) * (  t(Yh) %*% Ymh %*% Amhk - t(X%*%Bk) %*% Yh  )
# 
#       }
#     }
# 
#     if( mean( (A - Aold)^2 ) < eps.B ){
#       if(verbose) cat( "Error = ", mean( (A - Aold)^2 ), "\n" )
#       break
#     }
#   }
# 
#   return( A )
# }






# UpdateB.SGM <- function(X, Y, A, B, pvec, nrank, d, Lambda_ik, maxit.B, eps.B){
# 
#   pik <- as.numeric( table(pvec) )
# 
#   error <- 1000
#   it <- 1
#   while((it < maxit.B) & (error > eps.B)){
#     it <- it + 1
#     Bold <- B
# 
#     for( k in 1:nrank ){
# 
#       for( i in 1:d ){
# 
#         lambda_ik <- Lambda_ik[i,k] * sqrt(pik[i])
# 
#         pi <- which( pvec == i )
#         pmi <- which( pvec != i )
# 
#         Xi <- X[,pi]
#         Xmi <- X[,pmi]
# 
#         if( length(pmi) != 0 ){
#           Bmi <- B[pmi,]
#           Rmik <- (Y %*% A - Xmi %*% Bmi)[,k]
#           Sik <- t(Xi) %*% Rmik / lambda_ik
#         } else {
#           Rmik <- (Y %*% A)[,k]
#           Sik <- t(Xi) %*% Rmik / lambda_ik
#         }
# 
# 
#         if( norm( Sik, "2" ) < 1 ){
# 
#           B[pi,k] <- 0
# 
#         } else {
# 
#           if( norm( B[pi,k], "2" ) < 1e-10 ){
# 
#             B[pi,k] <- 0
# 
#           } else {
# 
#             SXiXi <- t(Xi) %*% Xi + diag(length(pi))*lambda_ik / norm( B[pi,k], "2" )
#             SXiXi_inv <- ginv(SXiXi)
#             # SXiXi_inv <- tryCatch(ginv(SXiXi), error = function(e) solve(SXiXi + 0.01 * diag(nrow(SXiXi))))
# 
#             B[pi,k] <- SXiXi_inv %*% t(Xi) %*% Rmik
# 
#           }
#         }
#       }
#     }
# 
#     error = norm(B - Bold, "F")
# 
#     if( error < eps.B ){
# 
#       # if(verbose) cat( "Error = ", mean( (B - Bold)^2 ), "\n" )
#       break
#     }
#   }
# 
#   # print(it)
#   # print(error)
# 
#   return( B )
# }




# UpdateB.BMD <- function(X, Y, A, B, pvec, nrank, maxit.B, eps.B, Lambda_ik, verbose=FALSE){
# 
#   eta_i <- NULL
#   for( i in 1:d ){
#     pi <- which( pvec == i )
#     eta_i[i] <- max( eigen(crossprod(X[,pi]))$values )
#   }
# 
#   pik <- as.numeric( table(pvec) )
# 
#   error <- 0
#   for( k in 1:nrank ){
#     lambda_k <- Lambda_ik[,k]
# 
#     bk <- B[,k]
# 
#     it <- 0
#     while( (it < maxit.B) ){
#       it <- it + 1
# 
#       bk_old <- bk
# 
#       for( i in 1:d ){
# 
#         pi <- which( pvec == i )
#         pmi <- which( pvec != i )
# 
#         lambda_ik <- lambda_k[i] * pik[i]
# 
#         Xi <- X[,pi]
#         bik <- bk[pi]
#         resid <- Y %*% A[,k] - X %*% bk
# 
#         Ui <- crossprod(Xi, resid)/n
#         gamma_i <- (1 + 1e-6)*eta_i[i]
#         Ri <- ( Ui + gamma_i * bik )
# 
#         bk[pi] <- (1/gamma_i) * Ri * max(0, 1 - lambda_ik/norm(Ri, "2"))
# 
#       }
# 
#       error <- error + mean( (bk - bk_old)^2 )
# 
#       if( error < eps.B ) break
#     }
# 
#     B[,k] <- bk
# 
#   }
# 
#   if(verbose) cat( "Error = ", mean( (B - Bold)^2 ), "\n" )
# 
#   return(B)
# }













UpdateA.ManifoldOptim <- function(X, Y, A=NULL, B, tau, eps=1e-10, Max_Iteration=100){
  # X <- data$X; Y <- data$Y; B <- data$B
  # eps=1e-10; Max_Iteration=100
  # tau=5
  
  n <- nrow(X); p <- ncol(X); q <- ncol(Y); r <- ncol(B)
  
  
  # tA <- function(A) { A = matrix(A, q, r); A[apply(A,1,norm,"2")<sqrt(eps),] <- 0; with(svd(A), tcrossprod(u,v)) }
  tA <- function(A) { A = matrix(A, q, r) }
  f <- function(A) { A <- tA(A); mean( (Y - X %*% B %*% t(A))^2 ) +
    tau * sum( apply(A, 1, function(xx) { sqrt( sum(xx^2) + eps ) - eps }) )
  }
  dF <- function(A) { A <- tA(A); As <- t( apply(A, 1, function(a) rep( 1/sqrt( sum(a^2) + eps), length(a) ) ) ); - 1/(n*q) * 2 * t(Y - X %*% B %*% t(A)) %*% (X %*% B) + tau * A * As }
  
  # f <- function(A) { A <- tA(A); mean( (Y - X %*% B %*% t(A))^2 ) + 
  #   tau * sum( sqrt( sum(A^2) + eps ) - eps )
  # }
  # dF <- function(A) { A <- tA(A); As <- 1/sqrt( sum(A^2) + eps); - 1/(n*q) * 2 * t(Y - X %*% B %*% t(A)) %*% (X %*% B) + tau * A * As }
  
  
  library(ManifoldOptim)
  mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
  prob <- new(mod$RProblem, f, dF)
  
  if( !is.null(A) ){
    A0 <- as.numeric( A )
  } else {
    A0 <- as.numeric( GPArotation::quartimax( with( svd(t(Y) %*% X %*% B), tcrossprod(u,v) ) )$loadings )
  }
  
  # A0 <- as.numeric(GPArotation::quartimax(orthonorm(matrix(rnorm(q*r), nrow=q, ncol=r)))$loadings)
  # A0 <- as.numeric(orthonorm(matrix(rnorm(q*r), nrow=q, ncol=r)))
  mani.params <- get.manifold.params(IsCheckParams = FALSE)
  solver.params <- get.solver.params(IsCheckParams = FALSE,
                                     isconvex = FALSE, 
                                     Tolerance = 1e-4,
                                     Max_Iteration = Max_Iteration)
  mani.defn <- get.stiefel.defn(q, r)
  
  res <- manifold.optim(prob, mani.defn, method = c("RTRSD", "LRTRSR1")[2],
                        mani.params = mani.params, 
                        solver.params = solver.params, 
                        x0 = A0)
  
  
  
  # out.tmp <- NULL
  # for( i in 1:50 ){
  #   methods <- c("LRBFGS","LRTRSR1","RBFGS","RBroydenFamily","RCG","RNewton","RSD","RTRNewton","RTRSD","RTRSR1","RWRBFGS")[-c(6,8)]
  #   ff <- lapply(methods, function(mm){
  #     res <- manifold.optim(prob, mani.defn, method = mm,
  #                           mani.params = mani.params, 
  #                           solver.params = solver.params, 
  #                           x0 = A0)
  #     res
  #   })
  #   
  #   tmp <- sapply(ff, function(res) c(angle=png.angle(data$A, tA(res$xopt))$max * 180/pi, time=res$elapsed) )
  #   colnames(tmp) <- methods
  #   out.tmp[[i]] <- tmp
  #   print( round(Reduce("+", out.tmp) / i, 3) )
  # }
  # LRBFGS LRTRSR1  RBFGS RBroydenFamily    RCG    RSD  RTRSD RTRSR1 RWRBFGS
  # angle 19.505  18.084 19.599         19.807 18.056 17.021 10.207 18.107  21.220
  # time   0.109   0.104  0.113          0.110  0.221  0.111  0.101  0.106   0.112
  # LRTRSR1; RTRSD
  
  
  Anew <- tA(res$xopt)
  # Anew <- ifelse(abs(Anew)<eps, 0, Anew)
  # Anew <- with(svd(Anew), tcrossprod(u,v))
  
  list(res=res, A=Anew)
}




# UpdateA.OptManifold <- function(X, Y, A, B, Delta, alpha=1e+3){
#   library(ManifoldOptim)
# 
#   if(FALSE){
#     alpha <- 100000
#   }
# 
#   p <- ncol(X);  q <- ncol(Y);  r <- ncol(A)
#   a0 <- as.numeric(A)
# 
#   tx <- function(x) { matrix(x, q, r) }
# 
#   YtX <- crossprod(Y, X)
#   F <- function(a){
#     A <- tx(a)
# 
#     # {
#     #   sum( A^4 )
#     # }
# 
#     # g <- 1/4 * sum( A^4 )
#     g <- sum( apply( A, 2, function(Ak) mean(Ak^4) + (mean(Ak^2))^2 ) )
# 
# 
#     - sum(diag( t(Y) %*% X %*% B %*% t(A) )) + delta * sum( abs(A) ) - alpha * g
# 
# 
#     # {
#     #
#     # }
# 
#   }
# 
#   dF <- function(a){
#     A <- tx(a)
# 
#     # {
#     #   4*A^3
#     # }
# 
#     # - t(Y) %*% X %*% B + delta * sign(A) - alpha * A^3
# 
#     g <- 4 * A^3 - 2 * A
#     # g <- A^3
# 
#     As <- 1/sqrt(A^2 + 1e-5)
#     AAs <- A * As
#     - t(Y) %*% X %*% B + delta * AAs - alpha * g
# 
#     # {
#     #   4 * A^3 - 2 * A
#     # }
# 
#   }
# 
#   mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
#   prob <- new(mod$RProblem, F, dF)
# 
#   mani.params <- get.manifold.params(IsCheckParams = FALSE)
#   solver.params <- get.solver.params(IsCheckParams = FALSE,
#                                      isconvex = FALSE, Max_Iteration = 100)
#   mani.defn <- get.stiefel.defn(q, r)
# 
#   res <- manifold.optim(prob, mani.defn, method = c("RCG", "RNewton", "RSD")[3],
#                         mani.params = mani.params, solver.params = solver.params, x0 = a0)
#   # print(res)
#   (tx(res$xopt))
#   data$A
#   Varimax((tx(res$xopt)))
# 
#   1/(q*r)*norm( abs(tx(res$xopt)) - abs(data$A), "F" )
# 
#   if(FALSE){
#     (tx(res$xopt)) %>% round(5)
#     crossprod(tx(res$xopt))
#     GPArotation::Varimax((tx(res$xopt)))
# 
#     1/(q*r)*norm( abs(tx(res$xopt)) - abs(data$A), "F" )
#   }
# 
#   (tx(res$xopt))
# }




# UpdateA.OptManifold2 <- function(X, Y, A, B, Delta){
#   library(ManifoldOptim)
# 
#   p <- ncol(X);  q <- ncol(Y);  r <- ncol(A)
#   a0 <- as.numeric(A)
# 
#   tx <- function(x) { matrix(x, q, r) }
# 
# 
#   YtX <- crossprod(Y, X)
#   F <- function(a){
#     A <- tx(a)
# 
#     {
#       YA <- Y %*% A
#       XB <- X %*% B
#       0.5 * sum(diag( t(YA) %*% YA )) - sum(diag( t(XB) %*% YA )) + delta * sum( abs(A) )
#     }
# 
#     # {
#     #   # A B^T X^T Y
#     #   1/(n*r)*sum( diag( tcrossprod(tcrossprod(A, B), YtX) ) ) + sum( Delta * abs(A) )
#     # }
# 
#   }
# 
#   dF <- function(a){
#     A <- tx(a)
# 
#     {
#       YA <- Y %*% A
#       XB <- X %*% B
#       t(Y) %*% YA - t(Y) %*% XB + delta * sign(A)
#     }
# 
#     # {
#     #   As <- 1/sqrt(A^2 + 1e-5)
#     #   AAs <- A * As
#     #   1/(n*r)*- YtX %*% B + Delta * AAs
#     # }
# 
#   }
# 
#   mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
#   prob <- new(mod$RProblem, F, dF)
# 
#   mani.params <- get.manifold.params(IsCheckParams = FALSE)
#   solver.params <- get.solver.params(IsCheckParams = FALSE,
#                                      isconvex = FALSE, Max_Iteration = 100)
#   mani.defn <- get.stiefel.defn(q, r)
# 
#   res <- manifold.optim(prob, mani.defn, method = c("RCG", "RNewton", "RSD")[1],
#                         mani.params = mani.params, solver.params = solver.params, x0 = a0)
#   # print(res)
#   (tx(res$xopt))
# }


















# UpdateQ.1 <- function(X, Y, A, B, alpha, lambda, delta, pvec){
#   # X: n x p
#   # Y: n x q
#   # A: q x r
#   # B: p x r
#   # lambda: d x r
#   # delta: 1 x 1
# 
#   delta <- 1
#   alpha <- 0.1
# 
#   r <- ncol(A)
# 
#   if( all( B == 0 ) ) return( diag(r) )
# 
#   obj.Q <- function(A, B, Q, lambda, delta, pvec){
# 
#     A <- A %*% Q
#     B <- B %*% Q
# 
#     r <- ncol(A)
# 
#     b_ik.norm <- sapply( 1:r, function(k) tapply( B[,k], pvec, function(x) norm(as.matrix(x), "F") ) )
# 
#     sum(lambda*b_ik.norm) - delta/4 * sum(A^4)
#   }
# 
#   dF <- function(A, B, Q, lambda, delta, pvec){
# 
#     AQ <- A %*% Q
#     BQ <- B %*% Q
#     p <- nrow(B);  r <- ncol(A);  d <- max(pvec);  q <- nrow(A)
# 
# 
#     b_ik.norm <- sapply( 1:r, function(k) tapply( BQ[,k], pvec, function(x) max(1e-5, norm(as.matrix(x), "F") ) ) )
#     b_ik.sum <- sapply( 1:r, function(k) tapply( BQ[,k], pvec, function(x) sum(as.matrix(x)) ) )
# 
# 
#     # Grad.B <- t(b_ik.norm) %*% (lambda*b_ik.sum/b_ik.norm)
#     Grad.B <- - t(B) %*% apply( lambda*b_ik.sum/b_ik.norm, 2, function(x) as.numeric(rep(x, each=p/d)) )
# 
#     # varimax
#     Grad.A <-  t(A) %*% ( 1/q*(AQ^3 - AQ %*% diag(drop(rep(1, q) %*% AQ^2))/q) )
# 
#     # quartimax
#     # Grad.A <- - delta * t(A) %*% AQ^3
# 
#     Grad.A + Grad.B
#   }
# 
# 
#   Q <- diag(r)
#   for( i in 1:100 ){
#     Qold <- Q
# 
#     Q <- with( png.svd( Q + alpha * dF(A=A, B=B, Q=Q, lambda=lambda, delta=delta, pvec=pvec) ), tcrossprod(u,v) )
# 
#     # print( norm( Qold - Q, "F" ) )
# 
#     print( obj.Q(A=A, B=B, Q=Q, lambda=lambda, delta=delta, pvec=pvec) )
# 
#     if( norm( Qold - Q, "F" ) < 1e-10 ) break
#   }
# 
#   # (A %*% Q)|> (function(.) ifelse( abs(.) > 1e-10, ., 0) )()
#   #
# 
#   Q
# }



# UpdateQ.2 <- function(X, Y, A, B, alpha, lambda, delta, pvec){
#   # X: n x p
#   # Y: n x q
#   # A: q x r
#   # B: p x r
#   # lambda: d x r
#   # delta: 1 x 1
# 
#   # delta <- 1
#   # alpha <- 0.1
# 
#   r <- ncol(A)
# 
#   if( all( B == 0 ) ) return( diag(r) )
# 
#   obj.Q <- function(X, Y, A, B, Q, lambda, delta, pvec){
# 
#     A <- A %*% Q
#     B <- B %*% Q
# 
#     r <- ncol(A)
# 
#     1/2 * norm( Y - X%*%B%*%t(A%*%Q), "F" ) - delta/4 * sum(A^4)
#   }
# 
#   dF <- function(X, Y, A, B, Q, lambda, delta, pvec){
# 
#     AQ <- A %*% Q
#     BQ <- B %*% Q
#     p <- nrow(B);  r <- ncol(A);  d <- max(pvec);  q <- nrow(A)
# 
# 
#     Grad.Y <- - t(Y %*% A) %*% (X %*% B)
# 
#     # varimax
#     # Grad.A <- - t(A) %*% ( 1/q*(AQ^3 - AQ %*% diag(drop(rep(1, q) %*% AQ^2))/q) )
# 
#     # quartimax
#     Grad.A <- - delta * t(A) %*% AQ^3
# 
#     Grad.A + Grad.Y
#   }
# 
# 
#   Q <- diag(r)
#   for( i in 1:100 ){
#     Qold <- Q
# 
#     Q <- with( png.svd( Q - alpha * dF(X=X, Y=Y, A=A, B=B, Q=Q, lambda=lambda, delta=delta, pvec=pvec) ), tcrossprod(u,v) )
# 
#     # print( norm( Qold - Q, "F" ) )
# 
#     # print( obj.Q(X=X, Y=Y, A=A, B=B, Q=Q, lambda=lambda, delta=delta, pvec=pvec) )
# 
#     if( norm( Qold - Q, "F" ) < 1e-10 ) break
#   }
# 
#   # (A %*% Q)|> (function(.) ifelse( abs(.) > 1e-10, ., 0) )()
#   #
# 
#   Q
# }







# A -----------------------------------------------------------------------


# UpdateA.1 <- function(X, Y, A, B, alpha, lambda, delta, pvec){
#   # delta <- 1
#   # alpha <- 10
# 
#   n <- nrow(Y); p <- ncol(X); q <- ncol(Y); r <- ncol(A)
# 
#   obj <- function(X, Y, A, B, delta){
#     Syx <- crossprod(Y, X)
# 
#     -sum(diag(Syx %*% tcrossprod(B, A))) - delta * sum(abs(A))/4
# 
#     norm(Y - X%*%B%*%t(A), "F") - delta * sum(A^4)/4
#   }
# 
#   dF <- function(X, Y, A, B, delta){
#     Syx <- crossprod(Y, X)
#     - Syx %*% B - delta * A^3
#   }
# 
#   # A <- A0
#   for( i in 1:100 ){
#     Aold <- A
# 
#     A <- with( png.svd( A - alpha * dF(X, Y, A, B, delta) ), tcrossprod(u,v) )
# 
#     # print( norm( Aold - A, "F" ) )
#     #
#     # print(obj(X, Y, A, B, delta))
# 
#     if( norm( Aold - A, "F" ) < 1e-10 ) break
#   }
# 
# 
#   # A |> (function(.) ifelse( abs(.) > 1e-10, ., 0) )()
#   # A0 |> (function(.) ifelse( abs(.) > 1e-10, ., 0) )()
# 
#   # norm(Y - X%*%B%*%t(A), "F")
#   # norm(Y - X%*%B%*%t(A0), "F")
#   #
#   # obj(X, Y, A, B, delta)
#   # obj(X, Y, A0, B, delta)
# 
#   A
# }


# UpdateA.stiefel1 <- function(X, Y, A, B, alpha=1, delta=1){
#   if(FALSE){
#     B <- data$B
#   }
# 
#   q <- ncol(Y)
#   r <- ncol(A)
# 
#   epsilon <- 1e-4
# 
#   # A <- data$A %>% {. + matrix(rnorm(prod(dim(.))), nrow(.), ncol(.))*0.1} %>% {svd(.)$u}
# 
#   # delta <- n*p*q*1000 #200000
# 
#   YtX <- crossprod(Y, X)
# 
#   f <- function(A){
#     As <- 1/sqrt(A^2 + 1e-50)
#     -sum(diag( t(Y) %*% X %*% B %*% t(A) )) + delta * sum( As )
#   }
# 
#   df <- function(A){
#     As <- 1/sqrt(A^2 + 1e-50)
#     - t(Y) %*% X %*% B + delta * A * As
#   }
# 
#   Anew <- optStiefel(f, df, Vinit=A,
#                                # method="curvilinear",
#                                # method="bb",
#                                # searchParams=list(rho=0.5, eta=0.1),
#                                # searchParams=list(rho1=0.1, rho2=0.9, tau=0.1),
#                                # maxIters = 1e+4,
#                                # maxLineSearchIters = 1e+4,
#                                tol=1e-15)
# 
#   ifelse( abs(Anew) < sqrt(1e-40), 0, Anew )
# 
# }



# UpdateA.stiefel2 <- function(X, Y, A, B, delta=1){
#   if(FALSE){
#     B <- data$B
#   }
# 
#   q <- ncol(Y)
#   r <- ncol(A)
# 
#   epsilon <- 1e-4
# 
#   # A <- data$A %>% {. + matrix(rnorm(prod(dim(.))), nrow(.), ncol(.))*0.1} %>% {svd(.)$u}
# 
#   # delta <- n*p*q*1000 #200000
# 
#   YtX <- crossprod(Y, X)
# 
#   f <- function(A){
#     -sum(diag( t(Y) %*% X %*% B %*% t(A) )) - delta * sum( A^4 )/4
#   }
# 
#   df <- function(A){
#     - t(Y) %*% X %*% B - delta * A^3
#   }
# 
#   Anew <- optStiefel(f, df, Vinit=A,
#                                # method="curvilinear",
#                                # method="bb",
#                                # searchParams=list(rho=0.5, eta=0.1),
#                                # searchParams=list(rho1=0.1, rho2=0.9, tau=0.1),
#                                # maxIters = 1e+4,
#                                # maxLineSearchIters = 1e+4,
#                                tol=1e-15)
# 
#   Anew
# 
# }




# UpdateA.stiefel2 <- function(X, Y, A, B, delta, alpha=1){
#   library(ManifoldOptim)
#
#   if(FALSE){
#     alpha <- 100000
#   }
#
#   p <- ncol(X);  q <- ncol(Y);  r <- ncol(A)
#   a0 <- as.numeric(A)
#
#   tx <- function(x) { matrix(x, q, r) }
#
#   f <- function(a){
#     A <- tx(a)
#     - sum(diag( t(Y) %*% X %*% B %*% t(A) )) + delta * sum( abs(A) )
#   }
#
#   df <- function(a){
#     A <- tx(a)
#     As <- 1/sqrt(A^2 + 1e-35)
#     AAs <- A * As
#     - t(Y) %*% X %*% B + delta * AAs
#   }
#
#   mod <- Module("ManifoldOptim_module", PACKAGE = "ManifoldOptim")
#   prob <- new(mod$RProblem, f, df)
#
#   mani.params <- get.manifold.params(IsCheckParams = FALSE)
#   solver.params <- get.solver.params(IsCheckParams = FALSE,
#                                      isconvex = TRUE, Max_Iteration = 100)
#   mani.defn <- get.stiefel.defn(q, r)
#
#   res <- manifold.optim(prob, mani.defn, method = c("RCG", "RNewton", "RSD")[3],
#                         mani.params = mani.params, solver.params = solver.params, x0 = a0)
#   # print(res)
#   (tx(res$xopt))
# }








