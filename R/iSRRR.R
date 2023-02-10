#' @name iSRRR
#' @title iSRRR: Integrative Sparse Reduced Rank Regression
#' @description
#' The iSRRR solves the constrained optimization problem with sparse regularization:
#' 
#' \eqn{ \frac{1}{2n} \| Y - X B A^T \|_F^2 +  \lambda \sum \sum \| b_{ik} \|_2}  subject to \eqn{A \in \mathcal{T}(v), ~ A \in \mathcal{O}_s(q,r) }
#'
#' @param X An \code{n} by \code{p} design matrix.
#' @param Y An \code{n} by \code{q} response matrix.
#' @param pvec A numeric vector of the number of variable for each dataset.
#'
#' @param nrank The rank of matrices to be estimated.
#'
#' @param cutoff Hard-thresholding parameter. Default is 0 (no thresholding).
#'
#' @param params A list with \code{lambda.seq}=NULL, \code{nlambda}=5, \code{lambda.factor}=1e-4, \code{log.scale}=TRUE, and \code{group.size}=FALSE.
#'
#' @param control A list with \code{best}=FALSE, \code{rot.method}="quartimax" \code{maxit.B}=3e8, \code{eps.B}=1e-8, \code{maxit.mse}=50, \code{eps.mse}=1e-6, \code{X.scale}=c("group", "each"), \code{Y.scale}=FALSE, \code{verbose}=FALSE, and \code{threads}=1.
#' @param trueB A list contains X with dimension \code{n} by \code{p} and Y with dimension \code{n} by \code{q}.
#' 
#' @param use.gglasso Default is TRUE.
#' 
#' 
#' @return A list with output objects
#' 
#' 
#' @author Kipoong Kim <statpng@snu.ac.kr>, Sungkyu Jung <sungkyu@snu.ac.kr>
#' 
#' @references
#' No reference.
#'
### #' @keywords "Data integration" "Multi-source data" "Orthogonal rotation" "Reduced-rank regression" "Structural learning"
#' @examples 
#' 
#' if(FALSE){
#' 
#' fit <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec,
#' nrank = data$nrank,
#' nu = c(0:4*0.1),
#' params=list(nlambda=10, lambda.factor=1e-4),
#' control=list(best=TRUE, maxit.mse=20, eps.mse=1e-8,
             #' maxit.B=6e6, eps.B=1e-6),
#' trueA = data$A, trueB = data$B )
#' }
#'
#' @importFrom gglasso gglasso
#' @import Rcpp
#' @importFrom MASS ginv
#' @useDynLib iSRRR
#'
#' @export arcov
#' 
#' @export iSRRR
iSRRR <- function(X, Y, pvec, nrank, nu, params=NULL, control=NULL, trueA=NULL, trueB=NULL, methodA="hard"){

  if(FALSE){
    fit <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec,
                  nrank = data$nrank,
                  nu = c(0:4*0.1),
                  methodA="hard",
                  params=list(nlambda=10, lambda.factor=1e-5),
                  control=list(best=TRUE, maxit.mse=200, eps.mse=1e-12,
                  maxit.B=3e6, eps.B=1e-8),
                  trueA = data$A, trueB = data$B )
  }
  
  
  if(FALSE){
    X = data$X; Y = data$Y; pvec = data$pvec;
    nrank = data$nrank;
    nu = (0:4*0.1);
    params=list(nlambda=10, lambda.factor=1e-4);
    control=list(best=TRUE, maxit.mse=100, eps.mse=1e-12, 
                 maxit.B=3e6, eps.B=1e-6);
    trueA = data$A; trueB = data$B
    methodA <- "hard"
    
    X = X
    Y = Y
    pvec = pvec
    nrank = 2
    nu = c(0:4*0.1)[4]
    methodA="hard"
    params=list(nlambda=10, lambda.factor=1e-4)
    control=list(best=TRUE, maxit.mse=200, eps.mse=1e-10,
                 maxit.B=1e4, eps.B=1e-4)
    trueA=NULL
    trueB=NULL
    
    
    
    X = data$X
    Y = data$Y
    pvec = data$pvec
    nrank = data$nrank
    nu = c(0,0.1,0.2)[1]
    params=list(nlambda=20, lambda.factor=1e-4)
    control=list(best=TRUE, maxit.mse=200, eps.mse=1e-10, 
                 maxit.B = 3e8, eps.B = 1e-10) 
    trueA=data$A
    trueB=data$B
  }
  
  
  
  if( length(pvec) != ncol(X) ) stop("pvec == ncol(X)")
  
  
  
  fit.default <- default.params()
  params0 <- fit.default$params
  control0 <- fit.default$control

  params <- append(params, params0[ which( !names(params0) %in% names(params) ) ])
  control <- append(control, control0[ which( !names(control0) %in% names(control) ) ])




  {

    lambda.seq <- params$lambda.seq
    nlambda <- params$nlambda
    lambda.factor <- params$lambda.factor
    log.scale <- params$log.scale
    group.size <- params$group.size

    best <- control$best
    rot.method <- control$rot.method
    maxit.B <- control$maxit.B
    eps.B <- control$eps.B
    maxit.mse <- control$maxit.mse
    eps.mse <- control$eps.mse
    X.scale <- control$X.scale
    Y.scale <- control$Y.scale
    use.gglasso <- control$use.gglasso
    verbose <- control$verbose
    threads <- control$threads

  }



  n <- nrow(Y);  p <- length(pvec);
  d <- max(pvec);  q <- ncol(Y)


  
  fit.init <- get.init(X=X, Y=Y, pvec=pvec, nrank=nrank, lambda.seq=lambda.seq, lambda.factor=lambda.factor, nlambda=nlambda, log.scale=log.scale, group.size=group.size)


  A0 <- fit.init$A0
  B0 <- fit.init$B0

  params$lambda.seq <- fit.init$lambda.seq
  params$nlambda <- fit.init$nlambda

  out.final <- NULL
  for( inu in 1:length(nu) ){
    
    nu_i <- nu[inu]
    
    start <- proc.time()
    out <- fit.iSRRR(
      X = X,
      Y = Y,
      A0 = A0,
      B0 = B0,
      pvec = pvec,
      nrank = nrank,
      nu = nu_i,
      methodA = methodA,
      params = params,
      control = control,
      trueB = trueB, 
      use.gglasso=use.gglasso
    )
    end <- proc.time()
    
    print(paste0(inu, " / ", length(nu), ":"))
    print(end-start)
    
    
    SSE.list <- lapply( out, function(x) mean( (Y-X%*%x$C)^2 ) )
    attr(SSE.list, "npqdr") <- c(n=nrow(X), p=ncol(X), q=ncol(Y), d=max(pvec), r=nrank)
    
    
    out.final[[inu]] <- list(
      # X = X, Y = Y,
      pvec=pvec, nrank=nrank,
      A = lapply( out, function(x) x$A ),
      B = lapply( out, function(x) x$B ),
      SSE = SSE.list,
      mse.path = lapply( out, function(x) x$mse.path ),
      diff.mse = lapply( out, function(x) x$diff ),
      simplexity = lapply( out, function(x) x$simplexity ),
      iter.mse = lapply( out, function(x) x$it.mse ),
      nu = nu_i,
      params = params,
      control = control,
      trueA = trueA, trueB = trueB, 
      time.lambda = lapply( out, function(x) x$time.lambda ),
      time = end-start )
    
  }
    
  out.final

}




#' @export get.init
get.init <- function(X, Y, pvec, nrank, lambda.seq, lambda.factor, nlambda, log.scale, group.size){
  d <- max(pvec)

  C.init <- crossprod(with( png.svd(X), u %*% diag(1/d) %*% t(v) ), Y)

  ## equivalently
  # Sxx <- crossprod(X)
  # Syx <- crossprod(Y, X)
  # Sxx_inv <- MASS::ginv(Sxx)
  # C.init <- tcrossprod(Sxx_inv, Syx)

  coefSVD <- png.svd(C.init)
  coefSVD$u <- coefSVD$u[, 1:nrank, drop = FALSE]
  coefSVD$d <- coefSVD$d[1:nrank]
  coefSVD$v <- coefSVD$v[, 1:nrank, drop = FALSE]

  A0 <- coefSVD$v
  B0 <- with(coefSVD, u %*% diag(d, nrank, nrank))


  if( is.null(lambda.seq) ){
    lambda.seq <- lambda_sequence(X, Y, A0, pvec=pvec, nlambda=nlambda, lambda.factor = lambda.factor, log.scale=log.scale, group.size=group.size)
  }

  nlambda <- length(lambda.seq)


  attr(B0, "pvec") <- as.numeric(table(pvec))

  list(A0=A0, B0=B0, lambda.seq=lambda.seq, nlambda=nlambda)
}




#' @export default.params
default.params <- function(){


  params0 <- list(
    lambda.seq=NULL,
    nlambda=10,
    lambda.factor=1e-3, # 0.0001
    log.scale=TRUE,
    group.size=TRUE
  )

  control0 = list(
    best = FALSE,
    rot.method = "quartimax", # c("varimax", "quartimax", "oblimin", "oblimax", "infomax", "CF", "McCammon", "simplimax", "tandem", "entropy"),
    maxit.B = 1e6,
    eps.B = 1e-6,
    maxit.mse = 20,
    eps.mse = 1e-6,
    X.scale = FALSE, # c("group", "each"),
    Y.scale = FALSE, # TRUE, FALSE
    use.gglasso=TRUE,
    verbose=FALSE,
    threads=1
  )

  list(params=params0, control=control0)
}







fit.iSRRR <- function(X,Y,A0,B0,pvec,nrank,nu,methodA,params,control,trueB=NULL, use.gglasso=TRUE){


  n <- nrow(Y);  p <- length(pvec);
  d <- max(pvec);  q <- ncol(Y)
  r <- nrank

  {

    lambda.seq <- params$lambda.seq
    nlambda <- params$nlambda
    lambda.factor <- params$lambda.factor
    log.scale <- params$log.scale
    group.size <- params$group.size


    best <- control$best
    rot.method <- control$rot.method
    maxit.B <- control$maxit.B
    eps.B <- control$eps.B
    maxit.mse <- control$maxit.mse
    eps.mse <- control$eps.mse
    X.scale <- control$X.scale
    Y.scale <- control$Y.scale
    verbose <- control$verbose
    threads <- control$threads

  }




  pik <- replicate(nrank, as.numeric(table(pvec)))
  nlambda <- length(lambda.seq)




  # iSRRR -----------------------------------------------------------------
  
  out <- NULL
  for( ll in 1:nlambda ){
    
    start <- proc.time()
    
    lambda <- lambda.seq[ll]
    # if(group.size){
    #   lambda <- lambda * sqrt(pik)[,1]
    # }
    
    it.mse <- 1
    
    mse.path <- matrix(NA, nrow=maxit.mse+1, ncol=2)
    colnames(mse.path) <- c("mse", "mse+B")
    diff <- matrix(NA, nrow=maxit.mse+1, ncol=4)
    colnames(diff) <- c("Y", "A", "B", "C")
    diff[1,] <- eps.mse*2
    
    simplexity <- matrix(NA, nrow=maxit.mse+1, ncol=2)
    colnames(simplexity) <- c("quartimax", "varimax")
    
    minimum <- 9e+5
    min.count <- 0
    while((it.mse < maxit.mse+1)&( diff[it.mse,"C"] > eps.mse)){
      # print(it.mse)
      it.mse = it.mse+1
      
      if(it.mse == 2){
        A <- A0
        B <- B0
      }
      
      
      # Old A, B ----------------------------------------------------------------
      {
        Aold <- A
        Bold <- B
      }
      
      
      
      
      if(methodA == "hard"){
        
        A <- with( png.svd( crossprod(Y, X) %*% B ), tcrossprod( u, v ) )
        # A <- ifelse( abs(A) > cutoff, A, 0 )
        A <- t( apply(A, 1, function(a){
          if( all( norm(a, "2") < nu ) ){
            rep(0, length(a))
          } else {
            a
          }
        }) )
        A <- with( png.svd(A), tcrossprod(u,v) )
        
      } else if(methodA == "RSD"){
        
        if( all( B == 0 ) ){
          A <- rbind( diag(r), matrix(0, q-r, r))
        } else {
          A <- UpdateA.ManifoldOptim(X, Y, A, B, nu, Max_Iteration = 50)$A
        }
        
      }
      
      
      Q <- Update.Rotation(A, rot.method = rot.method)$rotation
      
      A <- A %*% Q
      B <- B %*% Q
      
      
      
      
      if(verbose) print( paste0("(1) ", mean( (tcrossprod(B,A) - Cold)^2 )) )
      
      
      
      if(use.gglasso){
        
        fit.gglasso <- NULL
        for(k in 1:nrank){
          
          if(group.size){
            fit <- try( gglasso(X, Y %*% A[,k], group = pvec, lambda = lambda, loss = "ls", eps = eps.B, maxit = maxit.B, pf = sqrt(table(pvec))), silent = TRUE )
          } else {
            fit <- try( gglasso(X, Y %*% A[,k], group = pvec, lambda = lambda, loss = "ls", eps = eps.B, maxit = maxit.B, pf = rep(1, length(table(pvec))) ), silent = TRUE )
          }
          
          
          cc = 0
          maxit.B2 <- maxit.B
          while( class(fit)[1] == "try-error" ){
            cc = cc + 1
            maxit.B2 = maxit.B2 + 1e5
            
            if(group.size){
              fit <- try( gglasso(X, Y %*% A[,k], group = pvec, lambda = lambda, loss = "ls", eps = eps.B, maxit = maxit.B2, pf = sqrt(table(pvec))), silent = TRUE )
            } else {
              fit <- try( gglasso(X, Y %*% A[,k], group = pvec, lambda = lambda, loss = "ls", eps = eps.B, maxit = maxit.B2, pf = rep(1, length(table(pvec))) ), silent = TRUE )
            }
            
            if( cc == 100 ) stop("check the gglasso part or increase the maxit.B !");
          }
          
          fit.gglasso[[k]] <- fit
        }
        B <- do.call( "cbind", lapply(fit.gglasso, function(x) x$beta) )
        
      } else {
        
        fit.B <- UpdateB_BMD(X = X, Y = Y, A = A, B = B,
                             pvec = pvec, nrank = nrank, d = d,
                             lam = lambda, maxit = maxit.B, eps = eps.B, threads = threads)
        
        B <- fit.B$B
        
      }
      
      
      if(verbose) print( paste0("(2) ", mean( (tcrossprod(B,A) - Cold)^2 )) )
      
      
      
      B.norm <- function(B, pvec, group.size=FALSE){
        
        d <- length( table(pvec) )
        
        out <- NULL
        for( i in 1:d){
          out1 <- NULL
          for( k in 1:ncol(B) ){
            if(group.size){
              out1 <- cbind( out1, norm( B[which( pvec == i ),k], "2" ) / sqrt(sum(pvec == i)) )
            } else {
              out1 <- cbind( out1, norm( B[which( pvec == i ),k], "2" ))
            }
            
          }
          out <- rbind( out, out1 )
        }
        
        out
      }
      
      
      
      {
        C <- tcrossprod(B, A)
        Cold <- tcrossprod(Bold, Aold)
        
        simplexity[it.mse, "quartimax"] <- crit.quartimax(A)
        simplexity[it.mse, "varimax"] <- crit.varimax(A)
        
        mse.path[it.mse, "mse"] <- mean( Y - X %*% tcrossprod(B, A) )
        mse.path[it.mse, "mse+B"] <- mean( Y - X %*% tcrossprod(B, A) ) + lambda * sum( B.norm(B, pvec) )
        
        diff[it.mse, "Y"] <- mean( (X %*% C - X %*% Cold)^2 )
        diff[it.mse, "A"] <- mean( (A - Aold)^2 )
        diff[it.mse, "B"] <- mean( (B - Bold)^2 )
        diff[it.mse, "C"] <- mean( (C - Cold)^2 )
      }
      
      
    }
    
    end <- proc.time()
    
    
    
    
    if(it.mse == maxit.mse){
      warning(paste("Algorithm didn't converge in ", it.mse, " iterations at lambda[", ll,"] !", sep = ""))
    }
    
    simplexity <- simplexity[1:it.mse,]
    diff <- diff[1:it.mse,][-c(1),]
    mse.path <- mse.path[1:it.mse,][-c(1),]
    
    if(best & !is.null(trueB) ){
      
      if( any( dim(B) != dim(trueB) ) ) stop("Dimension of B must be the same.")
      
      best.cols <- rmse.best(trueB, B, best = TRUE)$cols
      A <- A[,best.cols]
      B <- B[,best.cols]
      
    }
    
    
    attr(B, "pvec") <- as.numeric(table(pvec))
    
    
    
    out[[ll]] <- list(A=A, B=B, C=C,
                      mse.path=mse.path,
                      diff=diff,
                      simplexity=simplexity,
                      it.mse=it.mse,
                      time.lambda=end-start)
  }
  
  

  out

}







##### #' @export cv.iSRRR
# cv.iSRRR <- function(fit, nrank=NULL, nfold=5, foldid=NULL){
# 
#   if(FALSE){
# 
#     library(iSRRR)
# 
#     set.seed(1)
#     data <- simdata(typeA="sparse0.2", typeB="partial", n=50, d=3, q=10, pvec=rep(50,3), rvec=rep(2,3), es="3", simplify=TRUE, snr=1.0)
# 
# 
#     fit <- with(data, iSRRR( X = X, Y = Y, pvec = pvec, nrank = 6, cutoff = 0.2) )
# 
#     nfold <- 5
#     foldid <- NULL
# 
#   }
# 
# 
# 
#   X <- fit$X;  Y <- fit$Y
#   pvec <- fit$pvec
#   cutoff <- fit$cutoff
#   params <- fit$params;  control <- fit$control
#   lambda.seq <- fit$params$lambda.seq
# 
#   n <- nrow(Y);  q <- ncol(Y);
#   p <- ncol(X);  d <- max(pvec);
# 
# 
#   if(is.null(nrank)){
#     nrank <- 1:min(c(n,p,q))
#   }
# 
#   nlambda <- length(lambda.seq)
#   n.nrank <- length(nrank)
# 
# 
# 
# 
#   if( is.null(foldid) ){
#     foldid <- sample( rep_len(1:nfold, n) )
#   }
#   if( max(foldid) != nfold ) stop("check your foldid and nfold.")
# 
# 
#   cv.err <- array(NA, dim=c(nfold, nlambda, n.nrank))
#   for( i in 1:nfold ){
# 
# 
#     id.tr <- which(foldid!=i);  id.ts <- which(foldid==i)
#     X.tr <- X[id.tr,];  Y.tr <- Y[id.tr,]
#     X.ts <- X[id.ts,];  Y.ts <- Y[id.ts,]
# 
#     for( k in 1:n.nrank ){
#       fit.tr <- iSRRR(X=X.tr, Y=Y.tr, pvec=pvec, nrank=nrank[k], cutoff=cutoff, params=params, control=control, trueB=NULL)
# 
#       for( j in 1:nlambda ){
#         cv.err[i,j,k] <- norm( Y.ts - X.ts %*% fit.tr$C[[j]], "F" )^2
#       }
#     }
# 
#   }
# 
# 
#   cvm <- apply( cv.err, 2:3, mean )
#   cvsd <- apply( cv.err, 2:3, sd )
#   cvup <- cvm + cvsd
#   cvlo <- cvm - cvsd
# 
#   wh.min <- which( cvm == min(cvm), arr.ind=TRUE )
# 
#   wh.rank.min <- wh.min[,2]
# 
#   wh.lambda.min <- wh.min[,1]
#   wh.lambda.1se <- min( which( cvm[,wh.rank.min] <= cvup[wh.lambda.min, wh.rank.min] ) )
# 
#   rank.opt <- nrank[wh.rank.min]
#   lambda.min <- lambda.seq[wh.lambda.min]
#   lambda.1se <- lambda.seq[wh.lambda.1se]
# 
#   list(fit=fit, nfold=nfold, foldid=foldid, cvm=cvm, cvsd=cvsd, cvup=cvup, cvlo=cvlo, lambda.min=lambda.min, lambda.1se=lambda.1se, nrank=rank.opt, wh.min=wh.lambda.min, wh.1se=wh.lambda.1se, wh.nrank=wh.rank.min)
# }





















#' @export png.iSRRR.lambda
png.iSRRR.lambda <- function(fit, idx){
  fit$A <- fit$A[[idx]]
  fit$B <- fit$B[[idx]]
  fit$SSE <- fit$SSE[[idx]]
  fit$diff.mse <- fit$diff.mse[[idx]]
  fit$simplexity <- fit$simplexity[[idx]]
  fit$iter.mse <- fit$iter.mse[[idx]]
  fit$time.lambda <- fit$time.lambda[[idx]]

  fit
}






#' @export iSRRR.boot.selfreq
iSRRR.boot.selfreq <- function(fit.boot, wh.nu){
  
  d <- fit.boot[[1]][[1]]$pvec %>% max
  nlambda <- fit.boot[[1]][[1]]$params$nlambda
  
  
  Structure <- lapply(1:nlambda, function(x) rbind( 0, png.GetStructure(d) ))
  for( i in 1:length(fit.boot)){
    
    a <- lapply( fit.boot[[i]][[wh.nu]]$B, function(x) x %>% print.B2() %>% {ifelse(abs(.)>0,1,0)} %>% {t(t(.) %>% {.[!duplicated(.),,drop=F]})} )
    
    for( h in 1:nlambda ){
      
      aa <- a[[h]]
      
      for( k in 1:ncol(aa) ){
        wh <- which( apply( Structure[[h]][-1,], 2, function(x) all(x==aa[,k]) ) )
        Structure[[h]][1,wh] <- Structure[[h]][1,wh] + 1
      }
      
      
    }
    
  }
  
  
  lapply(Structure, function(STR){
    as.data.frame( t( rbind(STR[-1,], STR[1,]) ) )
  })
  
  
  
}





#' @export iSRRR.boot
iSRRR.boot <- function(fit, X, Y, nboot=100){
  
  lambda.seq <- fit[[1]]$params$lambda.seq
  nrank <- fit[[1]]$nrank
  pvec <- fit[[1]]$pvec
  nu <- sapply(fit, function(x) x$nu)
  params <- fit[[1]]$params
  control <- fit[[1]]$control
  
  fit.boot <- NULL
  for( i in 1:nboot ){
    print(i)
    idx <- sample(1:nrow(X), replace=TRUE)
    fit.new <- iSRRR( X = X[idx,], Y = Y[idx,], pvec = pvec,
                  nrank = nrank,
                  nu = nu,
                  params=params,
                  control=control )
    
    fit.boot[[i]] <- fit.new
  }
  
  invisible( fit.boot )
  
}






#' @export iSRRR.boot.selarray
iSRRR.boot.selarray <- function(fit.boot){
  
  n.nu <- length(fit.boot[[1]])
  d <- fit.boot[[1]][[1]]$pvec %>% max
  
  out.BIC <- lapply( fit.boot, function(x) x %>% png.BICmat(print.out=F) %>% {attr(., "wh.min")} %>% apply(1, function(y) paste0(y, collapse=",") ) )
  
  
  Structure <- png.GetStructure(d)
  
  sel.array <- lapply(seq_len(n.nu), function(inu){
    fit.selfreq <- iSRRR.boot.selfreq(fit.boot, wh.nu=inu)
    mat <- do.call("rbind", lapply(fit.selfreq, function(x) round(x[,ncol(x)], 5) ))
  }) %>% simplify2array()
  
  dimnames(sel.array) <- list(paste0("lam.", 1:dim(sel.array)[1]),
                              paste0("S.", 1:dim(sel.array)[2]),
                              paste0("nu.", 1:dim(sel.array)[3])
  )
  
  attr(sel.array, "structure") <- Structure
  attr(sel.array, "wh.bic") <- table(unlist(out.BIC))
  
  sel.array
  
}





#' @export iSRRR.best
iSRRR.best <- function(fit.tmp, type="acc.C"){
  # trueA <- fit.all$train$data3$A
  # trueB <- fit.all$train$data3$B
  trueA <- fit.tmp[[1]]$trueA
  trueB <- fit.tmp[[1]]$trueB
  trueC <- tcrossprod(trueB,trueA)
  
  n.nu <- length(fit.tmp)
  n.lam <- fit.tmp[[1]]$params$nlambda
  
  acc.A.mat <- matrix(NA, nrow=n.lam, ncol=n.nu)
  acc.B.mat <- matrix(NA, nrow=n.lam, ncol=n.nu)
  acc.C.mat <- matrix(NA, nrow=n.lam, ncol=n.nu)
  rmse.B.mat <- matrix(NA, nrow=n.lam, ncol=n.nu)
  for( inu in 1:n.nu ){
    A.seq <- fit.tmp[[inu]]$A
    B.seq <- fit.tmp[[inu]]$B
    for( ilam in 1:n.lam ){
      Ahat <- fit.tmp[[inu]]$A[[ilam]]
      Bhat <- fit.tmp[[inu]]$B[[ilam]]
      acc.A.mat[ilam,inu] <- acc(Ahat, trueA, perm=F)
      acc.B.mat[ilam,inu] <- acc(Bhat, trueB, perm=F)
      acc.C.mat[ilam,inu] <- acc(tcrossprod(Bhat, Ahat), trueC, perm=F)
      rmse.B.mat[ilam,inu] <- rmse.best(trueB, Bhat)$min
    }
  }
  
  if(type=="acc.C"){
    wh.best <- (acc.C.mat) %>% {which(.==max(.),arr.ind=TRUE)} %>% {.[nrow(.),]}
  } else if(type=="acc.B"){
    wh.best <- (acc.B.mat) %>% {which(.==max(.),arr.ind=TRUE)} %>% {.[nrow(.),]}
  } else if(type=="rmse.B"){
    wh.best <- (rmse.B.mat) %>% {which(.==min(.),arr.ind=TRUE)} %>% {.[nrow(.),]}
  }
  
  out <- print.iSRRR(fit.tmp, wh.best[1], wh.best[2])
  best.cols <- rmse.best(trueB, out$B)$cols
  attr(out, "best.cols") <- best.cols
  
  out
}
