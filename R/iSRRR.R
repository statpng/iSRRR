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
#' @param control A list with \code{best}=FALSE, \code{early.stop}=TRUE, \code{rot.method}="quartimax" \code{maxit.B}=3e8, \code{eps.B}=1e-8, \code{maxit.mse}=50, \code{eps.mse}=1e-6, \code{X.scale}=c("group", "each"), \code{Y.scale}=FALSE, \code{verbose}=FALSE, and \code{threads}=1.
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
#'
#'
#' @importFrom gglasso gglasso
#' @import Rcpp
#' @importFrom MASS ginv
#' @useDynLib iSRRR
#'
#' @export arcov
#' @export iSRRR
iSRRR <- function(X, Y, pvec, nrank, cutoff, params=NULL, control=NULL, trueB=NULL, use.gglasso=TRUE){


  if(FALSE){

    library(iSRRR)
    set.seed(1)
    n=50; d=4; q=50; p=50; r=6; es="1"; snr=0.5

    params <- list(nlambda=5, lambda.factor=1e-4, group.size=FALSE)
    control <- list(best=TRUE, maxit.mse=20, eps.mse=1e-4, maxit.B=1e4, eps.B=1e-4)

    data <- simdata.sgl(typeA="sparse0.0", typeB="individual", n=n, d=d, q=q, pvec=c(10,10,100,200), rvec=rep(1,r), es=es, simplify=TRUE, snr=snr, rho_X = 0.5, param.sgl=c(1.0, 1.0, 0.5, 0.2))
    
    
    data <- simdata2(typeA="sparse0.0", typeB="row", n=n, d=d, d0=3, q=q, pvec=rep(p,d), rvec=rep(1,r), es=es, simplify=TRUE, snr=snr, rho_X = 0.5)
    
    
    
    fit1 <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec, nrank = 6,
                  cutoff = 0.4,
                  params=list(nlambda=20, lambda.factor=1e-4, group.size=TRUE),
                  control=list(best=TRUE, maxit.mse=20, eps.mse=1e-4, maxit.B = 1e4, eps.B = 1e-4), trueB = data$B, use.gglasso=TRUE )

    
    
    fit2 <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec, nrank = 6,
                  cutoff = 0.4,
                  params=list(nlambda=20, lambda.factor=1e-5, group.size=FALSE),
                  control=list(best=TRUE, maxit.mse=20, eps.mse=1e-4, maxit.B = 1e4, eps.B = 1e-4), trueB = data$B, use.gglasso=TRUE )
    
    
    fit1 %>% png.IC()
    fit2 %>% png.IC()
    
    fit1$B[[7]] %>% print.B2()
    fit2$B[[7]] %>% print.B2()
    
    #
    
    fit1 %>% { .$B[[which.min(png.IC(.)[,"BIC"])]] } |> print.B()
    fit2 %>% { .$B[[which.min(png.IC(.)[,"BIC"])]] } |> print.B()
    
    #
    
    

    X=data$X; Y=data$Y; pvec=data$pvec;
    nrank=r;
    cutoff=0.2;
    params=params;
    control=control; trueB=data$B
    use.gglasso=TRUE

  }
  


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
    early.stop <- control$early.stop
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



  n <- nrow(Y);  p <- length(pvec);
  d <- max(pvec);  q <- ncol(Y)



  fit.init <- get.init(X=X, Y=Y, pvec=pvec, nrank=nrank, lambda.seq=lambda.seq,
                       lambda.factor=lambda.factor, nlambda=nlambda, log.scale=log.scale, group.size=group.size)


  A0 <- fit.init$A0
  B0 <- fit.init$B0

  params$lambda.seq <- fit.init$lambda.seq
  params$nlambda <- fit.init$nlambda


  start <- proc.time()
  out <- fit.iSRRR(
    X = X,
    Y = Y,
    A0 = A0,
    B0 = B0,
    pvec = pvec,
    nrank = nrank,
    cutoff = cutoff,
    params = params,
    control = control,
    trueB = trueB, use.gglasso=use.gglasso
  )
  end <- proc.time()


  SSE.list <- lapply( out, function(x) mean( (Y-X%*%x$C)^2 ) )
  attr(SSE.list, "npqdr") <- c(n=nrow(X), p=ncol(X), q=ncol(Y), d=max(pvec), r=nrank)


  list(
    # X = X, Y = Y,
    pvec=pvec, nrank=nrank,
    A = lapply( out, function(x) x$A ),
    B = lapply( out, function(x) x$B ),
    # C = lapply( out, function(x) x$C ),
    SSE = SSE.list,
    diff.mse = lapply( out, function(x) x$diff ),
    # diff.B = lapply( out, function(x) x$diff.B ),
    simplexity = lapply( out, function(x) x$simplexity ),
    # iter.B = lapply( out, function(x) x$it.B ),
    iter.mse = lapply( out, function(x) x$it.mse ),
    cutoff = cutoff,
    params = params,
    control = control,
    trueB = trueB, time = end-start )

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
    nlambda=5,
    lambda.factor=1e-4, # 0.0001
    log.scale=TRUE,
    group.size=FALSE
  )

  control0 = list(
    best = FALSE,
    early.stop = TRUE,
    rot.method = "quartimax", # c("varimax", "quartimax", "oblimin", "oblimax", "infomax", "CF", "McCammon", "simplimax", "tandem", "entropy"),
    maxit.B = 3e8,
    eps.B = 1e-8,
    maxit.mse = 50,
    eps.mse = 1e-6,
    X.scale = FALSE, # c("group", "each"),
    Y.scale = FALSE, # TRUE, FALSE
    verbose=FALSE,
    threads=1
  )

  list(params=params0, control=control0)
}







fit.iSRRR <- function(X,Y,A0,B0,pvec,nrank,cutoff,params,control,trueB=NULL, use.gglasso=TRUE){


  n <- nrow(Y);  p <- length(pvec);
  d <- max(pvec);  q <- ncol(Y)


  {

    lambda.seq <- params$lambda.seq
    nlambda <- params$nlambda
    lambda.factor <- params$lambda.factor
    log.scale <- params$log.scale
    group.size <- params$group.size


    best <- control$best
    early.stop <- control$early.stop
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




  # gglasso -----------------------------------------------------------------
  if(use.gglasso){

    out <- NULL
    for( ll in 1:nlambda ){

      lambda <- lambda.seq[ll]
      # if(group.size){
      #   lambda <- lambda * sqrt(pik)[,1]
      # }

      it.mse <- 1

      diff <- matrix(NA, nrow=maxit.mse, ncol=4)
      colnames(diff) <- c("Y", "A", "B", "C")
      diff[1,] <- eps.mse*2

      simplexity <- matrix(NA, nrow=maxit.mse, ncol=2)
      colnames(simplexity) <- c("quartimax", "varimax")

      minimum <- 9e+5
      min.count <- 0
      while((it.mse < maxit.mse)&( diff[it.mse,"C"] > eps.mse)){

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


        {
          Q <- Update.Rotation(A, rot.method = rot.method)$rotation

          A <- A %*% Q
          B <- B %*% Q
        }


        {
          A <- with( png.svd( crossprod(Y, X) %*% B ), tcrossprod( u, v ) )
          A <- ifelse( abs(A) > cutoff, A, 0 )
          A <- with(png.svd(A), tcrossprod(u,v))
        }

        {
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
        }



        {
          C <- tcrossprod(B, A)
          Cold <- tcrossprod(Bold, Aold)

          simplexity[it.mse, "quartimax"] <- crit.quartimax(A)
          simplexity[it.mse, "varimax"] <- crit.varimax(A)

          diff[it.mse, "Y"] <- mean( (X %*% C - X %*% Cold)^2 )
          diff[it.mse, "A"] <- mean( (A - Aold)^2 )
          diff[it.mse, "B"] <- mean( (B - Bold)^2 )
          diff[it.mse, "C"] <- mean( (C - Cold)^2 )
        }


      }

      if(it.mse == maxit.mse){
        warning(paste("Algorithm didn't converge in ", it.mse, " iterations at lambda[", ll,"] !", sep = ""))
      }

      simplexity <- simplexity[1:it.mse,]
      diff <- diff[1:it.mse,][-c(1:2),]


      if(best & !is.null(trueB) ){

        if( any( dim(B) != dim(trueB) ) ) stop("Dimension of B must be the same.")

        best.cols <- rmse.best(trueB, B, best = TRUE)$cols
        A <- A[,best.cols]
        B <- B[,best.cols]

        fit.sign <- rmse.sign.best(trueB, B)

        B <- fit.sign$est
        A <- A %*% diag(as.numeric(fit.sign$sign))

      }


      attr(B, "pvec") <- as.numeric(table(pvec))



      out[[ll]] <- list(A=A, B=B, C=C,
                        diff=diff,
                        simplexity=simplexity,
                        it.mse=it.mse )
    }



    # iSRRR -------------------------------------------------------------------------

  } else {


    out <- NULL
    for( ll in 1:nlambda ){

      lambda <- lambda.seq[ll]
      if(group.size){
        lambda <- lambda * sqrt(pik)
      }
      lambda_mat <- lambda
      # lambda_mat <- matrix(lambda,p,nrank)

      it.mse <- 1

      diff <- matrix(NA, nrow=maxit.mse, ncol=4)
      colnames(diff) <- c("Y", "A", "B", "C")
      diff[1,] <- eps.mse*2

      simplexity <- matrix(NA, nrow=maxit.mse, ncol=2)
      colnames(simplexity) <- c("quartimax", "varimax")

      minimum <- 9e+5
      min.count <- 0
      while((it.mse < maxit.mse)&( diff[it.mse,"C"] > eps.mse)){

        it.mse <- it.mse + 1

        if( it.mse == 2 ){
          A <- A0
          B <- B0
        }



        # Old A, B ----------------------------------------------------------------
        {
          Aold <- A
          Bold <- B
        }


        # Update Q ----------------------------------------------------------------
        if(verbose) start <- proc.time()

        {
          Q <- Update.Rotation(A, rot.method = rot.method)$rotation

          A <- A %*% Q
          B <- B %*% Q
        }

        if(verbose){
          print("Update Q")
          print( proc.time() - start )
        }


        # Update A ----------------------------------------------------------------
        if(verbose) start <- proc.time()

        {
          A <- with( png.svd( crossprod(Y, X) %*% B ), tcrossprod( u, v ) )
          A <- ifelse( abs(A) > cutoff, A, 0 )
          A <- with(png.svd(A), tcrossprod(u,v))
        }

        if(verbose){
          print("Update A")
          print( proc.time() - start )
        }



        # Update B ----------------------------------------------------------------
        if(verbose) start <- proc.time()

        {
          fit.B <- UpdateB_BMD(X = X, Y = Y, A = A, B = B,
                               pvec = pvec, nrank = nrank, d = d,
                               lam = lambda_mat, maxit = maxit.B, eps = eps.B, threads = threads)
          # fit.B <- UpdateB_SGM(X = X, Y = Y, pvec = pvec,
          #                      lam = lam_mat, d = d, A = A, B = B,
          #                      nrank = nrank, eps = eps.B, maxit = maxit.B)

          diff.B <- mapply(function(x,y) x[1:y], fit.B$diff, fit.B$iter)
          iter.B <- fit.B$iter
          B <- fit.B$B
        }

        if(verbose){
          print("Update B")
          print( proc.time() - start )
        }



        # Convergence ---------------------------------------------------------------------
        {
          C <- tcrossprod(B, A)
          Cold <- tcrossprod(Bold, Aold)

          simplexity[it.mse, "quartimax"] <- crit.quartimax(A)
          simplexity[it.mse, "varimax"] <- crit.varimax(A)

          diff[it.mse, "Y"] <- mean( (X %*% C - X %*% Cold)^2 )
          diff[it.mse, "A"] <- mean( (A - Aold)^2 )
          diff[it.mse, "B"] <- mean( (B - Bold)^2 )
          diff[it.mse, "C"] <- mean( (C - Cold)^2 )
        }


        if( early.stop ){

          minimum <- min(minimum, diff[it.mse, "C"])
          if( diff[it.mse, "C"] > minimum ){
            min.count <- min.count+1
            if(min.count == 100){
              if(verbose){
                print("early stopping!")
                # cat("iteration=",it.mse,"\n")
                # cat("minimum=",minimum,"\n")
                # cat("current MSE=",diff[it.mse, "C"],"\n")
              }
              break;
            }
          } else {
            min.count <- 0
          }

        }


        # .Machine$double.eps

      }

      if(it.mse == maxit.mse){
        warning(paste("Algorithm didn't converge in ", it.mse, " iterations at lambda[", ll,"] !", sep = ""))
      }

      simplexity <- simplexity[1:it.mse,]
      diff <- diff[1:it.mse,][-c(1:2),]


      if(best & !is.null(trueB) ){

        if( any( dim(B) != dim(trueB) ) ) stop("Dimension of B must be the same.")

        best.cols <- rmse.best(trueB, B, best = TRUE)$cols
        A <- A[,best.cols]
        B <- B[,best.cols]

        fit.sign <- rmse.sign.best(trueB, B)

        B <- fit.sign$est
        A <- A %*% diag(as.numeric(fit.sign$sign))

      }


      attr(B, "pvec") <- as.numeric(table(pvec))



      out[[ll]] <- list(A=A, B=B, C=C,
                        diff=diff,
                        simplexity=simplexity,
                        it.mse=it.mse,
                        it.B=iter.B,
                        diff.B=diff.B )
    }


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
  fit$A <- fit$A[idx]
  fit$B <- fit$B[idx]
  fit$SSE <- fit$SSE[idx]
  fit$diff.mse <- fit$diff.mse[idx]
  fit$simplexity <- fit$simplexity[idx]
  fit$iter.mse <- fit$iter.mse[idx]

  fit
}



