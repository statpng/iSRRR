#' @importFrom MASS ginv
#' @importFrom rrr rrr
#' @export png.fit.best
png.fit.best <- function(fit, true){
  # true <- fit$trueB |> print.B()
  
  acc <- unlist( lapply( fit$B, function(est) png.accuracy(true, est) ) )
  fit.best <- png.iSRRR.lambda(fit, which.max(acc))
  
  list(fit.best=fit.best, which=which.max(acc))
}



#' @export
check_constraints <- function(fit, cutoff.type="element-wise"){
  # cutoff.type = c("element-wise", "row-wise")[1]
  
  eps.mse <- fit$control$eps.mse
  
  # Convergence
  (wh.BIC <- fit %>% png.IC() %>% .$which.min %>% .["BIC"])
  
  out_convergence <- sapply(fit$diff.mse, function(xx) sum(matrix(xx,ncol=4)[,4] < eps.mse))
  
  # df1 <- matrix( fit$diff.mse[[wh.BIC-1]], ncol=4 )
  df1 <- NULL
  df2 <- matrix( fit$diff.mse[[wh.BIC]], ncol=4 )
  
  NROW <- max( nrow(df1), nrow(df2) )
  # par(mfrow=c(1,2))
  # plot_null <- range(c(df1[,4,drop=T], df2[,4,drop=T])) %>% {seq(min(.), max(.), length.out=NROW)}
  # plot( plot_null, col="white", xlab="Iterations", ylab="MSE of C" )
  # df1[,4,drop=T] %>% points(type="b")
  # df2[,4,drop=T] %>% points(type="b", col="red")
  plot( df2[,4,drop=T], type="b", xlab="Iterations", ylab="MSE of C" )
  
  fit$B[[wh.BIC]] %>% print.B2()
  
  if(!is.null(fit$trueA)) cat( "Angle of A=", png.angle(fit$trueA, fit$A[[wh.BIC]])$max, "\n" )
  if(!is.null(fit$trueB)) cat( "Angle of B=", png.angle(fit$trueB, fit$B[[wh.BIC]])$max, "\n" )
    
  # 180/pi*0.004
  
  
  A <- fit$A
  nu <- fit$nu
  
  # Orthogonality of A
  out_orth <- lapply( A, function(x){
    cp <- crossprod(x)
    cp0 <- ifelse(abs(cp)<1e-12, 0, cp)
    Matrix::isDiagonal( ifelse(abs(cp0)<1e-12, 0, cp0) )
  } )
  
  
  # Quartimax-simple of A
  check_quartimax <- function(A, nrep=1e4){
    r <- ncol(A)
    crit_quartimax <- sum( ( A )^4 )
    for( i in 2:nrep ){
      Q <- with( svd( matrix( runif(r*r, -1, 1), r, r ) ), tcrossprod(u, v) )
      crit_quartimax[i] <- sum( (A %*% Q)^4 )
    }
    mean( crit_quartimax[-1] > crit_quartimax[1] )
  }
  
  out_quartimax <- sapply( A, check_quartimax, nrep=1e3 )
  
  
  # Hard-thresholding of A
  out_threshold0 <- sapply(A, function(x){
    ifelse(abs(x)<1e-12, 0, x)
  })
  
  out_threshold <- sapply(A, function(x){
    x <- ifelse(abs(x)<1e-12, 0, x)
    x <- x[!apply(x, 1, function(xx) all(xx==0)),]
    mean( apply(x, 1, function(a) ( norm(a,"2") < nu )) )
    # mean( abs( x[x!=0] ) < nu )
  })
  
  
  OUT <- cbind( orthogonal=out_orth, quartimax=out_quartimax, convergence=out_convergence, threshold=out_threshold )
  attr(OUT, "wh.BIC") <- wh.BIC
  
  OUT
}


#' @export BICk
BICk <- function(X,Y,pvec){
  if(FALSE){
    
    ex2 <- function(seed=1){
      
      # rm(list=ls()); detach("package:iSRRR"); devtools::document()
      
      # library(iSRRR)
      
      
      set.seed(seed)
      data <- simdata(typeA="sparse0.2", typeB="random", n=20, d=5, q=20, pvec=rep(50,5), rvec=rep(1,5), es="2", simplify=TRUE, snr=1.0)
      
      X <- data$X
      Y <- data$Y
      pvec <- data$pvec
      nrank <- 6
      cutoff <- 0.2
      params = list(lambda.seq=1:10*0.1)
      control = list(best=TRUE)
      trueB = data$B
      
      example <- list(X=X, Y=Y, pvec=pvec, nrank=nrank, cutoff=cutoff,
                      params=params, control=control, trueB=trueB)
    }
    
    
    for( i in 11:50){
      data <- ex2(i)
      
      X <- data$X
      Y <- data$Y
      pvec <- data$pvec
      
      c(data$trueB |> print.B() |> png.duplicated() |> ncol(),
        BICk(X, Y, pvec)) |> print()
    }
    
    
    #
    
    get.rank.max <- function(X,Y){
      n <- nrow(X);  p <- ncol(X);  q <- ncol(Y);
      rank.max <- min(c(n,p,q))
      rank.max
    }
    
  }
  
  n <- nrow(X);  p <- ncol(X);  q <- ncol(Y);
  
  rmax <- min(c(n,p,q))
  
  fit.svd <- svd(X)
  reduce.svd <- function(fit.svd, rank){
    fit.svd$u <- fit.svd$u[,1:rank,drop=F]
    fit.svd$v <- fit.svd$v[,1:rank,drop=F]
    fit.svd$d <- fit.svd$d[1:rank]
    fit.svd
  }
  
  
  
  
  # fit.RankMax <- iSRRR(X=X, Y=Y, pvec=pvec, nrank=rmax, cutoff=0.0, params = list(lambda.seq=1e-10, omega=10000))
  # A <- fit.RankMax$A[[1]]
  # B <- fit.RankMax$B[[1]]
  
  # fit.rrr <- rrpack::rrr(Y,X)
  fit.rrr <- rrr(X,Y,rank=rmax,k=1e-10)
  
  BIC <- NULL
  for( r in 1:rmax ){
    # Cr <- tcrossprod( B[,1:r], A[,1:r] )
    # Cr <- with( reduce.svd(fit.svd, r), as.matrix(v) %*% diag(1/d,r,r) %*% t(u) ) %*% Y
    Cr <- t(fit.rrr$A[,1:r,drop=F] %*% fit.rrr$B[1:r,,drop=F])
    
    if(TRUE){
      const <- sqrt(r*p*log(p))
      hn <- log(n)
      penalty <- const*hn
      BIC[r] <- sum((Y - X %*% Cr)^2) + penalty
    }
    
    if(FALSE){
      E <- (Y - X %*% Cr)
      loss <- n*log(det(crossprod(E)))
      penalty <- r*( log2(n) + 23*log2(q) )
      BIC[r] <- loss + penalty
    }
    
  }
  
  BIC |> which.min()
  
}




bic <- function(SSE, npqdr, A, B, sigma=1, use.dfA=TRUE, logarithm=TRUE){
  
  if(FALSE){
    X <- data$X
    Y <- data$Y
    A <- fit$A[[2]]
    B <- fit$B[[2]]
    use.dfA <- TRUE
    
    fit.svd <- svd( Y %*% A )
    
    C <- with( svd(X), v %*% diag(1/d) %*% u %*% Y  )
    
    sigma <- sqrt( mean( ( Y - X %*% C )^2 ) )
    
  }
  
  
  pvec <- attr(B, "pvec")
  pvec <- rep(1:length(pvec), pvec)
  
  # n <- nrow(Y); q <- ncol(Y)
  # d <- max(pvec)
  # p <- ncol(X); r <- ncol(A)
  
  n <- npqdr["n"] |> as.numeric()
  p <- npqdr["p"] |> as.numeric()
  q <- npqdr["q"] |> as.numeric()
  d <- npqdr["q"] |> as.numeric()
  r <- npqdr["r"] |> as.numeric()
  
  
  npq <- n*p*q
  nq <- n*q
  pq <- p*q
  
  # xrank <- sum(svd(X)$d > 0)
  # sigma <- sqrt( mean( ( Y - X %*% tcrossprod(B, A) )^2 ) )
  
  # SSE <- sum( ( Y - X %*% tcrossprod(B, A) )^2 ) / nq / sigma^2
  
  # # lasso
  # df <- sum(A!=0)+sum(B!=0)
  
  # group lasso
  Bnew <- apply(B,2,function(bk) tapply(bk, pvec, function(bik) norm(bik, "2")) )
  
  dfA <- sum( A !=0)
  dfB <- sum( B !=0)
  dfC <- sum( tcrossprod(B,A) !=0)
  
  # {
  #   # X %*% MASS::ginv( crossprod(X) ) %*% t(X)
  #   YPY <- t(Y) %*% with( png.svd(X), tcrossprod(u, u) ) %*% Y
  #   ypyv <- eigen(YPY)$values
  #
  #   sl <- 0
  #   for( kkk in 1:r ){
  #     lk <- ypyv[kkk]
  #     for( lll in (kkk+1):q ){
  #       ll <- ypyv[lll]
  #       sl <- 2 * ll / (lk-ll)
  #     }
  #   }
  #
  #   df <- (p+q-r)*r + 2 * sl
  # }
  
  
  if(use.dfA){
    df <- ifelse(dfB == 0, 9999, dfA + dfB)
  } else {
    df <- dfC
  }
  
  # df <- dfA + dfB * xrank/d - r^2
  
  
  #cat("SSE=", SSE * nq, "\n")
  #cat("df=", df, "\n")
  #cat("LHS=", (SSE)/sigma^2, "\n")
  #cat("RHS=", log(nq) / nq * df, "\n")
  #cat("BIC=", (SSE)/sigma^2 + log(nq) / nq * df, "\n")
  #cat("---------------------------", "\n")
  
  if( logarithm ){
    list(
      BIC = log(SSE) + log(nq) / nq * df,
      BIC2 = log(SSE/nq) + log(nq) / nq * df,
      GIC = log(SSE) + log(log(nq)) * log(pq)/nq * df,
      AIC = log(SSE) + 2/nq * df,
      GCV = nq*SSE/(nq - df)^2,
      # GCV.Yuan = SSE/(nq - df),
      SSE=SSE, logSSE=log(SSE), df=df
    )
  } else {
    list(
      BIC = (SSE)/sigma^2 + log(nq) / nq * df,
      GIC = (SSE) + log(log(nq)) * log(pq)/nq * df,
      AIC = (SSE) + 2/nq * df,
      GCV = nq*SSE/(nq - df)^2,
      # GCV.Yuan = SSE/(nq - df),
      SSE=SSE, df=df
    )
  }
  
  
}




#' @export png.IC
#' 
#' 
png.IC <- function(fit, use.dfA=TRUE, sigma=1, logarithm=TRUE){
  if(FALSE){
    library(iSRRR)
    set.seed(1)
    n=50; d=4; q=50; p=50; r=6; es="1"; snr=0.5
    
    params <- list(nlambda=5, lambda.factor=1e-4)
    control <- list(best=TRUE, maxit.mse=20, eps.mse=1e-4, maxit.B=1e4, eps.B=1e-4)
    
    data <- simdata2(typeA="sparse0.4", typeB="row", n=n, d=d, d0=3, q=q, pvec=rep(p,d), rvec=rep(1,r), es=es, simplify=TRUE, snr=snr, rho_X = 0.5)
    
    fit <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec, nrank = 6,
                  cutoff = 0.4, 
                  params=list(nlambda=20, lambda.factor=1e-4),
                  control=list(best=TRUE, maxit.mse=20, eps.mse=1e-4, maxit.B = 1e4, eps.B = 1e-4), trueB = data$B, gglasso=TRUE )
    
    fit |> png.IC()
    
    format( object.size(fit), "Mb" )
    # "0.5 Mb"
    format( object.size( fit$B[[5]]%*%t(fit$A[[5]]) ), "Mb" )
    format( object.size( fit$B[[5]]%*%t(fit$A[[5]]) ), "Mb" )
    
    {
      use.dfA=TRUE
      sigma=1
      logarithm=FALSE
    }
    
  }
  
  params <- fit$params
  control <- fit$control
  SSE <- fit$SSE
  npqdr <- attributes(SSE)$npqdr
  A <- fit$A;  B <- fit$B
  
  out <- NULL
  for( i in 1:params$nlambda){
    # A <- fit$A[[i]];  B <- fit$B[[i]]
    
    out[[i]] <- do.call("c", bic(SSE[[i]],npqdr, A[[i]],B[[i]], use.dfA=use.dfA, sigma = sigma, logarithm=logarithm))
    
  }
  
  # apply(out,2,which.min)
  out.bic <- do.call("rbind", out)
  
  out.bic2 <- out.bic[, colnames(out.bic) != "df"]
  list( BIC = out.bic,
        which.min = apply(out.bic2, 2, which.min), 
        which.min.new = apply(out.bic2, 2, function(x) { ifelse(which.min(x)!=1, which.min(x), min( which( x == min(x[x>min(x)]) ) ) ) } ),
        which.min2 = apply(out.bic2, 2, function(x) min( which( x == min(x[x>min(x)]) ) ) ) )
  
}



#' @export png.angle
png.angle <- function(true, est){
  # The largest principal angle
  qt <- qr.Q(qr(true))
  qe <- qr.Q(qr(est))
  fit.svd <- svd( crossprod(qe, qt) )
  theta <- acos(fit.svd$d |> round(12))
  
  # theta[1] * 180 / pi
  list( max = theta[1] * 180 / pi, Grassmanian = norm( theta, "2" ) * 180 / pi )
}


#' @export png.measure
png.measure <- function(true, est){
  
  true.01 <- true |> (function(.) ifelse(abs(.)>0, 1, 0))()
  est.01 <- est |> (function(.) ifelse(abs(.)>0, 1, 0))()
  
  tp <- sum((true.01==est.01)[true.01==1])
  fp <- sum((true.01!=est.01)[true.01==0])
  tn <- sum((true.01==est.01)[true.01==0])
  fn <- sum((true.01!=est.01)[true.01==1])
  
  tb <- table(true.01, est.01)
  tpr <- tp / sum(true.01==1)
  fdr <- ifelse(sum(est.01==1)==0, 0, fp / sum(est.01==1))
  f1 <- (2*tp)/(2*tp + fp + fn)
  mcc <- ifelse( any(c((tp+fp),(tp+fn),(tn+fp),(tn+fn)) == 0), 0, (tp*tn - fp*fn)/(sqrt(tp+fp)*sqrt(tp+fn)*sqrt(tn+fp)*sqrt(tn+fn)) )
  acc <- (tp+tn)/(tp+fp+tn+fn)
  
  list(tb=tb, tpr=tpr, fdr=fdr, f1=f1, mcc=mcc, acc=acc)
}





#' @export print.iSRRR
print.iSRRR <- function(fit, Lambda, Nu){
  Ahat <- fit[[ Nu ]]$A[[ Lambda ]]
  Bhat <- fit[[ Nu ]]$B[[ Lambda ]]
  
  list(A=Ahat, B=Bhat)
}


#' @export png.BICmat
png.BICmat <- function(fit, logarithm=TRUE, print.out=TRUE){
  
  tab <- do.call("cbind", lapply( 1:length(fit), function(idx) png.IC(fit[[idx]], logarithm=logarithm)$BIC[,1] ))
  
  attr(tab, "wh.min") <- which( tab == min(tab), arr.ind=TRUE )
  
  wh.lambda <- attr(tab, "wh.min")[1,1]
  wh.nu <- attr(tab, "wh.min")[1,2]
  
  if(print.out){
    print( print.B2( fit[[wh.nu]]$B[[wh.lambda]] ) )
    print( round( fit[[wh.nu]]$A[[wh.lambda]], 3 ) )
    
    if( !is.null(dim(fit[[wh.nu]]$diff.mse[[wh.lambda]])) ){
      plot( fit[[wh.nu]]$diff.mse[[wh.lambda]][,4], type="b", xlab="Iterations", ylab="MSE of C" )
    }
    
  }
  
  trueA <- fit[[wh.nu]]$trueA
  trueB <- fit[[wh.nu]]$trueB
  
  if(!is.null(trueA)) cat( "angle of A =", png.angle( trueA, fit[[wh.nu]]$A[[wh.lambda]] )$max, "\n" )
  if(!is.null(trueB)) cat( "angle of B =", png.angle( trueB, fit[[wh.nu]]$B[[wh.lambda]] )$max, "\n" )  
  
  
  tab
  
}