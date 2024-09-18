#' @importFrom MASS ginv
#' @importFrom rrr rrr

#' @export png.fit.best
png.fit.best <- function(fit, true){
  # true <- fit$trueB |> print.B()
  
  acc <- unlist( lapply( fit$B, function(est) png.accuracy(true, est) ) )
  fit.best <- iSRRR.choose.lambda(fit, which.max(acc))
  
  list(fit.best=fit.best, which=which.max(acc))
}








#' @export png.IC
png.IC <- function(fit, use.dfA=TRUE, sigma=1, logarithm=TRUE){
  
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
  
  out.bic <- do.call("rbind", out)
  
  out.bic2 <- out.bic[, colnames(out.bic) != "df"]
  list( BIC = out.bic,
        which.min = apply(out.bic2, 2, which.min), 
        which.min.new = apply(out.bic2, 2, function(x) { ifelse(which.min(x)!=1, which.min(x), min( which( x == min(x[x>min(x)]) ) ) ) } ),
        which.min2 = apply(out.bic2, 2, function(x) min( which( x == min(x[x>min(x)]) ) ) ) )
  
}






#' @export png.BICmat
png.BICmat <- function(fit, logarithm=TRUE, print.out=TRUE){
  
  tab <- do.call("cbind", lapply( 1:length(fit), function(idx) iSRRR.nu.IC(fit[[idx]], logarithm=logarithm)$BIC[,1] ))
  
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
















# Final ---------------------------------------------------------------------







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






#' @importFrom rrr rrr
#' @export BICk
BICk <- function(Y, X){
  
  n <- nrow(X);  p <- ncol(X);  q <- ncol(Y);
  
  rmax <- min(c(n,p,q))
  
  fit.svd <- svd(X)
  # reduce.svd <- function(fit.svd, rank){
  #   fit.svd$u <- fit.svd$u[,1:rank,drop=F]
  #   fit.svd$v <- fit.svd$v[,1:rank,drop=F]
  #   fit.svd$d <- fit.svd$d[1:rank]
  #   fit.svd
  # }
  
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
  
  pvec <- attr(B, "pvec") %>% { rep(1:length(.), .) }
  
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
  
  dfA <- sum( A !=0)
  dfB <- sum( B !=0)
  dfC <- sum( tcrossprod(B,A) !=0)
  
  
  if(use.dfA){
    df <- ifelse(dfB == 0, 9999, dfA + dfB)
  } else {
    df <- dfC
  }
  
  
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




#' @export iSRRR.nu.IC
iSRRR.nu.IC <- function(fit, use.dfA=TRUE, sigma=1, logarithm=TRUE){

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
  
  out.bic <- do.call("rbind", out)
  
  out.bic2 <- out.bic[, colnames(out.bic) != "df"]
  list( BIC = out.bic,
        which.min = apply(out.bic2, 2, which.min), 
        which.min.new = apply(out.bic2, 2, function(x) { ifelse(which.min(x)!=1, which.min(x), min( which( x == min(x[x>min(x)]) ) ) ) } ),
        which.min2 = apply(out.bic2, 2, function(x) min( which( x == min(x[x>min(x)]) ) ) ) )
  
}



#' @export iSRRR.BICmat
iSRRR.BICmat <- function(fit, logarithm=TRUE, print.out=TRUE){
  
  tab <- do.call("cbind", lapply( 1:length(fit), function(idx) iSRRR.nu.IC(fit[[idx]], logarithm=logarithm)$BIC[,1] ))
  
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





#' @export iSRRR.check.constraints
iSRRR.check.constraints <- function(fit){
  
  eps.mse <- fit$control$eps.mse
  
  # Convergence
  (wh.BIC <- fit %>% iSRRR.nu.IC() %>% .$which.min %>% .["BIC"])
  out_convergence <- sapply(fit$diff.mse, function(xx) sum(matrix(xx,ncol=4)[,4] < eps.mse))
  
  df1 <- NULL
  df2 <- matrix( fit$diff.mse[[wh.BIC]], ncol=4 )
  
  NROW <- max( nrow(df1), nrow(df2) )
  
  plot( df2[,4,drop=T], type="b", xlab="Iterations", ylab="MSE of C" )
  
  fit$B[[wh.BIC]] %>% print.B2()
  
  if(!is.null(fit$trueA)) cat( "Angle of A=", png.angle(fit$trueA, fit$A[[wh.BIC]])$max, "\n" )
  if(!is.null(fit$trueB)) cat( "Angle of B=", png.angle(fit$trueB, fit$B[[wh.BIC]])$max, "\n" )
  
  
  
  A <- fit$A
  nu <- fit$nu
  
  # Orthogonality of A
  out_orth <- lapply( A, function(x){
    cp <- crossprod(x)
    cp0 <- ifelse(abs(cp)<1e-12, 0, cp)
    isDiagonal( ifelse(abs(cp0)<1e-12, 0, cp0) )
  } )
  
  
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
