#' @importFrom gtools mixedsort
#' @export print.B
print.B <- function(B, pvec=NULL, type=c("norm", "head", "prop"), top=1, rowwise=FALSE){
  
  if(FALSE){
    B <- fit$B[[4]]
  }
  
  if( is.null(pvec) ){
    pvec2 <- attr(B, "pvec")
    pvec <- rep( 1:length(pvec2), pvec2 )
  }
  d <- max(pvec)
  
  
  out <- NULL
  for( i in 1:d){
    
    if(!rowwise){
      out1 <- NULL
      for( k in 1:ncol(B) ){
        
        if( type == "head" ){
          out1 <- cbind( out1, B[which( pvec == i )[1:top], k] )
        } else if(type == "norm"){
          out1 <- cbind( out1, norm( B[which( pvec == i ), k], "2") )
        } else if(type == "prop"){
          out1 <- cbind( out1, mean( abs( B[which( pvec == i ), k] ) > 1e-10 ) )
        }
        
      }
    } else {
      if( type == "head" ){
        out1 <- B[which( pvec == i )[1:top], ]
      } else if(type == "norm"){
        out1 <- norm( B[which( pvec == i ), ], "2")
      } else if(type == "prop"){
        out1 <- mean( apply(B[which( pvec == i ), ], 1, function(x) abs( x ) > 1e-10 ) )
      }
    }
    
    
    out <- rbind( out, out1 )
    
  }
  
  round( out, 4)
}


#' @export print.B2
print.B2 <- function(B, top=1, group.size=FALSE){
  
  if(FALSE){
    B <- fit$B[[4]]
  }
  
  pvec <- attr(B, "pvec")
  d <- length(pvec)
  pvec2 <- rep( 1:d, pvec )
  
  out <- NULL
  for( i in 1:d){
    out1 <- NULL
    for( k in 1:ncol(B) ){
      if(group.size){
        out1 <- cbind( out1, norm( B[which( pvec2 == i ),k], "2" ) / sqrt(sum(pvec2 == i)) )
      } else {
        out1 <- cbind( out1, norm( B[which( pvec2 == i ),k], "2" ))
      }
      
    }
    out <- rbind( out, out1 )
  }
  
  round( out, 4)
}



#' @export png.duplicated
png.duplicated <- function(A){
  a <- A
  a <- a[,!apply( a, 2, function(x) all(x==0) )]
  a.binary <- a |> (function(.) ifelse(.!=0, 1, 0))() |> as.data.frame()
  out <- a[,which(!duplicated(as.list(a.binary)))] |> as.matrix()
  colnames(out) <- NULL
  out
}



#' @export png.df.rename
png.df.rename <- function(df, names){
  if(FALSE){
    names <- c("n", "p", "q", "r", "d", "snr")
  }
  
  wh.names <- which(colnames(df) %in% names)
  for(j in 1:length(wh.names) ){
    df[wh.names[j]] <- paste0(names[j], "=", unlist(df[wh.names[j]]))
    
    xj <- as.factor(unlist(df[wh.names[j]])) 
    xj <- factor(xj, 
                 levels=mixedsort( levels(unique(xj)) ),
                 labels=mixedsort( levels(unique(xj)) ) )
    
    df[wh.names[j]] <- xj
  }
  
  df
}


#' @export png.list2str
png.list2str <- function(params){
  p1 <- lapply(params, function(x) ifelse(length(x)>=10, paste0(x[c(1,length(x))], collapse = "-"), paste0(x, collapse = "+") ) )
  
  p2 <- NULL
  for( i in 1:length(p1)){
    p2[[i]] <- paste0(names(p1[i]), "=", p1[[i]])
  }
  
  p3 <- paste0( p2, collapse="__" )
  
  p3
}


#' @export png.assign.list
png.assign.list <- function(params){
  if(FALSE){
    params = list(a=1, b=2)
  }
  
  list2env(params, envir = .GlobalEnv)
  
  NAMES <- names(params)
  LIST <- params
  
  invisible(
    lapply(seq_along(params), 
         function(x) {
           assign(NAMES[x], LIST[[x]], envir=.GlobalEnv)
         }
    )
  )
  
}


#' @export png.lambda
png.lambda <- function(X, Y, pvec, nlambda = 10, lambda.factor = 0.05){
  
  if(FALSE){
    Yk <- Y %*% A[,1]
    nlambda <- 10
    lambda.factor = 0.05
    
    # gglasso::gglasso(X, Yk, group=pvec, nlambda = 11, lambda.factor = 0.05)
    # ll <- lambda_sequence(X, Yk, pvec, nlambda = 11, lambda.factor = 0.05)
    # lamfix(ll)
  }

  
  n <- nrow(Y);  p <- length(pvec);
  d <- max(pvec);  q <- ncol(Y)
  
  C.init <- crossprod(with( png.svd(X), u %*% diag(1/d) %*% t(v) ), Y)
  
  r <- 1
  
  coefSVD <- png.svd(C.init)
  coefSVD$u <- coefSVD$u[, 1:r, drop = FALSE]
  coefSVD$d <- coefSVD$d[1:r]
  coefSVD$v <- coefSVD$v[, 1:r, drop = FALSE]
  
  A0 <- coefSVD$v
  B0 <- with(coefSVD, u %*% diag(d, r, r))
  
  lambda.seq <- lambda_sequence(X, Y, A0, pvec=pvec, nlambda=nlambda, lambda.factor = lambda.factor)
  
  
  list( lambda.seq=lambda.seq, nlambda=nlambda, lambda.factor=lambda.factor)
  
}




#' @export lambda_sequence
lambda_sequence <- function(X, Y, A, pvec, nlambda = 10, lambda.factor = 0.05, log.scale=FALSE, group.size=FALSE){

  if(FALSE){
    {
      
      data1 <- simdata(typeA="sparse0.2", typeB="partial", n=10, d=3, q=10, pvec=rep(5,3), rvec=rep(2,3), es="3", simplify=TRUE, snr=1.0)
      data2 <- simdata(typeA="sparse0.2", typeB="partial", n=500, d=3, q=10, pvec=rep(5,3), rvec=rep(2,3), es="3", simplify=TRUE, snr=1.0)
      
      data <- data1
      l1 <- get.init(data$X, data$Y, data$pvec, ncol(data$B), NULL, 1e-3, 0, 5)$lambda.seq
      fit1 <- iSRRR(data$X, data$Y, data$pvec, ncol(data$B), 0.2, params=list(lambda.seq=l1), control=NULL, trueB=data$B)
      
      data <- data2
      l2 <- get.init(data$X, data$Y, data$pvec, ncol(data$B), NULL, 1e-3, 0, 5)$lambda.seq
      fit2 <- iSRRR(data$X, data$Y, data$pvec, ncol(data$B), 0.2, params=list(lambda.seq=l2), control=NULL, trueB=data$B)
      
      fit1$B[[2]] |> print.B()
      fit2$B[[2]] |> print.B()
      
      #######################
      
      res <- res1
      
      X <- res$fit$X
      Y <- res$fit$Y
      pvec <- res$fit$pvec
      nrank <- res$fit$nrank
      lambda.seq <- NULL
      lambda.factor <- res$fit$params$lambda.factor
      gamma <- res$fit$params$gamma
      nlambda <- res$fit$params$nlambda
      
      
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
      
      A <- coefSVD$v
      # A0 <- coefSVD$v
      # B0 <- with(coefSVD, u %*% diag(d, nrank, nrank))
      
    }
    
    
    Yk <- Y %*% A[,1]
    nlambda <- 10
    lambda.factor = 0.05

    # gglasso::gglasso(X, Yk, group=pvec, nlambda = 11, lambda.factor = 0.05)
    # ll <- lambda_sequence(X, Yk, pvec, nlambda = 11, lambda.factor = 0.05)
    # lamfix(ll)
  }

  
  
  
  
  
  big <- 9.9e30
  bs <- ifelse( group.size, as.integer(as.numeric(table(pvec))), 1 )
  d <- max(pvec)
  nobs <- nrow(X)
  
  out <- NULL
  for( k in 1:ncol(A) ){

    Yk <- Y %*% A[,k]

    ga <- NULL
    for( i in 1:d ){
      Syx <- crossprod(Yk, X[, which(pvec==i)]) / nobs
      ga[i] <- norm( Syx, "F" )
    }

    out <- cbind(out, ga / sqrt(bs))
  }



  if( log.scale ){
    
    flmin <- max(1e-6, lambda.factor)
    alf <- flmin^(1.0/(nlambda-1))
    
    al <- NULL
    for( i in 1:nlambda ){
      if( i == 1 ){
        al[i] <- big
      } else if ( i == 2 ) {
        al[i] = max(0, max(out) ) * alf
      } else {
        al[i] <- al[i-1] * alf
      }
    }
    
    
    
    lamfix <- function(lam) {
      llam <- log(lam)
      lam[1] <- exp(2 * llam[2] - llam[3])
      lam
    }
    
    
    
    lamfix( al )
    
  } else {
    
    al <- NULL
    
    al[1] <- big
    al[2:nlambda] <- seq(max(out), max(out)*lambda.factor, length.out=nlambda-1)
    
    al
    
  }
  

}





mat.order.nonzero <- function(B, pvec){

  # nonzero.mat <- ifelse( B!=0, 1, 0 )
  nonzero.mat <- ifelse( B!=0, 1, 0 )

  out <- apply(nonzero.mat, 2, sum)
  for( i in 1:max(pvec) ){
    out <- out + apply(nonzero.mat[which(pvec == i),], 2, sum)*1/i
  }
  B[,order(out, decreasing=TRUE)]
}




# initialize <- function(X, Y, pvec, nrank, eta=0.1, method="JASA"){
# 
#   p <- length(pvec)
# 
# 
# 

# 
#   if( method == "JASA2" ){
#     C0 <- NULL
#     for(i in 1:max(pvec)){
#       Xi <- X[,which( pvec == i )]
#       pi <- sum( pvec == i )
# 
#       Sxx <- crossprod(Xi)
#       Syx <- crossprod(Y, Xi)
#       Sxx_inv <- tryCatch(ginv(Sxx), error = function(e) solve(Sxx + 0.1 * diag(pi)))
#       # Sxx_inv <- solve(Sxx + eta * diag(pi))
# 
#       C0 <- rbind(C0, tcrossprod(Sxx_inv, Syx))
#     }
# 
#     coefSVD <- svd(C0, nrank, nrank)
#     coefSVD$u <- coefSVD$u[, 1:nrank, drop = FALSE]
#     coefSVD$d <- diag(coefSVD$d[1:nrank])
#     coefSVD$v <- coefSVD$v[, 1:nrank, drop = FALSE]
# 
#     B0 <- coefSVD$u #%*% coefSVD$d
#     A0 <- coefSVD$v
# 
# 
# 
#   } else if( method == "JASA" ){
# 
#     eta <- sqrt( log(p) / n )
#     Sxx <- crossprod(X)
#     Syx <- crossprod(Y, X)
#     Sxx_inv <- solve(Sxx + eta * diag(p))
# 
#     C0 <- tcrossprod(Sxx_inv, Syx)
# 
#     coefSVD <- svd(C0, nrank, nrank)
#     coefSVD$u <- coefSVD$u[, 1:nrank, drop = FALSE]
#     # coefSVD$d <- diag(coefSVD$d[1:nrank])
#     coefSVD$d <- diag(coefSVD$d[1:nrank]^(-1))
#     coefSVD$v <- coefSVD$v[, 1:nrank, drop = FALSE]
# 
#     # B0 <- coefSVD$u
#     A0 <- coefSVD$v %*% coefSVD$d
# 
# 
# 
# 
#   } else if( method == "ridge" ){
# 
#     eta <- sqrt( log(p) / n )
#     Sxx <- crossprod(X)
#     Syx <- crossprod(Y, X)
#     Sxx_inv <- solve(Sxx + eta * diag(p))
# 
#     C0 <- tcrossprod(Sxx_inv, Syx)
# 
#     coefSVD <- svd(C0, nrank, nrank)
#     coefSVD$u <- coefSVD$u[, 1:nrank, drop = FALSE]
#     # coefSVD$d <- diag(coefSVD$d[1:nrank])
#     coefSVD$d <- diag(coefSVD$d[1:nrank]^(-1))
#     coefSVD$v <- coefSVD$v[, 1:nrank, drop = FALSE]
# 
#     # B0 <- coefSVD$u
#     A0 <- coefSVD$v %*% coefSVD$d
# 
# 
#   } else if( method == "MGLasso" ){
# 
#     for(k in 1:nrank){
#       gglassogglasso( X, Y%*%A[,k], pvec, lambda = 0.000001 )
#     }
# 
#     Sxx <- crossprod(X)
#     Syx <- crossprod(Y, X)
#     Sxx_inv <- solve(Sxx + eta * diag(p))
# 
#     C0 <- tcrossprod(Sxx_inv, Syx)
# 
#     coefSVD <- svd(C0, nrank, nrank)
#     coefSVD$u <- coefSVD$u[, 1:nrank, drop = FALSE]
#     coefSVD$d <- diag(coefSVD$d[1:nrank])
#     coefSVD$v <- coefSVD$v[, 1:nrank, drop = FALSE]
# 
#     B0 <- coefSVD$u
#     A0 <- coefSVD$v
# 
#   }
# 
# 
#   list(X=X, Y=Y, A0=A0, pvec=pvec)
# 
# }






X.standardize <- function(X, pvec){

  d <- max(pvec)

  X <- scale(X, center = TRUE)

  for( i in 1:d ){
    FNorm <- norm( X[ , which(pvec == i) ], "F" )
    X[ , which(pvec == i) ] <- X[ , which(pvec == i) ] / FNorm
  }

  X
}






print.sim <- function(B, pvec=NULL){
  # pvec

  if( is.null(pvec) ){
    d <- nrow(B) / 10
    pp <- rep(10, d)
  } else {
    pp <- pvec
  }


  selrow <- sort( unlist( lapply(1:2, function(x) c(0, cumsum(pp)[-length(pp)]) + (x) ) ) )
  round( B[selrow, ], 2)
}



DIFF <- function(fit.sRRR){

  lambda.seq <- fit.sRRR$lambda.seq

  out <- NULL
  for( i in 1:length(lambda.seq) ){
    PATH <- fit.sRRR$A.path[[i]]
    diff.A <- sapply( 1:(length(PATH)-1), function(i) norm( PATH[[i+1]] - PATH[[i]], "F" ) )
    PATH <- fit.sRRR$B.path[[i]]
    diff.B <- sapply( 1:(length(PATH)-1), function(i) norm( PATH[[i+1]] - PATH[[i]], "F" ) )
    PATH <- fit.sRRR$C.path[[i]]
    diff.C <- sapply( 1:(length(PATH)-1), function(i) norm( PATH[[i+1]] - PATH[[i]], "F" ) )
    out[[i]] <- cbind.data.frame(A=diff.A, B=diff.B, C=diff.C)
  }


  out
}






AB2Rot <- function(A, B, type=c("A", "B"), rotation="varimax"){
  if(FALSE){
    A <- mydata$A
    B <- mydata$B
  }

  if( type == "A" ){
    rot.mat <- Update.Rotation(A, rotation)$rotation
  } else if( type == "B" ){
    rot.mat <- Update.Rotation(B, rotation)$rotation
  }

  out <- list( AQ = A %*% rot.mat,
               BQ = B %*% rot.mat )

  out
}



Fit2Rot <- function(fit, type, rotation="varimax"){
  if(FALSE){
    fit <- fit1$eta1
  }

  A <- fit$A
  B <- fit$B

  out <- NULL
  for( i in 1:length(A) ){

    if( type == "A" ){
      rot.mat <- Update.Rotation(A[[i]], rotation)$rotation
    } else if( type == "B" ){
      rot.mat <- Update.Rotation(B[[i]], rotation)$rotation
    }

    out[[i]] <- list( AQ = A[[i]] %*% rot.mat,
                      BQ = B[[i]] %*% rot.mat )
  }

  out
}








png.accuracy <- function(true, est){
  if(FALSE){
    true <- Data$A
    est <- fit3.2.1$A[[4]]
  }

  true.01 <- ifelse(true==0, 0, 1)
  est.01 <- ifelse(est==0, 0, 1)

  mean( true.01 == est.01 )
}





acc <- function(X, true, perm=TRUE, rowwise=FALSE){
  X <- as.matrix(X)
  
  mat2binary <- function(mat, eps=1e-10){
    ifelse( abs(mat) > eps, 1, 0 )
  }
  
  if(rowwise){
    X2 <- apply( mat2binary(X), 1, function(x) any(x==1) )
    true2 <- apply( mat2binary(true), 1, function(x) any(x==1) )
    
    out <- mean( X2 == true2 )
    
    return(out)
  }
  
  if( perm ){
    out <- perm.acc(X, true)
  } else {
    out <- mean( mat2binary(true) == mat2binary(X) )
  }
  
  out
}


#' @export acc.best
acc.best <- function(true, est, best = TRUE){
  
  A <- true
  B <- est
  
  permutations <- function(n){
    if(n==1){
      return(matrix(1))
    } else {
      sp <- permutations(n-1)
      p <- nrow(sp)
      A <- matrix(nrow=n*p,ncol=n)
      for(i in 1:n){
        A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
      }
      return(A)
    }
  }
  
  if( best ){
    A <- abs(A);  B <- abs(B)
    
    r <- ncol(B)
    
    perms <- permutations(r)
    
    maximum <- NULL
    # minimum <- 1000
    for( i in 1:nrow(perms) ){
      maximum[i] <- acc(A, B[,perms[i,]], perm=F, rowwise=F)
    }
    
    out <- list( max = max(maximum),
                 cols = unlist(perms[which.max(maximum),]),
                 B = B[,unlist(perms[which.max(maximum),])] )
    
    
    # out <- min(maximum)
    
  } else {
    out <- acc(A, B, perm=F, rowwise=F)
  }
  
  
  out
}






#' @export rmse.best
rmse.best <- function(true, est, best = TRUE){

  A <- true
  B <- est

  permutations <- function(n){
    if(n==1){
      return(matrix(1))
    } else {
      sp <- permutations(n-1)
      p <- nrow(sp)
      A <- matrix(nrow=n*p,ncol=n)
      for(i in 1:n){
        A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
      }
      return(A)
    }
  }

  if( best ){
    A <- abs(A);  B <- abs(B)

    r <- ncol(B)

    perms <- permutations(r)

    minimum <- NULL
    # minimum <- 1000
    for( i in 1:nrow(perms) ){
      # minimum <- min(minimum, sqrt( mean( (A - B[,perms[i,]])^2 ) ))
      minimum[i] <- mean( (A - B[,perms[i,]])^2 )
    }

    out <- list( min = min(minimum),
                 cols = unlist(perms[which.min(minimum),]),
                 B = B[,unlist(perms[which.min(minimum),])] )


    # out <- min(minimum)

  } else {
    out <- mean( (A - B)^2 )
  }


  out
}





rmse.sign.best <- function(true, est){
  n <- nrow(est);  p <- ncol(est)
  
  GRID <- expand.grid(lapply(1:p, function(x) 0:1))
  fnorm <- NULL
  for( iii in 1:nrow(GRID) ){
    est.new <- est %*% diag( as.numeric( ifelse( GRID[iii,] == 1, 1, -1 ) ) )
    fnorm[iii] <- norm(true - est.new, "F")
  }
  
  # wh.min <- which( fnorm == min(fnorm) )
  wh.min <- which.min( fnorm )
  
  if( length(wh.min) == 1 ) opt.sign <- GRID[wh.min,]
  
  est.new <- est %*% diag( as.numeric( ifelse( opt.sign == 1, 1, -1 ) ) )
  
  list( est = est.new, sign = ifelse( opt.sign == 1, 1, -1 ), min=min(fnorm) )
}




crit.varimax <- function(L){
  QL <- sweep(L^2,2,colMeans(L^2),"-")
  sqrt(sum(diag(crossprod(QL))))^2

}

crit.quartimax <- function(L){
  sum(diag(crossprod(L^2)))
}

soft <- function(A, x){
  ifelse( abs(A) > x, A-sign(A-x)*x, 0 )
}











#' @export png.GetStructure
png.GetStructure <- function(d){
  out <- matrix(0, d, 2^d-1)
  
  count <- 0
  for( i in 1:d ){
    Combn <- combn(d,i)
    for( j in 1:ncol(Combn) ){
      count <- count+1
      out[Combn[,j],count] <- 1
    }
  }
  
  cbind(0, out)
}






#' @export permutations
permutations <- function(n){
  if(n==1){
    return(matrix(1))
  } else {
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}
