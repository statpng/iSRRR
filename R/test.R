#' @export test0
test0 <- function(){
  # devtools::load_all()
  
  library(tidyverse)
  library(iSRRR)
  
  n.seq <- c(10,30,50,100)
  p.seq <- c(100, 500)
  q.seq <- c(50)
  snr.seq <- c(0.5, 1.0)
  
  GRID <- expand.grid(n=n.seq, p=p.seq, q=q.seq, snr=snr.seq)
  
  II <- 1
  
  nthread <- 2
  AA <- matrix( 1:nrow(GRID), ncol=nthread, byrow=FALSE )[,II]
  for( i in AA ){
    
    # n=50; p=500; q=20; d=4; rvec=rep(2,3); snr=1
    # seq_typeA <- c("sparse0.2", "BlockDiagonal")[2]
    # seq_typeB <- c("row", "individual", "partial", "all")
    
    n=GRID[i,"n"]; p=GRID[i,"p"]; q=GRID[i,"q"]; snr=GRID[i,"snr"]
    d=4; rvec=rep(2,3)
    
    ( title1 <- paste0("i=",i,";", paste0(paste0( c("n","p","q","snr"), "=", GRID[i,] ), collapse=";")) )
    
    set.seed(1)
    data1 <- simdata2(typeA="quartimax", typeB="individual", n=n, d=d, q=q, p=p, rvec=rvec, snr=snr)
    set.seed(2)
    data2 <- simdata2(typeA="quartimax", typeB="partial", n=n, d=d, q=q, p=p, rvec=rvec, snr=snr)
    set.seed(3)
    data3 <- simdata2(typeA="quartimax", typeB="all", n=n, d=d, q=q, p=p, rvec=rvec, snr=snr)
    
    
    fit_rsd <- fit_hard <- NULL
    for( idx in 1:3 ){
      
      data <- list(data1, data2, data3)[[idx]]
      # fit_rsd[[idx]] <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec,
      #                   nrank = data$nrank,
      #                   nu = c(0, 0.2, 0.5, 0.7, 1.0),
      #                   typeA="RSD",
      #                   params=list(nlambda=10, lambda.factor=1e-4),
      #                   control=list(best=TRUE, maxit.mse=50, eps.mse=1e-10, 
      #                                maxit.B=1e6, eps.B=1e-6), 
      #                   use.gglasso=TRUE,
      #                   trueA = data$A, trueB = data$B,
      #                   rsd.it = 50 )
      # 
      # save(fit_rsd[[idx]], file=paste0("./fit_rsd - ", c("row", "ind", "partial", "all")[idx], ";",title1,".RData"))
      
      
      fit <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec,
                         nrank = data$nrank,
                         nu = (0:4*0.1),
                         params=list(nlambda=10, lambda.factor=1e-4),
                         control=list(best=TRUE, maxit.mse=20, eps.mse=1e-6, 
                                      maxit.B=1e6, eps.B=1e-6),
                         trueA = data$A, trueB = data$B )
      
      
      
      save(fit_hard[[idx]], file=paste0("./fit_hard - ", c("row", "ind", "partial", "all")[idx], ";",title1,".RData"))
      
    }
    
    save(fit_rsd, fit_hard, file=paste0("./fit - ", c("row", "ind", "partial", "all")[idx], ";",title1,".RData"))
    
  }
  
  
  
  
  fit0 <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec,
                nrank = data$nrank,
                nu = (0:4*0.1),
                params=list(lambda.seq=c(0,0.5) ),
                control=list(best=TRUE, maxit.mse=50, eps.mse=1e-12, 
                             maxit.B=1e6, eps.B=1e-6),
                trueA = data$A, trueB = data$B )
  
  
  
  
  # lapply(fit_hard, check_constraints)
  # check_constraints(fit_hard[[1]])
  # fit_hard[[1]]$diff.mse[[1]]
  # fit_hard[[1]]$B[[9]] %>% print.B2()
  
}






  
test1 <- function(){ 

  Algorithm <- c("RSD", "hard")[1]
  
  nu.seq <- list( c(0, 0.2, 0.5), 
                  c(0, 0.2, 0.4) )[[1]]
  seq_rsd.it <- c(20, 50)
  seq_typeA <- c("sparse0.2", "BlockDiagonal")[2]
  seq_typeB <- c("row", "individual", "partial", "all")
  
  GRID <- expand.grid(seq_rsd.it, seq_typeA, seq_typeB)
  
  i=13
  i=5
  for( i in 1:nrow(GRID) ){
    print( paste0( i, " / ", nrow(GRID) ) )
    
    (rsd.it <- GRID[i,1])
    (typeA <- GRID[i,2])
    (typeB <- GRID[i,3])
    
    # r <- switch(as.character(typeB), "partial" = 3)
    
    set.seed(1)
    data <- simdata2(typeA=typeA, typeB=typeB, n=n, d=d, q=q, pvec=rep(p,d), rvec=rep(1,r), d0=3, es="1", simplify=TRUE, snr=snr, rho_X=0.5)
    
    
    fit2 <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec,
                   nrank = r,
                   nu = 0:4*0.1,
                   typeA="hard",
                   params=list(nlambda=10, lambda.factor=1e-4),
                   control=list(best=TRUE, maxit.mse=1000, eps.mse=1e-12, 
                                maxit.B=3e6, eps.B=1e-6),
                   use.gglasso=TRUE, 
                   trueA = data$A, trueB = data$B )
    save(fit2, file="./RSD_fit2.RData")
    
    
    
    fit3_1 <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec,
                     nrank = r,
                     nu = c(0, 0.2, 0.5, 0.7, 1.0),
                     typeA="RSD",
                     params=list(nlambda=10, lambda.factor=1e-4),
                     control=list(best=TRUE, maxit.mse=50, eps.mse=1e-10, 
                                  maxit.B=1e6, eps.B=1e-6), 
                     use.gglasso=TRUE,
                     trueA = data$A, trueB = data$B,
                     rsd.it = 50 )
    save(fit3_1, file="./RSD_fit3_1.RData")
    
    fit3_2 <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec,
                     nrank = r,
                     nu = c(0, 0.2, 0.5, 0.7, 1.0),
                     typeA="RSD",
                     params=list(nlambda=10, lambda.factor=1e-4),
                     control=list(best=TRUE, maxit.mse=50, eps.mse=1e-10, 
                                  maxit.B=1e6, eps.B=1e-6), 
                     use.gglasso=TRUE,
                     trueA = data$A, trueB = data$B,
                     rsd.it = 100 )
    save(fit3_2, file="./RSD_fit3_2.RData")
    
    fit3_3 <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec,
                     nrank = r,
                     nu = c(0, 0.2, 0.5, 0.7, 1.0),
                     typeA="RSD",
                     params=list(nlambda=10, lambda.factor=1e-4),
                     control=list(best=TRUE, maxit.mse=100, eps.mse=1e-12, 
                                  maxit.B=1e6, eps.B=1e-6), 
                     use.gglasso=TRUE,
                     trueA = data$A, trueB = data$B,
                     rsd.it = 50 )
    save(fit3_3, file="./RSD_fit3_3.RData")
    
    
    fit1[[1]]$time
    fit2[[1]]$time
    fit3[[1]]$time
    
    fit2[[1]] %>% check_constraints()
    fit2[[2]] %>% check_constraints()
    fit2[[3]] %>% check_constraints()
    fit2[[4]] %>% check_constraints()
    fit2[[5]] %>% check_constraints()
    
    
    fit3_1[[1]] %>% check_constraints()
    fit3_1[[2]] %>% check_constraints()
    fit3_1[[3]] %>% check_constraints()
    fit3_1[[4]] %>% check_constraints()
    fit3_1[[5]] %>% check_constraints()
    
    
    fit3[[1]] %>% check_constraints()
    
    
    fit2[[4]]$A[[4]] %>% round(3)
    
    fit2[[1]] %>% check_constraints()
    fit3[[1]] %>% check_constraints()
    
    fit2[[1]] %>% png.IC() %>% .$BIC %>% .[,1] %>% min()
    fit2[[2]] %>% png.IC() %>% .$BIC %>% .[,1] %>% min()
    fit2[[3]] %>% png.IC() %>% .$BIC %>% .[,1] %>% min()
    fit2[[4]] %>% png.IC() %>% .$BIC %>% .[,1] %>% min()
    fit2[[5]] %>% png.IC() %>% .$BIC %>% .[,1] %>% min()
    
    fit3_1[[1]] %>% png.IC() %>% .$BIC %>% .[,1] %>% min()
    fit3_1[[2]] %>% png.IC() %>% .$BIC %>% .[,1] %>% min()
    fit3_1[[3]] %>% png.IC() %>% .$BIC %>% .[,1] %>% min()
    fit3_1[[4]] %>% png.IC() %>% .$BIC %>% .[,1] %>% min()
    fit3_1[[5]] %>% png.IC() %>% .$BIC %>% .[,1] %>% min()
    
    fit3[[2]]$B[[4]] %>% print.B2()
    
    
    fit2[[1]]$B[[7]] %>% print.B2()
    fit2[[2]]$B[[7]] %>% print.B2()
    fit2[[3]]$B[[7]] %>% print.B2()
    
    fit1[[1]]$diff.mse[[5]]
    fit2[[5]]$diff.mse[[2]][,4]
    
    fit2[[1]]$B[[3]] %>% print.B2()
    
    
    fit.rssvd <- rrpack::rssvd(data$Y, data$X, 6, control=list(maxit=500))
    
    res.path <- lapply(1:100, function(i){
      fit.rssvd$U.path[,i,] %*% 
        diag(fit.rssvd$D.path[i,]) %*%
        t(fit.rssvd$V.path[,i,])
    })
    
    res.diff <- NULL
    for( i in 1:99 ){
      res.diff[i] <- mean( (res.path[[i]] - res.path[[i+1]])^2 )
    }
    
    rev(res.diff) %>% plot(type="b")
    #
    
    
    fit2[[1]]$diff.mse[[2]][,4] %>% plot(type="b")
    fit2[[1]]$diff.mse[[3]][,4] %>% plot(type="b")
    
    fit3[[1]]$diff.mse[[4]][,4] %>% plot(type="b")
    
    
    fit1[[2]] %>% check_constraints()
    fit1[[3]] %>% check_constraints()
    
    fit1[[1]]$B[[5]] %>% print.B2()
    fit0[[2]]$B[[5]] %>% print.B2()
    fit0[[3]]$B[[5]] %>% print.B2()
    
    fit1[[1]]$A[[5]] %>% round(3)
    fit0[[2]]$A[[5]] %>% round(3)
    fit0[[3]]$A[[5]] %>% round(3)
    
    fit0[[1]]$diff.mse[[3]][,4] %>% plot(type="b")
    fit0[[2]]$diff.mse[[3]][,4] %>% plot(type="b")
    fit0[[3]]$diff.mse[[3]][,4] %>% plot(type="b")
    
    
    
    
    
    
    
    print("---------------------------")
    print("---------------------------")
    
    print("1 ---------------------------")
    print( check_constraints(fit[[1]]) )
    print("2 ---------------------------")
    print( check_constraints(fit[[2]]) )
    print("3 ---------------------------")
    print( check_constraints(fit[[3]]) )
    print("4 ---------------------------")
    print( check_constraints(fit[[4]]) )
    print("5 ---------------------------")
    print( check_constraints(fit[[5]]) )
    print("---------------------------")
    print("---------------------------")
    
    BIC.mat <- sapply(fit, function(x) round(png.IC(x)$BIC[,1],4) )
    colnames(BIC.mat) <- nu.seq
    rownames(BIC.mat) <- paste0("lam=", format(fit[[1]]$params$lambda.seq, digit=2, scientific = TRUE))
    print(BIC.mat)
    
    print("---------------------------")
    print("---------------------------")
    
    
    save(data, fit, BIC.mat, file=sprintf("./iSRRR - %s - (%s, %s); nu=(%s); rsd.it=%s.RData", Algorithm, typeA, typeB, paste0(nu.seq,collapse=","), rsd.it) )
    
  }
  
  
  
}




test2 <- function(){
  devtools::load_all()
  
  n=50; p=50; q=20; d=4; r=6;
  params <- list(nlambda=10, lambda.factor=1e-4)
  control <- list(best=TRUE, maxit.mse=100, eps.mse=1e-10, maxit.B=1e8, eps.B=1e-8)
  
  Algorithm <- c("RSD", "hard")[2]
  
  nu.seq <- list( c(0, 0.2, 0.5, 0.8, 1.0), 
                  c(0, 0.1, 0.2, 0.3, 0.4) )[[2]][c(1,3,5)]
  seq_rsd.it <- c(100, 500)[1]
  seq_typeA <- c("sparse0.0", "BlockDiagonal")
  seq_typeB <- c("row", "individual", "partial", "all")
  seq_snr <- c(0.5, 1)
  
  GRID <- expand.grid(seq_rsd.it, seq_typeA, seq_typeB, seq_snr)
  
  # i=7
  for( i in 1:nrow(GRID) ){
    print( paste0( i, " / ", nrow(GRID) ) )
    
    (rsd.it <- GRID[i,1])
    (typeA <- GRID[i,2])
    (typeB <- GRID[i,3])
    (snr <- GRID[i,4])
    
    set.seed(1)
    data <- simdata2(typeA=typeA, typeB=typeB, n=n, d=d, q=q, pvec=rep(p,d), rvec=rep(1,r), d0=3, es="1", simplify=TRUE, snr=snr, rho_X=0.5)
    
    # data <- simdata3_final(typeA=typeA, typeB=typeB, n=n, p=p, q=q, snr=snr, nu.A=0.05)
    
    #
    fit0 <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec,
                   nrank = r,
                   # nu = c(0,0.2,0.4),
                   nu = nu.seq[1],
                   typeA=Algorithm,
                   params=params,
                   control=control, use.gglasso=TRUE, 
                   rsd.it = rsd.it )
    
    
    
    fit0[[1]] %>% check_constraints()
    fit0[[2]] %>% check_constraints()
    fit0[[3]] %>% check_constraints()
    
    fit0[[1]]$B[[4]] %>% print.B2()
    fit0[[2]]$B[[5]] %>% print.B2()
    fit0[[3]]$B[[5]] %>% print.B2()
    
    fit0[[1]]$A[[5]] %>% round(3)
    fit0[[2]]$A[[5]] %>% round(3)
    fit0[[3]]$A[[5]] %>% round(3)
    
    fit0[[1]]$diff.mse[[3]][,4] %>% plot(type="b")
    fit0[[2]]$diff.mse[[3]][,4] %>% plot(type="b")
    fit0[[3]]$diff.mse[[3]][,4] %>% plot(type="b")
    
    
    
    fit0[[1]] %>% png.IC()
    fit0[[1]]$B[[10]] %>% print.B2()
    
    #
    
    
    for( ir in 1:4 ){
      fit[[ir]] <- iSRRR( X = data$X, Y = data$Y, pvec = data$pvec, 
                          nrank = ir+1,
                          # nu = c(0,0.2,0.4),
                          nu = nu.seq,
                          typeA=Algorithm,
                          params=params,
                          control=control,
                          trueA=data$A, use.gglasso=TRUE, rsd.it = rsd.it )
    }
    
    
    print("---------------------------")
    print("---------------------------")
    
    print("1 ---------------------------")
    print( check_constraints(fit[[1]]) )
    print("2 ---------------------------")
    print( check_constraints(fit[[2]]) )
    print("3 ---------------------------")
    print( check_constraints(fit[[3]]) )
    print("4 ---------------------------")
    print( check_constraints(fit[[4]]) )
    print("5 ---------------------------")
    print( check_constraints(fit[[5]]) )
    print("---------------------------")
    print("---------------------------")
    
    get.BICmat <- function(fit.new){
      BIC.mat <- sapply(fit.new, function(x) round(png.IC(x)$BIC[,1],4) )
      colnames(BIC.mat) <- nu.seq
      rownames(BIC.mat) <- paste0("lam=", format(fit.new[[1]]$params$lambda.seq, digit=2, scientific = TRUE))
      BIC.mat
    }
    
    get.BICmat(fit[[1]]) %>% min()
    get.BICmat(fit[[2]]) %>% min()
    get.BICmat(fit[[3]]) %>% min()
    get.BICmat(fit[[4]]) %>% min()
    
    print("---------------------------")
    print("---------------------------")
    
    
    save(data, fit, file=sprintf("./iSRRR - %s - (%s, %s); nu=(%s); rsd.it=%s.RData", Algorithm, typeA, typeB, paste0(nu.seq,collapse=","), rsd.it) )
    
  }
  
  
  
}







RCG <- function(){
  
  set.seed(1234)
  
  n=50; p=50; q=20; r=6; 
  
  tau <- 10
  
  
  
  
  out.rcg <- png.RCG(tau=2)
  out.rcg$A %>% round(5)
  data$A %>% round(5)
  
  cbind(
    apply(out.rcg$A, 1, norm, "2") %>% round(3) %>% order(decreasing=TRUE),
    apply(data$A, 1, norm, "2") %>% round(3) %>% order(decreasing=TRUE)
  )
  
  cbind(
    apply(out.rcg$A, 1, norm, "2") %>% round(3),
    apply(data$A, 1, norm, "2") %>% round(3)
  )
  
  
  out.rcg %>% GPArotation::quartimax() %>% {round(.$loadings,5)}
  data$A %>% GPArotation::quartimax()
  
  
  
  data$A %>% {cbind(round(.,5), round(sqrt(rowSums(.^2)),5) ) } %>% {ifelse(.<0.35,0,.)}
  png.RCG(tau=1) %>% round(5) %>% apply(1,norm,"2") %>% order(decreasing=TRUE)
  
  data$A %>% apply(1,norm,"2") %>% order(decreasing=TRUE)
  
  
  head(tx(res$xopt))
  
  
  
  
}