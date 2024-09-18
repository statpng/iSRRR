test <- function(){
  # devtools::load_all()
  
  library(dplyr)
  library(iSRRR)
  
  data(iSRRR.data)
  data(iSRRR.boot)
  
  X.list2 <- lapply(X.list, function(x) scale(x[,1:100],scale=F) %>% {./norm(.,"F")} )
  Y <- scale(Y.list$GI50)
  
  sel.list <- lapply(X.list2, function(X) mSIS(X,Y))
  
  X.list.new <- lapply(1:length(X.list2), function(i) X.list2[[i]][,sel.list[[i]]])
  
  X.list.new %>% str
  
  X <- do.call("cbind", X.list.new)
  Y <- Y
  pvec <- sapply(X.list.new, ncol) %>% {rep(1:length(.), .)}
  
  #
  # STRS::STRS(Y, X)
  # iSRRR::BICk(Y, X)
  
  nrank=2
  #
  
  fit <- iSRRR(X, Y, pvec=pvec, nrank=nrank, nu=0.1, control=list(maxit=5))
  
  iSRRR.lambda(X, Y, pvec=pvec, nrank=nrank)$lambda.seq
  fit[[1]]$param$lambda.seq
  
  #
  
  
  cor(Y) %>% corrplot(method="square", tl.pos="n", cl.pos="b", cl.ratio=0.06)
  
  
  library(dplyr)
  data(iSRRR.fit)
  C1 <- iSRRR.choose.LambdaNu(fit, 3, 1) %>% with(tcrossprod(B,A))
  C2 <- iSRRR.choose.LambdaNu(fit, 3, 2) %>% with(tcrossprod(B,A))
  C3 <- iSRRR.choose.LambdaNu(fit, 3, 3) %>% with(tcrossprod(B,A))
  C.List <- list(C1,C2,C3)
  
  iSRRR.plot.C.list(C.List)
  
  
  iSRRR.check.constraints(fit[[4]])
  #
  
  
  fit %>% png.BICmat(print.out=F)
  AB.opt <- print.iSRRR( fit, 2, 3 )
  AB.opt$B %>% print.B2()
  AB.opt$A[,2] %>% round(3) %>% {which(.!=0)}
  AB.opt$A[abs(AB.opt$A)<1e-8] <- 0
  
  Chat <- with(AB.opt, tcrossprod(B, A))
  
  
  
  #
  iSRRR.C.heatmap(Chat, legend.bar.width=20)
  ggsave(filename = "NCI60.plot2.C.rank2.(lam,nu)=(2,3).eps", width=5, height=3.5)
  #
  
  
  
  data(iSRRR.fit)
  fit %>% png.BICmat(print.out=F)
  
  print.iSRRR(fit, 3, 3)$B %>% print.B2()
  
  
  
  data(iSRRR.boot)
  boot.array <- iSRRR::iSRRR.boot.selarray(out.boot)
  boot.array/100
  boot.array[2,,3]/100
  
  SelFreq.Chat <- summary.boot.Chat(out.boot)
  iSRRR.SelFreq.C.heatmap(SelFreq.Chat)
  ggsave(filename="NCI60.plot.SelFreq.Chat.rank2.eps", width=5, height=2.5)
  #
  
  
}








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
  
  
  
}







