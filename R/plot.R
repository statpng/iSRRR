#' @import ggplot2
#' @import plyr
#' @import plot.matrix
#' 
#' @export png.iSRRR.plot
png.iSRRR.plot <- function(List, filename=NULL, print=TRUE, order=TRUE, order.ref=1, order.type=c("col", "row", "both", "none"), width=10, height=20, legend.bar.width=40, legend.bar.height=0.3){

  if(FALSE){
    filename=NULL
    print=FALSE
    order=TRUE
    order.ref=1
    order.type="col"
    width=10
    height=20
    legend.bar.width=40
    legend.bar.height=0.3
  }

  # B.list <- list(iSRRR0 = iSRRR.fit$B[[1]],
  #                iSRRR0.1 = iSRRR.fit$B[[2]],
  #                iSRRR0.2 = iSRRR.fit$B[[3]],
  #                srrr = fit.others$out$srrr$B,
  #                # rssvd = with(fit.others5$out$rssvd,U%*%D),
  #                sofar = with(fit.others$out$sofar,U%*%D),
  #                secure = with(fit.others$out$secure,U%*%D) )

  if(is.null(filename)) filename <- deparse(substitute(List))

  Names.List <- names(List)
  B.list <- List

  Array <- simplify2array(List)
  dimnames(Array) <- list(NULL, NULL, Names.List)

  ord_row <- hclust( dist(Array[,,order.ref]) )$order
  ord_col <- hclust( dist(t(Array[,,order.ref])) )$order

  if( order.type == "both" ){
    DF <- adply( Array[ord_row,ord_col,], 1:3 )
  } else if ( order.type == "col" ){
    DF <- adply( Array[,ord_col,], 1:3 )
  } else if ( order.type == "row" ){
    DF <- adply( Array[ord_row,,], 1:3 )
  } else {
    DF <- adply( Array[,,], 1:3 )
  }


  colnames(DF) <- c("x", "y", "method", "value")
  DF$x <- as.numeric(as.character(DF$x))
  DF$y <- as.numeric(as.character(DF$y))
  DF$method <- factor(DF$method,
                        levels = c("iSRRR0.0", "iSRRR0.1", "iSRRR0.2"),
                        labels = c( bquote(nu == 0.0), bquote(nu == 0.1), bquote(nu == 0.2) ) )

  # DF |> head()

  #
  #
  #

  DF.new <- DF

  for( idx in 1:nlevels(DF.new$method) ){
    x <- DF.new[ DF$method==levels(DF.new$method)[idx], ]$value

    x[ x > quantile(x, 0.99)] <- max(x) - (max(x) - quantile(x, 0.99))*0.9
    x[ x < quantile(x, 0.01)] <- min(x) + (quantile(x, 0.01) - min(x))*0.9

    DF.new[ DF.new$method==levels(DF.new$method)[idx], ]$value <- x
  }

  #

  plt <- cbind.data.frame(DF.new, w=2, h=2) |>
    # dplyr::filter( grepl("iSRRR",method) ) |>
    ggplot(aes(x=y,y=max(x)-x+1,fill=value)) +
    # ggplot(aes(y=x,x=y,fill=ifelse(abs(value)>1e-10,1,0)))+
    # geom_tile(color="gray80", size=0.001) +
    geom_tile(color=NA, size=0.001) +
    scale_fill_gradient2(low = "#AB3027",
                         mid = "white",
                         # mid = "#D2FBC5",
                         # mid = "#FFFFC2",
                         high = "#0F1BBF") +
    # scale_fill_gradientn(
    #   colours = c("#AB3027", "red", "white", "blue", "#0F1BBF"),
    #   values = c(-10, -3, 0, 3, 10)
    # ) +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw(base_size = 18) +
    facet_grid(~method, labeller = label_parsed) +
    # coord_fixed( ratio=1 ) +
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank() ) +
    guides(fill = guide_colorbar(title = "",
                                 label.position = "bottom",
                                 title.position = "left",
                                 title.vjust = 1,
                                 frame.colour = "black",
                                 barwidth = legend.bar.width,
                                 barheight = legend.bar.height))

  if(print){
    ggsave(filename=sprintf("%s.pdf", filename), plt, width=width, height=height)
  } else {
    print( plt )
  }

}




#' @export png.matplot
png.matplot <- function(x, position=NULL, label=NULL, cols=NULL, lty=1, lwd=2, legend.print=TRUE, ...){
  # x = fit$diff[[idx]]
  # sel.cols <- c(1,4)

  if(is.null(position)) position <- "topright"
  if(is.null(label)) label <- colnames(x)
  if(is.null(cols)) cols <- seq_len(ncol(x))

  matplot(x[,cols], type="l", col=cols, lty=lty, lwd=lwd, xlab=NA, ylab=NA, ...)
  if(legend.print){
    legend(position, label[cols], col=cols, lty=lty, lwd=lwd)
  }


}





#' @export png.plot.matrix
png.plot.matrix <- function(mat, binarize=FALSE, const=1, nbreaks=10, main=NA, ...){

  if(FALSE){
    mat <- fit$B[[4]]
    type="binary"
    nbreaks=10
    main=NA
    ... <- NULL
  }



  n <- nrow(mat)
  p <- ncol(mat)



  # if(type=="binary"){
    # mat <- ifelse( abs(mat) == 0 , 0, 1)
  # }
  sd.mat <- sd( as.numeric(mat) )
  mat <- ifelse( mat > quantile(mat, 0.75) + const*sd.mat, max(mat), mat)
  mat <- ifelse( mat < quantile(mat, 0.25) - const*sd.mat, min(mat), mat)

  MAX <- max( abs( min(mat) ), abs( max(mat) ) )
  # rg <- seq(-MAX, MAX, length.out=nbreaks) |> round(2)
  rg <- seq(0, MAX, length.out=nbreaks) |> round(2)
  # rg <- exp( seq( log(1e-10), log(MAX), length.out=10 ) ) |> round(2)



  # mat.new <- which( ifelse( abs(mat) == 0 , 0, 1) == 1, arr.ind=TRUE )
  # mat.new <- mat.new[,c(2,1)]
  #
  #
  # plot(mat.new[,1], mat.new[,2], pch=15,
  #      # xaxt="n",
  #      # yaxt="n",
  #      ylim=rev(range(mat.new[,2])),
  #      xlab=NA, ylab=NA, cex=cex)



  if(binarize){
    mat <- ifelse( mat!=0, 1, 0 )
    rg <- c(0,1)
  }

  mat <- ifelse( abs(mat) == 0 , NA, mat)


  
  # plot( mat, breaks=rg, col=colorRampPalette(c("red", "white", "blue")), main=main, xlab="Factors", ylab=NA, fmt.key="%.2f", na.col="white", border=NA, ...) #(nbreaks), ...)

  plot( abs(mat), breaks=rg, col=colorRampPalette(c("white", "black")), main=main, xlab="Factors", ylab=NA, fmt.key="%.2f", na.col="white", border=NA, key=NULL, axis.row=NULL, ...) #(nbreaks), ...)


}





png.plot.matrix.combined <- function(fit, title=NULL){
  if(is.null(title)) title <-deparse(substitute(fit))

  pdf(file=paste0(title, ".pdf"), width=7, height=3)
  par(mar=c(5.1, 1.1, 4.1, 0.1))
  par(mai=c(0.8,0.5,0.3,0.4))
  par(omi=c(0.1,0.1,0.1,0.4))
  par(mfrow=c(1,2))
  png.plot.matrix(mat=fit$B, nbreaks=11, main="B", cex.axis=0.8)
  png.plot.matrix(mat=fit$A, nbreaks=11, main="A", axis.row=NULL, cex.axis=0.8)
  dev.off()
}





plot.png.sRRR <- function(fit, type="B"){

  if(FALSE){
    fit <- fit.Rot.Lasso1.0
    type <- "B"
  }

  Data <- fit$data

  

  acc <- NULL
  for( i in 1:length(fit$B) ){
    # B.best <- rmse.best(Data$B, fit$B[[i]], best=TRUE)$B
    acc[i] <- png.accuracy(Data$B, fit$B[[i]])
  }

  opt.idx <- which.max( acc )
  opt.rmse <- rmse.best(Data$B, fit$B[[opt.idx]], best=FALSE)

  fit.opt <- list( A=fit$A[[opt.idx]], B=fit$B[[opt.idx]] )
  title <- paste0(deparse(substitute(fit)), "(rmse=", round(opt.rmse,4), ")")


  png.plot.matrix.combined(fit.opt, title=title)

}




#' @export png.image
png.image <- function(mat, binary=FALSE, col.type=2){
  if(FALSE){
    mat <- fit[[1]]$B[[3]][66:75,1:10]
    pvec <- rep(1:3, rep(10,3))
  }

  # mat <- abs(mat)

  colnames(mat) <- paste0("F", 1:ncol(mat))
  mat <- t(mat)
  pvec0 <- attr(mat, "pvec")
  pvec <- rep(1:length(pvec0), pvec0)



  if(binary){
    mat <- ifelse( abs(mat)>0, 1, 0 )
  } else {
    mat[mat>0] <- 1
    mat[mat<0] <- -1
    # mat <- abs(mat)
    # mat <- png.cut( mat, 9 )
    # mat <- mat-mean(mat)
    # mat <- mat - median(as.numeric(names(table(mat))))
  }

  # mat <- mat - min(mat) / (max(mat)-min(mat))
  # mat <- png.trim(mat)

  # image(1:nrow(mat), 1:ncol(mat), pmin(mat, 1), axes=F, col=c("white", "black"), xlab = "", ylab = "")


  # hcl.pals()
  col <- switch(col.type,
                `1` = c("black", "white", "black"),
                `2` = c("khaki1", "white", "midnightblue"),
                `3` = c("darkorange", "white", "midnightblue"),
                `4` = hcl.colors(3, "PuOr")
  )


  image(1:nrow(mat), (1:ncol(mat)), mat[,rev(1:ncol(mat))], axes=F,
        col = col,
        xlab = "", ylab = "")
  axis(1, at = 1:nrow(mat), labels=rownames(mat), tick=F)
  axis(2, at = cumsum( table(pvec) |> (function(.) c(.[1]/2,.[-1]) )() ), labels=names(table(pvec)), tick=F)


}

