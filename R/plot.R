#' @import ggplot2
#' @importFrom plyr adply
#' @import plot.matrix
#' 
#' @export iSRRR.plot.C.list
iSRRR.plot.C.list <- function(List, filename=NULL, print=TRUE, order=TRUE, order.ref=1, order.type=c("none", "col", "row", "both"), width=10, height=20, legend.bar.width=40, legend.bar.height=0.3){
  
  if( length(order.type)>1 ) order.type <- "none"

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
  # DF$method <- factor(DF$method,
  #                       levels = c("iSRRR0.0", "iSRRR0.1", "iSRRR0.2"),
  #                       labels = c( bquote(nu == 0.0), bquote(nu == 0.1), bquote(nu == 0.2) ) )

  

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



#' @importFrom corrplot corrplot
#' @export iSRRR.Y.corrplot
iSRRR.Y.corrplot <- function(Y){
  corrplot::corrplot(Y, method="square", tl.pos="n", cl.pos="b", cl.ratio=0.06)
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




















#' @export Chat.remove.outlier
Chat.remove.outlier <- function(Chat){
  x <- Chat
  x[ x > quantile(x, 0.99)] <- max(x) - (max(x) - quantile(x, 0.99))*0.5
  x[ x < quantile(x, 0.01)] <- min(x) + (quantile(x, 0.01) - min(x))*0.5
  x
}


#' @export iSRRR.C.heatmap
iSRRR.C.heatmap <- function(Chat, 
                            legend.bar.width=20, 
                            legend.bar.height=0.3){
  dimnames(Chat) <- NULL
  
  Chat.df <- adply( as.array(Chat), 1:2 )
  colnames(Chat.df) <- c("x", "y", "value")
  Chat.df$x <- as.numeric(as.character(Chat.df$x))
  Chat.df$y <- as.numeric(as.character(Chat.df$y))
  
  
  DF.new <- Chat.df
  DF.new$value <- Chat.remove.outlier(DF.new$value)
  
  
  # legend.bar.width <- 30
  # legend.bar.height <- 0.3
  
  
  plt <- DF.new %>% 
    cbind.data.frame(w=2, h=2) |>
    # dplyr::filter( grepl("ssrrr",method) ) |>
    ggplot(aes(x=y,y=max(x)-x+1,fill=value)) +
    # ggplot(aes(y=x,x=y,fill=ifelse(abs(value)>1e-10,1,0)))+
    # geom_tile(color="grey80", size=0.001) +
    geom_tile(color=NA, size=0.001) +
    # scale_fill_gradient2(low = "#AB3027",
    #                      mid = "white",
    #                      # mid = "#D2FBC5",
    #                      # mid = "#FFFFC2",
    #                      high = "#0F1BBF") +
    # scale_fill_binned(type="viridis") +
    # scale_fill_steps(low = "red",
    #                      # mid = "white",
    #                      high = "blue") +
    # scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0) + 
    scale_fill_gradient2(low = "red",
                         mid = "white",
                         high = "blue") +
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          axis.ticks = element_blank(),
          # panel.border = element_blank(),
          axis.text = element_blank() ) +
    guides(fill = guide_colorbar(title = "",
                                 label.position = "bottom",
                                 title.position = "left",
                                 title.vjust = 1,
                                 frame.colour = "black",
                                 barwidth = legend.bar.width,
                                 barheight = legend.bar.height))
  
  
  plt +
    annotate("text", x = -2.8, y = sum(c(558, 476, 409, 502/2)), label = "Exp", angle=90) +
    annotate("text", x = -2.8, y = sum(c(558, 476, 409/2)), label = "Mut", angle=90) +
    annotate("text", x = -2.8, y = sum(c(558, 476/2)), label = "Met", angle=90) +
    annotate("text", x = -2.8, y = sum(c(558/2)), label = "Cop", angle=90) +
    theme(plot.margin = margin(.6,.6,.6,.6, "cm")) +
    coord_cartesian(xlim=c(-1,53), clip = "off")
  
  # dim(Y)
  # # [1] 57 53
  # X.list %>% names
  # # [1] "Exp" "Mut" "Met" "Cop"
  # # 502 409 476 558 
  
  
}




#' @export summary.boot.Chat
summary.boot.Chat <- function(fit.boot, wh.tune=c(2,3)){
  
  opt.lam <- wh.tune[1]
  opt.nu <- wh.tune[2]
  nboot <- length(fit.boot)
  
  for( i.boot in 1:nboot ){
    A.opt <- fit.boot[[i.boot]][[opt.nu]]$A[[opt.lam]]
    A.opt[abs(A.opt)<1e-15] <- 0
    B.opt <- fit.boot[[i.boot]][[opt.nu]]$B[[opt.lam]]
    Chat <- tcrossprod(B.opt, A.opt)
    
    Chat2 <- ifelse( abs(Chat)==0, 0, 1 )
    
    if( i.boot == 1 ){
      out <- Chat2
    } else {
      out <- out+Chat2
    }
  }
  
  
  return( out/100 )
  # image(out/100, col = gray.colors(10))
  # heatmap(out/100, Rowv=NA, Colv=NA)
  
}



#' @export iSRRR.SelFreq.C.heatmap
iSRRR.SelFreq.C.heatmap <- function(SelFreq.Chat){
  
  dimnames(SelFreq.Chat) <- NULL
  SelFreq.Chat.df <- adply( as.matrix(SelFreq.Chat), 1:2 )
  # heatmap( SelFreq.Chat, Rowv=NA, Colv=NA, col=rev(gray.colors(10)) )
  colnames(SelFreq.Chat.df) <- c("x", "y", "value")
  SelFreq.Chat.df$x <- as.numeric(as.character(SelFreq.Chat.df$x))
  SelFreq.Chat.df$y <- as.numeric(as.character(SelFreq.Chat.df$y))
  
  legend.bar.width=20
  legend.bar.height=0.2
  plt <- SelFreq.Chat.df %>% 
    cbind.data.frame(w=2, h=2) |>
    ggplot(aes(x=y,y=max(x)-x+1,fill=value)) +
    geom_tile(color=NA, size=0) +
    # scale_fill_gradient2(low = "white",
    #                      high = "black") +
    # scale_fill_gradient2(low = "#AB3027",
    #                      mid = "white",
    #                      # mid = "#D2FBC5",
    #                      # mid = "#FFFFC2",
    #                      high = "#0F1BBF") +
    # scale_fill_binned(type="viridis") +
    # scale_colour_stepsn(colours = terrain.colors(10)) +
    scale_fill_steps(breaks=c(0:10*0.1),
                     low = "white",
                     high = "black") +
    # scale_fill_continuous_divergingx(palette = 'RdBu', mid = 0) + 
    scale_x_discrete(expand = c(0,0)) +
    scale_y_discrete(expand = c(0,0)) +
    theme_bw(base_size = 10) +
    theme(legend.position = "bottom",
          axis.title = element_blank(),
          axis.ticks = element_blank(), 
          # panel.border = element_blank(),
          axis.text = element_blank() ) +
    guides(fill = guide_colorbar(title = "",
                                 label.position = "bottom",
                                 title.position = "left",
                                 title.vjust = 1,
                                 frame.colour = "black",
                                 barwidth = legend.bar.width,
                                 barheight = legend.bar.height))
  
  plt +
    annotate("text", x = -1.2, y = sum(c(558, 476, 409, 502/2)), label = "Exp", angle=90) +
    annotate("text", x = -1.2, y = sum(c(558, 476, 409/2)), label = "Mut", angle=90) +
    annotate("text", x = -1.2, y = sum(c(558, 476/2)), label = "Met", angle=90) +
    annotate("text", x = -1.2, y = sum(c(558/2)), label = "Cop", angle=90) +
    # theme(plot.margin = margin(.6,.6,.6,.6, "cm")) +
    theme(plot.margin = margin(.3,.3,.3,.8, "cm")) +
    coord_cartesian(xlim=c(0.5,53.5), clip = "off")
  
}
