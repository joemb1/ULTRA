
#' ULTRA
#' 
#' Fits prediction model for all features and selected features.
#' 
#' @param data dataframe containing z-scaled continuous features and single binary outcome
#' @param namebin name of the binary outcome variable
#' @return Displays a combined cowplot. Note the individual plots can be accessed from the object returned by the ULTRA function.
#' @export
ULTRA <- function(data, namebin){
  
  #sum(sim_data[,selectbinary]==1)/N
  
  #sim_data[,namebin] <- as.factor(sim_data[,namebin])
  #ggduo(sim_data,   mapping = aes_string(color = namebin))
  
  ### univariate  regression per var
  
  sim_data <- data  
  
  zvars <-  colnames(sim_data)[colnames(sim_data)!=namebin]
  
  bvar <- as.factor(sim_data[,namebin])
  
  unistore <- data.frame(row.names = zvars, var=zvars, 
                         beta=rep(NA, length(zvars)),
                         pval=rep(NA, length(zvars)) )
  
  
  for(v in zvars){
    
    form <- as.formula(paste(namebin, "~", v))
    fit <- glm(data=sim_data, form, family="binomial")
    
    
    unistore[v, "beta"] <- fit$coefficients[2]
    unistore[v, "pval"] <- summary(fit)$coefficients[2,4]
    
  }
  
  
  head(unistore[order(unistore$pval),])
  
  
  order <- unistore[order(unistore$pval),]
  
  order$AUC <- NA
  
  order$iter_p <- NA
  
  order$toadd <- NA
  
  order$order <- 1:length(zvars)
  
  order$fracnewinfo <- NA
  
  platt0 <- glm(bvar ~ 1, family="binomial" )
  pred0 <- predict(platt0, type="response")
  
  var0 <- var(pred0)
  
  loglik0 <- sum(sim_data[,namebin] * log(pred0)  + (1- sim_data[,namebin])*log(1-pred0) )
  loglik0
  
  for(i in 1:length(zvars)){
    
    score <- as.matrix(sim_data[,order$var[1:i]]) %*% order$beta[1:i]
    
    order$AUC[i] <- pROC::auc(pROC::roc(bvar ~ as.numeric(score)))
    
    ## PLATT SCALING
    platt2 <- glm(bvar ~ score, family="binomial" )
    pred2 <- predict(platt2, type="response")
    
    var1 <- var(pred2)
    order$fracnewinfo[i] <- 1 - var0/var1
    
    
    loglik1 <- sum(sim_data[,namebin] * log(pred2)  + (1- sim_data[,namebin])*log(1-pred2) )
    loglik1
    
    LRT <- -2*(loglik0 - loglik1)
    LRT
    
    pvalLRT <- 1 - pchisq(LRT, df=1)
    pvalLRT
    order$iter_p[i] <-  pvalLRT
    
    ### add by LRT
    if(pvalLRT < 0.05){
      
      loglik0 <- loglik1
      order$toadd[i] <- order$var[i]
    }
    
    ### add by fraction new info (add if greater than 5%)
    if(1 - var0/var1 > 0.05){
      
      var0 <- var1
      order$toadd[i] <- order$var[i]
    }
  }
  
  
  ### calculate AUC based on FEATURE SELECTION USING LRT iteration
  select <- unique(na.omit(order$toadd ))
  
  selectdf <- order[order$var %in% select, c("var","beta")]
  
  selectdf$AUC <- NA
  
  for(s in 1:length(select)){
    
    score <- as.matrix(sim_data[,selectdf$var[1:s]]) %*% selectdf$beta[1:s]
    
    selectdf$AUC[s] <- pROC::auc(pROC::roc(bvar ~ as.numeric(score)))
    
  }
  
  
  ### PLATT SCALED FINAL MODELS
  num_signalfeats <- which.max(order$AUC)
  
  score <- as.matrix(sim_data[,order$var[1:num_signalfeats]]) %*% order$beta[1:num_signalfeats]
  
  platt_signal <- glm(bvar ~ score, family="binomial" )
  
  
  score <- as.matrix(sim_data[,selectdf$var]) %*% selectdf$beta
  
  platt_select <- glm(bvar ~ score, family="binomial" )
  
  
  
  order$varf <- as.factor(order$var)
  order$varf <- factor(order$varf, levels=order$var )
  
  
  selectdf$varf <- as.factor(selectdf$var)
  selectdf$varf <- factor(selectdf$varf, levels=selectdf$var )
  
  plot0 <- ggplot2::ggplot(order, ggplot2::aes(x=varf, y=AUC))+
    ggplot2::geom_point()+
    
    ggplot2::geom_line( size=1, group=1, alpha=0.2) +
    
    ggplot2::geom_point(data=selectdf, col="red", 
                        position = ggplot2::position_nudge(x = 0.4) )+
    ggplot2::geom_line(data=selectdf, color='red', size=1, group=1, alpha=0.2,
                       position = ggplot2::position_nudge(x = 0.4)) +
    ggplot2::theme_bw()+
    ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                   axis.title.x = ggplot2::element_blank())
  
  
  end <- which(order$var == tail(select,1) ) + 4
  
  plot2df <- order[1:end , ]
  
  plot <- ggplot2::ggplot(plot2df, ggplot2::aes(x=varf, y=AUC))+
    ggplot2::geom_point(position = ggplot2::position_nudge(x = -0.1))+
    
    ggplot2::geom_line( size=1, group=1, alpha=0.2,
                        position = ggplot2::position_nudge(x = -0.1)) +
    
    ggplot2::geom_point(data=selectdf, col="red", 
                        position = ggplot2::position_nudge(x = 0.1) )+
    ggplot2::geom_line(data=selectdf, color='red', size=1, group=1, alpha=0.2,
                       position = ggplot2::position_nudge(x = 0.1)) +
    ggplot2::theme_bw()+
    ggplot2::xlab("Feature")
  
  
  return(list(allfeats=order, 
              selectfeats=selectdf,
              platt_select=platt_select,
              platt_signal=platt_signal,
              allplot=plot0,
              selectplot=plot))
  
}


#' ULTRAplot
#' 
#' Plots the AUC from the prediction model for all features and selected features.
#' @param results The object returned by the ULTRA function.
#' @return Visual output displaying a combined cowplot. Note the individual plots can be accessed from the object returned by the ULTRA function.
#' @export
ULTRAplot <- function(results){
  
  num_signalfeats <- which.max(results$allfeats$AUC)
  num_selectfeat <- nrow(results$selectfeats)
  
  plot0 <- results$allplot +
    ggplot2::ggtitle(paste("Number of features with signal:", num_signalfeats,
                           "<br><span style = 'color:red;'>Number of features after selection"  
                           , num_selectfeat, "</span>"  ))+
    ggplot2::theme(plot.title = ggtext::element_markdown())
  
  print( cowplot::plot_grid(plot0, results$selectplot, ncol=1) )
  
}





