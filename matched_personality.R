# indices of matrix in descending order
library(tidyverse)
library(dplyr)
library(glmnet)

sorted_indices<-function( mtx ){
  v <- as.vector(mtx)
  w <- sort( v, decreasing=TRUE)  
  out <- list()
  for (k in 1:length(w)){
    if (w[k]> 1){
    idx <- which(mtx == w[k],arr.ind = T)
    out[[k]]<-list(x=idx[1],y=idx[2],z=w[k])
    }
  }
  out
}


paired_rows_df <- function( pairs, df ){
  outdf <- NULL
  n <- nrow(pairs)
  for ( j in 1:n ){
    x <- pairs$x[j]
    y <- pairs$y[j]
    z <- pairs$z[j]
    #print(paste('coords',x,y))
    first <- df[ df$type == x, ]
    second <- df[ df$type == y, ]
    #print(first)
    #print(second)
    if ( length(first) > 0 && length(second) > 0){
      # hack just put in the target row
        new_row <- bind_rows( unlist(first) - unlist(second) )
        names(new_row)<-names(df)
        new_row["target"] <- z
        #print(typeof(new_row))
        #print(length(new_row))
        print(dim(outdf))
        if (is.null(outdf)){
          outdf<-as_tibble(new_row)
        } else {
            outdf <- outdf %>% add_row( new_row )
        }
    }
  }
  print(dim(outdf))
  outdf
}

# we have pairs_match and ftm
pairs_match <- sorted_indices( data.matrix(countMatrix))
centroid_and_marsat <- paired_rows_df( 
  pairs_match %>% map_df(as_tibble), 
  ftm)


#########################################
# Fit linear models and penalized ones
tcm <- exp(-abs(centroid_and_marsat))
# Fit models here on tcm

tcm_predictors <- data.matrix(tcm[,!names(tcm) %in% c("type","target")])
tcm_response <- data.matrix(tcm[,"target"])

X <- tcm_predictors
y <- tcm_response

glm_mod<-cv.glmnet(scale(X), y, family="gaussian", standardize=F, type.measure = "mse",
                   nfolds=5,alpha=.5)
c<-coef(glm_mod,s='lambda.min',exact=TRUE)
inds<-which(c!=0)
variables<-row.names(c)[inds]
variables<-variables[variables %in% '(Intercept)']
#rsq = 1 - glm_mod$cvm/var(y)
#plot(glm_mod$lambda,rsq)
varImp <- function(object, lambda = NULL, ...) {
  
  ## skipping a few lines
  
  beta <- predict(object, s = lambda, type = "coef")
  if(is.list(beta)) {
    out <- do.call("cbind", lapply(beta, function(x) x[,1]))
    out <- as.data.frame(out, stringsAsFactors = TRUE)
  } else out <- data.frame(Overall = beta[,1])
  out <- abs(out[rownames(out) != "(Intercept)",,drop = FALSE])
  out
}
