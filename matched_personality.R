# indices of matrix in descending order
library(tidyverse)
library(dplyr)

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
pairs_match <- sorted_indices( countMatrix)
centroid_and_marsat <- paired_rows_df( 
  pairs_match %>% map_df(as_tibble), 
  ftm)