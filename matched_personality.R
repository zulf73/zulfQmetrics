# indices of matrix in descending order
library(tidyverse)

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


# assume each list has same names 
# put them on rows
list_to_df <- function( L ){
  
}

paired_rows_df <- function( pairs, df ){
  outdf <- data.frame()
  n <- nrow(pairs)
  for ( j in 1:n ){
    x <- pairs$x[j]
    y <- pairs$y[j]
    #print(paste('coords',x,y))
    first <- df[ which( df$type == x), ]
    second <- df[ which( df$type == y), ]
    #print(first)
    #print(second)
    if ( length(first) > 0 && length(second) > 0){
    new_row <- cbind( unlist(first), 
                      unlist(second) ) 
    tryCatch(
      { 
        outdf <- rbind( outdf, new_row )
      }, error = function (e){} )
    }
  }
  # hack just put in the target row
  outdf["target"] <- pairs$z
  outdf
}

# we have pairs_match and ftm
pairs_match <- sorted_indices( countMatrix)
centroid_and_marsat <- paired_rows_df( 
  pairs_match %>% map_df(as_tibble), 
  ftm)
