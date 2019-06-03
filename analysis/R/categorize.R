bootstrapping <- function(x,n,k,f){
  y <- NULL
  for (idx in 1:k){
    y <- c(y,jitter(sample(x,n,replace = T),factor=f))
  }
  return(y)
}

get_ct <- function(x,n){
  d <- density(x)
  
  x2 <- d$x
  y2 <- d$y
  
  each_area = sum(y2)/n
  tmp <- 0
  idx <- 0
  res <- NULL
  for (y in y2){
    tmp <- tmp+y
    if(tmp > each_area){
      res <- c(res,x2[idx])
      tmp <- 0
    }
    idx <- idx+1
  }
  return(res)
}

A_categorize <- function(x,n_ct,resamp_cnt=10000,noise_factor=10,plot=T){
  bts <- bootstrapping(x,length(x),resamp_cnt,noise_factor)
  ctz <- get_ct(bts,n_ct)
  if(plot){
    hist(bts, breaks=100, probability = T, main = 'Recommanded Categorization')
    lines(density(bts), col='red', lwd=2)
    abline(v=ctz, col='blue')
  }
  return (ctz)
}

