#' Function for performing latin-grid sampling
#' @param x numeric vector
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
sample_between_intervals <- function(x){
  intervals <- sort(unique(x))
  if(length(intervals) >= length(x)/2){
    x
  } else {
    sapply(x,function(y){
      ## if on the max, generate a uniform between a randomly-selected interval
      if (y == max(intervals)){
        s <- sample(1:(length(intervals)-1),1)
        runif(1,intervals[s],intervals[s+1])
        ## otherwise, generate uniform between obs and next interval
      } else {
        runif(1,y,min(intervals[intervals > y]))
      }
    })
  }
}


#' Function for centering and scaling the transformed metric
#' @param x numeric vector
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
center_and_scale <- function(x){
  xx <- na.omit(x)
  xx <- xx[is.finite(xx)]
  xx <- xx[!is.nan(xx)]
  if(length(xx) > 1){
    meann <- mean(xx[is.finite(xx)],na.rm = T)
    sdd <- sd(xx[is.finite(xx)],na.rm = T)
    if(sdd > 0){
      return(list((x-meann)/sdd,
                  meann,
                  sdd
      ))
    } else {
      return(list((x-meann),
                  meann,
                  1))
    }
  } else {
    return(list(x,
                0,
                1))
  }
}

#' Function useful for turning cubic expansions into knot expansions
#' @param mat matrix
#' @param Q matrix
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
#' @export
knot_expand2 <- function(mat,Q){
  mat[,rep(1:ncol(mat),ncol(cbind(Q)))] * 
    cbind(Q)[,rep(1:ncol(cbind(Q)),each = ncol(mat))]
}

#' Variance stabilizing/symmetry-inducing transformations
#' @param x numeric vector
#' @param type character
#' @param pwr_lmbda optional numeric
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
#' @export
transf <- function(x,type,pwr_lmbda = NULL,...){
  #pwr_lmbda <- 0
  if(type == 'distance'){
    ## use box-cox
    if(pwr_lmbda != 0){
      (x^pwr_lmbda+1)/pwr_lmbda
      # x^pwr_lmbda
    } else {
      log(x)
    }
  } else if(type == 'proportion'){
    ## use arcsin
    asin(sqrt(x))
  } else if(type == 'interval'){
    ## use hyperbolic-tangent transformation
    x <- ifelse(x <= -1,-0.99999,x)
    x <- ifelse(x >= 1,0.99999,x)
    x <- ifelse(is.nan(x),0,x)
    atanh(x)
  } else {
    return(...)
  }
}

#' variance stablizing/symmetry-inducing back-transformations
#' @param x numeric vector
#' @param type character
#' @param pwr_lmbda optional numeric
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
#' @export
un_transf <- function(x,type,pwr_lmbda = NULL,...){
  if(type == 'distance'){
    if(pwr_lmbda != 0){
      (x * pwr_lmbda - 1)^(1 / pwr_lmbda)
    } else {
      exp(x)
    }
  } else if(type == 'proportion'){
    (sin(x))^2
  } else if(type == 'interval'){
    tanh(x)
  } else {
    return(...)
  }
}

#' Assign a knot group (for quantiles q) to a vector x
#' @param x numeric vector
#' @param q numeric vector
#' @param rev logical
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
#' @export
get_group <- function(x,q,rev = F){
  if(length(q) > 1){
    rowsums <- rowSums(sapply(2:length(q),function(qq){
      qq*(x <= q[qq] &
            x > q[qq-1])
    }))
  } else {
    rowsums <- 0
  }
  rowsums + 
    ifelse(x <= q[1],
           1,
           0) +
    ifelse(x > q[length(q)],
           length(q)+1,
           0)
}

#' Generate inverse-gamma variables, with error handling
#' @param n numeric
#' @param a numeric
#' @param b numeric
#' @param s numeric
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
#' @export
my_rinvgamma <- function(n, a, b, s){
  try <- rgamma(n, 
                a, 
                rate = b)
  t <- try({
    if(is.na(try)){
      return(s)
    } else {
      if(is.nan(try)){
        return(s)
      } else {
        if(!is.finite(try) | try < 1e-16){
          return(s)
        } else {
          return(1/try)
        }
      }
    }},silent = T)
  if(class(t) == 'try-error'){
    return(s)
  }
}


#' Fix NAs
#' @param x numeric
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
#' @export
fixna <- function(x){
  if(any(na.omit(c(is.na(x),
                   is.nan(x),
                   !is.finite(x))))){
    rep(0,length(x))
  } else {
    x
  }
}

#' First attempting a normal inversion, before resorting to the M.P.G.I.
#' @param x matrix
#' @param tol numeric
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
#' @export
solvetry <- function(x, 
                     tol = sqrt(.Machine$double.eps)){
  t <- try(solve(x),silent = T)
  if(any(class(t) == 'try-error')){
    svd <- svd(x)
    d_inv <- ifelse(svd$d <= tol,
                    0,
                    1/svd$d)
    (svd$v) %*% 
      (d_inv * 
         t(svd$u))
  } else {
    t
  }
}

## Apply is.na, is.nan operation to all elements in a vector
#' @param x numeric vector
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
#' @export
is.nanf <- function(x){
  sapply(x,function(xx)!is.finite(xx) | is.nan(xx) | is.na(xx))
}

# Take the l2-norm
#' @param x numeric vector
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
#' @export
l2norm2 <- function(x,other_vec = 0)sqrt(mean(x-other_vec)^2)

#' Take a first derivative, with respect to a variable
#' @param dat matrix 
#' @param var character
#' @param stddevs numeric vector
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
#' @export
first_deriv <- function(dat,var,stdevs){
  cols <- grep(paste0(var,"_degree="),colnames(dat))
  which_sing <- which(colnames(dat) == paste0(var,"_degree=1"))
  cols <- cols[cols != which_sing]
  dat_deriv <- dat
  dat_deriv[,-cols] <- 0
  substrings <- sapply(colnames(dat)[c(cols,which_sing)],
                       function(cn)as.numeric(substr(cn,nchar(cn),nchar(cn))))
  names(substrings) <- colnames(dat)[c(cols,which_sing)]
  denom <- (dat[,which_sing]*stdevs[which_sing])
  tmp <- sapply(c(cols,which_sing),function(j){
    colname <- colnames(dat_deriv)[j]
    ## if interaction, simply divide by x = axy / x = ay
    if(grepl("xXx",colname)){
      ## = observed value divided by the variable 
      dat[,j]/denom
      ## if polynomial term, multiply by degree, divide by variable = ax^(a-1) = a * (x^a / x)
    } else {
      dgr <- substrings[colname]
      dgr*dat[,j]/denom
    }
  })
  dat_deriv[,c(cols,which_sing)] <- tmp
  dat_deriv[is.nan(dat_deriv)] <- 0
  dat_deriv[is.na(dat_deriv)] <- 0
  dat_deriv[!is.finite(dat_deriv)] <- 0
  dat_deriv
}

#' Take a second derivative for only interactions
#' @param dat matrix 
#' @param var character
#' @param stddevs numeric vector
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
#' @export
interact2_deriv <- function(dat,var,stdevs){
  cols <- grep(paste0(var,"_degree="),colnames(dat))
  which_sing <- which(colnames(dat) == paste0(var,"_degree=1"))
  cols <- cols[cols != which_sing]
  dat_deriv <- dat
  dat_deriv[,-cols] <- 0
  substrings <- sapply(colnames(dat)[c(cols,which_sing)],
                       function(cn)as.numeric(substr(cn,nchar(cn),nchar(cn))))
  names(substrings) <- colnames(dat)[c(cols,which_sing)]
  cols <- cols[grepl("xXx",cols)]
  tmp <- sapply(c(cols,which_sing),function(j){
    colname <- colnames(dat_deriv)[j]
    if(grepl('xxYYzz',colname)){
      catvar <- unlist(strsplit(colname,'xxYYzz'))[1]
      dat[,catvar] / stdevs[j]
    } else {
      rep(1 / stdevs[j],nrow(dat_deriv))
    }
  })
  dat_deriv[,c(cols,which_sing)] <- tmp
  dat_deriv[is.nan(dat_deriv)] <- 0
  dat_deriv[is.na(dat_deriv)] <- 0
  dat_deriv[!is.finite(dat_deriv)] <- 0
  dat_deriv
}

#' Compute second derivatives, no interactions
#' @param dat matrix 
#' @param var character
#' @param stddevs numeric vector
#'
#' @details for internal use with set_optim_thompson function
#'
#' @return Object of arbitrary class.
#'
#' @examples
#' \donttest{
#'   # No Examples Shown
#' }
#' @export
second_deriv <- function(dat,var,stdevs){
  cols <- grep(paste0(var,"_degree="),colnames(dat))
  which_sing <- which(colnames(dat) == paste0(var,"_degree=1"))
  cols <- cols[cols != which_sing]
  dat_deriv <- dat
  dat_deriv[,-cols] <- 0
  cols <- cols[!grepl("xXx",cols)]
  substrings <- sapply(colnames(dat)[c(cols,which_sing)],
                       function(cn)as.numeric(substr(cn,nchar(cn),nchar(cn))))
  names(substrings) <- colnames(dat)[c(cols,which_sing)]
  denom <- (dat[,which_sing]*stdevs[which_sing])^2
  tmp <- sapply(c(cols,which_sing),function(j){
    dgr <- substrings[colnames(dat_deriv)[j]]
    if(dgr <= 1){
      rep(0,nrow(dat))
    } else {
      dgr*(dgr-1)*dat[,j]/denom
    }
  })
  dat_deriv[,c(cols,which_sing)] <- tmp
  dat_deriv[is.nan(dat_deriv)] <- 0
  dat_deriv[is.na(dat_deriv)] <- 0
  dat_deriv[!is.finite(dat_deriv)] <- 0
  dat_deriv
}