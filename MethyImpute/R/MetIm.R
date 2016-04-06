#' Impute missing methylation data in sequence and microarray to their union dimensions
#'
#' @param sequence is the sequencing-based methylation measurement, usually only for a small number of individuals but almost all CPG sites
#' @param microarray is the array-based methylation measurement, usually only for a small number of CpG sites, but almost all individuals.
#' @param lambda is a parameter controlling for the shrinkage
#' @param use.mvalue is a paramter determine the option if MVN will be calculated in the M-value scale. Default=TRUE
#' @param cvfold is a fold times of cross validation.
#' @param partition is partition.
#' @param cut is cut.
#' @param qc_frac is qc_frac
#' @param overlap is overlap.
#' @return a MetIM object, (not fully implemented yet)
#' @export
#' @examples
#' data(dat1)
#' data(dat2)
#' qc_frac = 1-5/dim(dat2)[2]
#' Seq <- MetIm(sequence = t(dat1), microarray = t(dat2), lambda=0.3, cut=10, cvfold=0, use.mvalue = T, qc_frac = qc_frac)
MetIm <- function(sequence, microarray = NULL ,lambda = NULL, use.mvalue = NULL, cvfold = NULL, partition = NULL, cut = NULL, qc_frac = NULL, overlap = NULL, ...){
  library(MASS)
  # Input data sets sequence & microarray: row and column names must be numeric !!!
  
  # this option determines whether or not to use the M-value scale
  if (is.null(use.mvalue)) {
    use.mvalue = TRUE;
  }
  if(is.null(microarray)){
    res <- seqIM(sequence, cvfold = cvfold , lambda = lambda, partition = partition, cut = cut, overlap = overlap)
    res <- as.matrix(res)
    
  }
  
  if(!is.null(microarray)){
    # (1) create data frames of the union of the dimensions of arraya and seq
    A <- as.matrix(sequence)
    B <- as.matrix(microarray)
    
    if (use.mvalue) {
      C <- Beta_M(A)
      D <- Beta_M(B)
    } else {
      C <- A
      D <- B
    }
    
    
    
    ROW1 <- as.numeric(rownames(C))
    ROW2 <- as.numeric(rownames(D))
    COL1 <- as.numeric(colnames(C))
    COL2 <- as.numeric(colnames(D))
    
    ROW <- sort(union(ROW1,ROW2))
    COL <- sort(union(COL1,COL2))
    nrow <- length(ROW)
    ncol <- length(COL)
    
    C_star <- matrix(NA,nrow,ncol)
    colnames(C_star) <- COL
    rownames(C_star) <- ROW
    D_star <- C_star
    
    id1.row <- match(ROW1,ROW)
    id1.col <- match(COL1,COL)
    id2.row <- match(ROW2,ROW)
    id2.col <- match(COL2,COL)	
    
    C_star[id1.row,id1.col] <- as.matrix(C)
    D_star[id2.row,id2.col] <- as.matrix(D)
    
    # par(mfrow=c(2,2))
    # image(C); image(C_star); image(D); image(D_star)
    
    # (2) Sequence to Array transformation
    # Mapping the sequence and array data by using sample linear model
    beta <- FindMap(X = D_star,Y = C_star)
    
    # (3) Predict transformed array data
    fD <- NULL
    n <- dim(beta)[2]
    for(i in 1:n) {fD <- rbind(fD, beta[2,i]*D_star[,i] + beta[1,i])}
    fD = t(fD)
    
    # (3) Merge fD into C_star
    to_fill = which(is.na(C_star)| !is.na(fD))
    C_star[to_fill] <- fD[to_fill]
    
    if (use.mvalue) {
      A_star = M_Beta(C_star)
    } else {
      A_star = C_star
    }
    
    # (4) transfer array data (predicted) from M value to Beta value
    #   Imputation within A_star
    # has to use the A_star because yet no solution for the Inf problem
    res1 <- seqIM(A_star, cvfold = cvfold , lambda = lambda, partition = partition, cut = cut, qc_frac = qc_frac, overlap = overlap)
    res1 <- as.matrix(res1)
  }
  return(res1)
}

  
  
  ##########################################################################################################
  ####
  #### Imputation for Sequence Data (without training and testing) using CV for --HIGH missing rate---
  ####
  ##########################################################################################################
  
  seqIM <- function(data, cvfold = NULL, lambda = NULL, partition = NULL, cut = NULL, qc_frac = NULL, overlap = NULL,...){
    
    if (is.null(qc_frac)) qc_frac = 0.5
    
    dat1 <- as.matrix(data)
    # qc
    m <- dim(dat1)[2]
    num_na <- apply(is.na(dat1),1,sum)
    qc_out <- which(num_na > qc_frac*m)
    if(length(qc_out) > 0) {
      dat1.new <- dat1[,-qc_out]
    } else {
      dat1.new = dat1
    }
    
    
    # imputation
    Im_seq <- seqIm_low(data = dat1.new, cvfold = cvfold, lambda = lambda, partition = partition, cut = cut, overlap = overlap) 
    
    # renew the sequence data after imputation
    n <- dim(dat1)[2]
    OK_id <- setdiff(1:n,qc_out)
    dat1[,OK_id] <- Im_seq
    #Im_seq <- NULL
    
    
    # impute the high missing rate columns
    # this is only to extrapolate with average values currently
    # may need to be just skipped.
    if (0) {
      
    n <- length(qc_out)
    m <- dim(dat1)[2]
    id <- setdiff(1:m, qc_out)
    
    for(i in 1:n){
      id_NA <- qc_out[i]
      id1 <- id2 <- 0
      
      ok <- which(id < id_NA)
      if(length(ok)>0){
        ok <- max(ok)
        id1 <- id[ok]
      }
      ok <- which(id > id_NA)
      if(length(ok)>0){
        ok <- min(ok)
        id2 <- id[ok]
      }
      
      if(id1>0 & id2>0) A <- 0.5*(dat1[,id1] + dat1[,id2])
      if(id1<=0 & id2>0) A <-  dat1[,id2]
      if(id1>0 & id2<=0) A <- dat1[,id1]
      
      OK <- !complete.cases(dat1[,id_NA]) # missing points for that column	
      dat1[OK,id_NA] <- A[OK]
    }
    
    }
    
    return(dat1)
  }
  
  
  ################################################################################################
  ####
  #### Imputation for Sequence Data (with training and testing)
  ####
  ################################################################################################
  
  
  
  MetIm1 <- function(train, test,lambda = NULL, partition = NULL, cut = NULL, overlap = NULL, ...){
    # partition: a logistic variable, FALSE or True
    # overlap: a logistic variable, FALSE or True
    
    m <- dim(train)[1]
    l <- dim(train)[2] 
    n <- dim(test)[1]
    l1 <- dim(test)[2]
    if(l != l1)stop("reference methylation levels and to-be-filled methylation levels MUST have the same number of columns")
    
    if(is.null(cut))cut <- 100
    if(cut==0)cut <- l
    if(l <= cut)partition <- FALSE
    if(is.null(partition))partition <- TRUE
    if(is.null(overlap))overlap <- TRUE
    
    # small data
    if(!partition)res <- MI(train = train, test = test, lambda = lambda)
    
    # big data (no overlap)
    if(partition & !overlap){   
      num <- floor(l/cut)
      res_mu <- NULL
      res_s2 <- NULL
      for(i in 1:num){
        train_0 <- train[,((i-1)*cut+1):(i*cut)]
        test_0 <- test[,((i-1)*cut+1):(i*cut)]
        RES <- MI(train = train_0, test = test_0, lambda = lambda)
        res_mu <- cbind(res_mu, as.matrix(RES$out.y2))
        res_s2 <- cbind(res_s2, as.matrix(RES$out.s2))
      }
      if(l > num*cut){
        train_0 <- train[,(i*cut+1):l]
        test_0 <- test[,(i*cut+1):l]
        RES <- MI(train = train_0, test = test_0, lambda = lambda)
        res_mu <- cbind(res_mu, as.matrix(RES$out.y2))
        res_s2 <- cbind(res_s2, as.matrix(RES$out.s2))
      }
      res <- list(out.y2 = res_mu, out.s2 = res_s2)
    }
    
    # big data (overlap)
    if(partition & overlap){   
      num <- floor(l/cut)
      res_mu <- NULL
      res_s2 <- NULL
      for(i in 1:num){
        train_0 <- train[,((i-1)*cut+1):(i*cut)]
        test_0 <- test[,((i-1)*cut+1):(i*cut)]
        RES <- MI(train = train_0, test = test_0, lambda = lambda)
        res_mu <- cbind(res_mu, as.matrix(RES$out.y2))
        res_s2 <- cbind(res_s2, as.matrix(RES$out.s2))
      }
      if(l > num*cut){
        train_0 <- train[,(i*cut+1):l]
        test_0 <- test[,(i*cut+1):l]
        RES <- MI(train = train_0, test = test_0, lambda = lambda)
        res_mu <- cbind(res_mu, as.matrix(RES$out.y2))
        res_s2 <- cbind(res_s2, as.matrix(RES$out.s2))
      }
      
      
      over <- floor(cut/2)
      res_mu_o <- NULL
      res_s2_o <- NULL
      
      train_0 <- train[,1:(over-1)]
      test_0 <- test[,1:(over-1)]
      RES <- MI(train = train_0, test = test_0, lambda = lambda)
      res_mu_o <- cbind(res_mu_o, as.matrix(RES$out.y2))
      res_s2_o <- cbind(res_s2_o, as.matrix(RES$out.s2))
      
      
      if(l-over >= cut){
        
        num <- floor((l-over)/cut)
        for(i in 1:num){
          train_0 <- train[,((i-1)*cut+over):(i*cut+over-1)]
          test_0 <- test[,((i-1)*cut+over):(i*cut+over-1)]
          RES <- MI(train = train_0, test = test_0, lambda = lambda)
          res_mu_o <- cbind(res_mu_o, as.matrix(RES$out.y2))
          res_s2_o <- cbind(res_s2_o, as.matrix(RES$out.s2))
        }
        if(l-over > num*cut){
          train_0 <- train[,(i*cut+over):l]
          test_0 <- test[,(i*cut+over):l]
          RES <- MI(train = train_0, test = test_0, lambda = lambda)
          res_mu_o <- cbind(res_mu_o, as.matrix(RES$out.y2))
          res_s2_o <- cbind(res_s2_o, as.matrix(RES$out.s2))
        }
        
      }
      
      if(l-over < cut){    
        train_0 <- train[,over:l]
        test_0 <- test[,over:l]
        RES <- MI(train = train_0, test = test_0, lambda = lambda)
        res_mu_o <- cbind(res_mu_o, as.matrix(RES$out.y2))
        res_s2_o <- cbind(res_s2_o, as.matrix(RES$out.s2))
      }
      
      
      res_mu <- (res_mu + res_mu_o)/2
      res_s2 <- (res_s2 + res_s2_o)/2
      res_mu_o <- NULL
      res_s2_o <- NULL    
      
      res <- list(out.y2 = res_mu, out.s2 = res_s2)
    }
    
    
    return(res)
  }
  

########## GIBBS Sampling

MetIm.Gibbs <- function(train, test,lambda = NULL, cut = NULL, ...){
  # partition: a logistic variable, FALSE or True
  # overlap: a logistic variable, FALSE or True
  # train and test are all pre-filled
  # output: refined test using Gibbs sampling
  
  m <- dim(train)[1] # <- rows are individuals
  l <- dim(train)[2] # <- cols are sites
  n <- dim(test)[1]
  l1 <- dim(test)[2]
  if(l != l1)stop("reference methylation levels and to-be-filled methylation levels MUST have the same number of columns")
  # if(m != n)stop("reference methylation levels and to-be-filled methylation levels MUST have the same number of columns")
  
  if(is.null(cut))cut <- min(100,l)
  if(cut==0)cut <- l
  cut = min(cut,l)
  
  res_mu= test
  eta = 0.05 # the learning rate
  ITER = 5
  Burnin= 0
  trace.res_mu = array(NA, dim=c((ITER-Burnin),m,l))
  res_s2=NULL
  for (iter in (1:ITER))
  {
    for (the_l in sample(1:l))
    # for (the_l in 1:l)
      {
      # res_mu[,the_l] = NA;
      # cut off the window surrounding the_l
      d = cut / 2;
      left  = max(the_l-d,1)
      right = min(the_l+d,l)
      
      RES <- MI.1(train = train[,left:right], test = res_mu[,left:right], u = the_l-left+1, lambda = lambda)
      res_mu[,left:right] = (1-eta) * res_mu[,left:right] + eta * RES$out.y2   
    }
    if (iter > Burnin) {
      trace.res_mu[iter-Burnin,,] =  res_mu
    }
  }
  
  # may find a most efficient way to summing over the 3rd dimension
  # but for now we just use a simple one
  out.y2 = colSums(trace.res_mu, dims = 1) / (ITER-Burnin)
  out.s2 = NULL
  res <- list(out.y2 = out.y2, out.s2 = out.s2)
  return(res)
}



  
  MI <- function(train, test, lambda = NULL, ...){
    # train: m-by-l matrix of reference methylation levels
    # test:  n-by-l matrix of to-be-filled methylation levels
    #        unobserved columns are marked with NAs
    # output: fill up NAs in test via MVN
    # fill up default opt
    if(is.null(lambda))lambda <- 0.3
    
    m <- dim(train)[1]
    l <- dim(train)[2]
    
    n <- dim(test)[1]
    l1 <- dim(test)[2]
    
    o <- which(apply(!is.na(test),2,sum)>0)
    unobserved <- setdiff((1:l),o)
    
    # simplistic non-NA moment estimate for now
    mu <- as.matrix(colMeans(train, na.rm = TRUE, dims = 1)) # 3-31-2015 modified "as.matrix"
    sigma2 <- cov(train,use = "pairwise.complete.obs")
    diag_lambda = max(diag(sigma2), na.rm=T)*lambda # make lambda relative 
    s2 <- sigma2 + diag(diag_lambda,dim(sigma2)[1],dim(sigma2)[2]) # add shrinkage
    
    # NAs <- is.na(diag(s2))
    not.NAs <- which(!is.na(diag(s2)))
    
    old.train <- train
    train <- train[,not.NAs]
    old.test <- test
    test  <- test[,not.NAs]
    old.mu <- mu
    mu <- as.matrix(mu[not.NAs])
    old.s2 <- s2
    s2 <- s2[not.NAs, not.NAs]
    
    # estimate the CpGs for the sample
    y <- as.matrix(test) # intialized as a copy of the input
    out.y2 <- as.matrix(test)
    out.s2 <- as.matrix(test) # zero matrix, size of test
    
    # filling up one sample at a time, this loop can be time-consuming
    for(i in sample(1:n)){# i = 1
      u <- which(is.na(test[i,]))
      o <- which(!is.na(test[i,]))
      
      # the central algorithm
      if(length(o)==0) {
        mu_cond <- mu[u,1]
        s2_cond <- s2[u,u]
      }
      if(length(o)>0){    
        s2_o_inv <- ginv(s2[o,o])
        mu_cond <- mu[u,1] + s2[u,o] %*% s2_o_inv %*% (y[i,o]-mu[o,1])
        s2_cond <- s2[u,u] - s2[u,o] %*% s2_o_inv %*% s2[o,u]
      }
      out.y2[i,u] <- t(mu_cond)
      out.s2[i,u] <- diag(s2_cond)
    }
    
    o.y2 <- array(NA, dim(old.test))
    o.y2[,not.NAs] <- out.y2
    o.s2 <- array(NA, dim(old.test))
    o.s2[,not.NAs] <- out.s2
    res <- list(out.y2 = o.y2, out.s2 = o.s2)
    return(res)
  }
  
  
MI.1 <- function(train, test, u, lambda = NULL, ...){
  # train: m-by-l matrix of reference methylation levels
  # test:  n-by-l matrix of to-be-filled methylation levels
  #        unobserved columns are marked with NAs
  # u     : the unobserved column to be imputed
  # output: fill up NAs in test via MVN
  # just over 1 column
  # fill up default opt
  if(is.null(lambda))lambda <- 0.3
  
  m <- dim(train)[1]
  l <- dim(train)[2]
  
  n <- dim(test)[1]
  l1 <- dim(test)[2]
  if (u > l) {
    1;
  }
  unobserved <- u
  o <- setdiff(1:l, u)
  
  # simplistic non-NA moment estimate for now
  mu <- as.matrix(colMeans(train, na.rm = TRUE, dims = 1)) # 3-31-2015 modified "as.matrix"
  sigma2 <- cov(train,use = "pairwise.complete.obs")
  diag_lambda = max(diag(sigma2), na.rm=T)*lambda # make lambda relative 
  s2 <- sigma2 + diag(diag_lambda,dim(sigma2)[1],dim(sigma2)[2]) # add shrinkage
  
  # NAs <- is.na(diag(s2))
  not.NAs <- which(!is.na(diag(s2)))
  
  old.train <- train
  train <- train[,not.NAs]
  old.test <- test
  test  <- test[,not.NAs]
  old.mu <- mu
  mu <- as.matrix(mu[not.NAs])
  old.s2 <- s2
  s2 <- s2[not.NAs, not.NAs]
  
  if (u %in% not.NAs) {
    u <- match(u,not.NAs)
    o <- setdiff(1:length(not.NAs), u)
  } else {
    res <- list(out.y2 = old.test, out.s2 = old.test) # out.s2 is just for convenience, should change later
    return(res)
    
  }
  
  
  # estimate the CpGs for the sample
  y <- as.matrix(test) # intialized as a copy of the input
  out.y2 <- as.matrix(test)
  out.s2 <- as.matrix(test) # zero matrix, size of test
  
  # filling up one sample at a time, this loop can be time-consuming
  # if ()
  s2_o_inv <- ginv(s2[o,o])
  s2_cond <- s2[u,u] - s2[u,o] %*% s2_o_inv %*% s2[o,u]
  for(i in sample(1:n)){# i = 1
    mu_cond <- mu[u,1] + s2[u,o] %*% s2_o_inv %*% (y[i,o]-mu[o,1])
    out.y2[i,u] <- t(mu_cond)
    out.s2[i,u] <- diag(s2_cond)
  }
  
  o.y2 <- array(NA, dim(old.test))
  o.y2[,not.NAs] <- out.y2
  o.s2 <- array(NA, dim(old.test))
  o.s2[,not.NAs] <- out.s2
  res <- list(out.y2 = o.y2, out.s2 = o.s2)
  return(res)
}


  
  
  
  
  ################################################################################################
  ################################################################################################
  ################################################################################################
  ################################################################################################
  ################################################################################################
  ################################################################################################
  ################################################################################################
  ################################################################################################
  ################################################################################################
  
  
  
  
  
  ##########################################################################################################
  ####
  #### Imputation for Sequence Data (without training and testing) using CV for --low missing rate---
  ####
  ##########################################################################################################
  
  seqIm_low <- function(data, cvfold = NULL, lambda = NULL, partition = NULL, cut = NULL, overlap = NULL){
    
    data <- as.matrix(data)
    
    if(is.null(cvfold))cvfold <- 3
    if(cvfold<=0) 
    {
      test <- train <- data
      RES <- MetIm1(train=train, test=test,lambda = lambda, partition = partition, cut = cut, overlap = overlap)
      if(1) { # Bleeding edge here! Trying to figure out the best way for Gibbs
      res.quick = as.matrix(RES$out.y2)
      cut.d = 100
      RES = MetIm.Gibbs(train, res.quick,lambda = lambda, cut = cut.d) # though this cut is different from the cut in MetIm1
      }
      res = as.matrix(RES$out.y2)
      test <- train <- RES <- NULL
      } else {
      
      N <- dim(data)[1]
      x <- sample(N,N,replace = FALSE) 
      ids <- split(x, factor(sort(rank(x)%%cvfold))) # divided to cvfold groups
      
      res <- NULL
      ID <- NULL
      for (i in 1:cvfold){
        ID <- c(ID,as.vector(unlist(ids[i])))
        test <- data[unlist(ids[i]),]
        train <- data[-unlist(ids[i]),]
        RES <- MetIm1(train=train, test=test,lambda = lambda, partition = partition, cut = cut, overlap = overlap)
        res <- rbind(res,as.matrix(RES$out.y2))
        test <- train <- RES <- NULL
      }
      
      IDs <- data.frame(ID,1:N)
      IDs.new <- IDs[order(ID),2] # sort by ID
      res <- res[IDs.new,] #is.matrix(res)
    }
    
    res[res < 0] <- 0 ### replace all values <0 in a matrix with 0
    res[res > 1] <- 1 ### replace all values >1 in a matrix with 1 
    
    return(res)
  }
  
  
  
  ################################################################################################
  ####
  #### Imputation for Array Data by mappings from Sequence to Array Data
  ####
  ################################################################################################
  
  arrayIm <- function(array, seq, ...){
    
    seq <- as.matrix(seq)
    array <- as.matrix(array)
    
    # (1) Transfer from Beta value to M value
    seq.x <- Beta_M(seq) #is.matrix(seq.x)
    x <- seq.x
    x[x == Inf] <- NA
    x[x == -Inf] <- NA
    
    array.y <- Beta_M(array)
    y <- array.y
    y[y == Inf] <- NA
    y[y == -Inf] <- NA
    
    # (2) Mapping the sequence and array data by using sample linear model
    beta <- RegMap(X = x,Y = y)
    
    # (3) Predict transformed array data
    y.hat <- NULL
    n <- dim(beta)[2]
    for(i in 1:n) y.hat <- cbind(y.hat, beta[1,i]*seq.x[,i] + beta[2,i])
    
    # (4) transfer array data (predicted) from M value to Beta value
    y.hat <- M_Beta(y.hat)
    
    
    # (5) array data imputation
    array.im <- array
    array.im[is.na(array.im)] <- 0 #replace all missing values in a matrix with 0
    
    n <- dim(y.hat)[2]
    for(i in 1:n){
      weight <- complete.cases(array[,i]) + 0
      A <- data.frame(weight, array.im[,i], (1-weight), y.hat[,i])	
      array.im[,i] <- A[,1]*A[,2] + A[,3]*A[,4]	
    }
    
    res <- list(array.im = array.im, array.predicted = y.hat)
    return(res)
  }
  
  
  
  
  
  
  ## Function:  Beta:[0,1] ----> M:(-Inf Inf)
  # need to deal with -Inf and Inf
  Beta_M <- function(data,...){
    data <- as.matrix(data)
    id <- which(is.na(as.vector(data)))
    Range <- range(as.vector(data)[-id])
    if(Range[1]<0 || Range[2]>1)stop("The range of values of input dataset musr be between 0 and 1 !")
    Res <- log(data/(1-data))
    # Res[Res==-Inf] = 
    return(Res)
  }
  
  
  ## Function:  M:(-Inf Inf) ----> Beta:[0,1] 
  M_Beta <- function(data,...){
    data <- as.matrix(data)
    Res <- 1/(1+exp(data)^(-1))
    Res[which(Res==-Inf)]=0
    Res[which(Res== Inf)]=1
    return(Res)
  }
  
  
  ## Function: Simple Regression Transformation
  
  RegMap <- function(X,Y,...){
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    D1 <- dim(X)
    D2 <- dim(Y)
    if(D1[1] != D2[1] || D1[2] != D2[2]) stop("The sizes of two input datasets must be same")
    
    n <- D1[2]
    Beta <- matrix(NA,2,n) # estimated parameters 
    for(i in 1:n){
      x <- X[,i]	
      y <- Y[,i]
      ok<- complete.cases(cbind(x,y))
      if(sum(ok+0)){
        reg <- lm(y~x)
        coeff <- as.vector(reg$coefficients)
        Beta[1,i] <- coeff[1]
        Beta[2,i] <- coeff[2]
      }
    }
    id0 <- which(is.na(Beta[1,]))
    id1 <- which(!is.na(Beta[1,]))
    n0 <- length(id0)
    if(n0 > 0){
      for(i in 1:n0){
        low <- which(id0[i] > id1)
        upp <- which(id0[i] < id1)
        
        L <- U <- matrix(0,2,1)
        if(length(low)) L <- matrix(Beta[,id1[max(low)]],2,1)
        if(length(upp)) U <- matrix(Beta[,id1[min(upp)]],2,1)
        Beta[,id0[i]] <- 0.5*L + 0.5*U
        
        if(!length(low)) Beta[,id0[i]] <- U
        if(!length(upp)) Beta[,id0[i]] <- L
        
        
      }
    }
    return(Beta)
  }
  
  
  
  FindMap <- function(X,Y,...){
    # X, Y are matrices. 
    # Find intersect
    # Find transformation by linear regression for each row
    # In the future, there can be other transformations
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    D1 <- dim(X)
    D2 <- dim(Y)
    if(D1[1] != D2[1] || D1[2] != D2[2]) stop("The sizes of two input datasets must be same")
    
    n <- D1[2]
    Beta <- matrix(NA,2,n) # estimated parameters 
    for(i in 1:n){
      x <- X[,i]	
      y <- Y[,i]
      x[which(x==Inf)] = NA;
      x[which(x==-Inf)] = NA;
      x[which(is.nan(x))] = NA;
      y[which(y==Inf)] = NA;
      y[which(y==-Inf)] = NA;
      y[which(is.nan(y))] = NA;
      
      ok.x<- complete.cases(x)
      ok.y<- complete.cases(y)
      if(sum(ok.y+0) && sum(ok.x+0)){
        df=data.frame(x=x, y=y)
        #     x.na = (is.na(x))
        #      x.ex = x[1-x.na]
        #      y.ex = y[1-x.na]
        if (i==52) 
        {
          1;
        }
        reg <- lm(y~x,df, na.action=na.omit)
        coeff <- as.vector(reg$coefficients)
        Beta[1,i] <- coeff[1]
        Beta[2,i] <- coeff[2]
      }
    }
    
    # fill in the na cases by extrapolation, if needed
    return(Beta)
  }
  
