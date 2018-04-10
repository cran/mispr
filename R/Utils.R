
sym <- function (x)
{  (x + t(x))/2 }

class1 <- function(x) class(x)[1]
#################################

var.select = function(xobs, yobs, eps = 1e-04, maxcor = 0.99 )
{
  sel.tune.func=function(keep,xobs,rmv,eps ) {
       k <- sum(keep)
       cx <- cor(xobs[, keep, drop = FALSE], use = "all.obs")
       eig <- eigen(cx, symmetric = TRUE)
       ncx <- cx
       if(!rmv){
    while (eig$values[k]/eig$values[1] < eps) {
    j <- (1:k)[order(abs(eig$vectors[, k]), decreasing = TRUE)[1]]
    keep[keep][j] <- FALSE    #Keep[j] <- FALSE
    ncx <- cx[keep, keep, drop = FALSE]
    k <- k - 1
    eig <- eigen(ncx)
                               }
                }else{
                while (eig$values[k]/eig$values[1] < eps) {
                j <- (1:k)[order(abs(eig$vectors[, k]), decreasing = TRUE)[1]]

                keep[keep][j] <- FALSE
                ncx <- cor(xobs[, keep, drop = FALSE], use = "all.obs")
                k <- k - 1
                eig <- eigen(ncx)
                                                       }
                     }
    return(keep)
    }
  keep <- unlist(apply(xobs, 2, var) > eps)
  keep[is.na(keep)] <- FALSE
  highcor <- suppressWarnings((unlist(apply(xobs, 2, cor, yobs)) < maxcor))
  keep <- keep & highcor
  if(any(!keep)) rmv=TRUE else rmv=FALSE
  if (all(!keep))  {print("All predictors are constant or have too high correlation. Only choice is intercept model"); return(NULL) }
  if (length(keep) == 1)  keep[1] <- TRUE


res.sel=sel.tune.func(keep,xobs,rmv,eps=eps)
if(dim(xobs)[1]<sum(res.sel)+3){
while(dim(xobs)[1]<sum(res.sel)+3){
eps=ifelse(eps>=0.1,eps+.1,sqrt(eps))
res.sel=sel.tune.func(keep,xobs,rmv,eps=eps)
}
}
return(res.sel)
}
########################################################################################################

### utility functions
na.exist <- function(X) any(is.na(X))

random.impute <- function(x.missing){
for(i in 1:dim(x.missing)[2]) {
which.na <- which(is.na(x.missing[,i]))
if(length(which.na)>0) x.missing[which.na,i] <- sample(x.missing[-which.na,i],length(which.na),replace=TRUE)
 }
 return(x.missing)
 }
##############################################
initialise <- function(x,method){
  if(method=="median"){
    for( j in 1:ncol(x) ) {
      xx <- x[,j]
      if(class(xx) == "numeric") {xx <- as.vector(impute(as.matrix(xx), "mean"))}
      if(class(xx) == "integer") {xx <- as.integer(impute(as.matrix(xx), "median"))}
      if(class(xx) == "factor")  {xx <- as.character(xx)
        xx[which(is.na(xx))] <-  names(which.max(table(xx)))
        xx <- as.factor(xx)}
      x[,j] <- xx
    }
  }else   if(method=="random"){x <- random.impute(x)} #else x <- invisible(kNN(x,imp_var=FALSE))

  return(x)
}
###################################
## switch function to automatically select methods

lm.glm <- function(xReg, x.select, pen, index,type, form, L2.fix, cv, maxL2=maxL2, track ) {
   switch(type,
      numeric = useLM(xReg, x.select, pen, index, L2.fix, cv, maxL2=maxL2,track=track),
      bin     = useBin(xReg, x.select, pen, index, L2.fix, cv, maxL2=maxL2,track=track)
              )
}
################################################################################
## COnverting data.frame into matrix
Data.matrix <- function(dataframe){
                            if(is.matrix(dataframe)){return(dataframe)}else
                            if(is.data.frame(dataframe)){
                            bin.col <- which(sapply(dataframe,class)=="factor")
                            x.matrix <- data.matrix(dataframe)
                            if(length(bin.col)>0) x.matrix[,bin.col] <- x.matrix[,bin.col]-1
                            return(x.matrix)
                                                        }else {print("input object must be matrix or data.frame"); return(0)}
                                  }  # End of xPart.standardize Function


################################################################################


Fitting_glm <- function(xyData, ndata, wy, form ){
  options(warn = 2)
  glmfit <- glm(form, data=xyData, family="gaussian")
  options(warn = 1)
  if(   any(!is.finite(c(summary(glmfit)$coefficients[,1:2])))    ) stop("NAs or NaNs produced during normnal glm fit")

    df.resid <- max( glmfit$df.residual , 1 )
    Mu <- glmfit$coefficients

      sigma.new <- sqrt(  sum(glmfit$residuals^2)/rchisq(1,df.resid)  )
      Fish.inv = summary(glmfit)$cov.unscaled
      beta.hat.new <- Mu + (t(chol(sym(Fish.inv))) %*% rnorm(ncol(Fish.inv))) * sigma.new           #rnorm(ncol(x.intercept))
      eta.new <- cbind(1,Data.matrix(ndata)) %*% beta.hat.new + rnorm(dim(ndata)[1]) * sigma.new
      eta.new <- as.vector(eta.new)
      return( eta.new  )

  } #end of Fitting_glm function

  #############################################################################################################

  Fitting_ridge <- function(xyData, ndata, wy, form, Lambda2 ){
     pen <- penalized(form, data=xyData, trace=FALSE, standardize=TRUE, epsilon=0.0001 , lambda2=Lambda2, model='linear' , maxiter=50)

      df <- max( (length(pen@fitted) - dim(xyData)[2]) , 1 )
      Mu <- c(pen@unpenalized, pen@penalized)
      sigma.new <- sqrt(  sum(pen@residuals^2)/rchisq(1,df)  )
      x.intercept <- as.matrix( cbind(1, Data.matrix(xyData[,colnames(xyData)!="y.imp",drop=FALSE]))  )
      xwx <-  crossprod(x.intercept)
      Fish <- xwx + Lambda2*diag(  c(0, rep(1,dim(xyData[,colnames(xyData)!="y.imp",drop=FALSE])[2]) )   )
      Fish.inv <- solve(    Fish     )
      cov.mat.unscaled =  Fish.inv %*% crossprod(x.intercept) %*% Fish.inv
      beta.hat.new <- Mu + (t(chol(sym(cov.mat.unscaled))) %*% rnorm(ncol(cov.mat.unscaled))) * sigma.new
      eta.new <- cbind(1,Data.matrix(ndata)) %*% beta.hat.new + rnorm(dim(ndata)[1]) * sigma.new
      return( as.vector(eta.new)  )
  } #end of Fitting_ridge function
  #################################################################################

useLM <- function(xReg,  x.select, pen, wy,  L2.fix, cv, maxL2,track ){

ypart <- xReg[, which(names(xReg) == "y.imp"),drop=FALSE]
xpart <- xReg[, which(names(xReg) != "y.imp"),drop=FALSE]

    if( x.select || (dim(xReg)[1]-dim(xReg)[2]-length(wy))<=0 ){
    select <- var.select(yobs=as.numeric(ypart[-wy,1]),
            xobs=Data.matrix( xpart[-wy, , drop=FALSE]  )  )
		xyData <- data.frame(ypart[-wy, ,drop=FALSE], xpart[-wy,select,drop=FALSE])
		ndata<-data.frame(xpart[wy,select,drop=FALSE])
                   } else
                   {
             xyData <- data.frame(ypart[-wy, ,drop=FALSE], xpart[-wy, , drop=FALSE])
             ndata<-data.frame(xpart[wy,,drop=FALSE])
                   }

    form <- names(xyData)[names(xyData)!="y.imp"]
    if(length(form)>0)  { form <- as.formula(paste("y.imp ~",paste(form,collapse="+")))  } else {form <- y~.}

           #if(track){ print(paste("formula used:",form)) }
           if(track){   print(form)  }


if(!pen & dim(xyData)[1]>dim(xyData)[2]){

Eta.new <- try( Fitting_glm(xyData=xyData, ndata=ndata, wy, form=form) , silent=TRUE )
        if(inherits(Eta.new,"try-error")==TRUE){

CV <- try( optL2(form,data=xyData, trace=FALSE, model='linear',fold=5, standardize=TRUE, epsilon=0.0001,minlambda2=0.01,maxlambda2=maxL2,maxiter=50) , silent=TRUE )#maxlambda2=1 or 100 or any
        if(inherits(CV,"try-error")==TRUE) {CV <- optL2(form,data=xyData, trace=FALSE,model='linear',fold=dim(xyData)[1], standardize=TRUE, epsilon=0.0001,minlambda2=0.01,maxlambda2=maxL2,maxiter=50)}
        Lambda2 <- CV$fullfit@lambda2
        Eta.new <- Fitting_ridge(xyData=xyData,ndata=ndata, wy=wy, form=form, Lambda2=Lambda2)
                                                }
       }else{
 if(cv){
  CV <- try( optL2(form,data=xyData, trace=FALSE,model='linear',fold=5, standardize=TRUE, epsilon=0.0001,minlambda2=0.01,maxlambda2=maxL2,maxiter=50) ,silent=TRUE )   #maxlambda2=1 or 100 or any
  if(inherits(CV,"try-error")==TRUE) {CV <- optL2(form,data=xyData, trace=FALSE,model='linear',fold=dim(xyData)[1], standardize=TRUE, epsilon=0.0001,minlambda2=0.01,maxlambda2=maxL2,maxiter=50) }
  Lambda2 <- CV$fullfit@lambda2
  Eta.new <- Fitting_ridge(xyData=xyData,ndata=ndata, wy=wy, form=form, Lambda2=Lambda2)
         }else{
         Lambda2 <- L2.fix
         Eta.new <- Fitting_ridge(xyData=xyData,ndata=ndata, wy=wy, form=form, Lambda2=Lambda2)
               }
          }
return(Eta.new)

}
#################################################################
  Fitting_glm.bin <- function(xyData, ndata, wy, form ){

    options(warn = 2)
    glmfit <- glm(form, data=xyData, family="binomial")
    options(warn = 1)
    if(   any(!is.finite(c(summary(glmfit)$coefficients[,1:2])))    ) stop("NAs or NaNs produced during logistic glm fit")
    ###
    se.check <- which( sqrt(diag(summary(glmfit)$cov.unscaled))  > 50 )
    if(length(se.check) >= length(glmfit$coefficients)/2) stop("Infinitely large SEs are produced")

      sigma.beta.hat <- t(chol(sym(summary(glmfit)$cov.unscaled)))                           # OR t(chol(sym(Fish.inv)))
      beta.hat.new <- glmfit$coefficients + sigma.beta.hat %*% rnorm(ncol(sigma.beta.hat))
      eta.new <- as.vector( cbind(1,Data.matrix(ndata)) %*% beta.hat.new )
      fit.new <- 1/(1+exp(-eta.new))
      vec <- (runif(length(fit.new)) <= fit.new)
      vec[vec] <- 1
      return(vec)
  } #end of Fitting_glm.bin function
################################################################################
  Fitting_bin.ridge <- function(xyData, ndata, wy, form, Lambda2 ){
    pen <- penalized(form, data=xyData, trace=FALSE, standardize=TRUE, epsilon=0.0001, lambda2=Lambda2, model='logistic', maxiter=50)
    Mu = c(pen@unpenalized, pen@penalized)

      W <- diag( pen@fitted * (1-pen@fitted) )

      x.intercept <-  cbind(1, Data.matrix(xyData[,colnames(xyData)!="y.imp",drop=FALSE]))
      xwx <- t(x.intercept) %*% W %*% x.intercept
      Fish <- xwx + Lambda2*diag(  c(0, rep(1,dim(xyData[,colnames(xyData)!="y.imp",drop=FALSE])[2]) )   )

      Fish.inv <- solve(    Fish     )
      cov.beta <- Fish.inv %*% xwx %*% Fish.inv
      sigma.beta.hat <-   t(chol(sym(cov.beta)))

      beta.hat.new <- Mu + sigma.beta.hat %*% rnorm(ncol(sigma.beta.hat))
      eta.new <- as.vector( cbind(1,Data.matrix(ndata)) %*% beta.hat.new )

      fit.new <- 1/(1+exp(-eta.new))
      vec <- (runif(length(fit.new)) <= fit.new)
      vec[vec] <- 1
      return(vec)

  } #end of Fitting_ridge function


#################################################################################

useBin <- function(xReg, x.select, pen,  wy,  L2.fix, cv, maxL2, track ){

ypart <- xReg[, which(names(xReg) == "y.imp"),drop=FALSE]
xpart <- xReg[, which(names(xReg) != "y.imp"),drop=FALSE]

    if( x.select || (dim(xReg)[1]-dim(xReg)[2]-length(wy))<=0 ){
    select <- var.select(yobs=as.numeric(ypart[-wy,1]),
            xobs=Data.matrix( xpart[-wy, , drop=FALSE]  )  )

		xyData <- data.frame(ypart[-wy, ,drop=FALSE], xpart[-wy,select,drop=FALSE])
		ndata<-data.frame(xpart[wy,select,drop=FALSE])
                   } else
                   {
             xyData <- data.frame(ypart[-wy, ,drop=FALSE], xpart[-wy, , drop=FALSE])
             ndata<-data.frame(xpart[wy,,drop=FALSE])
                   }

    form <- names(xyData)[names(xyData)!="y.imp"]
    if(length(form)>0)  { form <- as.formula(paste("y.imp ~",paste(form,collapse="+")))  } else {form <- y~.}

              if(track){   print(form)  }


if(!pen & dim(xyData)[1]>dim(xyData)[2]){

Eta.new <- try( Fitting_glm.bin(xyData=xyData, ndata=ndata, wy=wy, form=form) , silent=TRUE )
        if(inherits(Eta.new,"try-error")==TRUE){

        CV <- try( optL2(form,data=xyData, trace=FALSE,model='logistic',fold=5,standardize=TRUE, epsilon=0.0001,minlambda2=0.01,maxlambda2=maxL2,maxiter=50), silent=TRUE )
        if(inherits(CV,"try-error")==TRUE) {CV <- optL2(form,data=xyData, trace=FALSE,model='logistic',fold=dim(xyData)[1], standardize=TRUE, epsilon=0.0001,minlambda2=0.01,maxlambda2=maxL2,maxiter=50)}
        Lambda2 <- CV$fullfit@lambda2

        Eta.new <- Fitting_bin.ridge(xyData=xyData,ndata=ndata, wy=wy, form=form, Lambda2=Lambda2)

                                                }
           } else{
  if(cv){
  CV <- try( optL2(form,data=xyData, trace=FALSE,model='logistic',fold=5, standardize=TRUE, epsilon=0.0001,minlambda2=0.01,maxlambda2=maxL2,maxiter=50), silent=TRUE )
  if(inherits(CV,"try-error")==TRUE) {CV <- optL2(form,data=xyData, trace=FALSE,model='logistic',fold=dim(xyData)[1], standardize=TRUE, epsilon=0.0001,minlambda2=0.01,maxlambda2=maxL2,maxiter=50)}
  Lambda2 <- CV$fullfit@lambda2
  Eta.new <- Fitting_bin.ridge(xyData=xyData,ndata=ndata, wy=wy, form=form, Lambda2=Lambda2)
          }else{
           Lambda2 <- L2.fix
           Eta.new <- Fitting_bin.ridge(xyData=xyData,ndata=ndata, wy=wy, form=form, Lambda2=Lambda2)
               }
                     }
   return(factor( Eta.new, c(0, 1), levels(unlist(ypart))  ))
}
################################################################
