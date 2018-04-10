#'Multiple Imputation with Sequential Penalized Regression 
#'
#'Generates Multivariate Imputations using sequential regression with L2 penalization. 
#'
#'
#'
#'Generates multiple imputations for incomplete multivariate data by fitting a 
#'sequence of regression models using L2 penalty iteratively. Missing data
#'can occur in one or more variables of the data. In each step of the iteration,
#'ridge regression is fitted according to the distributional form of the missing
#'variable taken as a response. All other variables are taken as predictors. 
#'If some predictors are incomplete, the most #'recently generated imputations 
#'are used to complete the predictors before using them as a predictor.
#'
#'
#'    
#'
#' @param x A data frame or a matrix containing the incomplete data.  Missing
#'values are coded as \code{NA}.

#' @param x.select A Boolean flag. If \code{TRUE}, linearly dependent columns 
#' will be removed before fitting of each imputation model. If \code{FALSE}, 
#' the linearly dependent columns will be removed only when number of predictors
#' is greater than the sample size for fitting an imputation model. The default 
#' is \code{FALSE}.
#'
#' @param pen A Boolean flag. If \code{TRUE}, each imputation model will be
#'fitted with L2 penalty. If \code{FALSE}, maximum likelihood estimation
#'(MLE) will be used. However, if MLE fails, L2 penalty is used for
#' fitting the imputation model. The default is \code{FALSE}.
#'
#' @param maxit A scalar giving the number of iterations. The default is 5.
#' @param m Number of multiple imputations. The default is \code{m=5}.
#' @param track A Boolean flag. If \code{TRUE}, \code{mispr} will print
#' additional information about iterations on console. The default is 
#' \code{FALSE} for silent computation.
#'
#' @param init.method Method for initialization of missing values. 
#' \code{random} for filling \code{NA} in each column with a random sample from 
#' the observed values of that column. \code{median} for mean imputation.
#'
#' @param L2.fix Fixed value of ridge penalty (optional) to use for each 
#' imputation model. For default i.e., \code{NULL}, L2 penalty will be decided
#' with k-fold cross-validation.
#'
#'@param cv A Boolean flag. If \code{TRUE} that is default, optimal value of L2
#'penalty will be decided indepndently for each imputation model using 5-fold
#'cross-validation.
#'
#' @param maxL2 The maximum value of the tuning parameter for L2 penalization
#' to be used for optimizing the cross-validated likelihood. Default value is
#' $2^10$.

#' @return a list containing the number of imputed datasets, number of iterations used to obtain imputed data, list of multiply imputed datasets, and summary of missing values.

#' @author Faisal Maqbool Zahid 
#' \email{faisalmz99@yahoo.com}.

#' @references Zahid, F. M., and Heumann, C. (2018). Multiple imputation with 
#' sequential penalized regression. 
#' \emph{Statistical Methods in Medical Research}, 0962280218755574.
#'
#'
#' @export
#'
#' @examples
#' data(data1)
#' # Select a subset of data1 
#' x=data1[ , 1:10]
#' res1 = mispr(x)
#' # to get 3 multiply imputed datasets
#' res2 = mispr(x, m=3)
#'
#' @import MASS
#' @import e1071 
#' @import penalized  
#' @importFrom stats rnorm
#' @importFrom stats runif
#' @importFrom stats glm
#' @importFrom stats rchisq
#' @importFrom utils head
#' @importFrom stats as.formula
#' @importFrom stats cor
#' @importFrom stats var




mispr <- function(x, x.select=FALSE, pen=FALSE, maxit=5, m=5, track=FALSE, init.method="random", L2.fix=NULL, cv=TRUE, maxL2=2^10)
{

if(!is.data.frame(x)){
    if(is.matrix(x))
      x <- as.data.frame(x)
    else
      stop("data frame must be provided")
                       }
                       
if(!is.null(L2.fix)) {if( L2.fix<0 | !is.numeric(L2.fix)) stop("L2 penalty must be non.negative value")}
if(!is.null(L2.fix) & cv==TRUE) {
                      cat("cv option is set automatically to FALSE  since fixed value of L2 penalty is given")
                      cv <- FALSE
                                 }

types <- lapply(x,class1)        

  attributes(types)$names <-NULL      
  types <- unlist(types)      

    if(any(types=="character")){
    chrInd <- which(types=="character")
    warning("At least one character variable is converted into a factor")
    for(ind in chrInd){
      x[,ind] <- as.factor(x[,ind])
      types[ind] <- "factor"
    }
  }
################################################################################


  indFac <- which(types=="factor")
  for(ind in indFac){
    if(length(levels(x[,ind]))==2)
      types[ind] <- "binary"
    else if(length(levels(x[,ind]))>2)
      stop("factor with more than 2 levels detected and has not yet implemented!")
    else stop("factor with less than 2 levels detected!")
  }
################################################################################

missingSummary <- data.frame(types,apply(x,2,function(x)sum(is.na(x))))
colnames(missingSummary) <- c("type","#missing")

################################################################################

N <- n <- dim(x)[1]       #No. of observations
P <- dim(x)[2]            #No. of variables/columns in x
  ## error management:
  if(dim(x)[2] < 2) stop("Less than 2 variables included in x.")

if(!any(is.na(x))) print("No missings in x. Nothing to impute")                   
if(any(apply(x, 1, function(x) all(is.na(x))))) stop("Unit non-responses included in x.")   

  factors <- vector()
  for(i in 1:ncol(x)){
    factors <- c(factors,is.factor(x[,i]))
  }
################################################################################


  ## variables/columns numbers that include missings
  x.na <-vector()
  for(i in seq(P)){
    if(na.exist(x[,i]))
      x.na <- c(x.na,i)
                 }

  w2  <- is.na(x)

################################################################################

  # Function Initialization of missing values

# Test <- require(e1071)
#  if(!Test){
#    install.packages("e1071",dependencies=TRUE,repos="http://cran.us.r-project.org")
#    warning("Package e1071 is installed on the run")
#    require(e1071)
#  }
#  
#  Test <- require(penalized)
#  if(!Test){
#    install.packages("penalized",dependencies=TRUE,repos="http://cran.us.r-project.org")
#    warning("Package penalized is installed on the run")
#    require(penalized)
#  }
#  
#  Test <- require(MASS)
#  if(!Test){
#    install.packages("MASS",dependencies=TRUE,repos="http://cran.us.r-project.org")
#    warning("Package MASS is installed on the run")
#    require(MASS)
#  }

if(requireNamespace("e1071"))
if(requireNamespace("MASS"))
if(requireNamespace("penalized"))
  res =list()
#  res$data = x
  res$m = m
  res$iteration = maxit
  res$imputation <- list()
  res$missingSummary <- missingSummary
  
if(m>0){
  

  data.input = x
  for(imp in 1:m){
    x=data.input
      if(init.method=="median"){
      x.init <- initialise(x,method="median")
    } else if(init.method=="random"){
      x.init <- initialise(x,method="random")
    } else {stop("The value of init.method is misspecified, please choose any of two options i.e., random or median") }
    x=x.init
    rm(x.init)
    if(track) print(head(x))
   ##############

    for(iteration in 1:maxit){

      for(i in x.na){
        if(track){
          print(paste("iteration ",iteration, "; imputed dataset", imp  ))
         # if(Sys.info()[1] == "Windows") flush.console()
        }
        yPart <- x[, i, drop=FALSE]
        wy <- which(w2[,i])      
        xPart <- x[, -i, drop=FALSE]

        dataForReg <- data.frame(yPart,xPart)
        

        if( types[i]=="numeric" ){ 
          meth = "numeric"
        } else if( types[i]=="binary" ){
          meth = "bin"
          }

        if(length(wy) > 0){

          if(track){ if(meth=="bin") print(paste("fitting model with ", colnames(dataForReg)[1], " as a binary response y.imp" )) else print(paste("fitting model with ", colnames(dataForReg)[1], " as a continuous response y.imp"))}
          

          colnames(dataForReg)[1] <- "y.imp"
        
                     x[wy,i] <- lm.glm(xReg=dataForReg, x.select=x.select, pen=pen, index=wy, type=meth,   L2.fix=L2.fix, cv=cv, maxL2=maxL2, track=track)

        }
              }  
    } 
     res$imputation[[imp]] <- x
  } 
} else       
{
stop("value of m should be positive")

}    

return(res)

  }   # End of 'mispr' function 
##################################################################################################################
