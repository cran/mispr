% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/misprImpute.R
\name{mispr}
\alias{mispr}
\title{Multiple Imputation with Sequential Penalized Regression}
\usage{
mispr(x, x.select = FALSE, pen = FALSE, maxit = 5, m = 5,
  track = FALSE, init.method = "random", L2.fix = NULL, cv = TRUE,
  maxL2 = 2^10)
}
\arguments{
\item{x}{A data frame or a matrix containing the incomplete data.  Missing
values are coded as \code{NA}.}

\item{x.select}{A Boolean flag. If \code{TRUE}, linearly dependent columns 
will be removed before fitting of each imputation model. If \code{FALSE}, 
the linearly dependent columns will be removed only when number of predictors
is greater than the sample size for fitting an imputation model. The default 
is \code{FALSE}.}

\item{pen}{A Boolean flag. If \code{TRUE}, each imputation model will be
fitted with L2 penalty. If \code{FALSE}, maximum likelihood estimation
(MLE) will be used. However, if MLE fails, L2 penalty is used for
fitting the imputation model. The default is \code{FALSE}.}

\item{maxit}{A scalar giving the number of iterations. The default is 5.}

\item{m}{Number of multiple imputations. The default is \code{m=5}.}

\item{track}{A Boolean flag. If \code{TRUE}, \code{mispr} will print
additional information about iterations on console. The default is 
\code{FALSE} for silent computation.}

\item{init.method}{Method for initialization of missing values. 
\code{random} for filling \code{NA} in each column with a random sample from 
the observed values of that column. \code{median} for mean imputation.}

\item{L2.fix}{Fixed value of ridge penalty (optional) to use for each 
imputation model. For default i.e., \code{NULL}, L2 penalty will be decided
with k-fold cross-validation.}

\item{cv}{A Boolean flag. If \code{TRUE} that is default, optimal value of L2
penalty will be decided indepndently for each imputation model using 5-fold
cross-validation.}

\item{maxL2}{The maximum value of the tuning parameter for L2 penalization
to be used for optimizing the cross-validated likelihood. Default value is
$2^10$.}
}
\value{
a list containing the number of imputed datasets, number of iterations used to obtain imputed data, list of multiply imputed datasets, and summary of missing values.
}
\description{
Generates Multivariate Imputations using sequential regression with L2 penalization.
}
\details{
Generates multiple imputations for incomplete multivariate data by fitting a 
sequence of regression models using L2 penalty iteratively. Missing data
can occur in one or more variables of the data. In each step of the iteration,
ridge regression is fitted according to the distributional form of the missing
variable taken as a response. All other variables are taken as predictors. 
If some predictors are incomplete, the most #'recently generated imputations 
are used to complete the predictors before using them as a predictor.
}
\examples{
data(data1)
# Select a subset of data1 
x=data1[ , 1:10]
res1 = mispr(x)
# to get 3 multiply imputed datasets
res2 = mispr(x, m=3)

}
\references{
Zahid, F. M., and Heumann, C. (2018). Multiple imputation with 
sequential penalized regression. 
\emph{Statistical Methods in Medical Research}, 0962280218755574.
}
\author{
Faisal Maqbool Zahid 
\email{faisalmz99@yahoo.com}.
}
