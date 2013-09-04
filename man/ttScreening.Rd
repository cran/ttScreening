\name{ttScreening}
\alias{ttScreening}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
 A screening process built upon training and testing samples
}
\description{
 A screening process to filter out non-informative DNA methylation sites by applying (robust) regressions to training data, and the results are further examined using testing samples. Surrogate variables are included to account for unknown factors.
}
\usage{
ttScreening(y = y, x1 = x1, x2 = x2, type = c("numeric", "factor"), 
interaction = c(TRUE, FALSE), iterations = 100, cv.cutoff = 85, 
n.sv = NULL, train.alpha = 0.05, test.alpha = 0.1, FDR.alpha = 0.05, 
Bon.alpha = 0.05, percent = (2/3), method = c("TT", "FDR", "Bonferroni"), 
linear = c("robust", "ls"), vfilter = NULL, B = 5, numSVmethod = "be", 
rowname = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
    Data matrix of logit transformed DNA methylation measures (m by n, m subjects and n CpG sites). Each row represents DNA methylation measures of all CpG sites for one subject.
}
  \item{x1}{
  a covariate with dimension m by 1. It can be a categorical or numeric.
}
  \item{x2}{
   Categorical variable with dimension m by 1 or NULL. This is to fit the situation of evaluating conditional methQTL with x2 being SNPs or interaction between x1 and x2.
}
  \item{type}{
    Type of variable x1, "numeric" or "categorical".
}
  \item{interaction}{
     logic values ("TRUE" or "FALSE").
}
  \item{iterations}{
   Number of loops for the training/testing (TT) procedure. The default is 100.
}
  \item{cv.cutoff}{
    The minimum frequency required for a DNA methylation site to be treated as an informative site. After "iterations" iterations, the frequency of each DNA methylation being selected out of "iterations" iterations is recorded. The higher the frequency, the more likely the site is informative. The default is 85.
}
  \item{n.sv}{
  Number of surrogate variables. If NULL, the number of surrogate variables will be determined based on the data. The default is NULL.
}
  \item{train.alpha}{
     Significance level for training samples. The default is 0.05.
}
  \item{test.alpha}{
  Significance level for testing samples. The default is 0.1.
}
  \item{FDR.alpha}{
    False discovery rate. The default is 0.05. This is to fit the need of selecting variables based on FDR.
}
  \item{Bon.alpha}{
    Overall significance level by use of the Bonferroni method for mulitple testing correction. The default is 0.05. This is to fit the need of selecting variables based on the Bonferroni multiple testing correction.
}
  \item{percent}{
     Proportion of the full sample to be used for training. The default is 2/3.
}
  \item{method}{
    Method for filtering, "TT" (training/testing), "FDR", or "Bonferroni".
}
  \item{linear}{
     Linear regression methods, "robust" (robust regression) or "ls" (ordinary least squares).
}
  \item{vfilter}{
    The number of most variable CpG sites to use when building SVs, must be between
100 and the number of genes; Must be NULL or numeric (> 0), The default is NULL. 
 }
  \item{B}{
   Number of iterations in generating surrogate variables. The default is 5.
}
  \item{numSVmethod}{
   The method for determining the number of surrogate variables to use. The default is "be", the other method is "leek".
}
\item{rowname}{
 Optional, NULL or "TRUE". The default is NULL. If rownames are not already present with the data, the order in which the DNA methylation sites are listed will become the rowname. Surrogate variable estimates are formed based on the algorithms in Leek and Storey (2007).
}
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
 \item{train.cpg}{
Number of DNA methylation sites selected after the training step of each loop.
}
  \item{test.cpg}{
Number of DNA methylation sites selected after the testing step of each loop.
}
  \item{selection}{
Indicator matrix for the TT method after the testing step. The number of rows is the number of methylation sites, and the number of columns is the number of iterations. An entry of 1 indiates the selection of a site, and 0 otherwise.
}
  \item{pvalue.matrix}{
Matrix of p-values of the selected DNA methylation sites after the testing step. The number of rows is the number of methylation sites and the number of columns is the number of iterations. For methylation sites not selected, NA is listed.
}
  \item{TT.cpg}{
Final list of the DNA methylation sites by their original rownames selected from the TT method.
}
  \item{FDR.cpg}{
Final list of the DNA methylation sites by their original rownames selected from the FDR method.
}
  \item{Bon.cpg}{
Final list of the DNA methylation sites by their original rownames selected from the Bonferroni method.
}
  \item{TT.output}{
Dataframe containing the list of DNA methylation sites selected from the TT method and the respective coefficients and pvalues for the variables and SVs.
}
  \item{FDR.output}{
Dataframe containing the list of DNA methylation sites selected from the FDR method and the respective coefficients and pvalues for the variables and SVs.
}
  \item{Bon.output}{
Dataframe containing the list of DNA methylation sites selected from the Bonferroni method and the respective coefficients and pvalues for the variables and SVs.
}

}
\references{
Meredith Ray, Xin Tong, Hongmei Zhang, and Wilfred Karmaus. (2013) "DNA methylation sites screening with surrogate variables", unpublished manuscript.

Leek JT and Storey JD. (2007) "Capturing heterogeneity in gene expression studies by `Surrogate Variable Analysis'." PLoS Genetics, 3: e161.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
  \code{\link{sva}}
}
\examples{
\dontrun{
##### Simulated Data Analysis #######
x1= rnorm(250,2,15)
size1<-rmultinom(1,250,c(0.4,0.4,0.2))
x2= matrix(sample(c(rep(0,size1[1,]),
			rep(1,size1[2,]),rep(2,size1[3,])),replace=F),byrow=250,ncol=1)
sigma1<-diag(250)
beta0<-3.5
beta1<-2.5
beta2<-2.5
beta3<-0.5
error<-mvrnorm(1,mu=rep(0,11), Sigma=diag(11))
beta4<-rnorm(250,3,10)
beta5<-rnorm(1,3,2)
beta6<-rnorm(1,2,3)
constant<-rmultinom(1,250,c(0.35,0.65))
beta7<-sample(c(rep(0,constant[1,1]),rep(1,constant[2,1])))
y1<- mvrnorm(10,mu=(beta0 + beta1*x1 + beta2*x2 + beta3*(x1*x2)),Sigma=sigma1)
y2<- mvrnorm(1,mu=(beta0 + beta1*x1 + beta3*(x1*x2)),Sigma=sigma1)
y4<- mvrnorm(1989,mu=(beta0 + beta5*beta4 + beta6*beta4*beta7),Sigma=sigma1)
y=rbind(y1,y2,y4)+c(error,rep(0,1989))


runs<-ttScreening(y=y,x1=x1,x2=x2,type="numeric",interaction="TRUE",iterations=100,
cv.cutoff=85,n.sv=NULL,train.alpha=0.05,test.alpha=0.1, FDR.alpha=0.05,Bon.alpha=0.05,
percent=(2/3), method="TT",linear="robust",vfilter = NULL, B = 5, numSVmethod = "be",
rowname=NULL)

runs$TT.output
}
}