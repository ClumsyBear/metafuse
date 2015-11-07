\name{metafuse-package}
\alias{metafuse-package}
\docType{package}
\title{\packageTitle{metafuse}}
\description{\packageDescription{metafuse}}

\details{
Very simple to use. Accepts X, y, and sid (study ID) for regression models. 
\packageIndices{metafuse}
}

\author{
\packageAuthor{metafuse}
Maintainer: \packageMaintainer{metafuse}
}

\references{
Lu Tang, Peter X.K. Song. Fused Lasso Approach in Regression Coefficients Clustering - Learning Parameter Heterogeneity in Data Integration. In preparation.

Fei Wang. Development of Joint Estimating Equation Approaches to Merging Clustered or Longitudinal Datasets from Multiple Biomedical Studies. PhD Thesis, The University of Michigan, 2012. (A Biometrics paper is tentatively accepted. We will revise this a later time.)

Jiahua Chen and Zehua Chen. Extended bayesian information criteria for model selection with large model spaces. Biometrika, 95(3):759-771, 2008.

Xin Gao and Peter X-K Song. Composite likelihood bayesian information criteria for model selection in high-dimensional data. Journal of the American Statistical Association, 105(492):1531-1540, 2010.
}


\keyword{ package }

\examples{
########### Generate Data ########### 
n <- 200    # sample size in each study
K <- 10     # number of studies
p <- 3      # number of covariates in X (including intercept)
N <- n*K    # total sample size

# the coefficient matrix, used this to set desired heterogeneous pattern (depends on p and K)
beta0 <- matrix(c(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, # intercept
                  0.0,0.0,0.0,0.0,0.0,1.0,1.0,1.0,1.0,1.0, # beta_1, etc.
                  0.0,0.0,0.0,0.0,0.5,0.5,0.5,1.0,1.0,1.0), K, p)

# generate a data set, family=c("gaussian", "binomial", "poisson")
data <- datagenerator(n=n, beta0=beta0, family="gaussian", seed=123)

# prepare the input (y, X, studyID)
y       <- data$y
sid     <- data$group
X       <- data[,-c(1,ncol(data))]

########### metafuse runs ###########
# fuse slopes of X1 (it is heterogeneous with 2 groups)
metafuse(X=X, y=y, sid=sid, fuse.which=c(1), family="gaussian", intercept=TRUE, alpha=0,
          criterion="EBIC", verbose=TRUE, plots=TRUE, loglambda=TRUE)

# fuse slopes of X2 (it is heterogeneous with 3 groups)
metafuse(X=X, y=y, sid=sid, fuse.which=c(2), family="gaussian", intercept=TRUE, alpha=0,
          criterion="EBIC", verbose=TRUE, plots=TRUE, loglambda=TRUE)

# fuse all three covariates
metafuse(X=X, y=y, sid=sid, fuse.which=c(0,1,2), family="gaussian", intercept=TRUE, alpha=0,
          criterion="EBIC", verbose=TRUE, plots=TRUE, loglambda=TRUE)

# fuse all three covariates, with sparsity penalty
metafuse(X=X, y=y, sid=sid, fuse.which=c(0,1,2), family="gaussian", intercept=TRUE, alpha=1,
          criterion="EBIC", verbose=TRUE, plots=TRUE, loglambda=TRUE)

# fit metafuse at a given lambda
metafuse.l(X=X, y=y, sid=sid, fuse.which=c(0,1,2), family="gaussian", intercept=TRUE, 
          alpha=1, lambda=0.5)
}
