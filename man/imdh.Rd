\name{imdh}
\alias{imdh}
\title{
Fast clustering with minimum density hyperplanes.
}
\description{
Obtains a divisive hierarchical clustering model using minimum density hyperplane separators. The function takes only a single mandatory argument, a n x d matrix of n observations on d variables. Can also be used to update a clustering model based on new observations.
}
\usage{
imdh(X, sol = NULL, depth = 7, hmult = 1, q = 1/6, r = 1, epoch = 2, C = 1, alpha = .25, tmax = 100000000, scale = 1, t_init0 = NULL, refine = TRUE, forward = FALSE, k = NULL)
}
\arguments{
  \item{X}{Matrix of data to be clustered with observations row-wise and variables column-wise.}
  \item{sol}{(optional) A clustering solution obtained from the imdh function. If supplied the solution will be updated using the observations in X. If omitted then a new solution will be created.}
  \item{depth}{(optional) integer maximum depth of hierarchical clustering model. Default is 7, meaning a maximum of 128 = 2^7 clusters.}
  \item{hmult}{(optional) positive numeric bandwidth multiplier for smoothing kernel used to approximate stochastic gradients. Default is 1.}
  \item{q}{(optional) positive numeric giving the rate of decrease of the sequence of bandwidths. Default is 1/6 meaning decay rate is t^(-q) = t^(-1/6)}
  \item{r}{(optional) positive numeric giving the rate of decrease of the sequence of step-sizes. Default is 1 meaning a decay rate of t^(-r) = t^(-1).}
  \item{epoch}{(optional) integer number of passes through the data allowed. Default is 2. Because the clustering solution is particularly variable early on it is recommended that either at least 2 passes are used, or if 1 pass is used then a subsequent assignment be performed using function imdh_assign}
  \item{C}{(optional) outdated, don't change from default.}
  \item{alpha}{(optional) non-negative numeric}
  \item{tmax}{(optional) integer maximum effective sample size. Introduces forgetting to stop convergence if underlying distribution might be non-stationary.}
  \item{scale}{(optional) should the variables be scaled to have variance 1? 1 for "yes" and 0 for "no". Default is 1/yes.}
  \item{t_init0}{(optional) integer sample size to obtain initial estimate of mean and standard deviation in each variable. Default is number of variables/dimensions.}
  \item{refine}{(optional) should the complete binary tree be pruned? Default is TRUE.}
  \item{forward}{(optional) should the pruning be based on forward addition of nodes or backward removal of nodes? Default is FALSE.}
  \item{k}{(optional) number of clusters desired in solution. Default is 2^depth if refine is FALSE and automatically selected if refine is TRUE.}
}
\examples{
# load pendigits data set
data('pendigits')

# obtain clustering solution
imdh_model <- imdh(pendigits$x)

# plot projection of data onto first two projection vectors in model with colours according to clusters and true classes
par(mfrow = c(1, 2))
plot(pendigits$x\%*\%imdh_model$V, col = imdh_model$clusters)
plot(pendigits$x\%*\%imdh_model$V, col = pendigits$c+1)
}
