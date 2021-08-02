\name{imdh}
\alias{imdh}
\title{
Fast clustering with minimum density hyperplanes.
}
\description{
Obtains a divisive hierarchical clustering model using minimum density hyperplane separators. The function takes only a single mandatory argument, a n x d matrix of n observations on d variables. Can also be used to update a clustering model based on new observations.
}
\usage{
imdh(X, sol = NULL, depth = 8, hmult = 1, epoch = 2, C = 10, alpha = .1, tmax = 100000000, scale = 1, t_init0 = NULL, refine = TRUE, k = NULL)
}
\arguments{
  \item{X}{Matrix of data to be clustered with observations row-wise and variables column-wise.}
  \item{sol}{(optional) A clustering solution obtained from the imdh function. If supplied the solution will be updated using the observations in X. If omitted then a new solution will be created.}
  \item{depth}{(optional) integer maximum depth of hierarchical clustering model. Default is 8, meaning a maximum of 128 = 2^(8-1) clusters.}
  \item{hmult}{(optional) positive numeric bandwidth multiplier for smoothing kernel used to approximate stochastic gradients. Default is 1.}
  \item{epoch}{(optional) integer number of passes through the data allowed. Default is 2. Because the clustering solution is particularly variable early on it is recommended that either at least 2 passes are used, or if 1 pass is used then a subsequent assignment be performed using function imdh_assign}
  \item{C}{(optional) positive numeric used as a Lagrange type multiplier. Recommend leaving at default value of 10}
  \item{alpha}{(optional) non-negative numeric}
  \item{hnorm}{(optional) to make implementation compatible with existing R functions for bandwidth selection, which usually by default provide the bandwidth value for the Gaussian kernel. If hnorm is provided then fkde will convert this to the corresponding value for the kernel being used, based on the AMISE value relative to that of the Gaussian.}
}
\examples{
# generate a sample of size 1000000 from heavily skewed density
x <- dataGen(1000000, 5)
# obtain exact kernel estimate and a binned approximation on a grid of 1000 points
system.time(f_exact <- fkde(x, ngrid = 1000))
system.time(f_binned <- fkde(x, nbin = 1000))
# estimate the MISE
MISE(5, f_exact$x, f_exact$y)
MISE(5, f_binned$x, f_binned$y)
}
