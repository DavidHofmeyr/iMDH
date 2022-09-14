\name{imdh_assign}
\alias{imdh_assign}
\title{
Assign data to clusters using an imdh model.
}
\description{
Performs cluster assignments from an imdh model without updating model parameters as passing through data.
}
\usage{
imdh_assign(X, sol, refine = TRUE, forward = FALSE, k = NULL)
}
\arguments{
  \item{X}{Matrix of data to be clustered with observations row-wise and variables column-wise.}
  \item{sol}{A clustering solution obtained from the imdh function.}
  \item{refine}{(optional) should the complete binary tree be pruned? Default is TRUE.}
  \item{forward}{(optional) should the pruning be based on forward addition of nodes or backward removal of nodes? Default is FALSE.}
  \item{k}{(optional) number of clusters desired in solution. Default is 2^(depth-1) if refine is FALSE and automatically selected if refine is TRUE.}
}
\examples{
# load pendigits data set
data('pendigits')

# obtain clustering solution with a single pass through the data. The assignments made early on will be from the model prior to any sort of
# convergence and hence we expect relatively poor performance
imdh_model <- imdh(pendigits$x, epoch = 1)

# make assignments using the fixed model obtained above
imdh_model2 <- imdh_assign(pendigits$x, imdh_model)

# plot projection of data onto first two projection vectors in models with colours according to clusters and also true classes
par(mfrow = c(1, 2))
plot(pendigits$x\%*\%imdh_model$V, col = imdh_model$clusters)
plot(pendigits$x\%*\%imdh_model2$V, col = imdh_model2$clusters)
plot(pendigits$x\%*\%imdh_model$V, col = pendigits$c+1)
}
