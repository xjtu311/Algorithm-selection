\name{successes}
\alias{successes}
\title{
Success
}
\description{
Was the problem solved successfully using the chosen algorithm?
}
\usage{
successes(data, model, timeout, addCosts = NULL)
}
\arguments{
  \item{data}{the data used to induce the model. The same as given to
      \code{classify}, \code{classifyPairs}, \code{cluster} or
      \code{regression}.}
  \item{model}{the algorithm selection model. Can be either a model
      returned by one of the model-building functions or a function that returns
      predictions such as \code{vbs} or the predictor function of a trained
      model.}
  \item{timeout}{the timeout value to be multiplied by the penalization factor.
      If not specified, the maximum performance value of all algorithms on the
      entire data is used.}
  \item{addCosts}{whether to add feature costs. You should not need to set this
      manually, the default of \code{NULL} will have LLAMA figure out
      automatically depending on the model whether to add costs or not. This
      should always be true (the default) except for comparison algorithms (i.e.
      single best and virtual best).}
}
\details{
Returns \code{TRUE} if the chosen algorithm successfully solved the problem
instance, \code{FALSE} otherwise for each problem instance.

If feature costs have been given and \code{addCosts} is \code{TRUE}, the cost of
the used features or feature groups is added to the performance of the chosen
algorithm. The used features are determined by examining the the \code{features}
member of \code{data}, not the model. If after that the performance value is
above the timeout value, \code{FALSE} is assumed. If whether an algorithm was
successful is not determined by performance and feature costs, don't pass costs
when creating the LLAMA data frame.

If the model returns \code{NA} (e.g. because no algorithm solved the instance),
\code{FALSE} is returned as success.

\code{data} may contain a train/test partition or not. This makes a difference
when computing the successes for the single best algorithm. If no train/test
split is present, the single best algorithm is determined on the entire data. If
it is present, the single best algorithm is determined on each test partition.
That is, the single best is local to the partition and may vary across
partitions.
}
\value{
A list of the success values.
}
\author{
Lars Kotthoff
}
\seealso{
\code{\link{misclassificationPenalties}}, \code{\link{parscores}}
}
\examples{
if(Sys.getenv("RUN_EXPENSIVE") == "true") {
data(satsolvers)
folds = cvFolds(satsolvers)

model = classify(classifier=makeLearner("classif.J48"), data=folds)
sum(successes(folds, model))
}
}
\keyword{ models }
