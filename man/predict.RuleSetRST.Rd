\name{predict.RuleSetRST}
\alias{predict.RST}
\alias{predict.RuleSetRST}
\title{The predicting function for rule induction methods based on RST}
\usage{
\method{predict}{RuleSetRST}(object, newdata, votingScheme = X.bestFirst, ...)
}
\arguments{
  \item{object}{a \code{"RuleSetRST"} class resulted by
  \code{\link{RI.indiscernibilityBasedRules.RST}}.}

  \item{newdata}{a \code{"DecisionTable"} class containing
  a data frame or matrix (m x n) of data for the prediction
  process, where m is the number of instances and n is the
  number of input attributes. It should be noted that this
  data must have \code{colnames} on each attribute.}

  \item{votingScheme}{a function used to choose predicted
  values. The default value is \code{\link{X.bestFirst}}.}

  \item{...}{the other parameters.}
}
\value{
The predicted values.
}
\description{
It is a function used to obtain predicted values after
obtaining rules by using rule induction methods. We have
provided the functions
\code{\link{RI.indiscernibilityBasedRules.RST}} to generate
rules based on RST.
}
\examples{
##############################################
## Example: Classification Task
##############################################
data(RoughSetData)
decision.table <- RoughSetData$hiring.dt

#### Rule induction based on RST ####
## determine feature subset/reduct
reduct <- FS.permutation.heuristic.reduct.RST(decision.table,  permutation = NULL)

## generate rules
rules.rst <- RI.indiscernibilityBasedRules.RST(decision.table, reduct)

## predicting newdata
## in this case, we are using the same dataset as training data
res.1 <- predict(rules.rst, decision.table[, -ncol(decision.table)])
}
\seealso{
\code{\link{RI.GFRS.FRST}} and
\code{\link{RI.hybridFS.FRST}}.
}

