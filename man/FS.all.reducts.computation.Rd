\name{FS.all.reducts.computation}
\alias{FS.all.reducts.computation}
\title{The function for computing all reducts}
\usage{
  FS.all.reducts.computation(discernibilityMatrix)
}
\arguments{
  \item{discernibilityMatrix}{a
  \code{"DiscernibilityMatrix"} class representing the
  discernibility matrix of RST and FRST.

  See \code{\link{BC.discernibility.mat.RST}} and
  \code{\link{BC.discernibility.mat.FRST}}.}
}
\value{
  A class \code{"ReductSet"}.
}
\description{
  It is a function used for computation of all reducts from
  a discernibility matrix based on RST and FRST
}
\examples{
########################################################
## Example 1: Generate a single reduct and
##            new decision table using RST
########################################################
data(RoughSetData)
decision.table <- RoughSetData$hiring.dt

## build the decision-relation discernibility matrix
res.2 <- BC.discernibility.mat.RST(decision.table, range.object = NULL)

## generate single reduct
reduct <- FS.all.reducts.computation(res.2)

## generate new decision table
new.decTable <- SF.applyDecTable(decision.table, reduct, control = list(indx.reduct = 1))

##############################################################
## Example 2: Generate a single reduct and
##            new decision table using FRST
##############################################################
data(RoughSetData)
decision.table <- RoughSetData$hiring.dt

## build the decision-relation discernibility matrix
control.1 <- list(type.relation = c("crisp"),
                type.aggregation = c("crisp"),
                t.implicator = "lukasiewicz", type.LU = "implicator.tnorm")
res.1 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "standard.red",
                                    control = control.1)

## generate single reduct
reduct <- FS.all.reducts.computation(res.1)

## generate new decision table
new.decTable <- SF.applyDecTable(decision.table, reduct, control = list(indx.reduct = 1))
}

