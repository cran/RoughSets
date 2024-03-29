% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BasicRoughSets.R
\name{BC.boundary.reg.RST}
\alias{BC.boundary.reg.RST}
\title{Computation of a boundary region}
\usage{
BC.boundary.reg.RST(decision.table, roughset)
}
\arguments{
\item{decision.table}{an object inheriting from the \code{"DecisionTable"} class, which represents a decision system. 
See \code{\link{SF.asDecisionTable}}.}

\item{roughset}{an object inheriting from the \code{"LowerUpperApproximation"} class, which represents
lower and upper approximations of decision classes in the data. Such objects are typically produced by calling 
the \code{\link{BC.LU.approximation.RST}} function.}
}
\value{
An object of a class \code{"BoundaryRegion"} which is a list with the following components:
        \itemize{
          \item \code{boundary.reg}: an integer vector containing indices of data instances belonging 
                to the boundary region,
          \item \code{degree.dependency}: a numeric value giving the degree of dependency,
          \item \code{type.model}: a varacter vector identifying the utilized model. In this case, 
                it is \code{"RST"} which means the rough set theory.       
        }
}
\description{
This function implements a fundamental part of RST: computation of a boundary region and the
degree of dependency. This function can be used as a basic building block for development 
of other RST-based methods. A more detailed explanation of this notion can be found 
in \code{\link{Introduction-RoughSets}}.
}
\examples{
########################################################
data(RoughSetData)
hiring.data <- RoughSetData$hiring.dt

## We select a single attribute for computation of indiscernibility classes:
A <- c(2)

## compute the indiscernibility classes:
IND.A <- BC.IND.relation.RST(hiring.data, feature.set = A)

## compute the lower and upper approximation:
roughset <- BC.LU.approximation.RST(hiring.data, IND.A)

## get the boundary region:
pos.boundary = BC.boundary.reg.RST(hiring.data, roughset)
pos.boundary

}
\references{
Z. Pawlak, "Rough Sets", International Journal of Computer and Information Sciences, 
vol. 11, no. 5, p. 341 - 356 (1982).
}
\seealso{
\code{\link{BC.IND.relation.RST}}, \code{\link{BC.LU.approximation.RST}}, \code{\link{BC.LU.approximation.FRST}}
}
\author{
Dariusz Jankowski, Andrzej Janusz
}
