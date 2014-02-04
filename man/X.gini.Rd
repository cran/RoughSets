\name{X.gini}
\alias{X.gini}
\title{The gini-index gain measure function}
\usage{
  X.gini(decisionDistrib)
}
\arguments{
  \item{decisionDistrib}{a distribution of values on the
  decision attribute.}
}
\value{
  \code{qualityF}.
}
\description{
  It is an auxiliary function for the \code{qualityF}
  parameter in the
  \code{\link{FS.greedy.heuristic.reduct.RST}} and
  \code{\link{FS.greedy.heuristic.superreduct.RST}}
  functions. It is based on the \emph{gini} index as a
  measure of information (K. Stoffel and L. E. Raileanu,
  2000).
}
\references{
  K. Stoffel and L. E. Raileanu, "Selecting Optimal
  Split-Functions with Linear Threshold Unit Trees and
  Madaline-Style Networks", in: Research and Development in
  Intelligent Systems XVII, BCS Conference Series (2000).
}

