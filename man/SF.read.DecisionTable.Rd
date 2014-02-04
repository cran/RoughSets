\name{SF.read.DecisionTable}
\alias{SF.read.DecisionTable}
\title{The importing function}
\usage{
  SF.read.DecisionTable(filename, decision.attr = NULL,
    indx.nominal = NULL, ...)
}
\arguments{
  \item{filename}{a file name that contains objects/data.}

  \item{decision.attr}{a index position of decision
  attribute. When we ignore this parameter, the function
  will treat the data as an information system or
  newdata/testing data. In other words, we must define the
  index of decision attribute if we have the decision
  attribute in the dataset.}

  \item{indx.nominal}{a indexes of nominal values. It is
  used to define specific indexes as nominal attributes.}

  \item{...}{other parameters which are involved in
  \code{read.table} function. See
  \code{\link{read.table}}.}
}
\value{
  A \code{"DecisionTable"} class which is the standard
  decision table. See \code{\link{SF.asDecisionTable}}.
}
\description{
  It is a additional function used to import dataset from
  files and then construct them into the standard decision
  table.
}
\details{
  It is used to read a dataset from files with the
  following format: ".csv", ".txt", ".dat", and ".xls",
  into decision table format. The data should be in table
  format containing rows and columns without header, where
  every row represents every object/instance while columns
  represent attributes of the objects.

  The output of this function is a decision table which
  fulfills a standard format of RoughSets package (See
  \code{\link{SF.asDecisionTable}}). In order to construct
  the decision table, especially in determining description
  of attributes, this function uses simple heuristic
  techniques as follows: \itemize{ \item string values:
  they will be recognized as nominal/symbolic values. \item
  integer values: they will be recognized as
  continuous/real values. \item decimal/real values: they
  will be recognized as continuous/real values. \item
  indx.nominal: the attributes that are set in this
  parameter will be assigned as nominal values. }
}
\examples{
#############################################################
## Example 1: dataset saved in a file
#############################################################
## Let us assume we have the following data which has been already saved to the file "tes.dat"
data <- data.frame(c(0.12, 0.23, 0.24), c(1,3,2), c(10, 12, 18), c("a", "a", "b"), c(1, 1, 2))
\dontrun{write.table(data, file = "tes.dat", row.names = FALSE, col.names = FALSE,
                    fileEncoding ="")}

## Then we would generate decision table from tes.dat file.
## in this case, we want to define that second and third attributes are nominal and continuous,
## respectively.
\dontrun{decision.table <- SF.read.DecisionTable(filename = "tes.dat", decision.attr = 4,
                  indx.nominal = 2, sep= " ", col.names = c("v1", "v2", "v3", "v4", "o1"))}
}

