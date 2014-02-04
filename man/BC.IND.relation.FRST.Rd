\name{BC.IND.relation.FRST}
\alias{BC.IND.relation.FRST}
\title{The indiscernibility relation based on fuzzy rough set theory}
\usage{
  BC.IND.relation.FRST(decision.table, attributes = NULL,
    control = list())
}
\arguments{
  \item{decision.table}{a \code{"DecisionTable"} class
  representing a decision table. See
  \code{\link{SF.asDecisionTable}}.}

  \item{attributes}{a numerical vector expressing indexes
  of subset of attributes to be considered.  The default
  value is \code{NULL} which means that all conditional
  attributes will be considered.}

  \item{control}{a list of other parameters which involves
  the following parameters.  \itemize{ \item
  \code{type.relation}: a vector containing string values
  that express the type of the fuzzy relation and its
  equation. The default value is \code{type.relation =
  c("tolerance", "eq.1")}. See in the Section
  \code{Details}.  \item \code{type.aggregation}: a list
  expressing type of aggregation. The default value is
  \code{type.aggregation = c("t.tnorm", "lukasiewicz")}.
  See in the Section \code{Details}.  }}
}
\value{
  A class \code{"IndiscernibilityRelation"} which contains
  \itemize{ \item \code{IND.relation}: a matrix
  representing indiscernibility relation over all objects.
  \item \code{type.relation}: a vector representing the
  type of relation.  \item \code{type.aggregation}: a
  vector representing the type of aggregation operator.
  \item \code{type.model}: a string showing the type of
  model which is used. In this case it is \code{"FRST"}
  which means fuzzy rough set theory.  }
}
\description{
  This is a function used to implement a fundamental
  concept of FRST which is fuzzy indiscernibility
  relations. It is used for any fuzzy relation that
  determines the degree to which two objects are
  indiscernibility. The detailed description about basic
  concepts of FRST can be found in
  \code{\link{B.Introduction-FuzzyRoughSets}}.
}
\details{
  Briefly, the indiscernibility relation is a relation that
  shows a degree of similarity among the objects. For
  example, \eqn{R(x_i, x_j) = 0} means the object \eqn{x_i}
  is completely different from \eqn{x_j}, and if
  \eqn{R(x_i, x_j) = 1}, while between those values we
  consider a degree of similarity. To calculate this
  relation, several methods have been implemented in this
  function which are approaches based on fuzzy tolerance,
  equivalence and \eqn{T}-equivalence relations. The fuzzy
  tolerance relations proposed by (R. Jensen and Q. Shen,
  2009) include three equations while (Q. Hu, 2004)
  proposed five \eqn{T_{cos}}-transitive kernel functions
  as fuzzy \eqn{T}-equivalence relations. The simple
  algorithm of fuzzy equivalence relation is implemented as
  well. Furthermore, we facilitate users to define their
  own equation for similarity relation.

  To calculate a particular relation, we should pay
  attention to several components in the parameter
  \code{control}. The main component in the \code{control}
  parameter is \code{type.relation} that defines what kind
  of approach we are going to use. The detailed explanation
  about the parameters and their equations is as follows:
  \itemize{ \item \code{"tolerance"}: It refers to fuzzy
  tolerance relations proposed by (R. Jensen and Q. Shen,
  2009).  In order to represent the \code{"tolerance"}
  relation, we must set \code{type.relation} as follows.

  \code{type.relation = c("tolerance", <chosen equation>)}

  where the chosen equation called as \code{t.similarity}
  is one of \code{"eq.1"}, \code{"eq.2"}, and \code{"eq.3"}
  equations which have been explained in
  \code{\link{B.Introduction-FuzzyRoughSets}}.

  \item \code{"transitive.kernel"}: It refers to
  \eqn{T}-equivalence relation proposed by (Q. Hu et al,
  2004).  In order to represent the relation, we must set
  \code{type.relation} parameter as follows.

  \code{type.relation = c("transitive.kernel", <chosen
  equation>, <delta>)}

  where the chosen equation is one of five following
  equations (called \code{t.similarity}): \itemize{ \item
  \code{"gaussian"}: It means gaussian kernel which is
  \eqn{R_G(x,y) = \exp (-\frac{|x - y|^2}{\delta})} \item
  \code{"exponential"}: It means exponential kernel which
  is \eqn{R_E(x,y) = \exp(-\frac{|x - y|}{\delta})} \item
  \code{"rational"}: It means rational quadratic kernel
  which is \eqn{R_R(x,y) = 1 - \frac{|x - y|^2}{|x - y|^2 +
  \delta}} \item \code{"circular"}: It means circular
  kernel which is if \eqn{|x - y| < \delta}, \eqn{R_C(x,y)
  = \frac{2}{\pi}\arccos(\frac{|x - y|}{\delta}) -
  \frac{2}{\pi}\frac{|x - y|}{\delta}\sqrt{1 - (\frac{|x -
  y|}{\delta})^2}} \item \code{"spherical"}: It means
  spherical kernel which is if \eqn{|x - y| < \delta},
  \eqn{R_S(x,y) = 1 - \frac{3}{2}\frac{|x - y|}{\delta} +
  \frac{1}{2}(\frac{|x - y|}{\delta})^3} } and \code{delta}
  is a specified value.  For example: let us assume we are
  going to use \code{"transitive.kernel"} as the fuzzy
  relation, \code{"gaussian"} as its equation and the delta
  is 0.2. So, we assign the \code{type.relation} parameter
  as follows:

  \code{type.relation = c("transitive.kernel", "gaussian",
  0.2)}

  If we omit the \code{delta} parameter then we are using
  \code{"gaussian"} defined as \eqn{R_E(x,y) =
  \exp(-\frac{|x - y|}{2\sigma^2})}, where \eqn{sigma} is
  the variance.  Furthermore, when using this relation,
  usually we set

  \code{type.aggregation = c("t.tnorm", "t.cos")}.

  \item \code{"transitive.closure"}: It refers to
  similarity relation (also called fuzzy equivalence
  relation). In this package, we consider a simple
  algorithm to calculate this relation as follows.

  Input: a fuzzy relation R

  Output: a min-transitive fuzzy relation \eqn{R^m}

  Algorithm:

  1. For every x, y: compute

  \eqn{R'(x,y) = max(R(x,y), max_{z \in U}min(R(x,y),
  R(z,y)))}

  2. If \eqn{R' \not= R}, then \eqn{R \gets R'} and go to
  1, else \eqn{R^m \gets R'}

  For interested readers, other algorithms can be seen in
  (H. Naessens et al, 2002). Let \code{"eq.1"} be the
  \eqn{R} fuzzy relations, to define it as parameter is

  \code{type.relation = c("transitive.closure", "eq.1")}.
  Instead of \code{"eq.1"}, we can also use other equations
  that have been explained in \code{"tolerance"} and
  \code{"transitive.kernel"}.

  \item \code{"crisp"}: It uses crisp equality for all
  attributes and we set the parameter \code{type.relation =
  "crisp"}. In this case, we only have \eqn{R(x_i, x_j) =
  0} which means the object \eqn{x_i} is completely
  different from \eqn{x_j}, and otherwise if \eqn{R(x_i,
  x_j) = 1}.

  \item \code{"custom"}: this value means that we define
  our own equation for the indiscernibility relation.  The
  equation should be defined in parameter
  \code{FUN.relation}.

  \code{type.relation = c("custom", <FUN.relation>)}

  The customized function should consist of three arguments
  which are \code{decision.table}, \code{x}, and \code{y},
  where \code{x} and \code{y} represent two objects which
  we want to compare.  It should be noted that the user
  must ensure that the values of this equation are alwayas
  between 0 and 1.  Furthermore, the function must have a
  return value. An example can be seen in Section
  \code{Examples}. }

  Beside the above \code{type.relation}, we provide several
  options of values for the \code{type.aggregation}
  parameter. The following is a description about it.
  \itemize{ \item \code{type.aggregation = c("crisp")}: It
  uses crisp equality for all attributes.  \item
  \code{type.aggregation = c("t.tnorm", <t.tnorm
  operator>)}: It means we are using \code{"t.tnorm"}
  aggregation which is a triangular norm operator with a
  specified operator \code{t.tnorm} as follows: \itemize{
  \item \code{"min"}: standard t-norm which is
  \eqn{min(x_1, x_2)}.  \item \code{"hamacher"}: hamacher
  product which is \eqn{(x_1 * x_2)/(x_1 + x_2 - x_1 *
  x_2)}.  \item \code{"yager"}: yager class which is \eqn{1
  - min(1, ((1 - x_1) + (1 - x_2)))}.  \item
  \code{"product"}: product t-norm which is \eqn{(x_1 *
  x_2)}.  \item \code{"lukasiewicz"}: lukasiewicz's t-norm
  (default) which is \eqn{max(x_2 + x_1 - 1, 0)}.  \item
  \code{"t.cos"}: \eqn{T_{cos}}t-norm which is \eqn{max(x_1
  * x_2 - \sqrt{1 - x_1^2}\sqrt{1 - x_2^2, 0})}.  \item
  \code{FUN.tnorm}: It is a user-defined function. It has
  to have two arguments, for example:

  \code{FUN.tnorm <- function(left.side, right.side)}

  \code{if ((left.side + right.side) > 1)}

  \code{return(min(left.side, right.side))}

  \code{else return(0)} } The default value is
  \code{type.aggregation = c("t.tnorm", "lukasiewicz")}.
  \item \code{type.aggregation = c("custom", <FUN.agg>)}:
  It is quite similar to the previous one but using this
  option users have more possibilities to define their own
  functions. \code{<FUN.agg>} is a function having one
  argument representing data that is produced by fuzzy
  similarity equation calculation.  The data is a list of
  one or many matrices which depend on the number of
  considered attributes and has dimension number of object
  x number of object. For example:

  \code{FUN.agg <- function(data) return(Reduce("+",
  data)/length(data))}

  which is a function to calculate average along all
  attributes. Then, we can set \code{type.aggregation} as
  follows:

  \code{type.aggregation = c("general.custom", <FUN.agg>)}.
  Another example can be seen in Section \code{Examples}.
  As we mentioned before, this option gives more
  flexibility to define functions since we do not need to
  consider the associative rule.  }

  Furthermore, the use of this function has been
  illustrated in Section \code{Examples}. Finally, this
  function is important since it is a basic function needed
  by other functions, such as
  \code{\link{BC.LU.approximation.FRST}} and
  \code{\link{BC.positive.reg.FRST}} for calculating lower
  and upper approximation and determining positive regions.
}
\examples{
###########################################################
## Example 1: Dataset containing nominal values for
## all attributes.
###########################################################
## Decision table is represented as data frame
dt.ex1 <- data.frame(c(1,0,2,1,1,2,2,0), c(0, 1,0, 1,0,2,1,1),
                        c(2,1,0,0,2,0,1,1), c(2,1,1,2,0,1,1,0), c(0,2,1,2,1,1,2,1))
colnames(dt.ex1) <- c("aa", "bb", "cc", "dd", "ee")
decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 5)

## In this case, we only consider the second and third attributes.
attributes <- c(2, 3)

#### calculate fuzzy indiscernibility relation ####
## in this case, we are using "crisp" as a type of relation and type of aggregation
control.ind <- list(type.relation = c("crisp"), type.aggregation = c("crisp"))
IND <- BC.IND.relation.FRST(decision.table, attribute = attributes, control = control.ind)

###########################################################
## Example 2: Dataset containing real-valued attributes
###########################################################
dt.ex2 <- data.frame(c(-0.4, -0.4, -0.3, 0.3, 0.2, 0.2),
                     c(-0.3, 0.2, -0.4, -0.3, -0.3, 0),
				        c(-0.5, -0.1, -0.3, 0, 0, 0),
				        c("no", "yes", "no", "yes", "yes", "no"))
colnames(dt.ex2) <- c("a", "b", "c", "d")
desc.attributes <- list(a = c(-0.4, 0.3), b = c(-0.4, 0.2),
                    c = c(-0.5, 0), d = c("no", "yes"))
decision.table <- SF.asDecisionTable(dataset = dt.ex2, decision.attr = 4)

## in this case, we only consider the first and second attributes
attributes <- c(1, 2)

#### Calculate fuzzy indiscernibility relation ####
## in this case, we choose "tolerance" relation and "eq.1" as similarity equation
## and "lukasiewicz" as t-norm of type of aggregation
control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"),
                    type.relation = c("tolerance", "eq.1"))
IND.1 <- BC.IND.relation.FRST(decision.table, attributes = attributes,
                              control = control.ind)

#### Calculate fuzzy indiscernibility relation: gaussian kernel ####
control.ind <- list(type.aggregation = c("t.tnorm", "t.cos"),
                    type.relation = c("transitive.kernel", "gaussian", 0.2))
IND.2 <- BC.IND.relation.FRST(decision.table, attributes = attributes,
                              control = control.ind)

#### calculate fuzzy transitive closure ####
control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"),
                    type.relation = c("transitive.closure", "eq.1"))
IND.3 <- BC.IND.relation.FRST(decision.table, attributes = attributes,
                              control = control.ind)

#### Calculate fuzzy indiscernibility relation: using user-defined relation #####
## The customized function should have three arguments which are : decision.table
## and object x, and y.
## This following example shows that we have an equation for similarity equation:
## 1 - abs(x - y) where x and y are two objects that will be compared.
## In this case, we do not consider decision.table in calculation.
FUN.relation <- function(decision.table, x, y) {
           return(1 - (abs(x - y)))
       }
control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"),
                     type.relation = c("custom", FUN.relation))
IND.4 <- BC.IND.relation.FRST(decision.table, attributes = attributes,
                              control = control.ind)
#### The other example using "custom"
## In this case, we calculate aggregation as average of all objects
FUN.average <- function(data){
  	 return(Reduce("+", data)/length(data))
}
control.ind <- list(type.aggregation = c("custom", FUN.average),
                     type.relation = c("tolerance", "eq.1"))
IND.5 <- BC.IND.relation.FRST(decision.table, attributes = attributes,
                              control = control.ind)

#### Calculate fuzzy indiscernibility relation: using user-defined tnorms ####
FUN.tnorm <- function(left.side, right.side) {
               if ((left.side + right.side) > 1)
                   return(min(left.side, right.side))
               else return(0)}
control.ind <- list(type.aggregation = c("t.tnorm", FUN.tnorm),
                    type.relation = c("tolerance", "eq.1"))
IND.6 <- BC.IND.relation.FRST(decision.table, attributes = attributes,
                              control = control.ind)
##################################################################
## Example 3: Dataset containing continuous and nominal attributes
## Note. we only consider type.relation = c("tolerance", "eq.1")
## but other approaches have the same way.
##################################################################
data(RoughSetData)
decision.table <- RoughSetData$housing7.dt

## in this case, we only consider the attribute: 1, 2, 3, 4
attributes <- c(1,2,3,4)

#### Calculate fuzzy indiscernibility relation ####
control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"),
                    type.relation = c("tolerance", "eq.1"))
IND.6 <- BC.IND.relation.FRST(decision.table, attributes = attributes, control = control.ind)
}
\references{
  H. Naessens, H. De Meyer, and B. De Baets, "Algorithms
  for the Computation of T-Transitive Closures", IEEE
  Trans. on Fuzzy Systems, vol. 10, No. 4, p. 541 - 551
  (2002).

  R. Jensen and Q. Shen, "New Approaches to Fuzzy-Rough
  Feature Selection", IEEE Trans. on Fuzzy Systems, vol.
  19, no. 4, p. 824 - 838 (2009).

  Q. Hu, D. Yu, W. Pedrycz, and D. Chen, "Kernelized Fuzzy
  Rough Sets and Their Applications", IEEE Trans. Knowledge
  Data Eng., vol. 23, no. 11, p. 1649 - 1471 (2011).
}
\seealso{
  \code{\link{BC.LU.approximation.FRST}},
  \code{\link{BC.IND.relation.RST}},
  \code{\link{BC.LU.approximation.RST}},

  and \code{\link{BC.positive.reg.FRST}}
}

