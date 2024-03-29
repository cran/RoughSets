% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BasicFuzzyRoughSets.R
\name{BC.discernibility.mat.FRST}
\alias{BC.discernibility.mat.FRST}
\title{The decision-relative discernibility matrix based on fuzzy rough set theory}
\usage{
BC.discernibility.mat.FRST(
  decision.table,
  type.discernibility = "standard.red",
  control = list()
)
}
\arguments{
\item{decision.table}{a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}. 
It should be noted that this case only supports the nominal/symbolic decision attribute.}

\item{type.discernibility}{a string representing a type of discernibility. See in Section \code{Details}.}

\item{control}{a list of other parameters.
       \itemize{
       \item \code{type.relation}: a type of fuzzy indiscernibility relation. The default value is \code{type.relation = c("tolerance", "eq.1")}.

                See \code{\link{BC.IND.relation.FRST}}.
       \item \code{type.aggregation}: a type of aggregation operator. The default value is \code{type.aggregation = c("t.tnorm", "lukasiewicz")}.

               See \code{\link{BC.IND.relation.FRST}}.
       \item \code{t.implicator}: a type of implicator operator. The default value is \code{"lukasiewicz"}.

               See \code{\link{BC.LU.approximation.FRST}}.
       \item \code{type.LU}: a type of method of lower and upper approximations. The default value is \code{"implicator.tnorm"}.
       
               See \code{\link{BC.LU.approximation.FRST}}.
       \item \code{alpha.precision}: a numeric value representing a precision variable. It is used when using \code{"alpha.red"} as \code{type.discernibility}.
             The default value is 0.05.

                See \code{\link{BC.LU.approximation.FRST}}.  
       \item \code{show.discernibilityMatrix}: a boolean value determining whether the discernibility matrix will be shown or not (NULL). The default value is \code{FALSE}.
       \item \code{epsilon}: a numeric between 0 and 1 representing the \eqn{\epsilon} value on 

            \code{type.discernibility = "gaussian.red"}. It should be noted that when having nominal values on all attributes then \eqn{\epsilon} should be 0. 
            The default value is 0.
       \item \code{delta}: a numeric representing the \eqn{\delta} value on \code{"gaussian"} equations 

            (see \code{\link{BC.IND.relation.FRST}}). The default value is 2.
          \item \code{range.object}: a vector representing considered objects to construct the \code{k}-relative discernibility matrix. 
               The default value is \code{NULL} which means that we are using all objects in the decision table.
       }}
}
\value{
A class \code{"DiscernibilityMatrix"} containing the following components: 
\itemize{
\item \code{disc.mat}: a matrix showing the decision-relative discernibility matrix \eqn{M(\mathcal{A})} 
       which contains \eqn{n \times n} where \eqn{n} is the number of objects. It will be printed when choosing 

       \code{show.discernibilityMatrix = TRUE}.
\item \code{disc.list}: the decision-relative discernibility represented in a list.
\item \code{discernibility.type}: a string showing the chosen type of discernibility methods.
\item \code{type.model}: in this case, it is \code{"FRST"}.
}
}
\description{
This is a function that is used to build the decision-relative discernibility matrix based on FRST.
It is a matrix whose elements contain discernible attributes among pairs of objects. 
By means of this matrix, we are able to produce all decision reducts of the given decision system.
}
\details{
In this function, we provide several approaches in order to generate the decision-relative discernibility matrix. 
Theoretically, all reducts are found by constructing 
the matrix that contains elements showing discernible attributes among objects. 
The discernible attributes are determined by a specific condition which depends on the selected algorithm. A particular approach can be executed by selecting
a value of the parameter \code{type.discernibility}. The following shows the different 
values of the parameter \code{type.discernibility} corresponding approaches considered in this function.
\itemize{
\item \code{"standard.red"}: It is adopted from (Tsang et al, 2008)'s approach. 
The concept has been explained briefly in \code{\link{Introduction-FuzzyRoughSets}}. 
In order to use this algorithm, we assign the \code{control} parameter
       with the following components:

       \code{control = list(type.aggregation, type.relation, type.LU, t.implicator)}

       The detailed description of the components can be seen in \code{\link{BC.IND.relation.FRST}} and 

       \code{\link{BC.LU.approximation.FRST}}. Furthermore, in this case the authors suggest to use the "min" t-norm  
      (i.e., \code{type.aggregation = c("t.tnorm", "min")}) and the implicator operator "kleene_dienes" (i.e., \code{t.implicator = "kleene_dienes"}).
  
\item \code{"alpha.red"}: It is based on (Zhao et al, 2009)'s approach where all reductions will 
      be found by building an \eqn{\alpha}-discernibility matrix. This matrix contains elements which are defined by

      1) if \eqn{x_i} and \eqn{x_j} belong to different decision concept,

      \eqn{c_{ij} = \{R : \mathcal{T}(R(x_i, x_j), \lambda) \le \alpha \}},

      where \eqn{\lambda = (R_{\alpha} \downarrow A)(u)} which is lower approximation 
      of FVPRS (See \code{\link{BC.LU.approximation.FRST}}). 

      2) \eqn{c_{ij}={\oslash}}, otherwise.

      To generate the discernibility matrix based on this approach, we use the \code{control} parameter
      with the following components:

       \code{control = list(type.aggregation, type.relation, t.implicator, alpha.precision)} 

       where the lower approximation \eqn{\lambda} is fixed to \code{type.LU = "fvprs"}. The detailed
       description of the components can be seen in \code{\link{BC.IND.relation.FRST}} and \code{\link{BC.LU.approximation.FRST}}.
       Furthermore, in this case the authors suggest to use \eqn{\mathcal{T}}-similarity relation 

       e.g., \code{type.relation = c("transitive.closure", "eq.3")},

       the "lukasiewicz" t-norm  (i.e., \code{type.aggregation = c("t.tnorm", "lukasiewicz")}), and \code{alpha.precision} from 0 to 0.5.
   
\item \code{"gaussian.red"}: It is based on (Chen et al, 2011)'s approach. The discernibility matrix contains elements which are defined by: 
       
       1) if \eqn{x_i} and \eqn{x_j} belong to different decision concept, 

        \eqn{c_{ij}= \{R : R(x_i, x_j) \le \sqrt{1 - \lambda^2(x_i)}\}},

       where \eqn{\lambda = inf_{u \in U}\mathcal{I}_{cos}(R(x_i, u), A(u)) - \epsilon}. To generate fuzzy relation \eqn{R} , we use the fixed parameters as follows:

       \code{t.tnorm = "t.cos"} and \code{type.relation = c("transitive.kernel", "gaussian")}. 
       
       2) \eqn{c_{ij}={\oslash}}, otherwise.

       In this case, we need to define \code{control} parameter as follows.

       \code{control <- list(epsilon)}

      It should be noted that when having nominal values on all attributes then \code{epsilon} (\eqn{\epsilon}) should be 0. 

\item \code{"min.element"}: It is based on (Chen et al, 2012)'s approach where we only consider finding 
      the minimal element of the discernibility matrix by introducing the binary relation \eqn{DIS(R)} the relative discernibility relation 
     of conditional attribute \eqn{R} with respect to decision attribute \eqn{d}, which is computed as
     
     \eqn{DIS(R) = \{(x_i, x_j) \in U \times U: 1 - R(x_i, x_j) > \lambda_i, x_j \notin [x_i]_d\}},

     where \eqn{\lambda_i = (Sim(R) \downarrow [x_i]_d)(x_i)} with \eqn{Sim(R)} a fuzzy equivalence relation. 
     In other words, this algorithm does not need to build the discernibility matrix. 
     To generate the fuzzy relation \eqn{R} and lower approximation \eqn{\lambda}, we use the \code{control} parameter
       with the following components:

       \code{control = list(type.aggregation, type.relation, type.LU, t.implicator)}. 

       The detailed description of the components can be seen in \code{\link{BC.IND.relation.FRST}} and 

       \code{\link{BC.LU.approximation.FRST}}.
}
}
\examples{
#######################################################################
## Example 1: Constructing the decision-relative discernibility matrix
## In this case, we are using The simple Pima dataset containing 7 rows. 
#######################################################################
data(RoughSetData)
decision.table <- RoughSetData$pima7.dt 

## using "standard.red"
control.1 <- list(type.relation = c("tolerance", "eq.1"), 
                type.aggregation = c("t.tnorm", "min"), 
                t.implicator = "kleene_dienes", type.LU = "implicator.tnorm")
res.1 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "standard.red", 
                                    control = control.1)

## using "gaussian.red"
control.2 <- list(epsilon = 0)
res.2 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "gaussian.red",
                                    control = control.2)

## using "alpha.red"
control.3 <- list(type.relation = c("tolerance", "eq.1"), 
                type.aggregation = c("t.tnorm", "min"),
                t.implicator = "lukasiewicz", alpha.precision = 0.05)
res.3 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "alpha.red", 
                                    control = control.3)

## using "min.element"
control.4 <- list(type.relation = c("tolerance", "eq.1"), 
                type.aggregation = c("t.tnorm", "lukasiewicz"),
                t.implicator = "lukasiewicz", type.LU = "implicator.tnorm")
res.4 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "min.element", 
                                    control = control.4)

#######################################################################
## Example 2: Constructing the decision-relative discernibility matrix
## In this case, we are using the Hiring dataset containing nominal values
#######################################################################
\dontrun{data(RoughSetData)
decision.table <- RoughSetData$hiring.dt 

control.1 <- list(type.relation = c("crisp"), 
                type.aggregation = c("crisp"),
                t.implicator = "lukasiewicz", type.LU = "implicator.tnorm")
res.1 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "standard.red", 
                                    control = control.1)

control.2 <- list(epsilon = 0)
res.2 <- BC.discernibility.mat.FRST(decision.table, type.discernibility = "gaussian.red",
                                    control = control.2)}
}
\references{
D. Chen, L. Zhang, S. Zhao, Q. Hu, and P. Zhu, "A Novel Algorithm for Finding Reducts 
with Fuzzy Rough Sets", IEEE Trans. on Fuzzy Systems, vol. 20, no. 2, p. 385 - 389 (2012). 

D. G. Chen, Q. H. Hu, and Y. P. Yang, "Parameterized Attribute Reduction with
Gaussian Kernel Based Fuzzy Rough Sets", Information Sciences, vol. 181, no. 23, 
p. 5169 - 5179 (2011).

E. C. C. Tsang, D. G. Chen, D. S. Yeung, X. Z. Wang, and J. W. T. Lee, 
"Attributes Reduction Using Fuzzy Rough Sets", IEEE Trans. Fuzzy Syst., 
vol. 16, no. 5, p. 1130 - 1141 (2008).


S. Zhao, E. C. C. Tsang, and D. Chen, "The Model of Fuzzy Variable Precision Rough Sets",
IEEE Trans. on Fuzzy Systems, vol. 17, no. 2, p. 451 - 467 (2009).
}
\seealso{
\code{\link{BC.discernibility.mat.RST}}, \code{\link{BC.LU.approximation.RST}}, and \code{\link{BC.LU.approximation.FRST}}
}
\author{
Lala Septem Riza
}
