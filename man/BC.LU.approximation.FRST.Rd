% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BasicFuzzyRoughSets.R
\name{BC.LU.approximation.FRST}
\alias{BC.LU.approximation.FRST}
\title{The fuzzy lower and upper approximations based on fuzzy rough set theory}
\usage{
BC.LU.approximation.FRST(
  decision.table,
  IND.condAttr,
  IND.decAttr,
  type.LU = "implicator.tnorm",
  control = list()
)
}
\arguments{
\item{decision.table}{a \code{"DecisionTable"} class representing the decision table. See \code{\link{SF.asDecisionTable}}.}

\item{IND.condAttr}{a \code{"IndiscernibilityRelation"} class of the conditional attributes which is produced by \code{\link{BC.IND.relation.FRST}}.}

\item{IND.decAttr}{a \code{"IndiscernibilityRelation"} class of the decision attribute which is produced by \code{\link{BC.IND.relation.FRST}}.}

\item{type.LU}{a string representing a chosen method to calculate lower and upper approximations. See the explanation in Section \code{Details}.}

\item{control}{a list of other parameters. In order to understand how to express the \code{control} parameter, 
        see the explanation in Section \code{Details}. 
        The descriptions of those components and their values is as follows. 
        \itemize{ 
             \item \code{t.tnorm}: a type of triangular functions which have been explained 

                  in \code{\link{BC.IND.relation.FRST}}.
               \item \code{t.implicator}: a type of implicator functions.  
                   The following are values of this parameter:               
                        \itemize{
                        \item \code{"kleene_dienes"} means \eqn{max(1 - x_1, x_2)}.
                        \item \code{"lukasiewicz"} means \eqn{min(1 - x_1 + x_2, 1)}. It is the default value. 
                        \item \code{"zadeh"} means \eqn{max(1 - x_1, min(x_1, x_2))}.
                        \item \code{"gaines"} means \eqn{(x_1 <= x_2 ? 1 : x_2 / x_1)}.
                        \item \code{"godel"} means \eqn{(x_1 <= x_2 ? 1 : x_2)}.
                        \item \code{"kleene_dienes_lukasiewicz"} means \eqn{1 - x_1 + x_1 * x_2}.
                        \item \code{"mizumoto"} means \eqn{(1 - x_1 + x_1 * x_2)}.
                        \item \code{"dubois_prade"} means \eqn{(x_2 == 0 ? 1 - x_1 : (x_1 == 1 ? x_2 : 1))}.
                        }
                 Where we consider the following rule: \eqn{x_1 -> x_2}. 
               \item \code{q.some}: a vector of alpha and beta parameters of vaguely quantified rough set  
                 for quantifier \code{some}. The default value is \code{q.some = c(0.1, 0.6)}.
               \item \code{q.most}: a vector of alpha and beta parameters of vaguely quantified rough set 
                 for quantifier \code{most}. The default value is \code{q.most = c(0.2, 1)}.
               \item \code{alpha}: a numeric between 0 and 1 representing the threshold parameter of the fuzzy variable precision rough sets 
                  (FVPRS) (see Section \code{Details}). The default value is 0.05.
             \item \code{m.owa}: an integer number (\eqn{m}) which is used in the OWA fuzzy rough sets (see Section \code{Details}). 
                   
                  The default value is \code{m.owa = round(0.5 * ncol(decision.table))}.
             \item \code{w.owa}: a vector representing the weight vector in the OWA fuzzy rough sets (see Section \code{Details}).
                   The default value is \code{NULL}, which means we use the \code{m.owa} type.
             \item \code{type.rfrs}: a type of robust fuzzy rough sets which is one of the following methods:
                   \code{"k.trimmed.min"}, \code{"k.mean.min"}, \code{"k.median.min"}, \code{"k.trimmed.max"},
                   \code{"k.mean.max"}, and \code{"k.median.max"} (see Section \code{Details}). The default value is \code{"k.trimmed.min"}.
             \item \code{k.rfrs}: a number between 0 and the number of data which is used to define considered data on 
                   robust fuzzy rough sets (RFRS) (see Section \code{Details}). The default value is 
                   \code{k.rfrs = round(0.5*nrow(decision.table))}.
             \item \code{beta.quasi}: a number between 0 and 1 representing \eqn{\beta}-precision t-norms and t-conorms in \eqn{\beta}-PFRS.
                   The default value is 0.05.
       }}
}
\value{
A class \code{"LowerUpperApproximation"} representing fuzzy rough set (fuzzy lower and upper approximations). It contains the following components:
        \itemize{
         \item \code{fuzzy.lower}: a list showing the lower approximation classified 
               based on decision concepts for each index of objects. The value refers to
               the degree of objects included in the lower approximation.  
               In case the decision attribute is continuous, the result is in a data frame 
               with dimension (number of objects x number of objects) and the value on position \eqn{(i,j)} 
               shows the membership of object \eqn{i} to the lower approximation of the similarity class of object \eqn{j}.
         \item \code{fuzzy.upper}: a list showing the upper approximation classified 
               based on decision concepts for each index of objects. The value refers to
               the degree of objects included in the upper approximation. 
               In case the decision attribute is continuous values, the result is in data frame 
               with dimension (number of objects x number of objects) and the value on position \eqn{(i,j)} 
               shows the membership of object \eqn{i} to the upper approximation of the similarity class of object \eqn{j}.
         \item \code{type.LU}: a string representing the type of lower and upper approximation approaches.
         \item \code{type.model}: a string showing the type of model which is used. In this case, it is \code{"FRST"} which means fuzzy rough set theory.
         }
}
\description{
This is a function implementing a fundamental concept of FRST: fuzzy lower and upper approximations. 
Many options have been considered for determining lower and upper approximations, 
such as techniques based on implicator and t-norm functions proposed by 
(Radzikowska and Kerre, 2002).
}
\details{
Fuzzy lower and upper approximations as explained in \code{\link{Introduction-FuzzyRoughSets}} are used
to define to what extent the set of elements can be classified into a certain class strongly or weakly. We can perform various methods by choosing the parameter \code{type.LU}. 
The following is a list of all \code{type.LU} values: 
   \itemize{
   \item \code{"implicator.tnorm"}: It means implicator/t-norm based model proposed by (Radzikowska and Kerre, 2002). 
         The explanation has been given in \code{\link{Introduction-FuzzyRoughSets}}.
         Other parameters in \code{control} related with this approach are \code{t.tnorm} and \code{t.implicator}.
         In other words, when we are using \code{"implicator.tnorm"} as \code{type.LU}, 
         we should consider parameters \code{t.tnorm} and \code{t.implicator}.
         The possible values of these parameters can be seen in the description of parameters. 
   \item \code{"vqrs"}: It means vaguely quantified rough sets proposed by 
         (Cornelis et al, 2007). Basically, this concept proposed to replace fuzzy lower and upper approximations 
         based on Radzikowska and Kerre's technique (see \code{\link{Introduction-FuzzyRoughSets}})
         with the following equations, respectively. 

         \eqn{(R_{Q_u} \downarrow A)(y) = Q_u(\frac{|R_y \cap A|}{|R_y|})} 

         \eqn{(R_{Q_l} \uparrow A)(y) = Q_l(\frac{|R_y \cap A|}{|R_y|})}

         where the quantifier \eqn{Q_u} and \eqn{Q_l} represent the terms \code{most} and \code{some}. 

   \item \code{"owa"}: It refers to ordered weighted average based fuzzy rough sets. 
           This method was introduced by (Cornelis et al, 2010) and computes the approximations
           by an aggregation process proposed by (Yager, 1988). The OWA-based lower and upper approximations of 
           \eqn{A} under \eqn{R} with weight vectors \eqn{W_l} and \eqn{W_u} are defined as

          \eqn{(R \downarrow W_l A)(y) = OW A_{W_l}\langle I(R(x, y), A(y))\rangle}

          \eqn{(R \uparrow W_u A)(y) = OW A_{W_u}\langle T(R(x, y), A(y))\rangle}

          We provide two ways to define the weight vectors as follows:

          \itemize{
          \item \code{m.owa}: Let \eqn{m.owa} be \eqn{m} and \eqn{m \le n}, this model is defined by

          \eqn{W_l = <w_i^l> = w_{n+1-i}^l = \frac{2^{m-i}}{2^{m}-1}} for \eqn{i = 1,\ldots, m} and \eqn{0} for \eqn{i = m + 1, \ldots, n}

          \eqn{W_u = <w_i^u> = w_i^u = \frac{2^{m - i}}{2^{m} - 1}} for \eqn{i = 1, \ldots, m} and \eqn{0} for \eqn{i = m + 1, \ldots, n}

          where \eqn{n} is the number of data.

          \item \code{custom}: In this case, users define the own weight vector.
          It should be noted that the weight vectors \eqn{<w_i>} should satisfy \eqn{w_i \in [0, 1]} and 
          their summation is 1.
          }

   \item \code{"fvprs"}: It refers to fuzzy variable precision rough sets (FVPRS) introduced by 
           (Zhao et al, 2009). It is a combination between variable precision rough sets (VPRS)
            and FRST. This function implements the construction of lower and upper approximations as follows.
      
		  \eqn{(R_{\alpha} \downarrow A)(y) = inf_{A(y) \le \alpha} \mathcal{I}(R(x,y), \alpha) \wedge inf_{A(y) > \alpha} \mathcal{I}(R(x,y), A(y))}

          \eqn{(R_{\alpha} \uparrow A)(y) = sup_{A(y) \ge N(\alpha)} \mathcal{T}(R(x,y), N(\alpha)) \vee sup_{A(y) < N(\alpha)} \mathcal{T}(R(x,y), A(y))}
            
          where \eqn{\alpha}, \eqn{\mathcal{I}} and \eqn{\mathcal{T}} are the variable precision parameter, implicator and t-norm operators, respectively.
   \item \code{"rfrs"}: It refers to robust fuzzy rough sets (RFRS) proposed by (Hu et al, 2012). 
           This package provides six types of RFRS which are k-trimmed minimum, k-mean minimum, k-median minimum, 
           k-trimmed maximum, k-mean maximum, and k-median maximum. 
           Basically, these methods are a special case of ordered weighted average (OWA) where they consider 
           the weight vectors as follows.
           \itemize{
           \item \code{"k.trimmed.min"}: \eqn{w_i^l = 1} for \eqn{i = n - k} and \eqn{w_i^l = 0} otherwise.
           \item \code{"k.mean.min"}: \eqn{w_i^l = 1/k} for \eqn{i > n - k} and \eqn{w_i^l = 0} otherwise.
           \item \code{"k.median.min"}: \eqn{w_i^l = 1} if k odd, \eqn{i = n - (k-1)/2} and \eqn{w_i^l = 1/2} if k even, \eqn{i = n - k/2} 
                                and \eqn{w_i^l = 0} otherwise.
           \item \code{"k.trimmed.max"}: \eqn{w_i^u = 1} for \eqn{i = k + 1} and \eqn{w_i^u = 0} otherwise. 
           \item \code{"k.mean.max"}: \eqn{w_i^u = 1/k} for \eqn{i < k + 1} and \eqn{w_i^u = 0} otherwise.   
           \item \code{"k.median.max"}: \eqn{w_i^u = 1} if k odd, \eqn{i = (k + 1)/2} and \eqn{w_i^u = 1/2} if k even, \eqn{i = k/2 + 1} 
                                or \eqn{w_i^u = 0} otherwise.
           }

   \item \code{"beta.pfrs"}: It refers to \eqn{\beta}-precision fuzzy rough sets (\eqn{\beta}-PFRS) proposed by 
           (Salido and Murakami, 2003). This algorithm uses \eqn{\beta}-precision quasi-\eqn{\mathcal{T}}-norm and 
           \eqn{\beta}-precision quasi-\eqn{\mathcal{T}}-conorm. The following are the \eqn{\beta}-precision versions of fuzzy lower and upper approximations of a fuzzy set \eqn{A} in \eqn{U}

				\eqn{(R_B \downarrow A)(y) = T_{\beta_{x \in U}} \mathcal{I}(R_B(x,y), A(x))}

				\eqn{(R_B \uparrow A)(y) = S_{\beta_{x \in U}} \mathcal{T}(R_B(x,y), A(x))} 

           where \eqn{T_{\beta}} and \eqn{S_{\beta}} are \eqn{\beta}-precision quasi-\eqn{\mathcal{T}}-norm and \eqn{\beta}-precision quasi-\eqn{\mathcal{T}}-conorm.
           Given a t-norm \eqn{\mathcal{T}}, a t-conorm \eqn{S}, \eqn{\beta \in [0,1]} and \eqn{n \in N \setminus \{0, 1\}}, the corresponding 
           \eqn{\beta}-precision quasi-t-norm \eqn{T_{\beta}} and \eqn{\beta}-precision-\eqn{\mathcal{T}}-conorm \eqn{S_{\beta}} of order \eqn{n} are
           \eqn{[0,1]^n \to [0,1]} mappings such that for all \eqn{x = (x_1,...,x_n)} in \eqn{[0,1]^n},

            \eqn{T_{\beta}(x) = \mathcal{T}(y_1,...,y_{n-m})},

            \eqn{S_{\beta}(x) = \mathcal{T}(z_1,...,z_{n-p})},

            where \eqn{y_i} is the \eqn{i^{th}} greatest element of \eqn{x} and \eqn{z_i} is the \eqn{i^{th}} smallest element of \eqn{x}, and

            \eqn{m = max(i \in \{0,...,n\}|i \le (1-\beta)\sum_{j=1}^{n}x_j)},

            \eqn{p = max(i \in \{0,...,n\}|i \le (1-\beta)\sum_{j=1}^{n}(a - x_j))}.

           In this package we use \code{min} and \code{max} for \eqn{\mathcal{T}}-norm and \eqn{\mathcal{T}}-conorm, respectively. 

   \item \code{"custom"}: It refers to user-defined lower and upper approximations. An example can be seen in Section \code{Examples}.
   }
The parameter \code{type.LU}, which is explained above, is related with parameter \code{control}. 
In other words, when choosing a specific value of \code{type.LU}, we should take into account to set values of related components in \code{control}.
The components that are considered depend on what kind of lower and upper approximations are used. 
So, we do not need to assign all components for a particular approach but only components related with \code{type.LU}.
The following is a list showing the components of each approaches.
               \itemize{
               \item \code{type.LU = "implicator.tnorm"}: 

                \code{control <- list(t.implicator, t.tnorm)}

                \item \code{type.LU = "vqrs"}:

                \code{control <- list(q.some, q.most, type.aggregation, t.tnorm)}

                \item \code{type.LU = "owa"}:

                \code{control <- list(t.implicator, t.tnorm, m.owa)} 
                
                or

               \code{control <- list(t.implicator, t.tnorm, w.owa)} 

                \item \code{type.LU = "fvprs"}:

                \code{control <- list(t.implicator, t.tnorm, alpha)}

               \item \code{type.LU = "beta.pfrs"}: 

                \code{control <- list(t.implicator, t.tnorm, beta.quasi)}

                \item \code{type.LU = "rfrs"}:

                \code{control <- list(t.implicator, t.tnorm, type.rfrs, k.rfrs)}

                \item \code{type.LU = "custom"}:

                \code{control <- list(t.implicator, t.tnorm, FUN.lower, FUN.upper)}
               }
The description of the components can be seen in the \code{control} parameter.
In Section \code{Examples}, we provide two examples showing different cases which are
when we have to handle a nominal decision attribute and a continuous one. 

It should be noted that this function depends on \code{\link{BC.IND.relation.FRST}}
which is a function used to calculate the fuzzy indiscernibility relation as input data. 
So, it is obvious that before performing this function, users must execute \code{\link{BC.IND.relation.FRST}} first.
}
\examples{
###########################################################
## 1. Example: Decision table contains nominal decision attribute
## we are using the same dataset and indiscernibility 
## relation along this example.
###########################################################
dt.ex1 <- data.frame(c(-0.4, -0.4, -0.3, 0.3, 0.2, 0.2), 
                     c(-0.3, 0.2, -0.4, -0.3, -0.3, 0),
			        c(-0.5, -0.1, -0.3, 0, 0, 0),
			        c("no", "yes", "no", "yes", "yes", "no"))
colnames(dt.ex1) <- c("a", "b", "c", "d")
decision.table <- SF.asDecisionTable(dataset = dt.ex1, decision.attr = 4)

## let us consider the first and second attributes 
## only as conditional attributes
condAttr <- c(1, 2)

## let us consider the fourth attribute as decision attribute
decAttr <- c(4)

#### calculate fuzzy indiscernibility relation ####
control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
                    type.relation = c("tolerance", "eq.1"))
control.dec <- list(type.aggregation = c("crisp"), type.relation = "crisp")

## fuzzy indiscernibility relation of conditional attribute
IND.condAttr <- BC.IND.relation.FRST(decision.table, attributes = condAttr, 
                                     control = control.ind)

## fuzzy indiscernibility relation of decision attribute 
IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = decAttr, 
                                     control = control.dec)

#### Calculate fuzzy lower and upper approximation using type.LU : "implicator.tnorm" ####	 
control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz")
FRST.LU <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
              type.LU = "implicator.tnorm", control = control)

#### Calculate fuzzy lower and upper approximation using type.LU : "vqrs" ####	 
control <- list(q.some = c(0.1, 0.6), q.most = c(0.2, 1), t.tnorm = "lukasiewicz")
FRST.VQRS <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
              type.LU = "vqrs", control = control)

#### Calculate fuzzy lower and upper approximation using type.LU : "owa" ####	 
control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz", m.owa = 3) 
FRST.OWA.1 <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
              type.LU = "owa", control = control)

#### Calculate fuzzy lower and upper approximation using type.LU : 
#### "owa" with customized function 
#### In this case, we are using the same weight vector as
#### previous one with m.owa = 3
control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz", 
               w.owa =  c(0, 0, 0, 0.14, 0.29, 0.57)) 
FRST.OWA.2 <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
              type.LU = "owa", control = control)

#### Calculate fuzzy lower and upper approximation using type.LU : "fvprs" ####	 
control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz", alpha = 0.05)
FRST.fvprs <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
              type.LU = "fvprs", control = control)


#### Calculate fuzzy lower and upper approximation using type.LU : "rfrs" ####	 
control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz", 
                type.rfrs = "k.trimmed.min", k.rfrs = 0)
FRST.rfrs <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
              type.LU = "rfrs", control = control)

#### Calculate fuzzy lower and upper approximation using type.LU : "beta.pfrs" ####	 
control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz", beta.quasi = 0.05)
FRST.beta.pfrs <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
              type.LU = "beta.pfrs", control = control)

#### Calculate fuzzy lower and upper approximation using type.LU : "custom" ####	
## In this case, we calculate approximations randomly. 
f.lower <- function(x){
        return(min(runif(1, min = 0, max = 1) * x))	
}
f.upper <- function(x){
        return(max(runif(1, min = 0, max = 1) * x))
}
control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz", FUN.lower = f.lower, 
                FUN.upper = f.upper)
FRST.custom <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
              type.LU = "custom", control = control)

#### In this case, we use custom function for triangular norm and implicator operator
## For example, let us define our implicator and t-norm operator as follows.
imp.lower <- function(antecedent, consequent){
                 return(max(1 - antecedent, consequent))
              }
tnorm.upper <- function(x, y){
                return (x * y)
             } 
control.custom <- list(t.implicator = imp.lower, t.tnorm = tnorm.upper)
FRST.LU.custom <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
              type.LU = "implicator.tnorm", control = control.custom)

###########################################################
## 2. Example: Decision table contains a continuous decision attribute.
## It should be noted that in this example, we are using
## the same dataset and indiscernibility relation.
## We only show one method but for other approaches 
## the procedure is analogous to the previous example
###########################################################
## In this case, we are using housing dataset containing 7 objects
data(RoughSetData)
decision.table <- RoughSetData$housing7.dt

## let us consider the first and second conditional attributes only,
## and the decision attribute at 14.
cond.attributes <- c(1, 2)
dec.attributes <- c(14)
control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
               type.relation = c("tolerance", "eq.1"))
IND.condAttr <- BC.IND.relation.FRST(decision.table, attributes = cond.attributes, 
                                     control = control.ind) 
IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = dec.attributes, 
                                    control = control.ind) 

#### Calculate fuzzy lower and upper approximation using type.LU : "implicator.tnorm" ####	 
control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz")
FRST.LU <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
              type.LU = "implicator.tnorm", control = control)

}
\references{
A. M. Radzikowska and E. E. Kerre, "A Comparative Study of Fuzzy Rough Sets", 
Fuzzy Sets and Systems, vol. 126, p. 137 - 156 (2002). 

C. Cornelis, M. De Cock, and A. Radzikowska, "Vaguely Quantified Rough Sets",
Proceedings of 11th International Conference on Rough Sets, Fuzzy Sets,
Data Mining and Granular Computing (RSFDGrC2007), Lecture Notes in
Artificial Intelligence 4482, p. 87 - 94 (2007).

C. Cornelis, N. Verbiest, and R. Jensen, "Ordered Weighted Average Based Fuzzy
Rough Sets", Proceedings of the 5th International Conference on Rough Sets
and Knowledge Technology (RSKT 2010), p. 78 - 85 (2010).


J. M. F. Salido and S. Murakami, "Rough Set Analysis of a General Type of Fuzzy Data
Using Transitive Aggregations of Fuzzy Similarity Relations", 
Fuzzy Sets Syst., vol. 139, p. 635 - 660 (2003).


Q. Hu, L. Zhang, S. An, D. Zhang, and D. Yu, "On Robust Fuzzy Rough Set Models",
IEEE Trans. on Fuzzy Systems, vol. 20, no. 4, p. 636 - 651 (2012).

R. Jensen and Q. Shen,  
"New Approaches to Fuzzy-Rough Feature Selection", 
IEEE Trans. on Fuzzy Systems, vol. 19, no. 4,
p. 824 - 838 (2009).

R. R. Yager, "On Ordered Weighted Averaging Aggregation Operators in Multicriteria
Decision Making", IEEE Transactions on Systems, Man, and Cybernetics, vol. 18, p. 183 - 190 (1988).

S. Y. Zhao, E. C. C. Tsang, and D. G. Chen, 
"The Model of Fuzzy Variable Precision Rough Sets",
IEEE Trans. Fuzzy Systems, vol. 17, no. 2,
p. 451 - 467 (2009).
}
\seealso{
\code{\link{BC.IND.relation.RST}}, \code{\link{BC.LU.approximation.RST}}, 
          and \code{\link{BC.positive.reg.FRST}}
}
\author{
Lala Septem Riza
}
