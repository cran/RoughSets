#############################################################################
#
#  This file is a part of the R package "RoughSets".
#
#  Author: Lala Septem Riza and Andrzej Janusz
#  Supervisors: Chris Cornelis, Francisco Herrera, Dominik Slezak and Jose Manuel Benitez
#  Copyright (c):
#       DiCITS Lab, Sci2s group, DECSAI, University of Granada and
#       Institute of Mathematics, University of Warsaw
#
#  This package is free software: you can redistribute it and/or modify it under
#  the terms of the GNU General Public License as published by the Free Software
#  Foundation, either version 2 of the License, or (at your option) any later version.
#
#  This package is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
#  A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
#############################################################################
#' This part contains global explanations about the implementation and use of the \code{RoughSets} package. 
#' The package \code{RoughSets} attempts to provide a complete tool to model and analyze 
#' information systems based on rough set theory (RST) and fuzzy rough set theory (FRST). 
#' From fundamental point of view, this package allows to construct rough sets by defining lower and upper approximations. 
#' Furthermore, recent methods for tackling common tasks in data mining such as preprocessing processes (e.g. discretization, feature selection, 
#' and instance selection), rule induction, and predicting classes or decision values of new datasets 
#' are available in this package as well. 
#'
#' There are two main parts considered in this package which are RST and FRST. 
#' RST was introduced by (Z. Pawlak, 1982; Z. Pawlak, 1991) which provides sophisticated mathematical tools to model 
#' and analyze information systems that involve uncertainty and imprecision. By employing indiscernibility relation among objects, RST does not require 
#' additional parameters to extract information. 
#' The detailed explanation about the fundamental concepts of RST can be read in Section \code{\link{A.Introduction-RoughSets}}. Secondly, FRST, an extension of
#' RST, was introduced by (D. Dubois and H. Prade, 1990) as a combination between 
#' fuzzy sets proposed by (L. A. Zadeh, 1965) and RST. This concept allows to analyze continuous attributes without 
#' performing discretization on data first. The basic concepts of FRST
#' can be seen in \code{\link{B.Introduction-FuzzyRoughSets}}. 
#'
#' Based on the above concepts, many methods have been proposed and applied for dealing with several different domains. 
#' In order to solve the problems, the methods employ the indiscernibility relation and lower and upper approximation concepts. 
#' All methods that have been implemented in this package will be explained by grouping based on their domains. The following is
#' a list of domains considered in this package:
#' \itemize{
#' \item Basic concepts of RST and FRST: This part, we can divide into four different tasks which are
#'       indiscernibility relation, lower and upper approximation, positive region and discernibility matrix.
#'       All of those tasks have been explained briefly in Section \code{\link{A.Introduction-RoughSets}} and
#'      \code{\link{B.Introduction-FuzzyRoughSets}}.
#' \item Discretization: It is used to convert real valued attributes into nominal/symbolic ones in an information system. 
#'       In RST point of view, this task attempts to maintain the discernibility between objects. 
#' \item Feature selection: It is a process for finding a subset of features which have the same quality as the complete feature set. 
#'      In other words, its purpose is to select the significant features and eliminate the dispensible ones.
#'      It is a useful and necessary process when we are facing datasets containing large numbers of features. From RST and FRST perspective, 
#'      feature selection refers to searching superreducts and reducts. The detailed information about reducts can be read in 
#'      \code{\link{A.Introduction-RoughSets}} and \code{\link{B.Introduction-FuzzyRoughSets}}.
#' \item Instance selection: This process is aimed to remove noisy, superfluous, or inconsistent instances from training datasets but retain consistent ones. 
#'      Therefore, good accuracy of classification is achieved by removing instances which do not give positive contributions.  
#' \item Prediction/classification: This task is used to predict decision values of a new dataset (test data). 
#'      We consider implementing some methods to perform this task such as fuzzy-rough nearest neighbours approaches, etc.
#' \item Rule induction: This task refers to generate IF - THEN rules. The rule represents knowledge which is contained in a dataset. 
#'      One advantage of building rules is that naturally the model is easy to interpret. Then, predicted values can be determined by
#'      considering the rules.
#' }
#'
#' As we mentioned before, we have embedded many well-known algorithms or techniques for handling the above domains. The algorithms were considered
#' since experimentally it has been proven that they were able to tackle complex tasks. They are implemented as functions that were organized 
#' to work with the same data structures. So, users can perform various approaches for a particular task easily and then compare their results. 
#' In order to be recognized quickly, generally we have chosen the names of the functions with some conventions. The names contain three parts 
#' which are \code{prefix}, \code{suffix}, and \code{middle} names that are separated by a point. The following is a description of each
#' part.  
#' \itemize{
#' \item \code{prefix}: There are some different prefixes for names of functions expressing a kind of task to be performed. 
#'                      The function names with prefix \code{BC} refer to \emph{basic concepts} which means that the functions are created for 
#'                      implementing the basic concepts of RST and FRST. For instance, the function \code{\link{BC.IND.relation.RST}} 
#'                      is used to calculate the indiscernibility relation which is one of the basic concepts of RST. 
#'                      Generally, functions having the prefix \code{BC} are called by many other functions. 
#'                      While prefix \code{D} refers to \emph{discretization}, \code{FS}, \code{IS}, \code{RI}, and \code{C} refer to \emph{feature selection}, 
#'                      \emph{instance selection}, \emph{rule induction} and \emph{classifier} domains. Furthermore, \code{SF} and \code{X} mean that
#'                      functions are used as \emph{supporting functions} which are not related directly with RST and FRST and \emph{auxiliary} functions which are called as a parameter.
#' \item \code{suffix}: It is located at the end of names. There are two types available which are \code{RST} and \code{FRST}. \code{RST} represents \emph{rough set theory} 
#'                      while \code{FRST} shows that the function is applied to \emph{fuzzy rough set theory}. Additionally, some functions that do not have 
#'                      \code{RST} or \code{FRST} suffix are used for both theories.
#' \item \code{middle}: All other words in the middle of the names are used to express the name of a particular method/algorithm or functionality. 
#'                      In this case, it could consist of more than one word separated by points.
#' }
#' Other functions that have names not based on the above rules are S3 functions e.g. \code{summary} and \code{predict} which are 
#' used to summarize objects and predict new data, respectively. 
#' 
#' The \code{RoughSets} package contains two main groups which are implementations of algorithms based on rough set and fuzzy rough set theories. 
#' For each part, we have considered many algorithms/methods to be implemented. The complete description about the algorithms and their associated functions is 
#' illustrated briefly as follows.   
#' \enumerate{
#' \item \bold{The implementations of RST}: This part outlines some considered algorihtms/methods based on RST. 
#'              The approaches can be classified based on their tasks as follows:
#'             \enumerate{
#'             \item The basic concepts of RST: The following is a list showing tasks and their implementation as functions.
#'                       \itemize{
#'                       \item Indiscernibility relation: It is a relation determining whether two objects are indiscernible by some attributes. 
#'                             It has been implemented in \code{\link{BC.IND.relation.RST}}. 
#'                       \item Lower and upper approximations: These approximations show whether objects can classified with certainty or not. 
#'                             It has been implemented in \code{\link{BC.LU.approximation.RST}}. 
#'                       \item Positive region: It is used to determine objects that are included in positive region and the corresponding degree of dependency.
#'                             It has been implemented in \code{\link{BC.positive.reg.RST}}.   
#'                       \item Discernibility matrix: It is used to create discernibility matrix showing attributes that discern each pair of objects. 
#'                             It has been implemented in \code{\link{BC.discernibility.mat.RST}}.  
#'                       }
#'             \item Discretization: There are several methods that have been considered in the package. The methods are implemented in the following functions.  
#'                      \itemize{
#'                       \item \code{\link{D.max.discernibility.matrix.RST}}: It implements the maximal discernibility algorithm to choose the cut values
#'                              which discern the largest number of pairs of objects in the decision-relative discernibility matrix. 
#'                       \item \code{\link{D.local.discernibility.matrix.RST}}: It implements the local strategy 
#'                              which implements decision tree to calculate the quality of a cut (i.e. number of objects discerned by cut).
#'                       \item \code{\link{D.global.discernibility.heuristic.RST}}: It implements the global discernibility algorithm 
#'                              which is computing globally semi-optimal cuts using the maximum discernibility heuristic.
#'						 \item \code{\link{D.discretize.quantiles.RST}}: It is a function used for computing cuts of the "quantile-based" discretization into \eqn{n} intervals. 	
#'                       \item \code{\link{D.discretize.equal.intervals.RST}}: It is a function used for computing cuts of the "equal interval size" discretization into \eqn{n} intervals.					
#'                       }
#'      		The output of these methods is a list of cut values which are the values for converting real to nominal values. 
#'              So, we provide the function \code{\link{SF.applyDecTable}} that is used to generate a new decision table according to the cut values we got by discretization methods. 
#'              Additionally, we have implemented \code{\link{D.discretization.RST}} as a wrapper function collecting all methods considered to perform discretization tasks.
#'             \item Feature selection: The output of this task can be classified into the following three groups as follows.
#'                      \itemize{
#'                       \item Feature subset: It refers to a superreduct which is not necessarily minimal. In other words, the methods in this group 
#'                              might generate just a subset of attributes.
#'                             \itemize{
#'                                  \item QuickReduct algorithm: It has been implemented in \code{\link{FS.quickreduct.RST}}. 
#'                                  \item Superreduct generation based on some criteria: It is based on different criteria which are 
#'                                                    entropy, gini index, discernibility measure, size of positive region. 
#'                                                    It has been implemented in \code{\link{FS.greedy.heuristic.superreduct.RST}}.                                                
#'                             }
#'                             Furthermore, we provide a wrapper function \code{\link{FS.feature.subset.computation}} in order to give a user interface for many methods of RST and FRST. 
#'                       \item Reduct: The following are methods that produce a single decision reduct. 
#'                             \itemize{
#'                                  \item Reduct generation based on some criteria: It is based on different criteria which are 
#'                                                    entropy, gini index, discernibility measure, size of positive region. 
#'                                                    It has been implemented in \code{\link{FS.greedy.heuristic.reduct.RST}}.
#'                                  \item Permutation reduct: It is based on a permutation schema over all attributes. 
#'                                                     It has been implemented in \code{\link{FS.permutation.heuristic.reduct.RST}}.
#'                              }
#'                             Furthermore, we provide a wrapper function \code{\link{FS.reduct.computation}} in order to give a user interface toward many methods of RST and FRST. 
#'                       \item All reducts: In order to get all decision reducts, first we execute the \code{\link{BC.discernibility.mat.RST}} function for
#'                             constructing a decision-relative discernibility matrix. After obtaining the matrix, \code{\link{FS.all.reducts.computation}}
#'                             is called to get all reducts.   
#'                       }
#'                       The output of the above methods is a class containing a decision reduct/feature subset and other descriptions. 
#'                       For generating a new decision table according to the decision reduct, we provide the function \code{\link{SF.applyDecTable}}. 
#'            \item Rule induction: We have provided the function \code{\link{RI.indiscernibilityBasedRules.RST}} to generate rules. This function requires the output of 
#'                       feature selection functions for getting a superreduct/reduct. After obtaining the rules, 
#'                       in the predicting process, we execute \code{\link{predict.RuleSetRST}} considering our rules and given newdata/testing data.
#'              }
#' \item \bold{The implementations of FRST}: As in the \code{RST} part, this group contains several algorithms that can be classified into several groups based on their purpose.
#'           The following is a description of all methods that have been implemented in functions.
#'  \enumerate{
#'	  \item Basic concepts of FRST: The following is a list showing tasks and their implementation as functions.
#'    \itemize{
#'       \item Indiscernibility relations: the are fuzzy relations determining to which degree two objects are similar. 
#'             This package provides several types of relations which are implemented in a single function 
#'             called \code{\link{BC.IND.relation.FRST}}. In this package, we consider several types of relations e.g. 
#'             fuzzy equivalence, tolerance, and \eqn{T}-similarity relations. These relations can be chosen by
#'             assigning \code{type.relation}. Additionally, In this function, we provide several options to 
#'             calculate aggregation e.g. triangular norm operator (e.g. \code{"lukasiewicz"}, \code{"min"}, etc)
#'             and user-defined operator. 
#' 	   	 \item Lower and upper approximations: These approximations show to what extent objects can be classified with certainty or not. 
#'           This task has been implemented in 
#'
#'           \code{\link{BC.LU.approximation.FRST}}. There are many approaches available in this package that can be selected by assigning the parameter \code{type.LU}. 
#'           The considered methods are
#'           implication/t-norm, \eqn{\beta}-precision fuzzy rough sets (\eqn{\beta}-PFRS), vaquely quantified rough sets (VQRS), fuzzy variable precision rough sets (FVPRS), ordered weighted average (OWA),
#'           soft fuzzy rough sets (SFRS), and robust fuzzy rough sets (RFRS). Furthermore, we provide a facility which is \code{"custom"} where users can create their own approximations by 
#'           defining functions to calculate lower and upper approximations. Many options to calculate implicator and triangular norm are also available.
#'      \item Positive region: It is used to determine the membership degree of each object to the positive region and the corresponding degree of dependency. 
#'                   It has been implemented in \code{\link{BC.positive.reg.FRST}}. 
#'      \item Discernibility matrix: It is used to construct the decision-relative discernibility matrix. There are some approaches to construct the matrix,
#'                   e.g. based on standard approach, gaussian reduction, alpha reduction, and minimal element in discernibility matrix. They have been implemented 
#'                   in \code{\link{BC.discernibility.mat.FRST}}. 
#'    }
#'    \item Feature selection: Generally, for dealing with this problem we have considered two different kinds of approaches 
#'       which are methods employing heuristics to produce a near-optimal reduction 
#'       such as the fuzzy QuickReduct algorithm and approaches based on 
#'       the decision-relative discernibility matrix
#'       (see \code{\link{BC.discernibility.mat.FRST}}). In the context of their outputs, 
#'       one of the differences between them is on the type of produced reduct which are a single superreduct or a subset of features, 
#'       a reduct and all reducts. The following is a description of all above types.
#'       \itemize{
#'           \item Feature subset: It refers to methods which produce a superreduct which is not necessarily a reduct. In other words methods in this group 
#'                               might generate just a subset of attributes.
#'                 The following is a complete list of methods considered in this package. 
#'                 \itemize{
#'                  \item positive region based algorithms: It refers to 
#'                        positive regions, as a way to evaluate attributes to be selected. 
#'                        Furthermore, we provide several different measures based on the positive region in this function.  
#'                        All methods included in this part employ the QuickReduct algorithm to obtain selected features.
#'                        In order to choose a particular algorithm, we assign parameter \code{type.method} in \code{\link{FS.quickreduct.FRST}}. 
#'                        The following is a description of all values of the parameter \code{type.method}. 
#'                        \itemize{
#'                          \item \code{"beta.pfrs"}: It is \eqn{\beta}-PFRS-based feature selection which employs the degree of 
#'                                  the dependency as a criterion to select features and it uses fuzzy lower approximations based on  
#'                                  \eqn{\beta}-precision fuzzy rough sets. 
#'                          \item \code{"vqrs"}: It employs the degree of the dependency which is obtained by using
#'                                  VQRS. 
#'                          \item \code{"fuzzy.dependency"}: It refers to feature selection based on degree of dependency (\eqn{\gamma}) as
#'                                 a value to select features. 
#'                          \item fuzzy discernibility function approach: It applies the decision-relative discernibility matrix to the QuickReduct algorithm. 
#'                          \item \code{"fvprs"}: It employs the degree of 
#'                                  the dependency based on FVPRS. 
#'                          \item \code{"min.positive.reg"}: it refers to feature selection with fuzzy decision reducts (\eqn{\delta}) which considers the fuzzy positive region
#'                                 to select features. 
#'                          \item \code{"owa"}: It uses the dependency degree of approximations based on OWA to select features.
#'                          \item \code{"rfrs"}: It uses the dependency degree of
#'                           		approximations based on RFRS.
#'                            }
#'                \item boundary region based algorithm: This algorithm is based on the membership degree to the fuzzy boundary region which is determined by subtracting 
#'                      the values of fuzzy upper and lower approximations.
#'                      This algorithm has been implemented in \code{\link{FS.quickreduct.FRST}} 
#'                      by setting the parameter as \code{type.method = "fuzzy.boundary.reg"}.
#'               } 
#'               Furthermore, we provide a wrapper function \code{\link{FS.feature.subset.computation}} in order to give a user interface for many methods of RST and FRST. 
#'           \item Reduct: The following methods produce a single decision reduct.
#'                \itemize{
#'                 \item The near-optimal reduction based algorithm: 
#'                       It is an algorihtm proposed by Zhao et al, which results 
#'                       a near-optimal reduction by modifying the discernibility matrix. 
#'                       Additionally, this algorithm uses fuzzy variable precision rough sets (FVPRS) for
#'                       calculating the lower approximation.
#'                       It has been implemented in \code{\link{FS.nearOpt.fvprs.FRST}}. 
#'                }
#'                Furthermore, we provide a wrapper function \code{\link{FS.reduct.computation}} in order to provide a user interface toward many methods of RST and FRST. 
#'           \item All reducts:  In order to get all decision reducts, first we execute the \code{\link{BC.discernibility.mat.FRST}} function for
#'                             constructing a decision-relative discernibility matrix. After obtaining the matrix, \code{\link{FS.all.reducts.computation}}
#'                             is called to get all reducts.    
#'     }
#'     The output of the above methods is a class containing a decision reduct/feature subset and other descriptions. 
#'     For generating a new decision table according to the decision reduct, we provide the function \code{\link{SF.applyDecTable}}. 
#'
#'     \item Rule induction: As we mentioned before, rule induction is a task used to generate 
#'           rules representing knowledge of a decision table. Commonly, this process is called learning phase in machine learning.
#'           The following methods are considered to generate rules:
#'           \itemize{
#' 	    		 \item \code{RI.hybridFS.FRST}: It combines fuzzy-rough rule induction
#'                      and feature selection. See \code{\link{RI.hybridFS.FRST}}.
#'               \item \code{RI.GFRS.FRST}: It refers to rule induction based on generalized fuzzy rough sets (GFRS). See \code{\link{RI.GFRS.FRST}}.
#'           }
#'           After generating rules, we can use them to predict decision values of new data 
#'           by executing a predicting function which is \code{\link{predict.RuleSetFRST}}.
#'     \item Instance selection: The following functions select instances to improve accuracy by
#'           removing noisy, superflous or inconsistent ones from training datasets.
#'           \itemize{
#'             \item \code{IS.FRIS.FRST}: It was proposed by Jensen and Cornelis. It evaluates the degree of membership to the positive region of each instance.
#'             		  If an instance's membership degree is less than the threshold, then the instance can be removed.      
#'             		  See \code{\link{IS.FRIS.FRST}}.
#'             \item \code{IS.FRPS.FRST}: The fuzzy-rough prototype selection (FRPS) based on Verbiest, et al's method. See \code{\link{IS.FRPS.FRST}}.
#'           }
#'       We provide the function \code{\link{SF.applyDecTable}} that is used to generate a new decision table according to the output of instance selection functions. 
#'    \item Fuzzy-rough nearest neighbour approaches: This part provides nearest neighbour based methods for 
#'          predicting decision values/classes of new datasets. In other words, by supplying a decision table as training data
#'          we can predict decision values of new data at the same time.
#'          We have considered the following methods:
#'          \itemize{
#'            \item \code{C.FRNN.FRST}: The fuzzy-rough nearest neighbors based on Jensen and Cornelis' technique. 
#'                  See \code{\link{C.FRNN.FRST}}.
#'            \item \code{C.FRNN.O.FRST}: The fuzzy-rough ownership nearest neighbour algorithm based on Sarkar's method. See \code{\link{C.FRNN.O.FRST}}.
#'                  \item \code{C.POSNN.FRST}: The positive region based fuzzy-rough nearest neighbour algorithm based on Verbiest et al's technique. See \code{\link{C.POSNN.FRST}}.
#'          }
#' }
#' }
#' To get started with the package, the user can have a look at the examples included in
#' the documentation on each function. Additionally, to show general usage of the package briefly, 
#' we also provide some examples in this section. 
#' 
#' If you have problems using the package, find a bug, or have suggestions, 
#' please contact the package maintainer by email, instead of writing to the general R lists 
#' or to other internet forums and mailing lists.
#' 
#' Furthermore, there are many demos that ship with the package. To get a list of them, type:
#' 
#' \code{demo()}
#' 
#' Then, to start a demo, type \code{demo(<demo_name_here>)}. All the demos are present as 
#' R scripts in the package sources in the "demo" subdirectory.
#' 
#' Currently, there are the following demos available:
#' 
#' \itemize{
#' \item Basic concepts of RST and FRST:
#'
#' \code{demo(BasicConcept.RST)},
#' \code{demo(BasicConcept.FRST)},
#'
#' \code{demo(DiscernibilityMatrix.RST)},
#' \code{demo(DiscernibilityMatrix.FRST)}.
#'
#' \item Discretization based on RST:
#'
#' \code{demo(D.local.discernibility.matrix.RST)},
#' \code{demo(D.max.discernibility.matrix.RST)},
#'
#' \code{demo(D.global.discernibility.heuristic.RST)},
#' \code{demo(D.discretize.equal.intervals.RST)}, 
#'
#' \code{demo(D.discretize.quantiles.RST)}.
#'
#' \item Feature selection based on RST:
#'
#' \code{demo(FS.permutation.heuristic.reduct.RST)},
#' \code{demo(FS.greedy.heuristic.reduct.RST)},
#'
#' \code{demo(FS.greedy.heuristic.reduct.RST)},
#' \code{demo(FS.quickreduct.RST)}.
#'
#' \item Feature selection based on FRST:
#'
#' \code{demo(FS.QuickReduct.FRST.Ex1)},
#' \code{demo(FS.QuickReduct.FRST.Ex2)},
#'
#' \code{demo(FS.QuickReduct.FRST.Ex3)},
#' \code{demo(FS.QuickReduct.FRST.Ex4)},
#'
#' \code{demo(FS.QuickReduct.FRST.Ex5)},
#' \code{demo(FS.nearOpt.fvprs.FRST)}.
#' 
#' \item Instance selection based on FRST:
#'
#' \code{demo(IS.FRIS.FRST)}, 
#' \code{demo(IS.FRPS.FRST)}
#' 
#' \item Classification using the Iris dataset:
#' 
#' \code{demo(FRNN.O.iris)},
#' \code{demo(POSNN.iris)},
#' \code{demo(FRNN.iris)}.
#'
#' \item Rule induction based on RST:
#'
#' \code{demo(RI.indiscernibilityBasedRules.RST)}.
#'
#' \item Rule induction based on FRST:
#'
#' \code{demo(RI.classification.FRST)},
#' \code{demo(RI.regression.FRST)}.
#' }
#'
#' Some decision tables have been embedded in this package which can be seen in 
#' \code{\link{RoughSetData}}. 
#'
#' Finally, you may visit the package webpage \url{http://sci2s.ugr.es/dicits/software/RoughSets}, 
#' where we provide a more extensive introduction as well as additional explanations of 
#' the procedures.
#' 
#' @name RoughSets-package
#' @aliases RoughSets
#' @docType package
#' @title Getting started with the RoughSets package
#' @references 
#' D. Dubois and H. Prade, "Rough Fuzzy Sets and Fuzzy Rough Sets",
#' International Journal of General Systems, vol. 17, p. 91 - 209 (1990).
#'
#' L.A. Zadeh, "Fuzzy Sets",
#' Information and Control, vol. 8, p. 338 - 353 (1965).
#'
#' Z. Pawlak, "Rough Sets", 
#' International Journal of Computer and Information System, 
#' vol. 11, no. 5, p. 341 - 356 (1982).
#' 
#' Z. Pawlak, "Rough Sets: Theoretical Aspects of Reasoning About Data, System Theory, Knowledge Engineering and Problem Solving",
#' vol. 9, Kluwer Academic Publishers, Dordrecht, Netherlands (1991). 
#'
# @keywords rough sets fuzzy rough sets instance selection feature selection classification prediction
#' @author Lala Septem Riza \email{lala.s.riza@@decsai.ugr.es},
#' 
#' Andrzej Janusz \email{andrzejanusz@@gmail.com}, 
#' 
#' Chris Cornelis \email{chriscornelis@@decsai.ugr.es}, 
#'
#' Francisco Herrera \email{herrera@@decsai.ugr.es}, 
#'
#' Dominik Slezak \email{slezak@@mimuw.edu.pl},
#' 
#' and Jose Manuel Benitez \email{j.m.benitez@@decsai.ugr.es}
#' 
#' DiCITS Lab, SCI2S group, CITIC-UGR, DECSAI, University of Granada,
#' 
#' \url{http://dicits.ugr.es}, \url{http://sci2s.ugr.es}
#'
#' Institute of Mathematics, University of Warsaw.
#' 
#' @useDynLib RoughSets 
#' @examples
#' ##############################################################
#' ## A.1 Example: Basic concepts of rough set theory
#' ##############################################################
#' ## Using hiring data set, see RoughSetData doc
#' data(RoughSetData)
#' decision.table <- RoughSetData$hiring.dt 
#'
#' ## define considered attributes which are first, second and 
#' ## third attributes
#' attr.P <- c(1,2,3)
#' 
#' ####### compute indiscernibility relation #######
#' IND <- BC.IND.relation.RST(decision.table, attribute = attr.P)
#'
#' ####### compute lower and upper approximations #####
#' ## Let us define fourth index as the decision attribute
#' decision.attr <- c(4)
#' roughset <- BC.LU.approximation.RST(decision.table, IND, decision.attr)
#' 
#' ####### Determine regions ######
#' region.RST <- BC.positive.reg.RST(decision.table, roughset)
#'
#' ####### The decision-relative discernibility matrix and reduct #####
#' disc.mat <- BC.discernibility.mat.RST(decision.table, range.object = NULL)
#'
#' ##############################################################
#' ## A.2 Example: Basic concepts of fuzzy rough set theory
#' ##############################################################
#' ## Using pima7 data set, see RoughSetData doc
#' data(RoughSetData)
#' decision.table <- RoughSetData$pima7.dt 
#' 
#' ## In this case, let us consider the first and second attributes
#' conditional.attr <- c(1, 2)
#'
#' ## We are using "lukasiewicz" t-norm and "tolerance" relation 
#' ## with "eq.1" as fuzzy similarity equation
#' control.ind <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'                     type.relation = c("tolerance", "eq.1"))
#'
#' ## Compute fuzzy indiscernibility relation 
#' IND.condAttr <- BC.IND.relation.FRST(decision.table, attributes = conditional.attr, 
#'                             control = control.ind) 
#'
#' ## Compute fuzzy lower and upper approximation using type.LU : "implicator.tnorm" 
#' ## Define index of decision attribute
#' decision.attr = c(9)
#'
#' ## Compute fuzzy indiscernibility relation of decision attribute
#' ## We are using "crisp" for type of aggregation and type of relation
#' control.dec <- list(type.aggregation = c("crisp"), type.relation = "crisp")
#'
#' IND.decAttr <- BC.IND.relation.FRST(decision.table, attributes = decision.attr, 
#'                             control = control.dec) 
#' 
#' ## Define control parameter containing type of implicator and t-norm
#' control <- list(t.implicator = "lukasiewicz", t.tnorm = "lukasiewicz")
#'
#' ## Compute fuzzy lower and upper approximation
#' FRST.LU <- BC.LU.approximation.FRST(decision.table, IND.condAttr, IND.decAttr, 
#'               type.LU = "implicator.tnorm", control = control)
#' 
#' ## Determine fuzzy positive region and its degree of dependency
#' fuzzy.region <- BC.positive.reg.FRST(decision.table, FRST.LU)
#'
#' ###############################################################
#' ## B Example : Data analysis based on RST and FRST 
#' ## In this example, we are using wine dataset for both RST and FRST
#' ###############################################################
#' ## Using wine data set, see RoughSetData doc
#' \dontrun{data(RoughSetData)
#' wine.decTable <- SF.asDecisionTable(dataset = RoughSetData$wine.dt[6 : 178, ], 
#'                                     decision.attr = 14, indx.nominal = 14)
#' 
#' ## define first five instances as newdata or testing data
#' tst.wine <- SF.asDecisionTable(dataset = wine.decTable[1:5, -ncol(wine.decTable)])
#' 
#' ###############################################################
#' ## B.1 Example : Rough Set Theory 
#' ###############################################################
#' #####  DISCRETIZATION STEP 
#' ## In this example, we are using local strategy algorithm
#' cut.values.tra <- D.discretization.RST(wine.decTable, 
#'                   type.method = "global.discernibility") 
#'
#' ## generate new decision table
#' d.new.tra.rst <- SF.applyDecTable(wine.decTable, cut.values.tra)
#' d.new.tst.rst <- SF.applyDecTable(tst.wine, cut.values.tra)
#'
#' ##### FEATURE SELECTION STEP 
#' ## In this example, we are using permutation algorithm
#' ## which generates a single superreduct
#' red.rst <- FS.feature.subset.computation(d.new.tra.rst, 
#'                   method = "quickreduct.rst")
#' 
#' ## generate new decision table according to the reduct (optional)
#' fs.new.decTable.rst <- SF.applyDecTable(d.new.tra.rst, red.rst)
#'
#' ##### RULE INDUCTION 
#' ## In this case, we are using the original decision table,
#' ## we also can use the decision table resulting from feature selection						 
#' rules.rst <- RI.indiscernibilityBasedRules.RST(d.new.tra.rst, red.rst)
#'					
#' ## predicting newdata
#' res.1 <- predict(rules.rst, d.new.tst.rst)
#'
#' ###############################################################
#' ## B.2 Example : Fuzzy Rough Set Theory
#' ###############################################################
#' ## FEATURE SELECTION STEP
#' red.frst <- FS.feature.subset.computation(wine.decTable, 
#'                   method = "quickreduct.frst")
#'
#' fs.new.decTable.frst <- SF.applyDecTable(decision.table = wine.decTable, red.frst)
#'
#' ## INSTANCE SELECTION STEP 
#' ## It should be noted in this case we are using the decision table
#' ## that results from feature selection
#' is.object <- IS.FRIS.FRST(decision.table = fs.new.decTable.frst, control = 
#'                        list(threshold.tau = 0.4, alpha = 1))
#' is.new.decTable.frst <- SF.applyDecTable(decision.table = fs.new.decTable.frst, is.object)
#'
#' ## CLASSIFIER
#' ## Predict the testing data by using the FRNN.O method
#' control <- list(m = 2, type.membership = "gradual")
#' new.tst.wine <- SF.applyDecTable(decision.table = tst.wine, red.frst)
#'
#' ## Compute predictions
#' pred.class <- C.FRNN.O.FRST(decision.table = is.new.decTable.frst, newdata = new.tst.wine, 
#'                           control = control)
#'
#' ## RULE INDUCTION 
#' control <- list(type.aggregation = c("t.tnorm", "lukasiewicz"), 
#'            type.relation = c("tolerance", "eq.3"), t.implicator = "lukasiewicz")
#' rules.hybrid <- RI.hybridFS.FRST(decision.table = is.new.decTable.frst, control)
#'
#' ## predicting newdata
#' res.1 <- predict(rules.hybrid, tst.wine)}
NULL


