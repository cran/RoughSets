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
## it is used to represent rules in order to get the standard format
# @param rules original form of rules
# @param type.method a type of method
# @param decision.table a decision table
# @param t.similarity a type of similarity equation
# @param t.tnorm a type of t-tnorm
std.rules <- function(rules, type.method = "RI.hybridFS.FRST", decision.table, t.similarity = NULL, t.tnorm = NULL){
	objects <- decision.table
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	decision.attr <- attr(decision.table, "decision.attr")
	consequent.attr <- names(desc.attrs)[decision.attr]
	
	if (any(type.method == c("RI.hybridFS.FRST", "RI.GFRS.FRST"))){
		all.attr = colnames(objects)
		
		if (nominal.att[length(nominal.att)] == TRUE){
			type.task <- c("classification")
		}
		else {
			type.task <- c("regression")
		}
		duplicated.pt <- duplicated(rules$attributes)
		new.obj <- c()
		ii <- 1
		for (i in 1 : length(duplicated.pt)){
			if (duplicated.pt[i] == FALSE){
				new.obj <- append(new.obj, list(rules$objects[[i]]))
				ii <- ii + 1
			}
			else {
				new.obj[[ii - 1]] <- c(new.obj[[ii - 1]], rules$objects[[i]])
			}
		}
		
		if (type.method == "RI.hybridFS.FRST"){
			rules <- list(attributes = rules$attributes[which(duplicated.pt == FALSE)], objects = new.obj)	
			attrs.rules <- unique(unlist(rules$attributes))				
			rules <- ch.rules(rules, type.method = type.method, decision.table)
			
			## calculate variance and range of data
			res.varRange <- cal.var.range(objects, c(which(colnames(objects) %in% attrs.rules)), nominal.att)
			variance.data <- res.varRange$variance.data
			rownames(variance.data) <- NULL
			colnames(variance.data) <- attrs.rules
			range.data <- res.varRange$range.data	
			colnames(range.data) <- attrs.rules
			
			mod <- list(rules = rules, type.model = "FRST", type.method = type.method, type.task = type.task, 
						t.similarity = t.similarity, t.tnorm = t.tnorm, variance.data = variance.data, range.data = range.data, 
						antecedent.attr = attrs.rules, consequent.attr = consequent.attr, nominal.att = nominal.att[c(which(colnames(objects) %in% attrs.rules))])			
			class.mod <- ObjectFactory(mod, classname = "RuleSetFRST")
			return(class.mod)
		}
		else if (type.method == "RI.GFRS.FRST"){
			rules <- list(attributes = rules$attributes[which(duplicated.pt == FALSE)], objects = new.obj, 
			            threshold = rules$threshold[sort(unlist(new.obj))])
			attrs.rules <- unique(unlist(rules$attributes))		
			rules <- ch.rules(rules, type.method = type.method, decision.table)

			## calculate variance and range of data
			res.varRange <- cal.var.range(objects, c(which(colnames(objects) %in% attrs.rules)), nominal.att)
			variance.data <- res.varRange$variance.data
			rownames(variance.data) <- NULL
			colnames(variance.data) <- attrs.rules
			range.data <- res.varRange$range.data	
			colnames(range.data) <- attrs.rules
			
			mod <- list(rules = rules, type.model = "FRST", type.method = type.method,
			            type.task = type.task, t.similarity = t.similarity, t.tnorm = t.tnorm, variance.data = variance.data, range.data = range.data, 
						antecedent.attr = attrs.rules, consequent.attr = consequent.attr, nominal.att = nominal.att[c(which(colnames(objects) %in% attrs.rules))])
			class.mod <- ObjectFactory(mod, classname = "RuleSetFRST")
			return(class.mod)
		}
	}
}

## it is used to change rules in order to get the standard format
# @param rules original form of rules
# @param type.method a type of method
# @param decision.table
ch.rules <- function(rules, type.method = "RI.hybridFS.FRST", decision.table){
	tra.objects <- decision.table
	desc.attrs <- attr(decision.table, "desc.attrs")
	nominal.att <- attr(decision.table, "nominal.attrs")
	rule <- c()
		
	for (i in 1 : length(rules$attributes)){
		exit <- FALSE
		ii <- 1
		while (exit == FALSE) {
			## get value for antecedent and consequent of each rule
			anteVal <- tra.objects[rules$objects[[i]][ii], c(rules$attributes[[i]]), drop = FALSE]
			conqVal <- tra.objects[rules$objects[[i]][ii], ncol(tra.objects), drop = FALSE]
			if (is.null(rule)){
				tmpRule <- cbind(anteVal, conqVal)
				rownames(tmpRule) <- "value"
				rule <- list(tmpRule)
			}
			else {
				tmpRule <- cbind(anteVal, conqVal)
				rownames(tmpRule) <- "value"
				rule <- append(rule, list(tmpRule))
			}
			
			## checking stopping criteria
			if (ii == length(rules$objects[[i]])){
				exit <- TRUE
			}
			else {
				ii <- ii + 1
			}
		}
	}
	
	if (type.method == "RI.GFRS.FRST"){
		return(list(rules = rule, threshold = rules$threshold))
	} 
	else {
		return(list(rules = rule))
	}
}

#an auxiliary function for constructing rules from indiscernibility classes
laplaceEstimate <- function(rule, dataS, clsVec, uniqueCls, suppIdx = NULL) {
	if (is.null(suppIdx)) {
		if (length(rule$idx) > 1) {
			dataS <- do.call(paste, dataS[rule$idx])
		}
		else {
			dataS <- as.character(dataS[[rule$idx]])
		}
	
		tmpValues = paste(rule$values, collapse = " ")
		suppIdx = which(dataS == tmpValues)
	}
  
	clsFreqs = table(clsVec[suppIdx])
	maxIdx = which.max(clsFreqs)
	nOfCorrect = clsFreqs[maxIdx]
  
	rule$consequent = names(clsFreqs)[maxIdx]
	rule$support = as.integer(suppIdx)
	rule$laplace = (nOfCorrect + 1)/(length(suppIdx) + length(uniqueCls))
	return(rule)
}

#' It is an auxiliary function of the best-first voting strategy.
#'
#' @title the best-first voting strategy function
#' @param object a class "RuleSetRST". 
#' @param ruleValues values of rules.
#' @param consequents values on the consequent part.
#' @param majorityCls  a value representing majority class of decision attribute.
#' @param ... other parameters.
#' @return predicted values. 
#' @export
X.bestFirst <- function(object, ruleValues, consequents, majorityCls, ...) { 
	matchIdx <- which(ruleValues == object)
	if (length(matchIdx) > 0) {
		prediction <- consequents[matchIdx[1]]
	}
	else {
		prediction <- majorityCls
	}
  
	return(prediction)
}
