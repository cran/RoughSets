 data(RoughSetData)
 decision.table <- RoughSetData$pima7.dt 
 
 ## get reduct
 reduct.2 <- FS.nearOpt.fvprs.FRST(decision.table)

 ## get new decision table according to the reduct
 new.decTable <- SF.applyDecTable(decision.table, reduct.2)