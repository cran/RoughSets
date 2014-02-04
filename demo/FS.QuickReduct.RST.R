 library(RoughSets)
 data(RoughSetData)
 decision.table <- RoughSetData$hiring.dt 

 ## evaluate single reduct
 res.1 <- FS.quickreduct.RST(decision.table)
 
 ## generate new decision table according to the reduct
 new.decTable <- SF.applyDecTable(decision.table, res.1)