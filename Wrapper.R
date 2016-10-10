#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)
source("Utilities.R")
if (length(args) < 3) {stop("At least 3 arguments must be supplied (2 input files and a reference)", call.=FALSE)}
File1 = args[1]
File2 = args[2]
databaseFile = args[3]
OutputDir = ifelse(length(args) == 3, getwd(), args[4])
print(paste("The arguments were:", paste(File1, File2, databaseFile, OutputDir, collapse = ",")))
myAlign = parseNex(databaseFile)
prunedTypes = pruneTypes(myAlign, threshold = THRESH)
redAlign = prunedTypes[[1]]
clusterMap = prunedTypes[[2]]
curFiles = c(File1, File2)
myOutput = wrapperKnownNew(curFiles, redAlign, N = K, graph = GRAPH, unique = UNIQ, split = SPLIT, hidden = FALSE, 
                maxFrac = MMF, cMap = clusterMap, outDir = OutputDir)
foundFractions = myOutput[[2]]
initDir = getwd()
setwd(OutputDir)
dir.create("Results")
setwd("Results")
dir.create("Fractions")
setwd("Fractions")
Filename = paste0(longestCommonPrefix(File1, File2), "Fractions.csv")
write.table(foundFractions, file = Filename, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ",")
setwd(initDir)