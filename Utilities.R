### This functon checks if the package is present, if it's not then it gets installed
checkAndRequire = function(pkgName) {
  if (!require(pkgName, character.only = TRUE)) {
    install.packages(pkgName, dependencies = TRUE, repos = 'http://cran.rstudio.com/')
  }
  ### Importing the package
  require(pkgName, character.only = TRUE)
}
### USEFUL PACKAGES ###
checkAndRequire("ape")
checkAndRequire("expm")
checkAndRequire("igraph")
### USEFUL CONSTANTS ###
ALPHABET = c("a","c","g","t")
MISSING = -2  # code for missing K-mers
START = -1    # code for starting K-mers
LEFT = -1     # code for left trimming
NONE = 0      # code for no trimming
RIGHT = 1     # code for right trimming
LOCAL <<- 1   # number of local alignments used
### CODE OPTIONS ###
MMF = 0.5     # maximum missing fraction for a path
THRESH = 0.05 # minimum pairwise distance
K = 25        # length of the K-mer to use for the library
SPLIT = TRUE  # TRUE if reads become non-overlapping K-mers
GRAPH = TRUE  # TRUE if using de Bruijn graphs to preprocess
UNIQ = TRUE   # TRUE if only unique K-mers from library used
HIDDEN = TRUE # TRUE if reference strains are hidden, one at a time
NDEC = 1      # number of decimals; used for rounding
NEW = TRUE    # TRUE if the new mixtures are being used
CORRECT=FALSE # TRUE if read pairing needs to be corrected

### This function is a wrapper function for all the relevant files found in the directory.
fullWrapper = function(databaseFile="../DatasetW1.nex.txt", maxFrac, threshold, K, split, graph, unique, hidden) {
  myAlign = parseNex(databaseFile)
	prunedTypes = pruneTypes(myAlign, threshold = threshold)
	redAlign = prunedTypes[[1]]
	clusterMap = prunedTypes[[2]]
	allFiles = list.files(pattern = ".fq")
	allFiles = gsub(".txt", "", allFiles)
	allFiles = gsub("R1.fq", "", allFiles)
	allFiles = gsub("R2.fq", "", allFiles)
	allFiles = unique(allFiles)
	allFiles = allFiles[grep("Hidden", allFiles, invert = TRUE)]
	Output = vector("list", length(allFiles))
	names(Output) = allFiles
	for (curFile in allFiles) {
	  print(curFile)
	  Output[[curFile]] = wrapperKnownNew(curFile, redAlign, K, graph, unique, split, hidden, maxFrac, clusterMap)
	}
	curDir = tail(unlist(strsplit(getwd(), "/")), 1)
	save(Output, file = paste0(curDir, "FullOutput", ifelse(hidden, "Hidden", "Known"), "K", K, ".RData"))
	opts1 = paste0("ResultsFrac", maxFrac, "Threshold", threshold, "K", K) 
	opts2 = paste0(rep("Split", split), rep("Graph", graph), rep("Unique", unique), rep("Hidden", hidden))
	outputFile = paste0(opts1, unlist(strsplit(databaseFile,"\\."))[1], opts2, ".csv")
	Table = formatOutput(Output, outputFile, clusterMap, hidden = hidden)
	Table
}

### This function takes the output of type frequency analysis, converts it into a table and writes it into a file
formatOutput = function(output, outputFile, clusterMap = NULL, hidden = HIDDEN, nDec = NDEC) {
	N = length(output)
	if (hidden) {
	  DimNames = list(NULL, c("Type", "True (%)", "Precision (%)", "Recall (%)"))
	}
	else {
	  DimNames = list(NULL, c("Type", "True (%)", "Predicted (%)", "Error (%)"))
	}
	Table = matrix(NA, nrow = 0, ncol = 4, dimnames = DimNames)
	for (ind in 1:N) {
	  curMat = processSingleResult(output[[ind]], clusterMap, hidden, nDec)
		Table = rbind(Table, curMat)
	}
	write.table(Table, file = outputFile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ",")
	Table
}

### This function prepares the output for the result of processing a single file 
processSingleResult = function(curPiece, clusterMap, hidden, nDec) {
  curTruth = curPiece[[1]]
  curPred = curPiece[[2]]
  curCombined = union(names(curTruth), names(curPred))
  curNames = unique(mapEntities(curCombined, clusterMap))
  curMat = matrix(0,  length(curNames), 3, dimnames = list(curNames, NULL))
  curMat[mapEntities(names(curTruth), clusterMap), 1] = round(curTruth * 100, nDec)
  if (hidden) {
    curMat[colnames(curPred), 2] = round(curPred[1, colnames(curPred)] * 100, nDec)
    curMat[colnames(curPred), 3] = round(curPred[2, colnames(curPred)] * 100, nDec)
  }
  else {
    curMat[mapEntities(names(curPred), clusterMap), 2] = round(curPred * 100, nDec)
    curMat[, 3] = round(abs(curMat[ ,1] - curMat[ ,2]), nDec)
    curError = sum(curMat[,3])/2
  }
  curMat = cbind(rownames(curMat), curMat)
  curMat = curMat[order(curMat[,1]), ]
  curMat = rbind(curMat, c(rep("", 3), ifelse(hidden, "", curError)))
  curMat
}

### This function returns all the K-mers of each row as well as the unique K-mers of each row of the input table
getUniqueKmers = function(Table, K) {
  N = nrow(Table)
  Names = rownames(Table)
  KmerList = vector("list", N)
  names(KmerList) = Names
  uniqueKmers = vector("list", N)
  names(uniqueKmers) = Names
  uniqueMapping = vector("list", N)
  names(uniqueMapping) = Names
  for (ind in 1:N) {
    curKmers = getKmers(Table[ind,], K, TRUE)
    KmerList[[ind]] = curKmers
  }
  for (ind in 1:N) {
    curList = KmerList[[ind]]
    stopifnot(!any(duplicated(curList)))
    otherLists = unlist(KmerList[setdiff(1:N, ind)])
    curUnique = which(!curList %in% otherLists)
    uniqueKmers[[ind]] = curList[curUnique]
    uniqueMapping[[ind]] = curUnique
  }
  output = list(KmerList, uniqueKmers, uniqueMapping)
  output
}

### This function creates the reverse complement of a string
reverseComplement = function(string) {
	inds = match(string, ALPHABET)
	revCompInds = 5 - rev(inds)
	revCompString = ALPHABET[revCompInds]
	revCompString
}

### This function maps the given entities (normally characters) using the given map
### The default function is applied to all the entities that do not have a mapping
mapEntities = function(Entities, Map, default = identity) {
	if (is.null(Map)) {
	  return(Entities)
	}
  else {
    mappedEntities = Map[Entities]
    bad = is.na(mappedEntities)
    mappedEntities[bad] = sapply(Entities[bad], default)
    mappedEntities = unlist(mappedEntities)
    return(mappedEntities)
  }
}

### This function constructs the de Bruijn graph to determine the strains that are present in a given file pair 
getGoodStrainsGraph = function (F1, F2, redAlign, N, maxFrac) {
  RCF1 = t(apply(F1, 1, reverseComplement))
  RCF2 = t(apply(F2, 1, reverseComplement))
  G0 = constructDBGraph(rbind(F1, F2, RCF1, RCF2), N)
  Res = identifyAllPaths(redAlign, G0[[1]], G0[[3]], maxMissingFrac = maxFrac, count = FALSE)
  goodStrains = names(Res)[!sapply(Res, is.null)]
  goodStrains
}

### This function processes all reads in the given tables, and provides the appropriate output
processAllReads = function(F1, F2, goodStrains, redAlign, N, unique, split, Name = NULL, cMap = NULL) {
  if (!is.null(Name)) {
    curStrains = setdiff(goodStrains, Name)
    curQ = getUniqueKmers(redAlign[curStrains, , drop = FALSE], N)
    curResult = mapPairedReadsNew(F1, F2, curQ, unique = unique, split = split, hidden = Name, cMap = cMap)
  }
  else {
    Qred = getUniqueKmers(redAlign[goodStrains, , drop = FALSE], N)
    curResult = mapPairedReadsNew(F1, F2, Qred, unique = unique, split = split, hidden = NULL, cMap = cMap)
  }
  curResult
}

### This function writes the specified subset of unmapped reads using baseFilename, Name, K, ind in the filename 
writeUnmappedReads = function(PF, unmappedReads, Filename) {
  Reads  = PF[[1]][unmappedReads, , drop = FALSE]
  Scores = PF[[2]][unmappedReads]
  writeF(Reads, Scores, Filename)
  return()
}

### This function tests if the file with given name is present in the directory, and extends it if it is not
completeFilename = function(filename, completion) {
  if (!filename %in% list.files()) {
    filename = paste0(filename, completion)
    stopifnot(filename %in% list.files())
  }
  filename
}

### This function finds the longest common prefix of two strings
longestCommonPrefix = function(string1, string2, removeDirs = TRUE) {
  Str1 = unlist(strsplit(string1, ""))
  Str2 = unlist(strsplit(string2, ""))
  L = min(length(Str1), length(Str2))
  maxPos = min(which(Str1[1:L] != Str2[1:L])) - 1
  Str = Str1[1:min(L, maxPos)]
  string = paste0(Str, collapse = "")
  if (removeDirs) {
    string = tail(unlist(strsplit(string, "/")), 1)
  }
  string
}

### This function is a new wrapper function that takes reverse complementarity of paired-end reads into account.
wrapperKnownNew = function(baseFilename, redAlign, N, graph, unique, split, hidden, maxFrac, cMap = NULL, outDir = NULL) {
  if (length(baseFilename) == 2) {
    File1 = baseFilename[1]
    File2 = baseFilename[2]
    correctFractions = NULL
    baseFilename = longestCommonPrefix(File1, File2, removeDirs = TRUE)
    DOUBLE = TRUE
  } 
  else {
    File1 = completeFilename(paste0(baseFilename, "R1.fq"), ".txt")
    File2 = completeFilename(paste0(baseFilename, "R2.fq"), ".txt")
    correctFractions = parseFilenameFracs(baseFilename)
    DOUBLE = FALSE
  }
  PF1 = parseF(File1, keepNames = TRUE, keepScores = (hidden || DOUBLE))
  PF2 = parseF(File2, keepNames = TRUE, keepScores = (hidden || DOUBLE))
  goodStrains = rownames(redAlign)
  if (graph) {
    goodStrains = getGoodStrainsGraph(PF1[[1]], PF2[[1]], redAlign, N, maxFrac)
    if (is.null(goodStrains)) {stop("Error: this should never happen!")}
  }
  if (hidden) {
    trueNames = mapEntities(names(correctFractions), cMap)
    estimatedFractions = matrix(NA, 2, length(trueNames), dimnames = list(c("Precision", "Recall"), trueNames))
    for (Name in trueNames) {
      curResult = processAllReads(PF1[[1]], PF2[[1]], goodStrains, redAlign, N, unique, split, Name, cMap = cMap)
      estimatedFractions[,Name] = unlist(curResult[[3]])
      Filename1 = paste0(baseFilename, "Hidden", Name, "K", N, ".R", 1, ".fq")
      Filename2 = paste0(baseFilename, "Hidden", Name, "K", N, ".R", 2, ".fq")
      writeUnmappedReads(PF1, curResult[[4]], Filename1)
      writeUnmappedReads(PF2, curResult[[4]], Filename2)
    }
  }
  else {
    curResult = processAllReads(PF1[[1]], PF2[[1]], goodStrains, redAlign, N, unique, split, Name = NULL, cMap = cMap)
    estimatedFractions = curResult[[2]]
    if (DOUBLE) {
      if (!is.null(outDir)) {
        if (!dir.exists(outDir)) {
          dir.create(outDir)
        }
        outDir = ifelse(substr(outDir, nchar(outDir), nchar(outDir)) == "/", outDir, paste0(outDir, "/"))
      }
      Filename1 = paste0(outDir, "LeftOver", baseFilename, 1, ".fq")
      Filename2 = paste0(outDir, "LeftOver", baseFilename, 2, ".fq")
      writeUnmappedReads(PF1, curResult[[4]], Filename1)
      writeUnmappedReads(PF2, curResult[[4]], Filename2)
    }
  }
  output = list(correctFractions, estimatedFractions)
  output
}

### This function parses the name of the file to determine the correct fractions of the strains present in it
parseFilenameFracs = function(Filename) {
	Filename = unlist(strsplit(Filename,"\\_"))[1]
	Filename = unlist(strsplit(Filename,"\\."))
	Strains = Filename[1]
	Fractions = as.numeric(Filename[-1])
	Fractions = Fractions[!is.na(Fractions)]
	Fractions = Fractions/sum(Fractions)
	splitPos = which(unlist(strsplit(Strains, "")) %in% LETTERS)
	numStrains = length(splitPos)
	splitStrains = substr(rep(Strains, numStrains), splitPos, c(splitPos[-1] - 1, nchar(Strains)))
	if (numStrains == 1) {
	  Fractions = 1
	}
	stopifnot(length(Fractions) == numStrains)
	names(Fractions) = splitStrains
	Fractions 
}

### This function constructs a de Bruijn graph for a table of reads or sequences with a specified K
### Returns the graph whose vertices are labelled by the K-mers, and a list of paths, one per sequence
constructDBGraph = function(Table, K) {
	R = nrow(Table)
	L = ncol(Table)
	allKmers = vector("list", R)
	Paths = vector("list", R)
	for (ind in 1:R) {
		allKmers[[ind]] = getKmers(Table[ind,], K, requireUnique = FALSE)
	}
	uniqueKmers = sort(unique(unlist(allKmers)))
	N = length(uniqueKmers)
	G = graph.empty(n = N)
	Map = 1:N
	names(Map) = uniqueKmers
	V(G)$name = uniqueKmers
	for (ind in 1:R) {
		curKmers = allKmers[[ind]]
		curL = length(curKmers)
		curEdges = cbind(curKmers[-curL], curKmers[-1])
		G = add.edges(G, curEdges)
		Paths[[ind]] = Map[curKmers]
	}
	output = list(G, Paths, Map)
	output
}

### This function traces the path corresponding to a sequence through a given input graph
### using K-mers of given length and a map to match the graph's vertices to these K-mers.
### If more than maxMissingFrac of the K-mers are absent, the return value is NULL; 
### otherwise, the K-mers are labelled, and if count = TRUE, they are also counted.
identifyPath = function(sequence, Graph, Map, maxMissingFrac, count) {
  K = nchar(V(Graph)$name[1])
	curKmers = getKmers(sequence, K, requireUnique = FALSE)
	numKmers = length(curKmers)
	mappedKmers = Map[curKmers]
	missing = is.na(mappedKmers)
	mappedKmers = mappedKmers[!missing]
	if (sum(missing)/numKmers < maxMissingFrac) {
		L = length(mappedKmers)
		checkEdges = rep(1, L - 1)
		if (count) {
		  for (ind in 1:(L-1)) {
		    checkEdges[ind] = length(E(Graph)[mappedKmers[ind] %--% mappedKmers[ind + 1]])
		  }
		}
		output = rep(NA, numKmers)
		output[missing] = MISSING
		output[!missing] = c(START, checkEdges) 
		return(output)
	}
	else {
		return()
	}
}

### This function is the batch version of the previous function: multiple sequences (Table) get processed.
identifyAllPaths = function(Table, Graph, Map, maxMissingFrac, count) {
	R = nrow(Table)
	Result = vector("list", R)
	for (ind in 1:R) {
		curSeq = Table[ind,]
		curName = rownames(Table)[ind]
		Result[[curName]] = identifyPath(curSeq, Graph, Map, maxMissingFrac = maxMissingFrac, count = count)
	}
	Result
}

### This function maps an individual pair of reads based on a list of unique K-mers
### If split = TRUE, then each read is partitioned into *non-overlapping* K-mers
mapReadPair = function(read1, read2, Kmers, unique, split) {
  uniqueKmers = Kmers[[1 + unique]]
  M = length(uniqueKmers)
  N = nchar(uniqueKmers[[1]][1])
  cur1 = getKmers(read1, N, FALSE, split = split)
  cur2 = getKmers(read2, N, FALSE, split = split)
  cur1R = getKmers(reverseComplement(read1), N, FALSE, split = split)
  cur2R = getKmers(reverseComplement(read2), N, FALSE, split = split)
  mult1 = lapply(1:M, function(x) {which(uniqueKmers[[x]] %in% cur1)}) ### later on: do this with hash tables?!
  mult2 = lapply(1:M, function(x) {which(uniqueKmers[[x]] %in% cur2)})
  mult1R = lapply(1:M, function(x) {which(uniqueKmers[[x]] %in% cur1R)})
  mult2R = lapply(1:M, function(x) {which(uniqueKmers[[x]] %in% cur2R)})
  multFR = sapply(mult1, length) + sapply(mult2R, length)
  multRF = sapply(mult1R, length) + sapply(mult2, length) 
  votes = which(multFR + multRF > 0)
  if (length(votes) > 0) {
    result = mapAmbiguousRead(read1, read2, Kmers, mult1, mult2, mult1R, mult2R, multFR, multRF, votes, unique)
    return(ifelse(is.na(result), NA, names(uniqueKmers)[result]))
  }
  else {
    return(NA)
  }
}

### This function maps a read for which the mapping is ambiguous
mapAmbiguousRead = function(read1, read2, Kmers, mult1, mult2, mult1R, mult2R, multFR, multRF, votes, unique) {
  if (length(votes) == 1) {
    return(votes)
  }
  maxError = computeNumErrors(length(read1)) + computeNumErrors(length(read2))
  bestFR = max(multFR)
  bestRF = max(multRF)
  scores = rep(Inf, length(Kmers[[1]]))
  if (bestFR <= bestRF) {
    Opts = which(multFR == bestFR)
    scores[Opts] = getScores(read1, reverseComplement(read2), Kmers, Opts, mult1[Opts], mult2R[Opts], unique)
  }
  if (bestRF <= bestFR) { # note: in case of equality, get both FR an RF scores!
    Opts = which(multRF == bestRF)
    scores[Opts] = getScores(reverseComplement(read1), read2, Kmers, Opts, mult1R[Opts], mult2[Opts], unique)
  }
  scores[scores > maxError] = NA
  bestScore = min(scores, na.rm = TRUE)
  bestPos = which(scores == bestScore)
  if (length(bestPos) != 1) {
    return(NA)
  }
  return(bestPos)
}

### This function computes the number of errors to be tolerated in a read 
computeNumErrors = function(readSize, errorRate = 0.01, quantile = 0.99) {
  stopifnot(errorRate > 0 && errorRate < 1)
  stopifnot(quantile  > 0 && quantile  < 1)
  maxErrors = qbinom(quantile, readSize, errorRate) 
  maxErrors
}

### This function computes the mismatch score of an alignment of two input reads to each element of the Options
### The corresponding lists of anchors provide positions where each of the reads has an exactly matching K-mer 
getScores = function(read1, read2, Kmers, Options, anchors1, anchors2, unique) {
  L = length(Options)
  scores1 = rep(NA, L)
  scores2 = rep(NA, L)
  for (ind in 1:L) {
    curOption = Options[ind]
    curAnchors1 = anchors1[[ind]]
    curAnchors2 = anchors2[[ind]]
    curKmers = lapply(Kmers, function(x) {x[[curOption]]})
    curScores1 = getAnchoredScores(read1, curAnchors1, curKmers, unique)
    scores1[ind] = min(curScores1)
    curScores2 = getAnchoredScores(read2, curAnchors2, curKmers, unique)
    scores2[ind] = min(curScores2)
  }
  totalScores = scores1 + scores2
  totalScores
}

getAnchoredScores = function(read, anchors, Kmers, unique) {
  N = nchar(Kmers[[1]][1])
  splitRead = getKmers(read, N, FALSE, split = TRUE)
  L = length(splitRead)
  if (length(anchors) == 0) {
    fullRef = reconstructReference(Kmers[[1]])
    scores = localAlignment(fullRef, read, match = 0, mismatch = 1, insertion = Inf, deletion = Inf)[[2]]
  }
  else {
    M = length(anchors)
    scores = rep(0, M)
    for (ind in 1:M) {
      anchor = anchors[ind]
      anchorPos = which(splitRead == Kmers[[1 + unique]][anchor])
      initPos = ifelse(unique, Kmers[[3]][anchor], anchor)
      for (index in 1:L) {
        if (index != anchorPos) {
          curPos = getTrim(N * (index - anchorPos) + initPos, length(Kmers[[1]]))
          scores[ind] = scores[ind] + HammingDistance(splitRead[index], Kmers[[1]][curPos[[1]]], curPos[[2]])
        }
      }
      tailEnd = read[-(1:(L * N))]
      if (length(tailEnd) > 0) {
        curPos = getTrim(N * (L - anchorPos + 1) + initPos, length(Kmers[[1]]))
        scores[ind] = scores[ind] + HammingDistance(tailEnd, Kmers[[1]][curPos[[1]]], RIGHT)
      }
    }
  }
  scores
}

localAlignment = function(sequence1, sequence2, match, mismatch, insertion, deletion) {
  LOCAL <<- LOCAL + 1
  L1 = length(sequence1)
  L2 = length(sequence2)
  stopifnot(L1 >= L2)
  if (is.infinite(insertion) && is.infinite(deletion) && match == 0) {
    numShifts = L1 - L2
    allScores = rep(Inf, numShifts)
    for (ind in 1 + (0:numShifts)) {
      curShift = sequence1[ind + (0:(L2 - 1))]
      allScores[ind] = HammingDistance(curShift, sequence2, NONE) * mismatch
    }  
    bestPos = which.min(allScores)
    bestCost = min(allScores)
    return(list(bestPos, bestCost))
  }
  else {
    print("Error: this configurations has not yet been implemented!")
    return(NULL)
  }
}

reconstructReference = function(Kmers) {
  N = nchar(Kmers[1])
  L = length(Kmers)
  Q = floor((L - 1)/N)
  R = L - 1 - Q * N
  allKmers = Kmers[1 + N * (0:Q)]
  reference = paste0(allKmers, collapse = "")
  reference = unlist(strsplit(reference, ""))
  if (R > 0) {
    tailEnd = tail(unlist(strsplit(Kmers[L], "")), R)
    reference = c(reference, tailEnd)
  }
  reference
}

getTrim = function(position, length) {
  trim = NONE
  if (position < 1) {
    position = 1
    trim = LEFT
  }
  if (position > length) {
    position = length
    trim = RIGHT
  }
  output = list(position, trim)
  output
}

HammingDistance = function(string, charVector, trim) {
  splitString = unlist(strsplit(string, ""))
  if (trim == LEFT) {
    charVector = tail(charVector, length(splitString))
  }
  if (trim == RIGHT) {
    charVector = head(charVector, length(splitString))  
  }
  distance = sum(splitString != charVector, na.rm = TRUE)
  distance
}

### This function computes the performance of a mapper based on read names and mapping results
### If hidden is not NULL, the precision and recall compare hidden to the unmapped reads only!
computePerformance = function(readNames, Mapped, hidden = NULL, cMap = NULL) {
  Mapped = mapEntities(Mapped, cMap)
  trueStrains = extractStrains(readNames)
  trueStrains = mapEntities(trueStrains, cMap)
  Unmapped = which(is.na(Mapped))
  if (!is.null(hidden)) {
    hidden = mapEntities(hidden, cMap)
    goodFinds = sum(trueStrains[Unmapped] == hidden)
    Precision = goodFinds/length(Unmapped)
    Recall = goodFinds/sum(trueStrains == hidden)
  }
  else {
    if (length(Unmapped) > 0) {
      trueStrains = trueStrains[-Unmapped]
      Mapped = Mapped[-Unmapped]
    }
    Table = table(trueStrains, Mapped)
    goodNames = intersect(rownames(Table), colnames(Table))
    Precision = diag(Table)/(rowSums(Table)[goodNames])
    Recall = diag(Table)/(colSums(Table)[goodNames])
  }
  Performance = list(Precision, Recall)
  Performance
}

### This function maps files of paired-end reads, taking reverse complementation of read pairs into account
### If hidden is not NULL, this corresponds to the name of the hidden strain (and we use read names to compare)
### If verbose is TRUE, then the precision and recall of the mapping are also printed out.
mapPairedReadsNew = function(F1, F2, Kmers, unique, split, hidden = NULL, verbose = FALSE, cMap = NULL) {
  L = nrow(F1)
  R = ncol(F1)
  stopifnot(L == nrow(F2))
  Mapped = rep(NA, L)
  for (ind in 1:L) {
    Mapped[ind] = mapReadPair(F1[ind,], F2[ind,], Kmers, unique, split)
  }
  fracMapped = length(na.omit(Mapped))/L
  print(paste0("A total of ", round(fracMapped * 100, NDEC), "% of the reads got mapped"))
  Fracs = table(na.omit(Mapped))/length(na.omit(Mapped))
  output = list(fracMapped, Fracs)
  Performance = NULL
  if (!is.null(hidden) || verbose) {
    Performance = computePerformance(rownames(F1), Mapped, hidden = hidden, cMap = cMap)
    if (verbose) {
      print(Performance) 
    }
  }
  output = c(output, list(Performance, which(is.na(Mapped))))
  output
}

### This function determines the K-mers in a string for a given K 
### If requireUnique is TRUE it also checks that they are unique.
### If split = TRUE, only returns maximal non-overlapping K-mers.
getKmers = function(string, K, requireUnique, split = FALSE) {
	string = string[string %in% ALPHABET]
	L = length(string)
	if (split) {
	  numKmers = floor(L/K)
	  startInds = 1 + K * (0:(numKmers - 1))
	}
	else {
	  numKmers = L - K + 1
	  startInds = 1:numKmers
	}
	allKmers = rep(NA, numKmers)
	for (ind in 1:numKmers) {
	  allKmers[ind] = paste(string[startInds[ind] + (0:(K-1))], collapse = "")
	}
	if (requireUnique && length(unique(allKmers)) < numKmers) {
		print("Not all K-mers are unique")
		return(NULL)
	}
	allKmers
}

### This function reads a nex file and outputs a matrix with the alignment
parseNex = function(nexFile) {
	Lines = read.nexus.data(nexFile)
	Result = matrix(unlist(Lines), nrow = length(Lines), byrow = TRUE, dimnames = list(names(Lines), NULL))
	Result
} 

### This function takes in an alignment and leaves one representative of any group of lines within a threshold
### Returns the reduced alignment as well as a map to convert from the eliminated rows to their representative
### If tree = TRUE, it also returns the neighbor-joining tree based on the distance matrix of the pruned types
pruneTypes = function(alignment, threshold) {
	Distances = computeAllPDistances(alignment)
	G = graph_from_adjacency_matrix(Distances <= threshold)
	G = simplify(G)
	C = clusters(G)
	classes = which(C$csize > 1)
	L = length(classes)
	clustersFound = vector("list", L)
	names(clustersFound) = classes
	for (class in classes) {
		curClass = which(C$membership == class)
		curSubgraph = induced_subgraph(G, curClass)
		curNames = V(G)[curClass]$name
		maxSize = 2 * choose(length(curClass), 2)
		stopifnot(ecount(curSubgraph) == maxSize)
		clustersFound[[as.character(class)]] = curNames
	}
	clusterMap = convertToMap(clustersFound)
	uniqueReps = !duplicated(C$membership)
	reducedAlignment = alignment[uniqueReps,]
	output = list(reducedAlignment, clusterMap)
	output
}

### This function converts a list of vectors into a map taking each non-head vector element into its head element
convertToMap = function(ListOfVectors) {
	Lens = sapply(ListOfVectors, length)
	Sums = cumsum(c(1, Lens[-length(Lens)]))
	Map = rep(NA, sum(Lens))
	firsts = sapply(ListOfVectors, function(x) {x[1]})
	Names = unlist(ListOfVectors)
	Map = rep(firsts, Lens)
	names(Map) = Names
	Map = Map[-Sums]
	Map
}

### This function takes in an alignment and a molecular evolutionary model and returns a matrix of P-distances
computeAllPDistances = function(alignment, model = "GTR", ignoreGaps = TRUE) {
	L = nrow(alignment)
	Distances = matrix(0, L, L)
	Names = rownames(alignment)
	dimnames(Distances) = list(Names, Names)
	for (ind1 in 1:(L-1)) {
		seq1 = alignment[ind1,]
		for (ind2 in (ind1 + 1):L) {
			seq2 = alignment[ind2,]
			dist = computePDistance(seq1, seq2, model = model, ignoreGaps = ignoreGaps)
			Distances[ind1, ind2] = dist
			Distances[ind2, ind1] = dist
		}
	}
	Distances
}

### This function takes in two sequences and a molecular evolutionary model and returns a P-distance
computePDistance = function(sequence1, sequence2, model, ignoreGaps) {
	if (model == "GTR") {
		Table = table(sequence1, sequence2)
		freqs = table(c(sequence1, sequence2))
		if (ignoreGaps) {
			Table = Table[ALPHABET, ALPHABET]
			freqs = freqs[ALPHABET]
		}
		dist = distanceFormula(Table, freqs)
	}
	dist
}

### This function estimates the distance according to the general time-revesible model
### The equation is taken from Ziheng Yang, "Computational Molecular Evolution" (1.68)
distanceFormula = function(Table, freqs) {
		Table = (Table + t(Table))/2
		Table = Table/sum(Table)
		freqs = freqs/sum(freqs)
		Freqs = diag(freqs)
		invFreqs = diag(1/freqs)
		dist = -sum(diag(Freqs %*% logm(invFreqs %*% Table)))
		dist
}

### This function parses a fastq file, returning a row-named matrix of characters
### If keepNames is TRUE, then it also keeps a short version of each read's name!
### If keepScores is TRUE, then it also keeps quality scores as a separate output.
### It fails if the reads are not all of the same length (TO BE CHANGED LATER ON)
parseF = function(fastqFile, keepNames = FALSE, keepScores = FALSE) {
	Lines = readLines(fastqFile)
	N = length(Lines)
	stopifnot(N %% 4 == 0)
	Reads = Lines[4 * (1:(N/4)) - 2]
	stopifnot(all(nchar(Reads) == nchar(Reads[1])))
	if (keepNames) {
		ReadIDs = Lines[4 * (1:(N/4)) - 3]
		shortIDs = unlist(lapply(ReadIDs, function(x) {tail(unlist(strsplit(x,":")),1)}))
	}
	Reads = lapply(Reads, function(x) {unlist(strsplit(tolower(x),""))})
	Reads = matrix(unlist(Reads), N/4, byrow = TRUE)
	if (keepNames) {
		rownames(Reads) = shortIDs
	}
	output = list(Reads)
	if (keepScores) {
	  Scores = Lines[4 * (1:(N/4))]
	  if (keepNames) {
	    names(Scores) = shortIDs
	  }
	  output = c(output, list(Scores))
	}
	output
}

### This function writes a fastq file, based on a given matrix of reads and scores, into a specified filename
writeF = function(Reads, Scores, Filename) {
  Names = rownames(Reads)
  L = nrow(Reads)
  stopifnot(L == length(Scores))
  names(Reads) = NULL
  names(Scores) = NULL
  Lines = rep("+", 4 * L)
  Lines[4 * (1:L) - 3] = Names
  Lines[4 * (1:L) - 2] = apply(Reads, 1, function(x) {paste(x, collapse = "")})
  Lines[4 * (1:L)] = Scores
  fw = file(Filename, "w")
  writeLines(Lines, fw)
  close(fw)
}

### This function identifies the strain name for each read based on its name
extractStrains = function(readNames) {
  shortNames = unlist(strsplit(readNames,"-"))[2 * (1:length(readNames)) - 1]
  strains = gsub("@", "", shortNames)
  strains
}

### This function gets the actual fractions of reads from each strain in a file
getTrueFractions = function(filename) {
  parsedFile  = parseF(filename, keepNames = TRUE, keepScores = FALSE)[[1]]
  Table = table(extractStrains(rownames(parsedFile)))
  actualFracs = Table/sum(Table)
  actualFracs
}

### This function computes the total variation distance between 2 distributions
TVDistance = function(distribution1, distribution2) {
  Names1 = names(distribution1)
  Names2 = names(distribution2)
  Names12 = intersect(Names1, Names2)
  Names11 = setdiff(Names1, Names12)
  Names22 = setdiff(Names2, Names12)
  Sum0 = sum(abs(distribution1[Names12] - distribution2[Names12]))
  Sum1 = sum(distribution1[Names11])
  Sum2 = sum(distribution2[Names22])
  Dist = (Sum0 + Sum1 + Sum2)/2
  Dist
}

### This function computes the TV distance between stated and actual distribution of reference strains in a file
testFractions = function(filename, returnFracs = FALSE) {
  statedFracs = parseFilenameFracs(filename)
  actualFracs = getTrueFractions(filename)
  if (returnFracs) {
    output = list(statedFracs, actualFracs)
  }
  else {
    output = TVDistance(actualFracs, statedFracs)
  }
  output
}

### This function checks if the reads are properly paired; if not, overwrites File2 by a corrected version
correctReadPairing = function(File1, File2) {
  PF1 = parseF(File1, keepNames = TRUE, keepScores = TRUE)
  PF2 = parseF(File2, keepNames = TRUE, keepScores = TRUE)
  Names1 = rownames(PF1[[1]])
  Names2 = rownames(PF2[[1]])
  stopifnot(length(Names1) == length(Names2))
  shortNames1 = gsub("/1", "", Names1)
  shortNames2 = gsub("/2", "", Names2)
  stopifnot(all(sort(shortNames1) == sort(shortNames2)))
  if (any(shortNames1 != shortNames2)) {
    print(paste("Correcting the order in", File2, "relative to", File1))
    order = match(shortNames1, shortNames2)
    writeUnmappedReads(PF2, order, File2)
  }
  else {
    print(paste("The files", File1, "and", File2, "are consistently paired; no changes made!"))
  }
}

### This function preprocesses a directory to consistently order all the paired-end read files
preprocessDirectory = function(Directory) {
  initDir = getwd()
  setwd(Directory)
  allFiles = list.files(pattern = ".fq")
  allFiles = gsub("R1.fq", "", allFiles)
  allFiles = gsub("R2.fq", "", allFiles)
  allFiles = unique(allFiles)
  allFiles = allFiles[grep("Hidden", allFiles, invert = TRUE)]
  for (File in allFiles) {
    File1 = paste0(File, "R1.fq")
    File2 = paste0(File, "R2.fq")
    correctReadPairing(File1, File2)
  }
  setwd(initDir)
}

### This function computes the difference between stated and actual strain fractions for all files in a directory
testAllFractions = function(directory) {
  initDir = getwd()
  setwd(directory)
  allFiles = list.files(pattern = ".fq")
  L = length(allFiles)
  allDiffs = rep(NA, L)
  names(allDiffs) = allFiles
  for (file in allFiles) {
    allDiffs[file] = testFractions(file, returnFracs = FALSE)
  }
  setwd(initDir)
  allDiffs
}

### This function describes what we are actually executing
main = function() {
  allDirs = as.vector(outer(c("AB mixtures", "Complex mixtures"), c("", " error-free"), paste0))
  if (CORRECT) {
    for (directory in allDirs[c(2,4)]) {
      preprocessDirectory(directory)
    }
  }
  if (NEW) {
    setwd("Accurate mixtures/")
  }
  for (HIDDEN in c(FALSE, TRUE)) {
    for (K in c(15, 18, 25)) {
      print(c(HIDDEN, K))
      if (!(HIDDEN || NEW)) {
        allDirs = c("Individual Strains", allDirs)
      }
      for (directory in allDirs) {
        setwd(directory)
        print(directory)
        Res = fullWrapper(maxFrac = MMF, threshold = THRESH, K = K, split = SPLIT, graph = GRAPH, unique = UNIQ, hidden = HIDDEN)
        setwd("..")
      }
    }
  }
}