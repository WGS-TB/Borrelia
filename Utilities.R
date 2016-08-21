library(ape)
library(expm)
library(igraph)
library(phylobase)

ALPHABET = c("a","c","g","t")
MISSING = -2  # code for missing K-mers
START = -1    # code for starting K-mers

MMF = 0.5     # maximum missing fraction for a path
THRESHOLD = 0.05 # minimum pairwise distance
MINK = 15     # shortest K-mer to try for the library
MAXK = 30     # longest K-mer to try for the library
NDEC = 1      # number of decimals; used for rounding
SAVE = TRUE   # TRUE if coverage vectors are saved
BYREAD = TRUE # TRUE if the reads get assigned to types
TREE = FALSE  # TRUE if types are processed in "tree order"
SPLIT = TRUE  # TRUE if reads become non-overlapping K-mers
UNIQ = TRUE   # TRUE if only unique K-mers from library used
GRAPH = TRUE  # TRUE if using de Bruijn graphs to preprocess
HIDDEN = TRUE # TRUE if reference strains are hidden, one at a time
WRITE = TRUE  # TRUE if the unmapped reads are written to a file
SEP = ","     # separator used internally for type names

### This function is a wrapper function for finding the copy numbers of all repeats
wrapperRepeats = function(baseFilename, templateFile = "template.txt") {
	File1 = paste0(baseFilename, ".R1.fq.txt")
	File2 = paste0(baseFilename, ".R2.fq.txt")
	F1 = parseF(File1, keepNames = FALSE, keepScores = FALSE)
	F2 = parseF(File2, keepNames = FALSE, keepScores = FALSE)
	
}

### UNFINISHED FUNCTION!
myAlign = function(sequence1, sequence2, costMatrix = NULL, local = FALSE) {
	L1 = length(sequence1)
	L2 = length(sequence2)
	M = matrix(NA, L1, L2)
	for (ind1 in 1:L1) {
		for (ind2 in 1:L2) {
			### TO BE COMPLETED LATER!
		}
	}
}

### This function is a wrapper function for all the relevant files found in the directory.
fullWrapperKnown = function(databaseFile="../DatasetW1.nex.txt", threshold, byRead, tree, split, graph, unique, hidden, write) {
	myAlign = parseNex(databaseFile)
	prunedTypes = pruneTypes(myAlign, threshold = threshold, tree = tree)
	redAlign = prunedTypes[[1]]
	clusterMap = prunedTypes[[2]]
	Tree = NA
	if (tree) {
		Tree = prunedTypes[[3]]
	}
	allKmers = getUniqueKmers(redAlign)
	N = allKmers[[3]]
	allFiles = list.files(pattern = ".fq")
	allFiles = gsub(".txt", "", allFiles)
	allFiles = gsub("R1.fq", "", allFiles)
	allFiles = gsub("R2.fq", "", allFiles)
	allFiles = unique(allFiles)
	Output = vector("list", length(allFiles))
	names(Output) = allFiles
	for (curFile in allFiles) {
		print(curFile)
		# Output[[curFile]] = wrapperKnown(curFile, redAlign, N, byRead = byRead, tree = Tree)
	  Output[[curFile]] = wrapperKnownNew(curFile, redAlign, N, graph = graph, unique = unique, split = split, hidden = hidden, write = write, clusterMap = clusterMap)
	}
	save(Output, file = "FullOutput.RData")
	startFile = unlist(strsplit(databaseFile,"\\."))[1]
	opts = paste0(rep("ByRead",byRead),rep("Tree",tree),rep("Split",split),rep("Graph",graph),rep("Unique",unique),rep("Hidden",hidden)) 
	outputFile = paste0(startFile, "ResultsFrac", MMF, "Threshold", threshold, "K", N, opts, ".csv")
	Table = formatOutput(Output, outputFile, clusterMap, hidden = hidden)
	Table
}

### This function takes the output of a type frequency analysis, converts it into a table and writes it into a file
formatOutput = function(output, outputFile, clusterMap = NULL, hidden = HIDDEN, nDec = NDEC) {
	Lens = sapply(output, function(x) {length(x[[1]])})
	N = length(output)
	if (hidden) {
	  DimNames = list(NULL, c("Type", "True (%)", "Precision (%)", "Recall (%)"))
	}
	else {
	  DimNames = list(NULL, c("Type", "True (%)", "Predicted (%)", "Error (%)"))
	}
	Table = matrix(NA, nrow = sum(Lens) + N, ncol = 4, dimnames = DimNames)
	pos = 0
	for (ind in 1:N) {
		curPiece = output[[ind]]
		curTruth = curPiece[[1]]
		curPred = curPiece[[2]]
		curLen = Lens[ind]
		curMat = matrix(0, curLen, 3)
		curNames = names(curTruth)
		if (!is.null(clusterMap)) {
			curNames = mapEntities(curNames, clusterMap)
		}
		rownames(curMat) = curNames
		curMat[, 1] = round(curTruth * 100, nDec)
		if (hidden) {
		  curMat[colnames(curPred), 2] = round(curPred[1, colnames(curPred)] * 100, nDec)
		  curMat[colnames(curPred), 3] = round(curPred[2, colnames(curPred)] * 100, nDec)
		}
		else {
		  curMat[names(curPred), 2] = round(curPred * 100, nDec)
		  curMat[, 3] = round(abs(curMat[ ,1] - curMat[ ,2]), nDec)
		  curError = sum(curMat[,3])
		}
		curMat = cbind(rownames(curMat), curMat)
		curMat = curMat[order(curMat[,1]), ]
		curMat = rbind(curMat, c(rep("", 3), ifelse(hidden, "", curError)))
		Table[pos + (1:(curLen + 1)), ] = curMat
		pos = pos + curLen + 1
	}
	write.table(Table, file = outputFile, row.names = FALSE, col.names = TRUE, quote = FALSE, sep = ",")
	Table
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
	stopifnot(!is.null(Map))
	mappedEntities = Map[Entities]
	bad = is.na(mappedEntities)
	mappedEntities[bad] = sapply(Entities[bad], default)
	mappedEntities
}

### This function is a new wrapper function that takes reverse complementarity of paired-end reads into account.
wrapperKnownNew = function(baseFilename, redAlign, N, graph, unique, split, hidden, write, clusterMap = NULL) {
  correctFractions = parseFilenameFracs(baseFilename)
  File1 = paste0(baseFilename, "R1.fq")
  if (!File1 %in% list.files()) {
    File1 = paste0(File1, ".txt")
  }
  File2 = paste0(baseFilename, "R2.fq")
  if (!File2 %in% list.files()) {
    File2 = paste0(File2, ".txt")
  }
  PF1 = parseF(File1, keepNames = TRUE, keepScores = hidden)
  PF2 = parseF(File2, keepNames = TRUE, keepScores = hidden)
  if (hidden) {
    F1 = PF1[[1]]
    F2 = PF2[[1]]
  }
  else {
    F1 = PF1
    F2 = PF2
  }
  if (graph) {
    print("Constructing the de Bruijn graph")
    RCF1 = t(apply(F1, 1, reverseComplement))
    RCF2 = t(apply(F2, 1, reverseComplement))
    G0 = constructDBGraph(rbind(F1, F2, RCF1, RCF2), N)
    Res = identifyPaths(redAlign, G0[[1]], N, G0[[3]], remove = FALSE, maxMissingFrac = MMF, count = FALSE)
    goodStrains = names(Res)[!sapply(Res, is.null)]
  }
  else {
    goodStrains = rownames(redAlign)
  }
  if (hidden) {
    trueNames = names(correctFractions)
    if (!is.null(clusterMap)) {
      trueNames = mapEntities(trueNames, clusterMap) 
    }
    C = matrix(NA, 2, length(trueNames), dimnames = list(c("Precision", "Recall"), trueNames))
    for (Name in trueNames) {
      curStrains = setdiff(goodStrains, Name)
      curQ = getUniqueKmers(redAlign[curStrains, , drop = FALSE], N)
      curResult = mapPairedReadsNew(F1, F2, curQ[[1 + unique]], split = split, hidden = Name)
      C[,Name] = curResult[[3]]
      unmappedReads = curResult[[4]]
      if (write) {
        Reads1  = F1[unmappedReads, , drop = FALSE]
        Reads2  = F2[unmappedReads, , drop = FALSE]
        Scores1 = PF1[[2]][unmappedReads]
        Scores2 = PF2[[2]][unmappedReads]
        Filename1 = gsub(".R1.fq", paste0("Hidden", Name, ".R1.fq"), File1)
        Filename2 = gsub(".R2.fq", paste0("Hidden", Name, ".R2.fq"), File2)
        writeF(Reads1, Scores1, Filename1)
        writeF(Reads2, Scores2, Filename2)
      }
    }
  }
  else {
    if (length(goodStrains) > 1) {
      Qred = getUniqueKmers(redAlign[goodStrains, , drop = FALSE], N)
      # C = mapPairedReads(F1, F2, Qred[[2]], tree = redTree)[[2]]
      C = mapPairedReadsNew(F1, F2, Qred[[1 + unique]], split = split)[[2]]
    }
    else {
      C = 1
      names(C) = goodStrains[1]
    }
  }
  output = list(correctFractions, C)
  output
}

readRefFile = function(filename) {
	Lines = readLines(filename)
	L = length(Lines)
	M = max(nchar(Lines))
	Mat = matrix('-', L, M)
	for (ind in 1:L) {
		curLine = Lines[ind]
		Mat[ind, 1:nchar(curLine)] = unlist(strsplit(tolower(curLine), ""))
	}
	Mat
}

hashReference = function(refTable, K, allowRepeats = FALSE, duplicate = FALSE) {
	M = nrow(refTable)
	N = ncol(refTable) + duplicate * K
	maxSize = round(M * N / K)
	hash = new.env(hash = TRUE, size = maxSize)
	allKmers = vector("list", M)
	for (ind in 1:M) {
		curReference = refTable[ind,]
		if (duplicate) {
			curReference = c(curReference, curReference[1:(K - 1)])
		}
		allKmers[[ind]] = getKmers(curReference, K, requireUnique = (!allowRepeats))
	}
	uniqueKmers = unique(unlist(allKmers))
	numKmers = length(uniqueKmers)
	hashedList = vector("list", numKmers)
	for (ind in 1:numKmers) {
		hashedList[[ind]] = matrix(NA,2,0)
	}
	names(hashedList) = uniqueKmers
	for (ind in 1:M) {
		curKmers = allKmers[[ind]]
		if (allowRepeats) {
			for (index in 1:length(curKmers)) {
				curKmer = curKmers[index]
				hashedList[[curKmer]] = cbind(hashedList[[curKmer]], c(ind, index))
			}
		}
		else {
			allInds = match(curKmers, uniqueKmers)
			curLen = length(allInds)
			hashedList[allInds] = lapply(1:curLen, function(x) {cbind(hashedList[[allInds[x]]], c(ind, x))})
		}
	}
	hashedList
}

### This function computes the number of errors to be tolerated in a read 
computeNumErrors = function(readSize, errorRate = 0.01, quantile = 0.99) {
	stopifnot(errorRate > 0 && errorRate < 1)
	stopifnot(quantile  > 0 && quantile  < 1)
	maxErrors = qbinom(quantile, readSize, errorRate) 
	maxErrors
}

### This function returns the reference to which each read in a read file can map
### If allowHalf = TRUE then more errors are allowed (only half a read must map);
### if revComp = TRUE then every read is reverse complemented before being mapped
mapReads = function(readFile, refTable, allowHalf = FALSE, revComp = FALSE) {
	M = ncol(readFile)
	numErrors = computeNumErrors(M)
	K = floor(M / (numErrors + 1))
	print(paste("Using k-mers of length", K))
	refHash = hashReference(refTable, K, duplicate = allowHalf)
	splitIndices = rep(1:(numErrors + 1), each = K)
	effectiveLength = length(splitIndices)
	L = nrow(readFile)
	finalMaps = vector("list", L)
	print(paste("There are", L, "reads to process"))
	for (ind in 1:L) {
		if (ind %% 100000 == 0) {
			print(ind)
		}
		curRead = readFile[ind, ]
		if (revComp) {
			curRead = reverseComplement(curRead)
		}
		curRead = curRead[1:effectiveLength]
		splitRead = split(curRead, splitIndices)
		splitReads = sapply(splitRead, function(x) {paste(x, collapse = "")})
		allMaps = refHash[splitReads]
		candidates = unlist(sapply(allMaps, function(x) {x[1,]}))
		if (length(candidates) > 0) {
			best = findModes(candidates)
			Len = length(best)
			if (Len > 1) {
				print('Using positions to decide')
				curBests = rep(NA, Len)
				for (index in 1:Len) {
					curBests[index] = findMaxCompatible(allMaps, best[index], K)
				}
				best = best[curBests == max(curBests)]
			}
			finalMaps[[ind]] = lapply(allMaps, function(x) {x[, x[1,] %in% best]})
		}
	}
	finalMaps
}

### This function finds the most frequent elements in a given vector
findModes = function(Vector) {
	votes = table(Vector)
	best = names(votes)[votes == max(votes)]
	best = as(best, class(Vector))
	best
}

### This function computes the set of leaf partitions induced by each internal branch in a tree
### For instance, for the tree ((A,B),(C,D)) the partitions are AB|CD, A|B, B|A, CD|AB, C|D, D|C.
computeContrasts = function(tree, node = NULL) {
	if (is.null(node)) {
		node = rootNode(tree)
	}
	kids = children(tree, node)
	if (length(kids) == 0) {
		return(list())
	}
	stopifnot(length(kids) == 2)
	leftKid = kids[1]
	rightKid = kids[2]
	leftTips = names(descendants(tree, leftKid, "tips"))
	rightTips = names(descendants(tree, rightKid, "tips"))
	rootContrasts = list(list(leftTips, rightTips), list(rightTips, leftTips))
	contrasts = c(rootContrasts, computeContrasts(tree, leftKid), computeContrasts(tree, rightKid))
	contrasts
}

### This function parses the name of the file to determine the correct fractions
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

### This function finds the minimum K such that each row in a table contains a unique K-mer
getUniqueKmers = function(Table, minK = MINK, maxK = MAXK, contrastList = NULL, intersect = FALSE) {
	N = nrow(Table)
	output = NULL
	KmerList = vector("list", N)
	Names = rownames(Table)
	names(KmerList) = Names
	if (is.null(contrastList)) {
		contrastList = vector("list", N)
		names(contrastList) = Names
		for (ind in 1:N) {
			contrastList[[ind]] = list(ind, setdiff(1:N, ind))
		}
	}
	M = length(contrastList)
	uniqueKmers = vector("list", M)
	# names(uniqueKmers) = sapply(contrastList, function(x) {paste(x[[1]], collapse = SEP)})
	names(uniqueKmers) = names(contrastList)
	for (K in minK:maxK) {
		FOUND = TRUE
		for (ind in 1:N) {
			curKmers = getKmers(Table[ind,], K, TRUE)
			if (is.null(curKmers)) {
				break
			}
			else {
				KmerList[[ind]] = curKmers
			}
		}
		for (ind in 1:M) {
			curContrast = contrastList[[ind]]
			if (intersect) {
				curList = KmerList[[curContrast[[1]][1]]]
				if (length(curContrast[[1]]) > 1) {
					for (index in 2:length(curContrast[[1]])) {
						curList = intersect(curList, KmerList[[curContrast[[1]][index]]])
					}
				}
			}
			else {
				curList = unique(unlist(KmerList[curContrast[[1]]]))
			}
			otherLists = unique(unlist(KmerList[curContrast[[2]]]))
			curUnique = curList[!curList %in% otherLists]
			if (length(curUnique) == 0 && (!intersect || length(curContrast[[1]]) == 1)) {
				FOUND = FALSE
				break
			}
			else {
				uniqueKmers[[ind]] = curUnique
			}
		}
		if (FOUND) {
			output = list(KmerList, uniqueKmers, K)
			break
		}
	}
	output
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
### If remove = TRUE, the found path is removed. If more than maxMissingFrac of the K-mers
### are absent, the return value is NULL; otherwise, the K-mers are labelled and counted (unless count = FALSE).
identifyPath = function(sequence, Graph, K, Map, remove, maxMissingFrac, count) {
	curKmers = getKmers(sequence, K, requireUnique = FALSE)
	numKmers = length(curKmers)
	mappedKmers = Map[curKmers]
	missing = is.na(mappedKmers)
	mappedKmers = mappedKmers[!missing]
	if (sum(missing)/numKmers < maxMissingFrac) {
		L = length(mappedKmers)
		if (count) {
		  checkEdges = rep(0, L - 1)
		  for (ind in 1:(L-1)) {
		    checkEdges[ind] = length(E(Graph)[mappedKmers[ind] %--% mappedKmers[ind + 1]])
		  }
		  if (remove && all(checkEdges > 0)) {
		    Graph = delete.edges(Graph, cbind(mappedKmers[-L], mappedKmers[-1]))
		  }
		}
		else {
		  checkEdges = rep(1, L - 1)
		}
		output = rep(NA, numKmers)
		output[missing] = MISSING
		output[!missing] = c(START, checkEdges) 
		output
	}
	else {
		# print("Error: too few k-mers found in the graph")
		return()
	}
}

### The batch version of the previous function: multiple sequences (Table) get processed.
identifyPaths = function(Table, Graph, K, Map, remove, maxMissingFrac, count) {
	R = nrow(Table)
	Result = vector("list", R)
	for (ind in 1:R) {
		print(ind)
		curSeq = Table[ind,]
		curName = rownames(Table)[ind]
		Result[[curName]] = identifyPath(curSeq, Graph, K, Map, remove = remove, maxMissingFrac = maxMissingFrac, count = count)
	}
	Result
}

### This function works like the one below, but takes reverse complementation of paired-end reads into account
### If hidden is not NULL, this corresponds to the name of the hidden strain (and we use read names to compare)
mapPairedReadsNew = function(F1, F2, uniqueKmers, split, hidden = NULL) {
  L = nrow(F1)
  R = ncol(F1)
  stopifnot(L == nrow(F2))
  Maps = vector("list", L)
  M = length(uniqueKmers)
  N = nchar(uniqueKmers[[1]][1])
  Mapped = rep(FALSE, L)
  FLAG = FALSE
  for (ind in 1:L) {
    cur1 = getKmers(F1[ind,], N, FALSE, split = SPLIT)
    cur2 = getKmers(F2[ind,], N, FALSE, split = SPLIT)
    mult1 = sapply(1:M, function(x) {sum(uniqueKmers[[x]] %in% cur1)})
    mult2 = sapply(1:M, function(x) {sum(uniqueKmers[[x]] %in% cur2)})
    set1 = which(mult1 > 0)
    set2 = which(mult2 > 0)
    L1 = length(set1)
    L2 = length(set2)
    if (L1 < L2) {
      if (L1 > 0) {
        print("Warning: paired-end read mapped in the same orientation both times!")
        print(ind)
      }
      revCompRead = reverseComplement(F1[ind,])
    }
    else {
      if (L1 == L2) {
        if (L2 > 0) {
          print("Warning: paired-end read produces an equal number of mappings; checking multiplicities")
          if (L1 == 1) {
            M1 = mult1[set1]
            M2 = mult2[set2]
            if (M1 == M2) {
              print("Unable to determine the best mapping; reverse-complementing by default!")
            }
            else {
              print("Correct read to be reverse-complemented determined!")
              if (M2 > M1) {
                FLAG = TRUE
              }
            }
          }
          print(ind)
        }
      }
      else if (L2 > 0) {
        print("Warning: paired-end read mapped in the same orientation both times!")
        print(ind)
      }
      if (FLAG) {
        revCompRead = reverseComplement(F1[ind,])
        FLAG = FALSE
      }
      else {
        revCompRead = reverseComplement(F2[ind,])
      }
    }
    cur3 = getKmers(revCompRead, N, FALSE, split = SPLIT)
    mult3 = sapply(1:M, function(x) {sum(uniqueKmers[[x]] %in% cur3)})
    set3 = which(mult3 > 0) 
    votes = intersect(union(set1, set2), set3)
    if (length(votes) > 0) {
      if (length(votes) == 1) {
        Mapped[ind] = names(uniqueKmers)[votes]
      }
      else {
        print(paste("There are", length(votes), "equally plausible options"))
        print(ind)
      }
    }
  }
  fracMapped = mean(Mapped != FALSE)
  Unmapped = which(Mapped == FALSE)
  Mapped = Mapped[Mapped != FALSE]
  Fracs = table(Mapped)/length(Mapped)
  print(paste0("A total of ", round(fracMapped * 100, NDEC), "% of the reads got mapped"))
  output = list(fracMapped, Fracs)
  if (!is.null(hidden)) {
    trueStrains = extractStrains(rownames(F1))
    goodFinds = sum(trueStrains[Unmapped] == hidden)
    Precision = goodFinds/length(Unmapped)
    Recall = goodFinds/sum(trueStrains == hidden)
    output = c(output, list(c(Precision, Recall), Unmapped))
  }
  output
}

### This function gets the mean coverage for the unique K-mers identified in the specified good strains.
### If frac = TRUE, the result is given as a fraction out of 100%; if not, it is given as mean coverage.
getMeanCoverage = function(Result, uniqueKmers, goodStrains, save, baseFilename, tree = NULL) {
	redResult = Result[goodStrains]
	L = length(uniqueKmers[[2]])
	if (!is.null(tree)) {
		fullNames = names(uniqueKmers[[2]])
		splitNames = sapply(fullNames,  function(x) {unlist(strsplit(x, SEP))[1]})
	}
	else {
		fullNames = goodStrains
		splitNames = goodStrains
	}
	Matches = lapply(1:L, function(x) {match(uniqueKmers[[2]][[x]], uniqueKmers[[1]][[splitNames[x]]])})
	names(Matches) = fullNames
	fullCoverages = lapply(1:L, function(x) {redResult[[splitNames[x]]][Matches[[fullNames[x]]]]})
	names(fullCoverages) = fullNames
	if (save) {
		save(fullCoverages, file = paste0(baseFilename, "Coverages.RData"))
	}
	Coverages = sapply(fullCoverages, function(Z) {mu = mean(Z[Z >= 0]); mu})
	if (is.null(tree)) {
		names(Coverages) = goodStrains
	}
	else {
		L = length(Coverages)
		for (ind in 1:(L/2)) {
			curPair = Coverages[2 * (ind - 1) + 1:2]
			if (any(is.na(curPair))) {
				print(paste("Uninformative split at index", ind))
				print(paste(names(Coverages)[2 * (ind - 1) + 1:2], collapse = " vs. "))
			}
			Coverages[2 * (ind - 1) + 1:2] = curPair/sum(curPair)
		}
		Fracs = rep(1, length(goodStrains))
		names(Fracs) = goodStrains
		Edges = edges(tree)
		for (ind in 2:nrow(Edges)) { # skip the edge leading to the root!
			node = Edges[ind, 2]
			Labels = labels(tree)[descendants(tree, node, "tips")]
			fullLabel = paste(Labels, collapse = SEP)
			if (!any(is.na(Coverages[fullLabel]))) {
				Fracs[Labels] = Fracs[Labels] * Coverages[fullLabel]
			}
		}
		Coverages = Fracs
	}
	Coverages = Coverages[!is.na(Coverages)]
	Coverages = Coverages/sum(Coverages)
	Coverages
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
pruneTypes = function(alignment, threshold, tree) {
	Distances = computeAllPDistances(alignment)
	G = graph_from_adjacency_matrix(Distances <= threshold)
	G = simplify(G)
	C = clusters(G)
	classes = which(C$csize > 1)
	L = length(classes)
	clustersFound = vector("list", L)
	names(clustersFound) = classes
	print(paste("There are", L, "clusters in the data"))
	for (class in classes) {
		curClass = which(C$membership == class)
		curSubgraph = induced_subgraph(G, curClass)
		curNames = V(G)[curClass]$name
		# print(curNames)
		maxSize = choose(length(curClass), 2)
		if (ecount(curSubgraph) < maxSize) {
			print("Warning: the current cluster is not a clique!")
		}
		clustersFound[[as.character(class)]] = curNames
	}
	clusterMap = convertToMap(clustersFound)
	uniqueReps = !duplicated(C$membership)
	reducedAlignment = alignment[uniqueReps,]
	output = list(reducedAlignment, clusterMap)
	if (tree) {
		Tree = nj(Distances[uniqueReps , uniqueReps])
		rTree = root(Tree, outgroup = tail(rownames(reducedAlignment), 1), resolve.root = TRUE)
		output = c(output, list(rTree))
	}
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
		shortIDs = unlist(lapply(ReadIDs, function(x) {
			tail(unlist(strsplit(x,":")),1)
		}))
	}
	Reads = lapply(Reads, function(x) {unlist(strsplit(tolower(x),""))})
	Reads = matrix(unlist(Reads), N/4, byrow = TRUE)
	if (keepNames) {
		rownames(Reads) = shortIDs
	}
	if (keepScores) {
	  Scores = Lines[4 * (1:(N/4))]
	  if (keepNames) {
	    names(Scores) = shortIDs
	  }
	  Reads = list(Reads, Scores)
	}
	Reads
}

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

extractStrains = function(readNames) {
  shortNames = unlist(strsplit(readNames,"-"))[2 * (1:length(readNames)) - 1]
  strains = gsub("@", "", shortNames)
  strains
}

# setwd("Individual Strains/")
# ResInd = fullWrapperKnown(threshold=THRESHOLD, byRead=BYREAD, tree=TREE, split=SPLIT, graph=GRAPH, unique=UNIQ, hidden=HIDDEN, write=WRITE)
# setwd("..")
# setwd("AB mixtures/")
# ResAB = fullWrapperKnown(threshold=THRESHOLD, byRead=BYREAD, tree=TREE, split=SPLIT, graph=GRAPH, unique=UNIQ, hidden=HIDDEN, write=WRITE)
# setwd("..")
setwd("Complex mixtures/")
ResComp = fullWrapperKnown(threshold=THRESHOLD, byRead=BYREAD, tree=TREE, split=SPLIT, graph=GRAPH, unique=UNIQ, hidden=HIDDEN, write=WRITE)
setwd("..")