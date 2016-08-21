### This function identifies the sequence to which each paired read should map
mapPairedReads = function(F1, F2, uniqueKmers, tree = NULL) {
  L = nrow(F1)
  stopifnot(L == nrow(F2))
  Maps = vector("list", L)
  M = length(uniqueKmers)
  N = nchar(uniqueKmers[[1]][1])
  Fracs = rep(0, M)
  allNames = sapply(names(uniqueKmers), function(x) {unlist(strsplit(x, SEP))})
  redNames = unlist(allNames[sapply(allNames, length) == 1])
  names(Fracs) = redNames
  for (ind in 1:L) {
    cur1 = getKmers(F1[ind,], N, FALSE)
    cur2 = getKmers(F2[ind,], N, FALSE)
    if (!is.null(tree)) {
      node = rootNode(tree)
      kids = children(tree, node)
      cur = c()
      while (length(kids) > 0) {
        stopifnot(length(kids) == 2)
        leftLabels = labels(tree)[descendants(tree, kids[1], "tips")]
        leftLabel = paste(leftLabels, collapse = SEP)
        rightLabels = labels(tree)[descendants(tree, kids[2], "tips")]
        rightLabel = paste(rightLabels, collapse = SEP)
        leftBranch1 = any(uniqueKmers[[leftLabel]] %in% cur1)
        rightBranch1 = any(uniqueKmers[[rightLabel]] %in% cur1)
        leftBranch2 = any(uniqueKmers[[leftLabel]] %in% cur2)
        rightBranch2 = any(uniqueKmers[[rightLabel]] %in% cur2)
        anyLeft = leftBranch1 || leftBranch2
        anyRight = rightBranch1 || rightBranch2
        stopifnot(!(anyLeft && anyRight))
        if (anyLeft) {
          node = kids[1]
        }
        else if (anyRight) {
          node = kids[2]
        }
        else {
          break
        }
        kids = children(tree, node)
      }
      if (length(kids) == 0) {
        cur = labels(tree)[node]
      }
      Maps[[ind]] = cur
      Fracs[cur] = Fracs[cur] + 1
    }
    else {
      set1 = which(sapply(1:M, function(x) {any(uniqueKmers[[x]] %in% cur1)}))
      set2 = which(sapply(1:M, function(x) {any(uniqueKmers[[x]] %in% cur2)}))
      if (length(setdiff(set1, set2)) == 0) {
        cur = set2
      }
      else if (length(setdiff(set2, set1)) == 0) {
        cur = set1
      }
      else {
        cur = c()
      }
      stopifnot(length(cur) <= 1)
      Maps[[ind]] = cur
      Fracs[cur] = Fracs[cur] + 1
    }
  }
  print(paste(round(sum(Fracs) / L * 100, 2), "% of the reads got assigned"))
  Fracs = Fracs[Fracs > 0]
  Fracs = Fracs/sum(Fracs)
  output = list(Maps, Fracs)
  output
}

### This function is a wrapper function for a case of only known strains present, error-free reads.
wrapperKnown = function(baseFilename, redAlign, N, byRead, tree) {
  correctFractions = parseFilenameFracs(baseFilename)
  File1 = paste0(baseFilename, "R1.fq")
  if (!File1 %in% list.files()) {
    File1 = paste0(File1, ".txt")
  }
  File2 = paste0(baseFilename, "R2.fq")
  if (!File2 %in% list.files()) {
    File2 = paste0(File2, ".txt")
  }
  F1 = parseF(File1, keepNames = TRUE)
  F2 = parseF(File2, keepNames = TRUE)
  print("Constructing the de Bruijn graph")
  G0 = constructDBGraph(rbind(F1, F2), N)
  Res = identifyPaths(redAlign, G0[[1]], N, G0[[3]], remove = FALSE, maxMissingFrac = MMF, count = TRUE)
  goodStrains = names(Res)[!sapply(Res, is.null)]
  redTree = NULL
  contrastList = NULL
  if (!is.na(tree) && length(goodStrains) > 2) {
    redTree = subset(phylo4(tree), tips.include = goodStrains)
    contrastList = computeContrasts(redTree)
  }
  Qred = getUniqueKmers(redAlign[goodStrains, , drop = FALSE], N, contrastList = contrastList, intersect = (!is.na(tree) && !byRead))
  if (byRead) {
    C = mapPairedReads(F1, F2, Qred[[2]], tree = redTree)[[2]] 
  }
  else {
    C = getMeanCoverage(Res, Qred, goodStrains, save = SAVE, baseFilename = baseFilename, tree = redTree)
  }
  output = list(correctFractions, C)
  output
}

### Unfinished function!
findMaxCompatible = function(allMaps, index, K) {
  curPos = lapply(allMaps, function(x) {x[2, x[1,] == index]})
  goodPos = sapply(curPos, function(x) {length(x) > 0})
  curPos[!goodPos] = NA
  if (length(goodPos) > 1) {
    curPos = unlist(curPos)
    curBest = length(goodPos) ### THIS IS WRONG, FIX LATER!
    ### FIND LONGEST PROGRESSION WITH DIFFERENCE K!
  }
  else {
    curBest = 1
  }
  curBest
}

### This function tries to map reads based on a table of unique K-mers
mapReads = function(readTable, uniqueKmers, allKmers, maxError = 1) {
	K = nchar(uniqueKmers[[1]][1])
	N = length(uniqueKmers)
	R = nrow(readTable)
	L = ncol(readTable)
	numKmers = L - K + 1
	map = rep(0, R)
	names(map) = rownames(readTable)
	for (ind in 1:R) {
		curRead = readTable[ind,]
		curKmers = getKmers(curRead, K)
		if (is.null(curKmers)) {
			next
		}
		for (index in 1:N) {
			mappedKmers = which(curKmers %in% uniqueKmers[[index]])
			if (length(mappedKmers) > 0) {
				firstIndex = mappedKmers[1]
				firstMapped = curKmers[firstIndex]
				curSequence = allKmers[[index]]
				mapPos = which(curSequence == firstMapped)
				if (mapPos - firstIndex + numKmers <= length(curSequence)) {
					curMapped = (curSequence[mapPos - firstIndex + (1:numKmers)] == curKmers)
					if (sum(!curMapped) <= maxError * K) {
						map[ind] = index
					}
				}
				break
			}
		}
	}
	map
}

### This function takes in a reference sequence and a read and returns the best place to align a given read to it and a score
alignRead = function(refSequence, read, method = "window") {
	if (method == "window") {
		L0 = length(refSequence)
		L1 = length(read)
		deltaL = L0 - L1 + 1
		scores = rep(NA, deltaL)
		freqs = table(c(refSequence, read))
		freqs = freqs[ALPHABET]
		for (pos in 1:deltaL) {
			Table = table(refSequence[pos + (0:(L1 - 1))], read)
			Table = Table[ALPHABET, ALPHABET]
			scores[pos] = distanceFormula(Table, freqs)
		}
		bestScore = max(scores)
		bestPos = which.max(scores)
	}
	result = list(position = bestPos, score = bestScore)
}

### This function takes in a reference database and a read and returns the vector of supports for each sequence in the database
alignReadToDB = function(DB, read, threshold, ranking = "best", tol = 0.01) {
	M = nrow(DB)
	supports = rep(0, M)
	scores = rep(0, M)
	for (ind in 1:M) {
		curAlign = alignRead(DB[ind,], read)
		scores[ind] = curAlign$score
	}
	bestRefs = which(scores >= max(score) - tol)
	supports[bestRefs] = 1
	supports 
}
