# Meer Mustafa
# August 7, 2018


# chosen functions that I have developed over the past few years of using R


# create a function for generating guide windows -----
createGuideBasedWindows = function (inputDF, 
                                    inputDFValueToGroup,
                                    numberOfTotalGuidesToConsiderSurrounding,
                                    metricToGroupBy,
                                    considerGenomicSpanFlag,
                                    maximumGenomicSpanToConsider,
                                    inputDFPositiveStrandStart,
                                    inputDFNegativeStrandStart) {
  
  # empty vector to hold all grouped guide values
  metricCSWindows = c()
  
  # size of guide-windows
  sizeOfGuideWindows = c()
  
  for (guideIndex in 1:nrow(inputDF)) {
    
    # for guides that are at the beginning, set the lowerGuideIndex to the 1st guide in the screen
    if ((guideIndex - 10) < 0) {
      cat('guide index < 0!\n')
      
      # set the initial guide as the first guide in the tiling set
      lowerGuideIndex = 1
    }
    
    
    # if guides are at beginning and don't have 10 guides before than include only up to n guides before the current guide
    if ( (guideIndex - 10) >= 0) {
      cat('guide index > 0!\n')
      
      # set the initial guide as the first guide in the tiling set
      lowerGuideIndex = guideIndex - 10
    }
    
    # set upper guide index as n guides after the last
    upperGuideIndex = guideIndex + 10
    cat( (upperGuideIndex - lowerGuideIndex),'is the # of guides in the window.\n',length(inputDFValueToGroup [lowerGuideIndex:upperGuideIndex] ),'is the # of values considering. \n')
    
    
    # compute the genomic size that spans the n guide-windows
    sizeOfGuideWindows = append(sizeOfGuideWindows,
                                (inputDF$sgRNA_start[upperGuideIndex] - inputDF$sgRNA_start[lowerGuideIndex] )
    )
    
    
    #now we know the genomic span that the guides span, limit the genomic span to n bp e.g. 1000 bp
    if (considerGenomicSpanFlag) {
      
      # only consider guide groups
      if ( (inputDF$sgRNA_start[upperGuideIndex] - inputDF$sgRNA_start[lowerGuideIndex] )  <= 1000 ) {
        
        # group guides by some metric
        metricCSofGuidesInWindow = metricToGroupBy( inputDFValueToGroup[lowerGuideIndex:upperGuideIndex], na.rm = T)
        
        # append this metric to the running vector
        metricCSWindows = append(metricCSWindows, metricCSofGuidesInWindow)
        
      }
      
    } # end of genomic span flag
    
    
    
    # even if genomic span limitation is not implemented, aggregate the guides' CS by a metric and save it
    metricCSofGuidesInWindow = metricToGroupBy( inputDFValueToGroup[lowerGuideIndex:upperGuideIndex], na.rm = T)
    
    metricCSWindows = append(metricCSWindows, metricCSofGuidesInWindow)
    
  }
  
  # assign any output objects to global env
  assign(x = paste('metricCSWindows'),
         value = metricCSWindows, 
         envir = globalenv())
  cat('Saved metricCSWindows!\n')
  
  assign(x = paste('sizeOfGuideWindows'),
         value = sizeOfGuideWindows, 
         envir = globalenv())
  cat('Saved sizeOfGuideWindows!\n')
  
}

# call function and iterate through several metrics (e.g. mean)
#### mean of CS of 20 sgRNA windows


# loop through all the shufflings

# createGuideBasedWindows(inputDF = inputDF
#                         ,inputDFValueToGroup = inputDF[, ncol(inputDF)]
#                         ,metricToGroupBy = median
#                         ,considerGenomicSpanFlag = F
# )






# genomic window grouped ----
# group values by genomic windows 

# create a function with 4 parameters
overlapLibScreenCSWithFeatureSignal_IEoverlap2Ranges = function(rangeDF1, 
                                                                rangeDF2
                                                                , flankingbpFor1stRangeDF
                                                                , typeOfOverlap
                                                                , functionToApplyToOverlappedEntries
                                                                , columnNameInLibToMeasure
                                                                , writeOutFilesFlag = F
) {
  
  # load in dependencies
  library(GenomicRanges)
  library(data.table)
  
  # prevent scientific notation from working
  options(scipen=999)
  
  # assign range 1 as library and range 2 as feature signal
  inputDF = rangeDF1
  windowIntervalsDF = rangeDF2
  
  # convert to data table for overlapping
  libDT = as.data.table(inputDF)
  windowIntervalsDT = as.data.table(windowIntervalsDF)
  
  # add window index
  
  
  
  # flank the ranges if desired
  Prime5flank = (flankingbpFor1stRangeDF / 2) - 1
  Prime3flank = flankingbpFor1stRangeDF / 2
  
  libDT$sgRNA_start = libDT$sgRNA_start - Prime5flank
  libDT$sgRNA_end = libDT$sgRNA_end + Prime3flank
  
  # after creating windows, overlap with the intervals that have values to be kept
  # use foverlaps with option that desired intervals must be fully contained inside windows
  # 1st set keys to overlap by
  setkey(x = windowIntervalsDT, chr, start, end)
  setkey(x = libDT, chr, sgRNA_start, sgRNA_end)
  
  # compute complete overlaps of guides that fall entirely inside the genomic window
  libOverlappedWithWindows = foverlaps(
    windowIntervalsDT, 
    libDT, 
    type = 'within',
    nomatch = 0L # when there are no interval overlaps, output a 0 instead of NA****
    #, which = T # use which to get the index of x that overlapped with the index of y
  )
  # now this outputs X on the RIGHT, and Y on the LEFT
  cat("Samples overlapped with windows DF:...\n")
  print(head(libOverlappedWithWindows))
  print(tail(libOverlappedWithWindows))
  
  # save file to env + write out
  assign(x = paste0(as.list( sys.call() )[2],'_libOverlappedWithWindows_',as.list( sys.call() )[6],'bpGenomicWindows_',as.list( sys.call() )[8],'Transformation'),
         libOverlappedWithWindows,
         envir = globalenv())
  
  if(writeOutFilesFlag) {
    write.table(libOverlappedWithWindows,
                paste0(as.list( sys.call() )[2],'.TilingsgRNAs_overlapped_with',as.list( sys.call() )[6],'bpGenomicWindows_',as.list( sys.call() )[8],'Transformation.txt'),
                col.names = T, row.names = F, quote = F, sep = '\t')
  }
  
  # now consolidate the guides that fall into each the window index by applying some function to all samples inside the window index
  # apply the function onto the last BC (CS)
  lastBCname = colnames(libOverlappedWithWindows)[grep(columnNameInLibToMeasure,colnames(libOverlappedWithWindows)) 
                                                  [length(grep(columnNameInLibToMeasure, colnames( libOverlappedWithWindows)))]]
  
  allBCnames = colnames(libOverlappedWithWindows)[grep(columnNameInLibToMeasure,colnames(libOverlappedWithWindows))]
  
  formula1 = as.formula(paste0('libOverlappedWithWindows$', lastBCname ,' ~ libOverlappedWithWindows$UniqueID'))
  
  consolidatedLibOverlappedWithWindows = aggregate( formula1, 
                                                    data = libOverlappedWithWindows,
                                                    FUN = functionToApplyToOverlappedEntries)
  
  # add back the window coordinates so that it doesn't just contain the window index
  consolidatedLibOverlappedWithWindows = merge(consolidatedLibOverlappedWithWindows, 
                                               libDT, 
                                               by.x = 'libOverlappedWithWindows$UniqueID', 
                                               by.y = 'UniqueID', 
                                               all = F)
  
  cat("Consolidated DF:...\n")
  print(head(consolidatedLibOverlappedWithWindows))
  print(tail(consolidatedLibOverlappedWithWindows))
  
  # save file to env + write out
  assign(x = paste0(as.list( sys.call() )[2],'_Consolidated_CS_',as.list( sys.call() )[7],'FeatureMeasurement',as.list( sys.call() )[6],'Transformation'), 
         consolidatedLibOverlappedWithWindows,
         envir = globalenv())
  cat("Assigned",paste0(as.list( sys.call() )[2],'_Consolidated_CS_',as.list( sys.call() )[7],'FeatureMeasurement',as.list( sys.call() )[6],'Transformation'),"\n")
  
  if(writeOutFilesFlag) {
    write.table(consolidatedLibOverlappedWithWindows,
                paste0(as.list( sys.call() )[2],'_Consolidated_CS_',as.list( sys.call() )[6],'bpGenomicWindows_',as.list( sys.call() )[8],'Transformation.txt'),
                col.names = T, row.names = F, quote = F, sep = '\t')
  }
  
  
  
}




# genomic window grouped ----
# group guides by genomic windows 


# create a function with 4 parameters
# output should be chr + genomic intervals, the CS of that interval (possibly gone through a transformation)
groupGuidesByGenomicWindows = function (inputDF,
                                        chr, 
                                        windowStartValue,
                                        windowEndValue,
                                        windowSpan,
                                        windowShift,
                                        typeOfOverlap, 
                                        functionToApplyToOverlappedEntries
                                        ,flankingbpFor1stRangeDF
                                        ,columnNameInLibToMeasure
                                        ,columnNameToConsolidateBy
                                        ,writeOutFilesFlag=F
                                        ) {
  
  library(GenomicRanges)
  library(data.table)
  
  # prevent scientific notation from working
  options(scipen=999)
  
  # define window intervals
  start = seq(from = windowStartValue, to = windowEndValue, by = windowShift)
  end = seq(from = windowStartValue + windowSpan, to = windowEndValue + windowSpan, by = windowShift)
  
  # create chr names for overlapping later
  chrom = rep(paste(chr), times = length(start))
  
  # create a window index number using seq_len
  windowIndex = seq_len(length(start))
  
  # assemble the df of tiling windows
  windowIntervalsDF = data.frame(chr = chrom, start = start, end = end, width = c(end - start), windowIndex = windowIndex )
  print(head(windowIntervalsDF))
  
  # convert to data table for overlapping
  libDT = as.data.table(inputDF)
  windowIntervalsDT = as.data.table(windowIntervalsDF)
  
  # save the DF version of the windows
  windowIntervalsDF = as.data.frame(windowIntervalsDT)
  
  assign(x = paste0( 'windowIntervalsDF',as.list( sys.call() )[6],'bpGenomicWindows'), windowIntervalsDF, envir = globalenv())
  
  # after creating windows, overlap with the intervals that have values to be kept
  # use foverlaps with option that desired intervals must be fully contained inside windows
  # 1st set keys to overlap by
  setkey(x = windowIntervalsDT, chr, start, end)
  setkey(x = libDT, chr, start, end)
  
  # 2nd overlap
  # compute complete overlaps of guides that fall entirely inside the genomic window -----
  libOverlappedWithWindows = foverlaps(libDT, 
                                       windowIntervalsDT, 
                                       type = typeOfOverlap,
                                       nomatch = 0L # when there are no interval overlaps, output a 0 instead of NA
                                       #, which = T # use which to get the index of x that overlapped with the index of y
  )
  # now this outputs X on the RIGHT, and Y on the LEFT
  
  cat("Samples overlapped with windows DF:...\n")
  print(head(libOverlappedWithWindows))
  print(tail(libOverlappedWithWindows))
  
  # save file to env + write out
  assign(x = paste0(as.list( sys.call() )[2],'_libOverlappedWithWindows_',as.list( sys.call() )[6],'bpGenomicWindows_',as.list( sys.call() )[8],'Transformation'),
         libOverlappedWithWindows,
         envir = globalenv())
  
  if(writeOutFilesFlag) {
    write.table(libOverlappedWithWindows,
                paste0(as.list( sys.call() )[2],'.TilingsgRNAs_overlapped_with',as.list( sys.call() )[6],'bpGenomicWindows_',as.list( sys.call() )[8],'Transformation.txt'),
                col.names = T, row.names = F, quote = F, sep = '\t')
  }
  
  # now consolidate the guides that fall into each the window index by applying some function to all samples inside the window index
  # apply the function onto the last BC (CS)
  lastBCname = colnames(libOverlappedWithWindows)[grep(columnNameInLibToMeasure,colnames(libOverlappedWithWindows)) 
                                                  [length(grep(columnNameInLibToMeasure, colnames( libOverlappedWithWindows)))]]
  lastBCname
  allBCnames = colnames(libOverlappedWithWindows)[grep(columnNameInLibToMeasure,colnames(libOverlappedWithWindows))]
  
  formula1 = as.formula(paste0('libOverlappedWithWindows$', lastBCname ,' ~ libOverlappedWithWindows$',columnNameToConsolidateBy))
  formula1
  
  consolidatedLibOverlappedWithWindows = aggregate( formula1, 
                                                    data = libOverlappedWithWindows,
                                                    FUN = functionToApplyToOverlappedEntries)
  
  cat("Consolidated",lastBCname," based on",columnNameToConsolidateBy,"...\n")
  print(head(consolidatedLibOverlappedWithWindows))
  print(tail(consolidatedLibOverlappedWithWindows))
  
  # add back the window coordinates so that it doesn't just contain the window index
  consolidatedLibOverlappedWithWindows = merge(consolidatedLibOverlappedWithWindows, 
                                               windowIntervalsDT, 
                                               by.x = paste0('libOverlappedWithWindows$',columnNameToConsolidateBy), 
                                               by.y = columnNameToConsolidateBy, 
                                               all = F)
  
  cat("Added information to consolidated DF:...\n")
  print(head(consolidatedLibOverlappedWithWindows))
  print(tail(consolidatedLibOverlappedWithWindows))
  
  # save file to env + write out
  assign(x = paste0(as.list( sys.call() )[2],'_Consolidated_CS_',as.list( sys.call() )[7],'FeatureMeasurement',as.list( sys.call() )[6],'Transformation'), 
         consolidatedLibOverlappedWithWindows,
         envir = globalenv())
  cat("Assigned",paste0(as.list( sys.call() )[2],'_Consolidated_CS_',as.list( sys.call() )[7],'FeatureMeasurement',as.list( sys.call() )[6],'Transformation'),"\n")
  
  if(writeOutFilesFlag) {
    write.table(consolidatedLibOverlappedWithWindows,
                paste0(as.list( sys.call() )[2],'_Consolidated_CS_',as.list( sys.call() )[6],'bpGenomicWindows_',as.list( sys.call() )[8],'Transformation.txt'),
                col.names = T, row.names = F, quote = F, sep = '\t')
  }
  
  return(consolidatedLibOverlappedWithWindows)
  
}

windowSpan = 500
windowShift = 10

# call function
groupGuidesByGenomicWindows(inputDF = input,
                            chr = input$chr[1],
                            windowStartValue = input$start[1] - 180, # start a lil bit upstream so that 1st window contains ≥ 1 guide
                            windowEndValue = input$start[nrow(input)] + 180, # end a lil bit downstream so that last window contains ≥ 1 guide
                            windowSpan = 500, 
                            windowShift = 10,
                            transformationAppliedToWindow = median)


ls(pattern = 'Consolidated_CS_')
ls(pattern = '_libOverlappedWithWindows_')



# convert bdg to library size normalized bdg file using the bam file -----
# meant to be parallel with the ls *bam
# July 8, 2018


#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# specify character to match in bdg file names (e.g. H3K27ac or DHS)
featureNamesToAdd = args[2]


# function to convert .perbp BedGraph -> Sushi input (chr, start, end, pileup), normalize pileups
BedGraphToNormalizedBedGraph_usingBAM = function(featureStringPatternToMatch, positionOfNameInInputBedGraph, positionOfFeatureInInputBedGraph, nameToAddToOutputNormalizedBDG ) {
  
  cat('Working on feature',featureStringPatternToMatch,'...\n')
  
  # extract all file names in the directory
  AllFeatureFiles = list.files(path = '.', pattern = featureStringPatternToMatch)
  
  NameOfRootCellLine = sapply(strsplit(args[1], split = '.', fixed = T), function(x) x[positionOfNameInInputBedGraph] )
  cat('Reading in cell line',NameOfRootCellLine,'...\n')
  
  # pattern match the matching bam and bdg files using the cell line name
  BedGraphfile = AllFeatureFiles[grep('.perbp$', AllFeatureFiles)]
  
  # subset all bdgs to specific bdg
  currentBedGraphfilenames = BedGraphfile[grep(NameOfRootCellLine, BedGraphfile)]
  
  # if more than one bedgraph belonging to that cell line (e.g. MYCcov.perbp & ZNF479cov.perbp), loop through all
  for (currentBedGraphfilename in  currentBedGraphfilenames) {
    
    cat('Current cell line regex matched file is',currentBedGraphfilename,'. About to read in...\n')
    
    # read in bedgraph (chr, position, pileup)
    currentBedGraph = read.table(currentBedGraphfilename)
    
    # index the bam file
    system(paste("samtools index ",args[1],sep=" "))
    
    # normalize using BAM library size, extracted using samtools
    CurrentBAMLibrarySizeInMillionsOfReads =  ( as.numeric(system(paste("samtools view -F 0x04 -c ",
                                                                        args[1],
                                                                        sep=" "),
                                                                  intern = T))
                                                / 1000000 ) # divide by 1e6 to get Reads Per Million
    cat(CurrentBAMLibrarySizeInMillionsOfReads,'M reads...\n')
    
    # divide the reads per bp (R/bp) in the bedgraph by the reads per million (RPM) to get the RPM/bp
    currentBedGraph$V3 = round( (currentBedGraph$V3 / CurrentBAMLibrarySizeInMillionsOfReads), digits = 4 )
    
    # write out the normalized bedgraph file
    write.table(currentBedGraph, file = paste0(currentBedGraphfilename,'.',nameToAddToOutputNormalizedBDG), row.names = F, col.names = F, quote = F)
    
  }
  
  
}

# loop through all features to load into envir
# for (currentFeatureName in featureNamesToAdd) {
#   
#   # invoke the previously defined function that converts BedGraph -> Sushi DF
#   BedGraphToNormalizedBedGraph_usingBAM(featureStringPatternToMatch = currentFeatureName,
#                                         positionOfNameInInputBedGraph = 1,
#                                         positionOfFeatureInInputBedGraph = 4,
#                                         nameToAddToOutputNormalizedBDG = 'normalized_RPMbp'
#   )
#   
# }





# convert bedgraph (3 column) file to bigWig
a = read.table('exampleNormalizedRPMbedgraph.bedgraph')

Bdg2BigWig = function(inputBedGraph, BdgColumnWithSignal# , outputBigWigName
                      ,writeOutFlag = F) {
  # begin the new DF with the first row of the input BDG
  OutputBigWig = inputBedGraph[1, ]
  
  for (row in 2:nrow(inputBedGraph)) {
    
    # progress report
    if (row %% 100000 == 0) {
      percentDone = row/nrow(inputBedGraph)*100
      cat(percentDone,'% done ... \n')
    }
    
    # # if signal at row = 1 IS NOT different than row = 2 just keep that initial row
    # if ( inputBedGraph [row, BdgColumnWithSignal] == inputBedGraph [(row-1), BdgColumnWithSignal]) {
    #   
    #   # pass this if loop
    #   #next()
    #   
    # }
    
    # if signal at row = 1 IS different than row = 2 add that new row to the 
    if ( inputBedGraph [row, BdgColumnWithSignal] != inputBedGraph [(row-1), BdgColumnWithSignal]) {
      
      # add the current row that is diff. than the previous one
      OutputBigWig = rbind(OutputBigWig, inputBedGraph [row, ])
      
    }
  }
  
  # if flag = T, write out file
  if (writeOutFlag) {
    write.table(OutputBigWig, 
                paste0(as.list(sys.call()[2]),'.bigWig'),
                quote = F, col.names = F, row.names = F, sep = '\t'
                )
  }
  
  return(OutputBigWig)
}

b = Bdg2BigWig(a, 3, writeOutFlag = T)
head(b);nrow(b)



# bedgraph to bigwig using RLE ----
#### much faster than Bdg2BigWig*********

Bdg2BigWig_viaRLE = function (inputBDG, columnWithSignal, writeOutFlag = F) {
  
  # subset BDG to rows in which the difference between successive signal differences are not 0
  OutputBigWig = inputBDG[c(1, 1 + which(diff(inputBDG[,columnWithSignal]) != 0) ) , ]
  
  assign(x = paste0(deparse(substitute(inputBDG)),'bigWig'),
         value = OutputBigWig,
         envir = globalenv())
  cat('Assigned', paste0(deparse(substitute(inputBDG)),'bigWig'),'\n')
  
  if(writeOutFlag == 't') {
    write.table(OutputBigWig, 
                paste0(as.list(sys.call()[2]),'.bigWig'),
                quote = F, col.names = F, row.names = F, sep = '\t'
    )
    cat('Wrote out', paste0(as.list(sys.call()[2]),'.bigWig'),'\n')
  }
  
}



# uppercase only the first letter of the character string ----
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  output = as.vector(paste0(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" "))
  return(output)
  
}

name <- c("zip code", "state", "final count")

sapply(name, simpleCap)





