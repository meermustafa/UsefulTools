# Meer Mustafa
# August 7, 2018

# functions that I have developed over the past 2 years of using R




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

createGuideBasedWindows(inputDF = mycGuides
                        ,inputDFValueToGroup = mycGuides[, ncol(mycGuides)]
                        ,metricToGroupBy = median
                        ,considerGenomicSpanFlag = F
)

# save the metric to the DF
mycGuides = cbind(mycGuides, 
                  mean20guideCSWindows = metricCSWindows)
# head(mycGuides)







# genomic window grouped ----
# group guides by genomic windows 

library(GenomicRanges)

# convert MYC tiling library to data table
library(data.table)
libDT = as.data.table(mycGuides)
# head(libDT)

# create a function with 4 parameters
# output should be chr + genomic intervals, the CS of that interval (possibly gone through a transformation)
groupGuidesByGenomicWindows = function (chr, windowStartValue, windowEndValue, windowSpan, windowShift, transformationAppliedToWindow ) {
  
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
  
  # assemble the df
  windowIntervalsDF = data.frame(chr = chrom, start = start, end = end, width = c(end - start), windowIndex = windowIndex )
  print(head(windowIntervalsDF))
  
  # convert to data table for overlapping
  windowIntervalsDT = as.data.table(windowIntervalsDF)
  
  assign(x = 'windowIntervalsDT', windowIntervalsDT, envir = globalenv())
  
  # after creating windows, overlap with the intervals that have values to be kept
  # use foverlaps with option that desired intervals must be fully contained inside windows
  
  # 1st set keys to overlap by
  setkey(x = windowIntervalsDT, chr, start, end)
  setkey(x = libDT, chr, sgRNA_start, sgRNA_end)
  
  # 2nd overlap
  # compute complete overlaps of guides that fall entirely inside the genomic window -----
  libOverlappedWithWindows = foverlaps(libDT, windowIntervalsDT, type = 'within',
                                       nomatch = 0L # when there are no interval overlaps, output a 0 instead of NA
                                       #, which = T # use which to get the index of x that overlapped with the index of y
  )
  # now this outputs X on the RIGHT, and Y on the LEFT
  
  
  assign(x = 'libOverlappedWithWindows', libOverlappedWithWindows, envir = globalenv())
  
  
}

# call function
groupGuidesByGenomicWindows(chr = mycGuides$chr[1],
                            windowStartValue = mycGuides$sgRNA_start[1]-180, # start a lil bit upstream so that 1st window contains ≥ 1 guide
                            windowEndValue = mycGuides$sgRNA_start[nrow(mycGuides)]+180, # end a lil bit downstream so that last window contains ≥ 1 guide
                            windowSpan = 200, 
                            windowShift = 10)

# head(libOverlappedWithWindows, 50)



# save this large genomic window - to - sgRNA overlapped DF
write.table(libOverlappedWithWindows, paste0(CellLineName,'.',RepName,'.TilingsgRNAs_overlapped_withGenomicWindows.txt'), col.names = T, row.names = F, quote = F, sep = '\t')



# consolidate this overlaps by the window index and apply some function per window index
lastBCname = colnames(libOverlappedWithWindows)[grep("bc",colnames(libOverlappedWithWindows)) [length(grep("bc", colnames( libOverlappedWithWindows)))]]

formula1 = as.formula(paste0('libOverlappedWithWindows$',lastBCname  ,' ~ libOverlappedWithWindows$windowIndex'))

consolidatedLibOverlappedWithWindows = aggregate( formula1, 
                                                  data = libOverlappedWithWindows,
                                                  FUN = mean)





