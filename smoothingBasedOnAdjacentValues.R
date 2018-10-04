# function for generating item smoothed values
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
  
  HalfNumberOfTotalGuidesToConsiderSurrounding = numberOfTotalGuidesToConsiderSurrounding / 2
  
  for (guideIndex in 1:nrow(inputDF)) {
    
    # for guides that are at the beginning, set the lowerGuideIndex to the 1st guide in the screen
    if ((guideIndex - HalfNumberOfTotalGuidesToConsiderSurrounding) < 0) {
      cat('guide index < 0!\n')
      
      # set the initial guide as the first guide in the tiling set
      lowerGuideIndex = 1
    }
    
    
    # if guides are at beginning and don't have 10 guides before than include only up to n guides before the current guide
    if ( (guideIndex - HalfNumberOfTotalGuidesToConsiderSurrounding) >= 0) {
      cat('guide index > 0!\n')
      
      # set the initial guide as the first guide in the tiling set
      lowerGuideIndex = guideIndex - HalfNumberOfTotalGuidesToConsiderSurrounding
    }
    
    # set upper guide index as n guides after the last
    upperGuideIndex = guideIndex + HalfNumberOfTotalGuidesToConsiderSurrounding
    
    
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
    metricCSofGuidesInWindow = metricToGroupBy( inputDFValueToGroup[ lowerGuideIndex:upperGuideIndex ], na.rm = T)
    
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
# createGuideBasedWindows(inputDF = inputDF
#                         ,inputDFValueToGroup = inputDF[ , 1]
#                         ,numberOfTotalGuidesToConsiderSurrounding = 20
#                         ,metricToGroupBy = median
#                         ,considerGenomicSpanFlag = F
#                         )
