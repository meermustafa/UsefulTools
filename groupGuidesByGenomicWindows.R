# Meer Mustafa

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
