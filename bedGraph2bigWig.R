# bedGraph to bigWig using RLE

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