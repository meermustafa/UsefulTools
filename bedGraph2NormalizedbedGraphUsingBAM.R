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