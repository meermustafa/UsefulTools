# Meer Mustafa
# Nov 9, 2018


# Intersecting genomic annotations (5' UTR, CDS exons, introns, 3' UTR) with a supplied table to produce barplots of magnitude of overlap
# i.e. Where do my intervals fall in the genome? Do they overlap with specific sites?




library(TxDb.Hsapiens.UCSC.hg19.knownGene)
## list the contents that are loaded into memory
ls('package:TxDb.Hsapiens.UCSC.hg19.knownGene')
## show the db object that is loaded by calling it’s name
TxDb.Hsapiens.UCSC.hg19.knownGene

txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
seqlevels(txdb)

# If you then wanted to only set Chromosome 1 to be active you could do it like this:
# seqlevels(txdb) = "chr1"

# If you need to reset back to the original seqlevels (i.e. to the seqlevels stored in the db), then set the seqlevels to seqlevels0(txdb).
# seqlevels(txdb) <- seqlevels0(txdb)



# Retrieving data using the select method
columns(txdb)
allColumns = columns(txdb)
keytypes(txdb)

keys <- c("100033416", "100033417", "100033420")
select(txdb, keys = keys, columns="TXNAME", keytype="GENEID")

# select all info
select(txdb, keys=keys, columns=allColumns, keytype="GENEID")


# Retrieving data with select is useful, but sometimes it is more convenient to extract the result as GRanges objects.
# The functions transcripts, exons, and cds return the coordinate information as a GRanges object. 
# As an example, all transcripts present in a TxDb object can be obtained as follows:
# TRANSCRIPTS
GR <- transcripts(txdb)
GR[1:10]

GR <- transcripts(txdb, filter=list(tx_chrom = "chr8"))
GR
length(GR)


# PROMOTERS
# The promoters function computes a GRanges object that spans the promoter region around the transcription start site for the transcripts in a TxDb object.
# The upstream and downstream arguments define the number of bases upstream and downstream from the transcription start site that make up the promoter region.
PR <- promoters(txdb, upstream=2000, downstream=0, filter=list(tx_chrom = "chr8"))
PR

# EXONS & CDS
# The exons and cds functions can also be used in a similar fashion to retrive genomic coordinates for exons and coding sequences.
EX = exons(txdb, filter=list(tx_chrom = "chr8"))
EX[1:4]
CDS = cds(txdb, filter=list(tx_chrom = "chr8"))
CDS[1:10]


introns(txdb, filter=list(tx_chrom = "chr8"))
length(intronsByTranscript(txdb))

head(intronsByTranscript(txdb))







library(GenomicRanges)

a375 = read.table('~/Dropbox/SLab histone/CRISPR_MYC_library/Screens/AllScreens/A375.Rep1.MYC20sgRNAGroupedGuides.txt', header = T)
head(a375)

sdToTake = 4
a375enhancers = a375[ a375[,ncol(a375)] > (sd(a375[,ncol(a375)], na.rm = T) * sdToTake) , ]
a375repressors = a375[ a375[,ncol(a375)] < -(sd(a375[,ncol(a375)], na.rm = T) * sdToTake) , ]


introns = read.table('~/Downloads/refseq_genes_introns')
head(introns)
introns_chr_subset = introns[introns$V1 == 'chr8',]



# convert DF -> Granges
a375enhancers_gr = makeGRangesFromDataFrame(a375enhancers,
                         keep.extra.columns=T,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seqnames", "seqname",
                                          "chromosome", "chrom",
                                          "chr", "chromosome_name",
                                          "seqid", 
                                          'V1'),
                         start.field=c("start",'sgRNA_start', 'V2'),
                         end.field=c("end", "stop", "sgRNA_end", 'V3'),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)


a375repressors_gr = makeGRangesFromDataFrame(a375repressors,
                                            keep.extra.columns=T,
                                            ignore.strand=FALSE,
                                            seqinfo=NULL,
                                            seqnames.field=c("seqnames", "seqname",
                                                             "chromosome", "chrom",
                                                             "chr", "chromosome_name",
                                                             "seqid", 
                                                             'V1'),
                                            start.field=c("start",'sgRNA_start', 'V2'),
                                            end.field=c("end", "stop", "sgRNA_end", 'V3'),
                                            strand.field="strand",
                                            starts.in.df.are.0based=FALSE)


introns_chr_subset_gr = makeGRangesFromDataFrame(introns_chr_subset,
                                                  keep.extra.columns=T,
                                                  ignore.strand=FALSE,
                                                  seqinfo=NULL,
                                                  seqnames.field=c("seqnames", "seqname",
                                                                   "chromosome", "chrom",
                                                                   "chr", "chromosome_name",
                                                                   "seqid", 
                                                                   'V1'),
                                                  start.field=c("start",'sgRNA_start', 'V2'),
                                                  end.field=c("end", "stop", "sgRNA_end", 'V3'),
                                                  strand.field="strand",
                                                  starts.in.df.are.0based=FALSE)


# countOverlaps & findOverlaps return > 1 hit for the query if it overlaps multiple times in subject ******** -------
# compute overlaps of query into subject of 2 GRs
sum(countOverlaps(a375enhancers_gr,
             PR,
             #maxgap = -1L,
             minoverlap = 1L,
             type = c("any"), # type=c("any", "start", "end", "within", "equal"),
             ignore.strand = T
             )
)

countOverlaps(a375enhancers_gr,
              PR,
              #maxgap = -1L,
              minoverlap = 1L,
              type = c("any"), # type=c("any", "start", "end", "within", "equal"),
              ignore.strand = T
)

o = findOverlaps(a375enhancers_gr,
                  PR,
                  #maxgap = -1L,
                  minoverlap = 1L,
                  type = c("any"), # type=c("any", "start", "end", "within", "equal"),
                 #select= 'first', #c("all", "first", "last", "arbitrary"),
                  ignore.strand = T
)

o

# get indices of query and subject from the overlap
queryHits(o)
subjectHits(o)



# -------


# subsetByOverlaps returns the original query subsetted by boolean overlaps in subject

subsetByOverlaps(a375enhancers_gr,
             PR,
             #maxgap = -1L,
             minoverlap = 1L,
             type = c("any"), # type=c("any", "start", "end", "within", "equal"),
             #select = c("all"),
             ignore.strand = T
)

subsetByOverlaps(a375repressors_gr,
                 PR,
                 #maxgap = -1L,
                 minoverlap = 1L,
                 type = c("any"), # type=c("any", "start", "end", "within", "equal"),
                 #select = c("all"),
                 ignore.strand = T
)

# get count of # of overlaps
length(subsetByOverlaps(a375repressors_gr,
                 PR,
                 #maxgap = -1L,
                 minoverlap = 1L,
                 type = c("any"), # type=c("any", "start", "end", "within", "equal"),
                 #select = c("all"),
                 ignore.strand = T
                 )
       )

length(subsetByOverlaps(a375enhancers_gr,
                        introns_chr_subset_gr,
                        #maxgap = -1L,
                        minoverlap = 1L,
                        type = c("any"), # type=c("any", "start", "end", "within", "equal"),
                        #select = c("all"),
                        ignore.strand = T
                        )
       )






# final code w/ loops, functions, etc ------


