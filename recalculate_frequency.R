library("optparse")

option_list <- list(make_option(c("-d", "--directory"), type="character", default=NULL, help="Specify working directory", metavar="character"),
                    make_option(c("-f", "--filename"), type="character", default=NULL, help="Specify file to recalculate frequency.", 
                                metavar="character"));
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

path <- opt$directory
filename <- (opt$filename)
#path <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline/AS_files/"

allreads_filtered <- read.table(opt$filename, sep=',', header=TRUE)
allreads_filtered1 <- allreads_filtered
numsamples <- (length(colnames(allreads_filtered)) - 2) / 4
for (sample in c(1:numsamples)){
  freqcol <- (sample * 2) + 1
  countcol <- (sample * 2) + 2
  for (region in unique(allreads_filtered1$Region)){
    allreads_filtered1[which(allreads_filtered1[,1]==region),freqcol] <- 
      allreads_filtered1[which(allreads_filtered1[,1]==region),countcol] / sum(allreads_filtered1[which(allreads_filtered1[,1]==region),countcol],na.rm=TRUE) * 100
  }
  allreads_filtered1[,freqcol] <- trimws(format(allreads_filtered1[,freqcol], digits = 4, nsmall = 4))
}

#out_filename <- substr(filename,1,nchar(filename)-4)
#allreads_out <- paste(path,"/",out_filename,".csv",sep='')
write.csv(allreads_filtered1,file=filename,row.names=FALSE,quote=FALSE)
