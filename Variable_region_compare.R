list.of.packages <- c("ggplot2", "grid", "nplr", "plyr", "dplyr", "scales", "gridExtra", "RColorBrewer", "optparse","randomcoloR", "cowplot",
                      "tidyr", "tibble")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
suppressMessages(invisible(lapply(list.of.packages,library,character.only=T)))

option_list <- list(make_option(c("-p", "--path"), type="character", default=NULL, help="Path to allreads.csv", 
                                metavar="character"));
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

##Finding files in directory
path <- opt$path
#path <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline/testing"
allreads <- paste(path,"/allreads_filtered.csv", sep="")
allreads_filtered <- paste(path,"/allreads_filtered.csv", sep="")
metadata_csv <- paste(path,"/metadata.csv", sep="")
metadata = read.csv(metadata_csv)
sample_names = as.character(metadata$SampleName)
setwd(path)
  
alldata <- read.csv(allreads,header=TRUE,sep=",",stringsAsFactors = FALSE)

for (i in 1:(length(sample_names) - 1)) {
  for (j in (i+1):(length(sample_names))) {
    print(paste("Generating figure for ",sample_names[i]," and ",sample_names[j],"...",sep=""))
    rfcol <- paste("Ill_",sample_names[i],"_RelativeFreq",sep = "")
    rfcol2 <- paste("Ill_",sample_names[j],"_RelativeFreq",sep = "")
    commondfIllumina <- select(alldata,Region,Read,rfcol,rfcol2)
    commondfIllumina <- filter(commondfIllumina,!((commondfIllumina[[rfcol]] == 0) & (commondfIllumina[[rfcol2]] == 0)))
    sortedIllumina <- commondfIllumina[order(commondfIllumina$Region,-commondfIllumina[[rfcol]],-commondfIllumina[[rfcol2]]),]
    sortedIllumina <- gather(sortedIllumina,rfcol,rfcol2,key="Sample",value="Frequency")
    sortedIllumina$Sample[sortedIllumina$Sample == rfcol] <- sample_names[i]
    sortedIllumina$Sample[sortedIllumina$Sample == rfcol2] <- sample_names[j]
    
    myColors <- distinctColorPalette(1068)
    names(myColors) <- levels(sortedIllumina$Read)
    colScale <- scale_colour_manual(name = NULL, guide = FALSE, values = myColors)
    h <- ggplot(sortedIllumina[which(sortedIllumina$Frequency>0),]) + geom_point(aes(y = Frequency, x = Sample, color=Read)) + geom_line(aes(y = Frequency, x = Sample, group=Read, color=Read)) +  
      facet_wrap(~Region,nrow=1) + theme(legend.position = "none") + theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank()) +
      theme(axis.text.x=element_text(angle=45,hjust=1))
    h1 <- h + colScale
    suppressMessages(ggsave(paste(sample_names[i],"_vs_",sample_names[j],"_VariableRegions_DotLine.pdf",sep=""),
           path="Figures/Variable_Region_Comparisons",plot=h1,width=5,height=4,units="in"))
  }
}