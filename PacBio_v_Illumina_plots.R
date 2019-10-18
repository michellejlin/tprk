list.of.packages <- c("ggplot2", "grid", "nplr", "dplyr", "scales", "gridExtra", "RColorBrewer", "optparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

library(ggplot2)
library(grid)
library(nplr)
library(dplyr)
library(scales)
library(gridExtra)
library(RColorBrewer)
library("optparse")

cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

option_list <- list(make_option(c("-p", "--path"), type="character", default=NULL, help="Path to allreads.csv", 
                                metavar="character"),
                    make_option(c("-s", "--sample_name"), type="character", default=NULL, 
                                help="Name of sample for both Illumina and PacBio", metavar="character"));
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

##Finding files in directory
path <- opt$path
#path <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline/testing"
allreads <- paste(path,"/allreads.csv", sep="")
allreads_filtered <- paste(path,"/allreads_filtered.csv", sep="")
sample_name <- opt$sample_name
#sample_name <- "AS10"

alldata <- read.csv(allreads,header=TRUE,sep=",",stringsAsFactors = FALSE)
alldatafilt <- read.csv(allreads_filtered,header=TRUE,sep=",",stringsAsFactors = FALSE)
illrf_col <- (paste("Ill_",sample_name,"_RelativeFreq",sep=""))
pbrf_col <- (paste("PB_",sample_name,"_RelativeFreq",sep=""))
pbcount_col <- (paste("PB_",sample_name,"_Count",sep=""))
illcount_col <- (paste("Ill_",sample_name,"_Count",sep=""))

vs_relative_freq <- select(alldata,Region,Read,illrf_col, pbrf_col)
vs_relative_freq_filt <- select(alldatafilt,Region,Read,illrf_col, pbrf_col)

model <- lm(alldata[[pbrf_col]] ~ alldata[[illrf_col]])
anno <- paste("r^2 =",round(summary(model)$adj.r.squared,3))
model_filt <- lm(alldatafilt[[pbrf_col]] ~ alldatafilt[[illrf_col]])
anno_filt <- paste("r^2 =",round(summary(model_filt)$adj.r.squared,3))

plot <- ggplot(alldata,aes(x=alldata[[illrf_col]],y=alldata[[pbrf_col]],color=Region)) + 
  geom_point(size=2) + 
  geom_abline(linetype="dashed",color="grey") + theme_classic() + scale_colour_manual(values=cbbPalette) + ggtitle(sample_name) +
  xlim(0,100) + ylim(0,100) + annotate("text", x=85, y=2, size=3, label=anno, 3) + 
  theme(legend.position = "bottom",legend.box = "horizontal") + theme(legend.key.width = unit(0.1, "cm")) + 
  xlab(label="Illumina Frequency") + ylab(label="PacBio Frequency") + guides(color = guide_legend(nrow = 1))

plot_filt <- ggplot(alldatafilt,aes(x=alldatafilt[[illrf_col]],y=alldatafilt[[pbrf_col]],color=Region)) + 
  geom_point(size=2) + ggtitle(paste(sample_name,"Filtered")) + 
  geom_abline(linetype="dashed",color="grey") + theme_classic() + scale_colour_manual(values=cbbPalette) + 
  xlim(0,100) + ylim(0,100) + annotate("text", x=85, y=2, size=3, label=anno_filt, 3) + theme(legend.position = "none") + 
  xlab(label="Illumina Frequency") + ylab(label="PacBio Frequency") 

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

mylegend<-g_legend(plot)
figure <- grid.arrange(arrangeGrob(plot + theme(legend.position="none"),
                                   plot_filt + theme(legend.position="none"),
                                   nrow=1),
                       mylegend, nrow=2,heights=c(10, 1))
ggsave(paste(sample_name,"Illumina_vs_PacBio.pdf",sep="_"),path="Figures",plot=figure,width=6,height=3,units="in")


lim10plot <- ggplot(alldata, aes(x=alldata[[illrf_col]],y=alldata[[pbrf_col]],color=Region)) + xlim(0,10) + ylim(0,10) + 
  geom_point(size=2) + xlab(label="Illumina Frequency") + ylab(label="PacBio Frequency") + 
  geom_abline(linetype="dashed", color = "grey") + theme_classic() + ggtitle(sample_name) + 
  scale_color_brewer(palette="Set1") + annotate("text", x=8, y=0.5, label=anno, 3)
lim10plot_filt <- ggplot(alldatafilt, aes(x=alldatafilt[[illrf_col]],y=alldatafilt[[pbrf_col]],color=Region)) +  xlim(0,10) + ylim(0,10) + 
  geom_point(size=2) + xlab(label="Illumina Frequency") + ylab(label="PacBio Frequency") + 
  geom_abline(linetype="dashed", color = "grey") + theme_classic() + ggtitle(paste(sample_name,"Filtered")) + 
  scale_color_brewer(palette="Set1") + annotate("text", x=8, y=0.5, label=anno_filt, 3)
figurelim10 <- grid.arrange(arrangeGrob(lim10plot + theme(legend.position="none"),
                                   lim10plot_filt + theme(legend.position="none"), nrow=1),
                       mylegend, nrow=2,heights=c(10, 1))
ggsave(paste(sample_name,"Illumina_vs_PacBio_lim0-10.pdf",sep="_"),path="Figures",plot=figurelim10,width=6,height=3,units="in")
file.exists("Rplots.pdf")
file.remove("Rplots.pdf")
