# Currently this takes the path and one sample name. It will extract the sample's PacBio and Illumina columns 
# from the allreads.csv file and output a scatterplot comparing PacBio vs. Illumina. Also outputs a zoomed in 
# plot with range 0-10. 

list.of.packages <- c("ggplot2", "grid", "nplr", "dplyr", "scales", "gridExtra", "RColorBrewer", "optparse")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
suppressMessages(invisible(lapply(list.of.packages,library,character.only=T)))


cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

option_list <- list(make_option(c("-p", "--path"), type="character", default=NULL, help="Path to allreads.csv", 
                                metavar="character"),
                    make_option(c("-s", "--sample_name"), type="character", default=NULL, 
                                help="Name of sample for both Illumina and PacBio", metavar="character"));
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

path <- opt$path
sample_name <- opt$sample_name

#####

## To run this script manually in R, uncomment the following lines. You do not need to change the preceding lines of path and sample_name,
## but remember to recomment the lines if you want to run the script automatically in the pipeline.
## path refers to the current working directory where the files and allreads_csv is.
## sample_name refers to the sample you're currently working with for nicer-looking graphs. Make sure this name matches 
## what was used in the allreads.csv (i.e. The sample name must be AS10 if PB_AS10_RelativeFreq is the column.)

#path <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline/paper2_redo/pacbio_v_illumina/"
#sample_name <- "AS10"

## This script can also be run from the command line. 
## Usage: rscript \path\to\Pacbio_v_Illumina_plots.R -p [path] -s [sample_name]

#####
allreads <- paste(path,"/allreads.csv", sep="")
allreads_filtered <- paste(path,"/allreads_filtered.csv", sep="")

alldata <- read.csv(allreads,header=TRUE,sep=",",stringsAsFactors = FALSE)
alldatafilt <- read.csv(allreads_filtered,header=TRUE,sep=",",stringsAsFactors = FALSE)
# Grabs column names relevant to sample.
illrf_col <- (paste("Ill_",sample_name,"_RelativeFreq",sep=""))
pbrf_col <- (paste("PB_",sample_name,"_RelativeFreq",sep=""))
pbcount_col <- (paste("PB_",sample_name,"_Count",sep=""))
illcount_col <- (paste("Ill_",sample_name,"_Count",sep=""))
vs_relative_freq <- select(alldata,Region,Read,illrf_col, pbrf_col)
vs_relative_freq_filt <- select(alldatafilt,Region,Read,illrf_col, pbrf_col)

# Gets r^2 value.
model <- lm(alldata[[pbrf_col]] ~ alldata[[illrf_col]])
anno <- paste("r^2 =",round(summary(model)$adj.r.squared,3))
model_filt <- lm(alldatafilt[[pbrf_col]] ~ alldatafilt[[illrf_col]])
summary(model_filt)
anno_filt <- paste("r^2 =",round(summary(model_filt)$adj.r.squared,3))

# Plots for non-filtered.
plot <- ggplot(alldata,aes(x=alldata[[illrf_col]],y=alldata[[pbrf_col]],color=Region)) + 
  geom_point(size=1.5) + coord_equal() + 
  theme(plot.title = element_text(size = 10,vjust=-2,hjust=0.03)) + 
  geom_abline(linetype="dashed",color="grey") + theme_classic() + scale_colour_manual(values=cbbPalette) + ggtitle(sample_name) +
  xlim(0,100) + ylim(0,100) + annotate("text", x=85, y=2, size=3, label=anno, 3) + 
  theme(legend.position = "bottom",legend.box = "horizontal") + theme(legend.key.width = unit(0.1, "cm")) + 
  xlab(label="Illumina Frequency") + ylab(label="PacBio Frequency") + guides(color = guide_legend(nrow = 1)) + 
  theme(aspect.ratio = 1,axis.title.x = element_text(size=9),axis.title.y = element_text(size=9))

# Plots for filtered.
plot_filt <- ggplot(alldatafilt,aes(x=alldatafilt[[illrf_col]],y=alldatafilt[[pbrf_col]],color=Region)) +
  geom_abline(slope = model_filt$coefficients[2], intercept = model_filt$coefficients[1], col = "slategray1",
              linetype = "solid") + 
  geom_point(size=1.5) + ggtitle(paste(sample_name,"Filtered")) + coord_equal() + 
  theme(plot.title = element_text(size = 10,vjust=-2,hjust=0.03)) + 
  geom_abline(linetype="dashed",color="grey") + theme_classic() + scale_colour_manual(values=cbbPalette) + 
  xlim(0,100) + ylim(0,100) + annotate("text", x=80, y=2, size=3, label=anno_filt, 3) + theme(legend.position = "none") + 
  xlab(label="Illumina Frequency") + ylab(label="PacBio Frequency")  + 
  annotate("text", x=80, y=10, size=3, label=paste("slope = ",round(model_filt$coefficients[2],3),sep=""), 3) + 
  annotate("text", x=80, y=18, size=3, label=paste("intercept = ",round(model_filt$coefficients[1],3),sep=""), 3) + 
  theme(aspect.ratio = 1,axis.title.x = element_text(size=9),axis.title.y = element_text(size=9))

ggsave(paste(sample_name,"Illumina_vs_PacBio_r.pdf",sep="_"),path="Figures",plot=plot_filt,width=6,height=3,units="in")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# Arranges non-filtered and filtered plots side by side.
mylegend<-g_legend(plot)
figure <- grid.arrange(arrangeGrob(plot + theme(legend.position="none"),
                                   plot_filt + theme(legend.position="none"),
                                   nrow=1),
                       mylegend, nrow=2,heights=c(10, 1))
# Saves the plot in the Figures/ subdirectory.
ggsave(paste(sample_name,"Illumina_vs_PacBio.pdf",sep="_"),path="Figures",plot=figure,width=6,height=3,units="in")


# Creates plot zoomed in on range 0-10.
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
# Saves the plot in the Figures/ subdirectory.
ggsave(paste(sample_name,"Illumina_vs_PacBio_lim0-10.pdf",sep="_"),path="Figures",plot=figurelim10,width=6,height=3,units="in")

file.exists("Rplots.pdf")
file.remove("Rplots.pdf")
