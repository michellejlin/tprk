# Creates a ggtree of all the PacBio samples. The script takes the path and will loop through
# all the final fastas to generate the tree. 
# Currently this script does not root (until I figure out a way to automatically do that?)

if (!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager") }
BiocManager::install("treeio")

list.of.packages <- c("treeio","ggtree","stringr","Biostrings","phylobase","pegas","tidyverse","lubridate","ape",
                      "phytools","plyr","phangorn","RColorBrewer","dplyr","optparse")
suppressMessages(invisible(lapply(list.of.packages,library,character.only=T)))

option_list <- list(make_option(c("-d", "--directory"), type="character", default=NULL, help="Specify working directory", metavar="character"));
opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser)

path <- opt$directory

#####

## To run this script manually in R, uncomment the following lines. You do not need to change the preceding lines of path and script.dir,
## but remember to recomment the lines if you want to run the script automatically in the pipeline.
## path refers to the folder your metadata.csv and sequencing files (.fastq) are.

#path <- "/Users/uwvirongs/Documents/Michelle/tprk_pipeline/testing2"

## This script can also be run from the command line.
## Usage: rscript \path\to\PacBio2tree.R -d [directory]

#####

# Works from the fasta files in the folder PacBio_frequencies
setwd(paste(path,"/PacBio_frequencies", sep=""))

# Grabs sample names and PacBio files from the metadata file
metadata <- read.table(paste(path,"/metadata.csv", sep=''), sep=',', header=TRUE)
PacBio_fns <- c(as.character(metadata$PacBio))
sample_names <- c(as.character(metadata$SampleName))

## Function to convert fasta to dataframe
BString2df=function(BString){
  names = names(BString)
  sequences = paste(BString)
  df <- data.frame(names,sequences)
  return(df)
}

## Function to convert fasta to dataframe.  Assumes the two columns are [names,sequences].
df2BString=function(df){
  BString <- BStringSet(df$sequences)
  names(BString) = paste0(df$names)
  return(BString)
}

## Function to remove translated ORFs with stop codons
removeTruncatedORF=function(AAFile){
  AAFile[,2] <- as.character(AAFile[,2])
  for (num in c(1:length(AAFile[,2]))){
    if (!(regexpr("\\*", AAFile[num,2]) == nchar(AAFile[num,2]))){
      AAFile[num,] <- NA}
  }
  return(AAFile)
}

fasta_files <- list()
df_list <- list()
aa_list <- list()
df_aa_list <- list()
df_aa_filtered_list <- list()
allAA <- list()
allAAfilt <- list()

for (i in 1:length(PacBio_fns)) {
  fastafile_name <- paste((substr(PacBio_fns[i],1,nchar(PacBio_fns[i])-6)),".noprimers.filtered.RAD.nolines.fix.fasta",sep="")
  fastafile <- reverseComplement(readDNAStringSet(fastafile_name))
  names(fastafile) <- paste(sample_names[i],"_",names(fastafile),sep="")
  fasta_files <- c(fasta_files,fastafile)
  df_list <- c(df_list,BString2df(fastafile))
  amino_acids <- translate(fastafile, getGeneticCode("11"), no.init.codon=FALSE,
                           if.fuzzy.codon="error")
  aa_list <- c(aa_list, amino_acids)
  df_aa <- BString2df(amino_acids)
  
  ##Filter frequency
  df_aa <- mutate(df_aa,sample=sapply(strsplit(as.character(names),"_"),"[",1),
                       count=as.numeric(sapply(strsplit(as.character(names),"_"),"[",3)))
  df_aa <- mutate(df_aa,percentage=round(count / sum(count)*100,3))
  df_aa_filt <- filter(df_aa,percentage>=0.20)
  df_aa_list[[i]] <- df_aa
  df_aa_filtered_list[[i]] <- df_aa_filt
  if (i == 1) {
    allAA <- df_aa
    allAAfilt <- df_aa_filt
  } else {
    allAA <- rbind(allAA, df_aa)
    allAAfilt <- rbind(allAAfilt, df_aa_filt)
  }
}

allAA_fullORFs_df <- drop_na(removeTruncatedORF(allAA))
allAAfilt_fullORFs_df <- drop_na(removeTruncatedORF(allAAfilt))
write.table(allAAfilt_fullORFs_df,"Table_allAAfilt_fullORFs.tsv",sep='\t',row.names=FALSE)
allAA_fullORFs_BString <- df2BString(allAA_fullORFs_df)
allAAfilt_fullORFs_BString <- df2BString(allAAfilt_fullORFs_df)

AAoutfile <- "Isolates_aa_fullORFs.fasta"
AAoutfilefilt <- "Isolates_aa_filt_fullORFs.fasta"

writeXStringSet(allAA_fullORFs_BString, AAoutfile, append=FALSE,format="fasta")
writeXStringSet(allAAfilt_fullORFs_BString, paste(AAoutfilefilt,sep=""), append=FALSE,format="fasta")

## Need to include our Illumina check right here!  ###

AAaln_outfile <- paste0(substr(AAoutfile,1,nchar(AAoutfile)-5),"aln.fasta")
AAaln_outfile_filt <- paste0(substr(AAoutfilefilt,1,nchar(AAoutfilefilt)-5),"aln.fasta")

#mafft_command <- paste0("mafft --auto ",AAoutfile," > ",AAaln_outfile)
#system(mafft_command)

mafft_command <- paste0("mafft --auto ",AAoutfilefilt," > ",AAaln_outfile_filt)
system(mafft_command)

AAtree_outfile <-paste0(substr(AAoutfile,1,nchar(AAoutfile)-5),"aln.tree.nwk")
AAtree_outfile_filt <-paste0(substr(AAaln_outfile_filt,1,nchar(AAaln_outfile_filt)-5),"aln.tree.nwk")

#fasttree_command <- paste0("fasttree ",AAaln_outfile," > ",AAtree_outfile)
#system(fasttree_command)
fasttree_command <- paste0("fasttree ",AAaln_outfile_filt," > ",AAtree_outfile_filt)
system(fasttree_command)

##Read in FastTree newick file
tree <- read.newick(AAtree_outfile_filt)
mycolors <- colorRampPalette(brewer.pal(name="Dark2", n = 8))(length(sample_names)+4)
#consider...
#tree <- root(tree, which(tree$tip.label == "148B_seq86_220"))
tree <- reorder(tree,"postorder")
tree <- phangorn::midpoint(tree)
tree <- reorder(tree)
##Expects tip labels to contain the sample name followed by a "_" and will parse the sample name accordingly.
d2 = data.frame(taxa=tree$tip.label,sample=sapply(strsplit(tree$tip.label,"_"),"[",1),count=as.numeric(sapply(strsplit(tree$tip.label,"_"),"[",3)))
d2 <- d2 %>% group_by(sample) %>% mutate(percentage=round(count / sum(count)*100,3))
ggtree <- ggtree(tree) %<+% d2 + geom_tippoint(aes(color=sample,size=percentage),alpha=0.8,) + 
  theme(legend.position = c(0.1,0.85),legend.title = element_blank(),legend.key.width=unit(0.2,"cm"),legend.text=element_text(size=10)) + 
  scale_colour_manual(values=mycolors) + 
  geom_treescale(x=0.005,y=150,fontsize=3.5,offset=2) + guides(colour = guide_legend(override.aes = list(size=3))) + 
  scale_size(name="Percentage", breaks=c(0.2,0.5,2,5)) 
ggsave(path=paste(path,"/Figures",sep=""),filename="PacBio_Tree_Filtered.pdf",height = 6, width=4)


# ggtree(tree) %<+% d2 + geom_tippoint(aes(color=sample,size=count),alpha=0.8) + 
#   theme(legend.position = c(0.8,0.2),legend.title = element_blank(),legend.key.width=unit(0.2,"cm"),legend.text=element_text(size=10)) + 
#   scale_colour_brewer(palette="Dark2") + geom_treescale(fontsize=3.5,offset=5) + guides(colour = guide_legend(override.aes = list(size=3)))
# ggsave("Figure_PacBio_Tree.pdf",height = 6, width=4)



######### Sandbox below this line.  Pay no attention!!
# 
# install.packages("vegan")
# library(vegan)
# d2
# 
# diversity(d2[d2$sample == "148B"]$percentage)
# 
# diversity(d2[d2$sample == "148B",]$percentage)
# diversity(d2[d2$sample == "148B2",]$percentage)
# 
# percentage_148B <- d2[d2$sample == "148B",]$percentage
# percentage_148B2 <- d2[d2$sample == "148B2",]$percentage
# 
# View(percentage_148B)
# View(percentage_148B2)
# 
# diversity(d2[d2$sample == "148B",]$percentage)/log10(specnumber(d2[d2$sample == "148B",]$percentage))
# diversity(d2[d2$sample == "148B2",]$percentage)/log10(specnumber(d2[d2$sample == "148B2",]$percentage))
# write.table(d2[d2$percentage >= 0.2,]$taxa,"myfastanames.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
# ?write.table
# length(d2$percentage)
# d2[d2$percentage >= 0.2,]
# write.table(paste0(d2[d2$percentage >= 0.2,]$taxa,"_freq_",d2[d2$percentage >= 0.2,]$percentage),"myfastanamesfreq.txt",quote=FALSE,col.names=FALSE,row.names=FALSE)
# 
# d = data.frame(sample=sapply(strsplit(tree$tip.label,"_"),"[",1))
# d = data.frame(color=sample(c('red', 'green'), 750, replace=T))
# rownames(d) = tree$tip.label
# g2 = phylo4d(g1, d)
# rNodeData <- data.frame(randomTrait = nNodes(g1),
#                         color = sample(c('black'), nNodes(g1), replace=T),
#                         row.names = nodeId(g1, "internal"))
# 
# nodeData(g2) <- rNodeData
# ggtree(g2, aes(color=sample)) + scale_colour_brewer(palette="Dark2")
# 
# tree <- read.newick(AAtree_outfile)
# tree$data
# treetib <- as_tibble(tree)
# treetib$label
# shortsample <- sapply(strsplit(treetib$label,"_"),"[",1)
# treetib <- add_column(treetib,sample = shortsample)
# tree_phylo <- as.phylo(treetib)
# treetib$sample
# 
# 
# ggplot(tree_phylo, aes(x, y)) + geom_tree() + theme_tree() + geom_tippoint(aes(color=sample), size=1.5)
# 
# tree <- read.newick(AAtree_outfile)
# g1 = as(tree,'phylo4')
# g1 = as(tr, 'phylo4')
# shortsample <- sapply(strsplit(tree$tip.label,"_"),"[",1)
