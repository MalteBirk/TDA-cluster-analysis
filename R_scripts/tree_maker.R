# Maybe setenv for xml2
packages <- c("xml2","rvest","tidyverse", "BiocManager", "ggtree", "beeswarm", "ggbeeswarm", "ggplot2")
# Function to check whether package is installed
is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
}
for(x in packages){
  
  if(!is.element(x, installed.packages()[,1]))
    
  {install.packages(x, repos="http://cran.fhcrc.org")
    
  } else {print(paste(x, " library already installed"))}
  
}

library("tidyverse")
library("ggtree")
library("beeswarm")
library("ggbeeswarm")
library("ggplot2")
library("ape")
setwd(getwd())
# Saves all data frames and information about where to write PNG's
paste(getwd())
all_tree_files = list.files("../python_scripts", pattern = "_tree.wrangled")
all_tree_files_list = as.list(strsplit(all_tree_files, '\\s+'))
result_path <- readLines(paste("../python_scripts", "/R_information_file", sep = ""))
exclude_list <- readLines(paste("../python_scripts", "/exclude_list.txt", sep = ""))

for (name in all_tree_files_list)
  {
  figure_name = strsplit(name, split = "_")
  figure_name <- sapply(figure_name, function(x) x[1:2])
  figure_name = figure_name[1]
  tsv_file = paste("../python_scripts/",figure_name,"KnownProducer.tsv", sep = "")
  
  tipcategories = read.csv(tsv_file, 
                           sep = "\t",
                           col.names = c("species", "genus_name", "operon", "producer"),
                           header = FALSE, 
                           stringsAsFactors = FALSE)
  
  META = as.data.frame(tipcategories)
  rownames(META)=META$species
  tree <- read.tree(paste("../python_scripts/",name, sep = ""))
  tree=ape::root(tree,"Paracoccus_limosus_strain_JCM_17370_Scaffold1")
  plotted_tree = ggtree(tree) + geom_treescale()
  
  combined_plot = plotted_tree %<+% META + 
    geom_tiplab(aes(fill = factor(genus_name)),
                size = 0.8,
                linesize = 0.005,
                color = "black", # color for label font
                geom = "label",  # labels not text
                label.padding = unit(0.0001, "lines"), # amount of padding around the labels
                label.size = 0, align = TRUE, linetype = "dashed")
    #theme(legend.position = "none")
    #geom_tippoint(aes(color = operon), 
      #           size=0.8) +
    #geom_tiplab(aes(color = producer, align=T, linetype=NA, 
     #            size=0.8, offset=8, hjust=0.5))
    
    #geom_tiplab2(aes(color = operon), size = 0.3, position = "identity") + 
    #geom_tiplab2(aes(color = producer), size = 0.3, position = "identity")
  #+ # size of label border
   # geom_tippoint(aes(shape = operon, color = operon), size = 0.3, position = "identity") +
    # guides(shape = guide_legend(override.aes = list(size = 5)))
  #combined_plot
  
  heatmap_plot <- gheatmap(combined_plot, META[c(3,4)],width = .08, offset = 0.7  ,
           color = "black", colnames = F) + theme(legend.position = "none")

  ggsave(file=paste(result_path,figure_name,"_heatmap_identity_tree_Rfig",".png",sep = ""), heatmap_plot, limitsize = FALSE, width = 11, height = 12)

  ###################################################### 
  # For trimmed tree
  trimmed_tree <- drop.tip(tree, exclude_list)
  plotted_tree = ggtree(trimmed_tree) + geom_treescale()
  
  combined_plot = plotted_tree %<+% META + 
    geom_tiplab(aes(fill = factor(genus_name)),
                size = 1.5,
                linesize = 0.05,
                color = "black", # color for label font
                geom = "label",  # labels not text
                label.padding = unit(0.0001, "lines"), # amount of padding around the labels
                label.size = 0, align = TRUE, linetype = "dashed")
  heatmap_plot <- gheatmap(combined_plot, META[c(3,4)],width = .08, offset = 2  ,
                           color = "black", colnames = F) + theme(legend.position = "none")
  
  ggsave(file=paste(result_path,figure_name,"_trimmed_heatmap_identity_tree_Rfig",".png",sep = ""), heatmap_plot, limitsize = FALSE, width = 11, height = 12)
  
  }

# Identity plots
all_data_frames = list.files("../python_scripts", pattern = "_identity_operon_data_frame")
all_data_frames = as.list(strsplit(all_data_frames, '\\s+'))

beeswarm_plots = strsplit(result_path, split = "/")
p = sapply(beeswarm_plots, function(x) x[1:4])
path_to_figure_folder = paste(p[1],p[2],p[3],p[4], sep="/")
path_to_figure_folder = paste(path_to_figure_folder,"/beeswarms/", sep = "")
dir.create(path_to_figure_folder)

for (name in all_data_frames)
{
  figure_name = strsplit(name, split = "_")
  figure_name <- sapply(figure_name, function(x) x[1:2])
  figure_name = figure_name[1]
  tipcategories = read.csv(paste("../python_scripts/",name, sep = ""), 
                           sep = "\t",
                           col.names = c("species", "genus_name", "BLAST_identity", "cluster"),
                           header = FALSE, 
                           stringsAsFactors = FALSE)
  identity_frame = as.data.frame(tipcategories)
  #beeswarm(identity_frame["operon"] ~ identity_frame["identity"], pch = 19, corral = "wrap")
  
  beeswarm = ggplot(identity_frame, aes(x = cluster, y = BLAST_identity, color = genus_name)) +
    geom_beeswarm(cex = 0.55) + ggtitle(figure_name)
  ggsave(file=paste(path_to_figure_folder,figure_name,"_beeswarm",".png",sep = ""), beeswarm, limitsize = FALSE, width = 13, height = 8)
  
  #identity = identity_frame["identity"]
  #group = c("Complete")
  #beeswarm(identity)
}
