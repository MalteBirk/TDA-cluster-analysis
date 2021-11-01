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
#library("beeswarm")
library("ggbeeswarm")
library("ggplot2")
setwd(getwd())
# Saves all data frames and information about where to write PNG's
paste(getwd())
all_tree_files = list.files("../python_scripts", pattern = "_tree.wrangled")
all_tree_files_list = as.list(strsplit(all_tree_files, '\\s+'))
result_path <- readLines(paste("../python_scripts", "/R_information_file", sep = ""))

for (name in all_tree_files_list)
  {
  figure_name = strsplit(name, split = "_")
  figure_name <- sapply(figure_name, function(x) x[1:2])
  figure_name = figure_name[1]
  tsv_file = paste("../python_scripts/",figure_name,"KnownProducer.tsv", sep = "")
  
  tipcategories = read.csv(tsv_file, 
                           sep = "\t",
                           col.names = c("species", "genus_name", "operon", "TDA-producer"),
                           header = FALSE, 
                           stringsAsFactors = FALSE)
  
  META = as.data.frame(tipcategories)
  rownames(META)=META$species
  tree <- read.tree(paste("../python_scripts/",name, sep = ""))
  plotted_tree = ggtree(tree) + geom_treescale()
  
  combined_plot = plotted_tree %<+% META + 
    geom_tiplab(aes(fill = factor(genus_name)),
                size = 0.8,
                linesize = 0.005,
                color = "black", # color for label font
                geom = "label",  # labels not text
                label.padding = unit(0.0001, "lines"), # amount of padding around the labels
                label.size = 0)
  #+ # size of label border
   # geom_tippoint(aes(shape = operon, color = operon), size = 0.3, position = "identity") +
    # guides(shape = guide_legend(override.aes = list(size = 5)))

  heatmap_plot <- gheatmap(combined_plot, META[c(3,4)],width = .1, offset = 0.8,
           colnames=T, legend_title="Cluster", color = "black") +
    scale_x_ggtree()
  
  ggsave(file=paste(result_path,figure_name,"heatmap_identity_tree_Rfig",".png",sep = ""), heatmap_plot, limitsize = FALSE, width = 8, height = 8)
}

# Identity plots
all_data_frames = list.files("../python_scripts", pattern = "_identity_operon_data_frame")
all_data_frames = as.list(strsplit(all_data_frames, '\\s+'))
for (name in all_data_frames)
{
  tipcategories = read.csv(paste("../python_scripts/",name, sep = ""), 
                           sep = "\t",
                           col.names = c("species", "genus_name", "identity", "operon"),
                           header = FALSE, 
                           stringsAsFactors = FALSE)
  identity_frame = as.data.frame(tipcategories)
  ggplot(identity_frame, aes(x = operon, y = identity)) +
    geom_beeswarm()
  
  #identity = identity_frame["identity"]
  #group = c("Complete")
  #beeswarm(identity)
}
