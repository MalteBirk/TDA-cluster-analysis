# Maybe setenv for xml2
packages <- c("xml2","rvest","tidyverse", "BiocManager", "ggtree")
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
  tsv_file = paste("../python_scripts/",figure_name,".tsv", sep = "")
  
  tipcategories = read.csv(tsv_file, 
                           sep = "\t",
                           col.names = c("species", "genus_name"), 
                           header = FALSE, 
                           stringsAsFactors = FALSE)
  
  dd = as.data.frame(tipcategories)
  tree <- read.tree(paste("../python_scripts/",name, sep = ""))
  plotted_tree = ggtree(tree) + geom_treescale()
  
  combined_plot = plotted_tree %<+% dd + 
    geom_tiplab(aes(fill = factor(genus_name)),
                size = 0.8,
                linesize = 0.005,
                color = "black", # color for label font
                geom = "label",  # labels not text
                label.padding = unit(0.0001, "lines"), # amount of padding around the labels
                label.size = 0)# size of label border
  
  ggsave(file=paste(result_path,figure_name,"_identity_tree_Rfig",".png",sep = ""), combined_plot, limitsize = FALSE)
}

all_tree_files = list.files("../python_scripts", pattern = "_tree.tree")
all_tree_files_list = as.list(strsplit(all_tree_files, '\\s+'))
result_path <- readLines(paste("../python_scripts", "/R_information_file", sep = ""))

for (name in all_tree_files_list)
{
  figure_name = strsplit(name, split = "_")
  figure_name <- sapply(figure_name, function(x) x[1:2])
  figure_name = figure_name[1]
  tsv_file = paste("../python_scripts/",figure_name,".tsv", sep = "")
  tree <- read.tree(paste("../python_scripts/",name, sep = ""))
  plotted_tree = ggtree(tree) + geom_treescale() + geom_tiplab(size = 0.8)
  
  ggsave(file=paste(result_path,figure_name,"_no_colours_identity_tree_Rfig",".png",sep = ""), plotted_tree, limitsize = FALSE)
}
