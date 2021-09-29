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
list.files()
all_tree_files = list.files("../python_scripts", pattern = "_tree")
all_tree_files_list = as.list(strsplit(all_tree_files, '\\s+'))
result_path <- readLines(paste("../python_scripts", "/R_information_file", sep = ""))

for (name in all_tree_files_list)
  {
  #tree <- read.tree(textConnection(paste("../python_scripts/",name, sep = "")))
  tree <- read.tree(paste("../python_scripts/",name, sep = ""))
  plotted_tree = ggtree(tree) + geom_tiplab(size=1)
  figure_name = strsplit(name, split = "_")
  figure_name <- sapply(figure_name, function(x) x[1:2])
  figure_name = figure_name[1]
  ggsave(file=paste(result_path,figure_name,"_identity_tree_Rfig",".png",sep = ""), plotted_tree, limitsize = FALSE)
}
plotted_tree
