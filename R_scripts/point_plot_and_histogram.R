# In order to work with R 3.5
packageurl<-"https://cran.r-project.org/src/contrib/Archive/nloptr/nloptr_1.2.1.tar.gz"
if(!require(rstudioapi)){install.packages(packageurl, repos=NULL, type="source")}
packages <- c("ggplot2", "ggpubr")
# Function to check whether package is installed
is.installed <- function(mypkg){
  is.element(mypkg, installed.packages()[,1])
}
for(x in packages){
  
  if(!is.element(x, installed.packages()[,1]))
    
  {install.packages(x, repos="http://cran.fhcrc.org")
    
  } else {print(paste(x, " library already installed"))}
  
}

library("ggplot2")
library("ggpubr")
setwd(getwd())
# Saves all data frames and information about where to write PNG's
all_data_frames = list.files("../python_scripts", pattern = "data_frame")
all_data_frames_list = as.list(strsplit(all_data_frames, '\\s+'))
result_path <- readLines(paste("../python_scripts", "/R_information_file", sep = ""))

# Make GGplots of identities, and lengths.
for (name in all_data_frames_list){
  data_frame = read.table(paste("../python_scripts/",name,sep = ""),
                          sep="\t")
  colnames(data_frame)=c("FILENAME","ID","LEN", "GENUS", "SPECIES")
  
  point_plot = ggplot(data_frame, aes(ID, LEN)) +
    geom_point(aes(colour = GENUS)) + scale_color_brewer(palette="Dark2")
  
  ID_hist_plot = ggplot(data_frame, aes(x=ID, colour=GENUS)) +
    geom_histogram(fill="white", alpha=0.5, position="identity") + 
    scale_color_brewer(palette="Dark2")
  
  LEN_hist_plot = ggplot(data_frame, aes(x=LEN, colour=GENUS)) +
    geom_histogram(fill="white", alpha=0.5, position="identity") + 
    scale_color_brewer(palette="Dark2")
  
  g <- ggarrange(
    point_plot,                # First row with line plot
    # Second row with box and dot plots
    ggarrange(ID_hist_plot, LEN_hist_plot, ncol = 2, labels = c("B", "C")), 
    nrow = 2, 
    labels = "A"       # Label of the line plot
  ) 
  figure_name = strsplit(name, split = "_")
  figure_name <- sapply(figure_name, function(x) x[1:2])
  figure_name = paste(figure_name[1],figure_name[2],sep = "_")
  ggsave(file=paste(result_path,figure_name,"_identity_and_length_Rfig",".png",sep = ""), g) #saves g
}

