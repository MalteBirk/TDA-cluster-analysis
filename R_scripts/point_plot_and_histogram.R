# Code from
# https://vbaliga.github.io/verify-that-r-packages-are-installed-and-loaded/

packages = c("ggplot2", "rstudioapi")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)

library("rstudioapi")
library("ggplot2")
library("ggpubr")

setwd(dirname(getActiveDocumentContext()$path))

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
  ggsave(file=paste(result_path,name,".png",sep = ""), g) #saves g
}

