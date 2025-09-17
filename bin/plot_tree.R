#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

# A simple function to parse named arguments --key value
parse_args <- function(args) {
  arg_list <- list()
  i <- 1
  while (i <= length(args)) {
    if (startsWith(args[i], "--")) {
      key <- substring(args[i], 3)
      if ((i + 1) <= length(args) && !startsWith(args[i + 1], "--")) {
        arg_list[[key]] <- args[i + 1]
        i <- i + 2
      } else {
        arg_list[[key]] <- TRUE  # flag without value
        i <- i + 1
      }
    } else {
      i <- i + 1
    }
  }
  return(arg_list)
}

opt <- parse_args(args)

# Access parsed arguments
tree_file <- opt$tree_file
tree_type <- opt$tree_type
sample_id <- opt$sample_id

if (is.null(tree_file) || is.null(tree_type) || is.null(sample_id)) {
  stop("Missing required arguments: --tree_file, --tree_type, and --sample_id must be provided")
}

cat("Sample ID:", sample_id, "\n")
cat("Tree file:", tree_file, "\n")
cat("Tree type:", tree_type, "\n")

library(tidyverse)

library(ggtree)


library(seqinr)
library(Biostrings)
library(ggmsa)
library(ape)


tree <- read.tree(tree_file)

node_data <- as_tibble(tree)

pattern <- ifelse(grepl("core",tree_type),"core","ns5b")

highlighted <- node_data[grepl(pattern,node_data$label),]$label
tip_df <- data.frame(label = tree$tip.label,
                     highlight = ifelse(tree$tip.label == highlighted, "yes", "no"))

subtypes <- sapply(strsplit(highlighted, "\\|"), tail, 1)


clade_nodes = NULL


for (clade in subtypes){
  
  if(!grepl("^7",clade)){ #do not highlight 7 otherwise it would highlight the whole thing since the tree is rooted on 7.
  
  labels = node_data$label[grep(paste0("^",clade), node_data$label)]
  
  clade_nodes = c(clade_nodes, getMRCA(tree,c(labels,highlighted[grep(clade,highlighted)])))
  }

}

if(!is.null(clade_nodes)){
  p=ggtree(tree) %<+% tip_df + geom_tiplab(size = 2,aes(color = highlight, fontface = ifelse(highlight == "yes", "bold", "plain"),
                                                      size = ifelse(highlight == "yes", 6, 3)),show.legend = FALSE) +
    scale_color_manual(values = c("yes" = "red", "no" = "black")) + 
  
    geom_hilight(node = c(clade_nodes), fill = "red") 
}else{
  p=ggtree(tree) %<+% tip_df + geom_tiplab(size = 2,aes(color = highlight, fontface = ifelse(highlight == "yes", "bold", "plain"),
                                                        size = ifelse(highlight == "yes", 6, 3)),show.legend = FALSE) +
    scale_color_manual(values = c("yes" = "red", "no" = "black")) 
}

if(grepl("subtype",tree_type)){
  ggsave(paste0(sample_id,"_",tree_type,"_tree.png"), p, width = 8, height = 10, units = "in", dpi = 300)
}else{
  ggsave(paste0(sample_id,"_",tree_type,"_tree.png"), p, width = 8, height = 18, units = "in", dpi = 300)

}


