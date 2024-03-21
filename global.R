if(FALSE){
  install.packages("shiny")
}

library(shiny)

library(stringr)
library(ggplot2)
library(patchwork)
library(sqldf)
library(reshape2)
library(cowplot)


print("======= reading data to be cached in memory ================ ")


getID <- function(s){
  s <- str_split(s,";")[[1]]
#  print(s)
  w <- which(str_starts(s,"ID"))
 # print(w)
  if(length(w)>0){
    str_replace_all(s[w[1]],"ID=","")
  } else {
    ""
  }
}
getDESC <- function(s){
  s <- str_split(s,";")[[1]]
  #print(s)
  w <- which(str_starts(s,"description"))
  #print(w)
  if(length(w)>0){
    str_replace_all(s[w[1]],"description=","")
  } else {
    ""
  }
}
getNAME <- function(s){
  s <- str_split(s,";")[[1]]
  #print(s)
  w <- which(str_starts(s,"Name"))
  #print(w)
  if(length(w)>0){
    str_replace_all(s[w[1]],"Name=","")
  } else {
    ""
  }
}
gff <- read.table("PlasmoDB-67_PbergheiANKA.gff.gz",sep="\t")
gff <- gff[str_detect(gff$V9,"description"),]$V9
geneinfo <- data.frame(
  gene=sapply(gff, getID),
  geneDesc=sapply(gff, getDESC),
  geneName=sapply(gff, getNAME)
)
geneinfo$geneDesc <- str_replace_all(geneinfo$geneDesc,"%2C","")
geneinfo <- unique(geneinfo)
geneinfo


print("======= reading data to be cached in memory ================ ")

all_samplemeta <- readRDS("samplemeta.rds")

all_grstats <- readRDS("grstats.rds")

all_timecourses <- readRDS("timecourses.rds")

all_coverage_stat <- readRDS("coverage_stat.rds")

print("========== global done ================")
