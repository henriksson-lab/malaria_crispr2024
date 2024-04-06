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
  w <- which(str_starts(s,"ID"))
  if(length(w)>0){
    str_replace_all(s[w[1]],"ID=","")
  } else {
    ""
  }
}
getDESC <- function(s){
  s <- str_split(s,";")[[1]]
  w <- which(str_starts(s,"description"))
  if(length(w)>0){
    str_replace_all(s[w[1]],"description=","")
  } else {
    ""
  }
}
getNAME <- function(s){
  s <- str_split(s,";")[[1]]
  w <- which(str_starts(s,"Name"))
  if(length(w)>0){
    str_replace_all(s[w[1]],"Name=","")
  } else {
    ""
  }
}
gff <- read.table("PlasmoDB-67_PbergheiANKA.gff.gz",sep="\t", quote="")
gff <- gff[str_detect(gff$V9,"description"),]$V9
gff <- gff[which(str_detect(gff,"ID"))]
#length(gff)
geneinfo <- data.frame(
  gene=sapply(gff, getID),
  geneDesc=sapply(gff, getDESC),
  geneName=sapply(gff, getNAME)
)
geneinfo$geneDesc <- str_replace_all(geneinfo$geneDesc,"%2C","")
geneinfo <- unique(geneinfo)
geneinfo
#rownames(geneinfo) <- NULL
#head(geneinfo)
#sum(geneinfo$gene=="PBANKA_1360100")

#foo <- geneinfo[str_detect(rownames(geneinfo),"PBANKA_1360100"),]
#dim(foo)
#rownames(foo) <- NULL
#foo

print("======= reading data to be cached in memory ================ ")

all_samplemeta <- readRDS("samplemeta.rds")
all_grstats <- readRDS("grstats.rds")
all_timecourses <- readRDS("timecourses.rds")
all_coverage_stat <- readRDS("coverage_stat.rds")

print("========== global done ================")
