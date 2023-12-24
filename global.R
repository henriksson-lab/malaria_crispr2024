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

all_samplemeta <- readRDS("samplemeta.rds")

all_grstats <- readRDS("grstats.rds")

all_timecourses <- readRDS("timecourses.rds")

all_coverage_stat <- readRDS("coverage_stat.rds")

print("========== global done ================")
