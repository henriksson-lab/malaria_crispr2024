library(ggtext)
library(ggplot2)
library(grid)
library(egg)

all_samplemeta <- readRDS("/corgi/websites/malariascreenviewer/samplemeta.rds")
#all_grstats <- readRDS("grstats.rds")
#all_timecourses <- readRDS("timecourses.rds")
#all_coverage_stat <- readRDS("coverage_stat.rds")

################################################################################
############## Fig ????  Barchart, relative abundance, pools ###################
################################################################################

for(current_pool in c("cr_2024march_half1","cr_2024march_p1","cr_2024march_p12","cr_2024march_p2")){
  #current_pool <- "cr_2024march_half1"
  samplemeta <- all_samplemeta[[current_pool]]
  coverage_stat <- all_coverage_stat[[current_pool]]
  
  ligpools <- sqldf::sqldf("select ligationwell, sum(cnt) as cnt from coverage_stat group by ligationwell")
  ligpools <- ligpools[ligpools$ligationwell!="spikein",]
  ligpools$frac <- ligpools$cnt/sum(ligpools$cnt)
  #ligpools$ligationwell <- str_sub(ligpools$ligationwell,2)
  ggplot(ligpools, aes(ligationwell, frac*100)) + geom_bar(stat="identity") + 
    geom_hline(yintercept=100/nrow(ligpools), color="blue")+
    xlab("") + ylab("Fraction (%)")+
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  ggsave(sprintf("/corgi/websites/malariascreenviewer/plots_crispr/fraction_per_ligpool %s.pdf",current_pool), width = 3, height = 3)  
}

current_pool <- "cr_2023march_screen"

################################################################################
####################### Fig xXX. Range of fc per gRNA ############################
################################################################################

#listpools <- c("cr_2024march_half1")
curpool <- "cr_2024march_half1"
curpool <- "cr_2024march_p2"
#for(curpool in listpools){
  print(curpool)
  
  grstats <- all_grstats[[curpool]]$stats_per_grna$`NP BL6`

#  grstats <- grstats[grstats$gene !="PBANKA_0914900",] #bad gene
  
  grstats$i <- str_split_fixed(grstats$grna,"g",2)[,2]
  
  grstats$xlabcol <- "gray" #Other etc
  grstats$xlabcol[grstats$genecat == "Essential"] <- "red"
  grstats$xlabcol[grstats$genecat == "Dispensable"] <- "darkgreen"
  grstats$xlabcol[grstats$genecat == "Slow growers"] <- "blue"
  #Disp vara grön, Ess vara röd

  grstats$x.label <- paste("<span style = 'color: ",grstats$xlabcol,";'>",grstats$gene,"</span>", sep = "")

  ggplot(grstats, aes(x=x.label, y=fc, group=i)) + 
    geom_pointrange(data=grstats, aes(ymin=fc-sd, ymax=fc+sd,color=gene), 
                    position=position_dodge(.9)) + 
    theme_bw()+
    theme(legend.position = "none")+
    theme(axis.text.y = ggtext::element_markdown())+
    coord_flip() +
    xlab("")+
    ylab("RGR") 
  ggsave(sprintf("/corgi/websites/malariascreenviewer/plots_crispr/fc_pergrna %s.pdf",curpool), width = 10, height = 10)

#}



################################################################################
####################### Fig 4. Range of fc per gRNA ############################ -------- alternative
################################################################################
# 
# #listpools <- c("cr_2024march_half1")
# curpool <- "cr_2024march_half1"
# curpool <- "cr_2024march_p2"
# #for(curpool in listpools){
# print(curpool)
# 
# grstats <- all_grstats[[curpool]]$stats_per_grna$`NP BL6`
# 
# #  grstats <- grstats[grstats$gene !="PBANKA_0914900",] #bad gene
# 
# grstats$i <- str_split_fixed(grstats$grna,"g",2)[,2]
# 
# grstats$xlabcol <- "gray" #Other etc
# grstats$xlabcol[grstats$genecat == "Essential"] <- "red"
# grstats$xlabcol[grstats$genecat == "Dispensable"] <- "darkgreen"
# grstats$xlabcol[grstats$genecat == "Slow growers"] <- "blue"
# #Disp vara grön, Ess vara röd
# 
# grstats$x.label <- paste("<span style = 'color: ",grstats$xlabcol,";'>",grstats$gene,"</span>", sep = "")
# 
# ggplot(grstats, aes(x=gene, y=fc, group=i)) + 
#   geom_pointrange(data=grstats, aes(ymin=fc-sd, ymax=fc+sd,color=xlabcol, group=gene), 
#                   position=position_dodge(.9)) + 
#   theme_bw()+
#   theme(legend.position = "none")+
#   theme(axis.text.y = ggtext::element_markdown())+
#   coord_flip() +
#   xlab("")+
#   ylab("RGR") 
# ggsave(sprintf("/corgi/websites/malariascreenviewer/plots_crispr/fc_pergrna %s.pdf",curpool), width = 10, height = 10)
#   

################################################################################
################# Fig 5. Relative abundance over time, per gRNA ################
################################################################################



#listpools <- c("cr_2023march_screen")
#for(curpool in listpools){
  curpool <- "cr_2024march_half1"
  
  print(curpool)
  
  grstats <- all_timecourses[[curpool]]$`Count/AllCount`

  #Figure out where in a grid all plots should go
  geneinfo <- unique(grstats[,c("gene","genecat")])
  geneinfo <- geneinfo[order(geneinfo$genecat),]
  geneinfo$i <- 1:nrow(geneinfo)
  geneinfo <- merge(geneinfo,sqldf::sqldf("select min(i) as mini, genecat from geneinfo group by genecat"))
  geneinfo$grid_y <- geneinfo$i-geneinfo$mini+1
  geneinfo$grid_x <- as.integer(factor(geneinfo$genecat))

  #Figure out coloring based on gene category  
  geneinfo$xlabcol <- "gray" #Other etc
  geneinfo$xlabcol[geneinfo$genecat == "Essential"] <- "red"
  geneinfo$xlabcol[geneinfo$genecat == "Dispensable"] <- "darkgreen"
  geneinfo$xlabcol[geneinfo$genecat == "Slow"] <- "blue"
  
  pushViewport(viewport(layout = grid.layout(max(geneinfo$grid_y), max(geneinfo$grid_x))))
  for(curgene in unique(grstats$gene)){
    
    p1 <- ggplot(grstats[grstats$gene==curgene,], aes(x=day, y=y, color=grna, group=paste(grna,mouse_ref))) + 
      geom_line()+
      theme_bw()+
      theme(legend.position = "none")+
      theme(axis.title.x = element_text(colour = geneinfo$xlabcol[geneinfo$gene==curgene]))+
      xlab(curgene)+
      ylab("Rel.count")
    p1
    
    print(p1, vp = viewport(
      layout.pos.row = geneinfo$grid_y[geneinfo$gene==curgene], 
      layout.pos.col = geneinfo$grid_x[geneinfo$gene==curgene]))

    666
  }

#  ggsave(sprintf("/corgi/websites/malariascreenviewer/plots_crispr/gridpanel_per_gene %s.pdf",curpool), width = 10, height = 10)
  
#}

  
  
  
  
  
  
  


##################################################################################
# Fig xxxx. Normalized abundance over time, per gRNA; all of them, dirty format ##
##################################################################################


for(curpool in c("cr_2024march_half1","cr_2024march_p1","cr_2024march_p12","cr_2024march_p2")){
  #curpool <- "cr_2024march_half1"
  print(curpool)
  grstats <- all_timecourses[[curpool]]$`Count/ControlCount`
  list_all_plot <- list()
  for(curgene in unique(grstats$gene)){
    
    p1 <- ggplot(grstats[grstats$gene==curgene,], aes(x=day, y=y, linetype=grna, color=grna, group=paste(grna,mouse_ref))) + 
      geom_line()+
      scale_color_manual(values = c("black","darkgray"))+
      theme_bw()+
      theme(legend.position = "none")+
      theme(axis.title.x = element_text(colour = geneinfo$xlabcol[geneinfo$gene==curgene]))+
      xlab(curgene)+
      ylab("Count/Disp. gene count")
    p1
    list_all_plot[[curgene]] <- p1
  }
  ptot <- egg::ggarrange(plots=list_all_plot, ncol=1)
  ggsave(sprintf("/corgi/websites/malariascreenviewer/plots_crispr/all_lineplot %s.pdf",curpool), width = 3, height = 2*length(list_all_plot), plot = ptot, limitsize=FALSE)
}  


  
  
  
################################################################################
################# Fig xxx. Volcano plots #######################################
################################################################################

  
for(current_pool in c("cr_2024march_half1","cr_2024march_p1","cr_2024march_p12","cr_2024march_p2")){
  #current_pool <- "cr_2024march_half1"
  print(current_pool)
  grstats <- all_grstats[[current_pool]]
  thecond <- "NP BL6"
  theplot <- grstats$volcano[[thecond]]
  
  toplot$y <- 1/toplot$sd
  yname <- paste("inverse s.d.")#,thecond)
  toplot$genecat <- factor(toplot$genecat, levels=c("Dispensable","Essential","Slow growers","Other"))
  ggplot(toplot, aes(fc, y, label=gene, color=genecat)) + 
    geom_point() + 
    xlab("RGR") + #xlab(paste("FC",thecond)) + 
    ylab(yname) +
    scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "turquoise3"))+ #"Dispensible","Essential","Slow growers","Other"
    theme_bw() + 
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  ggsave(sprintf("/corgi/websites/malariascreenviewer/plots_crispr/volcano %s.pdf",current_pool), width = 5, height = 4)
  
  #https://sape.inf.usi.ch/quick-reference/ggplot2/colour
  
}  

  
  
  
  
  
  
  
################################################################################
################# Fig 5 foo ####################################################
################################################################################
  
  
  


if(FALSE){
  
  #BiocManager::install("ComplexHeatmap")
  
  
  
  
  
  
  #Johan would it be possible to 
  #plot how the genes included in the ½ plate with all known pheno act in each pool, without the unknowns and with the spike in controls highlighted in some way?
  
  
  names(all_timecourses)
  
  current_pool <- "cr_2024march_half1"
  grstats_units <- "Count/ControlCount"
  grstats <- all_timecourses[[current_pool]]
  grstats <- grstats[[grstats_units]]  
  print(head(grstats))
  
  
  plotTC(
    grstats,
    grstats_avg_grna, input$grstats_avg_mouse, input$grstats_avg_genotype, input$grstats_avg_treatment,
    input$grstats_gene, input$grstats_colorby
  )
  
  ######## gene name from plasmodb with gene symbol (ideally in hover)
  
  
  # sen ellen väljer 3 gener => fig 1, time courses
  
  #fig4: half-plate 
  #fig5: the rest of the plates  
}


  