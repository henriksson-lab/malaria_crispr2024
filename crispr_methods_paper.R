library(ggtext)


################################################################################
####################### Range of fc per gRNA ###################################
################################################################################

listpools <- c("cr_2023march_screen")
for(curpool in listpools){
  print(curpool)
  
  grstats <- all_grstats[[curpool]]$stats_per_grna$`NP BL6`

  grstats$i <- str_split_fixed(grstats$grna,"g",2)[,2]
  
  grstats$xlabcol <- "gray" #Other etc
  grstats$xlabcol[grstats$genecat == "Essential"] <- "red"
  grstats$xlabcol[grstats$genecat == "Dispensable"] <- "darkgreen"
  grstats$xlabcol[grstats$genecat == "Slow"] <- "blue"
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

}


################################################################################
################# Relative abundance over time, per gRNA #######################
################################################################################



listpools <- c("cr_2023march_screen")
for(curpool in listpools){
  print(curpool)
  
  grstats <- timecourses[[curpool]]$`Count/AllCount`

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

  ggsave(sprintf("/corgi/websites/malariascreenviewer/plots_crispr/gridpanel_per_gene %s.pdf",curpool), width = 10, height = 10)
  
}
