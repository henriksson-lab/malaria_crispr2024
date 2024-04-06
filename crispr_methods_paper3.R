library(ggtext)
library(ggplot2)
library(grid)
library(egg)

pools_renamed <- list(
  cr_2024march_half1="48p",
  cr_2024march_p1="96p",
  cr_2024march_p2="p2",
  cr_2024march_p12="192p"
)

all_samplemeta <- readRDS("/corgi/websites/malariascreenviewer/samplemeta.rds")
#all_grstats <- readRDS("grstats.rds")
#all_timecourses <- readRDS("timecourses.rds")
#all_coverage_stat <- readRDS("coverage_stat.rds")

################################################################################
############## Fig ????  Barchart, relative abundance, pools ###################
################################################################################

for(current_pool in c("cr_2024march_half1","cr_2024march_p1","cr_2024march_p12")){ #"cr_2024march_p2"
  #current_pool <- "cr_2024march_half1"
  samplemeta <- all_samplemeta[[current_pool]]
  coverage_stat <- all_coverage_stat[[current_pool]]

  
  #Per ligation pool
  if(FALSE){
    ligpools <- sqldf::sqldf("select ligationwell, sum(cnt) as cnt from coverage_stat group by ligationwell")
    ligpools <- ligpools[ligpools$ligationwell!="spikein",]
    ligpools$frac <- ligpools$cnt/sum(ligpools$cnt)
    #ligpools$ligationwell <- str_sub(ligpools$ligationwell,2)
    ggplot(ligpools, aes(ligationwell, frac*100)) + geom_bar(stat="identity") + 
      geom_hline(yintercept=100/nrow(ligpools), color="blue")+
      xlab("") + ylab("Fraction (%)")+
      theme_bw() + 
      theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  }
  
  
  #Per grna
  ligpools <- sqldf::sqldf("select gene, ligationwell, sum(cnt) as cnt from coverage_stat group by gene, ligationwell")
  ligpools <- ligpools[ligpools$ligationwell!="spikein",]
  ligpools$frac <- ligpools$cnt/sum(ligpools$cnt)
  ggplot(ligpools, aes(ligationwell, frac*100, fill=gene)) + geom_bar(stat="identity", position = "stack") + 
    geom_hline(yintercept=100/length(unique(coverage_stat$ligationwell)), color="blue")+
    xlab("") + ylab("Fraction (%)")+
    theme_bw() + 
    theme(legend.position = "none")+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
    

  ggsave(sprintf("/corgi/websites/malariascreenviewer/plots_crispr/fraction_per_ligpool %s.pdf",current_pool), width = 3, height = 3)  
}

#current_pool <- "cr_2023march_screen"

################################################################################
####################### Fig xXX. Scatter plot of sgRNA 1 vs 2 ##################
################################################################################

# listpools <- c(
#   "cr_2024march_half1",
#   "cr_2024march_p1",
#   "cr_2024march_p2",
#   "cr_2024march_p12"
# )
# curpool <- "cr_2024march_half1"
# curpool <- "cr_2024march_p2"
# curpool <- "cr_2024march_p12"
# #for(curpool in listpools){
#   print(curpool)
#   
#   grstats <- all_grstats[[curpool]]$stats_per_grna$`NP BL6`
# 
# #  grstats <- grstats[grstats$gene !="PBANKA_0914900",] #bad gene
#   
#   grstats$i <- str_split_fixed(grstats$grna,"g",2)[,2]
#   
#   grstats$xlabcol <- "gray" #Other etc
#   grstats$xlabcol[grstats$genecat == "Essential"] <- "red"
#   grstats$xlabcol[grstats$genecat == "Dispensable"] <- "darkgreen"
#   grstats$xlabcol[grstats$genecat == "Slow growers"] <- "blue"
#   #Disp vara grön, Ess vara röd
# 
#   grstats$x.label <- paste("<span style = 'color: ",grstats$xlabcol,";'>",grstats$gene,"</span>", sep = "")
# 
#   ggplot(grstats, aes(x=x.label, y=fc, group=i)) + 
#     geom_pointrange(data=grstats, aes(ymin=fc-sd, ymax=fc+sd,color=gene), 
#                     position=position_dodge(.9)) + 
#     theme_bw()+
#     theme(legend.position = "none")+
#     theme(axis.text.y = ggtext::element_markdown())+
#     coord_flip() +
#     xlab("")+
#     ylab("RGR") +
#     ylim(-5,5)
#   ggsave(sprintf("/corgi/websites/malariascreenviewer/plots_crispr/fc_pergrna %s.pdf",curpool), width = 10, height = 10)
# 
# #}
# ######## Alternative way
# ggplot(grstats, aes(x=x.label, y=fc, group=gene)) + 
#   geom_line(aes(color=genecat))+#data=grstats, aes(ymin=fc-sd, ymax=fc+sd,color=gene), 
#               #    position=position_dodge(.9)) + 
#   theme_bw()+
#   theme(legend.position = "none")+
#   theme(axis.text.y = ggtext::element_markdown())+
#   coord_flip() +
#   xlab("")+
#   ylab("RGR")# +
# #    ylim(-5,5)
#   

listpools <- c(
  "cr_2024march_half1",
  "cr_2024march_p1",
  "cr_2024march_p2"
#  "cr_2024march_p12"
)
listplot <- list()
for(curpool in listpools){
  print(curpool)
  grstats <- all_grstats[[curpool]]$stats_per_grna$`NP BL6`
  
  g1 <- grstats[str_ends(grstats$grna,"1"),]
  g2 <- grstats[str_ends(grstats$grna,"2"),]
  toplot <- merge(
    data.frame(gene=g1$gene, fc1=g1$fc, genecat=g1$genecat),
    data.frame(gene=g2$gene, fc2=g2$fc)
  )
  toplot$genecat <- factor(toplot$genecat, levels=c("Dispensable","Essential","Slow growers","Other"))
  onep <- ggplot(toplot) + 
    xlab(paste("FC sgRNA #1", curpool))+
    ylab("FC sgRNA #2")+
    theme_bw()+
    theme(legend.position = "none")+
    geom_smooth(method = "lm", aes(fc1,fc2),color="black")+
    geom_point(aes(fc1,fc2,color=genecat),size=3)  +
    scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "turquoise3")) #"Dispensible","Essential","Slow growers","Other"
  listplot[[curpool]] <- onep
}  
totp <- egg::ggarrange(plots=listplot)
totp
ggsave(plot=totp, "/corgi/websites/malariascreenviewer/plots_crispr/scatterplot_sgrna_fc.pdf", width = 10, height = 10)



################################################################################
############### Fig xxx Comparison of screen FCs, scatter plots ################
################################################################################
  

listpools <- c(
#  "cr_2024march_half1",
  "cr_2024march_p1",
  "cr_2024march_p2",
  "cr_2024march_p12"
)
listplot <- list()
for(curpool in listpools){
  print(curpool)

  grstats1 <- all_grstats[["cr_2024march_half1"]]$volcano$`NP BL6`
  grstats2 <- all_grstats[[curpool]]$volcano$`NP BL6`
  
  toplot <- merge(
    data.frame(
      gene=grstats1$gene,
      fc1=grstats1$fc,
      genecat=grstats1$genecat
    ),  
    data.frame(
      gene=grstats2$gene,
      fc2=grstats2$fc
    )
  )
  toplot$genecat <- factor(toplot$genecat, levels=c("Dispensable","Essential","Slow growers","Other"))
  onep <- ggplot(toplot, aes(fc1,fc2,color=genecat)) + 
    geom_point() + 
    xlab("FC half-pool") +
    ylab(paste("FC",curpool))+
    theme_bw()+
    theme(legend.position = "none")+
    geom_smooth(method = "lm", aes(fc1,fc2),color="black")+
    geom_point(aes(fc1,fc2,color=genecat),size=3)  +
    scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "turquoise3")) #"Dispensible","Essential","Slow growers","Other"
  listplot[[curpool]] <- onep
  
  
  cor(toplot$fc1,toplot$fc2) ####
}  
totp <- egg::ggarrange(plots=listplot)
totp

ggsave(plot=totp, "/corgi/websites/malariascreenviewer/plots_crispr/scatterplot_pool_reproducibility.pdf", width = 10, height = 10)




################################################################################
################## Fig xxx Composite analysis of FC ############################
################################################################################


listpools <- c(
  "cr_2024march_half1",
  "cr_2024march_p1",
  "cr_2024march_p2",
  "cr_2024march_p12"
)

#### Merge all pools
poolstats <- NULL
for(curpool in listpools){
  print(curpool)
  grstats <- all_grstats[[curpool]]$volcano$`NP BL6`
  grstats$pool <- curpool
  grstats$fc <- grstats$fc - mean(grstats$fc[grstats$genecat=="Dispensable"])
  #grstats$fc <- grstats$fc / -mean(grstats$fc[grstats$genecat=="Essential"])  #this will not work for p2
  poolstats <- rbind(poolstats, grstats)
}

avgpool <- sqldf::sqldf("select count(*) as cnt, sum(sd*sd) as totvar, avg(fc) as fc, gene, genecat from poolstats group by gene")
avgpool$sd <- sqrt(avgpool$totvar/avgpool$cnt)

#### Distribution of Essential
fc_ess <- avgpool$fc[avgpool$genecat=="Essential"]
mean_ess <- mean(fc_ess)
sd_ess <- sd(fc_ess)
dist_ess <- data.frame(
  x=seq(from=-4,to=2, by=0.01)
)
dist_ess$p <- dnorm(dist_ess$x,mean=mean_ess, sd=sd_ess)
dist_ess$type<-"Essential"


#### Distribution of Dispensable
fc_disp <- avgpool$fc[avgpool$genecat=="Dispensable"]
mean_disp <- mean(fc_disp)
sd_disp <- sd(fc_disp)
dist_disp <- data.frame(
  x=seq(from=-4,to=2, by=0.01)
)
dist_disp$p <- dnorm(dist_disp$x,mean=mean_disp, sd=sd_disp)
dist_disp$type<-"Dispensable"


#### Produce plot of all pools
avgpool$genecat <- factor(avgpool$genecat, levels=c("Dispensable","Essential","Slow growers","Other"))
ggplot(avgpool, aes(fc, 1/sd, color=genecat)) + 
  geom_point()+
  geom_line(data=rbind(dist_disp, dist_ess),aes(x,p,color=type))+
  xlab("FC") +
  theme_bw()+
  theme(legend.position = "none") +
  scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "turquoise3")) #"Dispensable","Essential","Slow growers","Other"

ggsave("/corgi/websites/malariascreenviewer/plots_crispr/composite_histograms.pdf", width = 10, height = 10)



################### ROC curve
#library(verification)
tovery <- avgpool[avgpool$genecat %in% c("Dispensable", "Essential"),]
tovery$obs <- tovery$genecat=="Dispensable"
#roc.plot(
#  tovery$obs,
#  tovery$fc
#)
ggplot(tovery, aes(d = obs, m = fc)) + geom_roc() 

ggsave("/corgi/websites/malariascreenviewer/plots_crispr/composite_roc.pdf", width = 10, height = 10)


################################################################################
################## Fig xxx ROC curve for each screen ###########################
################################################################################

#install.packages("plotROC")
library(plotROC)

listpools <- c(
  "cr_2024march_half1",
  "cr_2024march_p1",
  #"cr_2024march_p2", #no essentials!
  "cr_2024march_p12"
)
listplot <- list()
for(curpool in listpools){
  print(curpool)
  
  grstats <- all_grstats[[curpool]]$volcano$`NP BL6`
  tovery <- grstats[grstats$genecat %in% c("Dispensable", "Essential"),]
  tovery$obs <- (tovery$genecat=="Dispensable")+0

  onep <- ggplot(tovery, aes(d = obs, m = fc)) + geom_roc() + xlab(paste("pos fraction,",curpool))
  listplot[[curpool]] <- onep
}
totp <- egg::ggarrange(plots=listplot)
totp

ggsave(plot = totp, "/corgi/websites/malariascreenviewer/plots_crispr/roc_perscreen.pdf", width = 10, height = 10)


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
      #theme(axis.title.x = element_text(colour = geneinfo$xlabcol[geneinfo$gene==curgene]))+
      xlab(curgene)+
      ylab("Count/Disp. gene count")
    #p1
    list_all_plot[[curgene]] <- p1
  }
  ptot <- egg::ggarrange(plots=list_all_plot, ncol=1)
  ggsave(sprintf("/corgi/websites/malariascreenviewer/plots_crispr/all_lineplot %s.pdf",curpool), width = 3, height = 2*length(list_all_plot), plot = ptot, limitsize=FALSE)
}  


  
  
  
################################################################################
################# Fig xxx. "Volcano" plots #####################################
################################################################################

toplot_all <- NULL
for(current_pool in c("cr_2024march_half1","cr_2024march_p1","cr_2024march_p12","cr_2024march_p2")){
  #current_pool <- "cr_2024march_half1"
  print(current_pool)
  thecond <- "NP BL6"
  grstats <- all_grstats[[current_pool]]
  toplot <- grstats$volcano[[thecond]]
    print(dim(toplot))
  
  toplot$y <- 1/toplot$sd
  toplot$pool <- current_pool
  
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
  
  toplot_all <- rbind(toplot_all, toplot)
#  print(head(toplot))
  
}

toplot_all[toplot_all$gene=="PBANKA_0314200",]  
  
  
ggplot(toplot_all[toplot_all$y>4,], aes(fc, y, label=gene, color=genecat, group=gene)) + 
  geom_line() + 
#  geom_point() + 
  xlab("RGR") + #xlab(paste("FC",thecond)) + 
  ylab(yname) +
  scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "black"))+ #"Dispensible","Essential","Slow growers","Other"
  #scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "turquoise3"))+ #"Dispensible","Essential","Slow growers","Other"
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  



ggplot(toplot_all, aes(fc, y, label=gene, color=gene, group=gene)) + 
  geom_line() + 
  #  geom_point() + 
  xlab("RGR") + #xlab(paste("FC",thecond)) + 
  ylab(yname) +
#  scale_color_manual(values = c("chartreuse4", "red", "dodgerblue", "turquoise3"))+ #"Dispensible","Essential","Slow growers","Other"
  theme_bw() + 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))  

  
  

ggplot(toplot_all[toplot_all$gene=="PBANKA_0314200",], aes(fc, y, color=gene, group=gene)) + 
  geom_line(color="black")
  
table(toplot_all$gene)
  

################################################################################
################# Fig 5 foo ####################################################
################################################################################
  

if(FALSE){
  
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
}


  