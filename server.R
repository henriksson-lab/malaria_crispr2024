library(plotly)
library(Cairo)
options(shiny.usecairo=T)

if(FALSE){
  #To run this app
  library(shiny)
  runApp(".")
}



server <- function(input, output, session) {

  observeEvent(input$grstats_pool,{
    current_pool <- input$grstats_pool
    grstats <- all_grstats[[current_pool]]
    
    updateSelectizeInput(session, 'grstats_volcano', choices = names(grstats$volcano), server = TRUE)
    updateSelectizeInput(session, 'grstats_scatter', choices = names(grstats$scatterplot), server = TRUE)
    
    grstats <- all_timecourses[[current_pool]]
    all_gr_type <- names(grstats)
    updateSelectizeInput(session, 'grstats_gene', choices = c("",unique(grstats[[all_gr_type[1]]]$gene)), server = TRUE)
    updateSelectizeInput(session, 'grstats_units', choices = all_gr_type, server = TRUE)

  })
  


  ################################################################################
  ########### Sample metadata ####################################################
  ################################################################################

  output$plotSamplemetaUmap <- renderPlot(height=700, {

    current_pool <- input$samplemeta_pool
    print(current_pool)
    
    samplemeta <- all_samplemeta[[current_pool]]
    coverage_stat <- all_coverage_stat[[current_pool]]
    
    print(samplemeta)
    
    samplemeta$day <- sprintf("d%s", samplemeta$day)
    p1 <- ggplot(samplemeta, aes(umap1,umap2,color=mouse_ref))+geom_point()
    p2 <- ggplot(samplemeta, aes(umap1,umap2,color=day))+geom_point()
    p3 <- ggplot(samplemeta, aes(umap1,umap2,color=is_input))+geom_point()
    p4 <- ggplot(samplemeta, aes(umap1,umap2,color=genotype))+geom_point()
    p5 <- ggplot(samplemeta, aes(umap1,umap2,color=primed))+geom_point()
    p6 <- ggplot(samplemeta, aes(umap1,umap2,color=total_count))+geom_point()
    ptot <- p1/p2|p3/p4|p5/p6 #|p7
    
    
    coverage_stat <- coverage_stat[order(coverage_stat$cnt, decreasing = TRUE),]
    coverage_stat$grna <- factor(coverage_stat$grna, levels=coverage_stat$grna)
    covstatplot <- ggplot(coverage_stat, aes(grna,cnt,color=genecat)) + geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    wellplot <- ggplot(coverage_stat, aes(genewiz, ligationwell, color=log10(cnt))) + geom_point()
    
    
    ptot/(covstatplot|wellplot)
    
    

  })

  ################################################################################
  ########### For a single screen - volcano ######################################
  ################################################################################
  
  #TODO click on a gene to show time course beneath
  
  get_current_volcano <- function(){
    current_pool <- input$grstats_pool
    print(current_pool)
    grstats <- all_grstats[[current_pool]]
    thecond <- input$grstats_volcano
    
    if(thecond %in% names(grstats$volcano)){
      grstats$volcano[[thecond]]
    } else {
      data.frame()
    }
  }
  
  output$plot_grstats_volcano <- renderPlotly({
    
    #could use get_current_volcano  above
    
    current_pool <- input$grstats_pool
    print(current_pool)
    grstats <- all_grstats[[current_pool]]
    thecond <- input$grstats_volcano
    
    if(thecond %in% names(grstats$volcano)){
      toplot <- grstats$volcano[[thecond]]
      theplot <- ggplot(toplot, aes(fc, logp, label=gene, color=genecat)) + 
        geom_point(color="gray") + 
        geom_text() +
        xlab(paste("FC",thecond)) + 
        ylab(paste("-log10 pval",thecond))
    } else {
      print("missing condition")
      theplot <- ggplot() + theme_void()
    }
    theplot  %>% ggplotly(source="plot_grstats_volcano") %>% event_register("plotly_click")
  })
  
  
  observeEvent(
    eventExpr = event_data("plotly_click", source = "plot_grstats_volcano"),
    handlerExpr = {
      print("plotly_click")
      event_data <- event_data("plotly_click", source = "plot_grstats_volcano")
      print(event_data)
      clicked_gene <- get_current_volcano()$gene[event_data$pointNumber+1]  #plotly seems to do 0-indexing
      print(get_current_volcano())
      print(clicked_gene)
      updateSelectInput(session, "grstats_gene", selected = clicked_gene)
    }
  )  
  
  
################################################################################
########### Comparison of two screens - scatter or volcano #####################
################################################################################

  get_current_scatter <- function(){
    current_pool <- input$grstats_pool
    grstats <- all_grstats[[current_pool]]
    thecond <- input$grstats_scatter
    
    if(thecond %in% names(grstats$scatterplot)){
      grstats$scatterplot[[thecond]]
    } else {
      data.frame()
    }
  }
  
  
  
  output$plot_grstats_scatterplot <- renderPlotly({
    
    current_pool <- input$grstats_pool
    grstats <- all_grstats[[current_pool]]
    thecond <- input$grstats_scatter
    
    represent_as <- input$grstats_scatter_type  ##### possibly store xlab and ylab names
    
    
    if(thecond %in% names(grstats$scatterplot)){
      toplot <- grstats$scatterplot[[thecond]]
      thecond2 <- str_split_fixed(thecond," / ",2)  #hopefully works
      cond1 <- thecond2[1]
      cond2 <- thecond2[2]
      
      
      if(represent_as=="FC scatter plot"){
        
        fc_range <- range(c(toplot$fc1, toplot$fc2))
        theplot <- ggplot(toplot, aes(fc1,fc2, label=gene, color=genedesc)) + geom_point(color="gray") + geom_text()+#size=1) +
          xlab(paste("FC",cond1)) + ylab(paste("FC",cond2)) +
          xlim(fc_range[1], fc_range[2]) + ylim(fc_range[1], fc_range[2])
        
      } else {
        
        theplot <- ggplot(toplot, aes(diff_fc, diff_log_p, label=gene, color=genedesc)) + 
          geom_point(color="gray") + 
          geom_text() +
          xlab(paste("FC",thecond)) + 
          ylab(paste("-log10 pval",thecond))
        
      }
      
    } else {
      print("missing comparison cond")
      theplot <- ggplot() + theme_void()
    }
    theplot %>% ggplotly(source="plot_grstats_scatterplot") %>% event_register("plotly_click")
  })
  
  
  
  observeEvent(
    eventExpr = event_data("plotly_click", source = "plot_grstats_scatterplot"),
    handlerExpr = {
      event_data <- event_data("plotly_click", source = "plot_grstats_scatterplot")
      #print(event_data)
      clicked_gene <- get_current_scatter()$gene[event_data$pointNumber+1]  #plotly seems to do 0-indexing
      updateSelectInput(session, "grstats_gene", selected = clicked_gene)
      }
  )  
  
  

  ################################################################################
  ########### GRstats - timecourse ###############################################
  ################################################################################
  
  output$plot_grstats_tcplot <- renderPlotly({
    
    current_pool <- input$grstats_pool
    grstats <- all_timecourses[[current_pool]]
    
    ## Pick the right unit to show
    grstats <- grstats[[input$grstats_units]]  
    print(head(grstats))
    
    ########### Average together based on user input
    
    if(input$grstats_avg_grna){
      grstats <- sqldf::sqldf(
        "select day, avg(y) as y, gene, primed, genotype, mouse_ref from grstats group by mouse_ref, gene, day, primed, genotype")
      grstats$grna <- paste(grstats$gene,"*",sep="")
    }

    if(input$grstats_avg_mouse){
      grstats <- sqldf::sqldf(
        "select day, avg(y) as y, gene, grna, primed, genotype from grstats group by gene, grna, day, primed, genotype")
      grstats$mouse_ref <- "m*"
    }

    if(input$grstats_avg_genotype){
      grstats <- sqldf::sqldf(
        "select day, avg(y) as y, gene, grna, primed from grstats group by gene, grna, day, primed")
      grstats$primed <- "g*"
    }
    
    if(input$grstats_avg_treatment){
      grstats <- sqldf::sqldf(
        "select day, avg(y) as y, gene, grna, genotype from grstats group by gene, grna, day, genotype")
      grstats$primed <- "t*"
    }

    grstats$group <- paste(grstats$grna, grstats$mouse_ref, grstats$primed, grstats$genotype)
    
    ######## Only show one gene, optionally
    current_gene <- input$grstats_gene
    if(current_gene!=""){
      grstats <- grstats[grstats$gene==current_gene,,drop=FALSE]
    }

    
    ######## Decide coloring strategy
    grstats$colorby <- grstats$mouse_ref
    if(input$grstats_colorby=="Gene"){
      grstats$colorby <- grstats$gene
    }
    if(input$grstats_colorby=="Genotype"){
      grstats$colorby <- grstats$genotype
    }
    if(input$grstats_colorby=="Treatment"){
      grstats$colorby <- grstats$primed
    }
    if(input$grstats_colorby=="Genotype+Treatment"){
      grstats$colorby <- paste(grstats$genotype,grstats$primed)
    }
    if(input$grstats_colorby=="Genetic construct"){
      grstats$colorby <- paste(grstats$grna)
    }
    if(input$grstats_colorby=="Genotype+Treatment+Genetic construct"){
      grstats$colorby <- paste(grstats$genotype,grstats$primed, grstats$grna)
    }

    ######## The actual plotting
    ggplotly(ggplot(grstats,aes(x=day,y=y, group=group, color=colorby, text=group)) + 
               geom_line()+
               xlab("Day")+
               ggtitle(""), tooltip = c("x", "y", "color", "text", "group"))

  })
  
    
}

