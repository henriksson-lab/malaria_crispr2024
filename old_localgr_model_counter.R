
library(umap)




################################################################################
####################### Generate count matrices -- barseq ######################
################################################################################

listpools <- c("barseq_slowhires_2023dec")
for(curpool in listpools){
  
  print(curpool)
  
  # re.compile("GGCGG"+"(\w{8,16})"+"CTGAC")
  
  seqbefore <- "TAGTCGCAGTAGGCGG"
  
  allpooldir <- "/corgi/otherdataset/ellenbushell/crispr_pools"
  pooldir <- file.path(allpooldir, curpool)
  fastqdir <- file.path(pooldir, "fastq")
  bcfile <- file.path(pooldir, "used_bc.csv")
  countfile <- file.path(pooldir,"counts.RDS")
  
  frank_bc <- read.csv("/corgi/otherdataset/ellenbushell/barcode_to_gene_210920_FRANK.csv")
  frank_bc$sgrna <- str_split_fixed(frank_bc$gene,"\\|",2)[,1]
  frank_bc$seq <- str_to_upper(frank_bc$barcode)
  
  #Subset by the BCs expected here
  usedbc <- read.csv(bcfile,sep="\t")
  usedbc <- frank_bc[frank_bc$sgrna %in% usedbc$gene,]
  bclength <- str_length(usedbc$seq[1]) #  TCTTTTCCCAG

  #R1 needs RC
  usedbc$seq <- as.character(Biostrings::reverseComplement(Biostrings::DNAStringSet(usedbc$seq)))
  
  
  list_bclist <- list()
  for(onef in list.files(fastqdir)){
    print(onef)
    if(str_ends(onef,"R1_001.fastq.gz")){
      onep <- pipe(paste("zcat",file.path(fastqdir,onef),"| grep TAGTCGCAGTAGGCGG"))
      li <- readLines(onep)
      close(onep)
      
      #revcomp BC for barseq R1??
      
      bclist <- str_split_fixed(li, seqbefore,2)[,2]
      bclist <- data.frame(bc=str_sub(bclist,1,bclength))
      bclist <- sqldf::sqldf("select count(bc) as cnt, bc from bclist group by bc order by cnt desc")
      bclist$file <- onef
      
      list_bclist[[onef]] <- bclist
    }
  }
  counts <- do.call(rbind,list_bclist)
  counts <- reshape2::acast(counts, bc~file, value.var = "cnt", fill = 0)
  
  #
  counts <- counts[order(rowSums(counts), decreasing = TRUE),]
  counts <- counts[rownames(counts) %in% usedbc$seq,]
  colnames(counts) <- str_sub(colnames(counts),1,11)  #dangerous. split by _S instead?
  
  rownames(usedbc) <- usedbc$seq
  rownames(counts) <- usedbc[rownames(counts),]$sgrna
  
  saveRDS(counts, countfile)
}


################################################################################
####################### Generate count matrices -- CRISPR ######################
################################################################################

# "tags" = EB_barseq_slowpool_CRISPR

#listpools <- c("tags")
listpools <- c("2023march_screen_noD4")
listpools <- c("2023march_screen")
#listpools <- c("aug_p192","aug_p24","aug_p96")
listpools <- c("2023aug_p192","2023aug_p24","2023aug_p96","2023jan_tags", "2023march_screen", "2023march_screen_noD4")
listpools <- c("2023dec_ligpool")
for(curpool in listpools){

  print(curpool)

  seqbefore <- "CAATATTATT"
  
  allpooldir <- "/corgi/otherdataset/ellenbushell/crispr_pools"
  pooldir <- file.path(allpooldir, curpool)
  fastqdir <- file.path(pooldir, "fastq")
  bcfile <- file.path(pooldir, "bc.csv")
  countfile <- file.path(pooldir,"counts.RDS")

  usedbc <- read.csv(bcfile,sep="\t")
  bclength <- str_length(usedbc$seq[1])
  
  list_bclist <- list()
  for(onef in list.files(fastqdir)){
    print(onef)
    if(str_ends(onef,"R1_001.fastq.gz")){
      onep <- pipe(paste("zcat",file.path(fastqdir,onef),"| grep CAATATTATT"))
      li <- readLines(onep)
      close(onep)
      
      bclist <- str_split_fixed(li, seqbefore,2)[,2]
      bclist <- data.frame(bc=str_sub(bclist,1,bclength))
      bclist <- sqldf::sqldf("select count(bc) as cnt, bc from bclist group by bc order by cnt desc")
      bclist$file <- onef
      
      list_bclist[[onef]] <- bclist
    }
  }
  counts <- do.call(rbind,list_bclist)
  counts <- reshape2::acast(counts, bc~file, value.var = "cnt", fill = 0)
  
  #
  counts <- counts[order(rowSums(counts), decreasing = TRUE),]
  counts <- counts[rownames(counts) %in% usedbc$seq,]
  colnames(counts) <- str_sub(colnames(counts),1,11)

  rownames(usedbc) <- usedbc$seq
  rownames(counts) <- usedbc[rownames(counts),]$sgrna

  saveRDS(counts, countfile)
}


################################################################################
####################### Coverage stats #########################################
################################################################################

if(FALSE){
  
  
  
}



################################################################################
####################### Generate curves ########################################
################################################################################




####### Loop over all pools
timecourses <- list()
all_grstats_per_grna <- list()
all_grstats <- list()
list_samplemeta <- list()
all_coverage_stat <- list()
for(curpool in listpools){
  
  print(curpool)
  
  allpooldir <- "/corgi/otherdataset/ellenbushell/crispr_pools"
  pooldir <- file.path(allpooldir, curpool)
  countfile <- file.path(pooldir,"counts.RDS")
  samplemetafile <- file.path(pooldir,"sampleinfo.txt")
  controlmetafile <- file.path(pooldir,"list_control.csv")
  
  #### Read sample metadata
  samplemeta <- read.csv(samplemetafile, sep = "\t")[,1:2]
  colnames(samplemeta) <- c("sampleid","samplename")
  samplemeta$day <- str_sub(str_split_fixed(samplemeta$samplename, "_",5)[,4],2)
  samplemeta$is_input <- str_count(samplemeta$samplename,"input")>0
  samplemeta$day <- as.integer(samplemeta$day)
  samplemeta$mouse_ref <- str_split_fixed(samplemeta$samplename, "_",5)[,5]
  samplemeta$genotype <- str_split_fixed(samplemeta$samplename, "_",5)[,3] ##"wt"
  samplemeta$primed <- str_split_fixed(samplemeta$samplename, "_",5)[,2]
  
  if(sum(samplemeta$is_input)>0){
    samplemeta$day[samplemeta$is_input] <- NA 
    samplemeta$mouse_ref[samplemeta$is_input] <- NA 
    samplemeta$genotype[samplemeta$is_input] <- NA 
    samplemeta$primed[samplemeta$is_input] <- NA 
  }
  
    
  ### Read count table
  counts <- readRDS(countfile)

  #### Read info about the cloning
  cloningfile <- file.path(pooldir,"cloning.csv")
  if(file.exists(cloningfile)){
    allgeneconstructs <- read.csv(cloningfile,sep="\t")  #read.csv("/corgi/otherdataset/ellenbushell/crispr_geneinfo.csv",sep="\t")
    allgeneconstructs$gene <- str_split_fixed(allgeneconstructs$grna,"gRNA",2)[,1]
    allgeneconstructs$genecat[allgeneconstructs$genecat==""] <- "Other"
  } else {
    print("No cloning.csv -- constructing equivalent")
    
    allgeneconstructs <- data.frame(
      grna=rownames(counts),
      gene=rownames(counts)
    )
    allgeneconstructs$genecat <- "Other"
    allgeneconstructs$genewiz <- "NA"
    allgeneconstructs$ligationwell <- "NA"

    geneinfotable <- read.csv("/corgi/otherdataset/ellenbushell/gene_description.csv",sep="\t")
    allgeneconstructs$genecat[allgeneconstructs$gene %in% geneinfotable$gene[geneinfotable$genedesc=="Dispensable"]] <- "Dispensable"
    allgeneconstructs$genecat[allgeneconstructs$gene %in% geneinfotable$gene[geneinfotable$genedesc=="Slow"]] <- "Slow"
    
  }
  
  grna_dispensible <- allgeneconstructs$grna[allgeneconstructs$genecat=="Dispensable"]
  genes_dispensible <- allgeneconstructs$gene[allgeneconstructs$genecat=="Dispensable"]
  
  
  
  ### Extract "input" control sample
  input_sampleid <- samplemeta$sampleid[samplemeta$is_input]
  coverage_stat <- data.frame(
    grna=rownames(counts),
    cnt=counts[,input_sampleid])
  coverage_stat <- merge(allgeneconstructs,coverage_stat)
  all_coverage_stat[[curpool]] <- coverage_stat

  ### Make pseudocounts
  counts <- counts + 1

  #Gather total count
  rownames(samplemeta) <- samplemeta$sampleid
  samplemeta <- samplemeta[colnames(counts),]
  samplemeta$total_count <- colSums(counts)
  
  
  #Filter bad samples
  count_stats <- colSums(counts)
  bad_libs <- count_stats<2000  #was 10000
  if(sum(bad_libs)>0){
    print("bad libs! here are counts before")
    print(count_stats)
    print("Removing:")
    print(colnames(counts)[bad_libs])
    counts <- counts[,!bad_libs]
    #print("bad libs! here are counts after")
    #count_stats <- colSums(counts)
    #print(count_stats)
  }
  
  #Normalize each library by depth
  for(i in 1:ncol(counts)){
    counts[,i] <- counts[,i]/sum(counts[,i]) 
  }
  
  #Align samplemeta with counts. Compute UMAP
  rownames(samplemeta) <- samplemeta$sampleid
  samplemeta <- samplemeta[colnames(counts),]
  umap.settings <- umap.defaults
  umap.settings$n_neighbors <- min(umap.settings$n_neighbors, ncol(counts))
  cnt.umap <- umap(t(counts), config=umap.settings)
  samplemeta$umap1 <- cnt.umap$layout[,1]
  samplemeta$umap2 <- cnt.umap$layout[,2]
  
  #rowMeans(counts[rownames(counts) %in% grna_dispensible,])
  #rowMeans(counts[!(rownames(counts) %in% grna_dispensible),])
  #counts <- counts[rowMeans(counts)>1e-5,]
  
  ####### Merge metadata with counts.
  ####### remove input samples from TC
  longcnt <- reshape2::melt(counts)
  colnames(longcnt) <- c("grna","sampleid","cnt") 
  longcnt <- merge(longcnt, samplemeta[!is.na(samplemeta$day),])
  longcnt$gene <- str_split_fixed(longcnt$grna,"gRNA",2)[,1]

  
  
  ######### Coverage stats
  if(FALSE){
    
    longcnt
    head(longcnt)
    allgeneconstructs
    qc <- merge(longcnt, allgeneconstructs)
    
    #onesample <- "pHIT_ligpool_d1_m2C"
    #part <- qc[qc$samplename==onesample,]
    
  #  ggplot(qc[qc$samplename==onesample & qc$cnt>0.0001,], aes(paste(grna),cnt)) + geom_point()
    #ggplot(qc[qc$samplename==onesample & qc$cnt>0.0001,], aes(grna,cnt)) + geom_point()

    ggplot(qc[qc$cnt>0.0001,], aes(paste(mouse_ref,grna),cnt,color=samplename)) + geom_point()
    ggplot(qc[qc$cnt>0.0001,], aes(paste(mouse_ref,grna),cnt,color=genecat)) + geom_point()
    
    qc <- qc[qc$cnt>0.005,c("samplename","grna","gene","cnt","genewiz","ligationwell","genecat","mouse_ref")]
    write.csv(qc,"~/ellen_2023dec_ligpool.csv")
    
    #contamination well
    qc[qc$mouse_ref=="m1A",]

    table(qc[,c("ligationwell","mouse_ref")])
    
    
    
    table(qc$mouse_ref)
    #sqldf::sqldf()
    
    ggplot(qc[qc$cnt>0.0001,], aes(paste(grna),cnt,color=genecat)) + facet_grid(. ~ mouse_ref) + geom_point()
    
    
  }
  
  
  #Calculate GR for each grna
  relcnt <- merge(longcnt, data.frame(
      grna=longcnt$grna,
      mouse_ref=longcnt$mouse_ref,
      prevcnt=longcnt$cnt,
      day=longcnt$day+1
    ))
  relcnt$gr <- log2(relcnt$cnt/relcnt$prevcnt)
  
  
  do_control_compensation <- TRUE

  ##### Figure 
  if(file.exists(controlmetafile)){
    print("using list of controls-file")
    list_controls <- read.csv(controlmetafile)[,1]
  } else {
    print("using dispensable as controls")
    list_controls <- unique(relcnt$gene[relcnt$gene %in% genes_dispensible])
  }
  print(list_controls)
  
    
  ##### Normalize genes GR by dispensable genes to get RGR --- #################### not always present; then stick with GR
  controlcnt <- relcnt[relcnt$gene %in% list_controls,]
  if(nrow(controlcnt)==0){
    print("No controls - all genes will be used as the control")
    controlcnt <- relcnt
    do_control_compensation <- FALSE
  } else {
    print("using controls")
    print(unique(relcnt$gene[relcnt$gene %in% genes_dispensible]))
  }
  controlcnt <- as.data.frame(controlcnt %>%
    group_by(day, mouse_ref) %>%
    summarise_at(vars(gr), list(control_var = var, control_gr=mean)))
  relcnt_norm <- merge(relcnt, controlcnt)  
  
  if(do_control_compensation){
    relcnt_norm$rgr <- relcnt_norm$gr - relcnt_norm$control_gr
  } else {
    relcnt_norm$rgr <- relcnt_norm$gr 
  }
  if(any(relcnt_norm$grna %in% allgeneconstructs$grna)){
    relcnt_norm <- merge(allgeneconstructs, relcnt_norm)  ## Only relevant for pools where we got info on source wells
  } else {
    print("No plate info")
  }
  relcnt_norm <- merge(longcnt[,c("grna","sampleid","cnt")], relcnt_norm)
  timecourses[[curpool]] <- relcnt_norm

  
  
  
  ##### PER GRNA: Calculate FC vs control
  localvar <- as.data.frame(relcnt_norm %>%
                              group_by(mouse_ref, grna) %>%
                              summarise_at(vars(rgr), list(sample_var = var, sample_mean=mean)))
  localvar <- merge(localvar,controlcnt)
  #localvar$sample_var[is.na(localvar$sample_var)] <- 0 #Problematic that this happens; put in average instead? or vst
  if(do_control_compensation){
    localvar$totalvar <- localvar$sample_var + localvar$control_var
  } else {
    localvar$totalvar <- localvar$sample_var 
  }
  
  grna_var <- sqldf::sqldf("select sum(sample_var) as totalvar, count(sample_mean) as cnt, avg(sample_mean) as fc, grna from localvar group by grna")
  grna_var$sd <- sqrt(grna_var$totalvar/grna_var$cnt)
  grna_var$p <- pnorm(grna_var$fc,0,sd=grna_var$sd)
  grna_var$p[grna_var$p>0.5] <- 1 - grna_var$p[grna_var$p>0.5]
  grna_var$logp <- -log10(grna_var$p)
  
  grna_var <- merge(grna_var,allgeneconstructs) #maybe an outer join later? TODO
  all_grstats_per_grna[[curpool]] <- grna_var

  
  #localvar[localvar$gene=="PBANKA_0515000",]
  
  
  ##### PER GENE: Calculate FC vs control (average over days!)
  gene_var <- sqldf::sqldf("select sum(totalvar) as totalvar, count(grna) as cnt, avg(fc) as fc, grna from grna_var group by grna")
  gene_var$sd <- sqrt(gene_var$totalvar/gene_var$cnt)
  
  
#  gene_var <- sqldf::sqldf("select sum(sample_var) as totalvar, count(sample_mean) as cnt, avg(sample_mean) as fc, gene from localvar group by gene")
 # gene_var$sd <- sqrt(gene_var$totalvar/gene_var$cnt)
  gene_var$p <- pnorm(gene_var$fc,0,sd=gene_var$sd)
  gene_var$p[gene_var$p>0.5] <- 1 - gene_var$p[gene_var$p>0.5]
  gene_var$logp <- -log10(gene_var$p)
  gene_var <- merge(gene_var,allgeneconstructs) #maybe an outer join later? TODO
  list_all_grstats <- list()
  list_all_grstats[["vs control"]] <- gene_var
  # 
  # 
  # ##### PER GENE: Calculate FC vs control (average over days!)
  # localvar <- as.data.frame(relcnt_norm %>%
  #                               group_by(day, mouse_ref, gene) %>%
  #                               summarise_at(vars(rgr), list(sample_var = var, sample_mean=mean)))
  # localvar <- merge(localvar,controlcnt)
  # localvar$sample_var[is.na(localvar$sample_var)] <- 0 #Problematic that this happens; put in average instead? or vst
  # localvar$totalvar <- localvar$sample_var + localvar$control_var
  # 
  # genevar <- sqldf::sqldf("select sum(sample_var) as totalvar, count(sample_mean) as cnt, avg(sample_mean) as fc, gene from localvar group by gene")
  # genevar$sd <- sqrt(genevar$totalvar/genevar$cnt)
  # genevar$p <- pnorm(genevar$fc,0,sd=genevar$sd)
  # genevar$p[genevar$p>0.5] <- 1 - genevar$p[genevar$p>0.5]
  # genevar$logp <- -log10(genevar$p)
  # genevar <- merge(genevar,allgeneconstructs) #maybe an outer join later? TODO
  # list_all_grstats <- list()
  # list_all_grstats[["vs control"]] <- genevar
  # 
  

  all_grstats[[curpool]] <- list(
    volcano=list_all_grstats,
    scatterplot=list()
  )
  
  list_samplemeta[[curpool]] <- samplemeta
}

saveRDS(all_grstats, file="/corgi/websites/malariascreenviewer/grstats.rds")
saveRDS(timecourses, file="/corgi/websites/malariascreenviewer/timecourses.rds")
saveRDS(list_samplemeta, file="/corgi/websites/malariascreenviewer/samplemeta.rds")
saveRDS(all_coverage_stat, file="/corgi/websites/malariascreenviewer/coverage_stat.rds")

# coverage_stat <- coverage_stat[order(coverage_stat$cnt, decreasing = TRUE),]
# coverage_stat$grna <- factor(coverage_stat$grna, levels=coverage_stat$grna)
# ggplot(coverage_stat, aes(grna,cnt,color=genecat)) + geom_point() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))




################################################################################
################### 2023jan_tags -- supplemental figure ########################
################################################################################


grstats <- timecourses[["2023jan_tags"]]

if(FALSE){
  #Average over mice
  grstats <- sqldf::sqldf(
    "select day, avg(rgr) as rgr, avg(cnt) as cnt, gene, grna from grstats group by gene, grna, day")
  grstats$mouse_ref <- "m*"
}

grstats$group <- paste(grstats$grna, grstats$mouse_ref)
grstats$grna_index <- str_split_fixed(grstats$grna, "g",2)[,2]

ggplot(grstats,aes(x=day,y=cnt, group=group, color=gene, text=group)) + 
  geom_line()+
  geom_point(aes(shape=grna_index))+
  xlab("Day")+
  ylab("Rel. abundance")+
  ggtitle("")
ggsave("/home/mahogny/ellen/2023jan_tags.pdf")




################################################################################
############### 2023march_screen main figure -- volcano ########################
################################################################################

forpool <- "2023march_screen_noD4"
forpool <- "2023march_screen"

### why are not some genes shown in the volcano?
### control genes in particular

#toplot <- all_grstats$`2023march_screen_noD4`$volcano$`vs control`
toplot <- all_grstats[[forpool]]$volcano$`vs control`
##TODO
#all_grstats$aug_p192$volcano$`vs control`
#all_grstats_per_grna$`2023march_screen`$grna


#current_pool <- input$grstats_pool
#print(current_pool)
#grstats <- all_grstats[[current_pool]]
#thecond <- input$grstats_volcano

#if(thecond %in% names(grstats$volcano)){
#  toplot <- grstats$volcano[[thecond]]
toplot <- unique(toplot[,c("fc","logp","gene","genecat")])
toplot$genecat <- factor(toplot$genecat, levels = c("Essential", "Dispensable", "Slow"))

ggplot(unique(toplot[,c("fc","logp","gene","genecat")]), aes(fc, logp, label=gene, color=genecat)) + 
  geom_point(color="gray") + 
  ggrepel::geom_text_repel(max.overlaps = 100)+
#  geom_text() +
#  xlim(c(-3,3)) +
#  ylim(c(0,10)) +
  xlab(paste("FC")) + 
  ylab(paste("-log10 pval")) + 
  scale_color_manual(values=c("#bf212f", "#27b376","#264b96"))
  #scale_color_manual(values=c("darkred", "darkgreen","blue"))

ggsave("/home/mahogny/ellen/2023march_screen__volcano.pdf", width = 12, height = 6)

#Disp: vara grön,
#Ess: vara röd
#slow: blue
#


################################################################################
############### 2023march_screen main figure -- rel abundance ##################
################################################################################

forpool <- "2023march_screen_noD4"
forpool <- "2023march_screen"
gene2genecat <- unique(all_grstats[[forpool]]$volcano$`vs control`[,c("gene","genecat")])
grstats <- timecourses[[forpool]]

#Average over mice
if(FALSE){
  grstats <- sqldf::sqldf(
    "select day, avg(rgr) as rgr, avg(cnt) as cnt, gene, grna from grstats group by gene, grna, day")
  grstats$mouse_ref <- "m*"
}

grstats$group <- paste(grstats$grna, grstats$mouse_ref)
grstats$grna_index <- str_split_fixed(grstats$grna, "g",2)[,2]


#Loop over all conditions
gene2genecat <- gene2genecat[order(gene2genecat$genecat),]
gene2genecat$i <- 1:4 #note, reusing!
gene2genecat <- gene2genecat[order(gene2genecat$i,gene2genecat$gene),]
for(i in 1:nrow(gene2genecat)){
    onegene <- gene2genecat$gene[i]
    
    substats <- grstats[grstats$gene==onegene,]
    substats <- substats[substats$day!=4,]
    
    oneplot <- ggplot(substats, aes(x=day,y=cnt, group=group, color=grna_index, text=group)) + 
      geom_line()+
      geom_point(aes(shape=grna_index))+
      xlab(onegene)+
      ylab("Rel. abundance")+
      ggtitle("")+ 
      scale_x_continuous(breaks=min(substats$day):max(substats$day))+
      theme_bw()+
      theme(legend.position = "none")

    allplots[[i]] <- oneplot
}
ptot <- do.call(gridExtra::grid.arrange, allplots)
plot(ptot)
ggsave(plot = ptot, "/home/mahogny/ellen/2023march_screen__relabund.pdf", width = 10, height = 10)

# essential, dispensible, slow

#note: 2 grna missing here, and it is as expected


################################################################################
############### 2023march_screen main figure -- fold change ####################
################################################################################

forpool <- "2023march_screen_noD4"
forpool <- "2023march_screen"

gene2genecat <- unique(all_grstats[[forpool]]$volcano$`vs control`[,c("gene","genecat")])

grstats <- all_grstats_per_grna[[forpool]]  #$`2023march_screen`
grstats$i <- str_split_fixed(grstats$grna,"g",2)[,2]

ggplot(grstats, aes(x=gene, y=fc, group=i)) + 
  geom_pointrange(data=grstats, aes(ymin=fc-sd, ymax=fc+sd,color=gene), 
                  #group=gene, color=i,
                position=position_dodge(.9)) + 
  theme_bw()+
  theme(legend.position = "none")+
  coord_flip() +
  #ylim(c(-7,7)) + 
  xlab("")+
  ylab("RGR")
ggsave("/home/mahogny/ellen/2023march_screen__fc.pdf")









#Volcano:
#  Ta bort dålig gen från control, räkna om
#Färg, volcano
#Disp vara grön,
#Ess vara röd
#
#Symmetri X axis
#Gör söt och fin för MS
#Clicka på gene få rel abund över tid för endast den genen (gRNA och mus separat)


