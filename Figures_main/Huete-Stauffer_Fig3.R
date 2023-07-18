#######################################
#####          FIGURE 3         #######
#######################################
# Created by Tamara Huete Stauffer tamara.huete@gmail.com

# - Access data and process the sequence table as in Sequence_prep.R
# - The starting point of this script is the table 'merged.rarify'
# - This script creates 3 individual plots that are joined at the last step to reproduce Figure 3

#######################################
#####      LIBRARIES           #######
#######################################
source("Requirements.R")

#######################################
#####         SEQ PREP          #######
#######################################
# Run Sequence_prep.R to get merged.rarefy
source("Sequence_prep.R")
sequence_prep()
merged.rarefy <- read.csv("Huete-Stauffer_rarefied_sequence_and_taxa.csv",row.names=1)

#Divide merged table in TAX and OTU
merged.rarefy$ASV <- row.names(merged.rarefy)
TAX <-  merged.rarefy[,c(75,68:74)]
OTU <- merged.rarefy[,c(75,1:67)]

#######################################
#####       FIG.3 A,B,C         #######
#######################################
# function to apply to each depth layer
subfig3 <- function(OTU, metadata, layer){
  row.names(metadata) <- metadata$X
  if(colnames(OTU)[1]=='ASV'){
    OTU <- OTU[2: ncol(OTU)]
  }
  if(layer == 'surface'){
    OTUs <- OTU[which(regmatches(colnames(OTU),regexpr("[_][0-9][0-9][0-9][_]",colnames(OTU)))=="_005_")] 
    metadata <- metadata[which(regmatches(row.names(metadata),regexpr("[_][0-9][0-9][0-9][_]",row.names(metadata)))=="_005_"),]
    layer_name='SURF'
  }
  
  if(layer == 'dcm'){
    OTUs <- OTU[,c(which(as.numeric(gsub("[_]","",regmatches(names(OTU),regexpr("[_][0-9][0-9][0-9][_]",names(OTU)))))>=50 & as.numeric(gsub("[_]","",regmatches(names(OTU),regexpr("[_][0-9][0-9][0-9][_]",names(OTU)))))<95),which(as.numeric(gsub("[_]","",regmatches(names(OTU),regexpr("[_][0-9][0-9][0-9][_]",names(OTU)))))==108),which(as.numeric(gsub("[_]","",regmatches(names(OTU),regexpr("[_][0-9][0-9][0-9][_]",names(OTU)))))==43))] # dcm
    metadata <- metadata[c(which(as.numeric(gsub("[_]","",regmatches(row.names(metadata),regexpr("[_][0-9][0-9][0-9][_]",row.names(metadata)))))>=50 & as.numeric(gsub("[_]","",regmatches(row.names(metadata),regexpr("[_][0-9][0-9][0-9][_]",row.names(metadata)))))<95),which(as.numeric(gsub("[_]","",regmatches(row.names(metadata),regexpr("[_][0-9][0-9][0-9][_]",row.names(metadata)))))==108),which(as.numeric(gsub("[_]","",regmatches(row.names(metadata),regexpr("[_][0-9][0-9][0-9][_]",row.names(metadata)))))==43)),] #dcm
    layer_name='DCM'
  }
  
  if(layer == 'dsl'){
    OTUs <- OTU[,which(as.numeric(gsub("[_]","",regmatches(names(OTU),regexpr("[_][0-9][0-9][0-9][_]",names(OTU)))))>=450 & as.numeric(gsub("[_]","",regmatches(names(OTU),regexpr("[_][0-9][0-9][0-9][_]",names(OTU)))))<700)] # fish
    metadata <-  metadata[which(as.numeric(gsub("[_]","",regmatches(row.names(metadata),regexpr("[_][0-9][0-9][0-9][_]",row.names(metadata)))))>=450 & as.numeric(gsub("[_]","",regmatches(row.names(metadata),regexpr("[_][0-9][0-9][0-9][_]",row.names(metadata)))))<700),] #fish
row.names(metadata)
metadata <-  metadata[-16,]
  layer_name='DSL'
  }
    
  ## Transpose table
  OTUt_rarefy <- t(sapply(OTUs, as.numeric))
  colnames(OTUt_rarefy) <- row.names(OTUs)
  OTUt_rarefy_nocero <- OTUt_rarefy[,-(which(colSums(OTUt_rarefy) == 0))] # r
  
  ### PCoA
  set.seed(200)
  OTU.pcoa <- cmdscale(vegdist(OTUt_rarefy_nocero,method="bray"),k=2, eig=TRUE)    
  Eigenvalues <- eigenvals(OTU.pcoa) 
  Variance <- Eigenvalues / sum(Eigenvalues) 
  Variance1 <- 100 * signif(Variance[1], 2) 
  Variance2 <- 100 * signif(Variance[2], 2)
  
  ### in ggplot to merge with cowplot to abundance
  data4Plot <- as.data.frame(OTU.pcoa$points)
  data4Plot$sample <- row.names(OTU.pcoa$points)
  
  metadata$sample <- row.names(metadata)
  data4Plot <- data4Plot %>% inner_join(metadata, by="sample")
  
  #### Correlation of PCoA scores to metadata
  names(data4Plot)
  corrvar <- c(1,2,16:19,22:24,27,47:54,67)
  res2 <- cor.mtest(data4Plot[,corrvar], conf.level = .95,adjust.method="Bonferroni")
  corrplot(cor(data4Plot[,corrvar]),method="number",p.mat = res2$p, sig.level = .05, type="upper")
  
  ### organize levels of factors
  data4Plot$season <- factor(data4Plot$season, levels=c("winter","spring","summer"))
  colSea <- c("deepskyblue4","chartreuse3","firebrick")
  step=0.09
  
  layer.gg <- ggplot(data4Plot, aes(x =V1, y = V2))+
    geom_point(size=0.1)+
    stat_ellipse(level = 0.95, aes(x=V1,y=V2,color=season))+
    geom_point(size = 1.5,shape=21,aes(x =V1, y = V2, fill= season),data=data4Plot,stroke=0.3)+
    geom_point(size = 1,shape=19,aes(x =V1, y = V2,color=season),data=data4Plot,stroke=0.3)+
    scale_fill_manual(values=rep("gray60",5))+
    scale_color_manual(values=c("deepskyblue3","yellowgreen","brown3"))+
    theme_bw()+
    theme(
      text=element_text(size=8),
      axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=0.5,color="black",size=6),
      axis.text.y = element_text(color="black",size=6),
      panel.background = element_rect(fill = "transparent"), # bg of the panel
      plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
      panel.grid.minor = element_blank(), # get rid of minor grid
      panel.grid.major = element_blank(), # get rid of major grid
      legend.position = c(0.95,0.25),
      panel.border = element_rect(color="black")) +
    guides(color=FALSE,fill=FALSE)
  
  if(layer == 'surface') {
   layer.gg = layer.gg + 
     annotate("text", x=-0.5,y=-0.4,label=paste(layer_name," (",Variance1+Variance2,"%)",sep=""),cex=3, hjust=0)+
    #pval pairwise adonis
    annotate("text", x=-0.5,y=0.3,label="PERMANOVA",cex=2.5, hjust=0)+
    annotate("point", x=rep(c(-0.5,-0.45),3),y=c(rep(0.25,2),rep(0.21,2),rep(0.17,2)),col=c("deepskyblue3","yellowgreen","deepskyblue3","brown3","yellowgreen","brown3"),cex=1.6) +
    annotate("text", x=rep(-0.40),y=c(0.25,0.21,0.171),label=c("p.adj = 0.003","p.adj = 0.021","p.adj = 0.003"),cex=2, hjust=0)+
    xlab(label=bquote(paste("PCoA 1 (" ,.(format(Variance1, digits = 3)),"%) (+PAb)"))) +
    ylab(label=bquote(paste("PCoA 2 (" ,.(format(Variance2, digits = 3)),"%) (-Temp)")))+
    theme(plot.margin = unit(c(0.2,0.2,0.2,0.4), "lines"))+
    theme(axis.title.y = element_text(margin = margin(t = 0, r = -1, b = 0, l = 0)))
  }
  
  if(layer == 'dcm') {
   layer.gg = layer.gg + annotate("text", x=-1.3,y=-0.4,label=paste("DCM (",Variance1+Variance2,"%)",sep=""),cex=3, hjust=0)+
  annotate("text", x=-1.3,y=0.36,label="PERMANOVA",cex=2.5, hjust=0)+
  annotate("point", x=rep(c(-1.3,-1.2),3),y=c(rep(0.31,2),rep(0.27,2),rep(0.23,2)),col=c("deepskyblue3","yellowgreen","deepskyblue3","brown3","yellowgreen","brown3"),cex=1.6) +
  annotate("text", x=rep(-1.1),y=c(0.31,0.27,0.23),label=c("p.adj = 0.012","p.adj = 0.021","p.adj = 0.003"),cex=2, hjust=0)+
  xlab(label=bquote(paste("PCoA 1 (" ,.(format(Variance1, digits = 3)),"%) (-Temp, -%HNA)"))) +
  ylab(label=bquote(paste("PCoA 2 (" ,.(format(Variance2, digits = 3)),"%) (+",NO[3]^{"2-"},", +DOC)")))+
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.4), "lines"))+
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -1, b = 0, l = 0)))
}
  
  if(layer == 'dsl') {
   layer.gg = layer.gg + annotate("text", x=-0.23,y=-0.28,label=paste("DSL (",Variance1+Variance2,"%)",sep=""),cex=3, hjust=0)+
  annotate("text", x=-0.23,y=0.3,label="PERMANOVA",cex=2.5, hjust=0)+
  annotate("point", x=rep(c(-0.23,-0.2),3),y=c(rep(0.26,2),rep(0.225,2),rep(0.19,2)),col=c("deepskyblue3","yellowgreen","deepskyblue3","brown3","yellowgreen","brown3"),cex=1.6) +
  annotate("text", x=rep(-0.18),y=c(0.26,0.225,0.19),label=c("p.adj = 0.009","p.adj = 0.027","p.adj = 0.015"),cex=2, hjust=0)+
  # legend
  annotate("point", x=rep(0.15,3),y=c(0.3,0.265,0.23),col=c("deepskyblue3","yellowgreen","brown3"),cex=1.6) +
  annotate("text", x=rep(0.17),y=c(0.3,0.265,0.23),label=c("winter","spring","summer"),cex=2, hjust=0)+
  
  xlab(label=bquote(paste("PCoA 1 (" ,.(format(Variance1, digits = 3)),"%) (+PVol, +%HNA, -",NO[3]^{"2-"},")"))) +
  ylab(label=bquote(paste("PCoA 2 (" ,.(format(Variance2, digits = 3)),"%) (+PAb)")))+
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.4), "lines"))+ 
  theme(axis.title.y = element_text(margin = margin(t = 0, r = -1, b = 0, l = 0)))
}
  
  return(layer.gg)
  
}

surf.gg <- subfig3(OTU,metadata,'surface')
dcm.gg <- subfig3(OTU,metadata,'dcm')
dsl.gg <- subfig3(OTU,metadata,'dsl')

#######################################
#####   ARRANGE &   EXPORT      #######
#######################################
gg_F2 = plot_grid(surf.gg, dcm.gg,dsl.gg, labels = c('A','B','C'), nrow = 1,ncol=3)
ggsave("Huete-Stauffer_Figure3.pdf",gg_F2,height=60,width=169,units="mm",colormode="cmyk")