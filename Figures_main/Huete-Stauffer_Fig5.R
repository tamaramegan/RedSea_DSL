#######################################
#####          FIGURE 5         #######
#######################################
# Created by Tamara Huete Stauffer tamara.huete@gmail.com

# - Access data and process the sequence table as in Sequence_prep.R
# - The starting point of this script is the table 'merged.rarify'

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

OTU <- merged.rarefy[,c(1:67)]

#######################################
#####           FIG 5          #######
#######################################
### select ms, dsl, md
OTUsm <- OTU[,which(as.numeric(gsub("[_]","",regmatches(names(OTU),regexpr("[_][0-9][0-9][0-9][_]",names(OTU)))))>200)] ### MESOPELAGIC including FISH
OTUs <- OTUsm[,which(regexpr("CCF_2018",names(OTUsm))>0)] ### 2018 samples only
names(OTUs)
dim(OTUs) #8889 18

## Transpose table
OTUt_rarefy <- t(sapply(OTUs, as.numeric))

#Subsampling 
rowSums(OTUt_rarefy)
OTUt_rarefy_nocero <- OTUt_rarefy[,-(which(colSums(OTUt_rarefy) == 0))] # r
row.names(OTUt_rarefy_nocero) 


## Only spring meso
metadata <- read.csv("https://zenodo.org/record/5816123/files/Huete-Stauffer_Metadata.csv?download=1")
row.names(metadata) <- metadata$X
metadata <- metadata[row.names(metadata) %in% names(OTUs),]
row.names(metadata)

### organize levels of factors
metadata$layer1 <- droplevels(metadata$layer1)
metadata$layer1 <- factor(metadata$layer1, levels=c("MESOs","FISH","MESOd"))
with(metadata, levels(layer1))
colLy <- c("deepskyblue4","darkorange","black")


### PCoA
set.seed(200)
OTU.pcoa <- cmdscale(vegdist(OTUt_rarefy_nocero,method="bray"),k=2, eig=TRUE)    ### this worked ### CMD is principal coordinate analysis (PCoA of euclidean distances is PCA)
Eigenvalues <- eigenvals(OTU.pcoa) 
Variance <- Eigenvalues / sum(Eigenvalues) 
Variance1 <- 100 * signif(Variance[1], 2) ## 47%
Variance2 <- 100 * signif(Variance[2], 2) ## 14%
#MESO= V1=47  V2=14%  sum = 61%


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

##V1(X): - 0.91 Sal, + 0.87 PO4, + 0.7 NO3, 0.97 BA, -0.67 Bsize, -0.47%HNA, + 0.91 C1, + 0.83 C3
##V2(X): - 0.57 DOC

### organize levels of factors
data4Plot$season <- factor(data4Plot$layer1, levels=c("MESOs","FISH","MESOd"))
colLy <- c("deepskyblue4","darkorange","black")
step=0.09
y1=0.38

meso.gg <- ggplot(data4Plot, aes(x =V1, y = V2))+
  geom_point(size=0.1)+
  stat_ellipse(level = 0.95, aes(x=V1,y=V2,color=layer1))+
  #stat_ellipse(aes(x=V1,y=V2,fill=season,color=NULL),geom="polygon", level=0.95, alpha=0.2,type="t") +
  geom_point(size = 1.5,shape=21,aes(x =V1, y = V2, fill= season),data=data4Plot,stroke=0.3)+
  geom_point(size = 1,shape=19,aes(x =V1, y = V2,color=season),data=data4Plot,stroke=0.3)+
  scale_fill_manual(values=rep("gray60",5))+
  scale_color_manual(values=c(colors()[637],colors()[631],colors()[132]))+
  #scale_fill_manual(values=c("deepskyblue3","yellowgreen","brown3"))+
  theme_bw()+
  #scale_y_continuous(expand=c(0,0),limits=c(-0.4, 0.4))+
  #scale_x_continuous(expand=c(0,0),limits=c(-1, 0.5))+
  #theme(text=element_text(family="Arial",size=8)) +
  theme(
    text=element_text(size=8),
    axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=0.5,color="black",size=6),
    axis.text.y = element_text(color="black",size=6),
    panel.background = element_rect(fill = "transparent"), # bg of the panel
    plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
    panel.grid.minor = element_blank(), # get rid of minor grid
    panel.grid.major = element_blank(), # get rid of major grid
    #legend.background = "none", # get rid of legend bg
    #legend.box.background = element_rect(fill = "transparent") , #get rid of legend panel bg
    legend.position = c(0.5,0.2),
    panel.border = element_rect(color="black")) +
  guides(fill=FALSE, color=FALSE)+
  #legend
  annotate("text", x=0.18,y=0.37,label="MESOPELAGIC",cex=3, hjust=0)+
  annotate("point", x=rep(0.2,3),y=c(0.32,0.27,0.22),col=c(colors()[637],colors()[631],colors()[132]),cex=1.6) +
  annotate("text", x=rep(0.23,3),y=c(0.32,0.27,0.22),label=c("MS","DSL","MD"),cex=2, hjust=0)+

  #pval pairwise adonis
  annotate("text", x=0.2,y=-0.25,label="PERMANOVA",cex=2.5, hjust=0)+
  annotate("point", x=rep(c(0.2,0.23),3),y=c(rep(-0.3,2),rep(-0.35,2),rep(-0.4,2)),col=c(colors()[631],colors()[637],colors()[631],colors()[132],colors()[637],colors()[132]),cex=1.6) +
  annotate("text", x=rep(0.25),y=c(-0.3,-0.35,-0.4),label=c("p.adj = 0.036","p.adj = 0.009","p.adj = 0.012"),cex=2, hjust=0)+
  xlab(label=bquote(paste("PCoA 1 (" ,.(format(Variance1, digits = 3)),"%)"))) +
  ylab(label=bquote(paste("PCoA 2 (" ,.(format(Variance2, digits = 3)),"%)")))+
  theme(plot.margin = unit(c(0.2,0.2,0.2,0.2), "lines"))
meso.gg

meso.gg + ggsave("Huete-Stauffer_Figure5.pdf",width=86,height=86,units="mm",colormodel="cmyk")


### season significant  r = 0.52 PERMANOVA (permutational ANOVA) p=0.001
set.seed(300)
adonis_16S <- adonis(OTUt_rarefy_nocero ~ layer1, data = metadata, na.rm=TRUE)
adonis_16S # p < 0.01

library(pairwiseAdonis)
# Compare habitats 2 by 2 
set.seed(500) ### 500 all layeres significant from one another
pairwise.adonis(OTUt_rarefy_nocero,metadata$layer1)
#           pairs Df  SumsOfSqs  F.Model        R2 p.value p.adjusted sig
#1  MESOs vs FISH  1 0.09017404 3.158736 0.2400486   0.012      0.036   .
#2 MESOs vs MESOd  1 0.34660603 9.757165 0.4938545   0.004      0.012   .
#3  FISH vs MESOd  1 0.17375875 6.081590 0.3781710   0.003      0.009   *
set.seed(300) ### 300 MS-DSL =
pairwise.adonis(OTUt_rarefy_nocero,metadata$layer1)
#pairs Df  SumsOfSqs  F.Model        R2 p.value p.adjusted sig
#1  MESOs vs FISH  1 0.09017404 3.158736 0.2400486   0.023      0.069    
#2 MESOs vs MESOd  1 0.34660603 9.757165 0.4938545   0.004      0.012   .
#3  FISH vs MESOd  1 0.17375875 6.081590 0.3781710   0.005      0.015 