#######################################
#####          FIGURE 2         #######
#######################################
# Created by Tamara Huete Stauffer tamara.huete@gmail.com

# - Access data and process the sequence table as in Sequence_prep.R
# - The starting point of this script is the table 'merged.rarify'
# - This script creates 4 individual plots that are joined at the last step to reproduce Figure 2

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

# Divide merged table in TAX and OTU
names(merged.rarefy)
merged.rarefy$ASV <- row.names(merged.rarefy)
TAX <-  merged.rarefy[,c(75,68:74)]
OTU <- merged.rarefy[,c(75,1:67)]

#######################################
#####       METADATA PREP       #######
#######################################

metadata <- read.csv("https://zenodo.org/record/5816123/files/Huete-Stauffer_Metadata.csv?download=1")
metadata <- metadata[- c(which(metadata$X=="CCF_20180321_10_510_F"),which(metadata$layer==100)),]
metadata$newlayer <- metadata$layer
metadata$newlayer <- plyr::mapvalues(metadata$newlayer, from = c("250", "300","DEEP","surface","FISH"), to = c("MESOs", "MESOs","MESOd","5m","DSL"))
metadata$newlayer <-  factor(metadata$newlayer)
metadata$newlayer <-  factor(metadata$newlayer,levels=c("5m","DCM","MESOs","DSL","MESOd"))
metadata$season <- factor(metadata$season, levels=c("winter","spring","summer"))
metadata$newlabel <-  paste(metadata$NewStation,metadata$newlayer,sep=" ") 

metadata <- metadata[order(metadata$newlayer,metadata$NewStation,metadata$season,metadata$newlayer),]
metadata$goodorder <- c(1:nrow(metadata))

metadata1 <- metadata %>% 
  dplyr::select(X,goodorder,newlabel,newlayer,NewStation,season) %>% 
  mutate(sample = X,X=NULL)

#######################################
#####           FIG.2 A         #######
#######################################

## Divide Proteobacteria
TAX$PhylaClass <- ifelse(TAX$Phylum == "Proteobacteria",paste("Proteobacteria (",substr(TAX$Class,1,5),")",sep=""),as.character(TAX$Phylum))
head(TAX)

### Top 10 Phyla + divide Proteobacteria in Alpha, Beta, Gamma, Others
top10 <- OTU %>%  #### OTU rela ab 0-100 ## Mean made in table
  gather(sample,value=count,2:68)  %>%
  group_by(sample,ASV) %>%
  summarise(count = sum(count)) %>%
  mutate(relab = count/sum(count)*100) %>% 
  #group_by(sample) %>% summarise(srelab = sum(relab)) %>% ungroup() %>%## Check that all samples sum=100
  inner_join(metadata1,by='sample') %>%
  inner_join(TAX,by="ASV") %>%
  group_by(PhylaClass,sample) %>% # instead of Phylum
  dplyr::summarise(relab=sum(relab)) %>%
  dplyr::summarise(mrelab=mean(relab)) %>% ### rarefied 41 rows (top 10 + 31) /// counts not rar 42 rows
  rename(Phylum = PhylaClass) %>%
  top_n(12,mrelab) %>% mutate_if(is.factor, as.character)


### How many Other Phylum
other <- OTU %>%  #### OTU rela ab 0-100 ## Mean made in table
  #gather(sample,value=relab,2:77) %>% inner_join(metadata1,by='sample') %>% # if with the relab table
  gather(sample,value=count,2:68)  %>%# if with the count table
  group_by(sample,ASV) %>%
  summarise(count = sum(count)) %>%
  mutate(relab = count/sum(count)*100) %>% 
  #group_by(sample) %>% summarise(srelab = sum(relab)) %>% ungroup() %>%## Check that all samples sum=100
  inner_join(metadata1,by='sample') %>%
  filter(newlayer != 100) %>%
  inner_join(TAX,by="ASV") %>%
  group_by(Phylum,sample) %>% 
  dplyr::summarise(relab=sum(relab)) %>%
  dplyr::summarise(mrelab=mean(relab)) %>% 
  filter(!Phylum %in% top10$Phylum) 

  nrow(other) ### 30 Phylum other than the top 10

# Other phylum in each layer
OTU %>%  #### OTU rela ab 0-100 ## Mean made in table
  #gather(sample,value=relab,2:77) %>% inner_join(metadata1,by='sample') %>% # if with the relab table
  gather(sample,value=count,2:68)  %>%# if with the count table
  group_by(sample,ASV) %>%
  summarise(count = sum(count)) %>%
  mutate(relab = count/sum(count)*100) %>% 
  #group_by(sample) %>% summarise(srelab = sum(relab)) %>% ungroup() %>%## Check that all samples sum=100
  inner_join(metadata1,by='sample') %>%
  filter(newlayer != 100) %>%
  inner_join(TAX,by="ASV") %>%
  filter(!Phylum %in% top10$Phylum) %>% # for other 
  group_by(Phylum,newlayer) %>%
  dplyr::summarise(relab=sum(relab)) %>%
  filter(relab > 0) %>%
  group_by(newlayer) %>%
  dplyr::summarise(relab=n_distinct(Phylum)) # n_distinct = length(unique(Phylum))
  
  # Phyla per depth
  #layer phyla other
  #5m	    25    15			
  #DCM	  33    23			
  #MESOs	35    25			
  #DSL	  39    29			
  #MESOd	39    29
  #Total 10+30 =40 phyla
  #(14+22)/2 = 18
  #(24+28+28)/3 = 26.67
  
rbPal <- colorRampPalette(brewer.pal(n =10, name = "Paired")[c(1:6,9,10,7,8)]) 
labels.legend <- sub("_.*)","",top10$Phylum) # fix Marinimicorbia
  
OTU.plot <- OTU %>% 
    gather(sample,value=count,2:68) %>%  # if with the count table
    group_by(sample,ASV) %>%
    summarise(count = sum(count)) %>% 
    mutate(relab = count/sum(count)*100)%>% 
    #group_by(sample) %>% summarise(srelab = sum(relab)) %>% ungroup() %>%## Check that all samples sum=100
    inner_join(metadata1,by='sample') %>% 
    filter(newlayer != 100) %>% 
    inner_join(TAX,by="ASV") %>%
    dplyr::select(sample, PhylaClass,relab,newlabel,NewStation,goodorder,season,newlayer) %>% #### select only the columns I will use for the plot
    rename(Phylum = PhylaClass) %>%
  left_join(
    top10 %>% transmute(Phylum, topphylum = Phylum), 
    by = 'Phylum') %>% 
  replace_na(list('topphylum' = 'ZOther')) %>% # Anything that is not in the top 10 phyla will be asigned as Other
  # Sum over samples and phyla
  group_by(sample, topphylum,NewStation,newlabel,goodorder,season,newlayer) %>% 
  dplyr::summarise(relab = sum(relab)) %>% 
  ungroup() 
    #group_by(sample) %>%
    #dplyr::summarise(srelab=sum(relab)) %>%filter(srelab > 100) #check that all groups sum 100%
 
OTU.plot %>% names()
OTU.plot1 <- OTU.plot %>% group_by(topphylum,newlayer) %>%
  summarise(mrelab = mean(relab)) %>% ungroup() #%>%
  #group_by(newlayer) %>%
  #dplyr::summarise(trelab=sum(mrelab))  #check that all groups sum 100%

# Modify some colors
colors <- c(rbPal(length(unique(top10$Phylum))+1),"gray") 
colors[6] <- '#FFD700'
colors[3] <- '#66CDAA'
colors[7] <- '#FA8072'
colors[11] <- '#551A8B'
colors <- colors[-4]

relab1.gg <- OTU.plot1 %>% ggplot(aes(x = newlayer, y = mrelab, fill = topphylum)) +
  geom_bar(stat="identity",position="stack",width=0.9) +
  scale_fill_manual(values=colors, name = "Phylum", labels = c(labels.legend, "Other (30)"))+ 
  scale_y_continuous(expand=c(0,0), name = "mean abundance (%)")+
  scale_x_discrete(name="",labels =c("SURF","DCM","MS","DSL","MD"),expand=expand_scale(mult = c(0.11, 0.11), add = c(0.1, 0)))+
  theme_bw()+
  theme(panel.border = element_blank(),
        panel.grid.minor = element_blank(), # get rid of minor grid
        panel.grid.major = element_blank(),
        text=element_text(size=8),
        #axis.line = element_blank(),
        axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=0.5,color="black",size=6),
        axis.text.y = element_text(color="black",size=6),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        #panel.background = element_rect(fill = "transparent"), # bg of the panel
        #plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
        legend.background = element_rect(fill = "transparent"), # get rid of legend bg
        legend.box.background = element_rect(fill = "transparent") , #get rid of legend panel bg
        legend.margin = margin(0.1,0.1,0.1,0.1,"lines"),
        plot.margin = unit(c(0.5,0.5,0,0.5), "lines"), #top,right,bottom,left
        legend.box=NULL,
        legend.position="none",
        legend.key.size = unit(0.8, "lines"))
        #axis.title.y = element_text(margin = margin(t = 0, r = -1, b = 1, l = -1))) 
# relab1.gg

#######################################
#####           FIG.2 B         #######
#######################################
# log of Odds ratio
phyla.diff <- data.frame(phyla = unique(OTU.plot$topphylum), changeEM= NA, seEM=NA, changeME=NA, seME=NA, abE= NA, abM=NA, probE=NA, probM=NA,odds.ratioME=NA,log.odds.ratioME=NA,odds.ratioEM=NA,log.odds.ratioEM=NA ) ### percent difference between the epi and mesopelagic

for (i in unique(OTU.plot$topphylum)){
  Epi <- OTU.plot %>% filter(newlayer %in% c("5m","DCM") & topphylum == i) %>% unite(id,c(3,6),sep=" ") %>% group_by(id) %>% summarise(meanE = mean(relab))
  Meso <- OTU.plot %>% filter(newlayer %in% c("MESOs","DSL","MESOd") & topphylum == i) %>% unite(id,c(3,6),sep=" ") %>% group_by(id) %>% summarise(meanM = mean(relab))
  valsE <- Epi %>% inner_join(Meso,by="id") %>% group_by(id) %>% summarise(change=100-(meanM*100/meanE))%>% summarise(mean=mean(change), se=std.error(change)) # how much increase in the Epi
    phyla.diff[phyla.diff$phyla==i,]$changeEM <- valsE$mean
    phyla.diff[phyla.diff$phyla==i,]$seEM <- valsE$se
  valsM <- Epi %>% inner_join(Meso,by="id") %>% group_by(id) %>% summarise(change=100-(meanE*100/meanM))%>% summarise(mean=mean(change), se=std.error(change)) # how much increase in the Meso
    phyla.diff[phyla.diff$phyla==i,]$changeME <- valsM$mean
    phyla.diff[phyla.diff$phyla==i,]$seME <- valsM$se
    phyla.diff[phyla.diff$phyla==i,]$abE <- mean(Epi$meanE)
    phyla.diff[phyla.diff$phyla==i,]$abM <- mean(Meso$meanM)
    probE <- (mean(Epi$meanE)/100)/(1-(mean(Epi$meanE)/100))
    probM <- (mean(Meso$meanM)/100)/(1-(mean(Meso$meanM)/100))
    phyla.diff[phyla.diff$phyla==i,]$probE <- probE
    phyla.diff[phyla.diff$phyla==i,]$probM <- probM
    phyla.diff[phyla.diff$phyla==i,]$odds.ratioME <- probM/probE
    phyla.diff[phyla.diff$phyla==i,]$log.odds.ratioME <- log(probM/probE)
    phyla.diff[phyla.diff$phyla==i,]$odds.ratioEM <- probE/probM
    phyla.diff[phyla.diff$phyla==i,]$log.odds.ratioEM <- log(probE/probM)
    
}

phyla.diff$col <- colors
phyla.diff$phyla <- sub("_.*)","",phyla.diff$phyla) 
phyla.diff$phyla <-sub("ZOther","Other",phyla.diff$phyla)
phyla.diff <- phyla.diff %>% arrange(desc(log.odds.ratioME))
phyla.diff$phyla <- factor(phyla.diff$phyla, levels = phyla.diff$phyla)

logOR <- ggplot(data=phyla.diff,aes(x=phyla.diff$phyla, y=log.odds.ratioME,fill=phyla))+
  geom_bar(stat='identity') +
  coord_flip() +
  #scale_x_discrete(limits =sort(phyla.diff$phyla,decreasing=TRUE),expand=c(0,0), labels=sub("ZOther","Other",sort(phyla.diff$phyla,decreasing=TRUE)))+
  scale_fill_manual(values=phyla.diff$col)+
  scale_x_discrete(expand=expand_scale(mult = c(0.17, 0), add = c(0.1, 0)))+
  theme(axis.title.y=element_blank(), 
        #panel.border = element_blank(),
        panel.grid.minor = element_blank(), # get rid of minor grid
        panel.grid.major = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        axis.text.y = element_text(angle = 0, hjust = 1,vjust=0.5,color="black",size=6),
        axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=0.5,color="black",size=6),
        axis.line.x.bottom = element_line(),
        axis.line.y.left = element_line(),
        text=element_text(size=8),
        plot.margin = unit(c(0.5,0.5,0,1.5), "lines"), #top, right,bottom,left
        axis.title.x = element_text(margin = margin(t = 0, r = -1, b = 2, l = -1)),
        legend.position="none") +
  labs(y="log of odds ratio")+
  annotate("text", x=c(0,0),y=c(-3,1),label=c("EPI","MESO"),cex=2, hjust=0) +
   annotate("segment", x=-0.4,y=0,xend=-0.4,yend=-5.2,col="black",cex=0.3,arrow=arrow(length = unit(0.02, "npc")),linejoin = c('round', 'mitre', 'bevel'))+ #EPI line
annotate("segment", x=-0.4,y=0,xend=-0.4,yend=3,col="black",cex=0.3,arrow=arrow(length = unit(0.02, "npc")),linejoin = c('round', 'mitre', 'bevel'))+ #MESO line
  annotate("segment", x=-0.6,y=0,xend=-0.2,yend=0,col="black",cex=0.5)
  
# logOR

#######################################
#####           FIG.2 C         #######
#######################################
#alpha diversity
div <- as.data.frame(diversity(ASVt_rarefy_nocero,"shannon"))
names(div)[1] <- "shannon"
div2 <-data.frame(t(estimateR(ASVt_rarefy_nocero))) #SObs, Chao1, ACE

H <- diversity(ASVt_rarefy_nocero) # from vegan shannon
S <- specnumber(ASVt_rarefy_nocero) ## rowSums(BCI > 0) does the same...
J <- as.data.frame(H/log(S))
div.alpha <- cbind(div,div2,J)
names(div.alpha)[7] <- "J"

div1.gg <- div.alpha %>% 
  tibble::rownames_to_column('sample') %>%
  inner_join(metadata1, by ='sample') %>%
  ggplot(aes(x=newlayer, y = shannon))+
  geom_boxplot(lwd=0.5,fatten=1,outlier.size =-1,color="gray10",fill="gray80")+
  scale_x_discrete(labels=c("SURF","DCM","MS","DSL","MD"),name=NULL,expand=expand_scale(mult = c(0.12, 0.13), add = c(0.1, 0)))+ #This would add space to add Chao1 and J label
  scale_y_continuous(expand=c(0,0),limits = c(3.95,6.4))+
  theme_bw() +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=0.5,color="black",size=6),
        axis.text.y = element_text(color="black",hjust = 0.5,vjust=0.5,size=6),
        panel.background = element_rect(fill = "white"), # bg of the panel
        plot.background = element_rect(fill = "white", color = NA), # bg of the plot
        axis.title.x = element_text(vjust=1, hjust = 0.4),
        axis.title.y = element_text(margin = margin(t = 0, r = 7, b = 1, l = 0)))+
  ylab(label="Shannon index (H)")+
  theme(text=element_text(size=8)) +
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
        plot.margin = unit(c(0.5,0.5,0.7,0.5), "lines"))+ #top, right,bottom,left
  guides(color=FALSE,fill=FALSE)+
  annotate("text",x=c(1:5),y=rep(6.3,5),label=c("a","b","c,d","c,d","b,d"),cex=2,col="black") + # shannon, S.Obs, Chao1
  annotate("text", x=0.55, y =4.2, label="(n)", cex=2)+
  annotate("text",x=c(1:5),y=rep(4.2,5),label=c("590 ± 37","971 ± 36","1066 ± 35","1325 ± 48","1169 ± 77"),cex=2)+
   annotate("text", x=0.55, y =4.055, label="(J)", cex=2)+
 annotate("text",x=c(1:5),y=rep(4.055,5),label=c("0.73" ,"0.80","0.84","0.83","0.80"),cex=2)

# div1.gg

#######################################
#####           FIG.2 D         #######
#######################################

metadata <- metadata[metadata$X %in% names(OTU),]
reorder_idx <- match(colnames(OTU)[2:length(colnames(OTU))],metadata$X)
metadata <- metadata[reorder_idx,]

### organize levels of factors
metadata$layer1 <- factor(metadata$layer1)
metadata$layer1 <- factor(metadata$layer1, levels=c("SURF","DCM","MESOs","FISH","MESOd"))
with(metadata,levels(layer1))

### PCoA
set.seed(200)
OTU.pcoa <- cmdscale(vegdist(ASVt_rarefy_nocero,method="bray"),k=2, eig=TRUE)    ### CMD is principal coordinate analysis (PCoA of euclidean distances is PCA)
Eigenvalues <- eigenvals(OTU.pcoa) 
Variance <- Eigenvalues / sum(Eigenvalues) 
Variance1 <- 100 * signif(Variance[1], 2) 
Variance2 <- 100 * signif(Variance[2], 2)
#ALL= V1=61  V2=16%  sum = 77%

metadata2 <- scale(metadata[,c(13,15,16,21,34,37,47,49,64)])
res2 <- cor.mtest(metadata2, conf.level = .95,adjust.method="Bonferroni")
corrplot(cor(metadata2),method="number",p.mat = res2$p, sig.level = .05, type="upper")
env <- envfit(OTU.pcoa,as.data.frame(metadata2),permu=999,na.rm=TRUE,labels=c("DOC","Temp","Sal",expression(NO[3]^{"2-"}),"Fluor","DIN/DIP","PAb","PVol","%HNA")) #NO3, PO4, SiO3 super correlated

colLy <- c(colors()[76],colors()[50],colors()[637],colors()[631],colors()[132])
pchLy <- c(23,22,24,21,25)

dev.new()
par(#xpd = NA, # switch off clipping, necessary to always see axis labels
    bg = "transparent", # switch off background to avoid obscuring adjacent plots
    mar = c(1, 2, 0.2, 0), #b,l,t,r
    oma=c(0.5,0.5,0.4,0.5)) #b,l,t,r

plot(OTU.pcoa$points,bg = colLy[metadata$layer1],col="gray10",cex=0.6,pch=pchLy[metadata$layer1], xlab="",ylab="", xlim=c(-0.6,0.7),ylim=c(-0.6,0.6), xaxs="i",yaxs="i",las=1,axes=F, lwd=0.4);### cex is point size
box();
axis(1,at=c(seq(-0.6,0.6,0.3),0.7), labels = c(seq(-0.6,0.6,0.3),0.7), padj=-4,tck=-0.01,cex.axis=0.5);
axis(2,at=seq(-0.6,0.6,0.3), labels = FALSE,las=1,tck=-0.01,cex.axis=0.5,hadj=-4);
text(-0.65, seq(-0.6,0.6,0.3),  labels =c("-0.6","-0.3","0","0.3","0.6"), srt = 0, xpd = TRUE,cex=0.5,pos=2,offset=-0.2);
mtext(1,text=bquote(paste("PCoA 1 (" ,.(format(Variance1, digits = 3)),"%)")),line=0.65,cex=0.7,las=1);
text(-0.725,-0.22, srt=90, adj = 0, labels = bquote(paste("PCoA 2 (" ,.(format(Variance2, digits = 3)),"%)")), xpd = TRUE,cex=0.7);
plot(env,labels=NULL,col = "gray40",p.max=0.05,cex=0.001);

sig <- which(env$vectors$pvals < 0.05)

env.lab <- as.data.frame(env$vectors$arrows*sqrt(env$vectors$r))[sig,]
multy <- c(0.65,0.69,0.62,0.62,0.65,0.67,0.67,0.67)
text(env.lab$Dim1*multx,env.lab$Dim2*multy,labels=c("DOC","Temp","Sal",expression(NO[3]^{"2-"}),"Fluor","PAb","PVol","%HNA"),cex=0.6,col='gray20')
text(x=rep(env.lab$Dim1[4]*0.62,2), y=c(-0.15,0.15),c(expression(PO[4]^{"3-"}),expression(SiO[3]^{"2-"})),cex=0.6,col='gray20')
legend(0.5,-0.25,legend=c("SURF","DCM","MS","DSL","MD"),pch=pchLy,col='gray10',cex=0.6,pt.bg=colLy,bty="n",y.intersp=1.8, pt.lwd=0.4,x.intersp=1)

beta <- recordPlot()


#######################################
#####   ARRANGE &   EXPORT      #######
#######################################
#### With reduced relative abundance
gg_F2 = plot_grid(relab1.gg,logOR, div1.gg,beta, labels = c('A','B','C','D'), nrow = 2,ncol=2)
ggsave("Huete-Stauffer_Figure2.pdf",gg_F2,height=110,width=169,units="mm",colormodel="cmyk")