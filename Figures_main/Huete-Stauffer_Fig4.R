#######################################
#####          FIGURE 4         #######
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

# Divide merged table in TAX and OTU
merged.rarefy$ASV <- row.names(merged.rarefy)
TAX <-  merged.rarefy[,c(68:75)]
OTU <- merged.rarefy[,c(1:67)]
dim(OTU) #8889

## Transpose table
OTUt <- t(sapply(OTU, as.numeric))
colnames(OTUt) <- row.names(OTU)
OTUt[1:4,1:4]
dim(OTUt) #67 * 8889
OTUt_rarefy_nocero <- OTUt

library(labdsv)
#### cluster is numeric. Assing an number to each MESO layer. Only 2018 spring
unique(metadata$season)

metadata1 <-  metadata %>%  
  filter(layer1 == "FISH") %>% 
  dplyr::select(X,season) %>% 
  filter(X != "CCF_20180321_10_510_F") %>% #filter spring sample out
  dplyr::rename(sample = X)

### only fish samples
OTUt_rarefy_nocero1 <- OTUt_rarefy_nocero[row.names(OTUt_rarefy_nocero) %in% metadata1$sample,]
OTUt_rarefy_nocero2 <- OTUt_rarefy_nocero1[,-which(colSums(OTUt_rarefy_nocero1)<1)]


#### check that the order of the variables is good and the levels assigned from clust are good also
match_idx <- match(metadata1$sample,row.names(OTUt_rarefy_nocero2))
metadata1 <- metadata1[order(match_idx), ]
clust <- as.numeric(factor(metadata1$season)) ### 1=winter, 2=spring, 3=summer
cbind(row.names(OTUt_rarefy_nocero2),metadata1,clust)
nrow(OTUt_rarefy_nocero2) == length(clust)

#### super picky indval, if there is anything wrong, the whole R session crashes
indval.out <- indval(OTUt_rarefy_nocero2,clust)

indval.table <- data.frame(ASV=rep(names(indval.out$maxcls),times=3),
                           indval=c(indval.out$indval[[1]][1:length(indval.out$indval[[1]])],
                                    indval.out$indval[[2]][1:length(indval.out$indval[[1]])],
                                    indval.out$indval[[3]][1:length(indval.out$indval[[1]])]),
                           
                           relfreq=c(indval.out$relfrq[[1]][1:length(indval.out$relfrq[[1]])],
                                     indval.out$relfrq[[2]][1:length(indval.out$relfrq[[1]])],
                                     indval.out$relfrq[[3]][1:length(indval.out$relfrq[[1]])]),
                           
                           relab=c(indval.out$relabu[[1]][1:length(indval.out$relabu[[1]])],
                                   indval.out$relabu[[2]][1:length(indval.out$relabu[[1]])],
                                   indval.out$relabu[[3]][1:length(indval.out$relabu[[1]])]),
                           
                           maxcls=rep(indval.out$maxcls,times=3),
                           pval=rep(indval.out$pval,times=3),
                           FISHseason=c(rep("winter",times=length(indval.out$indval[[1]])),
                                    rep("spring", times=length(indval.out$indval[[2]])),
                                    rep("summer", times=length(indval.out$indval[[3]]))))


write.csv(indval.table,"indval_table_DSL.csv",row.names=FALSE)


#### Analysis start here
indval.table <- read.csv("indval_table_DSL.csv")

head(indval.table)  
indval.table.filter <- indval.table %>%
  filter(pval <= 0.05 & indval>0.5) 
indval.table.filter %>% nrow() #186
sp <- indval.table.filter %>% filter(FISHseason=="spring") ### (pval <= 0.05 & indval>0.5) 18 indicator species (sp10F out : 28 sp) (23 rarefy 82127)
sm <- indval.table.filter %>% filter(FISHseason=="summer") ### (pval <= 0.05 & indval>0.5) 110 indicator species (sp10F out : 96 sp) (87 rarefy 82127)
wt  <- indval.table.filter %>% filter(FISHseason=="winter")### (pval <= 0.05 & indval>0.5) 63 indicator species (sp10F out : 58 sp) (66 rarefy 82127)

TAX %>% filter(ASV %in% sp$ASV) # 26
TAX %>% filter(ASV %in% sm$ASV) # 93
TAX %>% filter(ASV %in% wt$ASV) # 67


### Use the rarefied table for the relative abundances
head(OTU) ### OTUs has singletons removed and 100 m samples
mini <-  metadata %>%  
  dplyr::select(X,season,layer1) %>% 
  filter(X != "CCF_20180321_10_510_F") %>% #filter spring sample out
  dplyr::rename(sample = X)

relab <- OTU %>% tibble::rownames_to_column("ASV") %>% 
  gather(sample,count, 1:ncol(OTU)+1) %>% 
  inner_join(mini,by="sample") %>% ### mini metadata (just season, layer, sample and station)
  filter(layer1=="FISH") %>%
  group_by(ASV,sample,season) %>%
  summarise(mcounts = mean(count)) %>% ungroup() %>%
  filter(mcounts>0) %>%
  group_by(sample,season) %>%
  mutate(relab = mcounts/sum(mcounts)*100) %>%
  ungroup() 

relab %>% group_by(sample) %>% dplyr::summarise(srelab =sum(relab)) ## check they all sum 100


###top 10 Order
Taxa <- indval.table %>%
  filter(pval <= 0.05 & indval>0.5) %>% 
  inner_join(TAX, by="ASV") %>% 
  inner_join(relab,by="ASV") %>%
  #unite(col=Class.up,8:9, sep="/",remove=FALSE) %>% 
  unite(col=Class.up,8:11, sep="/",remove=FALSE) %>%
  filter(relab.y > 0) %>%
  group_by(sample,Class.up,FISHseason) %>%
  dplyr::summarise(relab=sum(relab.y)) %>% ungroup() %>%
  group_by(Class.up) %>%
  dplyr::summarise(mrelab=mean(relab)) %>% ungroup() %>%
  top_n(15,mrelab)


rbPal <- colorRampPalette(brewer.pal(n =11, name = "Paired")[c(2:1,8:7,10:9,6:5,4:3,11)])

### mean relative abundace of the Indval taxa per season
indval.table$FISHseason <- factor(indval.table$FISHseason,levels=c("winter","spring","summer"))
indval.fig <- indval.table %>%
  filter(pval <= 0.05 & indval>0.5) %>% 
  inner_join(TAX, by="ASV") %>% 
  inner_join(relab,by="ASV") %>%
  unite(col=Class.up,8:11, sep="/",remove=FALSE) %>% 
  filter(FISHseason==season) %>%
  filter(relab.y > 0) %>%
  group_by(sample,Class.up,FISHseason) %>%
  dplyr::summarise(relab=sum(relab.y)) %>% ungroup() %>%
  group_by(Class.up,FISHseason) %>%
  dplyr::summarise(mrelab=mean(relab)) %>% ungroup() %>%
  left_join(
    Taxa %>% transmute(Class.up, toptax = Class.up), 
    by = 'Class.up') %>% # Anything that is not in the top 10 phyla will be asigned as Other
  replace_na(list('toptax' = 'Other')) 
  
##set names nicely
  indval.fig$toptax1 <- gsub("/"," - ",indval.fig$toptax)
  indval.fig$toptax1 <- gsub("Marinimicrobia_[(]SAR406_clade[)]","Marinimicrobia",indval.fig$toptax1)
  indval.fig$toptax1 <- gsub("NA","unclassified",indval.fig$toptax1)
  
  indval.fig %>% ggplot(aes(x=FISHseason,y=mrelab,fill=toptax1)) +
    geom_bar(stat="identity",position="stack") +
    scale_fill_manual(values = c(rbPal(length(unique(Taxa$Class.up))),"gray")) +
    scale_y_continuous(expand=c(0,0),limits=c(0,13)) +
    theme_bw()+
    guides(fill=guide_legend(title="Indicator species (Order)"),colour=FALSE)+
    ylab(label= "mean abundance (%)")+
    xlab(label= "season (DSL)")+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black")) +
    theme(strip.text.x = element_text(size = 14, colour = "black", angle = 0),
          legend.spacing = unit(0.3,"cm"),
          legend.key.size = unit(0.3,"cm"),
          legend.title=element_text(size=10),
          legend.text = element_text(color="gray20", size=7)) +
    annotate("text", x=c(1,2,3),y= c(11.8,3.2,5.7),label=c("10.4 ± 1.6", "1.8 ± 0.7","4.3 ± 0.5"),cex=2.8,color="gray30")+
    annotate("text", x=c(1,2,3),y= c(11.1,2.4,4.9),label=c("(67 ASV)", "(26 ASV)","(93 ASV)"),cex=2.8,color="gray30") +
   ggsave("Huete-Stauffer_Figure4.pdf", width=176,height=70,units="mm",colormodel="cmyk")


## relab in each season (get values)
indval.table %>%
  filter(pval <= 0.05 & indval>0.5) %>% 
  inner_join(TAX, by="ASV") %>% 
  inner_join(relab,by="ASV") %>%
  unite(col=Class.up,8:11, sep="/",remove=FALSE) %>% 
  filter(FISHseason==season) %>%
  filter(relab.y > 0) %>%
  group_by(sample,Class.up,FISHseason) %>%
  dplyr::summarise(relab=sum(relab.y)) %>% ungroup() %>%
  group_by(sample,FISHseason) %>%
  summarise(sumst = sum(relab)) %>%
  group_by(FISHseason)%>%
  summarise(mean= mean(sumst),se=std.error(sumst))
  