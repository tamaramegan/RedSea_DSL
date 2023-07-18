#######################################
#####          FIGURE 6         #######
#######################################
# Created by Tamara Huete Stauffer tamara.huete@gmail.com

# - Access data and process the sequence table as in Sequence_prep.R
# - The starting point of this script is the table 'merged.rarify'

#######################################
#####      LIBRARIES           #######
#######################################
source("Requirements.R")


#######################################
#####     SEQ & METADATA PREP   #######
#######################################
# Run Sequence_prep.R to get merged.rarefy
source("Sequence_prep.R")
sequence_prep()
merged.rarefy <- read.csv("Huete-Stauffer_rarefied_sequence_and_taxa.csv",row.names=1)

OTU <- merged.rarefy[,c(1:67)]
OTUt <- t(sapply(OTU, as.numeric))
colnames(OTUt) <- row.names(OTU)

### read metadata
metadata <- read.csv("https://zenodo.org/record/5816123/files/Huete-Stauffer_Metadata.csv?download=1")
metadata$sample <- metadata$X
names(metadata)
mini <- metadata[c(5,42,63,65)]
head(mini)

#######################################
#####      VENN DIAGRAM DATA    #######
#######################################
#### MESOPELAGIC
#### BY STATION 
#### ONLY SPRING

library(devtools)
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram3.r")
source_url("http://raw.github.com/nielshanson/mp_tutorial/master/downstream_analysis_r/code/venn_diagram4.r")

DSL <- as.data.frame(OTUt) %>% ### it is not DSL just all data #### rarefied, but all 0 counts eliminated (and before eliminaing 100, I checked for singletons)
  tibble::rownames_to_column('sample') %>%
  inner_join(mini, by ='sample') 

ncol(DSL)
DSL[1:4,c(1:4,8889:8893)]
 
#4WAY VENN
### For spring, if I do a 4way Venn diagram I will be able to see also how much meso contributes to sinking
venn_4way <- data.frame(ASV=NA, season=NA,station=NA,layer=NA)
i="spring"
for (a in unique(DSL[DSL$season==i,]$NewStation)){
    
    ## WINTER (singles out, not rarefied)
DSLw <- DSL %>% filter(season ==i & NewStation ==a & layer1 %in% c("SURF","DCM")) %>%
  tibble::column_to_rownames("sample") %>% select(-"season",- "layer1", -"NewStation") %>% colSums() 
wt.surf <- names(DSLw[which(DSLw>1)])
DSLd <- DSL %>% filter(season ==i & NewStation ==a & layer1 == "FISH") %>%
  tibble::column_to_rownames("sample") %>% select(-"season",- "layer1", -"NewStation") %>% colSums() 
wt.dcm <- names(DSLd[which(DSLd>1)])
DSLf <- DSL %>% filter(season ==i & NewStation ==a & layer1 %in% "MESOs") %>%
  tibble::column_to_rownames("sample") %>% select(-"season",- "layer1", -"NewStation") %>% colSums() 
wt.dsl <- names(DSLf[which(DSLf>1)])
DSLf <- DSL %>% filter(season ==i & NewStation ==a & layer1 %in% "MESOd") %>%
  tibble::column_to_rownames("sample") %>% select(-"season",- "layer1", -"NewStation") %>% colSums() 
wt.md <- names(DSLf[which(DSLf>1)])

###For some reason I can't save the graphs, but doing this I can
pdf(paste("Venn_diagram_4layers_",i,"_",a,".pdf",sep=""))
  Wsdf <- venn_diagram4(wt.surf,wt.dcm,wt.dsl,wt.md, "EPI","FISH","MESOs","MESOd", colors= c("cyan","darkgreen","orange","blue"))
dev.off()
  #str(Wsdf)
  
  venn_winter <- data.frame(
ASV=c(Wsdf$"EPI_only",Wsdf$"EPI_FISH",Wsdf$"EPI_MESOs", Wsdf$"EPI_MESOd", Wsdf$"FISH_only",Wsdf$"FISH_MESOs",Wsdf$"FISH_MESOd", Wsdf$"MESOs_only",Wsdf$"MESOs_MESOd",Wsdf$"MESOd_only",Wsdf$"EPI_FISH_MESOs",Wsdf$"FISH_MESOs_MESOd",Wsdf$"EPI_FISH_MESOd",Wsdf$"EPI_FISH_MESOs_MESOd",Wsdf$"EPI_MESOs_MESOd"),
season=i,
station=a,
layer=c(rep("EPI_only",length(Wsdf$"EPI_only")),
        rep("EPI_FISH",length(Wsdf$"EPI_FISH")),
        rep("EPI_MESOs",length(Wsdf$"EPI_MESOs")),
        rep("EPI_MESOd",length(Wsdf$"EPI_MESOd")),
        rep("FISH_only",length(Wsdf$"FISH_only")),
        rep("FISH_MESOs",length(Wsdf$"FISH_MESOs")),
        rep("FISH_MESOd",length(Wsdf$"FISH_MESOd")),
        rep("MESOs_only",length(Wsdf$"MESOs_only")),
        rep("MESOs_MESOd",length(Wsdf$"MESOs_MESOd")),
        rep("MESOd_only",length(Wsdf$"MESOd_only")),
        rep("EPI_FISH_MESOs",length(Wsdf$"EPI_FISH_MESOs")),
        rep("EPI_MESOs_MESOd",length(Wsdf$"EPI_MESOs_MESOd")),
        rep("FISH_MESOs_MESOd",length(Wsdf$"FISH_MESOs_MESOd")),
        rep("EPI_FISH_MESOd",length(Wsdf$"EPI_FISH_MESOd")),
        rep("EPI_FISH_MESOs_MESOd",length(Wsdf$"EPI_FISH_MESOs_MESOd"))))
        
    
  venn_4way <- rbind(venn_4way,venn_winter)
  }
write.csv(venn_4way,"Venn_4way_bystation.csv",row.names=FALSE)
venn_4way <- read.csv("Venn_4way_bystation.csv")

venn_4way  <- venn_4way[-1,]
head(venn_4way)

#Number of ASVs of each layer
total <- venn_4way %>% plyr::count(c("season","station"))
fishly <- unique(as.character(venn_4way$layer))[grep('FISH', unique(as.character(venn_4way$layer)))]
mesoly <- unique(as.character(venn_4way$layer))[grep('MESO', unique(as.character(venn_4way$layer)))]

total.fish <- venn_4way %>% filter(layer %in% fishly) %>% plyr::count(c("season","station")) 
onlyF <- venn_4way %>% filter(layer %in% c("FISH_only")) %>% plyr::count(c("season","station")) 
meso.shared <- venn_4way %>% filter(layer %in% c("FISH_MESOs_MESOd")) %>% plyr::count(c("season","station")) 
meso.s <- venn_4way %>% filter(layer %in% c("FISH_MESOs")) %>% plyr::count(c("season","station")) 
meso.d <- venn_4way %>% filter(layer %in% c("FISH_MESOd")) %>% plyr::count(c("season","station")) 
sinking <- venn_4way %>% filter(layer %in% c("EPI_FISH","EPI_FISH_MESOs","EPI_FISH_MESOd","EPI_FISH_MESOs_MESOd")) %>% plyr::count(c("season","station"))
total.meso <- venn_4way %>% filter(layer %in% unique(c(mesoly,fishly))) %>% plyr::count(c("season","station")) 
meso <- venn_4way %>% filter(layer %in% c("FISH_MESOs","FISH_MESOd","FISH_MESOs_MESOd")) %>% plyr::count(c("season","station"))

sdf.fish <- as.data.frame(total.fish)
colnames(sdf.fish)[3] <- 'total.fish'
### form venn 4way :)
sdf.fish$shared.meso.n <- ifelse(sdf.fish$season=="spring",meso.shared$freq, NA)
sdf.fish$shared.meso.per <- ifelse(sdf.fish$season=="spring",meso.shared$freq/total.fish$freq*100, NA)
sdf.fish$mesos.n <- ifelse(sdf.fish$season=="spring",meso.s$freq, NA)
sdf.fish$mesos.per <- ifelse(sdf.fish$season=="spring",meso.s$freq/total.fish$freq*100, NA)
sdf.fish$mesod.n <- ifelse(sdf.fish$season=="spring",meso.d$freq, NA)
sdf.fish$mesod.per <- ifelse(sdf.fish$season=="spring",meso.d$freq/total.fish$freq*100, NA)
sdf.fish$uniq.n <- ifelse(sdf.fish$season=="spring",onlyF$freq, NA)
sdf.fish$uniq.per <- ifelse(sdf.fish$season=="spring",onlyF$freq/total.fish$freq*100, NA)
sdf.fish$sink.n <- ifelse(sdf.fish$season=="spring",sinking$freq, NA)
sdf.fish$sink.per <- ifelse(sdf.fish$season=="spring",sinking$freq/total.fish$freq*100, NA)
sdf.fish$meso.n <- ifelse(sdf.fish$season=="spring",meso$freq, NA)
sdf.fish$meso.per <- ifelse(sdf.fish$season=="spring",meso$freq/total.fish$freq*100, NA)
sdf.fish$total.meso <- ifelse(sdf.fish$season=="spring",total.meso$freq, NA)

write.csv(sdf.fish,"div_partitioning_summary_table_rar_spring.csv",row.names=FALSE)

#######################################
#####        BAR GRAPH          #######
#######################################
div <- read.csv("div_partitioning_summary_table_rar_spring.csv")
head(div)
nrow(div)


#### The mesopelagic values
####  ASVs
ASVper <- data.frame(per=c(div$sink.per,div$uniq.per,div$shared.meso.per,div$mesos.per,div$mesod.per), season=rep(div$season, 5),partition =rep(c("sinking","unique","meso shared","meso shallow","meso deep"),each=5))
data2 = ASVper %>% filter(season=="spring") %>% group_by(season,partition) %>% 
  summarise(mper = mean(per), se=std.error(per)) 
data2$partition <- factor(data2$partition, levels=c("sinking","unique","meso shallow","meso deep","meso shared"))
data2 <- data2 %>% arrange(partition)

### ASVs contribution
# Meso common
sum(data2[data2$partition %in% c('meso shallow', 'meso deep', 'meso shared'),]$mper) #72.06366
# sinking
data2[data2$partition=='sinking',]$mper #6.9453
# unique
data2[data2$partition=='unique',]$mper #20.99104

### total ASVS in DSL
div %>% summarize(DSL= mean(total.fish), se =std.error(total.fish)) #1189.8 +- 119.0174

### Shared ASVs with higher ab in DSL (#from next section)
10/164 # 0.06097561% of meso shallow
45/160.4 #0.2805486% of meso deep
145/518.6 #0.2795989% de shared.meso


## ASVs DSL in MESO
fish.meso <- mean((div$uniq.n/div$total.meso)*100) ## 13.68%
fish.dsl <- mean((div$uniq.n/div$total.fish)*100) ## 20.99%
high.meso <- mean((10+45+145)/div$total.meso)*100 #10.335
# 13.68 + 10.335 = 24.015 % DSL influence in Meso


#(sinking, unique, meso shallow, meso deep, meso shared)
data3 <- data2$mper*c(1,1,0.06097561,0.2805486,0.2795989)
sum(data3[3:5]) # 17.16% (+21% = 38.16) # total  higher 

data4 <- c(data2$mper[1:2],data3[3],data2$mper[3]-data3[3],data3[4],data2$mper[4]-data3[4],data3[5],data2$mper[5]-data3[5])



###############################################
# Fig 6 with colors and Venn diagram values
df <- data.frame(mper=round(data4,1),partition=c("sinking","unique","mesos-dsl","mesos","mesod-dsl","meso deep","mesosh-dsl","meso shared"), ly=c("EPI","DSL","MS","MS","MD","MD","MESO","MESO"), colLy=as.character(c(1,2,6,3,6,4,6,5)))
df$ly <- factor(df$ly,levels=rev(c("EPI","MS","DSL","MD","MESO")))
a=-0.1
b=-0.8

ASV <- ggplot(df, aes(x=ly, y=mper, fill=colLy)) + 
  geom_bar(stat='identity') +
  #geom_errorbar(aes(ymin=mper-se, ymax=mper+se), width=.2,
  #               position=position_dodge(.9)) +
  coord_flip() +
  # colors epi, dsl, ms, md, dsl+
  #scale_fill_manual(values=c(colors()[516],colors()[504],colors()[617],colors()[564],"gray60",colors()[621]))+
  scale_fill_manual(values=c(colors()[104],colors()[503],colors()[400],colors()[128],"gray60",colors()[621]))+
  scale_y_continuous(limits=c(0,50), expand=expand_scale(mult = c(0.02, 0), add = c(0.02, 0)))+
  scale_x_discrete(expand=expand_scale(mult = c(0.1, 0), add = c(0.1, 0)),labels=rev(c("DSL+EPI","DSL+MS","DSL","DSL+MD","DSL+MS+MD")))+
  #theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        #panel.background = element_rect(fill = "transparent"),
        #plot.background = element_rect(fill = "transparent"),
        panel.background = element_rect(fill = "white"),
        #axis.line.x.bottom = element_line(color="white"),
        axis.line.x.bottom = element_line(),
        #axis.line.y.left = element_line(color="white"),
        axis.line.y.left = element_line(),
        axis.title.y=element_blank(), 
        axis.text.y = element_text(angle = 0, hjust = 1,vjust=0.5,color="black",size=6.5,margin=margin(0,1,0,0)),
        axis.text.x = element_text(angle = 0, hjust = 0.5,vjust=0.5,color="black",size=6.5),
        #text=element_text(size=8,color="white"),
        text=element_text(size=8),
        plot.margin = unit(c(0.5,0.5,0.3,0.5), "lines"), #top, right,bottom,left
        #axis.title.x = element_text(margin = margin(t = 0, r = -1, b = 2, l = -1)),
        axis.title.x = element_text(margin = unit(c(1, 0, 0, 0), "mm")),
        #axis.ticks = element_line(color="white"),
        legend.position="none") +
  labs(y="ASVs (%)")+
  #annotate("text",x=-10, y=30,label="Species (%)")+
  # (s,dsl,dsl-ms,ms,dsl-md,md,dsl-m,m)
annotate("text", x=(as.numeric(df$ly) + c(0,0,b,a,b,a,0,0)),y=(df$mper +c(2,0,2,2+0.4,2,2+3.7,0,13.4*2))*c(1.4,0.5,1.3,1.3,1.3,1.3,0.5,0.5),label=c(paste(df$mper[1:2], "%",sep=""),"",paste(df$mper[4], "%",sep=""),"",paste(df$mper[6:8], "%",sep="")),cex=2.4, hjust=0.5, col="black")+
  annotate("text", x=4.15,y=(df$mper[4]+2+0.4)*1.3,label=paste(df$mper[3], "%",sep=""),cex=2.4, hjust=0.5, col=colors()[621])+
annotate("text", x=2.15,y=(df$mper[6]+2+3.7)*1.3,label=paste(df$mper[5], "%",sep=""),cex=2.4, hjust=0.5, col=colors()[621]) 
  
ASV

#######################################
#####   ARRANGE &   EXPORT      #######
#######################################
# - Fig A: Venn diagram (image meade outside of R, saved as png)
# - Fig B: Bar graph
p2 <- cowplot::ggdraw() + cowplot::draw_image("Venn4way_Fig6A.png", scale = 0.95)
ASV_Venn <- cowplot::plot_grid(p2, ASV, labels = c('A','B'),nrow=1,ncol=2)
ggsave("Huete-Stauffer_Figure6.pdf",ASV_Venn,height=60,width=169,units="mm",colormode="cmyk")