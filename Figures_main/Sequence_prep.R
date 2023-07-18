###################################
#######   SEQUENENCE PREP    ######
###################################
# - Remove unused samples (samples taken at 100 m where excluded from the analysis)
# - Rarefication of remaining samples
# - Remove empty ASVs
# - Join taxa and remaining ASVs

sequence_prep <- function(){
# rRNA table has samples as columns and asvs as rows
rRNA <- read.csv("https://zenodo.org/record/5816123/files/Huete-Stauffer_sequence_count_and_taxonomy.csv?download=1")

# name first column ASV
colnames(rRNA)[1] <- "ASV"

# delete samples of 100 m and a bad sample
rRNA <-  rRNA[,-c(which(names(rRNA)=="CCF_20180321_10_510_F"),grep("_100_",names(rRNA)))]

# remove taxa columns and add ASV as row identifier
row.names(rRNA) <- rRNA$ASV
TAX <- rRNA[,69:75]
ASV <- rRNA[,2:68]

#Transpose table (samples as rows and asvs as columns)
ASVt <- t(sapply(ASV, as.numeric))

# rename cols
colnames(ASVt) <- row.names(ASV)

# Subsampling
# Find the sample with smallest number of  reads
minReads <- min(rowSums(ASVt)) #### Min reads is now 82130
set.seed(1) # make the aleatory sampling always the same

# Rarefication
ASVt_rarefy<-rrarefy(ASVt, minReads) # rarefication
# rowSums(ASVt_rarefy) # check all samples have the same count
ASVt_rarefy_nocero <- ASVt_rarefy[,-(which(colSums(ASVt_rarefy) == 0))] # remove ASVs with cero value (because we loose ASVs by normalizing)

# Match taxa to rarefied table
### revert transposed table to normal way and match ASVs that remain
ASVtt <- t(ASVt_rarefy_nocero)
merged.rarefy <- cbind(ASVtt,TAX[row.names(TAX) %in% row.names(ASVtt),])
# nrow(merged.rarefy) ### 9796 ASVs
# ncol(merged.rarefy) ### 74 cols (67 samples + 7 tax categories)

## remove no Phylum ASVs
# length(which(is.na(merged.rarefy$Phylum))) # 242 ASVs missing Phylum
notaxa <- which(is.na(merged.rarefy$Phylum)) # 242 no Phylum classification

## remove notaxa in ASV table
merged.rarefy <- merged.rarefy[-notaxa,]
# dim(merged.rarefy) # 9554 * 67 samples + 7 tax categories
# 9554/9796*100 # 97.5% of ASVs retained

# Correct singletons after rarefication
# length(which(rowSums(merged.rarefy[1:67]) < 2)) #665 singletons 
singles <- which(rowSums(merged.rarefy[1:67]) < 2)
merged.rarefy <- merged.rarefy[-singles,]
# dim(merged.rarefy) # 8889 * 67 samples + 7 phyla
# 8889/9796*100 # 90.74%

write.csv(merged.rarefy,"Huete-Stauffer_rarefied_sequence_and_taxa.csv",row.names=TRUE)
}