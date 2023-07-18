###################################
#######     DATA DOWNLOAD    ######
###################################

# The main tables of ASVs and associated metadata are stored in a zenodo repository

rRNA <- read.csv("https://zenodo.org/record/5816123/files/Huete-Stauffer_sequence_count_and_taxonomy.csv?download=1")
metadata <- read.csv("https://zenodo.org/record/5816123/files/Huete-Stauffer_Metadata.csv?download=1")


