############### SYNOPSIS ###################
# Function to read RNAseq expression file
# This function similar to "function_read_marray.R" expect for the INPUT and OUTPUT (.RData) files
# READ RNAseq: the RNAseq expression data NEEDs normalization: we do the processing of the RNAseq data (Winsorizing + log2)
# The ENSEMBL gene identifiers are CLEANED (duplicates removed), but THIS HAS NO EFFECT ON THE RNAseq data.
# OUTPUT: 2x.RData files:
# 1) "data_rnaseq_expression_unprocessed.RData"
# 2) "data_rnaseq_expression_processed.RData": winsorized and log2 transformed
# OBS: RNAseq has $gene_length set to a DUMMY VARIABLE (123)
############################################
library(dplyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext


rm(list=ls())

wd <- "/Users/pascaltimshel/Dropbox/0_Projects/p_EAv2/git/EAv2/src_BrainCloud"
setwd(wd)


########### SOURCING - **OBS**: HIGH RES! ###########
##source("function_def_stages_BrainSpan.R", echo=TRUE)

order.stages <- c("s1", "s2a", "s2b", "s3a", "s3b", "s4", "s5", "s6", "s7", "s8", "s9", "s10", "s11")

########### FILES ##############
path.data <- "/Users/pascaltimshel/Dropbox/0_Projects/p_EAv2/data_BrainCloud/processed_data"

file.expression_matrix <- paste0(path.data, "/expression_matrix.csv")
file.columns <- paste0(path.data, "/columns_metadata.csv")
file.rows <- paste0(path.data, "/rows_metadata.csv")
#** No gene length file


########### READ columns file ###########
df.columns <- read.csv(file.columns,h=T)
names(df.columns)
str(df.columns)
### *OBS* IMPORTANT: substitute underscore ("_") in donor_id with "." ###
df.columns$donor_id <- gsub("_", ".", df.columns$donor_id)

#### Add stage column
## SKIPPED
#### Sort factor levels of "stage"
df.columns$stage <- factor(df.columns$stage, levels(df.columns$stage)[match(order.stages, levels(df.columns$stage))])





########### READ row file ###########
df.rows <- read.csv(file.rows,h=F)
colnames(df.rows)
colnames(df.rows) <- c("ensembl_gene_id", "probeID")
dim(df.rows) # [1] 23354     2

########### READ gene_length file ###########
###*** NO GENE LENGTH FILE

########### READ AND MANIPULATE expression file ###########
### ** THIS TAKES SOME TIME ** ###
df.expression_matrix <- read.csv(file.expression_matrix,h=F) # HEADER FALSE.
dim(df.expression_matrix) # [1] 23354   269

############################# ** BrainCloud switch ** ###########################
df.expression_matrix.clean <- df.expression_matrix
df.rows.clean <- df.rows


############################# ** NEW 12/01/2014 **###########################
### Setting column names - must be done first!

#colnames(df.expression_matrix.clean) <- with(df.columns, paste(donor_id, structure_acronym, stage, age, gender, sep="_")) # NEW - added 
colnames(df.expression_matrix.clean) <- with(df.columns, paste(donor_id, structure, stage, round(age_raw), sex, sep="_")) # BrainCloud
colnames(df.expression_matrix.clean)


### Setting ensemblID
df.expression_matrix.clean$ensembl_gene_id <- df.rows.clean$ensembl_gene_id


str(df.expression_matrix.clean,list.len=Inf)


############################### MANIPULAION - ALL GENES - full ###########################
### Melting dataframe
df.expression_matrix.clean.melt <- melt(df.expression_matrix.clean, id=c("ensembl_gene_id"))
#df.expression_matrix.clean.melt <- melt(df.expression_matrix.clean, id=c("ensembl_gene_id", "gene_type"))
head(df.expression_matrix.clean.melt)

### Creating new variables from string
variable_split <- strsplit(as.character(df.expression_matrix.clean.melt$variable), "_") # THIS takes a LOOONG time
df.expression_matrix.clean.melt$donor_id <- as.factor(sapply(variable_split, "[[", 1))
df.expression_matrix.clean.melt$structure_acronym <- as.factor(sapply(variable_split, "[[", 2))
df.expression_matrix.clean.melt$stage <- as.factor(sapply(variable_split, "[[", 3))
### sorting stage levels
df.expression_matrix.clean.melt$stage <- with(df.expression_matrix.clean.melt, factor(stage, levels(stage)[match(order.stages, levels(stage))]))
head(df.expression_matrix.clean.melt)
############################# ** NEW 12/01/2014 
### adding age
#df.expression_matrix.clean.melt$age <- as.factor(sapply(variable_split, "[[", 4))
df.expression_matrix.clean.melt$age <- as.numeric(sapply(variable_split, "[[", 4)) # *EDIT VERSION for BrainCloud*
levels(df.expression_matrix.clean.melt$age) # check the order
# soring age factor levels - *OUTCOMMENTED BrainCloud*
#df.expression_matrix.clean.melt$age <- with(df.expression_matrix.clean.melt, factor(age, levels(age)[match(order.age, levels(age))]))
#levels(df.expression_matrix.clean.melt$age)

head(df.expression_matrix.clean.melt)
### adding gender
df.expression_matrix.clean.melt$gender <- as.factor(sapply(variable_split, "[[", 5)) 


############################# Setting pre and post natal stages ###########################
## adding stage_natal (prenatal vs post-natal)
df.expression_matrix.clean.melt$natal <- NA
df.expression_matrix.clean.melt[df.expression_matrix.clean.melt$stage %in% c("s1","s2a","s2b","s3a","s3b","s4","s5"), "natal"] <- "prenatal"
df.expression_matrix.clean.melt[df.expression_matrix.clean.melt$stage %in% c("s6","s7","s8","s9","s10","s11"), "natal"] <- "postnatal"
df.expression_matrix.clean.melt$natal <- as.factor(df.expression_matrix.clean.melt$natal)
str(df.expression_matrix.clean.melt)



df.expression_matrix.clean.unprocessed <- df.expression_matrix.clean
df.expression_matrix.clean.melt.unprocessed <- df.expression_matrix.clean.melt
################################# SAVING *UNPROCESSED* DATA ###################################
#save(df.expression_matrix.clean.unprocessed, df.expression_matrix.clean.melt.unprocessed, file="RData/data_braincloud_expression_unprocessed.RData")
#load(file="RData/data_braincloud_expression_unprocessed.RData")
#df.expression_matrix.clean <- df.expression_matrix.clean.unprocessed
#df.expression_matrix.clean.melt <- df.expression_matrix.clean.melt.unprocessed


################################# DIAGNOSING DATA - **BrainCloud Specific** #####################################

### Checking for NA
anyNA(df.expression_matrix.clean) # --> MUST BE FALSE


### Histogram of number of probes per ENSG gene
x <- df.expression_matrix.clean %>% 
  select(ensembl_gene_id) %>%
  group_by(ensembl_gene_id) %>%
  summarise(
    count = n()
  )
ggplot(x, aes(x=count)) + geom_histogram(binwidth=1)
ggsave(file="diagnostics_histogram_number_of_probes_per_ENSG_gene.pdf", w=12, h=8)

### Distribution of expression values
str(df.expression_matrix.clean.melt)
#x1 <- df.expression_matrix.clean.melt %>% filter(donor_id=="HB.16.14") # TEST
#p <- ggplot(df.expression_matrix.clean.melt, aes(x=value, color=donor_id)) + geom_density()
#p <- p + labs(title="Distribution of expression values over donor_id (color=donor_id)", x="Expression value")
#ggsave(file="diagnostics_donor_id_expression_values.pdf", w=12, h=8)


df.tmp.stats <- df.expression_matrix.clean.melt %>% 
  group_by(donor_id) %>%
  #group_by(donor_id, stage, gender, natal, structure_acronym) %>%
  summarise(
    #n.donor_id = n_distinct(donor_id),
    n.stage = n_distinct(stage),
    n.gender = n_distinct(gender),
    n.natal = n_distinct(natal)
  )
# We see that each donor is only represented in one stage. This makes sense, since they use post-mortem brains


################################# PROCESSING DATA - calculating MEDIAN over probes #####################################
# ** NEW DATA FRAME - copy
str(df.expression_matrix.clean.melt)

system.time(df.expression_matrix.clean.melt.transform <- df.expression_matrix.clean.melt %>% 
  group_by(donor_id, ensembl_gene_id) %>%
  summarise(
    variable = unique(variable),
    value = median(value),
    stage = unique(stage),
    age = unique(age),
    gender = unique(gender),
    natal = unique(natal)
  )
)
# elapsed --> 911.817 seconds

p <- ggplot(df.expression_matrix.clean.melt.transform, aes(x=value, color=donor_id)) + geom_density()
p <- p + labs(title="Distribution of expression values over donor_id (color=donor_id)", x="MEDIAN Expression value")
ggsave(file="diagnostics_donor_id_expression_values_AFTER_MEDIAN_NEW.pdf", w=12, h=8)


######## COPYING BACK AGAIN
df.expression_matrix.clean.melt <- df.expression_matrix.clean.melt.transform

################################# SAVING *PROCESSED* DATA ###################################
### Notice that the "unprocessed" and "processed" data has different variable names
save(df.expression_matrix.clean, df.expression_matrix.clean.melt, file="RData/data_braincloud_expression_processed.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt


#str(df.expression_matrix.clean)
#str(df.expression_matrix.clean.melt)

#levels(df.expression_matrix.clean.melt$age)
#levels(df.expression_matrix.clean.melt$stage)


