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
library(plyr)
library(ggplot2)
library(reshape2)
library(tools) # for file_path_sans_ext, path.expand()

rm(list=ls())

wd <- path.expand("~/Dropbox/0_Projects/p_EAv2/git/EAv2/src")
setwd(wd)

########### SOURCING - **OBS**: HIGH RES! ###########
#source("function_def_stages.R", echo=TRUE)
source("function_def_stages.corrected.R", echo=TRUE) # New march 2015

#### New: this will define the order of the age factor
#### I created this vector partial manually. 
##### I have TESTED it and it is correct FOR RNASEQ (needs validation in Microarray)
#### BE CAREFUL IF UPDATEING THE DATA SET (including new years)
order.age <- c("8 pcw","9 pcw","12 pcw","13 pcw","16 pcw","17 pcw","19 pcw","21 pcw","24 pcw","25 pcw","26 pcw","35 pcw","37 pcw",
               "4 mos","10 mos",
               "1 yrs","2 yrs","3 yrs","4 yrs","8 yrs",
               "11 yrs","13 yrs","15 yrs","18 yrs","19 yrs",
               "21 yrs","23 yrs","30 yrs","36 yrs","37 yrs","40 yrs")


########### FILES ##############
file.expression_matrix <- "../data/141031/rnaseq/expression_matrix.csv"
file.columns <- "../data/141031/rnaseq/columns_metadata.csv"
file.rows <- "../data/141031/rnaseq/rows_metadata.csv"
#** No gene length file


########### READ columns file ###########
df.columns <- read.csv(file.columns,h=T,row.names=1)

#### Add stage column
df.columns$stage <- as.factor(sapply(df.columns$age, function(x) {names(stages)[sapply(stages, function(stage) x %in% stage)]}))
#### Sort factor levels of "stage"
df.columns$stage <- factor(df.columns$stage, levels(df.columns$stage)[match(order.stages, levels(df.columns$stage))])





########### READ row file ###########
df.rows <- read.csv(file.rows,h=T,row.names=1)

########### READ gene_length file ###########
###*** NO GENE LENGTH FILE

########### READ AND MANIPULATE expression file ###########
### ** THIS TAKES SOME TIME ** ###
df.expression_matrix <- read.csv(file.expression_matrix,h=F,row.names=1) # HEADER FALSE
### Removing duplicates
df.expression_matrix.clean <- df.expression_matrix[!(duplicated(df.rows$ensembl_gene_id) | duplicated(df.rows$ensembl_gene_id, fromLast = TRUE)), ]
df.rows.clean <- df.rows[!(duplicated(df.rows$ensembl_gene_id) | duplicated(df.rows$ensembl_gene_id, fromLast = TRUE)), ]
# Print the sum of duplicates
n_duplicates <- sum(duplicated(df.rows$ensembl_gene_id) | duplicated(df.rows$ensembl_gene_id, fromLast = TRUE))
print(paste("Number of ENSEMBL Gene IDs (duplicates) removed:", n_duplicates))

############################# ** NEW 12/01/2014 **###########################
### Setting column names - must be done first!
#colnames(df.expression_matrix.clean) <- with(df.columns, paste(donor_id, structure_acronym, stage, sep="_")) # BEFORE 12/01/2014
colnames(df.expression_matrix.clean) <- with(df.columns, paste(donor_id, structure_acronym, stage, age, gender, sep="_")) # NEW - added 
colnames(df.expression_matrix.clean)



### Setting ensemblID
df.expression_matrix.clean$ensembl_gene_id <- df.rows.clean$ensembl_gene_id

### Setting gene_length - DUMMY VARIBLE FOR RNAseq
df.expression_matrix.clean$gene_length <- 123

str(df.expression_matrix.clean,list.len=Inf)


############################### MANIPULAION - ALL GENES - full ###########################
### Melting dataframe
df.expression_matrix.clean.melt <- melt(df.expression_matrix.clean, id=c("ensembl_gene_id", "gene_length"))
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
df.expression_matrix.clean.melt$age <- as.factor(sapply(variable_split, "[[", 4))
levels(df.expression_matrix.clean.melt$age) # check the order
# soring age factor levels
df.expression_matrix.clean.melt$age <- with(df.expression_matrix.clean.melt, factor(age, levels(age)[match(order.age, levels(age))]))
levels(df.expression_matrix.clean.melt$age)
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
#save(df.expression_matrix.clean.unprocessed, df.expression_matrix.clean.melt.unprocessed, file="RData/data_rnaseq_expression_unprocessed_full_res.corrected-03-23-2015.RData")
#load(file="RData/data_rnaseq_expression_unprocessed_full_res.RData")
#df.expression_matrix.clean <- df.expression_matrix.clean.unprocessed
#df.expression_matrix.clean.melt <- df.expression_matrix.clean.melt.unprocessed

################################# PROCESSING DATA #####################################
########## "Winsorizing" samples higher than 50
# **** NEW DATA FRAME
df.expression_matrix.rnaseq.clean.melt.transform <- df.expression_matrix.clean.melt
df.expression_matrix.rnaseq.clean.melt.transform[df.expression_matrix.rnaseq.clean.melt.transform$value > 50, "value"] <- 50
######### Log2(RPKM+1) transforming
df.expression_matrix.rnaseq.clean.melt.transform$value <- log2(df.expression_matrix.rnaseq.clean.melt.transform$value+1)
######## COPYING BACK AGAIN
df.expression_matrix.clean.melt <- df.expression_matrix.rnaseq.clean.melt.transform

################################# SAVING *PROCESSED* DATA (winsorizing, log2) ###################################
### Notice that the "unprocessed" and "processed" data has different variable names
#save(df.expression_matrix.clean, df.expression_matrix.clean.melt, file="RData/data_rnaseq_expression_processed_full_res-03-23-2015.RData") # df.expression_matrix.clean, df.expression_matrix.clean.melt


#str(df.expression_matrix.clean)
#str(df.expression_matrix.clean.melt)

#levels(df.expression_matrix.clean.melt$age)
#levels(df.expression_matrix.clean.melt$stage)


