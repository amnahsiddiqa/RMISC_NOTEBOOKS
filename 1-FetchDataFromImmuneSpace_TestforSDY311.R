#Amnah Siddiqa
#To Fetch all Gene expression matrices, HAI matrices and Biosample and HAI annotation files for raw data directory
#check how to configure settings for ImmuneSpaceR @ https://github.com/amnahsiddiqa/VaccinesCohortsMetaAnalysis/wiki/ConnectToImmuneSpaceR






# Remove Manually Extra Cohorts -------------------------------------------



#delete following files
# SDY269 LAIV
# SDY113LAIV
# SDY1119 t2Dold cohort 
# SDY305 IDTIV Cohort 






# Get data for SDY311 -----------------------------------------------


#//////////////////SDY311
#create Directory SDY311
library(limma) 
library(lumiHumanIDMapping) 
library(lumiHumanAll.db) 
library(lumi) 
setwd("/Volumes/UNTITLED/VaccinesCohortsMetaAnalysis/April/Data/Raw/SDY311")
#get Biosample File fof SDY311
get_biosample("SDY311")

# #Read in raw unnormalized Illumina expression data
#downloaded file from ImmPort "SDY311_EXP13635_microarray.tsv"
df.SDY311<- read.delim("/Volumes/UNTITLED/VaccinesCohortsMetaAnalysis/April/Data/Raw/SDY311/SDY311_EXP13635_microarray.703343.tsv", sep = '\t',
                     header = TRUE, row.names = 1, stringsAsFactors = FALSE)
colnames(df.SDY311)

setwd("/Volumes/UNTITLED/VaccinesCohortsMetaAnalysis/April/Data/Raw/SDY311")

#if it was eset i would check this like this 
#table(df.SDY311$QC)

#get the detection pvalues separately
pvalues <- df.SDY311[,grep('PVALUE', colnames(df.SDY311))]


# remove detection p-values,and nbeads from dataframe
df.SDY311 <- df.SDY311[,-grep('PVALUE$|SYMBOL$|NBEADS$', colnames(df.SDY311))]

# Group A Group B Group C 
# 19      26      28 

#kepp the probes for only those values for which we have detection p value in greater than at least 19 subjects 
#as it is the samllest number of individual group od participants
keep <- (rowSums(pvalues < 0.05) >= 19)
table(keep)


#get the matrix and rename columns and rows back to proper annotation 
df.SDY311.Filtered <-df.SDY311[keep,]
pvalues <- pvalues[keep,]
dim(df.SDY311.Filtered)


#log 2 and quantile normalize the data 
df.SDY311.Normalized<- data.matrix(log2(df.SDY311))
df.SDY311.Normalized<-preprocessCore::normalize.quantiles((df.SDY311.Normalized))%>%data.frame()
dim(df.SDY311.Normalized)

colnames(df.SDY311.Normalized)<-colnames(df.SDY311)
row.names(df.SDY311.Normalized)<-row.names(df.SDY311)
df.SDY311.Normalized$feature<-row.names(df.SDY311.Normalized)
colnames(df.SDY311.Normalized)

#get feature data now 
#downloaded file from ImmPort "SDY311_EXP13635_microarray.tsv"
df.fdata<- read.delim("/Volumes/UNTITLED/VaccinesCohortsMetaAnalysis/April/Data/Raw/SDY311/SDY311_EXP13635_microarray.703343.tsv", sep = '\t',
                       header = TRUE, row.names = 1, stringsAsFactors = FALSE)

df.fdata$feature<-row.names(df.fdata)
#keep these columns only
df.fdata<-df.fdata[,c("feature", "SYMBOL")]

df.fdata<-subset(df.fdata, !is.na(SYMBOL))

#now get the gene symbols back on normalized data 
df.SDY311.Normalized<- dplyr::left_join(df.SDY311.Normalized, df.fdata, by="feature")
library("illuminaHumanv4.db") #Get this library if you don't have
# #worked
testfnames<-row.names(df.SDY311.Normalized)
Fnames<-data.frame(Gene=unlist(mget(x = testfnames,envir = illuminaHumanv4SYMBOL)))

Fnames$PROBE_ID<-rownames(Fnames)# dim(test)
Fnames<-subset(Fnames, !is.na(Gene))#14296     2
names(Fnames)[2]<-c("ID_REF")


#Visual inspection 
#boxplot(df.SDY311.Normalized)[1:10]


#Check how many have one:many mapping 
dim(df.SDY311.Normalized)[1]-length(unique(df.SDY311.Normalized$SYMBOL))
#4320
#[1] 28766    75
length(unique(df.SDY311.Normalized$SYMBOL))
# 12840



#keep only the probes which are matched to single gene 

df.SDY311.Normalized$feature<-NULL
df.SDY311.Normalized<- aggregate(. ~ SYMBOL, data = df.SDY311.Normalized, mean)
colnames(df.SDY311.Normalized)

colnames(df.SDY311.Normalized)<-str_replace_all(colnames(df.SDY311.Normalized),"_SIGNAL","")

#get biosamples ids matched 
biosampleData<-read.delim("/Volumes/UNTITLED/VaccinesCohortsMetaAnalysis/April/Data/Raw/SDY311/SDY311-DR38_Subject_2_Illumina_BeadArray.txt")

#concat Arm group and Biosample Accession Id 
biosampleData$newid<-paste0(biosampleData$Biosample.Accession, "_",biosampleData$ARM.Name)

names(df.SDY311.Normalized) <- plyr::mapvalues(names(df.SDY311.Normalized), from = biosampleData$Expsample.Accession , to = biosampleData$newid)


#Groupa<-dm2%>%select(dm2, contains("SYMBOL|Group A"))


GroupA<-df.SDY311.Normalized%>% dplyr::select(matches("SYMBOL|Group A"))
colnames(GroupA)<-(str_replace_all(colnames(GroupA), "_Group A", ""))
colnames(GroupA)[1]<-c("Genes")
write_tsv(GroupA,path="/Volumes/UNTITLED/VaccinesCohortsMetaAnalysis/April/Data/Raw/SDY311/SDY311_GroupA_Expression_Matrices.tsv")

GroupB<-df.SDY311.Normalized%>% dplyr::select(matches("SYMBOL|Group B"))
colnames(GroupB)<-(str_replace_all(colnames(GroupB), "_Group B", ""))
colnames(GroupB)[1]<-c("Genes")
write_tsv(GroupB,path="/Volumes/UNTITLED/VaccinesCohortsMetaAnalysis/April/Data/Raw/SDY311/SDY311_GroupB_Expression_Matrices.tsv")

GroupC<-df.SDY311.Normalized%>% dplyr::select(matches("SYMBOL|Group C"))
colnames(GroupC)<-(str_replace_all(colnames(GroupC), "_Group C", ""))
colnames(GroupC)[1]<-c("Genes")
write_tsv(GroupC,path="/Volumes/UNTITLED/VaccinesCohortsMetaAnalysis/April/Data/Raw/SDY311/SDY311_GroupC_Expression_Matrices.tsv")


#downloaded from ImmuneSpace too
hai<-read.delim("/Volumes/UNTITLED/VaccinesCohortsMetaAnalysis/April/Data/Raw/SDY311/SDY311_hai_data_matrix.txt")
#write_tsv(hai,path="/Volumes/UNTITLED/VaccinesCohortsMetaAnalysis/notebooks/MarchAnalysis/APrilData/Data/Raw/SDY311/SDY311_hai_data_matrix.txt")


# Get Data for Anis Larbi -------------------------------------------------







# Checking few errors -----------------------------------------------------


#////Checking other studies issues
library(ImmuneSpaceR)
sdy1276 <- CreateConnection("SDY1276")
sdy1276$listDatasets()

df<-sdy1276$getDataset("hai")
#all <- CreateConnection("")

unique(df$virus)
#there were three 
#1] "B/Florida/4/2006"   "A/Brisbane/59/2007" "A/Brisbane/10/2007"
