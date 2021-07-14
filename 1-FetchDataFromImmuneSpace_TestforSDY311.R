#Amnah Siddiqa
#To Fetch all Gene expression matrices, HAI matrices and Biosample and HAI annotation files for raw data directory
#check how to configure settings for ImmuneSpaceR @ https://github.com/amnahsiddiqa/VaccinesCohortsMetaAnalysis/wiki/ConnectToImmuneSpaceR



#dependencies
library(ImmuneSpaceR)
library(Rlabkey)
library(dplyr)
library(Biobase)#-->for exprs function 
library(readr)#for write_tsv
library(forcats)
library(tidyr)
library(preprocessCore)
library(plyr)
library(stringi)
library(stringr)
# library(ggplot2) 
# library(reshape2)
 library(plyr)
 library(dplyr)
# library(magrittr)
# library(purrr)
# library(stringr)
# library(readr)
# library(tidyverse) 
# library(preprocessCore)



# Get data from ImmuneSpace -----------------------------------------------


rm( list = ls( envir = globalenv() ), envir = globalenv() ) 

#read the study names --> 25 that i selected for analysis 
#base::load(file="/notebooks/01-PrepareDataForAnalysis/studiesList.rda")
#these are all studies that we can get from Immunespace 
# studiesList<-structure(list(Immport.Name = c("SDY61", "SDY269", "SDY270", "SDY56", "SDY1119",
#                                               "SDY63", "SDY400","SDY404", "SDY520", "SDY640", 
#                                              "SDY212", "SDY312", "SDY112", "SDY315", "SDY67",
#                                              "SDY113", "SDY305",  "SDY296", "SDY301", 
#                                              "SDY1276", "SDY144", "SDY364", "SDY368", "SDY387",
#                                              "SDY372")), class = "data.frame", row.names = c(NA, -26L))

studiesList<-structure(list(Immport.Name = c("SDY1439", "SDY648", "SDY622", "SDY739", "SDY80")),
                                              class = "data.frame", row.names = c(NA, -5L))


studiesList

#"SDY314"reported to ImmPORT for possible error
#SDY311 not in ImmuneSpace ; will process Independantly
#SDY from Anis Larbi ; processed independently 
#

setwd("/Volumes/UNTITLED/VaccinesCohortsMetaAnalysis/April")
studies <- as.vector(studiesList$Immport.Name)
#studies<-studies[1:3]
SetPath=getwd()
mainDir <- paste0(SetPath,"/Data/Raw")

#create /Data/Rawdata directory if does not exist 
if(!dir.exists(mainDir) )dir.create( file.path(mainDir),recursive = TRUE)
  
   

#Function for fetching expression data from ImmuneSpace  
###############Expression Matrix Fetch Function####################
expression_fetch <- function(study) 
{
  #Input:
  # study name e.g. "SDY80"
  #Output:
  # list of expression files associated with data in tsv format
  #
  #Example:
  #
  #>study="SDY80"
  #>expression_fetch(study)
  
  tmp <- CreateConnection(study)
  fmat=tmp$listGEMatrices()
  if (!is.null(fmat))
  {
    flist=as.vector(fmat[,name])
    
    for( i in 1:length(flist)){
      outmat=paste0(flist[i],'_Expression_Matrices.tsv')
      #i have added arguments outputType and annotation because default was unknown and
      #i dont want to assume anything
      obj<-tmp$getGEMatrix(flist[i],outputType = "summary", annotation = "latest")
      #https://www.bioconductor.org/packages/release/bioc/vignettes/ImmuneSpaceR/inst/doc/getGEMatrix.html
      #By default, the returned ExpressionSets have probe names as features (or rows). However, multiple probes often match the same gene and merging experiments from different arrays is impossible at feature level. When they are available, the summary argument allows to return the matrices with gene symbols instead of probes. You should use currAnno set to TRUE to use the latest official gene symbols mapped for each probe, but you can also set this to FALSE to retrieve the original mappings from when the matrix was created.
      #getGEMAtrix --> getGEMatrix(matrixName = NULL, cohortType = NULL, outputType = "summary", annotation = "latest", reload = FALSE, verbose = FALSE)
      #outputType: A character. one of 'raw', 'normalized' or 'summary'. If 'raw', returns an expression matrix of non-normalized values by probe. 'normalized' returns normalized values by probe. 'summary' returns normalized values averaged by gene symbol.
      #annotation: A character. one of 'default', 'latest', or 'ImmSig'. Determines which feature annotation set (FAS) is used. 'default' uses the FAS from when the matrix was generated. latest' uses a recently updated FAS based on the original. 'ImmSig' is specific to studies involved in the ImmuneSignatures project and uses the annotation from when the meta-study's 
      mat<-data.frame(exprs(obj))
      mat_rowname=data.frame(cbind(gene=row.names(mat),mat),stringsAsFactors = FALSE)
      write_tsv(mat_rowname,path=outmat)
      
    }
  }else{print(paste0("No gene expression data in ", study))
  }
}


####################Biosample Fetch Function_SA#######################
#now I am using biosample file from latest release of Immport that i downloaded (SqlImmportData_12122020) 
#and have on local drive 
base::load("/Volumes/UNTITLED/VaccinesCohortsMetaAnalysis/April/biosample.rda")





get_biosample <- function(study) 
{
  #
  # Function to get biosample.txt files for an arbitrary ImmPort study
  # Example:
  #
  # > studies <- c("SDY80","SDY112")
  # > sapply(studies,get_biosample(x))
  #
  
  output <- biosample %>% filter(STUDY_ACCESSION == study)
  
  filename<-paste(study,"biosample.tsv",sep="_")
  print(filename)
  
  
  write_tsv(output,filename)
  msg <- paste0("Processed ",study," to filename: ",filename)
  print(msg)
  
}


####################Fetch HAI Data from ImmuneSPace #######################
# Note : same script can be used to fetch  ImmuneSpace defined types of data as 
# either "hai" or ""neut_ab_titer" but I am right now passing only "hai" in type
## antibody
#antibody_fetch <- function(study="SDY80") {
antibody_fetch <- function(study) {
  #
  # Take a study name and a data type and trues to fetch it
  #
  # INPUT: study - a name of a valid study
  #        
  #
  # OUTPUT: two output files:
  #         1) an annotation file prefixed with the study name
  #         2) a matrix file prefixed with the study name
  #
  # Example:
  # setwd("~/Downloads") # Change this to a working directory
  #
  # sdy135.hai  <- antibody_fetch(study="SDY315")
  # sdy135.neut <- antibody_fetch(study="SDY315")
  #
  # OR could use in a call to sapply as a loop
  #
  # studies <- c("SDY63", "SDY80", "SDY112", "SDY144", "SDY180")
  # neut <- sapply(studies,function(x) antibody_fetch(x)
  #
  
  #type=c("hai","neut_ab_titer")
  #type=("neut_ab_titer")
  #type=("hai")
  type="hai"
  tmp <- CreateConnection(study)
  df  <- tmp$getDataset(type)
  
  # Check to see if anything was returned for the desired data type
  
  if (nrow(df) > 0) {
    virus <- df$virus
    #    unique(virus)
    
    # So as with the fcs data we want to create a variable called observation_id that is based on a combination of
    # the participant id and time point.
    
    ids <- sapply(strsplit(df$`participant_id`,"\\."),`[`,1)
    observation_id <- gsub(" ","",
                           paste(ids,format(df$`study_time_collected`,nsmall=1),
                                 sep="_"))
    
    # Create a new data frame that has the observation_id as the first column.
    
    df.2 <- data.frame(cbind(observation_id=observation_id,df),stringsAsFactors = FALSE)
    colnames(df.2) <- tolower(colnames(df.2))
    
    # Pick out the desired row names as specified by the code sample
    # https://storage.cloud.google.com/immport/Shuzhao_example_SDY80/README.txt
    
    df.3 <- df.2[,c("observation_id","participant_id","age_reported","gender","race","cohort",
                    "study_time_collected","study_time_collected_unit")]
    
    # Get just the the non duplicated rows-->unique(df.3$observation_id)==? SA
    
    df.4 <- df.3[!duplicated(df.3),]
    
    # Write this out as the annotation file
    
    anno.out <- paste(study,type,sep="_")
    
    anno.out <- paste(anno.out,"annotation.tsv",sep="_")
    
    write_tsv(df.4,path=anno.out)
    
    # Now let's write out the matrix portion. We will use df.2 since it has all the data
    # Now, the code at https://storage.cloud.google.com/immport/Shuzhao_example_SDY80/README.txt
    # wants to collape the virus names into a normalized set of four names. We can use the forcats
    # package to simplify this.
    
    #changed all Steve script in here 
    df.2$virus <- as.factor(df.2$virus)
    
    
    #have removed the collapsing script of Steve ; have it in my original collection if needed 
    #remove punctuation marks back slashes and brackets from names and spaces as well and keep rest 
    df.2$hai <- str_replace_all(df.2$virus,"[:punct:]|[:space:]","")
    
    
    
    if (type == "neut_ab_titer") {
      #df.2$neut_ab_titer <- virus.tmp
      
      out <- df.2 %>% dplyr::group_by(neut_ab_titer,observation_id) %>%
        dplyr::summarize(med=median(value_preferred)) %>%
        tidyr::spread(observation_id,med)
      #remove columns with all NAS
      out<-out[, colMeans(is.na(out)) != 1]
      
    } else {
      #df.2$hai <- virus.tmp
      
      out <- df.2 %>% dplyr::group_by(hai,observation_id) %>%
        dplyr::summarize(med=median(value_preferred)) %>%
        tidyr::spread(observation_id,med)
      
      out<-out[, colMeans(is.na(out)) != 1]
      
      #remove columns with all NAS
    }
    
    
    matrix.out <- paste(study,type,sep="_")
    matrix.out <- paste(matrix.out,"data_matrix.txt",sep="_")
    write_tsv(out,matrix.out)
    #    return(out)
  } else {
    # If there is nothing to process keep going
    print(paste("Study",study,"has no data of type",type,sep=" "))
  }
}




####################Main Call#######################

for (i in 1:length(studies)){
  #its study name too so i am directly passing this in required pathnames and study name in arguments
  subDir=studies[i]
  #if study folder exist move into it and fetch data and set back wd
  if (file.exists(file.path(mainDir, subDir))) 
  {
    print(paste0(subDir," exists in ",mainDir," and is a directory"))
    setwd(file.path(mainDir, subDir))
    expression_fetch(subDir)
    get_biosample(subDir)
    antibody_fetch(subDir)
    setwd(mainDir)
  } 
  #else create one , move into it and fetch data and set back wd
  else
  {
    dir.create(file.path(mainDir, subDir))
    setwd(file.path(mainDir, subDir))
    expression_fetch(subDir)
    get_biosample(subDir)
    antibody_fetch(subDir)
    setwd(mainDir)
  }
  #setwd(file.path(mainDir, subDir))
}



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
