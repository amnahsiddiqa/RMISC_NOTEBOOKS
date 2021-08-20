# CreateHeatmap -----------------------------------------------------------
rm( list = ls( envir = globalenv() ), envir = globalenv() ) 
#file has Nas 

df.Activity<-read.delim("/Users/siddia/Google Drive/1-FluVaccinesStudyManuscript/Analysis/UnadjustedDataResults/OPLSAnalysis/FullLengthBTMactivities/OPLSR_BTMactivityFileforPlot.tsv",row.names="GeneSet")
df.Activity<-na.omit(df.Activity)

dim(df.Activity)
#So this is a GeneSets (rows) by samples(columns) matrix

