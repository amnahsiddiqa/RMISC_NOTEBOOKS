# CreateHeatmap -----------------------------------------------------------
rm( list = ls( envir = globalenv() ), envir = globalenv() ) 
#file has Nas 

df.Activity<-read.delim("/Users/siddia/Google Drive/1-FluVaccinesStudyManuscript/Analysis/UnadjustedDataResults/OPLSAnalysis/FullLengthBTMactivities/OPLSR_BTMactivityFileforPlot.tsv",row.names="GeneSet")
df.Activity<-na.omit(df.Activity)

dim(df.Activity)
#So this is a GeneSets (rows) by samples(columns) matrix

#checking last name of column 
colnames(df.Activity)[ncol(df.Activity)]
#chacking first name of column 
colnames(df.Activity)[1:2]

row.names(df.Activity)[1]


#cor function uses columns to compute correlations in a matrix ; threfore passing transposed matrix
# Pairwise correlation between Gene Sets (rows)
cols.cor <- cor(t(df.Activity), use = "pairwise.complete.obs", method = "pearson")

#get the dissimilarity matrix 
df.dist= as.dist(1 - cols.cor)


#Get the clustering done 
data.tree <- hclust(df.dist, method="complete")


#plot the tree/dendrogram 
#par(cex=0.3, mar=c(5, 8, 4, 1))
par( mar=c(5, 8, 4, 1))

plot(data.tree,labels=FALSE)
par(cex=1)
abline(h=0.5, col="red")
# par(cex=1)
# title(xlab="xlab", ylab="ylab", main="main")
# axis(2)

df.clusters<-cutree(data.tree , h = 0.5) %>%data.frame()
df.clusters$GeneSet<-row.names(df.clusters)
colnames(df.clusters)<-c("ClusterNum", "GeneSet")


#how many unique clusters are in there
length(unique(df.clusters$ClusterNum))
#304




#get the frequency of clusters
df.Freq<-data.frame(table(df.clusters$ClusterNum))

table(df.Freq$Freq)
#@amnahsiddiqa Please do it
#@amnahsiddiqa dont do it
# 1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  27  28  30  34  52 
# 132  54  26  12  17  13   6   4   4   3   5   2   6   1   2   1   2   2   1   1   2   2   1   2   1   1   1 



df.clusters<-cutree(data.tree , h = 0.4) %>%data.frame()
df.clusters$GeneSet<-row.names(df.clusters)
colnames(df.clusters)<-c("ClusterNum", "GeneSet")


#how many unique clusters are in there
length(unique(df.clusters$ClusterNum))
#304




#get the frequency of clusters
df.Freq<-data.frame(table(df.clusters$ClusterNum))

table(df.Freq$Freq)



df.clusters<-cutree(data.tree , h = 0.3) %>%data.frame()
df.clusters$GeneSet<-row.names(df.clusters)
colnames(df.clusters)<-c("ClusterNum", "GeneSet")


#how many unique clusters are in there
length(unique(df.clusters$ClusterNum))
#304




#get the frequency of clusters
df.Freq<-data.frame(table(df.clusters$ClusterNum))

table(df.Freq$Freq)








