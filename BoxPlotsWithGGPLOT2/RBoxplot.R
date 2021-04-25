
'
Amnah Siddiqa: Updated April,22,2021
Boxplots of Rank normalized Rsquared data for non Redundant clusters

'



library(ggplot2)



# Data Wrangling ----------------------------------------------------------

#write the file for direct 
df.QN.long<-read.delim("/Users/siddia/Google Drive/1-FluVaccinesStudyManuscript/Figure2/Figure2b/Input/PointFive_DChremoved_max_long.tsv")

#verify
length(unique(df.QN.long$GenesSet))
#86 Clusters 

#for getting the correct order of boxplots
load("/Users/siddia/Google Drive/1-FluVaccinesStudyManuscript/Figure2/Figure2b/Input/vec.rda")

df.QN.long$GenesSet<- factor(df.QN.long$GenesSet, levels=vec.order$GenesSet)




# Plot Visualization ------------------------------------------------------


#pdf(file="/Users/siddia/Google Drive/1-FluVaccinesStudyManuscript/Figure2/Figure2b/Output/Figure2b.pdf",
#    width=4.6,
#    height=8.2)


# Box plot
ggplot(df.QN.long, aes(x=Rsq, y=GenesSet))+
  geom_boxplot(fill="grey", coef=1e30, alpha=0.5)+                                        # boxplot fill color
  theme_set(theme_bw(base_size  = 6, base_family = "Times"))+                             # base font family and size
  theme(text=element_text(color = "black"),                                               # only axis titles font color
        #axis.title.y =element_text(color = "red"),                                       # for only y or x axis title 
        axis.text = element_text(color = "black"))+                                       # axis  font color
  ylab("Genes Sets")+                                                                     # y-axis labels 
  xlab (expression(paste("R"^"2")))+                                                      # x-axis label
  theme(axis.title = element_text(size = 8))                                              # Axis titles font size


#dev.off()


#saving a plot Command 
