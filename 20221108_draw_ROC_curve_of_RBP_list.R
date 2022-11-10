# this script is to draw ROC curve by each RBP research which listed in RBP2GO
# 2022/11/08 made

# activate package for drawing ROC curve
library(Epi)

# set function to find row that contains specify RBP name
find.row <-function(x,y,a){
  for (i in 1:length(x)) {
    g <-grep(x[i],y[,a])
    df <-y[g,]
  if(i==1){
    data <-df
  }else{
    data <-rbind(data,df)
  }
  }
  return(data)
}

# import RBP2GO list
# this talbe is located at "\\fsw-q02\okamura-lab\Files_related_to_M1_Projects\Hirota\RBP2GO_data"
setwd("C:/Rdata/RBP2GO_data")
RBP2GO <-read.table("RBP2GO_no_blank.txt",sep="\t",header = T,stringsAsFactors = F,quote = "")
RBP2GO <-RBP2GO[,c(1,2,11,15:57)]

# import table of RBP correspond to gene name in CCLE
# this talbe is located at "Dropbox/Okamura Lab share folder/Hirota/results_and_matterials/RBP2GO's_RBPs_contained_to_CCLE"
setwd("C:/Rdata/RBP2GO's_RBPs_contained_to_CCLE")
RBP.list <-read.table("all_RBP2GO_RBPs_which_merged_with_CCLE.txt",sep="\t",header = T,stringsAsFactors = F)

# import result of correlation analysis
# this result is located at "Dropbox/Okamura Lab share folder/Hirota/results_and_matterials/20220616_ROC_curve_for_determing_cutoff_of_cell_number"
setwd("C:/Rdata/ROC_curve_for_cutoff")
result <-read.table("combine_results_of_correlation_between_residual_and_RBP_exp.txt",sep="\t",header = T,stringsAsFactors = F)
result[,1] <-paste0(result[,1],"_",result[,2],"_",result[,3],"_vs_",result[,4])
r <-unique(result[,4])

# import list of physical interactions
# this list is located at "https://github.com/Ryosuke-Hirota/20221017_ROC_curve_with_list_of_functional_or_physical_interactions"
setwd("C:/Rdata/20221017_ROC_curve_for_cutoff_with_functional_interactions")
phy.list <-read.table("list_of_treiber_physical_interaction_between_RBP_and_miRNA.txt",sep="\t",header = T,stringsAsFactors = F)
phy.list <-phy.list[phy.list[,5]>=3,]
phy.list[,1] <-paste0(phy.list[,1],"_",phy.list[,2],"_",phy.list[,3],"_vs_",phy.list[,4])

# create new directory
setwd("C:/Rdata")
dir.create("20221107_ROC_curve_of_RBP_list")
setwd("C:/Rdata/20221107_ROC_curve_of_RBP_list")

# list number of rowsthat contain specify RBP name 
rn <-as.list(matrix(nrow =length(r),ncol = 1))

for (i in 1:length(r)) {
  data <-result[result[,4]==r[i],]
  n <-as.list(as.numeric(row.names(data)))
  rn[i] <-list(n)
}
names(rn) <-r

# annotate true combination to result
p <-match(phy.list[,1],result[,1])
p <-na.omit(p)

result[p,10] <-1
result[-p,10] <-0
colnames(result)[10] <-"physical_interaction"

# draw ROC curve
for (i in 4:46) {
  # extract RBP by each research
  df <-RBP2GO[,c(1:3,i)]
  df <-df[df[,4]=="X",]
  
  # extract RBP name in CCLE
  RBP <-unique(df[,2])
  RBP.df <-find.row(RBP,RBP.list,3)
  CCLE.RBP <-unique(RBP.df[,2])
  
  # extract number of rows from list and remove unneccesary rows
  r <-as.numeric(unlist(rn[CCLE.RBP]))
  data <-result[r,]
  
  # draw ROC curve
  pdf(paste0("ROC_curve_of_",colnames(RBP2GO)[i],".pdf"))
  ROC(test=data$number_of_cell, stat=data$physical_interaction, plot="ROC")
  dev.off()
  
  # output extracted result
  write.table(data,paste0(colnames(RBP2GO)[i],".txt"),sep="\t",row.names = F,quote = F)
  }
