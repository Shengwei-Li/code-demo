

# GSVA + hierarchical clustering + heatmap -------------------------------------------------
rm(list = ls())
library(GSVA)
library(estimate)
library(pheatmap)

# Read in immune gene sets
genesets = read.csv("./toy data/immune signatures.csv", stringsAsFactors = FALSE,header = FALSE)
genesets = as.data.frame(t(genesets))
geneSet=c()  

# Remove the blank spaces from each row and turn the dataframe into a list
for(i in 1:dim(genesets)[1]){
  if(genesets[i,dim(genesets)[2]] == ""){geneset = list(as.character(genesets[i,2:length(genesets[i,])])[-which(genesets[i,2:length(genesets[i,])] == "")])}
  else geneset = list(as.character(genesets[i,2:length(genesets[i,])]))
  names(geneset) = genesets[i,1]
  geneSet = c(geneSet,geneset)
}


# Quantification of the activity of the immune gene set for patients in each data set using the ssgsea algorithm
project = c("TCGA","METABRIC","GSE24450","GSE11121","GSE2034")

for(i in 1:length(project)){
  rna_seq = read.csv(paste0("./toy data/",project[i],"_expr.csv"),stringsAsFactors = F,header = T,check.names = F,row.names = 1)
  rna_seq = as.matrix(rna_seq)
  res = gsva(rna_seq,geneSet,method="ssgsea",ssgsea.norm = TRUE,verbose = TRUE)
  colnames(res) = gsub(colnames(res),pattern=".",replacement="-",fixed = TRUE)
  res = as.data.frame(res)
  write.csv(res, paste0(project[i],"_ssgsea.csv"))
}

# Hierarchical clustering and heatmap visualisation

for(i in 1:length(project)){
  ssgsea = read.csv(paste0(project[i],"_ssgsea.csv"),header=F,encoding="UTF-8")
  rownames(ssgsea) = as.character(ssgsea[,1])
  ssgsea = ssgsea[,-1]
  colnames(ssgsea) = t(ssgsea[1,])
  ssgsea = ssgsea[-1,]
  ssgsea = as.matrix(ssgsea)
  
 # Converting expression values to numeric
  ssgsea = matrix(as.numeric(ssgsea),nrow=nrow(ssgsea),ncol=ncol(ssgsea),dimnames=list(rownames(ssgsea),colnames(ssgsea)))
  
 # Normalisation
  df=scale(t(ssgsea)) 
  d=dist(df,method = "euclidean")   # dist() - Calculate the Euclidean distance between samples
  sample.hc = hclust(d,method="ward.D2")
  sample.id <- cutree(sample.hc,3)   #k=3
 
  sample.id = as.data.frame(sample.id)
  sample.id = data.frame(X = row.names(sample.id),x = sample.id[,1])
  
  sample.id[,2] = paste("cluster",sample.id[,2],sep = "")
  anno_col = data.frame(cluster=factor(sample.id[,2]))
  rownames(anno_col)=as.character(sample.id[,1])
  ann_colors = list(cluster = c(cluster1 = "#80B1D3", cluster2="#FDB462",cluster3="#FB8072"))
  pdf(paste0(project[i],"_heatmap.pdf"),width=8,height=15)
  heatmap = pheatmap(ssgsea,scale = 'row',cellheight = 12,show_colnames = FALSE,color=colorRampPalette(c("blue2", "white", "red"))(20),legend=F,
                   clustering_distance_cols = "euclidean",cluster_rows=FALSE,annotation_col=anno_col, annotation_colors = ann_colors,
                   clustering_method = "ward.D2",cutree_cols=3)
  dev.off()
}

# Significance statistics -------------------------------------------------------------------

####################   Mann–Whitney U test (the data are non-normally distributed), Student’s t test (the data are normally distributed)

rm(list=ls())

project = c("METABRIC","TCGA","GSE24450","GSE2034","GSE11121")

compare_immune_activity = function(project,x,y){
  # Parameter x and parameter y are the results of hierarchical clustering
  if (missing(project) || missing(x) || missing(y) || !is.numeric(c(x,y)))
    stop("'parameter' is missing or incorrect")
  
  result = c()
  
  for(i in 1:length(project)){
    ki = read.table(paste0("./toy data/",project[i],"_scores.txt"),sep = "\t")
    rownames(ki) = gsub(rownames(ki),pattern = ".",replacement = "-",fixed = TRUE)
    cluster = read.csv(paste0(project[i],"_cluster.csv"))
    cluster[,1] = gsub(cluster[,1],pattern = ".",replacement = "-",fixed = TRUE)
    cordata = merge(cluster,ki,by.x="X",by.y="row.names")
    pv = c()
    for(f in 3:dim(cordata)[2]){
      #wilcox
      pvG = wilcox.test(as.numeric(cordata[,f])[which(cordata[,2] == x)],as.numeric(cordata[,f])[which(cordata[,2] == y)],alternative = "greater")$p.value
      pvL = wilcox.test(as.numeric(cordata[,f])[which(cordata[,2] == x)],as.numeric(cordata[,f])[which(cordata[,2] == y)],alternative = "less")$p.value
      pv2 = cbind(project[i],colnames(cordata)[f],pvG,pvL)
      pv = rbind(pv,pv2)
      rownames(pv) = NULL
    }
    ord = order(as.numeric(pv[,"pvL"]),decreasing=FALSE)
    pvalue1 = pv[ord,]
    fdr = rep(1,dim(pvalue1)[1]); for(n in 1:dim(pvalue1)[1]) fdr[n]<-as.numeric(pvalue1[n,4])*dim(pvalue1)[1]/n
    pv2 = cbind(pvalue1,fdr)
    colnames(pv2) = c("project","ESTIMATE","pvG","pvL","fdr")
    result=rbind(result,pv2)
  }
  result = as.data.frame(result)
  write.csv(result,paste0(x,"-",y,"_wilcox.csv"),row.names = FALSE)
}

### wilcox cluster1 vs cluster2
compare_immune_activity(project,1,2)


### wilcox cluster1 vs cluster3
compare_immune_activity(project,1,3)

### wilcox cluster2 vs cluster3
compare_immune_activity(project,2,3)


