
library(GEOquery)
library(dplyr)
library(VennDiagram)
library(metaMA)
setwd('~/geo_analysis')

############################## Load data (from GSE) ##############################
data1_our = getGEO('GSE83452')
#data2_our = getGEO('GSE126848')
data3_our = getGEO('GSE24807')

data1_r = read.table('/data/home/mickey24/signatureSearch/GSE83452_nodup.entrezid.txt', header=TRUE, row.names=1, check.names=FALSE)
data1_r = data.frame(data1_r)
data1_eset = data1_r %>% dplyr::select(!contains("entrezgene"))

#data2_r = read.table('/data/home/mickey24/signatureSearch/GSE126848_nodup.entrezid.txt', header=TRUE, row.names=1, check.names=FALSE)
#data2_r = data.frame(data2_r)
#data2_eset = data2_r %>% dplyr::select(!contains("entrezgene"))

data3_r = read.table('/data/home/mickey24/signatureSearch/GSE24807_nodup.entrezid.txt', header=TRUE, row.names=1, check.names=FALSE)
data3_r = data.frame(data3_r)
data3_eset = data3_r %>% dplyr::select(!contains("entrezgene"))

############################## Box plot for gene expression ##############################
par(mfrow=c(2,2))
boxplot(data1_eset,main="GSE83452",outline=FALSE)
#boxplot(data2_eset,main="GSE126848",outline=FALSE)
boxplot(data3_eset,main="GSE24807",outline=FALSE)

#data2_eset_log<-log2(data2_eset)
#boxplot(data2_eset_log,main="GSE126848_log2",outline=FALSE)

#data3_eset_log<-log2(data3_eset+0.001)
#boxplot(data3_eset_log,main="GSE24807_log2",outline=FALSE)

############################## Get condition (case/control) ##############################
# Condition
c1_our = as.numeric( pData(data1_our[[1]])["source_name_ch1"]=='NASH liver baseline' ) # NASH liver follow-up
#c2_our = as.numeric( pData(data2_our[[1]])["disease:ch1"]=='NASH' ) # NASH liver follow-up
c3_our = as.numeric(pData(data3_our[[1]])["disease state:ch1"]=='non-alcoholic steatohepatitis (NASH)')
#classes_our=list(c1_our,c2_our,c3_our)
classes_our=list(c1_our,c3_our)

############################## summation data for gene expression ##############################
#esets_our = list(data1_eset, data2_eset, data3_eset)
#esets_our = list(data1_eset, data2_eset_log, data3_eset)
esets_our = list(data1_eset, data3_eset)


############################## Load data (unigene) ##############################
data1_u = read.table('/data/home/mickey24/signatureSearch/GSE83452.unigene.csv', sep=',', header=TRUE)
result1 = data1_u %>%
  group_by(unigene) %>% summarize(entrezgene = paste(sort(unique(entrezgene)),collapse=", "))
uni_dict1 = list()
for(i in 1:nrow(result1) ) {
  key <- result1[i,1]$unigene
  value <- unlist(strsplit(result1[i,2]$entrezgene, ","))
  uni_dict1[[ key  ]] =  value
}
data2_u = read.table('/data/home/mickey24/signatureSearch/GSE126848.unigene.csv', sep=',', header=TRUE)
result2 = data2_u %>%
  group_by(unigene) %>% summarize(entrezgene = paste(sort(unique(entrezgene)),collapse=", "))
uni_dict2 = list()
for(i in 1:nrow(result2) ) {
  key <- result2[i,1]$unigene
  value <- unlist(strsplit(result2[i,2]$entrezgene, ","))
  uni_dict2[[ key  ]] =  value
}
data3_u = read.table('/data/home/mickey24/signatureSearch/GSE24807.unigene.csv', sep=',', header=TRUE)
result3 = data3_u %>%
  group_by(unigene) %>% summarize(entrezgene = paste(sort(unique(entrezgene)),collapse=", "))
uni_dict3 = list()
for(i in 1:nrow(result3) ) {
  key <- result3[i,1]$unigene
  value <- unlist(strsplit(result3[i,2]$entrezgene, ","))
  uni_dict3[[ key  ]] =  value
}

#conv_unigene_our = list(uni_dict1, uni_dict2, uni_dict3)
conv_unigene_our = list(uni_dict1, uni_dict3)

############################## Run metaMA ##############################
## run MA
res_our=pvalcombination(esets=esets_our,classes=classes_our)
length(res_our$Meta)
Hs.Meta_our=rownames(esets_our[[1]])[res_our$Meta]

write.table(res_our$Meta, "runMeta_result.csv", sep='\t', col.names=NA)

############################## Venn diagram ##############################
#venn.plot<-venn.diagram(x = list(study1=res_our$study1, study2=res_our$study2, study3=res_our$study3,
venn.plot<-venn.diagram(x = list(study1=res_our$study1, study2=res_our$study2,
                                 meta=res_our$Meta),
                        filename = NULL, col = "black",
                        fill = c("blue", "red","green"),
                       # fill = c("blue", "red", "purple","green"),
                        margin=0.05, alpha = 0.6)

jpeg("venn_jpeg.jpg")
grid.draw(venn.plot)
dev.off()
