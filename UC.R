setwd('D:\\UCPF')
library(GDCRNATools)
library(tidyverse)
library(caret)
library(ggplot2)
library(dplyr)
library(factoextra)
options(stringsAsFactors = F)
flm=read.table('GSE75214.lnc_mrna.exp.1.txt',row.names = 1,header=T,stringsAsFactors = F)
lprob=read.table('GSE75214.prob.1.lnc.1.txt',header=F)
head(lprob)
probl=lprob$V2
names(probl)=lprob$V1
exprl1=flm[rownames(flm)%in%names(probl),]
exprl1$gene=probl[rownames(exprl1)]
sal=na.omit(unique(exprl1$gene))
exprl2=data.frame()
#for (i in sal)
#{
#  dbs=na.omit(exprl1[exprl1$gene%in%c(i),])
#  dbx=dbs[,1:(ncol(dbs)-1)]
#  dbx1=t(data.frame(apply(dbx,2,median)))
#  rownames(dbx1)=i
#  exprl2=rbind(exprl2,dbx1)
#}
#exprl2
#write.table(exprl2,'GSE75214.lnc.expression.1.txt',row.names = T,col.names = T,sep = '\t')
exprl2=read.table('GSE75214.lnc.expression.1.txt',row.names = 1,header=T,sep='\t')

mprob=read.table('GSE75214.prob.1.pc.1.txt',header=F)
probm=mprob$V2
names(probm)=mprob$V1
exprm=flm[rownames(flm)%in%names(probm),]
exprm$gene=probm[rownames(exprm)]
sam=na.omit(unique(exprm$gene))
#exprm1=data.frame()
#for (i in sam)
#{
#  dbs=na.omit(exprm[exprm$gene%in%c(i),])
#  dbx=dbs[,1:(ncol(dbs)-1)]
#  dbx1=t(data.frame(apply(dbx,2,median)))
#  rownames(dbx1)=i
#  exprm1=rbind(exprm1,dbx1)
#}
#exprm1
#write.table(exprm1,'GSE75214.expression.1.txt',row.names = T,col.names = T,sep = '\t')
exprm1=read.table('GSE75214.expression.1.txt',row.names = 1,header=T,sep='\t')
fcl=read.table('GSE75214.lnc_mrna.exp.cli.1.txt',header=T)
slc=fcl$group
slc1=fcl$group1
names(slc)=fcl$sampleID
names(slc1)=fcl$sampleID
smc=slc


####plasm deg-lnc####
library(limma)
exprl2=exprl2[,names(slc)]
group=slc[colnames(exprl2)]
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
contrast.matrix<-makeContrasts("ulcerative_colitis-control",
                               levels = design)
contrast.matrix
ddlnc<-as.data.frame(t(exprl2))
res.pca <- prcomp(ddlnc, scale = TRUE)
fviz_eig(res.pca)
color=c("#00AFBB",  "#FC4E07")
fviz_pca_ind(res.pca,
             col.ind = factor(slc[rownames(ddlnc)]), # 颜色对应group信息
             palette = color,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "LncRNA Group",## Legend名称
             repel = TRUE
)

fit <- lmFit(exprl2,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEGlnc<-topTable(fit2, coef=1, n=Inf) %>% na.omit()

require(EnhancedVolcano)
EnhancedVolcano(DEGlnc,title="DEGlnc",
                lab = rownames(DEGlnc),
                x = "logFC",
                y = "P.Value",
                selectLab = rownames(DEGlnc)[1:5],
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,## pvalue阈值
                FCcutoff = 0.58,## FC cutoff
                xlim = c(-5,5),
                colAlpha = 1,
                legendPosition = "bottom",
                legendLabSize = 10,
                legendIconSize = 3.0)
write.table(DEGlnc,'GSE75214.DEG_lnc.txt',row.names = T,col.names = T,sep='\t')
nrow(DEGlnc[DEGlnc$P.Value<0.05&abs(DEGlnc$logFC)>0.58,])
DEGlnc1=DEGlnc[DEGlnc$P.Value<0.05&abs(DEGlnc$logFC)>0.58,]
write.table(DEGlnc1,'GSE75214.DEG_lnc.sig.txt',row.names = T,col.names = T,sep='\t')
index<-DEGlnc1 %>% as_tibble() %>% 
  mutate(genename=rownames(DEGlnc1)) %>% 
  dplyr::arrange(desc(logFC)) %>% 
  .$genename 

library(pheatmap)
exprl1s=exprl2[,names(group)]
index_matrix<-t(scale(t(exprl1s[index,])))##归一化
index_matrix[index_matrix>2]=2
index_matrix[index_matrix<(-2)]=-2
head(index_matrix)
anno=data.frame(group=group)
rownames(anno)=colnames(index_matrix)
pheatmap(index_matrix,
         show_colnames =F,
         show_rownames = T,
         cluster_cols = F, 
         annotation_col=anno)
###########DEG-mRNA###############
exprm1=exprm1[,names(smc)]
group=smc[colnames(exprm1)]
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
contrast.matrix<-makeContrasts("ulcerative_colitis-control",
                               levels = design)
contrast.matrix
ddmrna<-as.data.frame(t(exprm1))
res.pca <- prcomp(ddmrna, scale = TRUE)
fviz_eig(res.pca)
color=c("#00AFBB",  "#FC4E07")
fviz_pca_ind(res.pca,
             col.ind = factor(smc[rownames(ddmrna)]), # 颜色对应group信息
             palette = color,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "mRNA Group",## Legend名称
             repel = TRUE
)

fit <- lmFit(exprm1,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEGmRNA<-topTable(fit2, coef=1, n=Inf) %>% na.omit()

require(EnhancedVolcano)
EnhancedVolcano(DEGmRNA,title="DEGmRNA",
                lab = rownames(DEGmRNA),
                x = "logFC",
                y = "P.Value",
                selectLab = rownames(DEGmRNA)[1:5],
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,## pvalue阈值
                FCcutoff = 0.58,## FC cutoff
                xlim = c(-5,5),
                colAlpha = 1,
                legendPosition = "bottom",
                legendLabSize = 10,
                legendIconSize = 3.0)
write.table(DEGmRNA,'GSE75214.DEG_mRNA.txt',row.names = T,col.names = T,sep='\t')
nrow(DEGmRNA[DEGmRNA$P.Value<0.05&abs(DEGmRNA$logFC)>0.58,])
DEGmRNA1=DEGmRNA[DEGmRNA$P.Value<0.05&abs(DEGmRNA$logFC)>0.58,]
write.table(DEGmRNA1,'GSE75214.DEG_mRNA.sig.txt',row.names = T,col.names = T,sep='\t')
index<-DEGmRNA1 %>% as_tibble() %>% 
  mutate(genename=rownames(DEGmRNA1)) %>% 
  dplyr::arrange(desc(logFC)) %>% 
  .$genename 

library(pheatmap)
exprm1s=exprm1[,names(group)]
index_matrix<-t(scale(t(exprm1s[index,])))##归一化
index_matrix[index_matrix>2]=2
index_matrix[index_matrix<(-2)]=-2
head(index_matrix)
anno=data.frame(group=group)
rownames(anno)=colnames(index_matrix)
pheatmap(index_matrix,
         show_colnames =F,
         show_rownames = T,
         cluster_cols = F, 
         annotation_col=anno)
########################################plasm mrna-lncRna###################
dsln=exprl2[rownames(DEGlnc1),]
dsln1=as.data.frame(t(dsln))
dsmrna=exprm1[rownames(DEGmRNA1),]
dsmrna1=as.data.frame(t(dsmrna))
dlnm=data.frame()
for (i in 1:ncol(dsln1))
{
  for (j in 1:ncol(dsmrna1))
  {
    print(i)
    print(j)
    dtmp=data.frame()
    dtmp[1,1]=colnames(dsln1)[i]
    dtmp[1,2]=colnames(dsmrna1)[j]
    dtmp[1,3]=cor.test(as.numeric(dsln1[,i]),as.numeric(dsmrna1[,j]))$estimate
    dtmp[1,4]=cor.test(as.numeric(dsln1[,i]),as.numeric(dsmrna1[,j]))$p.value
    dlnm=rbind(dlnm,dtmp)
 }
}
nrow(dlnm)

colnames(dlnm)=c('LNCRNA','mRNA','R','pvalue')
head(dlnm)
write.table(dlnm,'correlation_DEGlncRNA-DEGmRNA.GSE75214.txt',row.names = F,col.names = T,sep='\t')
dlnm=read.table('correlation_DEGlncRNA-DEGmRNA.GSE75214.txt',header = T,sep = '\t')
dlnm$fdr=p.adjust(dlnm$pvalue,n=nrow(dlnm))
write.table(dlnm,'correlation_DEGlncRNA-DEGmRNA.complete.GSE75214.fdr.txt',col.names = T,sep='\t')
dlnm1=dlnm[dlnm$pvalue<0.05,]
write.table(dlnm1,'correlation_DEGlncRNA-DEGmRNA.complete.GSE75214.fdr.sig.txt',col.names = T,sep='\t')
dal=data.frame()
for (i in 1:nrow(DEGlnc1))
{
  dall=data.frame()
  dall[1,1]=rownames(DEGlnc1)[i]
  if(DEGlnc1[i,1]>0)
  {
    dall[1,2]='up'
  }else{
    dall[1,2]='down'
  }
  dall[1,3]='lnc'
  dal=rbind(dal,dall)
}
dal

dam=data.frame()
for (i in 1:nrow(DEGmRNA1))
{
  dall=data.frame()
  dall[1,1]=rownames(DEGmRNA1)[i]
  if(DEGmRNA1[i,1]>0)
  {
    dall[1,2]='up'
  }else{
    dall[1,2]='down'
  }
  dall[1,3]='mRNA'
  dam=rbind(dam,dall)
}
dam
dalm=rbind(dal,dam)
colnames(dalm)=c('node','class','group')
write.table(dalm,'GSE75214-node-lncrna-annotation.txt',col.names = T,sep = '\t')

#########################miRNA-mRNA##################


fcmi=read.table('GSE48959-GPL14613_series_matrix.miRNA.cli.1.txt',header=T,row.names = 1)
fcmm=read.table('GSE48959-GPL6244_series_matrix.mRNA.cli.1.txt',header=T,row.names = 1)
fcmi=fcmi[rownames(fcmi)%in%rownames(fcmm),]
fcmm=fcmm[rownames(fcmi),]
scsss=rownames(fcmi)
names(scsss)=fcmi$sampleID

exprxmi=read.table('GSE48959-GPL14613_series_matrix.miRNA.exp.1.txt',row.names = 1,header=T,sep='\t')
exprxmi=exprxmi[,names(scsss)]
colnames(exprxmi)=scsss[colnames(exprxmi)]

miprob=read.table('GSE48959-GPL14613_series_matrix.miRNA.prob.1.1.1.txt',header=F)
probmi=miprob$V1
names(probmi)=miprob$V4

exprxmi=exprxmi[rownames(exprxmi)%in%names(probmi),]
exprxmi$gene=probmi[rownames(exprxmi)]
sam=na.omit(unique(exprxmi$gene))
exprxmi1=data.frame()
for (i in sam)
{
  dbs=na.omit(exprxmi[exprxmi$gene%in%c(i),])
  dbx=dbs[,1:(ncol(dbs)-1)]
  dbx1=t(data.frame(apply(dbx,2,median)))
  rownames(dbx1)=i
  exprxmi1=rbind(exprxmi1,dbx1)
}
#exprxmi1
#write.table(exprxmi1,'GSE48959-GPL14613_series_matrix.miRNA.exp.1.txt',row.names = T,col.names = T,sep = '\t')
exprxmi1=read.table('GSE48959-GPL14613_series_matrix.miRNA.exp.1.txt',row.names = 1,header=T,sep='\t')
exprxmi1=exprxmi1[,names(scsss)]
colnames(exprxmi1)=scsss[colnames(exprxmi1)]

smic=fcmi$group
names(smic)=rownames(fcmi)

exprxmi1=exprxmi1[,names(smic)]
group=smic[colnames(exprxmi1)]
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
contrast.matrix<-makeContrasts("UC-control",
                               levels = design)
contrast.matrix
ddmirna<-as.data.frame(t(exprxmi1))
res.pca <- prcomp(ddmirna, scale = TRUE)
fviz_eig(res.pca)
color=c("#00AFBB",  "#FC4E07")
fviz_pca_ind(res.pca,
             col.ind = factor(smic[rownames(ddmirna)]), # 颜色对应group信息
             palette = color,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "miRNA Group",## Legend名称
             repel = TRUE
)

fit <- lmFit(exprxmi1,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEGmiRNA<-topTable(fit2, coef=1, n=Inf) %>% na.omit()

require(EnhancedVolcano)
EnhancedVolcano(DEGmiRNA,title="DEGmiRNA",
                lab = rownames(DEGmiRNA),
                x = "logFC",
                y = "P.Value",
                selectLab = rownames(DEGmiRNA)[1:5],
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,## pvalue阈值
                FCcutoff =0.58,## FC cutoff
                xlim = c(-5,5),
                colAlpha = 1,
                legendPosition = "bottom",
                legendLabSize = 10,
                legendIconSize = 3.0)
write.table(DEGmiRNA,'GSE48959-GPL14613_DEG_miRNA.txt',row.names = T,col.names = T,sep='\t')
DEGmiRNA1=DEGmiRNA[DEGmiRNA$P.Value<0.05&abs(DEGmiRNA$logFC)>0.58,]
write.table(DEGmiRNA1,'GSE48959-GPL14613_DEGmirna.def.genes.txt',row.names = T,col.names = T,sep='\t')
index<-DEGmiRNA1 %>% as_tibble() %>% 
  mutate(genename=rownames(DEGmiRNA1)) %>% 
  dplyr::arrange(desc(logFC)) %>% 
  .$genename 

library(pheatmap)
exprmi1s=exprxmi1[,names(sort(smic))]
index_matrix<-t(scale(t(exprmi1s[index,])))##归一化
index_matrix[index_matrix>1]=1
index_matrix[index_matrix<(-1)]=-1
head(index_matrix)
anno=data.frame(group=sort(group))
rownames(anno)=colnames(index_matrix)
pheatmap(index_matrix,
         show_colnames =F,
         show_rownames = T,
         cluster_cols = F, 
         annotation_col=anno)

#####deg_mRNA
exprmx=read.table('GSE48959-GPL6244_series_matrix.mRNA.exp.1.txt',row.names = 1,header=T,sep='\t')
mprob=read.table('GSE48959-GPL6244_series_matrix.prob.txt',header=F)
probm=mprob$V2
names(probm)=mprob$V1
exprmmx=exprmx[rownames(exprmx)%in%names(probm),]
exprmmx$gene=probm[rownames(exprmmx)]
sam=na.omit(unique(exprmmx$gene))
#exprmmx1=data.frame()
#for (i in sam)
#{
#  dbs=na.omit(exprmmx[exprmmx$gene%in%c(i),])
#  dbx=dbs[,1:(ncol(dbs)-1)]
#  dbx1=t(data.frame(apply(dbx,2,median)))
#  rownames(dbx1)=i
#  exprmmx1=rbind(exprmmx1,dbx1)
#}
#exprmmx1
#write.table(exprmmx1,'GSE48959-GPL6244.mrna.expression.1.txt',row.names = T,col.names = T,sep = '\t')
exprmmx1=read.table('GSE48959-GPL6244.mrna.expression.1.txt',row.names = 1,header=T,sep='\t')

scsssm=rownames(fcmm)
names(scsssm)=fcmm$sampleID

exprmmx1=exprmmx1[,names(scsssm)]
colnames(exprmmx1)=scsssm[colnames(exprmmx1)]

smmc=fcmm$group
names(smmc)=rownames(fcmm)
smmc1=fcmm$disease_activity_
names(smmc1)=rownames(fcmm)
exprmmx1=log2(exprmmx1[,names(smmc)])
group=smmc[colnames(exprmmx1)]
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
contrast.matrix<-makeContrasts("UC-control",
                               levels = design)
contrast.matrix
ddmrna<-as.data.frame(t(exprmmx1))
res.pca <- prcomp(ddmrna, scale = TRUE)
fviz_eig(res.pca)
color=c("#00AFBB",  "#FC4E07")
fviz_pca_ind(res.pca,
             col.ind = factor(smmc[rownames(ddmrna)]), # 颜色对应group信息
             palette = color,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "mRNA Group",## Legend名称
             repel = TRUE
)

fit <- lmFit(exprmmx1,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEGmmRNA<-topTable(fit2, coef=1, n=Inf) %>% na.omit()

require(EnhancedVolcano)
EnhancedVolcano(DEGmmRNA,title="DEGmRNA",
                lab = rownames(DEGmmRNA),
                x = "logFC",
                y = "P.Value",
                selectLab = rownames(DEGmmRNA)[1:5],
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,## pvalue阈值
                FCcutoff = 0.58,## FC cutoff
                xlim = c(-5,5),
                colAlpha = 1,
                legendPosition = "bottom",
                legendLabSize = 10,
                legendIconSize = 3.0)
write.table(DEGmmRNA,'GSE48959-GPL6244.DEG_mRNA.txt',row.names = T,col.names = T,sep='\t')
nrow(DEGmmRNA[DEGmmRNA$P.Value<0.05&abs(DEGmmRNA$logFC)>0.58,])
DEGmmRNA1=DEGmmRNA[DEGmmRNA$P.Value<0.05&abs(DEGmmRNA$logFC)>0.58,]
write.table(DEGmmRNA1,'GSE48959-GPL6244.DEG_mRNA.sig.txt',row.names = T,col.names = T,sep='\t')
index<-DEGmmRNA1 %>% as_tibble() %>% 
  mutate(genename=rownames(DEGmmRNA1)) %>% 
  dplyr::arrange(desc(logFC)) %>% 
  .$genename 

library(pheatmap)
exprm1s=exprmmx1[,names(group)]
index_matrix<-t(scale(t(exprm1s[index,])))##归一化
index_matrix[index_matrix>2]=2
index_matrix[index_matrix<(-2)]=-2
head(index_matrix)
anno=data.frame(group=group)
rownames(anno)=colnames(index_matrix)
pheatmap(index_matrix,
         show_colnames =F,
         show_rownames = T,
         cluster_cols = F, 
         annotation_col=anno)
#########################miRNA-mRNA
dsmi=exprxmi1[rownames(DEGmiRNA1),]
dsmi1=as.data.frame(t(dsmi))
dsmrna=exprmmx1[rownames(DEGmmRNA1),]
dsmrna1=as.data.frame(t(dsmrna))
dmim=data.frame()
for (i in 1:ncol(dsmi1))
{
  for (j in 1:ncol(dsmrna1))
  {
    print(i)
    print(j)
    dtmp=data.frame()
    dtmp[1,1]=colnames(dsmi1)[i]
    dtmp[1,2]=colnames(dsmrna1)[j]
    dtmp[1,3]=cor.test(as.numeric(dsmi1[,i]),as.numeric(dsmrna1[,j]))$estimate
    dtmp[1,4]=cor.test(as.numeric(dsmi1[,i]),as.numeric(dsmrna1[,j]))$p.value
    dmim=rbind(dmim,dtmp)
  }
}
nrow(dmim)

colnames(dmim)=c('MiRNA','mRNA','R','pvalue')
head(dmim)
write.table(dmim,'correlation_DEGmiRNA-DEGmRNA.GSE48959.txt',row.names = F,col.names = T,sep='\t')
dmim=read.table('correlation_DEGmiRNA-DEGmRNA.GSE48959.txt',header = T,sep = '\t')
dmim$fdr=p.adjust(dmim$pvalue,n=nrow(dmim))
write.table(dmim,'correlation_DEGlncRNA-DEGmRNA.complete.GSE48959.fdr.txt',col.names = T,sep='\t')
dmim1=dmim[dmim$pvalue<0.05&dmim$R<0,]
write.table(dmim1,'correlation_DEGlncRNA-DEGmRNA.complete.GSE48959.fdr.sig.txt',col.names = T,sep='\t')
rownames(dmim1)=paste(dmim1$MiRNA,dmim1$mRNA,sep='|')
f23=read.table('3database_mir_target_genes.2_db_gene_mir_predicted.txt',header = F)
rownames(f23)=paste(f23$V1,f23$V2,sep='|')
dmim1=dmim1[rownames(dmim1)%in%rownames(f23),]
write.table(dmim1,'correlation_DEGlncRNA-DEGmRNA.complete.GSE48959.fdr.sig.1.txt',col.names = T,sep='\t')
unique(dmim1$MiRNA)
unique(dmim1$mRNA)

dami=data.frame()
for (i in 1:nrow(DEGmiRNA1))
{
  dall=data.frame()
  dall[1,1]=rownames(DEGmiRNA1)[i]
  if(DEGmiRNA1[i,1]>0)
  {
    dall[1,2]='up'
  }else{
    dall[1,2]='down'
  }
  dall[1,3]='miRNA'
  dami=rbind(dami,dall)
}
dami

dam=data.frame()
for (i in 1:nrow(DEGmRNA1))
{
  dall=data.frame()
  dall[1,1]=rownames(DEGmRNA1)[i]
  if(DEGmRNA1[i,1]>0)
  {
    dall[1,2]='up'
  }else{
    dall[1,2]='down'
  }
  dall[1,3]='mRNA'
  dam=rbind(dam,dall)
}
dam
damim=rbind(dami,dam)
colnames(damim)=c('node','class','group')
write.table(damim,'GSE48959-node-miRna-annotation.txt',col.names = T,sep = '\t')

###############################miRNA-mrna-lncRNA####################
dmim11=dmim1[dmim1$mRNA%in%dlnm1$mRNA,]
colnames(dmim11)[1:2]=c('node1','node2')
dlnm11=dlnm1[dlnm1$mRNA%in%dmim1$mRNA,]
colnames(dlnm11)[1:2]=c('node1','node2')
dmimln=rbind(dmim11,dlnm11)
dmimln1=data.frame()
pxx=vector()
for (i in  1:nrow(dmimln))
{
  if(dmimln[i,3]>0)
  {
    pxx[i]=1
  }else{
      pxx[i]=0
    }
}
dmimln$positive=pxx
dmimln$R1=abs(dmimln$R)
write.table(dmimln,'UC-lnc-mrna-mirna.txt',row.names = F,col.names = T,sep='\t')
###################################mRNA ROC
#fox=read.table('oxsidtive.all.txt',header=F)
#ox=colnames(exprm1ss)[colnames(exprm1ss)%in%fox$V1]

exprm1ss=data.frame(t(exprm1[unique(dmimln$node2),]))
exprmmx1ss=data.frame(t(exprmmx1[unique(dmimln$node2),]))

#exprm1ss=data.frame(t(exprm1[ox,]))
#exprmmx1ss=data.frame(t(exprmmx1[ox,]))
exprm1ss$gr=smc[rownames(exprm1ss)]
exprm1ss$gr[exprm1ss$gr=='control']=0
exprm1ss$gr[exprm1ss$gr=='ulcerative_colitis']=1
exprm1ss$gr=factor(exprm1ss$gr)
exprmmx1ss$gr=smmc[rownames(exprmmx1ss)]
exprmmx1ss$gr[exprmmx1ss$gr=='control']=0
exprmmx1ss$gr[exprmmx1ss$gr=='UC']=1
exprmmx1ss$gr=factor(exprmmx1ss$gr)

library(lars)
data4=exprm1ss
x = as.matrix(data4[,1:(ncol(data4)-1)]) 
y = as.matrix(as.numeric(data4[,ncol(data4)]))
rownames(y)=rownames(x)
colnames(y)=c("gr")
lar1 <-lars(x,y,type = "lasso")
summary(lar1)
larcp=as.vector(lar1$Cp)
for (i in 1:length(larcp))
{
  if (larcp[i] == min(larcp))
  {
    n=i
  }
}#n
beta <- lar1$beta[n,]
beta
mgene <- names(beta[beta!=0])
beta=beta[mgene]

data4=exprmmx1ss
x = as.matrix(data4[,1:(ncol(data4)-1)]) 
y = as.matrix(as.numeric(data4[,ncol(data4)]))
rownames(y)=rownames(x)
colnames(y)=c("gr")
lar1 <-lars(x,y,type = "lasso")
summary(lar1)
larcp=as.vector(lar1$Cp)
for (i in 1:length(larcp))
{
  if (larcp[i] == min(larcp))
  {
    n=i
  }
}#n
beta <- lar1$beta[n,]
beta
mgene <- names(beta[beta!=0])
beta=beta[mgene]


m20=exprm1ss[exprm1ss$gr==0,]
m21=exprm1ss[exprm1ss$gr==1,]

#drx2=vector()
dcb=data.frame()
p=0.7
nxb=1
indm20 <- sample(2, nrow(m20), replace = TRUE, prob = c(p, 1-p))
if (length(indm20[indm20==1])<3||length(indm20[indm20==2])<3)
{
  while (length(indm20[indm20==1])<3||length(indm20[indm20==2])<3)
  {
    indm20 <- sample(2, nrow(m20), replace = TRUE, prob = c(p, 1-p))
    if (length(indm20[indm20==1])>=3&&length(indm20[indm20==2])>=3)
    {
      trainm20 <- data.frame(m20[indm20==1,])
      testm20 <- as.data.frame(m20[indm20==2,])
      break
    }
    
  }
}else
{
  trainm20 <- data.frame(m20[indm20==1,])
  testm20 <- as.data.frame(m20[indm20==2,])
}

indm21 <- sample(2, nrow(m21), replace = TRUE, prob = c(p, 1-p))
if (length(indm21[indm21==1])<3||length(indm21[indm21==2])<3)
{
  while (length(indm21[indm21==1])<3||length(indm21[indm21==2])<3)
  {
    indm21 <- sample(2, nrow(m21), replace = TRUE, prob = c(p, 1-p))
    if (length(indm21[indm21==1])>=3&&length(indm21[indm21==2])>=3)
    {
      trainm21 <- data.frame(m21[indm21==1,])
      testm21 <- as.data.frame(m21[indm21==2,])
      break
    }
    
  }
}else
{
  trainm21 <- data.frame(m21[indm21==1,])
  testm21 <- as.data.frame(m21[indm21==2,])
}


trainm22 <- rbind(trainm20,trainm21)
#trainm2 <- data5


testm22 <- rbind(testm20,testm21)
#testm2 <- data5


library('randomForest')
library('pROC')
rfm2 <- randomForest(gr~., data=trainm22, ntree = 500,mtry=3,importance=TRUE,proximity=TRUE)
p2m2 <- predict(rfm2, testm22)
#obs_p_ranm2 = data.frame(prob=p2m2,obs=testm2$cluster)
ran_rocm2 <- roc(testm22$gr,as.numeric(p2m2))
#drx2[nxb]=as.numeric(ran_rocm2$auc)
drx2=as.numeric(ran_rocm2$auc)
names(drx2)='AUC'
drf2=importance(rfm2) 
dip1x2=drf2[,(ncol(drf2)-1)]
dcu=dip1x2
dcu=data.frame(feature=names(rev(sort(dcu))),MeanDecreaseAccuracy=rev(sort(dcu)))
dcu1=dcu

ggroc(ran_rocm2)
#ROC只用于2分类，对多分类不适用
plot(ran_rocm2, print.auc=TRUE, auc.polygon=TRUE, 
     grid=c(0.1, 0.2),grid.col=c("green", "red"), 
     max.auc.polygon=TRUE,auc.polygon.col="skyblue", 
     print.thres=TRUE,
     main='随机森林模型ROC曲线,mtry=3,ntree=500')

dippx2=as.data.frame(c(dip1x2,drx2))
colnames(dippx2)=paste("round",as.character(nxb),sep='_')
dipx2=dippx2
for (i in 2:1000)
{
  nxb=i

  indm20 <- sample(2, nrow(m20), replace = TRUE, prob = c(p, 1-p))
  if (length(indm20[indm20==1])<3||length(indm20[indm20==2])<3)
  {
    while (length(indm20[indm20==1])<3||length(indm20[indm20==2])<3)
    {
      indm20 <- sample(2, nrow(m20), replace = TRUE, prob = c(p, 1-p))
      if (length(indm20[indm20==1])>=3&&length(indm20[indm20==2])>=3)
      {
        trainm20 <- data.frame(m20[indm20==1,])
        testm20 <- as.data.frame(m20[indm20==2,])
        break
      }
      
    }
  }else
  {
    trainm20 <- data.frame(m20[indm20==1,])
    testm20 <- as.data.frame(m20[indm20==2,])
  }
  
  indm21 <- sample(2, nrow(m21), replace = TRUE, prob = c(p, 1-p))
  if (length(indm21[indm21==1])<3||length(indm21[indm21==2])<3)
  {
    while (length(indm21[indm21==1])<3||length(indm21[indm21==2])<3)
    {
      indm21 <- sample(2, nrow(m21), replace = TRUE, prob = c(p, 1-p))
      if (length(indm21[indm21==1])>=3&&length(indm21[indm21==2])>=3)
      {
        trainm21 <- data.frame(m21[indm21==1,])
        testm21 <- as.data.frame(m21[indm21==2,])
        break
      }
      
    }
  }else
  {
    trainm21 <- data.frame(m21[indm21==1,])
    testm21 <- as.data.frame(m21[indm21==2,])
  }
  
  
  trainm22 <- rbind(trainm20,trainm21)
  #trainm2 <- data5
  
  
  testm22 <- rbind(testm20,testm21)
  #testm2 <- data5
  
  
  library('randomForest')
  library('pROC')
  rfm2 <- randomForest(gr~., data=trainm22, ntree = 500,mtry=3,importance=TRUE,proximity=TRUE)
  p2m2 <- predict(rfm2, testm22)
  #obs_p_ranm2 = data.frame(prob=p2m2,obs=testm2$cluster)
  ran_rocm2 <- roc(testm22$gr,as.numeric(p2m2))
  #drx2[nxb]=as.numeric(ran_rocm2$auc)
  drx2=as.numeric(ran_rocm2$auc)
  names(drx2)='AUC'
  drf2=importance(rfm2)
  dip1x2=drf2[,(ncol(drf2)-1)]
  dippx2=as.data.frame(c(dip1x2,drx2))
  colnames(dippx2)=paste("round",as.character(nxb),sep='_')
  dipx2=cbind(dipx2,dippx2)
  
}
dipx2$prob=rep(p,nrow(dipx2))
dcb=rbind(dcb,dipx2)
f4=data.frame(t(dcb))
da4=data.frame(AUC=f4[,c('AUC')])

p<-ggplot(data = da4) + geom_histogram(aes(x = AUC, y = ..density..),fill = "red", alpha = 0.2) +
  geom_density(aes(x= AUC, y = ..density..), col = "red")
p 
f5=f4
d=data.frame()
for (i in 1:(ncol(f5)-1))
{
  d1=data.frame(importance=f5[,i],factors=rep(colnames(f5)[i],nrow(f5)))
  d=rbind(d,d1)
}
d
sssi=names(sort(tapply(d$importance,d$factors,median),decreasing = T))
library(ggsignif)
library(ggpubr)
library(ggplot2)
p=ggplot(d,aes(x=factor(factors,levels = sssi),y=importance,fill=factors),color='#C3241F')+geom_boxplot()+
  geom_hline(yintercept = 0, color='red',linetype = 2)+
  #scale_colour_manual(values = rev(c('#FF6347','#1E90FF')))+
  #scale_fill_manual(values = rev(c('#C3241F','#195FAC')))+
  #stat_compare_means(aes(group = group),method='wilcox.test'
  #                   ,label="p.signif"
  #)+
  #ylim(c(0,0.2))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()
        ,panel.background = element_blank(),legend.position = 'none',
        axis.text.x = element_text( hjust = 1,vjust= 1,angle =30,size = 15),
        axis.line = element_line(colour = "black"))
p 



m20=exprmmx1ss[exprmmx1ss$gr==0,]
m21=exprmmx1ss[exprmmx1ss$gr==1,]
#drx2=vector()
dcb=data.frame()
p=0.7
nxb=1
indm20 <- sample(2, nrow(m20), replace = TRUE, prob = c(p, 1-p))
if (length(indm20[indm20==1])<3||length(indm20[indm20==2])<3)
{
  while (length(indm20[indm20==1])<3||length(indm20[indm20==2])<3)
  {
    indm20 <- sample(2, nrow(m20), replace = TRUE, prob = c(p, 1-p))
    if (length(indm20[indm20==1])>=3&&length(indm20[indm20==2])>=3)
    {
      trainm20 <- data.frame(m20[indm20==1,])
      testm20 <- as.data.frame(m20[indm20==2,])
      break
    }
    
  }
}else
{
  trainm20 <- data.frame(m20[indm20==1,])
  testm20 <- as.data.frame(m20[indm20==2,])
}

indm21 <- sample(2, nrow(m21), replace = TRUE, prob = c(p, 1-p))
if (length(indm21[indm21==1])<3||length(indm21[indm21==2])<3)
{
  while (length(indm21[indm21==1])<3||length(indm21[indm21==2])<3)
  {
    indm21 <- sample(2, nrow(m21), replace = TRUE, prob = c(p, 1-p))
    if (length(indm21[indm21==1])>=3&&length(indm21[indm21==2])>=3)
    {
      trainm21 <- data.frame(m21[indm21==1,])
      testm21 <- as.data.frame(m21[indm21==2,])
      break
    }
    
  }
}else
{
  trainm21 <- data.frame(m21[indm21==1,])
  testm21 <- as.data.frame(m21[indm21==2,])
}


trainm22 <- rbind(trainm20,trainm21)
#trainm2 <- data5


testm22 <- rbind(testm20,testm21)
#testm2 <- data5


library('randomForest')
library('pROC')
rfm2 <- randomForest(gr~., data=trainm22, ntree = 500,mtry=3,importance=TRUE,proximity=TRUE)
p2m2 <- predict(rfm2, testm22)
#obs_p_ranm2 = data.frame(prob=p2m2,obs=testm2$cluster)
ran_rocm2 <- roc(testm22$gr,as.numeric(p2m2))
#drx2[nxb]=as.numeric(ran_rocm2$auc)
drx2=as.numeric(ran_rocm2$auc)
names(drx2)='AUC'
drf2=importance(rfm2) 
dip1x2=drf2[,(ncol(drf2)-1)]
dcu=dip1x2
dcu=data.frame(feature=names(rev(sort(dcu))),MeanDecreaseAccuracy=rev(sort(dcu)))
dcu1=dcu

ggroc(ran_rocm2)
#ROC只用于2分类，对多分类不适用
plot(ran_rocm2, print.auc=TRUE, auc.polygon=TRUE, 
     grid=c(0.1, 0.2),grid.col=c("green", "red"), 
     max.auc.polygon=TRUE,auc.polygon.col="skyblue", 
     print.thres=TRUE,
     main='随机森林模型ROC曲线,mtry=3,ntree=500')

dippx2=as.data.frame(c(dip1x2,drx2))
colnames(dippx2)=paste("round",as.character(nxb),sep='_')
dipx2=dippx2
for (i in 2:1000)
{
  nxb=i
  
  indm20 <- sample(2, nrow(m20), replace = TRUE, prob = c(p, 1-p))
  if (length(indm20[indm20==1])<3||length(indm20[indm20==2])<3)
  {
    while (length(indm20[indm20==1])<3||length(indm20[indm20==2])<3)
    {
      indm20 <- sample(2, nrow(m20), replace = TRUE, prob = c(p, 1-p))
      if (length(indm20[indm20==1])>=3&&length(indm20[indm20==2])>=3)
      {
        trainm20 <- data.frame(m20[indm20==1,])
        testm20 <- as.data.frame(m20[indm20==2,])
        break
      }
      
    }
  }else
  {
    trainm20 <- data.frame(m20[indm20==1,])
    testm20 <- as.data.frame(m20[indm20==2,])
  }
  
  indm21 <- sample(2, nrow(m21), replace = TRUE, prob = c(p, 1-p))
  if (length(indm21[indm21==1])<3||length(indm21[indm21==2])<3)
  {
    while (length(indm21[indm21==1])<3||length(indm21[indm21==2])<3)
    {
      indm21 <- sample(2, nrow(m21), replace = TRUE, prob = c(p, 1-p))
      if (length(indm21[indm21==1])>=3&&length(indm21[indm21==2])>=3)
      {
        trainm21 <- data.frame(m21[indm21==1,])
        testm21 <- as.data.frame(m21[indm21==2,])
        break
      }
      
    }
  }else
  {
    trainm21 <- data.frame(m21[indm21==1,])
    testm21 <- as.data.frame(m21[indm21==2,])
  }
  
  
  trainm22 <- rbind(trainm20,trainm21)
  #trainm2 <- data5
  
  
  testm22 <- rbind(testm20,testm21)
  #testm2 <- data5
  
  
  library('randomForest')
  library('pROC')
  rfm2 <- randomForest(gr~., data=trainm22, ntree = 500,mtry=3,importance=TRUE,proximity=TRUE)
  p2m2 <- predict(rfm2, testm22)
  #obs_p_ranm2 = data.frame(prob=p2m2,obs=testm2$cluster)
  ran_rocm2 <- roc(testm22$gr,as.numeric(p2m2))
  #drx2[nxb]=as.numeric(ran_rocm2$auc)
  drx2=as.numeric(ran_rocm2$auc)
  names(drx2)='AUC'
  drf2=importance(rfm2)
  dip1x2=drf2[,(ncol(drf2)-1)]
  dippx2=as.data.frame(c(dip1x2,drx2))
  colnames(dippx2)=paste("round",as.character(nxb),sep='_')
  dipx2=cbind(dipx2,dippx2)
  
}
dipx2$prob=rep(p,nrow(dipx2))
dcb=rbind(dcb,dipx2)
f4=data.frame(t(dcb))
da4=data.frame(AUC=f4[,c('AUC')])

p<-ggplot(data = da4) + geom_histogram(aes(x = AUC, y = ..density..),fill = "red", alpha = 0.2) +
  geom_density(aes(x= AUC, y = ..density..), col = "red")
p 
f5=f4
d=data.frame()
for (i in 1:(ncol(f5)-1))
{
  d1=data.frame(importance=f5[,i],factors=rep(colnames(f5)[i],nrow(f5)))
  d=rbind(d,d1)
}
d
sssi=names(sort(tapply(d$importance,d$factors,median),decreasing = T))
library(ggsignif)
library(ggpubr)
library(ggplot2)
p=ggplot(d,aes(x=factor(factors,levels = sssi),y=importance,fill=factors),color='#C3241F')+geom_boxplot()+
  geom_hline(yintercept = 0, color='red',linetype = 2)+
  #scale_colour_manual(values = rev(c('#FF6347','#1E90FF')))+
  #scale_fill_manual(values = rev(c('#C3241F','#195FAC')))+
  #stat_compare_means(aes(group = group),method='wilcox.test'
  #                   ,label="p.signif"
  #)+
  #ylim(c(0,0.2))+
  theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank()
        ,panel.background = element_blank(),legend.position = 'none',
        axis.text.x = element_text( hjust = 1,vjust= 1,angle =30,size = 15),
        axis.line = element_line(colour = "black"))
p 


fg=read.table('PFgenes.pharmmapper_swisstargetpredict.stitch.1.txt',header=F)
#fg=read.table('immuGene.merged.txt',header=F)
#fg=rbind(fg1,fg2)
mit=fg$V1
ls1=unique(dmimln$node2)
ls2=mit

library (VennDiagram) 
l1=list()
l1[[1]]=ls1
l1[[2]]=ls2

names(l1)=c('markergene','oxsidtive')
T1<-venn.diagram(l1,
                 filename=NULL,
                 lwd=1,#圈线粗度
                 lty=1, #圈线类型
                 col=c("Gold","Turquoise4",'red','blue','purple')[1:length(l1)], #圈线颜色
                 fill=c("Gold","Turquoise4",'red','blue','purple')[1:length(l1)], #填充颜色
                 cat.col=c("Gold","Turquoise4",'red','blue','purple')[1:length(l1)],#A和B的颜色
                 cat.cex = 2,# A和B的大小
                 rotation.degree = 0,#旋转角度
                 main.cex = 2,#主标题大小
                 cex=1.5,#里面交集字的大小
                 alpha = 0.5,#透明度 
                 reverse=TRUE)
p <- ggplot()
p
grid.draw(T1)


oxgenes=ls1[ls1%in%ls2]
dmimln1=dmimln[dmimln$node2%in%oxgenes,]
write.table(dmimln1,'oxidative-cerna.net.txt',row.names = F,col.names = T,sep='\t')

rnaExpr1=t(exprm1ss[exprm1ss$gr==1,oxgenes])
rnaExpr2=as.matrix(rnaExpr1)
imExpr=rnaExpr2[rownames(rnaExpr2)%in%oxgenes,]
mads=apply(imExpr,1,mad)
imExpr=imExpr[rev(order(mads)),]
imExpr=na.omit(imExpr)
d = sweep(imExpr,1, apply(imExpr,1,median))

library(ConsensusClusterPlus)
title='D:\\UC\\merge'
results = ConsensusClusterPlus(d,maxK=10,reps=1000,pItem=0.8,pFeature=1,
                               title=title,clusterAlg="km",distance="euclidean",seed=1262118388.71279,plot="png")
#title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

nxn=2
icl = calcICL(results,title=title,plot="png")
results[[nxn]]$consensusClass
write.table(results[[nxn]]$consensusClass,'cluster-samples-lnc-mRNA.txt',sep='\t')


timmExpr=t(imExpr[rownames(imExpr)%in%oxgenes,])
timmExpr=na.omit(timmExpr)
pca <- prcomp(timmExpr, scale=TRUE,rank=4)#,retx=T)
pca.var <- pca$sdev^2 
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       cluster=factor(results[[nxn]]$consensusClass[rownames(pca$x)]))

mcolor=c("Gold",'blue',"Turquoise4","Firebrick")
library(vegan)
AN=anosim(timmExpr,factor(results[[nxn]]$consensusClass[rownames(pca$x)]),permutations = 999, distance = "bray")
ggplot(data=pca.data,aes(x=X,y=Y,color=cluster))+
  geom_point(size=3)+
  geom_text(x=max(pca.data$X)-1,y=max(pca.data$Y)-1,label=paste("R=",round(AN$statistic,3),"  ","P=",round(AN$signif,3),sep=''),color='black')+
  stat_ellipse()+
  theme_bw()+theme(panel.grid=element_blank())+
  scale_color_manual(values = c("Gold","Chocolate","Turquoise4","red",'blue','pink'))+
  xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
  ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))+
  theme(legend.position = 'top')

library(rgl)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
pcad <- prcomp(timmExpr)
pcad1 <- data.frame(pcad$x)
groups <-  results[[nxn]]$consensusClass[rownames(pcad1)] 
sco=c("Gold","Chocolate","Turquoise4","Firebrick",'blue')
groups[groups==1]="Gold"
groups[groups==2]="Chocolate"
groups[groups==3]="Turquoise4"
groups[groups==4]="red"
groups[groups==5]="blue"
#groups[groups==6]="pink"
groups=groups[rownames(pcad1)]
plot3d(pcad1[,1:3], # 取前三个主成分
       xlab="Comp.1", ylab="Comp.2", zlab="Comp.3", 
       col=groups, # 按groups填充颜色
       type="s", # 画球，'p' for points, 's' for spheres, 'l' for lines, 'h' for line segments 
       size=2, #球的大小
       lwd=2, box=T)

rgl.snapshot("PCA01.png")

slc2=slc1[colnames(imExpr)]

drgg=data.frame(slc2)
drgg$cluster=results[[nxn]]$consensusClass[rownames(drgg)]

library(ggplot2)
library(dplyr)
d2=data.frame(cluster=drgg$cluster,group=drgg$slc2)
rownames(d2)=rownames(drgg)
ds=data.frame()
for (i in 1:length(unique(as.character(d2$cluster))))
{
  
  dx=d2[d2$cluster==unique(as.character(d2$cluster))[i],]
  print(dx)
  counti=nrow(dx)#count in i
  countnoti=nrow(d2)-counti# count in noti
  dc=data.frame()
  for (j in 1:length(unique(as.vector(d2$group))))
  {
    #       i   noti
    #  j    x1  x3
    # notj  x2  x4
    countj=nrow(dx[dx$group==unique(as.vector(d2$group))[j],])# count j in cluster i  x1
    x1=countj
    countnotj=counti-countj#count not j in cluster  x2
    x2=countnotj
    totalj=nrow(d2[d2$group==unique(as.vector(d2$group))[j],])# total j in all
    coutjinnoti=totalj-countj#countj not in i x3
    x3=coutjinnoti
    countnotinotj=countnoti-coutjinnoti# count not j phase in not i cluster x4
    x4=countnotinotj
    #matrix(c(x1,x2,x3,x4),nrow=2)
    p=fisher.test(matrix(c(x1,x2,x3,x4),nrow=2),alternative='greater')$p.value
    dc[j,1]=unique(as.vector(d2$group))[j]
    dc[j,2]=unique(as.vector(d2$cluster))[i]
    dc[j,3]=x1
    dc[j,4]=x2
    dc[j,5]=x3
    dc[j,6]=x4
    dc[j,7]=p
  }
  colnames(dc)=c('group','Cluster','group_in_cluster','not_group_in_cluster','group_in_not_cluster','not_group_in_not_clusters','Pvalue')
  ds=rbind(ds,dc)
}
p=ds$Pvalue
pa=p.adjust(p,method = 'fdr',n=length(p))
ds$fdr=pa
ds

f3=d2
tf3=table(f3$cluster)
f31=as.data.frame(f3 %>%
                    group_by_at( vars( matches("cluster|group") ) ) %>%count())
f41=f31
myco1=c('#20B2AA','#FF8C00','#00BFFF','FireBrick','#9400D3','#1E90FF')
x=as.vector(unique(f41$cluster))
x1=as.vector(unique(f41$group))
d=data.frame()
for (i in 1:length(x))
{
  d1=f41[f41$cluster==x[i],]
  print(d1)
  
  for (j in 1:length(x1))
  {
    #print(x1[j])
    if (!(x1[j] %in% as.vector(unique(d1$group))))
    {
      print(x1[j])
      d2=data.frame()
      d2[1,1]=x[i]
      d2[1,2]=x1[j]
      d2[1,3]=0
      colnames(d2)=c('cluster','group', 'n')
      print(d2)
      d1=rbind(d1,d2)
    }
  }
  d=rbind(d,d1)
}
d
ggplot(d,aes(x=factor(cluster),y=n,fill=factor(group)))+geom_bar(stat = "identity",position=position_dodge(),width=0.5)+
  scale_fill_manual(values=myco1)+scale_colour_manual(values=myco1)+
  theme(axis.text.x = element_text(size = 11),
        legend.position="top",
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        ,panel.border = element_blank()
  )
write.table(ds,'risk_enriched_in_cluster.fisher-test.lnc-mrna.txt',col.names=T,row.names=F,sep='\t')


exprm1x=exprm1[,names(slc2)]
group=paste('cluster',results[[nxn]]$consensusClass[colnames(exprm1x)],sep='')
names(group)=colnames(exprm1x)
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
contrast.matrix<-makeContrasts("cluster1-cluster2",
                               levels = design)
contrast.matrix
ddmrna<-as.data.frame(t(exprm1x))
res.pca <- prcomp(ddmrna, scale = TRUE)
fviz_eig(res.pca)
color=c("#FC4E07","#00AFBB")
fviz_pca_ind(res.pca,
             col.ind = factor(group[rownames(ddmrna)]), # 颜色对应group信息
             palette = color,
             addEllipses = TRUE, # Concentration ellipses
             ellipse.type = "confidence",
             legend.title = "mRNA Group",## Legend名称
             repel = TRUE
)

fit <- lmFit(exprm1x,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEGmRNA<-topTable(fit2, coef=1, n=Inf) %>% na.omit()

require(EnhancedVolcano)
EnhancedVolcano(DEGmRNA,title="DEGmRNA",
                lab = rownames(DEGmRNA),
                x = "logFC",
                y = "P.Value",
                selectLab = rownames(DEGmRNA)[1:5],
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,## pvalue阈值
                FCcutoff = 1,## FC cutoff
                xlim = c(-5,5),
                colAlpha = 1,
                legendPosition = "bottom",
                legendLabSize = 10,
                legendIconSize = 3.0)
write.table(DEGmRNA,'GSE75214.activ.DEG_mRNA.cluster1-2.txt',row.names = T,col.names = T,sep='\t')
nrow(DEGmRNA[DEGmRNA$P.Value<0.05&abs(DEGmRNA$logFC)>1,])
DEGmRNA1=DEGmRNA[DEGmRNA$P.Value<0.05&abs(DEGmRNA$logFC)>1,]
write.table(DEGmRNA1,'GSE75214.activ.DEG_mRNA.cluster1-2.sig.txt',row.names = T,col.names = T,sep='\t')
index<-DEGmRNA1 %>% as_tibble() %>% 
  mutate(genename=rownames(DEGmRNA1)) %>% 
  dplyr::arrange(desc(logFC)) %>% 
  .$genename 

library(pheatmap)
exprm1s=exprm1x[,names(sort(group))]
index_matrix<-t(scale(t(exprm1s[index,])))##归一化
index_matrix[index_matrix>2]=2
index_matrix[index_matrix<(-2)]=-2
head(index_matrix)
anno=data.frame(group=sort(group))
rownames(anno)=colnames(index_matrix)
pheatmap(index_matrix,
         show_colnames =F,
         show_rownames = T,
         cluster_cols = F, 
         annotation_col=anno)

dsmi=exprm1x[rownames(DEGmRNA1),]
dsmi1=as.data.frame(t(dsmi))
dsmrna=exprm1x[oxgenes,rownames(dsmi1)]
dsmrna1=as.data.frame(t(dsmrna))
dmim=data.frame()
for (i in 1:ncol(dsmi1))
{
  for (j in 1:ncol(dsmrna1))
  {
    print(i)
    print(j)
    dtmp=data.frame()
    dtmp[1,1]=colnames(dsmi1)[i]
    dtmp[1,2]=colnames(dsmrna1)[j]
    dtmp[1,3]=cor.test(as.numeric(dsmi1[,i]),as.numeric(dsmrna1[,j]))$estimate
    dtmp[1,4]=cor.test(as.numeric(dsmi1[,i]),as.numeric(dsmrna1[,j]))$p.value
    dmim=rbind(dmim,dtmp)
  }
}
nrow(dmim)
colnames(dmim)=c('degmrna','oxditivegenes','R','pvalue')
dmim$fdr=p.adjust(dmim$pvalue,n=nrow(dmim))
dmim1=dmim[dmim$fdr<0.05,]
dmim1=dmim1[abs(dmim1$R)>0.8,]
write.table(dmim1,'correlation_DEGactive.cluster1-2.fdr.07.txt',row.names = F,col.names = T,sep='\t')
ac12s=DEGmRNA1[rownames(DEGmRNA1)%in%c(dmim1$degmrna,dmim1$oxditivegenes),]
write.table(ac12s,'DEGactive.cluster1-2.siggene.net.txt',row.names = T,col.names = T,sep='\t')

dam=data.frame()
for (i in 1:nrow(DEGmRNA1))
{
  dall=data.frame()
  dall[1,1]=rownames(DEGmRNA1)[i]
  if(DEGmRNA1[i,1]>0)
  {
    dall[1,2]='up'
  }else{
    dall[1,2]='down'
  }
  dall[1,3]='cluster1'
  dam=rbind(dam,dall)
}
dam
dal=data.frame()
dal[1,1]="CAV1"
dal[1,2]="even"
dal[1,3]="oxditive"
dalm=rbind(dal,dam)
colnames(dalm)=c('node','class','group')
write.table(dalm,'correlation_DEGactive.cluster1-2.fdr.07.annotation.txt',row.names = F,col.names = T,sep='\t')
######################################active-mirna-mrna###########
rnaExpr1=t(exprmmx1ss[exprmmx1ss$gr==1,oxgenes])
rnaExpr2=as.matrix(rnaExpr1)
imExpr=rnaExpr2[rownames(rnaExpr2)%in%oxgenes,]
mads=apply(imExpr,1,mad)
imExpr=imExpr[rev(order(mads)),]
imExpr=na.omit(imExpr)
d = sweep(imExpr,1, apply(imExpr,1,median))

library(ConsensusClusterPlus)
title='D:\\UC\\merge1'
results = ConsensusClusterPlus(d,maxK=10,reps=1000,pItem=0.8,pFeature=1,
 #                              title=title,clusterAlg="km",distance="euclidean",seed=1262118388.71279,plot="png")
title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

nxn=2
icl = calcICL(results,title=title,plot="png")
results[[nxn]]$consensusClass
write.table(results[[nxn]]$consensusClass,'cluster-samples-miRNA-mRNA.txt',sep='\t')


timmExpr=t(imExpr[rownames(imExpr)%in%oxgenes,])
timmExpr=na.omit(timmExpr)
pca <- prcomp(timmExpr, scale=TRUE,rank=4)#,retx=T)
pca.var <- pca$sdev^2 
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       cluster=factor(results[[nxn]]$consensusClass[rownames(pca$x)]))

mcolor=c("Gold",'blue',"Turquoise4","Firebrick")
library(vegan)
AN=anosim(timmExpr,factor(results[[nxn]]$consensusClass[rownames(pca$x)]),permutations = 999, distance = "bray")
ggplot(data=pca.data,aes(x=X,y=Y,color=cluster))+
  geom_point(size=3)+
  geom_text(x=max(pca.data$X)-1,y=max(pca.data$Y)-1,label=paste("R=",round(AN$statistic,3),"  ","P=",round(AN$signif,3),sep=''),color='black')+
  stat_ellipse()+
  theme_bw()+theme(panel.grid=element_blank())+
  scale_color_manual(values = c("Gold","Chocolate","Turquoise4","red",'blue','pink'))+
  xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
  ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))+
  theme(legend.position = 'top')

library(rgl)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
pcad <- prcomp(timmExpr)
pcad1 <- data.frame(pcad$x)
groups <-  results[[nxn]]$consensusClass[rownames(pcad1)] 
sco=c("Gold","Chocolate","Turquoise4","Firebrick",'blue')
groups[groups==1]="Gold"
groups[groups==2]="Chocolate"
groups[groups==3]="Turquoise4"
groups[groups==4]="red"
groups[groups==5]="blue"
#groups[groups==6]="pink"
groups=groups[rownames(pcad1)]
plot3d(pcad1[,1:3], # 取前三个主成分
       xlab="Comp.1", ylab="Comp.2", zlab="Comp.3", 
       col=groups, # 按groups填充颜色
       type="s", # 画球，'p' for points, 's' for spheres, 'l' for lines, 'h' for line segments 
       size=2, #球的大小
       lwd=2, box=T)

rgl.snapshot("PCA02.png")

smmc2=smmc1[colnames(imExpr)]

drgg=data.frame(smmc2)
drgg$cluster=results[[nxn]]$consensusClass[rownames(drgg)]

library(ggplot2)
library(dplyr)
d2=data.frame(cluster=drgg$cluster,group=drgg$smmc2)
rownames(d2)=rownames(drgg)
ds=data.frame()
for (i in 1:length(unique(as.character(d2$cluster))))
{
  
  dx=d2[d2$cluster==unique(as.character(d2$cluster))[i],]
  print(dx)
  counti=nrow(dx)#count in i
  countnoti=nrow(d2)-counti# count in noti
  dc=data.frame()
  for (j in 1:length(unique(as.vector(d2$group))))
  {
    #       i   noti
    #  j    x1  x3
    # notj  x2  x4
    countj=nrow(dx[dx$group==unique(as.vector(d2$group))[j],])# count j in cluster i  x1
    x1=countj
    countnotj=counti-countj#count not j in cluster  x2
    x2=countnotj
    totalj=nrow(d2[d2$group==unique(as.vector(d2$group))[j],])# total j in all
    coutjinnoti=totalj-countj#countj not in i x3
    x3=coutjinnoti
    countnotinotj=countnoti-coutjinnoti# count not j phase in not i cluster x4
    x4=countnotinotj
    #matrix(c(x1,x2,x3,x4),nrow=2)
    p=fisher.test(matrix(c(x1,x2,x3,x4),nrow=2),alternative='greater')$p.value
    dc[j,1]=unique(as.vector(d2$group))[j]
    dc[j,2]=unique(as.vector(d2$cluster))[i]
    dc[j,3]=x1
    dc[j,4]=x2
    dc[j,5]=x3
    dc[j,6]=x4
    dc[j,7]=p
  }
  colnames(dc)=c('group','Cluster','group_in_cluster','not_group_in_cluster','group_in_not_cluster','not_group_in_not_clusters','Pvalue')
  ds=rbind(ds,dc)
}
p=ds$Pvalue
pa=p.adjust(p,method = 'fdr',n=length(p))
ds$fdr=pa
ds

f3=d2
tf3=table(f3$cluster)
f31=as.data.frame(f3 %>%
                    group_by_at( vars( matches("cluster|group") ) ) %>%count())
f41=f31
myco1=c('#20B2AA','#FF8C00','#00BFFF','FireBrick','#9400D3','#1E90FF')
x=as.vector(unique(f41$cluster))
x1=as.vector(unique(f41$group))
d=data.frame()
for (i in 1:length(x))
{
  d1=f41[f41$cluster==x[i],]
  print(d1)
  
  for (j in 1:length(x1))
  {
    #print(x1[j])
    if (!(x1[j] %in% as.vector(unique(d1$group))))
    {
      print(x1[j])
      d2=data.frame()
      d2[1,1]=x[i]
      d2[1,2]=x1[j]
      d2[1,3]=0
      colnames(d2)=c('cluster','group', 'n')
      print(d2)
      d1=rbind(d1,d2)
    }
  }
  d=rbind(d,d1)
}
d
ggplot(d,aes(x=factor(cluster),y=n,fill=factor(group)))+geom_bar(stat = "identity",position=position_dodge(),width=0.5)+
  scale_fill_manual(values=myco1)+scale_colour_manual(values=myco1)+
  theme(axis.text.x = element_text(size = 11),
        legend.position="top",
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        ,panel.border = element_blank()
  )
write.table(ds,'risk_enriched_in_cluster.fisher-test.mir-mrna.txt',col.names=T,row.names=F,sep='\t')

###################################therapy-response################

exprmx=read.table('GSE114603-GPL14951_series_matrix.miRNA_mRNA.exp.1.txt',row.names = 1,header=T,sep='\t')
mprob=read.table('GSE114603-GPL14951_series_matrix.miRNA_mRNA.prob.1.pc.txt',header=F)
probm=mprob$V2
names(probm)=mprob$V1
exprmmx=exprmx[rownames(exprmx)%in%names(probm),]
exprmmx$gene=probm[rownames(exprmmx)]
sam=na.omit(unique(exprmmx$gene))
#exprmmx1=data.frame()
#for (i in sam)
#{
#  dbs=na.omit(exprmmx[exprmmx$gene%in%c(i),])
#  dbx=dbs[,1:(ncol(dbs)-1)]
#  dbx1=t(data.frame(apply(dbx,2,median)))
#  rownames(dbx1)=i
#  exprmmx1=rbind(exprmmx1,dbx1)
#}
#exprmmx1
#write.table(exprmmx1,'GSE114603-GPL14951.mrna.expression.1.txt',row.names = T,col.names = T,sep = '\t')
exprmmx1=read.table('GSE114603-GPL14951.mrna.expression.1.txt',row.names = 1,header=T,sep='\t')

fcmi=read.table('GSE114603-GPL14951_series_matrix.miRNA_mRNA.cli.2.txt',header=T,row.names = 1)
scsss=fcmi$group2
names(scsss)=rownames(fcmi)
exprmmx11=as.matrix(exprmmx1[,names(scsss)])
imExpr=exprmmx11[rownames(exprmmx11)%in%oxgenes,]
mads=apply(imExpr,1,mad)
imExpr=imExpr[rev(order(mads)),]
imExpr=na.omit(imExpr)
d = sweep(imExpr,1, apply(imExpr,1,median))

library(ConsensusClusterPlus)
title='D:\\UC\\merge2'
results = ConsensusClusterPlus(d,maxK=10,reps=1000,pItem=0.8,pFeature=1,
                               #                              title=title,clusterAlg="km",distance="euclidean",seed=1262118388.71279,plot="png")
                               title=title,clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

nxn=3
icl = calcICL(results,title=title,plot="png")
results[[nxn]]$consensusClass
write.table(results[[nxn]]$consensusClass,'GSE114603-GPL14951.cluster-samples-miRNA-mRNA.txt',sep='\t')


timmExpr=t(imExpr[rownames(imExpr)%in%oxgenes,])
timmExpr=na.omit(timmExpr)
pca <- prcomp(timmExpr, scale=TRUE,rank=4)#,retx=T)
pca.var <- pca$sdev^2 
pca.var.per <- round(pca.var/sum(pca.var)*100, 1)
pca.data <- data.frame(Sample=rownames(pca$x),
                       X=pca$x[,1],
                       Y=pca$x[,2],
                       cluster=factor(results[[nxn]]$consensusClass[rownames(pca$x)]))

mcolor=c("Gold",'blue',"Turquoise4","Firebrick")
library(vegan)
AN=anosim(timmExpr,factor(results[[nxn]]$consensusClass[rownames(pca$x)]),permutations = 999, distance = "bray")
ggplot(data=pca.data,aes(x=X,y=Y,color=cluster))+
  geom_point(size=3)+
  geom_text(x=max(pca.data$X)-1,y=max(pca.data$Y)-1,label=paste("R=",round(AN$statistic,3),"  ","P=",round(AN$signif,3),sep=''),color='black')+
  stat_ellipse()+
  theme_bw()+theme(panel.grid=element_blank())+
  scale_color_manual(values = c("Gold","Chocolate","Turquoise4","red",'blue','pink'))+
  xlab(paste("PC1(",pca.var.per[1],"%","variance)",sep=""))+
  ylab(paste("PC2(",pca.var.per[2],"%","variance)",sep=""))+
  theme(legend.position = 'top')

library(rgl)
Sys.setenv(LANGUAGE = "en") #显示英文报错信息
options(stringsAsFactors = FALSE) #禁止chr转成factor
pcad <- prcomp(timmExpr)
pcad1 <- data.frame(pcad$x)
groups <-  results[[nxn]]$consensusClass[rownames(pcad1)] 
sco=c("Gold","Chocolate","Turquoise4","Firebrick",'blue')
groups[groups==1]="Gold"
groups[groups==2]="Chocolate"
groups[groups==3]="Turquoise4"
groups[groups==4]="red"
groups[groups==5]="blue"
#groups[groups==6]="pink"
groups=groups[rownames(pcad1)]
plot3d(pcad1[,1:3], # 取前三个主成分
       xlab="Comp.1", ylab="Comp.2", zlab="Comp.3", 
       col=groups, # 按groups填充颜色
       type="s", # 画球，'p' for points, 's' for spheres, 'l' for lines, 'h' for line segments 
       size=2, #球的大小
       lwd=2, box=T)

rgl.snapshot("PCA02.png")

scsss=scsss[colnames(imExpr)]

drgg=data.frame(scsss)
drgg$cluster=results[[nxn]]$consensusClass[rownames(drgg)]

library(ggplot2)
library(dplyr)
d2=data.frame(cluster=drgg$cluster,group=drgg$scsss)
rownames(d2)=rownames(drgg)
ds=data.frame()
for (i in 1:length(unique(as.character(d2$cluster))))
{
  
  dx=d2[d2$cluster==unique(as.character(d2$cluster))[i],]
  print(dx)
  counti=nrow(dx)#count in i
  countnoti=nrow(d2)-counti# count in noti
  dc=data.frame()
  for (j in 1:length(unique(as.vector(d2$group))))
  {
    #       i   noti
    #  j    x1  x3
    # notj  x2  x4
    countj=nrow(dx[dx$group==unique(as.vector(d2$group))[j],])# count j in cluster i  x1
    x1=countj
    countnotj=counti-countj#count not j in cluster  x2
    x2=countnotj
    totalj=nrow(d2[d2$group==unique(as.vector(d2$group))[j],])# total j in all
    coutjinnoti=totalj-countj#countj not in i x3
    x3=coutjinnoti
    countnotinotj=countnoti-coutjinnoti# count not j phase in not i cluster x4
    x4=countnotinotj
    #matrix(c(x1,x2,x3,x4),nrow=2)
    p=fisher.test(matrix(c(x1,x2,x3,x4),nrow=2),alternative='greater')$p.value
    dc[j,1]=unique(as.vector(d2$group))[j]
    dc[j,2]=unique(as.vector(d2$cluster))[i]
    dc[j,3]=x1
    dc[j,4]=x2
    dc[j,5]=x3
    dc[j,6]=x4
    dc[j,7]=p
  }
  colnames(dc)=c('group','Cluster','group_in_cluster','not_group_in_cluster','group_in_not_cluster','not_group_in_not_clusters','Pvalue')
  ds=rbind(ds,dc)
}
p=ds$Pvalue
pa=p.adjust(p,method = 'fdr',n=length(p))
ds$fdr=pa
ds

f3=d2
tf3=table(f3$cluster)
f31=as.data.frame(f3 %>%
                    group_by_at( vars( matches("cluster|group") ) ) %>%count())
f41=f31
myco1=c('#20B2AA','#FF8C00','#00BFFF','FireBrick','#9400D3','#1E90FF')
x=as.vector(unique(f41$cluster))
x1=as.vector(unique(f41$group))
d=data.frame()
for (i in 1:length(x))
{
  d1=f41[f41$cluster==x[i],]
  print(d1)
  
  for (j in 1:length(x1))
  {
    #print(x1[j])
    if (!(x1[j] %in% as.vector(unique(d1$group))))
    {
      print(x1[j])
      d2=data.frame()
      d2[1,1]=x[i]
      d2[1,2]=x1[j]
      d2[1,3]=0
      colnames(d2)=c('cluster','group', 'n')
      print(d2)
      d1=rbind(d1,d2)
    }
  }
  d=rbind(d,d1)
}
d
ggplot(d,aes(x=factor(cluster),y=n,fill=factor(group)))+geom_bar(stat = "identity",position=position_dodge(),width=0.5)+
  scale_fill_manual(values=myco1)+scale_colour_manual(values=myco1)+
  theme(axis.text.x = element_text(size = 11),
        legend.position="top",
        panel.grid.major =element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank()
        ,panel.border = element_blank()
  )
write.table(ds,'risk_enriched_in_cluster.fisher-test.therapy.txt',col.names=T,row.names=F,sep='\t')

scs=results[[nxn]]$consensusClass
scsg=scs[scs%in%c(1,2,3)]
exprm1psl=exprmmx1[,names(scsg)]
scsd=data.frame(scsg)
scs1=scsg
scs1[scs1==1]="y"
scs1[scs1!="y"]="n"
scsd$scs1=scs1[rownames(scsd)]
group=scsd$scs1
names(group)=rownames(scsd)

design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
contrast.matrix<-makeContrasts("y-n",
                               levels = design)
contrast.matrix
ddmrna<-as.data.frame(t(exprm1psl))

fit <- lmFit(exprm1psl,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEGmRNAclu1<-topTable(fit2, coef=1, n=Inf) %>% na.omit()

require(EnhancedVolcano)
EnhancedVolcano(DEGmRNAclu1,title="Represent cluster 1",
                lab = rownames(DEGmRNAclu1),
                x = "logFC",
                y = "P.Value",
                selectLab = rownames(DEGmRNAclu1)[1:5],
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,## pvalue阈值
                FCcutoff = 0.58,## FC cutoff
                xlim = c(-5,5),
                colAlpha = 1,
                legendPosition = "bottom",
                legendLabSize = 10,
                legendIconSize = 3.0)
DEGmRNAclu11=DEGmRNAclu1[DEGmRNAclu1$P.Value<0.05&DEGmRNAclu1$logFC>0.58,]
nrow(DEGmRNAclu11)
rpg1=rownames(DEGmRNAclu11)
write.table(DEGmRNAclu11,'DEG_mRNA_cluster1.txt',row.names = T,col.names = T,sep='\t')

scs2=scsg
scs2[scs2==2]="y"
scs2[scs2!="y"]="n"
scsd$scs2=scs2[rownames(scsd)]
group=scsd$scs2
names(group)=rownames(scsd)

design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
contrast.matrix<-makeContrasts("y-n",
                               levels = design)
contrast.matrix
ddmrna<-as.data.frame(t(exprm1psl))

fit <- lmFit(exprm1psl,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEGmRNAclu2<-topTable(fit2, coef=1, n=Inf) %>% na.omit()

require(EnhancedVolcano)
EnhancedVolcano(DEGmRNAclu2,title="Represent cluster 1",
                lab = rownames(DEGmRNAclu1),
                x = "logFC",
                y = "P.Value",
                selectLab = rownames(DEGmRNAclu1)[1:5],
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,## pvalue阈值
                FCcutoff = 0.58,## FC cutoff
                xlim = c(-5,5),
                colAlpha = 1,
                legendPosition = "bottom",
                legendLabSize = 10,
                legendIconSize = 3.0)
DEGmRNAclu12=DEGmRNAclu2[DEGmRNAclu2$P.Value<0.05&DEGmRNAclu2$logFC>0.58,]
nrow(DEGmRNAclu12)
rpg2=rownames(DEGmRNAclu12)
write.table(DEGmRNAclu12,'DEG_mRNA_cluster2.txt',row.names = T,col.names = T,sep='\t')


scs3=scsg
scs3[scs3==3]="y"
scs3[scs3!="y"]="n"
scsd$scs3=scs3[rownames(scsd)]
group=scsd$scs3
names(group)=rownames(scsd)

design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
contrast.matrix<-makeContrasts("y-n",
                               levels = design)
contrast.matrix
ddmrna<-as.data.frame(t(exprm1psl))

fit <- lmFit(exprm1psl,design)
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2) 
DEGmRNAclu3<-topTable(fit2, coef=1, n=Inf) %>% na.omit()

require(EnhancedVolcano)
EnhancedVolcano(DEGmRNAclu3,title="Represent cluster 1",
                lab = rownames(DEGmRNAclu3),
                x = "logFC",
                y = "P.Value",
                selectLab = rownames(DEGmRNAclu1)[1:5],
                xlab = bquote(~Log[2]~ "fold change"),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,## pvalue阈值
                FCcutoff = 0.58,## FC cutoff
                xlim = c(-5,5),
                colAlpha = 1,
                legendPosition = "bottom",
                legendLabSize = 10,
                legendIconSize = 3.0)
DEGmRNAclu13=DEGmRNAclu3[DEGmRNAclu3$P.Value<0.05&DEGmRNAclu3$logFC>0.58,]
nrow(DEGmRNAclu13)
write.table(DEGmRNAclu13,'DEG_mRNA_cluster3.txt',row.names = T,col.names = T,sep='\t')

rpg1=data.frame(rg=rownames(DEGmRNAclu11))
rpg1$cluster=rep('cluster1',nrow(rpg1))

rpg2=data.frame(rg=rownames(DEGmRNAclu12))
rpg2$cluster=rep('cluster2',nrow(rpg2))

rpg3=data.frame(rg=rownames(DEGmRNAclu13))
rpg3$cluster=rep('cluster3',nrow(rpg3))

ox=data.frame(rg=oxgenes)
ox$cluster=rep('oxgene',nrow(ox))

rpgs=rbind(rpg1,rpg2,rpg3,ox)
write.table(rpgs,'therapy-representive-cluster-genes.txt',col.names = T,sep='\t')
dbb=data.frame()
for (x in unique(rpgs$rg))
{
  dx=rpgs[rpgs$rg==x,]
  db=data.frame()
  db[1,1]=x
  db[1,2]=0
  db[1,3]=0
  db[1,4]=0
  db[1,5]=0
  colnames(db)=c('gene','cluster1','cluster2','cluster3','oxgene')
  sv=rep(1,nrow(dx))
  names(sv)=dx$cluster
  db[1,names(sv)]=sv
  dbb=rbind(db,dbb)
}
head(dbb)
nrow(dbb)

rpg123=unique(c(rownames(DEGmRNAclu13),rownames(DEGmRNAclu12),rownames(DEGmRNAclu11)))
dsmi=exprm1psl[rpg123,]
dsmi1=as.data.frame(t(dsmi))
dsmrna=exprm1psl[oxgenes,rownames(dsmi1)]
dsmrna1=as.data.frame(t(dsmrna))
dmim=data.frame()
for (i in 1:ncol(dsmi1))
{
  for (j in 1:ncol(dsmrna1))
  {
    print(i)
    print(j)
    dtmp=data.frame()
    dtmp[1,1]=colnames(dsmi1)[i]
    dtmp[1,2]=colnames(dsmrna1)[j]
    dtmp[1,3]=cor.test(as.numeric(dsmi1[,i]),as.numeric(dsmrna1[,j]))$estimate
    dtmp[1,4]=cor.test(as.numeric(dsmi1[,i]),as.numeric(dsmrna1[,j]))$p.value
    dmim=rbind(dmim,dtmp)
  }
}
nrow(dmim)
colnames(dmim)=c('rgg','oxditivegenes','R','pvalue')
dmim$fdr=p.adjust(dmim$pvalue,n=nrow(dmim))
dmim1=dmim[dmim$fdr<0.05,]
dmim1=dmim1[abs(dmim1$R)>0.8,]
write.table(dmim1,'correlation_resisteance.cluster1-2-2.fdr.08.txt',row.names = F,col.names = T,sep='\t')

rggn=rpgs[rpgs$rg%in%(c(dmim1$rgg,dmim1$oxditivegenes)),]
write.table(rggn,'therapy-representive-cluster-genes..signet.txt',col.names = T,sep='\t')

c1=rggn[rggn$cluster=='cluster1',]$rg
c2=rggn[rggn$cluster=='cluster2',]$rg
c3=rggn[rggn$cluster=='cluster3',]$rg
oxx=rggn[rggn$cluster=='oxgene',]$rg

library (VennDiagram) 
l1=list()
l1[[1]]=c1
l1[[2]]=c2
l1[[3]]=c3
l1[[4]]=oxx
names(l1)=c('cluster1','cluster2','cluster3','oxgenes')
T1<-venn.diagram(l1,
                 filename=NULL,
                 lwd=1,#圈线粗度
                 lty=1, #圈线类型
                 col=c("Gold","Turquoise4",'red','blue','purple')[1:length(l1)], #圈线颜色
                 fill=c("Gold","Turquoise4",'red','blue','purple')[1:length(l1)], #填充颜色
                 cat.col=c("Gold","Turquoise4",'red','blue','purple')[1:length(l1)],#A和B的颜色
                 cat.cex = 2,# A和B的大小
                 rotation.degree = 0,#旋转角度
                 main.cex = 2,#主标题大小
                 cex=1.5,#里面交集字的大小
                 alpha = 0.5,#透明度 
                 reverse=TRUE)
p <- ggplot()
p
grid.draw(T1)














