library(ggplot2)
library(phyloseq)
library(tidyverse)
library(doParallel)
library(multcomp) 
library(stringr) 
library(dplyr)
library(seqinr)
library(ape)
library(magrittr)
library(vegan)
library(tibble)
library(stringr)
library(ape)
library(ggtree)
library(GUniFrac)
library(phangorn)

theme_set(theme_classic())

######################################################
#########第一种方法：RF——比较拓扑图的相似性###########
######################################################


##########################封装函数RFmeasures############################
##用法：RFmeasures(your_host_tree,your_symbiont_dendrogram,nRandom=a_random_number)

RFmeasures=function(HostTree,SymbiontDendro,nRandom=100)
{
  obs=multiRF(HostTree,as.phylo(SymbiontDendro))#phangorn
  #shuffling tips ("null model")
  
  nullSymbiontDendro=lapply(1:nRandom,shuffleTips,as.phylo(SymbiontDendro))
  
  class(nullSymbiontDendro)="multiPhylo"
  
  null=sapply(1:nRandom,FUN=function(x){multiRF(HostTree,nullSymbiontDendro[[x]])})
  
  #new (random) trees ("random model")
  tips=as.phylo(SymbiontDendro)$tip.label
  randomSymbiontDendro=rmtree(N=nRandom,n=length(tips),tip.label = tips)
  
  rd=sapply(1:nRandom,FUN=function(x){multiRF(HostTree,randomSymbiontDendro[[x]])})
  pval=sum(obs>null)/(nRandom+1)
  pvalRd=sum(obs>rd)/(nRandom+1)
  res=matrix(NA,ncol=2,nrow=2,dimnames = list(c("stat","pval"),c("RF","RFrd")))
  res["stat",]=rep(obs,2)
  
  res["pval","RF"]=pval
  
  res["pval","RFrd"]=pvalRd
  return(res)
}
##########################封装函数multiRF############################
multiRF<-function(tree1,tree2){
  trees=list(tree1,tree2)
  class(trees)="multiPhylo"
  if(class(trees)!="multiPhylo")
    stop("trees should be an object of class \"multiPhylo\"")
  N<-length(trees)
  RF<-matrix(0,N,N)
  if(any(sapply(unclass(trees),is.rooted))){
    #cat("Some trees are rooted. Unrooting all trees.\n")
    trees<-lapply(unclass(trees),unroot)
  }
  foo<-function(pp) lapply(pp,function(x,pp)
    sort(attr(pp,"labels")[x]),pp=pp)
  xx<-lapply(unclass(trees),function(x) foo(prop.part(x)))
  for(i in 1:(N-1)) for(j in (i+1):N)
    RF[i,j]<-RF[j,i]<-2*sum(!xx[[i]]%in%xx[[j]])
  RF=RF[1,2]
  RF=RF / (Nnode(tree1) + Nnode(tree2) - 2)
  RF
}


asv<- read.table("coral_V4_OTU_unrarefy.txt",header = T, row.names = 1, check.names = F,sep="\t",strip.white = T,comment.char ='')


meta <- read.delim('Endo-sample-table-coral-V4.txt',header = T, row.names = 1, check.names = F,strip.white = T,sep="\t",comment.char ='')
meta$"sample-id"<-rownames(meta)  #增加列名至最后1列


asv_tree_bray_crutis<-read.tree("coral_V4_den_grouped_bc_upgma.tre") 
asv_tree_unweighted_UniFrac<-read.tree("coral_V4_den_FI-tree_trim_uwuf_upgma.tre") 
asv_tree_weighted_UniFrac<-read.tree("coral_V4_den_FI-tree_trim_wuf_upgma.tre") 


# Format taxonomy table
tax<- read.delim('taxonomy-coral-V4.txt',header = T, row.names = 1, check.names = F,strip.white = T,sep="\t",comment.char ='',stringsAsFactors = F)
tax$Confidence <- NULL 
# Combine to make phyloseq object  


## Import host tree and format ##########
host_tree <- read.tree("coral_host_V4_Phylo.tre") 
plot(host_tree) 
# Edit tip labels 
#tip_labs <- host_tree$tip.label

# 去除叶节点标签上的引号
#host_tree$tip.label <- gsub("'", "", host_tree$tip.label,fixed = TRUE)
# 将叶节点标签中的下划线替换为空格
#host_tree$tip.label <- gsub("_", " ", host_tree$tip.label)
#host_tree
# 输出处理后的树文件
#write.tree(host_tree, "coral_host_modify.nwk")

asv_tree_bray_crutis$tip.label
asv_tree_unweighted_UniFrac$tip.label
asv_tree_weighted_UniFrac$tip.label

host_tree$tip.label


# 获取host_tree$tip.label和micro$tip.label的交集
intersection <- intersect(host_tree$tip.label, asv_tree_bray_crutis$tip.label)

asv_tree_bray_crutis<-keep.tip(asv_tree_bray_crutis,intersection)
asv_tree_unweighted_UniFrac<-keep.tip(asv_tree_unweighted_UniFrac,intersection)
asv_tree_weighted_UniFrac<-keep.tip(asv_tree_weighted_UniFrac,intersection)
host_tree<-keep.tip(host_tree,intersection)

print(is.element(host_tree$tip.label,asv_tree_bray_crutis$tip.label))####is.element(x,y)表示x元素是否在y中

print(is.binary(host_tree))#确保宿主树是二进制


##############检验宿主树和微生物树状图是否相同#################
source("PhyloSymbiosis_Functions.R")

#调用多线程计算
doParallel::registerDoParallel(90) 
RF_s_bray_crutis<-RFmeasures(host_tree,asv_tree_bray_crutis,nRandom = 9999)
RF_s_unweighted_UniFrac<-RFmeasures(host_tree,asv_tree_unweighted_UniFrac,nRandom = 9999)
RF_s_weighted_UniFrac<-RFmeasures(host_tree,asv_tree_weighted_UniFrac,nRandom = 9999)


write.table(RF_s_bray_crutis, "RF_s_bray_crutis.txt")
write.table(RF_s_unweighted_UniFrac, "RF_s_unweighted_UniFrac.txt")
write.table(RF_s_weighted_UniFrac,"RF_s_weighted_UniFrac.txt")


######################################################
#########第二种方法：Mantel test #####################
######################################################
# create distance matrix of host phylogenetic distances
host_tree_mat <- cophenetic.phylo(host_tree) 


# Create BC and UniFrac distance matrices of ASV table (collapsed by host species in qiime2)

# Import and format table
asv<-read.csv("coral_V4_OTU_unrarefy-min-2-grouped-table.tsv", sep = '\t', row.names = 1, header = T, skip=1,strip.white = T,comment.char = "",)
colnames(asv)
asv_t <- t(asv)
# Edit row names
row.names(asv_t) <- gsub("\\.", "_", row.names(asv_t))

# Check matching row names
print(is.element(row.names(asv_t), row.names(host_tree_mat)))

# align rows between dataframes
asv_t<- asv_t[match(row.names(host_tree_mat), rownames(asv_t)),] 
row.names(asv_t)
# standardise data
asv_wis <- wisconsin(asv_t) # more approriate for Bray-Curtis, but check both


asv_r <- rrarefy(asv_t, min(apply(asv_t, 1, sum))) # Rarefy more appropriate for UniFrac

# Test for Bray-Curtis correlation
bc_dist <- vegdist(asv_wis,  method = "bray")
bc_mat  <- as.matrix(bc_dist)

mantel_bray_curtis<- mantel(bc_mat, host_tree_mat, method = "pearson", permutations = 9999)

# Check stats
mantel_bray_curtis

# For UniFrac
asv_tree <- read.tree("coral_asv_rooted_V4_unrarefy.nwk") # import bacteria phylogeny

#使tree和ASV表对应上
otu_names <- colnames(asv_r)
tree_tip_labels <- asv_tree$tip.label

# Find the matching elements
matching_elements <- intersect(otu_names, tree_tip_labels)

# Subset the OTU table and tree to include only the matching elements
asv_r<- asv_r[,matching_elements]
asv_tree<- drop.tip(asv_tree, tree_tip_labels[!tree_tip_labels %in% matching_elements])


UF<- GUniFrac(asv_r, asv_tree, alpha = c(0, 0.5,1))

wUF<- UF$unifracs[, , "d_1"] # weighted
uwUF<- UF$unifracs[, , "d_UW"] # Unweighted


# Test for correlation
mantel_wUF<- mantel(wUF, host_tree_mat, method = "pearson", permutations = 9999)
mantel_uwUF <- mantel(uwUF, host_tree_mat, method = "pearson", permutations = 9999)

# Check stats
mantel_wUF
mantel_uwUF
mantel_bray_curtis
# Plot using Plot Phylosymbiosis script

output1 <- capture.output(print(mantel_wUF))
output2 <- capture.output(print(mantel_uwUF))
output3 <- capture.output(print(mantel_bray_curtis))
# 将捕获的输出写入文本文件
writeLines(output1, "mantel_wUF.txt")
writeLines(output2, "mantel_uwUF.txt")
writeLines(output3, "mantel_bray_curtis.txt")



###########利用cophylopot绘制Face to face系统发育进化树###########
library(phytools)
tree1 <- host_tree
tree2<-asv_tree_bray_crutis #bray_curtis树形图
#tree2<-asv_tree_unweighted_UniFrac
#tree2<-asv_tree_weighted_UniFrac

#creation of the association matrix:
#creation of the association matrix:
association <- cbind(tree1$tip.label, tree2$tip.label)
association
as.data.frame(association)
colnames(association)<-c("tree1$tip.label","tree2$tip.label")

cophyloplot(tree1, tree2, assoc = association, use.edge.length = FALSE,
                 length.line =0,space =50, gap = 6, font=3, lwd=0.1, lty=2) 

#保存最终的珊瑚宿主树
ggtree(tree1)+geom_tiplab()+geom_treescale(x=0,y=0)+xlim(NA,0.4)
ggsave("coral_V4_host_final.pdf")

#use.edge.length=TRUE,表示使用边长信息来绘制树形图
#length.line表示连接两个树形图节点对应关系的线的长度
#space表示两个树形图之间的水平间距
#gap表示节点对应关系线的垂直间距
plot(cophylo)
ggsave("coral_phylosymbiosis.tiff")

#plot with rotations
## Not run: 
#cophylo<- cophyloplot(tree1, tree2, assoc = association,use.edge.length=TRUE,length.line = 0, space = 500, gap = 6,rotate=TRUE)

#####设置连接线的颜色
library(randomcoloR)
palette<-distinctColorPalette(nrow(association))
pie(rep(1,nrow(association)),col=palette)
palette<-make.transparent(palette,0.5)

par(lend=3)
plot(cophylo,link.col=palette,link.lwd=4,link.type="curved",
     link.lty="solid",fsize=c(0.8,0.8))
#link.col: 用于指定连接线的颜色，palette 是一个颜色向量，每个线段使用不同的颜色。
#link.lwd: 用于指定连接线的线宽，设置为 4。
#link.type: 用于指定连接线的类型，这里使用 "curved" 表示曲线连接线。
#link.lty: 用于指定连接线的线型，这里使用 "solid" 表示实线连接线。
#fsize: 用于调整节点标签的大小，设置为 (0.8, 0.8)。
#另外，通过 par(lend = 3) 设置了线段末端的样式，其中 lend = 3 表示线段末端为圆形。


obj<-cophylo(tree1, tree2, assoc = association, print=TRUE)
plot(obj,link.type="curved",link.lwd=2,link.lty="solid",
     link.col=make.transparent("blue",0.25),fsize=0.8)
