###############Ecosystem###########

library(UpSetR)
load("dataset.RData")
# 筛选出总reads数大于等于100的ASV
dataset$otu_table <- dataset$otu_table[rowSums(dataset$otu_table) >= 100, ]

dataset$tidy_dataset()

tmp <- dataset$merge_samples(use_group = "Ecosystem")
tmp

t1 <- trans_venn$new(dataset = tmp)
t1$data_details
t1$data_summary
write.table(t1$data_details,file="Ecosystem/Ecosystem_data_details.txt",sep ="\t",quote = F, row.names = T)
write.table(t1$data_summary,file="Ecosystem/Ecosystem_data_summary.txt",sep ="\t",quote = F, row.names = T)

df<-tmp$otu_table
# 使用逻辑运算将大于 0 的值替换为 1
df[df > 0] <- 1
df
save(df,file="Ecosystem/df.RData")


pdf("Ecosystem/Ecosystem_upset.pdf", width=4,height = 3)
upset(
  df,
  #nsets = 5, #显示集合的数量
  nintersects = NA, #要绘制的最多交集数量，NA则全部绘制
  sets = c( "Symbiotic","non-reef-sediment","non-reef-seawater","coral reef-sediment", "coral reef-seawater"),
  keep.order = TRUE,
  order.by = "freq",
  decreasing = TRUE,
  set_size.show = T, #显示每个集合的总元素数(在左侧条形图中)
  main.bar.color = '#2CA02CFF', #上方柱形图颜色#FFBB78FF
  sets.bar.color = c("#66a61e", "#9467BD", "#e6ab02", "#e377c2", "#729ece"), #左方条形图颜色
  matrix.color = '#2CA02CFF', #交点颜色#2CA02CFF
  point.size = 3.3, #交点大小
  line.size = 0.8, #交点连线粗细
  shade.color = 'grey', #交集矩阵中阴影行的颜色
  shade.alpha = 0.2, #阴影行的不透明度
  matrix.dot.alpha = 0.7, #交集矩阵中空交点的不透明度
  mb.ratio = c(0.7, 0.3),#上方柱形图和下方交集矩阵的占比
  set_size.scale_max=6000)#设置左边柱状图的最大刻度值)
dev.off()




###########host_phylum##############
dataset5<-clone(dataset)
dataset5$sample_table <- subset(dataset5$sample_table,Ecosystem =="Symbiotic")
dataset5$sample_table<-dataset5$sample_table[dataset5$sample_table$host_phylum != 0, ]
dataset5$tidy_dataset()

tmp <- dataset5$merge_samples(use_group = "host_phylum")
tmp

t1 <- trans_venn$new(dataset = tmp)
t1$data_details
t1$data_summary
write.table(t1$data_details,file="host_phylum/host_phylum_data_details.txt",sep ="\t",quote = F, row.names = T)
write.table(t1$data_summary,file="host_phylum/host_phylum_data_summary.txt",sep ="\t",quote = F, row.names = T)

df<-tmp$otu_table
# 假设 df 是你的数据框
# 使用逻辑运算将大于 0 的值替换为 1
df[df > 0] <- 1
df
save(df,file="host_phylum/df.RData")


pdf("host_phylum/host_phylum_upset.pdf",width=,height = 6)
upset(
  df,
  #nsets = 5, #显示集合的数量
  nintersects = NA, #要绘制的最多交集数量，NA则全部绘制
  sets = c("Annelida","Arthropoda", "Chordata","Cnidaria","Echinodermata","Mollusca","Platyhelminthes","Porifera"),
  keep.order = TRUE,
  order.by = "freq",
  decreasing = TRUE,
  set_size.show = T, #显示每个集合的总元素数(在左侧条形图中)
  main.bar.color = '#2CA02CFF', #上方柱形图颜色#FFBB78FF
  sets.bar.color = c("#1F78B4", "#729ECE","#FF7F0E","#E6AB02","#1B9E77","#66A61E","#9467BD","#E377C2"), #左方条形图颜色
  matrix.color = '#2CA02CFF', #交点颜色#2CA02CFF
  point.size = 3.3, #交点大小
  line.size = 0.8, #交点连线粗细
  shade.color = 'grey', #交集矩阵中阴影行的颜色
  shade.alpha = 0.2, #阴影行的不透明度
  matrix.dot.alpha = 0.7, #交集矩阵中空交点的不透明度
  mb.ratio = c(0.7, 0.3),#上方柱形图和下方交集矩阵的占比
  set_size.scale_max=5000)
dev.off()



###########host_class##############
dataset6<-clone(dataset)
dataset6$sample_table <- subset(dataset6$sample_table,Ecosystem =="Symbiotic")
dataset6$sample_table<-dataset6$sample_table[dataset6$sample_table$host_class != 0, ]
dataset6$tidy_dataset()
tmp <- dataset6$merge_samples(use_group = "host_class")
tmp

t1 <- trans_venn$new(dataset = tmp)
t1$data_details
t1$data_summary
write.table(t1$data_details,file="host_class/host_class_data_details.txt",sep ="\t",quote = F, row.names = T)
write.table(t1$data_summary,file="host_class/host_class_data_summary.txt",sep ="\t",quote = F, row.names = T)

df<-tmp$otu_table
# 假设 df 是你的数据框
# 使用逻辑运算将大于 0 的值替换为 1
df[df > 0] <- 1
df
save(df,file="host_class/df.RData")

sets<-c("Anthozoa","Polychaeta","Cephalopoda","Demospongiae","Ascidiacea","Echinoidea","Asteroidea","Thaliacea", "Homoscleromorpha","Bivalvia", 
        "Calcarea","Scyphozoa","Gastropoda","Holothuroidea","Malacostraca","Ophiuroidea","Pelmatozoa","Rhabditophora","Thecostraca","Hydrozoa",        
        "Hexactinellida")
sci_colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728",
                "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2",
                "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf") 


pdf("upset/host_class/host_class_upset.pdf",width=28,height =20)
upset(
  df,
  #nsets = 5, #显示集合的数量
  nintersects = NA, #要绘制的最多交集数量，NA则全部绘制
  sets<-c("Anthozoa", "Ascidiacea", "Asteroidea", "Bivalvia", "Calcarea", "Cephalopoda", "Demospongiae",
          "Echinoidea", "Gastropoda" , "Holothuroidea", "Homoscleromorpha", "Hydrozoa", "Malacostraca",
          "Ophiuroidea", "Pelmatozoa", "Polychaeta", "Rhabditophora","Scyphozoa","Thecostraca"),
  keep.order = TRUE,
  order.by = "freq",
  decreasing = TRUE,
  set_size.show = T, #显示每个集合的总元素数(在左侧条形图中)
  main.bar.color = '#2CA02CFF', #上方柱形图颜色#FFBB78FF
  sets.bar.color = sci_colors, #左方条形图颜色
  matrix.color = '#2CA02CFF', #交点颜色#2CA02CFF
  point.size = 3.3, #交点大小
  line.size = 0.8, #交点连线粗细
  shade.color = 'grey', #交集矩阵中阴影行的颜色
  shade.alpha = 0.2, #阴影行的不透明度
  matrix.dot.alpha = 0.7, #交集矩阵中空交点的不透明度
  mb.ratio = c(0.7, 0.3),#上方柱形图和下方交集矩阵的占比
  set_size.scale_max=6000)
dev.off


