library(microeco)
library(magrittr)
library(ggplot2)
library(treeio)
require(ggtree)
library(mecodev)
require(tidyr)
require(dplyr)
require(tidytree)

sci_colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728",
                "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2",
                "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5", "#393b79",
                "#5254a3", "#6b6ecf", "#9c9ede", "#637939", "#8ca252", "#b5cf6b", "#cedb9c",
                "#8c6d31", "#bd9e39") 


# 定义要生成的编号数量
seed(123)
num_ids <- 60

# 生成随机的 RGB 颜色
generate_random_colors <- function(n) {
  sample_colors <- sample(0:255, n * 3, replace = TRUE)
  matrix(sample_colors, ncol = 3)
}

# 生成颜色差异大的 RGB 颜色
colors <- generate_random_colors(num_ids)

# 将 RGB 颜色转换为十六进制颜色代码
rgb_to_hex <- function(r, g, b) {
  sprintf("#%02X%02X%02X", r, g, b)
}

# 将颜色矩阵转换为十六进制颜色代码
color_codes <- apply(colors, 1, function(color) {
  rgb_to_hex(color[1], color[2], color[3])
})

# 打印生成的颜色代码
print(color_codes)
random_colors <- color_codes #随机产生50个颜色编号
###########################dataset_rarefy#################
load("dataset.RData")
t1 <- trans_rarefy$new(dataset, alphadiv = "Shannon", depth = c(0,50,100,200,500,1000,2000,4000,6000,8000,10000))
save(t1,file="rarefy/dataset_rarefy.RData")
t1$plot_rarefy(color = "host_type", show_point = FALSE, add_fitting = FALSE,color_values =sci_colors)
t1$plot_rarefy(color = "host_type", show_point = FALSE, add_fitting = TRUE,color_values =sci_colors)
p<-t1$plot_rarefy(color = "Ecosystem", show_point = FALSE)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(p,file='rarefy/dataset.pdf',width=8,height = 6)


###########################dataset5_rarefy#################
load("dataset5.RData")
t1 <- trans_rarefy$new(dataset5, alphadiv = "Shannon", depth = c(0,500,1000,2000,4000,6000,8000,10000))
save(t1,file="rarefy/dataset5_rarefy.RData")
p<-t1$plot_rarefy(color = "host_phylum", color_values =sci_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
ggsave(p,file='rarefy/dataset5.pdf',width=12,height = 8)


t1 <- trans_rarefy$new(dataset5, alphadiv = "Shannon", depth = c(0,50,100,200,500,800,1000,2000))
save(t1,file="rarefy/dataset5_rarefy_1.RData")
p<-t1$plot_rarefy(color = "host_phylum", color_values =sci_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
ggsave(p,file='rarefy/dataset5_1.pdf',width=12,height = 8)

t1 <- trans_rarefy$new(dataset5, alphadiv = "Shannon", depth = c(0,50,100,200,500))
save(t1,file="rarefy/dataset5_rarefy_2.RData")
p<-t1$plot_rarefy(color = "host_phylum", color_values =sci_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
print(p)
ggsave(p,file='rarefy/dataset5_2.pdf',width=12,height = 8)





#########################dataset17_rarefy########################
load("dataset17.RData")
t1 <- trans_rarefy$new(dataset17, alphadiv = "Shannon", depth = c(0,50,100,200,400,500,800,1000,2000))
p<-t1$plot_rarefy(color = "Ecosystem")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
save(t1,file="rarefy/dataset17_rarefy.RData")
ggsave(p,file='rarefy/dataset17.pdf',width=8,height = 6)

t1 <- trans_rarefy$new(dataset17, alphadiv = "Shannon", depth = c(0,50,100,200,400,500))
p<-t1$plot_rarefy(color = "Ecosystem")+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
save(t1,file="rarefy/dataset17_rarefy1.RData")
ggsave(p,file='rarefy/dataset17_1.pdf',width=8,height = 6)


#########################dataset11_rarefy########################
load("dataset11.RData")
t1 <- trans_rarefy$new(dataset11, alphadiv = "Shannon", depth = c(0,50,100,200,400,500,800,1000,2000))
p<-t1$plot_rarefy(color = "host_family",color_values =random_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
save(t1,file="rarefy/dataset11_rarefy_family.RData")
ggsave(p,file='rarefy/dataset11_family.pdf',width=8,height = 6)

t1 <- trans_rarefy$new(dataset11, alphadiv = "Shannon", depth = c(0,50,100,200,300,400,500,800))
p<-t1$plot_rarefy(color = "host_family",color_values =random_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
save(t1,file="rarefy/dataset11_rarefy_family_1.RData")
ggsave(p,file='rarefy/dataset11_family_1.pdf',width=8,height = 6)





load("dataset11.RData")
t1 <- trans_rarefy$new(dataset11, alphadiv = "Shannon", depth = c(0,50,100,200,400,500,800,1000,2000))
p<-t1$plot_rarefy(color = "host_order",color_values =sci_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
save(t1,file="rarefy/dataset11_rarefy_order.RData")
ggsave(p,file='rarefy/dataset11_order.pdf',width=8,height = 6)

t1 <- trans_rarefy$new(dataset11, alphadiv = "Shannon", depth = c(0,50,100,200,300,400,500,800))
p<-t1$plot_rarefy(color = "host_order",color_values =sci_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
save(t1,file="rarefy/dataset11_rarefy_order_1.RData")
ggsave(p,file='rarefy/dataset11_order_1.pdf',width=8,height = 6)




#########################dataset12_rarefy########################
load("dataset12.RData")
t1 <- trans_rarefy$new(dataset12, alphadiv = "Shannon", depth = c(0,50,100,200,300,400,500,600,800))
p<-t1$plot_rarefy(color = "host_family",color_values =random_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
save(t1,file="rarefy/dataset12_rarefy_family.RData")
ggsave(p,file='rarefy/dataset12_family.pdf',width=8,height = 6)


load("dataset12.RData")
t1 <- trans_rarefy$new(dataset12, alphadiv = "Shannon", depth = c(0,50,100,200,300,400,500,600,800))
p<-t1$plot_rarefy(color = "host_order",color_values =sci_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
save(t1,file="rarefy/dataset12_rarefy_order.RData")
ggsave(p,file='rarefy/dataset12_order.pdf',width=8,height = 6)


