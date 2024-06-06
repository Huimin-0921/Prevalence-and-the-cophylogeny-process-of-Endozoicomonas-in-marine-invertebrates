library(ggplot2)
library(microeco)
library(dplyr)
library(ggpubr)
theme_set(theme_classic())

#load("dataset.RData")
#dataset13<-clone(dataset)
#dataset13$sample_table<-subset(dataset13$sample_table,host_phylum == "Cnidaria")
#dataset13$tidy_dataset()
#dataset13$rarefy_samples(sample.size = 1000)
#dataset13$cal_alphadiv(PD = FALSE)
#dataset13$save_alphadiv(dirpath = "dataset13_1000reads_alpha_diversity")
#dataset13$cal_abund()
#dataset13$save_abund(dirpath = "dataset13_1000reads_taxa_abund")
#metadata<-dataset$sample_table

########################coral######################################
load("dataset13.RData")

a1<- read.table('coral-metadata_a-diversity-mya.txt',sep="\t", header=TRUE,check.names = F)
df1 <- a1[, c("sample-id","lat", "Shannon"), drop = FALSE]
df1 <- df1[!(is.na(df1$lat) | df1$lat == "unknown"), ]
#删除df的列为unknown

cor1<-cor.test(abs(as.numeric(df1$lat)), as.numeric(as.character(df1$Shannon)),
         alternative = "two.sided", method = "pearson", conf.level = 0.95,group=host_order)

plot1 <-ggplot(df1, aes(x=abs(as.numeric(as.character(lat))), y=as.numeric(as.character(Shannon))#, color=host_class
)) +
    geom_point(shape=16,colour="#E7B800",fill="#E7B800",size=3,stroke= .5,alpha=0.25) +
    geom_smooth(method=lm,   # Add linear regression lines
                se=TRUE,linetype=2,color="black",fill="#808080") +
  
    xlab("latitude")+#scale_color_manual(values = c("#ff7400","#009999","#3914AF")) +  # 自定义颜色映射
    ylab("Alpha diversity\nShannon")+
    stat_cor(method='pearson',aes(x=abs(as.numeric(lat)), y=as.numeric(Shannon)))+
  scale_x_continuous(limits = c(0, 70), breaks = c(0, 10, 20, 30, 40, 50,60,70)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3, 4)) 
  
plot1 
save(plot1,file="data/Shannon_lat/plot1.RData")
ggsave(plot1,file='data/Shannon_lat/Cnidaria.pdf',width=8,height = 4)
ggsave(plot1,file='data/Shannon_lat/Cnidaria.tiff',width=8,height = 4,dpi=300)


########################sponge######################################
load("dataset14.RData")

a2<- read.table('sponge-metadata_a-diversity-mya.txt',sep="\t", header=TRUE,check.names = F)
df2 <- a2[, c("sample-id","lat", "Shannon"), drop = FALSE]
df2 <- df2[!(is.na(df2$lat) | df2$lat == "unknown"), ]
#删除df的列为unknown

cor2<-cor.test(abs(as.numeric(df2$lat)), as.numeric(as.character(df2$Shannon)),
         alternative = "two.sided", method = "pearson", conf.level = 0.95,group=host_order)

plot2 <-ggplot(df2, aes(x=abs(as.numeric(as.character(lat))), y=as.numeric(as.character(Shannon))#, color=host_class
)) +
    geom_point(shape=16,colour="#42b983",fill="#42b983",size=3,stroke= .5,alpha=0.25) +
    geom_smooth(method=lm,   # Add linear regression lines
                se=TRUE,linetype=2,color="black",fill="#808080") +##6495ED
  
    xlab("Divergence time (lat)")+#scale_color_manual(values = c("#ff7400","#009999","#3914AF")) +  # 自定义颜色映射
    ylab("Alpha diversity\nShannon")+
    stat_cor(method='pearson',aes(x=abs(as.numeric(lat)), y=as.numeric(Shannon)))+
  scale_x_continuous(limits = c(0, 70), breaks = c(0, 10, 20, 30, 40, 50,60,70)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3, 4)) 
plot2 
save(plot2,file="data/Shannon_lat/plot2.RData")
ggsave(plot2,file='data/Shannon_lat/Porifera.pdf',width=8,height = 4)
ggsave(plot2,file='data/Shannon_lat/Porifera.tiff',width=8,height = 4,dpi=300)

p<-plot1+plot2
library(gridExtra)
p <- grid.arrange(plot1, plot2, ncol = 1)

df1$Group <- "Cnidaria"
df2$Group<-"Porifera"
data<-rbind(df1,df2)


p4 <- ggplot(data, aes(y = Shannon, x = abs(data$lat),color=Group)) + 
  geom_point(shape = 16, aes(color = Group, fill = Group), size = 3, stroke = 1, alpha = 0.25) + 
  geom_smooth(method = lm, se = TRUE, color = "black",fill="#808080",linetype = 2,aes(group=Group)) + 
  labs(x = "latitude", y ="Alpha diversity\nShannon") + 
  stat_cor(method = 'pearson', aes(group = Group, x = lat, y = Shannon), color = "black") +
  scale_x_continuous(limits = c(0, 80), breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80)) +
  scale_y_continuous(limits = c(0, 4), breaks = c(0, 1, 2, 3, 4)) +
  scale_color_manual(values = c("#E7B800", "#42b983")) +
  scale_fill_manual(values = c("#E7B800", "#42b983")) 

p4

save(p4,file="data/Shannon_lat/p4.RData")

ggsave(p4,file='data/Shannon_lat/Cnidaria_Porifera1.pdf',width=8,height = 4,dpi=300)
ggsave(p4,file='data/Shannon_lat/Cnidaria_Porifera.tiff',width=8,height = 4,dpi=300)

