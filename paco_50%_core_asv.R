##############format_data###############

library(ggplot2)
library(phyloseq)
library(tidyverse)
library(paco)
library(doParallel)
library(future) 
library(future.apply) 
library(multcomp) 
library(stringr) 
library(dplyr)
library(seqinr)
library(ape)
library(magrittr)
library(randomForest)
library(vegan)
library(tibble)
library(stringr)
library(ape)
library(stringr)
library(microeco)
library(file2meco)

theme_set(theme_classic())

sci_colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728",
                "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2",
                "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5", "#393b79",
                "#5254a3", "#6b6ecf", "#9c9ede", "#637939", "#8ca252", "#b5cf6b", "#cedb9c",
                "#8c6d31", "#bd9e39") 



#physeq <- meco2phyloseq(dataset_S_sw_rarefy)
#physeq

# Set core limit 
core_cutoff = 0.5     ########核心菌群的设计

# Get list of host species to subset core ASVs  
species_list<- as.character(unique(sample_data(physeq_rare)$host_species))
species_list
# Filter core microbes——50%
#Coral 
core_coral<- filter_taxa(physeq,function(x) sum(x > 0) / length(x) > core_cutoff, TRUE) 
print(core_coral) 


# Create a separate phyloseq object for each host species with core OTUs filtered 
lapply(names(core_coral), function(x) assign(x, core_coral[[x]], envir = .GlobalEnv)) # assigns a 



###########cophylogeny_PACo#############
## Import host tree and format ---- 

# Import tree 

host <- read.tree("coral_host_f.nwk") 

# Edit tip labels 
tip_labs <- host$tip.label

# 去除叶节点标签上的引号
host$tip.label <- gsub("'", "", host$tip.label,fixed = TRUE)
# 将叶节点标签中的下划线替换为空格
host$tip.label <- gsub("_", " ", host$tip.label)
host
# 输出处理后的树文件
#setwd("/home/juhuimin/202305/coral_V4/region/unrarefy")
#write.tree(host, "Cophylogeny_50%_core_asv/coral_host_modify.nwk")


###找到coral宿主物种和抽平后ASV表都有的host_species
all_host_species_keep<-intersect(host$tip.label,unique(sample_data(physeq_rare)$host_species))
all_host_species_keep
write.csv(all_host_species_keep,"coral_host_species_keep.csv")


#####host
coral_host<-keep.tip(host,all_host_species_keep)
write.tree(coral_host,"coral_tree_final.tre")

####ASV
physeq_f<- prune_samples(core_coral@sam_data$host %in% coral_host$tip.label, core_coral) 
write.tree(physeq_f@phy_tree,"physeq_f@phy_tree_final.tre")


######存在与否矩阵#########
otu <- t(as.data.frame(as.matrix(otu_table(physeq_f)))) %>% 
  apply(2, function(x) ifelse(x > 0, 1, 0)) %>% as.matrix 
otu
HP<-otu
write.table(HP,file="HP.csv",sep =",",quote = F, row.names = T)



############# Cophylogeny and specificity between cryptic coral species (Pocillopora spp.) at Mo’orea and their symbionts (Symbiodiniaceae)###########################################
### Establish PACo function: adjustment prior to procustes analysis

#Read in association data matrix
HP <- read.csv("HP.csv", header = T) #这里的HP是0和1的存在与不存在表
HP[] <- lapply(HP, function(x) type.convert(as.character(x)))
HP
sample_table <- physeq_f@sam_data
sample_table$"sample-id"<-row.names(sample_table)
# 保留列名和host_species列
sample<- subset(sample_table, select = c('sample.id', "host"))
sample<-data.frame(sample)
class(sample)
HP$sample.id<-rownames(HP)
merged_df <- merge(HP, sample, by = "sample.id")
merged_df <-merged_df[,-1]
# 使用aggregate函数重新整合数据框，并以host_species为新的列名
aggregated_df <- aggregate(. ~ host, data = merged_df, FUN = mean)
# 假设第一列为列名
rownames(aggregated_df) <- aggregated_df[, 1]

# 移除第一列
aggregated_df <- aggregated_df[, -1]


aggregated_df

HP <- as.data.frame(as.matrix(aggregated_df)) %>% 
  apply(2, function(x) ifelse(x > 0, 1, 0)) %>% as.matrix
HP<-as.data.frame(HP)
# Replace spaces with underscores in the row names of HP
rownames(HP) <- gsub(" ", "_", rownames(HP))

HP[] <- lapply(HP, function(x) type.convert(as.character(x)))
head(HP)
write.table(HP,file="HP_final.csv",sep =",",quote = F,row.names = T)



### Establish PACo function: adjustment prior to procustes analysis
PACo <- function(H.dist, P.dist, HP.bin)
{ HP.bin <- which(HP.bin > 0, arr.in=TRUE)
H.PCo <- pcoa(H.dist, correction="cailliez") #Performs PCo of Host distances 
P.PCo <- pcoa(P.dist, correction="cailliez") #Performs PCo of Parasite distances
if (is.null(H.PCo$vectors.cor)==TRUE) H.PCo <- H.PCo$vectors else
  H.PCo <- H.PCo$vectors.cor      # returns corrected pcoord 
if (is.null(P.PCo$vectors.cor)==TRUE) P.PCo <- P.PCo$vectors else
  P.PCo <- P.PCo$vectors.cor
H.PCo <- H.PCo[HP.bin[,1],]  #adjust Host PCo vectors 
P.PCo <- P.PCo[HP.bin[,2],]  #adjust Symbiont PCo vectors
list (H.PCo = H.PCo, P.PCo = P.PCo)}


#Read in association data matrix
HP <- read.csv("HP_final.csv", header =T)

rownames(HP)<-HP$host

HP$host <- NULL



### 2. DATA INPUT
#2.1 Phylogenetic trees:

#Read host tree
TreeH <- read.tree("coral_host_f.nwk")

#Read host tree
TreeP <- read.tree("phy_tree.nwk")

host.D <- cophenetic(TreeH)
para.D <- cophenetic(TreeP)

#Sort host and parasite taxa in distance matrices to match the HP matrix:
host.D <- host.D[rownames(HP), rownames(HP)]
para.D <- para.D[colnames(HP), colnames(HP)]

#D <- prepare_paco_data(host.D, para.D , HP)

#D <- add_pcoord(D,correction = "cailliez")#cailliez             D<-add_pcoord(D)

#D <- PACo(D, nperm=100000, seed=8765, method="quasiswap",symmetric=TRUE,, shuffled = TRUE)#quasiswap

#D$H
#D$P
#D$HP
#D$H_PCo
#D$P_PCo
#D$correction
#D$note
#D$proc
#D$gof
#D$method
#D$symmetric

#print(D$gof)

#D
#D <- paco_links(D)
#D$jackknife
#residuals_paco(D$proc)
#NLinks = sum(HP)

#save(D,file = "D.RData")
#z_final <- PACo(paco_data, nperm = 10000, seed = 99, symmetric = FALSE, shuffled = TRUE)
#z_links <- paco_links(z)
#residuals_paco(z_links)
#hist(z_final$shuffled, xlim=c(2, 3.2), las = 1, col = "lightblue", xlab = "Sum of squared residuals", main = '')
#abline(v = z_final$gof$ss, col = "blue", lwd = 3)








### 3. APPLY PACo FUNCTION  
PACo.fit <- PACo(host.D, para.D, HP)
HP.proc <- procrustes(PACo.fit$H.PCo, PACo.fit$P.PCo) #Procrustes Ordination 
HP.proc
NLinks = sum(HP) #Number of H-P links; needed for further computations
NLinks
save(HP.proc,file="HP.proc.RData")
HP_proc <- capture.output(print(HP.proc))
# 将捕获的输出写入文本文件
writeLines(HP_proc, "HP_proc.txt")
output <- capture.output(print(NLinks))
# 将捕获的输出写入文本文件
writeLines(output, "Nlinks.txt")


#3.1 Goodness-of-fit-test --- ## Takes a long time to run
m2.obs <- HP.proc$ss #observed sum of squares
N.perm = 100000 #set number of permutations for testing
P.value = 0
set.seed(8765) ### use this option to obtain reproducible randomizations
for (n in c(1:N.perm))
{ 
  if (NLinks <= nrow(HP) | NLinks <= ncol(HP)) 	#control statement to avoid all symbionts being associated to a single host 
  {	flag2 <- TRUE 
  while (flag2 == TRUE)	{ 
    HP.perm <- t(apply(HP,1,sample))
    if(any(colSums(HP.perm) == NLinks)) flag2 <- TRUE else flag2 <- FALSE
  }  
  } else { HP.perm <- t(apply(HP,1,sample))} #permutes each HP row independently
  PACo.perm <- PACo(host.D, para.D, HP.perm)
  m2.perm <- procrustes(PACo.perm$H.PCo, PACo.perm$P.PCo)$ss #randomized sum of squares
  if (m2.perm <= m2.obs)
  {P.value = P.value + 1} 
}
P.value <- P.value/N.perm
cat(" The observed m2 is ", m2.obs, "\n", "P-value = ", P.value, " based on ", N.perm," permutations.")

output <- capture.output(cat(" The observed m2 is ", m2.obs, "\n", "P-value = ", P.value, " based on ", N.perm," permutations."))
# 将捕获的输出写入文本文件
writeLines(output, "paco_results.txt")

#3.2 Contribution of individual links
HP.ones <- which(HP > 0, arr.in=TRUE)
SQres.jackn <- matrix(rep(NA, NLinks**2), NLinks)# empty matrix of jackknifed squared residuals
colnames (SQres.jackn) <- paste(rownames(HP.proc$X),rownames(HP.proc$Yrot), sep="-") #colnames identify the H-P link
t.critical = qt(0.975,NLinks-1) #Needed to compute 95% confidence intervals.
for(i in c(1:NLinks)) #PACo setting the ith link = 0
{HP.ind <- HP
HP.ind[HP.ones[i,1],HP.ones[i,2]]=0
PACo.ind <- PACo(host.D, para.D, HP.ind)
Proc.ind <- procrustes(PACo.ind$H.PCo, PACo.ind$P.PCo) 
res.Proc.ind <- c(residuals(Proc.ind))
res.Proc.ind <- append (res.Proc.ind, NA, after= i-1)
SQres.jackn [i, ] <- res.Proc.ind	#Append residuals to matrix of jackknifed squared residuals
} 
SQres.jackn <- SQres.jackn**2 #Jackknifed residuals are squared
SQres <- (residuals (HP.proc))**2 # Vector of original square residuals
#jackknife calculations:
SQres.jackn <- SQres.jackn*(-(NLinks-1))
SQres <- SQres*NLinks
SQres.jackn <- t(apply(SQres.jackn, 1, "+", SQres)) #apply jackknife function to matrix
phi.mean <- apply(SQres.jackn, 2, mean, na.rm = TRUE) #mean jackknife estimate per link
phi.UCI <- apply(SQres.jackn, 2, sd, na.rm = TRUE) #standard deviation of estimates
phi.UCI <- phi.mean + t.critical * phi.UCI/sqrt(NLinks) #upper 95% confidence interval


#convert vectors into dataframes
phi.UCI.df <- as.data.frame(phi.UCI)
phi.mean.df <- as.data.frame(phi.mean)
phi.UCI.df <- cbind(UCI_links = rownames(phi.UCI.df), phi.UCI.df)
rownames(phi.UCI.df) <- 1:nrow(phi.UCI.df)
phi.mean.df <- cbind(mean_links = rownames(phi.mean.df), phi.mean.df)
rownames(phi.mean.df) <- 1:nrow(phi.mean.df)

#Merge dataframes
phi.UCI.mean = merge(phi.UCI.df, phi.mean.df, by.x=c("UCI_links"), by.y=c("mean_links"))
save(phi.UCI.mean,file="phi.UCI.mean.RData")
write.table(phi.UCI.mean,file="phi.UCI.mean.txt",sep ="\t",quote = F, row.names = T)



#Add new column for host
phi.UCI.mean$host <- phi.UCI.mean$UCI_links
phi.UCI.mean$host <- gsub("\\-.*","",phi.UCI.mean$host)


#Add new column for ASV clades
phi.UCI.mean$ASV <- gsub(".*-","",phi.UCI.mean$UCI_links)

#Find mean of phi.mean
median(phi.UCI.mean$phi.mean) # median = 1.195144e-05
phi.UCI.mean


phi.UCI.mean$host  <- gsub("_", " ", phi.UCI.mean$host )

phi.UCI.mean$host <- factor(phi.UCI.mean$host ,levels=c("Dipsastraea sp2.", "Dipsastraea speciosa", "Favites halicora", "Hydnophora bonsai",
                                                        "Hydnophora exesa", "Leptastrea transversa", "Lobophyllia hemprichii", "Pavona decussata",
                                                        "Pavona duerdeni", "Pavona frondifera", "Platygyra daedalea", "Platygyra sp1.", "Plesiastrea versipora",
                                                        "Porites lutea"))

phi.UCI.mean

p1<-ggplot(phi.UCI.mean, aes(x=reorder(UCI_links, phi.mean), y=phi.mean, fill=host))+
  geom_bar(stat='identity')+
  geom_hline(yintercept=median(phi.UCI.mean$phi.mean), linetype='dashed', col = 'red')+
  theme_minimal(base_size=16)+
  geom_errorbar(aes(ymin=phi.mean, ymax=phi.UCI), width=.1)+ 
  scale_fill_manual(name="host",
                    values =sci_colors)+
  coord_flip()+
  theme(axis.text=element_text(size=2))+
  labs(title="",
       x ="ASV - host", y = "Squared residuals")
p1
ggsave(p1,file="Squared residuals.pdf")



