library(paco)

sci_colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728",
                "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2",
                "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5", "#393b79",
                "#5254a3", "#6b6ecf", "#9c9ede", "#637939", "#8ca252", "#b5cf6b", "#cedb9c",
                "#8c6d31", "#bd9e39","#E6550D","#D74E9C","#5AAA35","#5BB6F0","#F849B6","#F38BB6",
                "#82ABF2","#D085C2","#13BC97","#6DF4F6","#FEBE4C","#B884EF","#AEDF14","#F0758A",
                "#6276D9","#4CEE24")



# residuals
res = residuals_paco(D$proc) %>% as.data.frame 
colnames(res) = 'residuals'
res$"link"<-rownames(res)
#Formatting output
jackknife<-D$jackknife%>% as.data.frame 
colnames(jackknife) = 'jackknife'
jackknife$"link"<-rownames(jackknife)
phi.UCI.mean<-merge(res,jackknife,by="link")
phi.UCI.mean = phi.UCI.mean %>%
  mutate(comparison =phi.UCI.mean$link ) %>%
  separate(comparison, c('host', 'ASV'), sep='-') 
phi.UCI.mean


#Find mean of phi.mean
median(phi.UCI.mean$residuals) # median = 0.08116533
phi.UCI.mean


phi.UCI.mean$host  <- gsub("_", " ", phi.UCI.mean$host )

phi.UCI.mean$host <- factor(phi.UCI.mean$host ,levels=c(
  "Acropora cervicornis", "Acropora cytherea", "Acropora elseyi", "Acropora gemmifera",
  "Acropora hyacinthus", "Acropora intermedia", "Acropora latistella", "Acropora millepora",
  "Acropora nasuta", "Agaricia undata", "Astrea curta", "Ctenactis crassa",
  "Cyphastrea microphthalma", "Cyphastrea serailia", "Echinopora horrida", "Echinopora lamellosa",
  "Favites abdita", "Galaxea astreata", "Goniastrea edwardsi", "Goniastrea favulus",
  "Herpolitha limax", "Hydnophora exesa", "Isopora palifera", "Leptastrea transversa",
  "Lobophyllia hemprichii", "Merulina ampliata", "Montipora aequituberculata", "Montipora capitata",
  "Montipora efflorescens", "Orbicella annularis", "Pavona cactus", "Pavona varians",
  "Platygyra daedalea", "Platygyra sinensis", "Pocillopora verrucosa", "Porites astreoides",
  "Porites compressa", "Porites cylindrica", "Porites lobata", "Porites porites",
  "Psammocora digitata", "Seriatopora caliendrum", "Seriatopora hystrix"
))

phi.UCI.mean

save(phi.UCI.mean,file="phi.UCI.mean.RData")

p1<-ggplot(phi.UCI.mean, aes(x=reorder(link, residuals), y=residuals, fill=host))+
  geom_bar(stat='identity')+
  geom_hline(yintercept=median(phi.UCI.mean$residuals), linetype='dashed', col = 'red')+
  theme_minimal(base_size=16)+
  geom_errorbar(aes(ymin=residuals, ymax=jackknife), width=.1)+ 
  scale_fill_manual(name="host",
                    values =sci_colors)+
  coord_flip()+
  theme(axis.text=element_text(size=4))+
  labs(title="",
       x ="host - ASV", y = "Squared residuals")
p1
ggsave(p1,file="Squared residuals_图例.pdf",width = 8, height =20,units = "cm")

p2<-p1+ theme(legend.position = "none")
ggsave(p2,file="Squared residuals.pdf",width = 8, height =20,units = "cm")



#箱型图
phi.UCI.mean
load("dataset.RData")
taxonomy<-dataset$tax_table
taxonomy$ASV<-row.names(taxonomy)
res<-merge(phi.UCI.mean,taxonomy,by="ASV")
unique(res$host)
meta<-dataset_S_sw_rarefy$sample_table
sort(unique(meta$host))
res<-merge(res,meta,by="host")
res
write.table(res,file="phi.UCI.mean.txt",sep ="\t",quote = F, row.names = T)


p <- ggplot(res, aes(host, residuals)) +
  geom_boxplot() +
  labs(x='host', y='Residuals') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1)
  )
p 
ggsave(p,file="p.pdf")

res$fourth <-(res$residuals ^ (1/4)) # Transform data
resids.lm <- lm(fourth ~ASV , res)
plot(resids.lm)

summary(resids.lm) # summarise linear model.
anova(resids.lm) # for anova
output <- capture.output(anova(resids.lm))
# 将捕获的输出写入文本文件
writeLines(output, "anova(resids.lm)_res_ASV.txt")


# Anova on residuals (host group)
res$fourth <-(res$residuals ^ (1/4)) # Transform data
resids.lm <- lm(fourth ~host , res)
plot(resids.lm)

summary(resids.lm) # summarise linear model.
anova(resids.lm) # for anova
output <- capture.output(anova(resids.lm))
# 将捕获的输出写入文本文件
writeLines(output, "anova(resids.lm)_res_host.txt")


host_link_status<-as.data.frame(summary(res$host))
host_link_status <- host_link_status %>%
  rename(host_link = `summary(res$host)`)
host_link_status$host<-row.names(host_link_status)
write.table(host_link_status,file="host_link_status.txt",sep ="\t",quote = F, row.names = T)


# Create a bar plot
bar_plot <- ggplot(host_link_status, aes(x = host, y = host_link, fill = host)) +
  geom_bar(stat = "identity") +
  labs(x = 'Host', y = 'Host Link') +
  scale_fill_manual(values = sci_colors) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  scale_y_continuous(breaks = seq(0, 9, by = 1))

# Save or display the plot
ggsave("bar_plot_host_link_status_图例.pdf", bar_plot, width = 8, height = 6)
bar_plot<-bar_plot+ theme(legend.position = "none")
ggsave("bar_plot_host_link_status.pdf", bar_plot, width = 8, height = 6)



#post hoc test
glht.resids.sum <- summary(glht(resids.lm, linfct = mcp(host = "Tukey")), test = adjusted("bonferroni"))#这一行进行多重比较校正的事后检验
resids.groups <- cld(glht.resids.sum)
resids.groups.df <- fortify(resids.groups)#将 resids.groups 转化为数据框，存储在 resids.groups.df 中。后续你可以查看该数据框的结构
colnames(resids.groups.df) <- c("host", "letters")

ymax <- tapply(res$residuals, res$host, max)
resids.groups.df$Ymax <- ymax # add to plot above
resids.groups.df

# plot by host species
phylum <- ggplot(res) +
  geom_boxplot(aes(host, residuals, fill = host), outlier.shape=NA) +
  geom_jitter(aes(host, residuals,colour = host), position=position_jitter(width=.002, height=0)) +#
  labs(x='Species', y='Residuals') +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle=90, hjust=1)
  )+scale_color_manual(values = sci_colors)+
  scale_fill_manual(values = sci_colors) 
phylum
ggsave(phylum,file="phylum.pdf", width = 16, height = 6)

