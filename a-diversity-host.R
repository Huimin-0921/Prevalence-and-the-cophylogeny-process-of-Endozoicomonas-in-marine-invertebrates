load("dataset5.RData")

sci_colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728",
                "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2",
                "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5", "#393b79",
                "#5254a3", "#6b6ecf", "#9c9ede", "#637939", "#8ca252", "#b5cf6b", "#cedb9c",
                "#8c6d31", "#bd9e39") 
dataset5$sample_sums() %>% range
dataset5$rarefy_samples(sample.size = 100)
dataset5$tidy_dataset()
dataset5$cal_alphadiv(PD = FALSE)

######################host_phylum#########################
t1 <- trans_alpha$new(dataset = dataset5, group = "host_phylum")# return t1$alpha_stat
t1$data_stat
t1$data_alpha
write.table(t1$data_stat, file="alpha_diversity_dataset5_host_phylum/data_stat.csv",sep =",",quote = F, row.names = T)
write.table(t1$data_alpha, file="alpha_diversity_dataset5_host_phylum/data_alpha.csv",sep =",",quote = F, row.names = T)

#使用KW_dunn检验
t1$cal_diff(method = "KW_dunn",measure="Shannon")
# return t1$res_alpha_diff
t1$res_diff
write.table(t1$res_diff, file="alpha_diversity_dataset5_host_phylum/α-diff_Shannon.csv",sep =",",quote = F, row.names = T)
p1<-t1$plot_alpha(measure = "Shannon", add_sig_text_size = 6, boxplot_add = "jitter", order_x_mean = TRUE,color_values =sci_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(p1,file='alpha_diversity_dataset5_host_phylum/Shannon_host_phylum.pdf',width=8,height = 6)
save(p1,file="alpha_diversity_dataset5_host_phylum/p1.RData")

######################host_class#########################
load("dataset6.RData")

sci_colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728",
                "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2",
                "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5", "#393b79",
                "#5254a3", "#6b6ecf", "#9c9ede", "#637939", "#8ca252", "#b5cf6b", "#cedb9c",
                "#8c6d31", "#bd9e39") 
dataset6$sample_sums() %>% range
dataset6$rarefy_samples(sample.size = 100)
dataset6$tidy_dataset()
dataset6$cal_alphadiv(PD = FALSE)
t1 <- trans_alpha$new(dataset = dataset6, group = "host_class")# return t1$alpha_stat
t1$data_stat
t1$data_alpha
write.table(t1$data_stat, file="alpha_diversity_dataset6_host_class/data_stat.csv",sep =",",quote = F, row.names = T)
write.table(t1$data_alpha, file="alpha_diversity_dataset6_host_class/data_alpha.csv",sep =",",quote = F, row.names = T)

#使用KW_dunn检验
t1$cal_diff(method = "KW_dunn",measure="Shannon")
# return t1$res_alpha_diff
t1$res_diff
write.table(t1$res_diff, file="alpha_diversity_dataset6_host_class/α-diff_Shannon.csv",sep =",",quote = F, row.names = T)
p2<-t1$plot_alpha(measure = "Shannon", add_sig_text_size = 6, boxplot_add = "jitter", order_x_mean = TRUE,color_values =sci_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(p2,file='alpha_diversity_dataset6_host_class/Shannon_host_class.pdf',width=20,height = 16)
save(p2,file="alpha_diversity_dataset6_host_class/p2.RData")





######################host_order#########################
load("dataset16.RData")

sci_colors <- c("#1f77b4", "#aec7e8", "#ff7f0e", "#ffbb78", "#2ca02c", "#98df8a", "#d62728",
                "#ff9896", "#9467bd", "#c5b0d5", "#8c564b", "#c49c94", "#e377c2", "#f7b6d2",
                "#7f7f7f", "#c7c7c7", "#bcbd22", "#dbdb8d", "#17becf", "#9edae5", "#393b79",
                "#5254a3", "#6b6ecf", "#9c9ede", "#637939", "#8ca252", "#b5cf6b", "#cedb9c",
                "#8c6d31", "#bd9e39") 
dataset16$sample_sums() %>% range
dataset16$rarefy_samples(sample.size = 100)
dataset16$tidy_dataset()
dataset16$cal_alphadiv(PD = FALSE)
t1 <- trans_alpha$new(dataset = dataset16, group = "host_order")# return t1$alpha_stat
t1$data_stat
t1$data_alpha
write.table(t1$data_stat, file="alpha_diversity_dataset16_host_order/data_stat.csv",sep =",",quote = F, row.names = T)
write.table(t1$data_alpha, file="alpha_diversity_dataset16_host_order/data_alpha.csv",sep =",",quote = F, row.names = T)

#使用KW_dunn检验
t1$cal_diff(method = "KW_dunn",measure="Shannon")
# return t1$res_alpha_diff
t1$res_diff
write.table(t1$res_diff, file="alpha_diversity_dataset16_host_order/α-diff_Shannon.csv",sep =",",quote = F, row.names = T)
p3<-t1$plot_alpha(measure = "Shannon", add_sig_text_size = 6, boxplot_add = "jitter", order_x_mean = TRUE,color_values =sci_colors)+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(p3,file='alpha_diversity_dataset16_host_order/Shannon_host_order.pdf',width=28,height = 16)
save(p3,file="alpha_diversity_dataset16_host_order/p3.RData")
