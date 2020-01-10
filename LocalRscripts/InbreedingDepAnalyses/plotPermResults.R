#set libraries
library(tidyverse)
library(cowplot)

#read data in
set.seed(4000)
dfPlink = read.table("~/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/allPerms_plink.txt", col.names = c("count")) %>%
  mutate(method= "plink") %>%
  sample_n(4000) 

dfVcf = read.table("~/DogProject_Jaz/LocalRscripts/InbreedingDepAnalyses/allPerms.txt", col.names = c("count")) %>%
  mutate(method= "vcftools") 

mergeDF = rbind.data.frame(dfPlink, dfVcf)

#plot 
inset = ggplot(dfVcf, aes(x=count)) +
  geom_density(fill="#3288BD") + 
  labs(y = "Density", x=paste("Number of genes", "without ROH", sep = "\n")) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 14), 
                   axis.text.y = element_text(size = 14), 
                   plot.title=element_text(size=14, face = "bold", hjust=0.5), 
                   axis.title=element_text(size=16),
                   legend.title=element_text(size=16), 
                   legend.text=element_text(size=14))

vcfTools = ggplot(dfVcf, aes(x=count)) +
  geom_density(fill="#3288BD") + 
  geom_vline(xintercept = 27, col="#D53E4F") +
  labs(y = "", x="", title = "vcftools") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title=element_text(size=14, face = "bold", hjust=0.5), 
        axis.title=element_text(size=18),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16))

plot1 = ggdraw() +
  draw_plot(vcfTools) +
  draw_plot(inset, x = 0.6, y = .64, width = .3, height = .3)

##plots plink

inset2 = ggplot(dfPlink, aes(x=count)) +
  geom_density(fill="#FDAE61") + 
  labs(y = "Density", x=paste("Number of genes", "without ROH", sep = "\n")) + 
  theme_bw() +
  theme(axis.text.x = element_text(size = 14), 
        axis.text.y = element_text(size = 14), 
        plot.title=element_text(size=14, face = "bold", hjust=0.5), 
        axis.title=element_text(size=16),
        legend.title=element_text(size=16), 
        legend.text=element_text(size=14))

plink = ggplot(dfPlink, aes(x=count)) +
  geom_density(fill="#FDAE61") + 
  geom_vline(xintercept = 76, col="#D53E4F") +
  labs(y = "", x="", title = "plink") + 
  theme_bw() + 
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), 
        plot.title=element_text(size=14, face = "bold", hjust=0.5), 
        axis.title=element_text(size=18),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=16))

plot2 = ggdraw() +
  draw_plot(plink) +
  draw_plot(inset2, x = 0.6, y = .64, width = .3, height = .3)


#put it all together
figure = ggarrange(plot1, plot2, nrow = 1)
annotate_figure(figure,
                bottom = text_grob("Number of genes without ROH", size = 20,face = "bold"),
                left = text_grob("Density", size = 20, rot = 90, face = "bold")
)
