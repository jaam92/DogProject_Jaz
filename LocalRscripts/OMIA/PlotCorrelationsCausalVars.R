#####Load Libraries
source("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/R_rainclouds.R")
source("~/Documents/DogProject_Jaz/LocalRscripts/OMIA/SummarizeData.R")
library(ggpubr)
library(ggrepel)
library(randomcoloR)

######Plot Linear Regression Function###
ggplotRegression = function (fit) {
  ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) + 
    geom_point(size = 2) + 
    stat_smooth( method = 'lm', col = "blue") +  
    theme_bw() + 
    labs(title = bquote(R^2== ~.(signif(summary(fit)$adj.r.squared, 5))~"&"~"p"==~.(signif(summary(fit)$coef[2,4], 5))))
  
}
######Plot Causal Vars Fxn######
######Labelling the top 5% of Scores and ~1% of Causal Var Counts
plotCausal = function(dataFrame, scoreCutOff, xaxisLabel, indicator){
  if (indicator == "ROH") {
    newDF = dataFrame %>% 
      rename("NormPopScore" = NormPopScore_ROH)
  }else{
    newDF = 
      newDF = dataFrame %>% 
      rename("NormPopScore" = NormPopScore_IBD)
  }
  Causualplot = ggplot(newDF, aes(x=newDF$NormPopScore, y=newDF$CausalVars)) + 
  geom_point(aes(colour = cut(newDF$CausalVars, c(-Inf, 0, 1, 5, 20))),
               size = 3) + 
  scale_color_manual(name = "Count Causal Variants", 
                       values = c("(-Inf,0]" = "black","(0,1]" = "yellow", "(1,5]" = "orange", "(5,20]" = "red"),
                       labels = c("0","1", "1 < variants <= 5", "5 < variants <= 20")) + 
  geom_text_repel(aes(label=ifelse(newDF$CausalVars >= 10 | newDF$NormPopScore > scoreCutOff, as.character(Population),'')), size = 6) + 
  theme_bw() +
  theme(plot.title=element_text(size = 18, face = "bold", hjust= 0.5), 
          axis.text.x = element_text(size = 24, vjust= 1, hjust= 0.5), 
          axis.text.y = element_text(size = 24), 
          axis.title=element_text(size= 24), 
          legend.title=element_text(size= 24), 
          legend.text=element_text(size= 18)) +
  labs(x=paste(xaxisLabel), y="Count Causal Variants")
  
  return(Causualplot)
}
######Plot Causal Vars with Correlation Fxn######
######Labelling the top 5% of Scores and ~1% of Causal Var Counts
plotCausalCorrs = function(regModel, dataFrame, varOfInterest, scoreCutOff, xaxisLabel){
  ggplotRegression(regModel)   + 
  geom_point(aes(colour = cut(dataFrame$CausalVars, c(-Inf, 0, 1, 5, 20))),
               size = 3) + 
  scale_color_manual(name = "Count Causal Variants", 
                       values = c("(-Inf,0]" = "black","(0,1]" = "yellow", "(1,5]" = "orange", "(5,20]" = "red"),
                       labels = c("0","1", "1 < variants <= 5", "5 < variants <= 20")) + 
  geom_text_repel(aes(label=ifelse(dataFrame$CausalVars >= 10 | dataFrame[,varOfInterest] > scoreCutOff, as.character(dataFrame$Population),'')), size = 6) + 
  theme_bw() +
  theme(plot.title=element_text(size = 18, face = "bold", hjust= 0.5), 
          axis.text.x = element_text(size = 24, vjust= 1, hjust= 0.5), 
          axis.text.y = element_text(size = 24), 
          axis.title=element_text(size= 24), 
          legend.title=element_text(size= 24), 
          legend.text=element_text(size= 18)) +
  labs(x=paste(xaxisLabel), y="Count Causal Variants") 
}
######Raincloud plotting fxns######
plotRainClouds = function(dataFrame, ylabTitle){
  dataFrame$Bin = cut(dataFrame$CausalVars, c(-Inf, 0, 1, 5, 20))
  ggplot(dataFrame, aes(x=Bin, y= NormPopScore, fill = Bin, colour = Bin)) +
  geom_flat_violin(size=2,position = position_nudge(x = .25, y = 0),adjust =2, trim = FALSE) +
  geom_point(aes(x = as.numeric(Bin)-.15, y = NormPopScore, colour = Bin), position = position_jitter(width = .05), size = 1, shape = 20) +
  geom_boxplot(aes(as.numeric(Bin), y = NormPopScore, fill = Bin), outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
  coord_flip() +
  guides(fill = FALSE, colour = FALSE) +
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  scale_x_discrete(labels=c("(-Inf,0]" = "0", "(0,1]" = "1", "(1,5]" = "1 < variants <= 5", "(5,20]" = "5 < variants <= 20"))+ 
    theme_bw() +
    theme(plot.title=element_text(size = 18, face = "bold", hjust= 0.5), 
          axis.text.x = element_text(size = 24, vjust= 1, hjust= 0.5), 
          axis.text.y = element_text(size = 24), 
          axis.title=element_text(size= 24), 
          legend.title=element_text(size= 24), 
          legend.text=element_text(size= 18)) +
    labs(y=paste(ylabTitle), x="Count Causal Variants") 
  }

######without Correlation
pROHScore = plotCausal(comboDF_noWolves, 200, "ROH Score in Mb (Normalized)", "ROH")
pIBDScore = plotCausal(comboDF_noWolves, 900, "IBD Score in Mb (Normalized)", "IBD")

#####Plot with correlation
plotPopularityCausVars = plotCausalCorrs(corrPopularitycausVars, PopularityDF, "OverallPopularityRank", 70, "Overall Popularity") 

plotFinalROHScoresCausVars = plotCausalCorrs(corrROHScorecausVars, comboDF_noWolves, "NormPopScore_ROH", 200, "ROH Score in Mb (Normalized)") 

plotFinalIBDScoresCausVars = plotCausalCorrs(corrIBDScorecausVars, comboDF_noWolves, "NormPopScore_IBD", 900, "IBD Score in Mb (Normalized)") 

#####Multiplot scores and Causals
OMIAplots = ggarrange(plotPopularityCausVars + theme(axis.title.y = element_blank()),
          ggarrange(plotFinalROHScoresCausVars + theme(axis.title.y = element_blank()), 
                    plotFinalIBDScoresCausVars + theme(axis.title.y = element_blank()), 
                    ncol = 2, 
                    labels = c("B", "C"), 
                    legend = "none"), # Second row with ROH and IBD plots
          nrow = 2, 
          legend = "none",
          labels = "A"                                        
) 

OMIAplots_addAxes = annotate_figure(OMIAplots, 
                                    left = text_grob("Count Causal Variants", 
                                                     size = 24, 
                                                     face = "bold", 
                                                     rot = 90))

ROHvsIBDCausals = ggarrange(plotFinalROHScoresCausVars + theme(axis.title.y = element_blank()), 
                            plotFinalIBDScoresCausVars + theme(axis.title.y = element_blank()),
                            nrow = 1, 
                            labels = c("A", "B"),
                            common.legend = T,
                            legend = "right")

ROHvsIBDCausals_addAxes = annotate_figure(ROHvsIBDCausals, 
                left = text_grob("Count Causal Variants", 
                                 size = 24, 
                                 face = "bold", 
                                 rot = 90))
print(OMIAplots_addAxes)
print(ROHvsIBDCausals_addAxes)
###Plot the ROH vs IBD relationship
plotROHvsIBD = ggplotRegression(corrROHvsIBD) + 
  labs(x="IBD Score (Mb)", y = "ROH Score (Mb)") +
  theme_bw() +
  theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), 
        axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=18), 
        legend.position = "bottom") 

plotROHvsIBDScores = ggplotRegression(corrROHScorevsIBDScore) + 
  labs(x="IBD Score(Mb) Normalized ", y = "ROH Score(Mb) Normalized") +
  theme_bw() +
  theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), 
        axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=18), 
        legend.position = "bottom")

ggarrange(plotROHvsIBD, 
          plotROHvsIBDScores,
          ncol = 2, 
          labels = c("A", "B"))


#color by clade
colourCount_pop = length(unique(comboDF$Clade)) 
palette = distinctColorPalette(colourCount_pop)
plotNormROHScorevsNormIBDScore = ggplotRegression(corrROHScorevsIBDScore) +
  geom_point(aes(colour=comboDF$Clade), size=3) + 
  scale_colour_manual(name= "Clade", values = palette, na.value="grey") +
  geom_text_repel(data=subset(comboDF, NormPopScore_ROH >= 200 | NormPopScore_IBD > 900), aes(label=paste(Population,",",CausalVars)), size = 6) + 
  labs(x="IBD Score in Mb (Normalized)", y="ROH Score in Mb (Normalized)") +
  theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), 
        axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), 
        axis.text.y = element_text(size  = 24), 
        axis.title=element_text(size=24),
        legend.title=element_text(size=24), 
        legend.text=element_text(size=18), 
        legend.position = "bottom") + 
  guides(colour = guide_legend(nrow = 6))

print(plotNormROHScorevsNormIBDScore)

###Raincloud plots
RCROH = plotRainClouds(FinalROHScores, "ROH Score in Mb (Normalized)")

RCIBD = plotRainClouds(FinalIBDScores, "IBD Score in Mb (Normalized)")

ggarrange(RCROH, RCIBD + xlab(NULL))

####Multiplot popularity with ROH and IBD Scores
ggarrange(ggplotRegression(corrPopularityROHScore) + 
            labs(x="Breed Popularity", y = "ROH Score(Mb) Normalized"), 
          ggplotRegression(corrPopularityIBDScore) + 
            labs(x="Breed Popularity", y = "IBD Score(Mb) Normalized"),
          ncol = 2, 
          labels = c("A", "B"))