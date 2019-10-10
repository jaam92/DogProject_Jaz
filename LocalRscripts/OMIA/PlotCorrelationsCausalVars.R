#####Load Libraries
source("~/DogProject_Jaz/LocalRscripts/OMIA/R_rainclouds.R")
source("~/DogProject_Jaz/LocalRscripts/OMIA/SummarizeData.R")
library(cowplot)
library(ggplot2)
library(grid)
library(gridExtra)
library(gtable)
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

######Plot Causal Vars Fxn ######
######Labelling the top 5% of Scores and ~1% of Causal Var Counts
plotCausal = function(dataFrame, scoreCutOff, xaxisLabel){
  ggplot(dataFrame, aes(x=dataFrame$NormPopScore, y=dataFrame$CausalVars)) + 
  geom_point(aes(colour = cut(dataFrame$CausalVars, c(-Inf, 0, 1, 5, 20))),
               size = 3) + 
  scale_color_manual(name = "Count Causal Variants", 
                       values = c("(-Inf,0]" = "black","(0,1]" = "yellow", "(1,5]" = "orange", "(5,20]" = "red"),
                       labels = c("0","1", "1 < variants <= 5", "5 < variants <= 20")) + 
  geom_text_repel(aes(label=ifelse(dataFrame$CausalVars >= 10 | dataFrame$NormPopScore > scoreCutOff, as.character(Population),''))) + 
  theme_bw() +
  theme(plot.title=element_text(size = 18, face = "bold", hjust= 0.5), 
          axis.text.x = element_text(size = 24, vjust= 1, hjust= 0.5), 
          axis.text.y = element_text(size = 24), 
          axis.title=element_text(size= 24), 
          legend.title=element_text(size= 24), 
          legend.text=element_text(size= 18)) +
  labs(x=paste(xaxisLabel), y="Count Causal Variants") 
}

######Plot Causal Vars with Correlation Fxn ######
######Labelling the top 5% of Scores and ~1% of Causal Var Counts
plotCausalCorrs = function(regModel, dataFrame, varOfInterest, scoreCutOff, xaxisLabel){
  selectVarofInterestCol = enquo(varOfInterest)
  ggplotRegression(regModel) + 
  ggplot(dataFrame, aes(x=dataFrame$NormPopScore, y=dataFrame$CausalVars)) + 
  geom_point(aes(colour = cut(dataFrame$CausalVars, c(-Inf, 0, 1, 5, 20))),
               size = 3) + 
    scale_color_manual(name = "Count Causal Variants", 
                       values = c("(-Inf,0]" = "black","(0,1]" = "yellow", "(1,5]" = "orange", "(5,20]" = "red"),
                       labels = c("0","1", "1 < variants <= 5", "5 < variants <= 20")) + 
  geom_text_repel(aes(label=ifelse(dataFrame$CausalVars >= 10 | dataFrame$selectVarofInterestCol > scoreCutOff, as.character(Population),''))) + 
  theme_bw() +
  theme(plot.title=element_text(size = 18, face = "bold", hjust= 0.5), 
          axis.text.x = element_text(size = 24, vjust= 1, hjust= 0.5), 
          axis.text.y = element_text(size = 24), 
          axis.title=element_text(size= 24), 
          legend.title=element_text(size= 24), 
          legend.text=element_text(size= 18)) +
  labs(x=paste(xaxisLabel), y="Count Causal Variants") 
}

####Plot without Correlation
pROHScore = plotCausal(FinalROHScores, 200, "ROH Score in Mb (Normalized)")

pIBDScore = plotCausal(FinalROHScores, 900, "IBD Score in Mb (Normalized)")

#####Multiplot scores and Causals
AlternativeScoresCausal = plot_grid(pROHScore + theme(legend.position="none"), 
                                    pIBDScore + theme(legend.position="none", axis.title.y = element_blank()),
                                    nrow = 1, align = "hv")

AlternativeFinalPlotScoresCausal= plot_grid(AlternativeScoresCausal, legendScores, rel_widths = c(2, 0.4))

#####Plot with correlation
plotPopularityCausVars = ggplotRegression(corrPopularitycausVars)  +  
  geom_point(aes(colour = cut(PopularityDF$CausalVars, c(-Inf, 0, 1, 5, 20))), size = 3) + 
  scale_color_manual(name = "Count Causal Variants", values = c("(-Inf,0]" = "black","(0,1]" = "yellow", "(1,5]" = "orange", "(5,20]" = "red"),labels = c("0","1", "1 < variants <= 5", "5 < variants <= 20")) + 
  geom_text_repel(data=subset(PopularityDF, CausalVars >= 10 | OverallPopularityRank >= 160), aes(label=Population)) + 
  labs(x="Ranked Overall Popularity", y="Count Causal Variants") 

plotFinalROHScoresCausVars = ggplotRegression(corrROHScorecausVars) + geom_point(aes(colour = cut(FinalROHScores$CausalVars, c(-Inf, 0, 1, 5, 20))), size = 3) + scale_color_manual(name = "Count Causal Variants", values = c("(-Inf,0]" = "black","(0,1]" = "yellow", "(1,5]" = "orange", "(5,20]" = "red"),labels = c("0","1", "1 < variants <= 5", "5 < variants <= 20")) + geom_text_repel(data=subset(FinalROHScores, CausalVars >= 10 | NormPopScore >= 200), aes(label=Population)) + labs(x="ROH Score in Mb (Normalized)", y="Count Causal Variants") + theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=18)) 

plotFinalIBDScoresCausVars = ggplotRegression(corrIBDScorecausVars) + geom_point(aes(colour = cut(FinalIBDScores$CausalVars, c(-Inf, 0, 1, 5, 20))), size = 3) + scale_color_manual(name = "Count Causal Variants", values = c("(-Inf,0]" = "black","(0,1]" = "yellow", "(1,5]" = "orange", "(5,20]" = "red"),labels = c("0","1", "1 < variants <= 5", "5 < variants <= 20")) + geom_text_repel(data=subset(FinalIBDScores, CausalVars >= 10 | NormPopScore > 900), aes(label=Population)) + labs(x="IBD Score in Mb (Normalized)", y="Count Causal Variants") + theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=18))

plotMeanROHCountCausVars = ggplotRegression(corrROHCountcausVars) + geom_point(aes(colour = cut(FinalROHScores$CausalVars, c(-Inf, 0, 1, 5, 20))), size = 3) + scale_color_manual(name = "Count Causal Variants", values = c("(-Inf,0]" = "black","(0,1]" = "yellow", "(1,5]" = "orange", "(5,20]" = "red"),labels = c("0","1", "1 < variants <= 5", "5 < variants <= 20")) + geom_text_repel(data=subset(FinalROHScores, CausalVars >= 10 | MeanROHperIndivCount > 141.5), aes(label=Population)) + labs(x="Mean ROH (Mb) per Individual", y="Count Causal Variants") + theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=18))

plotMeanIBDCountCausVars = ggplotRegression(corrIBDCountcausVars) + geom_point(aes(colour = cut(FinalIBDScores$CausalVars, c(-Inf, 0, 1, 5, 20))), size = 3) + scale_color_manual(name = "Count Causal Variants", values = c("(-Inf,0]" = "black","(0,1]" = "yellow", "(1,5]" = "orange", "(5,20]" = "red"),labels = c("0","1", "1 < variants <= 5", "5 < variants <= 20")) + geom_text_repel(data=subset(FinalIBDScores, CausalVars >= 10 | MeanIBDperIndivCount > 3720.8), aes(label=Population))  + labs(x="Mean IBD (Mb) per Individual", y="Count Causal Variants") + theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=18))

plotROHvsIBD = ggplotRegression(corrROHvsIBD) + labs(x="IBD Score (Mb)", y = "ROH Score (Mb)") + theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=18)) 

#expand color palette of choice to have number of colors equal to number of clades
colourCount_pop = length(unique(comboDF$Clade)) 
palette = distinctColorPalette(colourCount_pop)
plotNormROHScorevsNormIBDScore = ggplotRegression(corrROHScorevsIBDScore) + geom_point(aes(colour=comboDF$Clade), size=3) + scale_colour_manual(name= "Clade", values = palette, na.value="grey") + geom_text_repel(data=subset(comboDF, NormPopScore_ROH >= 200 | NormPopScore_IBD > 900), aes(label=paste(Population,",",CausalVars))) + labs(x="IBD Score in Mb (Normalized)", y="ROH Score in Mb (Normalized)") + theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24)) + theme(legend.title=element_text(size=24), legend.text=element_text(size=18), legend.position = "bottom") + guides(colour = guide_legend(nrow = 6))

###Raincloud plots by bin
FinalROHScores$Bin = cut(FinalROHScores$CausalVars, c(-Inf, 0, 1, 5, 20))
FinalIBDScores$Bin = cut(FinalIBDScores$CausalVars, c(-Inf, 0, 1, 5, 20))

#plotROH
RCROH = ggplot(FinalROHScores,aes(x=Bin,y=NormPopScore, fill = Bin, colour = Bin))+
  geom_flat_violin(size=2,position = position_nudge(x = .25, y = 0),adjust =2, trim = FALSE)+
  geom_point(aes(x = as.numeric(Bin)-.15, y = NormPopScore, colour = Bin),position = position_jitter(width = .05), size = 1, shape = 20)+
  geom_boxplot(aes(as.numeric(Bin), y = NormPopScore, fill = Bin),outlier.shape = NA, alpha = .5, width = .1, colour = "black") +
  labs(y="ROH Score in Mb (Normalized)", x="Count Causal Variants")+coord_flip()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  scale_x_discrete(labels=c("(-Inf,0]" = "0","(0,1]" = "1", "(1,5]" = "1 < variants <= 5", "(5,20]" = "5 < variants <= 20"))+ theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24))

RCIBD = ggplot(FinalIBDScores,aes(x=Bin,y=NormPopScore, fill = Bin, colour = Bin))+
  geom_flat_violin(size=2,position = position_nudge(x = .25, y = 0),adjust =2, trim = FALSE)+
  geom_point(aes(x = as.numeric(Bin)-.15, y = NormPopScore, colour = Bin),position = position_jitter(width = .05), size = 1, shape = 20)+
  geom_boxplot(aes(as.numeric(Bin), y = NormPopScore, fill = Bin),outlier.shape = NA, alpha = .5, width = .1, colour = "black")+
  labs(y="IBD Score in Mb (Normalized)", x="Count Causal Variants")+coord_flip()+guides(fill = FALSE, colour = FALSE) +
  scale_colour_brewer(palette = "Dark2")+
  scale_fill_brewer(palette = "Dark2")+
  scale_x_discrete(labels=c("(-Inf,0]" = "0","(0,1]" = "1", "(1,5]" = "1 < variants <= 5", "(5,20]" = "5 < variants <= 20"))+ theme(plot.title=element_text(size =18, face = "bold", hjust=0.5), axis.text.x = element_text(size  = 24, vjust=1, hjust=0.5), axis.text.y = element_text(size  = 24), axis.title=element_text(size=24))

plot_grid(RCROH, RCIBD + xlab(NULL))


#####Multiplot scores and Causals
plotScoresCausal = plot_grid(plotFinalROHScoresCausVars + theme(legend.position="none"), plotFinalIBDScoresCausVars + theme(legend.position="none", axis.title.y = element_blank()), nrow = 1, align = "hv")
legendScores = get_legend(plotFinalROHScoresCausVars) #pull universal legend
FinalPlotScoresCausal= plot_grid(plotScoresCausal, legendScores, rel_widths = c(2, 0.4))

plot_grid(AlternativeFinalPlotScoresCausal, plotNormROHScorevsNormIBDScore, ncol = 1) 
plot_grid(FinalPlotScoresCausal, plotNormROHScorevsNormIBDScore, ncol = 1) 

plot_grid(plotFinalROHScoresCausVars + theme(legend.position="none"), plotFinalIBDScoresCausVars + theme(legend.position="none", axis.title.y = element_blank()), plotPopularityCausVars + theme(legend.position="none"), nrow = 2, legendScores)
