##Codes to plot the PCA results

library(ggplot2)

raw<-read.csv(file="PCA_result_with_group.csv",header=T)

ggplot(raw,aes(x=PC2,y=PC3,col=POPULATION,fill=POPULATION,shape=GROUP))+
	geom_point(alpha=raw$ALPHA)+
	geom_segment(aes(x = -0.05, y = -0.04, xend = -0.05, yend = 0.04),color="black")+
	geom_segment(aes(x = -0.05, y = -0.07, xend = 0.05, yend = -0.07),color="black")+
	theme_bw()+
	theme(panel.border=element_blank(),axis.line=element_blank(),panel.grid =element_blank())
