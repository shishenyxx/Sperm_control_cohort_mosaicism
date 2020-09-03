library(ggplot2)

raw<-read.delim(file="summary_percentage_of_coverage_with_group.txt",header=T)

ggplot(raw,aes(x=DEPTH,y=RATIO,col=TISSUE,group=ID))+
	geom_line()+
	scale_colour_manual(values = c("#ff5b00","#107010"))+
	theme_classic()
