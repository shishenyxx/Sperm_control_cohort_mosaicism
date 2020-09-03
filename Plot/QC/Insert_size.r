
library(ggplot2)

raw<-read.delim(file="summary_insert_size.txt",header=T)

raw2<-raw[raw$cohort=="Old"|raw$cohort=="Young",]


ggplot(raw2,aes(x=size,y=count,col=group,group=name))+
	geom_line()+
	scale_colour_manual(values = c("#ff5b00","#107010"))+
	theme_classic()
