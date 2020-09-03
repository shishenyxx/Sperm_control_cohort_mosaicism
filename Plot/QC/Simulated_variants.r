library(ggplot2)

raw<-read.delim(file="sensitivity.txt",header=T)

raw2<-raw[raw$AF!=50,]

bks=c(1,2,3,4,5,10,15,20,25)

ggplot(raw2,aes(x=sqrt(AF),y=Sensitivity,col=factor(DEPTH),group=DEPTH))+
	geom_point(shape=2)+
	geom_line()+
	scale_x_continuous(breaks=sqrt(bks),labels=bks)+
	ylim(0,1)+
	theme_classic()
