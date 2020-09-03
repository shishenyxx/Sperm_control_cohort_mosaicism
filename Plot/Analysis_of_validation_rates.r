library(ggplot2)

raw<-read.csv(file="Detected_numbers_YA_and_AA.csv",header=T)

ggplot(raw)+
	geom_smooth(aes(x=AGE,y=BL_ND,col="Sperm"),method='lm',formula=y~x,col="red")+
	geom_smooth(aes(x=AGE,y=1-Pos_BL,col="Blood"),method='lm',formula=y~x,col="blue")+
	geom_point(aes(x=AGE,y=BL_ND,col="Sperm"),alpha=0.5,size=2,col="red")+
	geom_point(aes(x=AGE,y=1-Pos_BL,col="Blood"),alpha=0.5,size=2,col="blue")+
	ylim(0,1)+
	theme_bw()
