#data: visitors
data1<-read.csv("/Users/alysonmak/Downloads/visitors_rawdata.csv")
samplesize=60
a<-nrow(data1)
b<-a1-(samplesize-1)
x1<-data1[b:a,]
c1<-as.numeric(gsub(",", "", x1$Unique.Visits))
mean1<-mean(c1)
sd1<-sd(c1)
Unique.Visits<-(c1-mean1)/sd1
data1<-matrix(Unique.Visits,ncol=1)
colnames(data1)="x"
path1="/Users/alysonmak/Desktop/visitors.csv"
write.csv(data1,path1,row.names=FALSE)

#data: sales
data2<-read.csv("/Users/alysonmak/Downloads/sales_rawdata.csv")
data2<-matrix(data2[1:samplesize,]$ProductP4,ncol=1)
colnames(data2)="x"
path2="/Users/alysonmak/Desktop/sales.csv"
write.csv(data2,path2,row.names=FALSE)