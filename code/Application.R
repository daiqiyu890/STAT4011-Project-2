#Step 1: data pre-processing----------------------------------------------------
#visitors
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

#sales
data2<-read.csv("/Users/alysonmak/Downloads/sales_rawdata.csv")
data2<-matrix(data2[1:samplesize,]$ProductP4,ncol=1)
colnames(data2)="x"


#Step 2: Model Selection--------------------------------------------------------


#Step 3: Parameter Estimation---------------------------------------------------

#Step 4: Forecast---------------------------------------------------------------


#Step 4: Plot-------------------------------------------------------------------
