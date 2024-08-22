#Latin square design
a<-read.csv("c://5.1cake.csv", header=TRUE)
a ; attach(a)
cake <- as.factor(a$cake)
mixer <- as.factor(a$mixer)
oven <- as.factor(a$oven)

by(y,cake,FUN=function(x){ c(mean(x), sd(x))})
by(y,mixer,FUN=function(x){ c(mean(x), sd(x))})
by(y,oven,FUN=function(x){ c(mean(x), sd(x))})

a1 <- lm(y~oven+mixer+cake) # 순서: row+col+trt
anova(a1)
pairwise.t.test(y,cake,p.adjust="none")

par(mfrow=c(1,3))
  plot(y~oven+mixer+cake)
  
par(mfrow=c(2,2))
  plot(a1)
  
detach(a)

#2반복 이요인설계
press<-C(100,100,150,150,100,100,150,150)
temp<-c(50,50,50,50,100,100,100,100)
y<-c(9,11,14,16,11,13,14,20)
press<-as.factor(press)
temp<-as.factor(temp)
dd<-aov(y~press*temp)
summary(dd);
model.tables(dd)

#부분교락법
da<-read.csv("c://예제7.9.csv", header=T)
da
attach(da)
rep<-as.factor(Rep) #반복
block<-as.factor(block)
A<-as.factor(A)
B<-as.factor(B)
C<-as.factor(C)

par(mfrow=c(1,2))
boxplot(y~rep, sub="rep") ; boxplot(y~block, sub="block")

op=par(mfrow=c(1,3))
boxplot(y~A, sub="factor A"); boxplot(y~B, sub="factor B"); boxplot(y~C, sub="factor C")
par(op)

L<-aov(y~rep+block+A+B+C)
anova(L)
model.tables(L, type="means") # 모든 요인 수준평균
model.tables(L, type="effects") # 모든 요인 수준효과

L1<-aov(y~rep+block+A+B+C+A:B+A:C)
anova(L1)
model.tables(L1, type="means") # 모든 요인 수준평균
model.tables(L1, type="effects") # 모든 요인 수준효과

#부분요인설계
frac<-read.csv("C://그림 8.3데이터.csv", header=T)
frac
attach(frac)
A<-as.factor(A)
B<-as.factor(B)
C<-as.factor(C)
D<-as.factor(D)

op=par(mfrow=c(1,4))
boxplot(y~A, sub="A") ; boxplot(y~B, sub="B")
boxplot(y~C, sub="C") ; boxplot(y~D, sub="D")

tapply(y,A,mean) ; tapply(y,B,mean)
tapply(y,C,mean) ; tapply(y,D,mean)

faov1 <- aov(y~A+B+C+D)
summary(faov1)
model.tables(faov1) # model parameters