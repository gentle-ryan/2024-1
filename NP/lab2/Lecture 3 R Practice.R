
##########################################
########## Lecture 3 R Practice ##########
##########################################

library(ggplot2)
library(dplyr)
# install.packages("agricolae")
library(agricolae)
library(DescTools)


### (1) Permuational F-test

### Example: 주말 오후 1-3시에 서울 잠실 일대 공유 전동킥보드 세 브랜드(킥고잉, 라임, 고고씽)의 이용 시간을 조사하였다.

kick <- read.csv("kickboard.csv")
head(kick)

group_by(kick, brand) %>%
  summarise(count=n(), mean=mean(time), sd=sd(time))

ggplot(kick, aes(x=time, color=brand, fill=brand)) +
  geom_density(alpha=0.5) + geom_rug()

summary(aov(time ~ brand, data=kick))


N <- nrow(kick); K <- 3 # N에 총데이터 수 k는 그룹수
group.length <- aggregate(kick$time, by=list(kick$brand), length)[,2]
group.mean <- aggregate(kick$time, by=list(kick$brand), mean)[,2]

SSTotal <- var(kick$time) * (length(kick$time)-1)
SST <- sum(group.length * (group.mean - mean(kick$time))^2)
SSE <- SSTotal - SST

F.obs <- (SST/(K-1))/(SSE/(N-K)); F.obs

# 데이터가 정규성을 만족하지 않으니 퍼뮤테이션!
R = 999 # number of permutations
F.perm <- vector(length=R)
for (r in 1:R) {
  perm <- sample(1:N, replace=FALSE)
  group.length <- aggregate(kick$time, by=list(kick$brand[perm]), length)[,2]
  group.mean <- aggregate(kick$time, by=list(kick$brand[perm]), mean)[,2]
  
  SSTotal <- var(kick$time) * (length(kick$time)-1)
  SST <- sum(group.length * (group.mean - mean(kick$time))^2)
  SSE <- SSTotal - SST
  F.perm[r] <- (SST/(K-1))/(SSE/(N-K))
}

p.value <- (sum(F.perm>=F.obs)+1)/(R+1); p.value

data.frame(F.perm=F.perm) %>% ggplot(aes(x=F.perm)) +
  geom_histogram(aes(y=..density..), binwidth=0.1, color="white") +
  geom_vline(xintercept=F.obs, color='red')



### (2) Kruskal-Wallis test: without ties

### Example: 떡갈나무 16그루에 대해 연간 밑둥 직경의 증가량을 측정하였다. 각 나무는 위치한 부지의 경사도에 따라 분류되었다. (1: 경사 급함 \~ 5: 경사 없음)

tree <- read.csv("tree.csv")
head(tree)

ggplot(tree, aes(x=as.factor(slope), y=increment)) + geom_boxplot()

kruskal.test(increment ~ slope, data=tree) #카이스퀘어 근사


tree.rank <- tree %>% mutate(rank=rank(increment))
head(tree.rank)

group.length <- aggregate(rank ~ slope, data=tree.rank, length)[,2]
group.rank.mean <- aggregate(rank ~ slope, data=tree.rank, mean)[,2]

N <- sum(group.length); K <- length(group.length)
KW.obs <- ( 12/(N*(N+1))) * sum( group.length * (group.rank.mean - (N+1)/2)^2 )
KW.obs

p.value <- pchisq(KW.obs, df=K-1, lower.tail=F) #카이스퀘어 근사를 했으니까 pchisq
p.value



### (3) Jonckheere-Terpstra Test: ordered alternative

tree$slope <- factor(tree$slope, ordered=TRUE) #ordered factor로 변환
# exact test
DescTools::JonckheereTerpstraTest(increment ~ slope, data=tree, alternative="increasing")



### (4) Multiple testing: controlling FWER
#20개의 sample로 부터 20개의 귀무가설(ml = 0)을 테스트
### The following simulations shows that setting $\alpha=0.05$ for each test leads to undesired results.

L=20 # number of tests
n=100 # size of samples
R=1000
err.vec <- vector(length=R)
for (r in 1:R) {
  fam.err <- FALSE
  x=matrix(rnorm(n*L), nrow=n)
  for (i in 1:L) {
    if (t.test(x[,i], mu=0)$p.value < 0.05) { #h0가 참이라는 가설하에 test
      fam.err <- TRUE #false rejection이 일어나면 fwer를 true로
      break
    }
  }
  err.vec[r] <- fam.err
}

(fam.err.rate <- sum(err.vec)/R)

1 - (1-0.05)^20

### Need to use a smaller significance level $\alpha\prime$ for each pairwise test!



### (5) Multiple comparisons: Bonferroni adjustment

### Example: Percentages of clay in soil samples from four locations. (Populations are normally distributed with equal variances.)

clay <- read.csv("clay.csv")
head(clay)

group_by(clay, location) %>%
  summarise(count=n(), mean=mean(pct), sd=sd(pct))

summary(aov(pct ~ as.factor(location), data=clay))


pct.all <- matrix( c(clay$pct[clay$location==1], 
                     clay$pct[clay$location==2],
                     clay$pct[clay$location==3],
                     clay$pct[clay$location==4]), 
                   byrow=FALSE, ncol=4 )

pvals <- matrix(NA, nrow=4, ncol=4)
for (i in 1:4) {
  j <- i+1
  while(j<=4) {
    pvals[i,j] <- t.test(pct.all[,i], pct.all[,j], var.equal=TRUE)$p.value
    j <- j+1}
}

pvals # raw p-values

pvals[is.na(pvals)==FALSE] <- p.adjust(pvals[is.na(pvals)==FALSE], method="bonferroni")
pvals # Bonferroni-adjusted p-values

t.test(clay$pct[clay$location==2], clay$pct[clay$location==4], var.equal=TRUE)$p.value*6

(bonf.critical.value <- qt(0.025/6, df=2*6-2, lower.tail=F))



### (6) Multiple comparisons: Fisher's protected LSD

N <- nrow(clay); K <-4
clay$location <- as.factor(clay$location)
model <- aov(pct ~ location, data=clay)

MSE <- deviance(model)/df.residual(model)

(lsd.critical.value <- qt(0.025, df=N-K, lower.tail=F) * sqrt(MSE*(2/6)))

(agricolae::LSD.test(model, "location", p.adj="none"))

#유의하지 않은 것들끼리 하나의 그룹으로 묶임(?)

### (7) Multiple comprisons: Tukey's HSD

(out <- TukeyHSD(model, "location", ordered=FALSE))

Tij.obs <- abs(out$location[,1]); Tij.obs

# Permutational HSD
R = 2000
Q.perm <- vector(length=R)
for (r in 1:R) {
  clay.perm <- data.frame(pct = clay$pct, location=clay$location[sample(1:N, replace=FALSE)])
  model <- aov(pct ~ location, data=clay.perm)
  MSE <- deviance(model)/df.residual(model)
  group.mean <- aggregate(pct ~ location, data=clay.perm, mean)[,2]
  Tij <- abs(outer(group.mean, group.mean, "-"))/ sqrt(MSE/6)
  Q.perm[r] <- max(Tij)
} 

( p.values <- unlist(purrr::map(Tij.obs, function(x) {sum(Q.perm>=x)/R})) )

data.frame(Q.perm=Q.perm) %>% ggplot(aes(Q.perm)) +
  geom_histogram(aes(y=..density..), binwidth=0.1, color="white", alpha=0.6) +
  geom_vline(xintercept=Tij.obs, color="red") +
  geom_text(data=data.frame(x=Tij.obs, label=names(Tij.obs)),
            aes(x=x, y=0.1, label=label), color="red")