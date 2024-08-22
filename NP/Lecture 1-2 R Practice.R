
############################################
########## Lecture 1-2 R Practice ##########
############################################

# install.packages("tidyr")
# install.packages("ggplot2")
 install.packages("fGarch") # Try: options("install.lock"=FALSE), if Error occurs
install.packages("gtools")
# install.packages("dplyr")
install.packages("DescTools")

library(tidyr)
library(ggplot2)
library(fGarch)
library(gtools)
library(dplyr)
library(DescTools)


##### Lecture 1: One-Sample Methods #####

## Example: Sodium contents (in mg) of 40 servings of a food product
x <- read.csv("sodium.csv")[,1]
x

n <- length(x); n

mean(x); sd(x)


## (1) T-test: ???균의 검???
## H_0: mean(x) = 75 vs. H_a: mean(x) != 75

g <- data.frame(x=x) %>% ggplot(aes(x=x))
g + geom_histogram(aes(y=..density..), color="white", bins=30)

qqnorm(x); qqline(x, col=2)

t.test(x, alternative="two.sided", mu=75)

## Consider a sample of size 40 from a skewed Normal dist.
y <- fGarch::rsnorm(40, mean=mean(x), sd=sd(x), xi=1.5)

par(mfrow=c(1,2))
hist(x, main="Histogram of Sodium Data (x)")
hist(y, main="Histogram of Hypothetical Data (y)")
par(mfrow=c(1,1))

## Power of T-test with right-skewed data
mu.a <- 75 + 1
alpha=0.05
R=1000 # number of replications
rejections <- numeric()
for (r in 1:R) {
  y <- rsnorm(40, mean=mu.a, sd=sd(x), xi=1.5)
  rejections[r] <- t.test(y, alternative="two.sided", mu=75)$p.value < alpha
}
mean(rejections)



## (2) Binomial (sign) test: 중앙값의 검???
## H_0: median(x) = 75 vs. H_a: median(x) > 75

median(x)

signs <- sign(x-75)
S <- sum(signs > 0); S

## 만약 H_a가 참이??????, S ~ B(n,p) with p > 0.5
binom.test(S, n=n, p=0.5, alternative="greater") # exact test using Binomial dist.`
pbinom(S-1, size=n, prob=0.5, lower.tail=FALSE)

## 95% CI (X_(a), X_(b)) for median(x)
x.ordered <- sort(x)
# at least a obs. must fall less than median(x)
a <- qbinom(0.025, size=n, prob=0.5)
# at most (b-1) obs. must fall less than median(x)
b <- qbinom(0.975, size=n, prob=0.5) + 1
ci <- c(x.ordered[a], x.ordered[b]); ci


prop.test(S, n=n, p=0.5, alternative="greater", correct=FALSE) # normal approx. test
# prop.test(S, n=n, p=0.5, alternative="greater", correct=TRUE)

z <- (S-0.5*n)/sqrt(0.25*n)
pnorm(z, 0, 1, lower.tail=FALSE)

## Approx. 95% CI (X_(a), X_(b)) for median(x)
a.temp <- 0.5*n - qnorm(0.975, 0, 1) * sqrt(0.25*n)
a <- round(a.temp)
b.temp <- (0.5*n + 1) + qnorm(0.975, 0, 1) * sqrt(0.25*n)
b <- round(b.temp)
ci <- c(x.ordered[a], x.ordered[b]); ci





##### Lecture 2: Two-Sample Methods #####

## Example: 체체???리의 ?????? 종들 ??? ?????? ?????? 종을 9마리, 강에 ?????? 종을 6마리 ???????????? 각각??? ????????? 길이??? ?????????.
flies <- read.csv("tsetse_flies.csv")
head(flies)

m <- length(flies$length[flies$label=="forest"]); m
n <- length(flies$length[flies$label=="river"]); n
N <- m + n

ggplot(flies, aes(x=label, y=length)) + 
  geom_boxplot() + 
  geom_dotplot(binaxis="y", stackdir="center", dotsize=.5)



## (1): Permutation test: ???치모?????? 검???
## x = ?????? ?????? 종의 ????????? 길이/ y = 강에 ?????? 종의 ????????? 길이
## Assume $F_x(z) = F_y(z-Delta)
## H_0: Delta = 0 vs. H_a: Delta > 0

D.observed <- with(flies, mean(length[label=="forest"])-mean(length[label=="river"]))
D.observed

all.perms <- gtools::combinations(N, m) # Consider all permutations
head(all.perms)

D <- numeric()
for (i in 1:nrow(all.perms)) {
  length.forest.permuted <- flies$length[all.perms[i,]]
  length.river.permuted <- flies$length[-all.perms[i,]]
  D[i] <- mean(length.forest.permuted) - mean(length.river.permuted)
}

ggplot(data.frame(D), aes(D)) + 
  geom_histogram(binwidth = 0.005, color="white") +
  geom_vline(xintercept=D.observed, color='red')

p.value <- sum(D >= D.observed)/nrow(all.perms) # exact P-value
p.value


nrow(all.perms)

B <- 1000 # Consider only 1000 permuations
set.seed(0318)
D <- numeric()
for (i in 1:B) {
  rand.id <- sample(1:N, m, replace=FALSE)
  length.forest.permuted <- flies$length[rand.id]
  length.river.permuted <- flies$length[-rand.id]
  D[i] <- mean(length.forest.permuted) - mean(length.river.permuted)
}

p.value <- 2*(sum(D >= D.observed))/B # approx. P-value
p.value



## (2) Wilcoxon-Mann-Whitney test: ???치모?????? 검???
## H_0: Delta = 0 vs. H_a: Delta > 0

flies.rank <- flies %>% mutate(rank=rank(length))
flies.rank

W1.observed <- with(flies.rank, sum(rank[label=="forest"])) # Wilcoxon rank sum
W1.observed

U.observed <- W1.observed - m*(m+1)/2 # Mann-Whitney U
U.observed

# Caution: dwilcox, pwilcox??? Mann-Whitney U ???계량 기??? (???콕슨 ???????????? ??????)
data.frame(x=0:m*n, density=dwilcox(0:m*n, m, n)) %>%
  ggplot(aes(x=x, y=density)) + geom_line() + 
  geom_vline(xintercept=U.observed, color="red")

pwilcox(U.observed, m=m, n=n, lower.tail=FALSE) #pvalue

x <- with(flies, length[label=="forest"])
y <- with(flies, length[label=="river"])
x

# Caution: wilcox.test??? W ???계량 = U ???계량???
wilcox.test(x, y, alternative="greater") #alternative -> x?? ???? ?? Ŭ ???̴?
# wilcox.test(length ~ label, data=flies)

## 95% Confidence Interval for Delta
pwd <- sort(outer(x, y, "-")) # pairwise difference, Xi - Yj
pwd

a <- qwilcox(0.025, m=m, n=n) + 1 #??��???? ??ũ???? ?ƴ϶? ????Ʈ?? ??Ÿ????
b <- qwilcox(0.025, m=m, n=n, lower.tail=F) + 1
ci <- c(pwd[a], pwd[b]); ci

median(pwd) # Hodges-Lehammn estimate

wilcox.test(x, y, conf.int = TRUE)
DescTools::HodgesLehmann(x, y, conf.level=0.95)
# None of these are exact as there are ties



## (3) Siegel-Tukey test: 척도모수??? ?????? 검??? (equal location)

DescTools::SiegelTukeyTest(x,y)
#equal location?϶? ?? ???? ?л??? ??��?? ??��


## (4) Ansari-Bradley test: 척도모수??? ?????? 검??? ((un)equal location)

ansari.test(x,y) # equal location 
ansari.test(x - median(x), y - median(y)) # unequal location

# Note: We can use F-test when data are normally distributed
var.test(x,y)



## (5) Permutation test on deviances ((un)equal location)

x.mad <- mean(abs(x-median(x))) # mean absolute deviation
y.mad <- mean(abs(y-median(y)))

# ratio of mads
rmd.observed <- x.mad/y.mad; rmd.observed

B <- 1000
rmd <- numeric()
for (i in 1:B){
  ids = sample(1:N, m, replace=FALSE)
  permuted.x = with(flies, length[ids])
  permuted.y = with(flies, length[-ids])
  permuted.x.mad <- mean(abs(permuted.x - median(permuted.x)))
  permuted.y.mad <- mean(abs(permuted.y - median(permuted.y)))
  rmd[i] <- permuted.x.mad / permuted.y.mad
} 

ggplot(data.frame(rmd), aes(rmd)) + 
  geom_histogram(binwidth = 0.05, color="white") +
  geom_vline(xintercept=rmd.observed, color='red')

p.value <- 2*min( c(sum(rmd >= rmd.observed)/B, sum(rmd <= rmd.observed)/B) )
p.value



## (6) Kolmogorov-Smirnov test
## H_0: F_x = F_y vs. H_a: F_x != F_y

plot(ecdf(x), col="blue", xlim=c(1, 1.7))
lines(ecdf(y), col="red")

ks.test(x,y)

