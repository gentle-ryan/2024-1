## ---- Regression Analysis and Lab. --- Lab 2 ---------------------------------


## -----------------------------------------------------------------------------
deliver <- read.csv("ex3.1.csv")
colnames(deliver) <- c("obs", "time", "case", "distance")
head(deliver)
plot(deliver[,2:4]) # scatterplot matrix


## -----------------------------------------------------------------------------
lmfit.full <- lm(time ~ case + distance, data=deliver)
summary(lmfit.full)


## -----------------------------------------------------------------------------
names(lmfit.full)
lmfit.full$residuals # e
lmfit.full$fitted.values # y^

names(summary(lmfit.full))
summary(lmfit.full)$r.squared # R^2
predict(lmfit.full) # y.hat; or use [lmfit.full$fitted.values]


## -----------------------------------------------------------------------------
lmfit.full <- lm(time ~ case + distance, data=deliver)
lmfit.reduced <- lm(time ~ case, data=deliver)
anova(lmfit.reduced, lmfit.full) # is distance significant? 


## -----------------------------------------------------------------------------
anova(lm(time ~ I(case + distance), data = deliver),
      lm(time ~ case + distance, data = deliver)) # reject null


## -----------------------------------------------------------------------------
lmfit.full <- lm(time ~ case + distance, data=deliver)
lmfit.reduce <- lm(time ~ 1, data = deliver)
anova(lmfit.reduce, lmfit.full)


## -----------------------------------------------------------------------------
confint(lmfit.full) # 95% confidence interval


## -----------------------------------------------------------------------------
predict(lmfit.full, data.frame(case = c(5, 10, 15), distance=c(700, 500, 300)))
predict(lmfit.full, data.frame(case = c(5, 10, 15), distance=c(700, 500, 300)), interval="confidence")
predict(lmfit.full, data.frame(case = c(5, 10, 15), distance=c(700, 500, 300)), interval="predict")


## -----------------------------------------------------------------------------
windmill <- read.csv("windmill.csv")
fit1 <- lm(output ~ wind, data = windmill) # (1)
fit2 <- lm(output ~ I(1/wind), data = windmill) # (2)
fit3 <- lm(output ~ wind + I(wind^2), data = windmill) # (3)
par(mfrow=c(1,3))
plot(fit1, which = 1)
plot(fit2, which = 1)
plot(fit3, which = 1)
par(mfrow=c(1,1))
plot(output ~ wind, data = windmill)
abline(fit1)
lines(windmill$wind, predict(fit2), col = "red")
lines(windmill$wind, predict(fit3), col = "blue")
legend("bottomright", c("linear", "reciprocal", "quadratic"), lty=1, col=c("black", "red", "blue"))


rocket <- read.csv("rocket.csv")
fit1 <- lm(x ~ I(x^2), data = rocket)
1/(1-summary(fit1)$r.squared)

fit2 <- lm(I(x-mean(x))~I((x-mean(x))^2), data = rocket)
1/(1-summary(fit2)$r.squared)

## -----------------------------------------------------------------------------
## lm(output ~ I(1/wind), data = windmill, weight = c(w1, w2, [...], wn))


## -----------------------------------------------------------------------------
hardwood <- read.csv("hardwood.csv")
lmfit1 <- lm(strength ~ hardwood, data=hardwood)
lmfit2 <- lm(strength ~ hardwood + I(hardwood^2), data=hardwood)
lmfit3 <- lm(strength ~ hardwood + I(hardwood^2) + I(hardwood^3), data=hardwood)
lmfit4 <- lm(strength ~ hardwood + I(hardwood^2) + I(hardwood^3) + I(hardwood^4), data=hardwood)
lmfit5 <- lm(strength ~ poly(hardwood, 4, raw = TRUE), data=hardwood)
lmfit6 <- lm(strength ~ poly(hardwood, 4), data=hardwood)
# model.matrix(lmfit4) # matrix X
# model.matrix(lmfit5)
# model.matrix(lmfit6)
model.matrix(lmfit6) |> apply(2, \(col) sum(col^2))
data.frame(coefficients(lmfit4), coefficients(lmfit5), coefficients(lmfit6))
data.frame(predict(lmfit4), predict(lmfit5), predict(lmfit6))
plot(hardwood$hardwood, hardwood$strength, xlab="hardwood", ylab="strength")
lines(hardwood$hardwood, lmfit1$fitted.values, col=2, lwd=1.5)
lines(hardwood$hardwood, lmfit2$fitted.values, col=3, lwd=1.5)
lines(hardwood$hardwood, lmfit3$fitted.values, col=4, lwd=1.5)
lines(hardwood$hardwood, lmfit4$fitted.values, col=5, lwd=1.5)
legend("bottomright", c("1st order", "2nd order", "3rd order", "4th order"), col=c(2,3,4,5), lty=1)


## -----------------------------------------------------------------------------
library(readxl)
VD <- read_xls("Voltage_Drop.xls")
colnames(VD) <- c("obs","time","voltage")
head(VD)
hinge <- \(x,t) return((x - t)*(x > t)) # returns (x-t)_+
fit.pl <- lm(voltage ~ time + hinge(time, 6.5)  + hinge(time, 13), data=VD)
# model.matrix(fit.pl)
plot(voltage ~ time, data = VD)
lines(VD$time, fit.pl$fitted.values, col="blue", lwd=1.5)


## -----------------------------------------------------------------------------
lmfit.full <- lm(time ~ case + distance, data=deliver)
# ?influence.measures
e <- residuals(lmfit.full) # residuals
d <- residuals(lmfit.full)/summary(lmfit.full)$sigma # standardized residuals
r <- rstandard(lmfit.full) # *Studentized* residuals
t <- rstudent(lmfit.full)  # R-Student
press <- rstandard(lmfit.full, type="predictive") # PRESS residuals
data.frame(residuals=e,
           standardized=d,
           studentized=r,
           Rstudent=t,
           PRESS = press)


## -----------------------------------------------------------------------------
# ?plot.lm
plot(lmfit.full, which = 1)


## -----------------------------------------------------------------------------
plot(lmfit.full, which = 2)


## -----------------------------------------------------------------------------
plot(1:25, residuals(lmfit.full), 
     xlab = "Time", ylab = "Residuals")


## -----------------------------------------------------------------------------
#install.packages("lmtest")
library(lmtest)
dwtest(time ~ case + distance, data=deliver, alternative="two.sided")


## -----------------------------------------------------------------------------
x <- c(1.0, 1.0, 2.0, 3.3, 3.3, 4.0, 4.0, 4.0, 4.7, 5.0, 5.6, 5.6, 5.6, 6.0, 6.0, 6.5, 6.9)
y <- c(10.84, 9.30, 16.35, 22.88, 24.35, 24.56, 25.86, 29.16, 24.59, 22.25, 25.90, 27.20, 25.61, 25.45, 26.56, 21.03, 21.46)
plot(x, y)
anova(lm(y ~ 1), lm(y ~ x), lm(y ~ factor(x)))


## -----------------------------------------------------------------------------
hatvalues(lmfit.full) # h_ii
par(mfrow = c(1,2))
plot(hatvalues(lmfit.full))
abline(h = 2*3/nrow(deliver), lty = 2)
plot(deliver$case, deliver$distance, pch = NA)
text(deliver$case, deliver$distance, labels = 1:25)
which(hatvalues(lmfit.full) > 2*3/nrow(deliver))


## -----------------------------------------------------------------------------
cooks.distance(lmfit.full)
which(cooks.distance(lmfit.full) > 1) # Decision rule: D_i > 1 -> i is a influence point
dffits(lmfit.full) # |DFFITS| large -> influence point


## -----------------------------------------------------------------------------
library(readxl)
TL <- read_xls("Tool_Life.xls")
colnames(TL) <- c("obs","hours","rpm","type")
head(TL)
tail(TL)
TL$type <- factor(TL$type)
lmfit.TL <- lm(hours ~ rpm + type, data = TL)
summary(lmfit.TL)


## -----------------------------------------------------------------------------
TL$type |> levels()
TL$type2 <- factor(TL$type, levels = c("B", "A"))
lm(hours ~ rpm + type2, data = TL) |> summary()

plot(hours ~ rpm, data=TL, col=as.integer(type), 
     xlab="rpm", ylab="hours", pch=20)
abline(a=coefficients(lmfit.TL)[1], b=coefficients(lmfit.TL)[2], col=1, lwd=2)
abline(a=sum(coefficients(lmfit.TL)[c(1,3)]), 
       b=coefficients(lmfit.TL)[2], col=2, lwd=2)


## -----------------------------------------------------------------------------
# (1) and (2) are equivalent
lm(hours ~ rpm*type, data = TL) # (1)
lm(hours ~ rpm + type + rpm:type, data = TL) # (2)


## -----------------------------------------------------------------------------
## ### Example
## x1 <- c(1,1,0,0,-1,-1)
## x2 <- c(0,0,1,1,-1,-1)
## lm(y ~ x1 + x2) |> summary()

