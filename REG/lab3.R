## ---- Regression Analysis and Lab. --- Lab 3 ---------------------------------


## -----------------------------------------------------------------------------
acetylene <- read.csv("acetylene.csv")
colnames(acetylene) <- c("Obs", "P", "Temp", "H2", "Cont") 
lmfit <- lm(P ~ poly(scale(Temp), scale(H2), scale(Cont), degree = 2, raw = TRUE), 
            data = acetylene)
summary(lmfit)


## -----------------------------------------------------------------------------
X <- model.matrix(lmfit)
XPX <- cor(X[,-1]) 
XPX
X2 <- scale(X[,-1])/sqrt(15)
t(X2) %*% X2 ## equivalent to `cor(X[,-1])`.


## -----------------------------------------------------------------------------
# install.packages("rms")
library(rms)
vif(lmfit)


## -----------------------------------------------------------------------------
XPX |> solve() |> diag()


## -----------------------------------------------------------------------------
kappa(XPX, exact = TRUE)


## -----------------------------------------------------------------------------
e <- eigen(XPX)$values # eigenvalues
max(e) / min(e)


## -----------------------------------------------------------------------------
# install.packages("glmnet")
library(glmnet)
X <- model.matrix(lmfit)[,-1]
y <- acetylene$P
lmfit.ridge <- glmnet(x = X, y = y, alpha = 0) # `lambda`s are automatically chosen
lmfit.ridge |> coef()


## -----------------------------------------------------------------------------
lmfit.ridge$lambda # 100 lambda values are automatically chosen


## -----------------------------------------------------------------------------
# ridge trace
# ?plot.glmnet
plot(lmfit.ridge, xvar = "lambda", label = TRUE)
abline(v = log(lmfit.ridge$lambda[60]), lty = 2)


## -----------------------------------------------------------------------------
lmfit.ridge$lambda[60]
coef(lmfit.ridge)[, 60]


## -----------------------------------------------------------------------------
# ?predict.glm
# `s` does *not* have to be included in lmfit.ridge$lambda !
predict(lmfit.ridge, s = 45, type ="coefficients") # estimation of beta's when lambda = 45
predict(lmfit.ridge, s = 45, newx = X) # predictions at observed points with the model fitted at lambda = 45


## -----------------------------------------------------------------------------
hald <- read.csv("hald.csv")
lmfit <- lm(y ~ ., data = hald) # y ~ x1 + x2 + x3 + x4


## -----------------------------------------------------------------------------
summary(lmfit)$r.squared # R^2
summary(lmfit)$adj.r.squared # R^2_Adj
summary(lmfit)$sigma^2 # MSE
AIC(lmfit) # AIC
BIC(lmfit) # BIC


## -----------------------------------------------------------------------------
# install.packages("leaps")
library(leaps)
reg <- regsubsets(y ~ ., data = hald, method = "exhaustive")
reg.summary <- summary(reg)
data.frame(reg.summary$outmat,
           R2 = reg.summary$rsq,
           R2Adj = reg.summary$adjr2,
           SSE = reg.summary$rss,
           Cp = reg.summary$cp,
           BIC = reg.summary$bic)


## -----------------------------------------------------------------------------
coef(reg, 3) # estimated regression coefficients when k=3


## -----------------------------------------------------------------------------
# ?plot.regsubsets
# Plots a table of models showing which variables are in each model. 
# The models are ordered by the specified model selection statistic (`scale` argument).
plot(reg, scale="bic")
plot(reg, scale="Cp")
plot(reg, scale="adjr2")
plot(reg, scale="r2")


## -----------------------------------------------------------------------------
reg <- regsubsets(y ~ ., data = hald, method = "forward")
reg.summary <- summary(reg)
data.frame(reg.summary$outmat,
           R2 = reg.summary$rsq,
           R2Adj = reg.summary$adjr2,
           SSE = reg.summary$rss,
           Cp = reg.summary$cp,
           BIC = reg.summary$bic)
plot(reg) # different from the all possible regressions


## -----------------------------------------------------------------------------
reg <- regsubsets(y ~ ., data = hald, method = "backward")
reg.summary <- summary(reg)
data.frame(reg.summary$outmat,
           R2 = reg.summary$rsq,
           R2Adj = reg.summary$adjr2,
           SSE = reg.summary$rss,
           Cp = reg.summary$cp,
           BIC = reg.summary$bic)
plot(reg) 


## -----------------------------------------------------------------------------
library(readxl)
Puromycin = read_xls("Puromycin.xls")
plot(y ~ x, data = Puromycin)


## -----------------------------------------------------------------------------
lmfit <- lm(1/y ~ I(1/x), data = Puromycin)
par(mfrow = c(1,2))
plot(1/y ~ I(1/x), data = Puromycin)
abline(lmfit)
plot(y ~ x, data = Puromycin)
lines(Puromycin$x, 1/fitted.values(lmfit))


## -----------------------------------------------------------------------------
lmfit <- nls(y ~ I(t1*x/(t2+x)), data = Puromycin, 
             start = list(t1=205, t2=0.08))
summary(lmfit)
plot(y ~ x, data = Puromycin)
lines(Puromycin$x, fitted.values(lmfit))


## -----------------------------------------------------------------------------
pneumoconiosis <- read.csv("pneumoconiosis.csv")
colnames(pneumoconiosis) <- c("x","pneum","miners","y")
glmfit <- glm(cbind(pneum, miners-pneum) ~ x, family = binomial, data = pneumoconiosis)
summary(glmfit)


## -----------------------------------------------------------------------------
success <- sum(pneumoconiosis$pneum)
failure <- sum(pneumoconiosis$miners) - success
y.new <- rep(c(1,0),c(success, failure))
x.new <- c(rep(pneumoconiosis$x, pneumoconiosis$pneum),
           rep(pneumoconiosis$x, pneumoconiosis$miners - pneumoconiosis$pneum))
glm(y.new ~ x.new, family = binomial, data = pneumoconiosis)


## -----------------------------------------------------------------------------
glm(y ~ x, family = binomial, weights = miners, data = pneumoconiosis)


## -----------------------------------------------------------------------------
coef(glmfit)
fitted.values(glmfit)
residuals(glmfit)


## -----------------------------------------------------------------------------
predict(glmfit)
model.matrix(glmfit) %*% coef(glmfit) 
predict(glmfit, type = "response")

