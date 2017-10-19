library(rpart)
require(broom)

#wine <- read.csv("C:/Users/kcaudle/Documents/SDSMT/Teaching/Spring 2017/Brenna Mollet/DelayedGreedy/winequality.csv",header=T)
attach(wine)

gmodel <- rpart(quality~.,data=wine, maxcompete=0, usesurrogate=0)
cpvalg <- gmodel$cptable[which.min(gmodel$cptable[,"xerror"]),"CP"]
gmodel <- prune(gmodel,cp=cpvalg)
SSg <- tidy(gmodel$cptable)
r2g <- 1-tail(SSg$rel.error,n=1)

dgmodel <- rpart(quality~.,data=wine,delayed=1, maxcompete=0) 
cpvaldg <- dgmodel$cptable[which.min(dgmodel$cptable[,"xerror"]),"CP"]
dgmodel <- prune(dgmodel,cp=cpvaldg)
SSdg <- tidy(dgmodel$cptable)
r2dg <- 1-tail(SSdg$rel.error,n=1)

r2g
r2dg
