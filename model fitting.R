## read in wine data files ###
rw <- read.csv('data/winequality-red.csv', header = T, sep = ';')
ww <- read.csv('data/winequality-white.csv', header = T, sep = ';')
## combine datasets
rw$type <- 'red'
ww$type <- 'white'
aw <- data.frame(rbind(rw, ww)) ## create full dataset with all wines

## look at quality distributions
table(rw$quality)
table(ww$quality)
table(aw$quality)
## histograms
hist(rw$quality, breaks = 0:10)
hist(ww$quality, breaks = 0:10)
## test normality
shapiro.test(rw$quality)
shapiro.test(ww$quality)

## look at correlations between numeric predictors

cor(ww)
pairs(ww)
pairs(rw)
names(aw)
par(mar = c(8,4,4,4))

round(cor(aw[-c(length(names(aw)), length(names(aw))-1)]),2)

heatmap(cor(aw[-c(length(names(aw)), length(names(aw))-1)]),
        scale = 'column',
        col = cm.colors(256),
        Colv = NA, Rowv = NA)

summary(ww)
summary(rw)

## create dummy vars for quality
mm <- data.frame(model.matrix(~ factor(quality) + 0, aw))
aw2 <- cbind(aw, mm)
names(aw2)
cols <- cbind(aw2$factor.quality.3, aw2$factor.quality.4, aw2$factor.quality.5, 
              aw2$factor.quality.6, aw2$factor.quality.7, aw2$factor.quality.8, 
              aw2$factor.quality.9)
## fit cumulative logit model with proportional odds
library(VGAM)
fit.po <- vglm(cols ~ fixed.acidity + volatile.acidity + citric.acid +
                 residual.sugar + chlorides + free.sulfur.dioxide +
                 total.sulfur.dioxide + density + pH +
                 sulphates + alcohol + factor(type),
               family=cumulative(parallel=TRUE), data=aw2)


fit.po.sc <- vglm(cols ~ scale(fixed.acidity) + scale(volatile.acidity) + scale(citric.acid) +
                 scale(residual.sugar) + scale(chlorides) + scale(free.sulfur.dioxide) +
                 scale(total.sulfur.dioxide) + scale(density) + scale(pH) +
                 scale(sulphates) + scale(alcohol) + factor(type),
                family=cumulative(parallel=TRUE), data=aw2)

fit.po2 <- vglm(cols ~ fixed.acidity + volatile.acidity +
                   residual.sugar + free.sulfur.dioxide +
                   total.sulfur.dioxide + density + pH +
                   sulphates + alcohol + factor(type),
            family=cumulative(parallel=TRUE), data=aw2)
summary(fit.po)
summary(fit.po.sc)

pchisq(deviance(fit.po2)-deviance(fit.po),df=df.residual(fit.po2)-df.residual(fit.po),lower.tail=FALSE)
## diagnostic check against model without proportional odds
fit.npo <- vglm(cols ~ scale(fixed.acidity) + scale(volatile.acidity) + scale(citric.acid) +
                  scale(residual.sugar) + scale(chlorides) + scale(free.sulfur.dioxide) +
                  scale(total.sulfur.dioxide) + scale(density) + scale(pH) +
                  scale(sulphates) + scale(alcohol),
            family=cumulative(), data=aw2) ## model returns error, probably because there are too many parameters
summary(fit.npo) 
pchisq(deviance(fit.po)-deviance(fit.npo),df=df.residual(fit.po)-df.residual(fit.npo),lower.tail=FALSE)

## fit log-log link model
fit.ll <- vglm(cols ~ fixed.acidity + volatile.acidity + citric.acid +
                 residual.sugar + chlorides + free.sulfur.dioxide +
                 total.sulfur.dioxide + density + pH +
                 sulphates + alcohol,
               family=cumulative(link = cloglog,parallel=TRUE), data=aw2)
summary(fit.ll)


## more complex models are unable to fit
## change to three categories y to simplify

aw$quality.bucket <- ifelse(aw$quality <= 4, 'Low',
                            ifelse(aw$quality <= 6, 'Med', 'High'))
table(aw$quality.bucket)
mm <- data.frame(model.matrix(~ factor(quality.bucket) + 0, aw))
aw2 <- cbind(aw, mm)
names(aw2)
cols <- cbind(aw2$factor.quality.bucket.Low,
              aw2$factor.quality.bucket.Med,
              aw2$factor.quality.bucket.High)

fit.po <- vglm(cols ~ fixed.acidity + volatile.acidity + citric.acid +
                 residual.sugar + chlorides + free.sulfur.dioxide +
                 total.sulfur.dioxide + density + pH +
                 sulphates + alcohol,
               family=cumulative(parallel=TRUE), data=aw2)
summary(fit.po)
## diagnostic check against model without proportional odds
fit.npo <- vglm(cols ~ fixed.acidity + volatile.acidity + citric.acid +
                  residual.sugar + chlorides + free.sulfur.dioxide +
                  total.sulfur.dioxide + density + pH +
                  sulphates + alcohol,
                family=cumulative(link = probit), data=aw2) ## model returns error, probably because there are too many parameters
summary(fit.npo) 

## fit log-log link model
fit.ll <- vglm(cols ~ fixed.acidity + volatile.acidity + citric.acid +
                 residual.sugar + chlorides + free.sulfur.dioxide +
                 total.sulfur.dioxide + density + pH +
                 sulphates + alcohol,
               family=cumulative(link = cloglog,parallel=TRUE), data=aw2)
summary(fit.ll)
?vglm
## diagnostic check comparing link functions by deviance, both have the same df
deviance(fit.po)
deviance(fit.ll)
deviance(fit.npo)


##Harry's code
##Compute the probability of Y=3 using cumulative logit model
rw <- winequality_red
ww <- winequality_white
rw$type <- 'red'
ww$type <- 'white'
aw <- data.frame(rbind(rw, ww))
par(mar = c(8,4,4,4))
round(cor(aw[-c(length(names(aw)), length(names(aw))-1)]),2)
mm <- data.frame(model.matrix(~ factor(quality) + 0, aw))
aw2 <- cbind(aw, mm)
cols <- cbind(aw2$factor.quality.3, aw2$factor.quality.4, aw2$factor.quality.5, 
              aw2$factor.quality.6, aw2$factor.quality.7, aw2$factor.quality.8, 
              aw2$factor.quality.9)
library(VGAM)
fit.po <- vglm(cols ~ fixed.acidity + volatile.acidity + 
                 residual.sugar + chlorides + free.sulfur.dioxide +
                 total.sulfur.dioxide + density + pH +
                 sulphates + alcohol + factor(type),
               family=cumulative(parallel=TRUE), data=aw2)
##Extract the coefficients of the cumulative logit model
betahat<-coef(fit.po)[7:17]
##Transform the list of estimated coefficients into a vector
names(betahat)<-NULL
##Combine the datasets
rw1<-winequality_red
ww1<-winequality_white
aw1<-data.frame(rbind(rw1,ww1)) 
aw1$quality<-NULL
aw1$citric.acid<-NULL
##Compute column means
cms<-colMeans(aw1)
names(cms)<-NULL
cms[11]<-0
##Use the formula to calculate the probability
sum(betahat*cms)+coef(fit.po)[1]
exp(sum(betahat*cms)+coef(fit.po)[1])
prob3<-exp(sum(betahat*cms)+coef(fit.po)[1])/(1+exp(sum(betahat*cms)+coef(fit.po)[1]))

##Conduct two-sample t-test to test the 
##difference between the means of these two types of wine
rw <- winequality_red
ww <- winequality_white
##Sample sizes
nr<-length(rw$quality)
nw<-length(ww$quality)
##Compute the vecotr of the sample mean of each column
meanr<-colMeans(rw)
meanw<-colMeans(ww)
##Compute the vecotr of the sample variance of each column
varr<-sapply(rw, var)
varw<-sapply(ww, var)
##Compute the vecotr of the statistic of each column
t<-(meanr-meanw)/sqrt(varr/nr+varw/nw)
##Define the significance level
alpha<-0.05
##Degrees of freedom for each explanatory variable
dfv<-(varr/nr+varw/nw)^2/((varr/nr)^2/(nr-1)+(varw/nw)^2/(nw-1))
##Compute the vector of the critical value of each explanatory variable
talpha<-qt(1-alpha/2, df=dfv)
##Compare the statistics with the critical values
abs(t)>talpha


##Build a separate model for each dataset
rw <- winequality_red
ww <- winequality_white

par(mar = c(8,4,4,4))
round(cor(rw[-c(length(names(rw)), length(names(rw))-1)]),2)
round(cor(ww[-c(length(names(ww)), length(names(ww))-1)]),2)
##Some necessary changes to the datasets
mmr <- data.frame(model.matrix(~ factor(quality) + 0, rw))
mmw <- data.frame(model.matrix(~ factor(quality) + 0, ww))
rw2 <- cbind(rw, mmr)
ww2 <- cbind(ww,mmw)
rw2<-data.frame(rw2)
ww2<-data.frame(ww2)

rcols <- cbind(rw2$factor.quality.3, rw2$factor.quality.4, rw2$factor.quality.5, 
              rw2$factor.quality.6, rw2$factor.quality.7, rw2$factor.quality.8)
wcols<-cbind(ww2$factor.quality.3, ww2$factor.quality.4, ww2$factor.quality.5, 
             ww2$factor.quality.6, ww2$factor.quality.7, ww2$factor.quality.8,ww2$factor.quality.9)
library(VGAM)
##Fit the model for the red wine
rfit.po <- vglm(rcols ~ fixed.acidity + volatile.acidity + citric.acid +
                 residual.sugar + chlorides + free.sulfur.dioxide +
                 total.sulfur.dioxide + density + pH +
                 sulphates + alcohol,
               family=cumulative(parallel=TRUE), data=rw2)

##Fit the model for the white wine
wfit.po<-vglm(wcols ~ fixed.acidity + volatile.acidity + citric.acid +
                residual.sugar + chlorides + free.sulfur.dioxide +
                total.sulfur.dioxide + density + pH +
                sulphates + alcohol,
              family=cumulative(parallel=TRUE), data=ww2)
summary(rfit.po)
summary(wfit.po)
