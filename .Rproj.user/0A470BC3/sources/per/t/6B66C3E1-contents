#data modelling-----
#influence

x <- c(2,3,3,3,4)
y <- c(2,3,2,1,2)
par(mfrow=c(1,2))
plot(x,y,xlim=c(0,8),ylim=c(0,8))

#outlier influence
x1 <- c(x,7)
y1 <- c(y,6)
plot(x1,y1,xlim=c(0,8),ylim=c(0,8))
abline(lm(y1~x1),col='blue')

reg <- lm(y1~x1)
summary(reg)

influence.measures(reg)
#outlier is influential since it has a star

influence.measures(reg)$is.inf

#OR 
lm.influence(reg):
  lm.influence(reg)

lm(y1[-6]~x1[-6])
summary.aov(lm(y1[-6]~x1[-6]))
#summary of anova

data <- read.table("regression.txt",header=T)
attach(data)
names(data)
model <-lm(growth~tannin)

logLik(model) 

AIC(model) 

#influence.measures(model)

graze <- read.table("ipomopsis.txt",header=T)
attach(graze)
names(graze)
#AIC as a measure of the fit of a model 
model.1 <- lm(Fruit~Grazing*Root)
model.2 <- lm(Fruit~Grazing+Root)

AIC(model.1, model.2)

#Model checking in R ----
decay <- read.table("Decay.txt",header=T)
attach(decay) 
names(decay)
model <- lm(amount~time) 

par(mfrow=c(1,2))
plot(model) 

#Extracting information from model objects
#We often want to extract material from fitted models (e.g. slopes, residuals or p values) and there are three different ways of doing this:
 # by name, e.g. coef(model);
model <- lm(growth~tannin)
summary(model)
coef(model)# gives the parameter estimates (‘coefficients’) for the intercept (a) and slope (b)
fitted(model)# gives the fitted values (yˆ = a + bx) of the model (its predictions) for each value of the explanatory variables
resid(model)# gives the residuals (y minus fitted values) for the nine data points
vcov(model)# extracts the variance–covariance matrix of the model’s parameters

#with list subscripts, e.g. summary(model)[[3]];
#The two model summary objects summary.aov(model) and summary.lm(model) are lists 
#With many components
summary.aov(model)# gives columns of the ANOVA table
summary.aov(model)[[1]][1] # can be used to extracted one column at a time
#Here is summary.lm: which can follow similar application trend
#summary(model)
#using $ to name the component, e.g. model$resid. 
#To extract model components one can use the $ symbol
model$coef

comp <- read.table("competition.txt",header=T)
attach(comp)
names(comp) 

model <- lm(biomass~clipping) 
summary.lm(model)
summary.aov(model) 

means <- tapply(biomass,clipping,mean)
means

#Scale-dependent correlations----
prodct<-read.table("productivity.txt", header=TRUE)
attach(prodct)
names(prodct)
head(prodct)
plot(prodct$x, prodct$y, pch=21, col="blue", bg="green",  xlab="Productivity", ylab="Mammal species")
cor.test(prodct$x,prodct$y)

formula<-prodct$x ~ prodct$y | prodct$f
coplot(formula, prodct,col = 'blue')

#Bootstrap-----
data <- read.table("skewdata.txt",header=T)
attach(data)
names(data)
head(data)

sample(values,replace=T)

ms <- numeric(10000)
for (i in 1:10000){
  ms[i] <- mean(sample(values,replace=T))}

quantile(ms,c(0.025,0.975))

mean(values)-quantile(ms,c(0.025,0.975))

1.96*sqrt(var(values)/length(values))

library(boot)
#boot(data, statistic, R)
mymean <- function(values,i) mean(values[i])
#The key point is that we write mean(values[i]) not mean(values)
#Now we can run the bootstrap for 10 000 iterations
myboot <- boot(values,mymean,R=10000)
myboot
#boot confidence interval
boot.ci(myboot)

#PCA----
#kmd <- read.table("kmeansdata.txt",header=T)
#Partitioning
kmd <- read.table("kmeansdata.txt",header=T)
attach(kmd)
names(kmd)
par(mfrow=c(2,2))
plot(kmd$x,kmd$y,pch=16)

plot(kmd$x,kmd$y,col=group,pch=16)
model <- kmeans(data.frame(x,y),6)
plot(kmd$x,kmd$y,col=model[[1]])
model <- kmeans(data.frame(x,y),4)
plot(kmd$x,kmd$y,col=model[[1]])
par(mfrow=c(1,1))

#Hierarchical
taxa <- read.table("taxon.txt",header=T)
attach(taxa)
names(taxa)
pairs(taxa)
kmeans(taxa,4)

plot(hclust(dist(taxa)),main="")

library(MASS)
model <- lda(Taxon~.,taxa)
plot(model,col=rep(1:4,each=30))
