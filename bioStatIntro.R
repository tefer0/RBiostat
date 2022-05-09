#data structures-----
dbl_var <- c(1,2.5,4.5)
int_var <- c(1L,6L,10L)
log_var <- c(TRUE,FALSE,T,F)
chr_var <- c('These are','some strings')

dbl_var

int_var

log_var

chr_var

a <- c(1,c(2,c(3,4)))

a
b <- c(1,2,3,4)

typeof(a)
typeof(int_var)
typeof(dbl_var)
typeof(chr_var)

is.integer(int_var)

d <- str(c('a',1))

d

x <- c(FALSE,FALSE,TRUE)
as.numeric(x)
sum(x)
mean(x)
#lists-----
x1 <- list(1:3,'a',c(TRUE,FALSE,TRUE),c(2.3,5.9))
x1

x2 <- list(list(list(list)))

str(x1)
x2
str(x2)
is.recursive(x2)

x3 <- list(list(1,2),c(3,4))
y <- c(list(1,2),c(3,4))

str(x3)
str(y)

sum(x3)
sum(y)

#mtcars dataset----
mtcars
is.list(mtcars)
is.data.frame(mtcars)

#linear model-----
#independent variables predict(explanatory) dependant variables
#dependent variables are predicted(respond)
mod <- lm(mpg~wt, data = mtcars)
mod
is.list(mod)

#attributes-----
y1 <- 1:10
y1
attr(y1,"my_attribute") <- "This is a vector"
attr(y1,"my_attribute")
str(attributes(y1))
attributes(y1)

#atributes are lost when a vector is modified
attributes(y1[1])
attributes(sum(y1))

#NULL and NA ------
y3 <- c(a=1,2,3)
names(y3)
v <- c(1,2,3)
names(v) <- c('a','b','c')
names(v)
z <- c(1,2,3)
names(z)
#all names missing NULL
#NA not applicable (not available yet)

#factors-----
x4 <- factor(c('a','b','b','a'))
class(x4)
levels(x4)
x4

x4[2] <- 'c'
#[,] for subsetting before comma row, then after col, (for assignment)

c(factor('a'),factor('b'))
#R<4.1 gives 1 1 for true true, R>4.0 gives the factors a,b and levels

sex_char <- c('m','m','m')
sex_factor <- factor(sex_char,levels=c('m','f'))

levels(sex_char)
sex_factor

#tables----
table(sex_char)

table(sex_factor)

is.data.frame(table(sex_factor))
is.data.frame(table(sex_char))
#table and dataframe are different algorithms hence FALSE

#csv-----
z1 <- read.csv(text='value\n12\n1\n.\n9')
z1

typeof(z1$value)

as.double(z1$value)

class(z1$value)

as.double(as.character(z1$value))

#matrix----
a1 <- matrix(1:6, ncol=3,nrow=2)
a1          

#array----
b1 <- array(1:12,c(2,3,2))
b1

#dimmensions----
xy <- c(1:6)
dim(xy) <- c(3,2);xy
dim(xy) <- c(2,3);xy

xy[,1];xy

a1
length(xy)
nrow(xy)
ncol(xy)

rownames(xy) <- c('A','B')
rownames(xy)
xy
colnames(xy) <- c('A','B','C')
colnames(xy)

length(xy)
b1
length(b1)
dim(b1)

dimnames(b1) <- list(c('one','two'),c('a','b','c'),c('A','B'))

#dataframe-----
df <- data.frame(x=1:3, y=c('a','b','c'), stringsAsFactors = FALSE)
str(df)
typeof(df)
class(df)
is.data.frame(df)
is.table(df)

cbind(df, data.frame(z=3:1)) #cbind column bind

rbind(df, data.frame(x=10, y='z'))

gdata <- data.frame(cbind(a=1:2, b=c('a','b')))
str(gdata)

good <- data.frame(a=1:2, b=c('a','b'), stringsAsFactors = FALSE)
str(good)

df1 <- data.frame(x=1:3)
str(df1)

df1$y <- list(1:2,1:3,1:4) #error dif row nos
data.frame(df1$y)

df2 <- data.frame(x=1:3, y=list(1:2,1:3,1:4)) #error diff row nos

df3 <- data.frame(x=1:3, y=I(list(1:2,1:3,1:4))) #I for AsIs
str(df3)

dfm <- data.frame(x=1:3, y=I(matrix(1:9,nrow=3))) 
str(dfm)
dfm[2,'y']                  

#in built fxns----
mean()
max()
sum()
colMeans()
rowSums()
cummax()
pmin()

#probability fxns----
#normal, poisson, binomial, factorial, etc

#1. plots----
par(mfrow=c(1,1)) #par for drawing parameters
x <- 0:6
plot(x,factorial(x),type='s',main = 'factorial x', log='y')
#main is title, type is label, 's' for standard

?par(mfrow)
#for data visualisation

#2.choose----
choose(8,3)
#out of 8 choose 3, 56 ways of choosing, binomial 
plot(0:8,choose(8,0:8),type='s',main='binomial coefficients')

#continous probability----
#bionomial, chisq, fishers f, gamma,exponential,students t test, normal,
#poisson 
binom.test()
chisq.test()
gamma()
exp()
t.test()
norm()
poisson.test()

curve(pnorm(x,-3,3))
#curve of probability of normal distribution

arrows(-1,0,1,pnorm(-1),col = 'red') #col for color, pnorm probability 
arrows(-1,pnorm(-1),-3,pnorm(-1),col = 'green')
?arrows

curve(dnorm(x),-3,3) #dnorm density normal

par(mfrow=c(2,2)) #2 * 2 no of graphs on page
x <- seq(-3,3,0.01) #seq sequence generation
y <- exp(-abs(x)) #exp exponential, abs absolute nos
plot(x,y,type = 'l', main = 'x', col = 'green') #type l for line

y <- exp(-abs(x)^2)
plot(x,y,type = 'l', main = 'x^2', col = 'yellow')

y <- exp(-abs(x)^3)
plot(x,y,type = 'l', main = 'x^3', col = 'blue')

y <- exp(-abs(x)^8)
plot(x,y,type = 'l', main = 'x^8', col = 'red')

#probability density----
pnorm(-1.25)

pnorm(1.875)

1-pnorm(1.875)

pnorm(1.25)-pnorm(0.625)

x <- seq(-3,3,0.01)
z <- seq(-3,-1.25,0.01)
p <- dnorm(z)
z <- c(z,-1.25,-3)
p <- c(p,min(p),min(p))
plot(x,dnorm(x),type = 'l',xaxt='n',ylab='probability density',xlab = 'height')
axis(1,at=-3:3,labels = c('146','154','162','170','178','186','192'))
polygon(z,p,col='red') #polygon to insert region in curve, axis to add x axis

z <- seq(-0.635,1.25,0.01)
p <- dnorm(z)
z <- c(z,1.25,-0.635)
p <- c(p,0,0)
plot(x,dnorm(x))
axis(1,at=-3:3,labels = c('146','154','162','170','178','186','192'))
polygon(z,p,col='red')

#uniformly distrbuted----
?par(mfrow)

par(mfrom=c(1,1))
hist(runif(10000)*10,main = 'title') #10000 random numbers repeated 10 times

xv <- seq(0,10,0.1) #sequence of numbers from 0 to 10, probability 0.1
yv <- dnorm(xv,mean = 4.998581,sd=1.28996)*5000
lines(xv,yv) #plot line on histogram

par(mfrow=c(2,2))
hist(sample(1:6,replace = T,10000),breaks = 0.5:6.5,main = 'title',xlab = 'one die')
#replace means you replace the die face after casting once to avoid same no

a <- sample(1:6,replace = T,10000)
b <- sample(1:6,replace = T,10000)
hist(a+b,breaks = 1.5:12.5,main = 'title2',xlab = 'two die')
#get triangular disbn centered around 7

c <- sample(1:6,replace = T,10000)
hist(a+b+c,breaks = 2.5:18.5,main = 'title3',xlab = 'three die')
#get belly shaped disbn (binomial)

d <- sample(1:6,replace = T,10000)
e <- sample(1:6,replace = T,10000)
hist(a+b+c+d+e,breaks = 4.5:30.5,main = 'title4',xlab = 'five die')
#get normal disbn as no of die increases

mean(a+b+c+d+e)
sd(a+b+c+d+e)
lines(seq(1,30,0.1),dnorm(seq(1,30,0.1),17.5937,3.837668)*10000)

#comparing data with a normal disbn----
#shapiro-wilks test and Kolmogorov–Smirnov test to test for normality

#par(mfrow=c(1,1))
#importing data from file
cars1 <- read.table('cars',header = T)
attach(cars1)
names(cars1)
mean(speed)
max(cars$speed)
min(speed)
max(dist)
min(dist)
var(cars)
var(speed)
length(speed)
length(dist)
plot(cars1)
#positive correlation graph
par(mfrow=c(1,2))
hist(speed,breaks = -0.5:25.5,col = 'green',main = 'cars1')
hist(speed,breaks = max(speed),min(speed),col = 'blue',main = 'cars2')
hist(dist,breaks = max(dist),min(dist),col = 'red',main = 'cars3')
#max before min
lines(seq(0,26,1),length(speed)*dnorm(seq(0,26,1),mean(speed),sqrt(var(speed))))
lines(seq(02,120,10),length(dist)*dnorm(seq(02,120,10),mean(dist),sqrt(var(dist))))

shapiro.test(cars$speed)
shapiro.test(cars$dist)

#data transformation----
par(mfrow=c(1,1))
fishes <- read.table('fishes.txt', header = T)
attach(fishes)
names(fishes)
mean(mass)
max(mass)
min(mass)
hist(mass,breaks = -0.5:16.5, col = 'red',main = 'fishes')
lines(seq(0,16,0.1),length(mass)*dnorm(seq(0,16,0.1),mean(mass),sqrt(var(mass))))

shapiro.test(mass)
#wald test for normality null: normally distributed
#alternative: not normally distributed
#when p<alpha we reject null, hence data not normally distributed.
#try to normalise data. by transformation using log
logfishes <- log10(mass)
par(mfrow=c(1,1))
max(logfishes)
min(logfishes)
hist(logfishes,breaks = max(logfishes),min(logfishes), col = 'blue', main = 'new')
lines(log10(seq(0,16,0.1)),length(logfishes)*dnorm(log10(seq(0,16,0.1)),mean(logfishes),sqrt(var(logfishes))))

shapiro.test(logfishes)
#p>0.05, we fail to reject null, log10 data is normally distributed
#we can now do a (one sample) simple (student) t-test
#for two samples we use two sample t test (mean)
#for three samples or more we use f test (ANOVA)
#for logical data, binomial tests like chi2

t.test(logfishes)
#if data is not normally distributed we use wilkocson test(non-parametric)
#for one sample
#alt ANOVA is Kruscawallis test

#maximum likelihood normal dist----
#the z test for std normal dist
par(mfrow=c(2,2))
curve(dnorm,-3,3,xlab = 'z',ylab = 'probability density',main='density')
#for density distribution, area under the curve
curve(pnorm,-3,3,xlab = 'z',ylab = 'probability', main='probability')
#for probability distribution s shaped
curve(qnorm,0.1,xlab = 'p',ylab = 'quantile (z)', main='quantiles')
#for quantile values
y <- rnorm(1000)

#1000 random numbers, mean = 0, sd = 1 
hist(y,xlab = 'z',ylab = 'frequency',main = 'random numbers')

#test stat vs critical stat----
1-pchisq(14.3,9) #calculated value

qchisq(0.95,9) #q for quatile probabilities at 95% confidence interval, df=9.

1-pf(2.85,8,12) #probability value
#f distbn, var 2.85, df numerator 8, df deno 12

qt(0.975,10) #qt quantile value
#students t test, 2 tailed 5%/2 = 0.025

#chi2,x disb----
#mean = df
par(mfrow=c(1,2))
x <- seq(0,30,0.25)
plot(x,pchisq(x,3,7.25), type = 'l',ylab='p(x)',xlab = 'x')

plot(x,pchisq(x,5,10),type = 'l',ylab = 'p(x)',xlab = 'x')

#var = 10.2 , df =8, interval alpha2 quantiles, area under curve
8*10.2/qchisq(0.975,8)
8*10.2/qchisq(0.025,8)

#fisher's test----
qf(0.95,2,18)
x <- seq(0.05,4,0.05)
#df2,df18
plot(x,df(x,2,18),type = 'l',ylab = 'f(x)',xlab = 'x')
plot(x,df(x,6,18),type = 'l',ylab = 'f(x)',xlab = 'x')
#plot 2 increase in sample size

par(mfrow=c(1,1))
df <- seq(1,30,0.1)
plot(df,qf(0.95,df,30),type = 'l',ylab = 'critical F',)
lines(df,qf(0.95,df,10),lty=2) #lty - line type 2 broken line

x <- seq(0.01,3,0.01)
plot(x,df(x,1,10),type = 'l',ylim = c(0,1),ylab = 'f(x)')
lines(x,df(x,2,10),lty=6,col='red')
lines(x,df(x,5,10),lty=2,col='green')
lines(x,df(x,30,10),lty=3,col='blue')
legend(2,0.9,c('1','2','5','30'),col=(1:4),lty=c(1,6,2,3),title = 'numerator d.f')

#students t distribution----
curve((1+x^2)^(-0.5),-3,3,ylab = 't(x)',col='red')

plot(1:30,qt(0.975,1:30),ylim=c(0,12),type='l',ylab='student
t value',xlab='d.f',col='red')
abline(h=2,lty=2,col='green')

xvs <- seq(-4,4,0.01)
plot(xvs,dnorm(xvs),type='l',lty=26,ylab='prob density',xlab='deviates',col='red')
lines(xvs,dt(xvs,df=5),col='blue') #density of t disbn
#two tail alpha/2 hence 0.05/2. so 0.975
qt(0.975,5)

#data summary ----
data <- read.table('temp.txt',header = T)
attach(data)
names(data)
par(mfrow=c(1,1))
plot(temps)
boxplot(temps)
hist(temps, main='')
summary(temps)

par(mfrow=c(1,1))
qqnorm(temps) #quatile fxn
qqline(temps,lty=2) #should tuch head and tail
#our data is not normally distributed

x <- exp(rnorm(30))
shapiro.test(x)
shapiro.test(temps)
#pvalue < 0.05, non parametric, not normally distributed

log(temps)
par(mfrow=c(1,1))
qqnorm(log(temps))
qqline(log(temps),lty=2)
shapiro.test(log(temps))

#after log transformation data still not normaly disbt
#hence we use non-parametric alt test
wilcox.test(temps)
#V is test statistic
#H0: equal, H1: not equal, p<&, reject null
#observations highly significant difference, V=12246, p<&

#one sample, H0: data points equal
#2 sample, H0: sample means equal
#more than 2, H0: sample variances equal

#2 sample -----
fdata <- read.table('f.test.data.txt',header = T)
attach(fdata)
names(fdata)
#assumptions variances are equal
var(gardenB)
var(gardenC)
#t.test(gardenB,gardenC) two sample t.test
f.ratio <- var(gardenC)/var(gardenB)
f.ratio
#critical f=4.026, cal f=10.66667, reject H0.
#because critical F < calculated F
2*(1-pf(f.ratio,9,9))
#two tailed nature of the test
#shortcut. f test compare two vars----
var.test(gardenB,gardenC)

#t test----
ttest.data <- read.table('gardens.txt',header = T)
attach(ttest.data)
names(ttest.data)

effect <- c(gardenA, gardenB)
label <- factor(c(rep('A',10),rep('B',10)))
boxplot(effect~label,notch=F,xlab='Garden',ylab = 'ozone',col='Red')

s2A <- var(gardenA)
s2B <- var(gardenB)
(mean(gardenA)-mean(gardenB))/sqrt(s2A/10+s2B/10)
qt(0.975,18)
#df = (10-1) + (10-1)
t.test(gardenA,gardenB)

#wilcoxon----
ozone <- c(gardenA,gardenB)
label <- c(rep('A',10),rep('B',10))
combined.ranks <- rank(ozone)
combined.ranks

wilcox.test(gardenA,gardenB)

#paired samples ----
stream <- read.table('streams.txt',header = T)
attach(stream)
names(stream)
t.test(down,up)
#use paired=True
t.test(down,up,paired = TRUE)
#now gives significant pvalue
diff <- down-up
t.test(diff)
#this also works, diff shows correlation

#signed test----
binom.test(1,9)

sign.test <- function(x,y) 
{
if(length(x)!=length(y))stop('The two variables must be same lenght')
d <- x-y
binom.test(sum(d>0),length(d))
}
#sign.test(gardenA,gardenB)

#binomial ----
#vec1 success f,m
#vec2 total f,m
prop.test(c(4,196),c(40,3270))
#pvalue > alpha hence no evidence for positive discrimination

#chi-sq test----
qchisq(0.95,1)

count <- matrix(c(38,14,11,51),nrow = 2)
count
chisq.test(count)
chisq.test(count,correct = FALSE)
#more accurate
chisq.test(count,correct=F)$expected

chisq.test(c(10,3,2,6))
chisq.test(c(10,3,2,6),p=c(0.2,0.2,0.3,0.3))

die <- ceiling(runif(100,0,6))
table(die)

#contigency tables, fishers exact test----
#small expected freq
factorial(8)*factorial(12)*factorial(10)*factorial(10)/(factorial(6)*factorial(2)*factorial(4)*factorial(8)*factorial(20))
x <- as.matrix(c(6,4,2,8))
dim(x) <- c(2,2)

table <- read.table('fisher.txt', header=TRUE)
attach(table)
names(table)
head(table)
fisher.test(tree,nests)
#uses odds ratio

data <- read.table('twosample.txt',header = T)
attach(data)
names(data)
plot(x,y,pch=21,col='Red',bg='orange')
plot(data$x,data$y,pch=21,col='Red',bg='orange')
var(data$x)
var(data$y)
cov(data$x,data$y)
var(data$x,data$y)
var(data$x,data$y)/sqrt(var(data$x)*var(data$y))
cor(data$x,data$y)
