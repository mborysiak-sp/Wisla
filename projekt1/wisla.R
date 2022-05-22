library(tidyr)
library(gamlss)
library(maps)
library(fitdistrplus)

coords <- read.csv("coord233_alt.csv")
city <- coords[coords$place=="WISLA",]
station_id = city['station']$station

pdf("mapa.pdf")
poland <- map('world', 'poland', fill=T, col='gray');
points(city[c('lon', 'lat')], pch=19, col=2);
dev.off()

reload_file <- FALSE

path_to_files <- "temp.stations-all"
L <- as.list(list.files(path=path_to_files, pattern = ".*[6-8].csv"))
L <- paste0(path_to_files, "//", L)
data0 <- lapply(L, read.csv)

n <- length(data0); n
i <- 1
data0[[1]]
x <- data0[[i]]$X249180230 
head(x)
max10 <- c()
datetime <- c()

for(i in i:n) {
  max10 <- c(max10, data0[[i]]$X249180230)
  datetime <- c(datetime, as.character(data0[[i]]$datetime))
}

x <- c()
for(i in 1:n) {
  sti <- colnames(data0[[i]])
  x[i] <- "X249180230" %in% sti
}

sum(x)
which(x==TRUE)

library(tidyr)
max10 <- data.frame(date = as.Date(datetime), max10=max10)
head(max10)
tail(max10)
rownames(max10) <- c()

max10 <- separate(max10, date, c("year", "mth", "day"), convert=TRUE)
head(max10)

max10 <- data.frame(datetime=datetime, max10)

library(gamlss)
t1 <- Sys.time()
fit <- fitDist(max10$max10, type="realline")
t2 <- Sys.time()
t2-t1
save(max10, fit, file="Wisla.Rdata")

pdf("gestosc chyba.pdf")
hist(max10$max10, prob=TRUE)
dev.off()
fit$family
fit$fits    
fit$parameters  
mu <- fit$mu
sigma <- fit$sigma
nu <- fit$nu
tau <- fit$tau
par(mfrow=c(2,2))

pdf("gestosc.pdf")
hist(max10$max10, prob=TRUE,xlab=NA)
curve(dSEP1(x,mu,sigma,nu,tau),add=T,col=2)
dev.off()

pdf("kwantylkwantyl.pdf")
alpha=ppoints(100)
dev.off()

kwantyle_teo <- qSEP1(alpha,mu,sigma,nu,tau)
kwantyle_emp <- quantile(max10$max10,alpha,na.rm=TRUE)

pdf("kwantyle emp-teo.pdf")
plot(kwantyle_emp,kwantyle_teo)
abline(a=0,b=1,col=2)
dev.off()

X <- as.numeric(na.omit(max10$max10))

fSEP1 <- fitdist(X, "SEP1", start =list(mu=mu,sigma=sigma,nu=nu,tau=tau))
pdf("fSEP1.pdf")
plot(fSEP1)
dev.off()

pdf("dystrybuanta emp-teo.pdf")
plot(ecdf(max10$max10))
curve(pSEP1(x,mu,sigma,nu,tau), xlim=c(-10,35),col=2,add=TRUE)
dev.off()

x20 <- 1-(1/(20*92*24*6))
result20 <- qSEP1(x20,mu,sigma,nu,tau); result20
x50 <- 1-(1/(50*92*24*6))
result50 <- qSEP1(x50,mu,sigma,nu,tau); result50

data <- max10  

data <- max10[,5]
hist(data)
data <- data[!is.na(data)]

library(evir)

b <- 6*24*92
fit1 <- evir::gev(data,b)  
fit1

parGEV <- fit1$par.ests; parGEV

xi <- parGEV[[1]]
sigma <- parGEV[[2]]
mu <- parGEV[[3]]

Max <- fit1$data


pdf("wykres diagnostyczny 1.pdf")
hist(Max,prob=TRUE)
curve(evir::dgev(x,xi,mu,sigma),col=2,add=TRUE) 
dev.off()

pdf("wykres diagnostyczny 2.pdf")
plot(ecdf(Max))
curve(evir::pgev(x,xi,mu,sigma),col=2,add=TRUE) 
dev.off()

kwantyle.emp <- quantile(Max,ppoints(100))
kwantyle.teo <- evir::qgev(ppoints(100),xi,mu,sigma)

pdf("kwantyle gev.pdf")
plot(kwantyle.emp,kwantyle.teo)
abline(a=0,b=1,col=2)
dev.off()

fit2 <- ismev::gev.fit(Max)
fit2$mle

pdf("wykresy diagnostyczne fit2.pdf")
ismev::gev.diag(fit2)
dev.off()

fit3 <- fExtremes::gevFit(data,b)
fit3

pdf("wykresy diagnostyczne fit3.pdf")
summary(fit3)
dev.off()


X <- max10[,5]
u <- 27
X <- na.omit(X)
pdf("wykresy rozrzutu.pdf")
par(mfrow=c(2,1))
plot(X,type="h")
abline(h=u,lwd=2,col='red')   
dev.off()
Y=X[X>u]-u
plot(Y,type='h')

library(ismev)
library(evir)

fitGPD=ismev::gpd.fit(X,u)

xi=fitGPD$mle[[2]]
beta=fitGPD$mle[[1]]
xi; beta
pdf("diag.pdf")
ismev::gpd.diag(fitGPD)
dev.off()
fitGPD=gpd(X,u)

x20 <- 1-(1/(20*92*24*6))
evir::riskmeasures(fitGPD, x20)[2]

x50 <- 1-(1/(50*92*24*6))
evir::riskmeasures(fitGPD, x50)[2]
