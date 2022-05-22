library(tidyr)
library(gamlss)
library(maps)
library(fitdistrplus)

coords <- read.csv("coord233_alt.csv")
city <- coords[coords$place=="WISLA",]
station_id = city['station']$station

reload_file <- FALSE

path_to_files <- "temp.stations-all"
L <- as.list(list.files(path=path_to_files, pattern = ".*((([0][3])|([0][4])|([0][5])).csv)$"))
L <- paste0(path_to_files, "//", L)
data0 <- lapply(L, read.csv)

n <- length(data0); n
i <- 1
x <- data0[[i]]$X249180230 
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
save(max10, fit, file="WislaJesien.Rdata")

mu <- fit$mu
sigma <- fit$sigma
nu <- fit$nu
tau <- fit$tau

X <- as.numeric(na.omit(max10$max10))

print("1 metoda")

x20 <- 1-(1/(20*90*24*6))
result20 <- qSEP1(x20,mu,sigma,nu,tau); result20
x50 <- 1-(1/(50*90*24*6))
result50 <- qSEP1(x50,mu,sigma,nu,tau); result50

data <- max10  #wstaw 10-minutowe maksima dla Twojej stacji i wykonaj dalsze analizy

data <- max10[,5]
data <- data[!is.na(data)]

library(evir)

#dla uproszczenia przyjmijmy ze w lecie jest 6*24*90 obserwacji
b <- 6*24*90
fit1 <- evir::gev(data,b)

parGEV <- fit1$par.ests
xi <- parGEV[[1]]
sigma <- parGEV[[2]]
mu <- parGEV[[3]]
Max <- fit1$data
fit2 <- ismev::gev.fit(Max)
print("2 metoda")

print("X20 lato GEV")
x20 = evir::rlevel.gev(fit1, k.blocks = 20); x20
print("X50 lato GEV")
x50 = evir::rlevel.gev(fit1, k.blocks = 50); x50

### 3 sposob
#-=-=-=-=-=-=-=-=- metoda 3. Rozkład GPD. Metoda POT -=-=-=-=-=-=-=-=-=-=-=

X <- max10[,5]
u <- 27
X <- na.omit(X)

#nadwyzki nad prog u
Y=X[X>u]-u

#estymujemy parametry rozkladu GPD
#gpd(dane,u) - podajemy prog lub liczbe nadwyzek
library(ismev)
library(evir)

fitGPD=ismev::gpd.fit(X,u)   #u=kwantyl 90%, 
#można podać liczbę nadwyżek, tutaj 0.10*length(X)

#wyestymowane parametry rozkladu GPD
xi=fitGPD$mle[[2]]
beta=fitGPD$mle[[1]]

ismev::gpd.diag(fitGPD)

fitGPD <-gpd(X,u)
print("3 metoda")
x20 <- 1-(1/(20*90*24*6))
riskmeasures(fitGPD, x20)[2]

x50 <- 1-(1/(50*90*24*6))
riskmeasures(fitGPD, x50)[2]
