library(tidyr)
library(gamlss)
library(maps)
library(fitdistrplus)
library(ismev)
library(evir)

coords <- read.csv("coord233_alt.csv")
city <- coords[coords$place=="WISLA",]
station_id = city['station']$station

pdf("mapa.pdf")
poland <- map('world', 'poland', fill=T, col='gray');
points(city[c('lon', 'lat')], pch=19, col=2);
dev.off()

path_to_files <- "temp.stations-all"
lato <- as.list(list.files(path=path_to_files, pattern = ".*[6-8].csv"))
lato <- paste0(path_to_files, "//", lato)
data_lato <- lapply(lato, read.csv)

jesien <- as.list(list.files(path=path_to_files, pattern = ".*[9-11].csv"))
jesien <- paste0(path_to_files, "//", jesien)
data_jesien <- lapply(jesien, read.csv)

zima <- as.list(list.files(path=path_to_files, pattern = ".*[12-2].csv"))
zima <- paste0(path_to_files, "//", zima)
data_zima <- lapply(zima, read.csv)

wiosna <- as.list(list.files(path=path_to_files, pattern = ".*[3-5].csv"))
wiosna <- paste0(path_to_files, "//", wiosna)
data_wiosna <- lapply(wiosna, read.csv)

### Lato
n_lato <- length(data_lato)
x_lato <- data_lato[[1]]$X249180230 
max10_lato <- c()
datetime_lato <- c()

for(i in 1:n_lato) {
  max10_lato <- c(max10_lato, data_lato[[i]]$X249180230)
  datetime_lato <- c(datetime_lato, as.character(data_lato[[i]]$datetime_lato))
}

x_lato <- c()
for(i_lato in 1:n_lato) {
  sti_lato <- colnames(data_lato[[i]])
  x_lato[i] <- "X249180230" %in% sti_lato
}


max10_lato <- data.frame(date = as.Date(datetime_lato), max10_lato=max10_lato)

max10_lato <- separate(max10_lato, date, c("year", "mth", "day"), convert=TRUE)

max10_lato <- data.frame(datetime_lato=datetime_lato, max10_lato)

fit_lato <- fitDist(max10_lato$max10_lato, type="realline")

save(max10_lato, fit_lato, file="Wisla_lato.Rdata")


mu_lato <- fit_lato$mu
sigma_lato <- fit_lato$sigma
nu_lato <- fit_lato$nu
tau_lato <- fit_lato$tau

print("X20 lato GAMLSS")
x20_lato <- 1-(1/(20*92*24*6))
result20_lato <- qSEP1(x20,mu_lato,sigma_lato,nu_lato,tau_lato); result20_lato
print("X50 lato GAMLSS")
x50_lato <- 1-(1/(50*92*24*6))
result50_lato <- qSEP1(x50,mu_lato,sigma_lato,nu_lato,tau_lato); result50_lato

### 2 sposob
#-=-=-=-=-=-=-=-=- metoda 2. Rozkład GEV -=-=-=-=-=-=-=-=-=-=-=

b_lato <- 6*24*92
fit_lato1 <- evir::gev(data_lato,b_lato)  

Max_lato <- fit_lato1$data_lato

fit_lato2 <- ismev::gev.fit_lato(Max_lato)

print("X20 lato GEV")
x20_lato = evir::rlevel.gev(fit_lato2, k.blocks = 20)
print("X50 lato GEV")
x50_lato = evir::rlevel.gev(fit_lato2, k.blocks = 50)

### 3 sposob
#-=-=-=-=-=-=-=-=- metoda 3. Rozkład GPD. Metoda POT -=-=-=-=-=-=-=-=-=-=-=

X_lato <- max10_lato[,5]
u_lato <- 27
X_lato <- na.omit(X_lato)

fit_latoGPD=ismev::gpd.fit(X_lato,u_lato)   #u=kwantyl 90%, 

print("X20 lato GPD")
x20_lato <- 1-(1/(20*92*24*6))
evir::riskmeasures(fit_latoGPD, x20_lato)[2]
print("X50 lato GPD")
x50_lato <- 1-(1/(50*92*24*6))
evir::riskmeasures(fit_latoGPD, x50_lato)[2]

### jesien
n_jesien <- length(data_jesien)
x_jesien <- data_jesien[[1]]$X249180230 
max10_jesien <- c()
datetime_jesien <- c()

for(i in i_jesien:n_jesien) {
  max10_jesien <- c(max10_jesien, data_jesien[[i]]$X249180230)
  datetime_jesien <- c(datetime_jesien, as.character(data_jesien[[i]]$datetime_jesien))
}

x_jesien <- c()
for(i in 1:n_jesien) {
  sti_jesien <- colnames(x_jesien[[i]])
  x_jesien[i_jesien] <- "X249180230" %in% sti_jesien
}


max10_jesien <- data.frame(date = as.Date(datetime_jesien), max10_jesien=max10_jesien)

max10_jesien <- separate(max10_jesien, date, c("year", "mth", "day"), convert=TRUE)

max10_jesien <- data.frame(datetime_jesien=datetime_jesien, max10_jesien)

fit_jesien <- fitDist(max10_jesien$max10_jesien, type="realline")

save(max10_jesien, fit_jesien, file="Wisla_jesien.Rdata")


mu_jesien <- fit_jesien$mu
sigma_jesien <- fit_jesien$sigma
nu_jesien <- fit_jesien$nu
tau_jesien <- fit_jesien$tau

print("X20 jesien GAMLSS")
x20_jesien <- 1-(1/(20*92*24*6))
result20_jesien <- qSEP1(x20_jesien,mu_jesien,sigma_jesien,nu_jesien,tau_jesien); result20_jesien
print("X50 jesien GAMLSS")
x50_jesien <- 1-(1/(50*92*24*6))
result50_jesien <- qSEP1(x50_jesien,mu_jesien,sigma_jesien,nu_jesien,tau_jesien); result50_jesien

### 2 sposob
#-=-=-=-=-=-=-=-=- metoda 2. Rozkład GEV -=-=-=-=-=-=-=-=-=-=-=

b_jesien <- 6*24*92
fit_jesien1 <- evir::gev(data_jesien,b_jesien)  
fit_jesien1

Max_jesien <- fit_jesien1$data_jesien

fit_jesien2 <- ismev::gev.fit_jesien(Max_jesien)

print("X20 jesien GEV")
x20_jesien = evir::rlevel.gev(fit_jesien2, k.blocks = 20)
print("X50 jesien GEV")
x50_jesien = evir::rlevel.gev(fit_jesien2, k.blocks = 50)

### 3 sposob
#-=-=-=-=-=-=-=-=- metoda 3. Rozkład GPD. Metoda POT -=-=-=-=-=-=-=-=-=-=-=

X_jesien <- max10_jesien[,5]
u_jesien <- 27
X_jesien <- na.omit(X_jesien)

fit_jesienGPD=ismev::gpd.fit(X_jesien,u_jesien)   #u=kwantyl 90%, 

print("X20 jesien GPD")
x20_jesien <- 1-(1/(20*92*24*6))
evir::riskmeasures(fit_jesienGPD, x20_jesien)[2]
print("X50 jesien GPD")
x50_jesien <- 1-(1/(50*92*24*6))
evir::riskmeasures(fit_jesienGPD, x50_jesien)[2]

### zima
n_zima <- length(data_zima)
x_zima <- data_zima[[1]]$X249180230 
max10_zima <- c()
datetime_zima <- c()

for(i in i_zima:n_zima) {
  max10_zima <- c(max10_zima, data_zima[[i]]$X249180230)
  datetime_zima <- c(datetime_zima, as.character(data_zima[[i]]$datetime_zima))
}

x_zima <- c()
for(i in 1:n_zima) {
  sti_zima <- colnames(data__zima[[i]])
  x_zima[i] <- "X249180230" %in% sti_zima
}


max10_zima <- data.frame(date = as.Date(datetime_zima), max10_zima=max10_zima)

max10_zima <- separate(max10_zima, date, c("year", "mth", "day"), convert=TRUE)

max10_zima <- data.frame(datetime_zima=datetime_zima, max10_zima)

fit_zima <- fitDist(max10_zima$max10_zima, type="realline")

save(max10_zima, fit_zima, file="Wisla_zima.Rdata")


mu_zima <- fit_zima$mu
sigma_zima <- fit_zima$sigma
nu_zima <- fit_zima$nu
tau_zima <- fit_zima$tau

print("X20 zima GAMLSS")
x20_zima <- 1-(1/(20*92*24*6))
result20_zima <- qSEP1(x20,mu_zima,sigma_zima,nu_zima,tau_zima); result20_zima
print("X50 zima GAMLSS")
x50_zima <- 1-(1/(50*92*24*6))
result50_zima <- qSEP1(x50,mu_zima,sigma_zima,nu_zima,tau_zima); result50_zima

### 2 sposob
#-=-=-=-=-=-=-=-=- metoda 2. Rozkład GEV -=-=-=-=-=-=-=-=-=-=-=

b_zima <- 6*24*92
fit_zima1 <- evir::gev(data_zima,b_zima)  

Max_zima <- fit_zima1$data_zima

fit_zima2 <- ismev::gev.fit_zima(Max_zima)

print("X20 zima GEV")
x20_zima = evir::rlevel.gev(fit_zima2, k.blocks = 20)
print("X50 zima GEV")
x50_zima = evir::rlevel.gev(fit_zima2, k.blocks = 50)

### 3 sposob
#-=-=-=-=-=-=-=-=- metoda 3. Rozkład GPD. Metoda POT -=-=-=-=-=-=-=-=-=-=-=

X_zima <- max10_zima[,5]
u_zima <- 27
X_zima <- na.omit(X_zima)

fit_zimaGPD=ismev::gpd.fit(X_zima,u_zima)   #u=kwantyl 90%, 

print("X20 zima GPD")
x20_zima <- 1-(1/(20*92*24*6))
evir::riskmeasures(fit_zimaGPD, x20_zima)[2]
print("X50 zima GPD")
x50_zima <- 1-(1/(50*92*24*6))
evir::riskmeasures(fit_zimaGPD, x50_zima)[2]

### wiosna
n_wiosna <- length(data_wiosna)
x_wiosna <- data_wiosna[[1]]$X249180230 
max10_wiosna <- c()
datetime_wiosna <- c()

for(i in i_wiosna:n_wiosna) {
  max10_wiosna <- c(max10_wiosna, data_wiosna[[i]]$X249180230)
  datetime_wiosna <- c(datetime_wiosna, as.character(data_wiosna[[i]]$datetime_wiosna))
}

x_wiosna <- c()
for(i in 1:n_wiosna) {
  sti_wiosna <- colnames(data_wiosna[[i]])
  x_wiosna[i] <- "X249180230" %in% sti_wiosna
}


max10_wiosna <- data.frame(date = as.Date(datetime_wiosna), max10_wiosna=max10_wiosna)

max10_wiosna <- separate(max10_wiosna, date, c("year", "mth", "day"), convert=TRUE)

max10_wiosna <- data.frame(datetime_wiosna=datetime_wiosna, max10_wiosna)

fit_wiosna <- fitDist(max10_wiosna$max10_wiosna, type="realline")

save(max10_wiosna, fit_wiosna, file="Wisla_wiosna.Rdata")


mu_wiosna <- fit_wiosna$mu
sigma_wiosna <- fit_wiosna$sigma
nu_wiosna <- fit_wiosna$nu
tau_wiosna <- fit_wiosna$tau

print("X20 wiosna GAMLSS")
x20_wiosna <- 1-(1/(20*92*24*6))
result20_wiosna <- qSEP1(x20,mu_wiosna,sigma_wiosna,nu_wiosna,tau_wiosna); result20_wiosna
print("X50 wiosna GAMLSS")
x50_wiosna <- 1-(1/(50*92*24*6))
result50_wiosna <- qSEP1(x50,mu_wiosna,sigma_wiosna,nu_wiosna,tau_wiosna); result50_wiosna

### 2 sposob
#-=-=-=-=-=-=-=-=- metoda 2. Rozkład GEV -=-=-=-=-=-=-=-=-=-=-=

b_wiosna <- 6*24*92
fit_wiosna1 <- evir::gev(data_wiosna,b_wiosna)  
fit_wiosna1

Max_wiosna <- fit_wiosna1$data_wiosna

fit_wiosna2 <- ismev::gev.fit_wiosna(Max_wiosna)

print("X20 wiosna GEV")
x20_wiosna = evir::rlevel.gev(fit_wiosna2, k.blocks = 20)
print("X50 wiosna GEV")
x50_wiosna = evir::rlevel.gev(fit_wiosna2, k.blocks = 50)

### 3 sposob
#-=-=-=-=-=-=-=-=- metoda 3. Rozkład GPD. Metoda POT -=-=-=-=-=-=-=-=-=-=-=

X_wiosna <- max10_wiosna[,5]
u_wiosna <- 27
X_wiosna <- na.omit(X_wiosna)

fit_wiosnaGPD=ismev::gpd.fit(X_wiosna,u_wiosna)   #u=kwantyl 90%, 

print("X20 wiosna GPD")
x20_wiosna <- 1-(1/(20*92*24*6))
evir::riskmeasures(fit_wiosnaGPD, x20_wiosna)[2]
print("X50 wiosna GPD")
x50_wiosna <- 1-(1/(50*92*24*6))
evir::riskmeasures(fit_wiosnaGPD, x50_wiosna)[2]
