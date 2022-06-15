library(tidyr)

coords = read.csv("D:/IT/Nauka/R/DATA/coord233_alt.csv")
path_to_files <- "D:/IT/Nauka/R/DATA/wind_stations_all"

L <- as.list(list.files(path=path_to_files, pattern = ".*8.csv"))
L <- paste0(path_to_files,"\\",L)
dataWind <- lapply(L,read.csv)
path_to_files <- "D:/IT/Nauka/R/DATA/temp.stations-all"

L <- as.list(list.files(path=path_to_files, pattern = ".*8.csv"))
L <- paste0(path_to_files,"\\",L)
data0 <- lapply(L,read.csv)
n <- length(data0)
newData <- list()
missing <- c()
min_temp <- 0
max_temp <- 40
data_max <- 6 * 24 * 31 * 11
for(i in 1:nrow(coords)){
  station_temps <- c()
  station_winds <- c()
  for(x in 1:n){
    if(exists(coords[i,]$station, data0[[x]])){
      station_temps <- c(station_temps, get(coords[i,]$station, data0[[x]]))
    }
    if(exists(coords[i,]$station, dataWind[[x]])){
      station_winds <- c(station_winds, get(coords[i,]$station, dataWind[[x]]))
    }
  }
  if(length(station_temps) == data_max && length(station_winds) == data_max){
    station_temps[station_temps > max_temp] = NA
    station_temps[station_temps < min_temp] = NA
    station_winds[station_winds < 0] = NA
    newData[[length(newData)+1]] <- list(station_temps, station_winds, coords[i,])
    missing <- c(missing, sum(is.na(station_temps)))
  }
}
median_missing <- median(missing)

for(i in 1:length(missing)){
  if(missing[i] > median_missing){
    newData <- newData[-i]
  }
}
datetime <- c()

for(i in 1:n){
  datetime <- c(datetime,as.character(data0[[i]]$datetime))
}
library("dplyr")
formatted <- list()
for(i in 1: length(newData)){
  formatted[[length(formatted) + 1]] <- list(
    newData[[i]][[3]],
    data.frame(date=as.Date(datetime),maxd_temp=newData[[i]][[1]],maxd_wind=newData[[i]][[2]]))
  rownames(formatted[[i]][[2]]) <- c()
  formatted[[i]][[2]] <- separate(
    formatted[[i]][[2]],
    date,
    c("year","mth","day"),
    convert=TRUE
  )
  formatted[[i]][[2]] <- data.frame(datetime=datetime,formatted[[i]][[2]])
  formatted[[i]][[2]] <- formatted[[i]][[2]][complete.cases(formatted[[i]][[2]]),]
  tmp_max_wind <- formatted[[i]][[2]] %>% group_by(year, mth, day) %>% slice(which.max(maxd_wind))
  formatted[[i]][[2]] <- formatted[[i]][[2]] %>% group_by(year, mth, day) %>% slice(which.max(maxd_temp))
  formatted[[i]][[2]]$maxd_wind <- tmp_max_wind$maxd_wind
  print(round(100 * i / length(newData), 1))
}
save(formatted, file="D:/IT/Nauka/R/DATA/Daily_Max_August.Rdata")
library(copula)
library(VineCopula)
library(ggplot2)
library(ggExtra)
library(MASS)

cop_stats <-c()

for(station in formatted){
  X1 <- station[[2]]$maxd_temp
  X2 <- station[[2]]$maxd_wind
  X <- data.frame(X1,X2)
  X <- X[complete.cases(X),]
  V <- pobs(X)
  colnames(V) <- c("V1","V2")
  cop.npar <- BiCopSelect(V[,1],V[,2], selectioncrit="AIC", se=TRUE)
  norm.cop <- BiCop(family = cop.npar$family, par = cop.npar$par, par2= cop.npar$par2)
  cop_stats <- c(cop_stats, norm.cop$familyname)
  N <- 30 * 11
  temp <- c()
  wind <- c()
  for (i in 1 : 100) {
    V <- BiCopSim(N,norm.cop)
    Z.npar <- cbind(quantile(X$X1, V[,1], na.rm=TRUE), quantile(X$X2, V[,2]))
    max_blocks <- c()
    for(k in 1 : 11){
      block_temp <- c();block_wind <- c();
      from <- 1 + (k-1) * 30
      to <- k * 30
      block_temp <- Z.npar[ from : to, 1]
      block_wind <- Z.npar[ from : to, 2]
      temp <- c(temp, max(block_temp))
      wind <- c(wind, max(block_wind))
      
    }
    
  }
  
  
  
  
}
occurence <- data.frame(table(unlist(cop_stats)))
colnames(occurence) <- c("Rodzaj kopuly", "liczba wystapien")
#
#
#to dalej nie dotyczy projektu 
#
#













X <- data.frame(X1,X2)
X <- X[complete.cases(X),]
#estymujemy parametry
fit.exp=fitdistr(X$X1,"exponential")
fit.norm=fitdistr(X$X2,"normal")

#parametry wyestymowanych rozkladow i qq-ploty (wykresy diagnostyczne na dobroc dopasowania rozkladu)
lambda=round(fit.exp$estimate,2)
par=round(fit.norm$estimate,2)
lambda; par

alpha <- ppoints(100)
X1_emp <- quantile(X1,alpha)
X2_emp <- quantile(X2,alpha)
X1_teo <- qexp(alpha,lambda)
X2_teo <- qnorm(alpha,par[1],par[2])


#--- B. podejscie nieparametryczne 
#pseudo-obserwacje nieparametryczne(korzystamy z funkcji pobs w library(copula))
V <- pobs(X)
colnames(V) <- c("V1","V2")
#===== w bibliotece VineCopula za pomoca funkcji 'BiCopSelect'
#===== wybieramy najlepsza, z zaimplementowanych tam kopul (kryterium wyboru AIC)
#library(VineCopula)
cop.npar <- BiCopSelect(V[,1],V[,2], selectioncrit="AIC", se=TRUE)     #na podstawie obserwavji nieparametrycznych

#===== wyniki na wykresach (gestosci, kontury, dane vs. wygenerowane) 
#parametry kazdej z kopul potrzebne dalej do wykresow

theta.npar <-  cop.npar$par

#tworzymy obiekty 'copula' dla kazdej z kopul
norm.cop <- BiCop(family = cop.npar$family, par = cop.npar$par, par2= cop.npar$par2)
norm.cop

#--- wykres konturowy (z rozkladami brzegowymi N(0,1))
#Zobacz CC, Def.3.11, p.59

contour(norm.cop)

#Wygenerujmy, dla porownania, proby licznosci danych: N=1000, 
#z obydwu kopul i wykorzystamy je do utworzenia prob z rozkladu F=C(F1,F2)
N <- 1000

#--- proby z kopul
V <- BiCopSim(N,norm.cop)
V2 <-BiCopSim(N, cop)
u1=U[,1]; u2=U[,2]
v1=V[,1]; v2=V[,2]

#--- proby z rokladow F=C(F1,F2)
X.npar <- cbind(quantile(X1,v1),quantile(X2,v2))            #podejscie nieparametryczne
Z.npar <- cbind(quantile(X$X1, v1, na.rm=TRUE), quantile(X$X2, v2, na.rm=TRUE))
par(mfrow=c(2,3))
plot(X); plot(X.par); plot(X.npar)
plot(U); plot(V)
plot(Z.npar)
#mozemy  wykorzystac funkcje 'BiCopCompare' do selekcji kopul
#wedlug wartosci AIC, BIC czy LogLik
#Mozemy obejrzec wygenerowane probki z kopul przy roznych rozkladach brzegowych
#------------------------------------------------------------------------------
BiCopCompare(U[,1], U[,2])
BiCopCompare(V[,1], V[,2])

#===========
#Przyklad 3. Generowanie z kopul: Claytona, Franka, Gumbela (w bibliotece 'copula')
#===========
#proby z kopul Claytona 
#---------------
par(mfrow=c(3,3))

theta <- c(seq(-1,0,by=0.2),seq(1,5,by=2))
theta

for(i in theta){
  uC <- rCopula(copula=claytonCopula(i), n=1000)
  plot(uC)
  tau=round(i/(i+2),2)
  text(0.3,0, i, col=2)
  #text(0.4,0, tau, col=2)
}

#proby z kopul Gumbela
#--------------------
theta <- seq(1,5,by=0.5)

par(mfrow=c(3,3))

for(i in theta){
  uC <- rCopula(copula=gumbelCopula(i), n=1000)
  plot(uC)
  text(0.8,0, i, col=2)
}

#proby z kopul Franka
#--------------------
theta <- seq(-16,16,by=4)

par(mfrow=c(3,3))

for(i in theta){
  uC <- rCopula(copula=frankCopula(i), n=1000)
  plot(uC)
  text(0.8,0, i, col=2)
}


#proby z kopul Joe
#--------------------
theta <- seq(1,9,by=1)

par(mfrow=c(3,3))

for(t in theta){
  
  uC <- rCopula(1000,joeCopula(t))
  plot(uC)
  text(0.8,0, t, col=2)
}


#Przyklad 4 -- wspolczynniki:Pearsona, Kendalla i Spearmana
#==========
?cor
#--- A wspolczynniki dla zestawu ,,dane''
#Z <- pobierz plik dane.csv z PE (dane do Przykladu 5 w skrypcie kopuly01.R)
#Z <- X z Przykladu 3
Z <- X
sapply(list("pearson", "kendall", "spearman"), function(x) cor(X,method=x)[1,2])

#--- B wsp. Pearsona zalezy od rozkladow brzegowych
#Wykorzystamy dane z Przykladu 1 - kopula gaussowska zlozona z roznymi rozkladami brzegowymi
norm.cop <-  normalCopula(0.7)

N= 1000

set.seed=1010
U <- rCopula(N,norm.cop)                           #rozklady brzegowe U(01)
X <- cbind(qexp(U[,1]),qexp(U[,2]))                #exp(1)
Y <- cbind(qnorm(U[,1],1,2),qexp(U[,2],1/4))       #N(0,1), Exp(1/4)

df <- data.frame(U=U, X=X,Y=Y)
head(df)

#dane na wykresach
p1 <- ggplot(df, aes(x=U.1, y=U.2))+  geom_point()
p2 <- ggplot(df, aes(x=X.1, y=X.2))+  geom_point()
p3 <- ggplot(df, aes(x=Y.1, y=Y.2))+  geom_point()
p1.hist <- ggMarginal(p1, type="histogram")
p2.hist <- ggMarginal(p2, type="histogram")
p3.hist <- ggMarginal(p3, type="histogram")

cowplot::plot_grid(p1.hist,p2.hist, p3.hist,ncol = 1,nrow = 3)

df <- list(U=U, X=X,Y=Y)
sapply(df,function(x) cor(x, method="pearson")[1,2])
sapply(df,function(x) cor(x, method="kendall")[1,2])
sapply(df,function(x) cor(x, method="spearman")[1,2])


#Przyklad 5 (wiatry Ustka-Leba)
#==========
#LU <- pobierz dane z pliku ,,LebaUstka_2015_2018_winter.Rdata"
#LU <- LU <- load(file="sciezka dostepu/LebaUstka_2015_2018_winter.Rdata")
LU

head(Leba); head(Ustka)
tail(Leba); tail(Ustka)
dim(Leba); dim(Ustka)

wind <- data.frame(w_Leba=Leba$w_sp,w_Ustka=Ustka$w_sp)
head(wind)
dim(wind)

X <- wind[complete.cases(wind),]  #wybieramy wiersze z pelnymi danymi
dim(X)

#===== dobor kopuly metoda  nieparametryczna 
#pseudo-obserwacje nieparametryczne
U <- pobs(X)
colnames(U) <- c("u","v")
Y <- qnorm(U)  #'normalna' normalizacja
colnames(Y) <- c("y1","y2")

#wykresy rozrzutu
df <- data.frame(X,U,Y)
head(df)

p1 <- ggplot(df, aes(w_Leba,w_Ustka))+  geom_point()
p2 <- ggplot(df, aes(u,v))+  geom_point()
p3 <- ggplot(df, aes(y1,y2))+  geom_point()
p1.hist <- ggMarginal(p1, type="histogram")
p2.hist <- ggMarginal(p2, type="histogram")
p3.hist <- ggMarginal(p3, type="histogram")

cowplot::plot_grid(p1.hist,p2.hist,p3.hist,ncol = 1,nrow = 3)

#--- dobor kopuly
t1 <- Sys.time()
cop.npar <- BiCopSelect(U[,1],U[,2])
t2 <- Sys.time()
t2-t1 #ok. 30s

cop.npar
cop.npar$family

#--- inne kopuly - porownanie, sortujemy wzgledem AIC, BIC
t1 <- Sys.time()
comp.npar <- BiCopEstList(U[,1],U[,2])
t2 <- Sys.time() 
t2-t1 #30s

comp.npar

AIC3.npar <- head(comp.npar$summary[order(comp.npar$summary$AIC),],3) #204,2,1
BIC3.npar <- head(comp.npar$summary[order(comp.npar$summary$BIC),],3) #204,1,2
logLik3.npar <- head(comp.npar$summary[order(comp.npar$summary$logLik,decreasing = TRUE),],3)#204,2,1
AIC3.npar; BIC3.npar; logLik3.npar

#===== dobor kopuly metoda parametryczna 
#tworzymy pseudo-obserwacje parametryczne

#korzystajac z funkcji fitDist() z biblioteki 'gamlss' dobieramy 
#najlepszy z rozkladow zaimplementowanych w gamlss.family
library(gamlss)

t1 <- Sys.time()
Fits <- lapply(1:2,function(i) fitDist(X[,i]))
t2 <- Sys.time()
t2-t1 #ok. 2 min

#gdyby z jakichs powodow estymacja wyzej nie poszla
#wyniki sa na PE, w pliku LebaUstka_gamlss.Rdata
#plik zawiera liste Fits z wynikami
#fit <- load(file="sciezka dostepu/LebaUstka_gamlss.Rdata")
#fit

#rozklady posortowane wedlug wartosci AIC
Fits[[1]]$fits  #Leba
Fits[[2]]$fits   #Ustka

#Przeksztalcamy AIC (exp((AICmin???AICi)/2)), aby zobaczyc na ile istotnie kolejne modele róznia sie od pierwszego
#Zobacz: How to use AIC in practice, https://en.wikipedia.org/wiki/Akaike_information_criterion
par(mfrow=c(2,1))
for(i in 1:2){
  plot(exp((Fits[[i]]$fits[[1]]- Fits[[i]]$fits)/2))
}

#dopasowane rozklady
lapply(1:2,function(i) Fits[[i]]$fits)
sapply(1:2,function(i) Fits[[i]]$family)
sapply(1:2,function(i) Fits[[i]]$family[[1]])


#histogramy i QQ-ploty
par(mfrow=c(2,2))
hist(X$w_Leba,prob=TRUE)
curve(dSEP4(x, mu=Fits[[1]]$mu,
            sigma= Fits[[1]]$sigma,
            nu = Fits[[1]]$nu,
            tau = Fits[[1]]$tau),col=2,add=TRUE)
hist(X$w_Ustka,prob=TRUE)
curve(dRG(x,mu=Fits[[2]]$mu,sigma= Fits[[2]]$sigma),col=2,add=TRUE)

alpha <- ppoints(100)
X1emp <- quantile(X$w_Leba,alpha)
X1teo <- qSEP4(alpha, mu=Fits[[1]]$mu,
               sigma= Fits[[1]]$sigma,
               nu = Fits[[1]]$nu,
               tau = Fits[[1]]$tau)
X2emp <- quantile(X$w_Ustka,alpha)
X2teo <- qRG(alpha, mu=Fits[[2]]$mu,
             sigma= Fits[[2]]$sigma)


plot(X1emp,X1teo)
abline(a=0,b=1,col=2)
plot(X2emp,X2teo)
abline(a=0,b=1,col=2)


#==== pseudo-obserwacje parametryczne
V <- cbind(pSEP4(X[,1], mu=Fits[[1]]$mu,
                 sigma= Fits[[1]]$sigma,
                 nu = Fits[[1]]$nu,
                 tau = Fits[[1]]$tau),
           pRG(X[,2], mu=Fits[[2]]$mu,
               sigma= Fits[[2]]$sigma))

colnames(V) <- c("v1","v2")
Y <- qnorm(V)
colnames(Y) <- c("y1","y2")

#wykresy rozrzutu
df <- data.frame(X,V,Y)
head(df)

#p1 <- ggplot(df, aes(w_Leba,w_Ustka))+  geom_point()
p4 <- ggplot(df, aes(v1,v2))+  geom_point()
p5 <- ggplot(df, aes(y1,y2))+  geom_point()
p4.hist <- ggMarginal(p4, type="histogram")
p5.hist <- ggMarginal(p5, type="histogram")

cowplot::plot_grid(p1.hist,p4.hist,p5.hist,ncol = 1,nrow = 3)

#--- dobor kopuly
t1 <- Sys.time()
cop.par <- BiCopSelect(V[,1],V[,2])
t2 <- Sys.time()
t2-t1 #ok. 40s

cop.par
cop.par$family

#--- inne kopuly - porownanie, sortujemy wzgledem AIC, BIC
t1 <- Sys.time()
comp.par <- BiCopEstList(V[,1],V[,2])
t2 <- Sys.time() 
t2-t1 #30s

comp.par

#Pierwsze trzy ,,najlepsze kopuly'' 
AIC3.par <- head(comp.par$summary[order(comp.par$summary$AIC),],3) #204,2,1
BIC3.par <- head(comp.par$summary[order(comp.par$summary$BIC),],3) #204,1,2
logLik3.par <- head(comp.par$summary[order(comp.par$summary$logLik,decreasing = TRUE),],3) #204,2,1
AIC3.par; BIC3.par; logLik3.par

#----
cowplot::plot_grid(p1.hist,p2.hist,p4.hist,ncol = 3,nrow = 1)

#Przyklad 6 (Diagnostyka)
#=========
#1. porownanie gestosci empirycznej i teoretycznej kopuly
p6 <- BiCopKDE(U[,1],U[,2],type = "surface")
p7 <- plot(cop.npar)
p8 <- plot(cop.par)

cowplot::plot_grid(p6,p7,p8,ncol = 1,nrow = 3)

#2. --- porownanie konturu 'empirycznego' i teoretycznego (kopuly)
par(mfrow=c(4,1))
BiCopKDE(U[,1],U[,2],type = "contour")
BiCopKDE(V[,1],V[,2],type = "contour")
contour(cop.npar)
contour(cop.par)

#kontur 'empiryczny' inaczej
UU <- as.copuladata(U)
VV <- as.copuladata(V)
pairs(UU)
pairs(VV)

#3. --- proby wygenerowane z kopul
N <- dim(X)[1]; N

sem.npar <- BiCopSim(N,cop.npar)
sem.par <-BiCopSim(N,cop.par)

sem <- data.frame(sem.npar=sem.npar,sem.par=sem.par)

p9 <- ggplot(sem, aes(sem.npar.1,sem.npar.2))+  geom_point()
p10 <- ggplot(sem, aes(sem.par.1,sem.par.2))+  geom_point()
p9.hist <- ggMarginal(p9, type="histogram")
p10.hist <- ggMarginal(p10, type="histogram")

cowplot::plot_grid(p9.hist,p10.hist,ncol = 1,nrow = 2)

#4. -- proby z rozkladow F=C(F1,F2)
u1 <- sem.npar[,1];  u2 <- sem.npar[,2]
v1 <- sem.par[,1];  v2 <- sem.par[,2]

x1.npar <- quantile(X$w_Leba,u1)
x2.npar <- quantile(X$w_Ustka,u2)
x1.par <- qSEP4(v1, mu=Fits[[1]]$mu,
                sigma= Fits[[1]]$sigma,
                nu = Fits[[1]]$nu,
                tau = Fits[[1]]$tau)
x2.par <- qRG(v2,mu=Fits[[2]]$mu,
              sigma= Fits[[2]]$sigma)

sem.F <- data.frame(x1.npar=x1.npar,x2.npar =x2.npar,
                    x1.par=x1.par,x2.par =x2.par)
head(sem.F)

p11 <- ggplot(sem.F, aes(x1.npar,x2.npar))+  geom_point()
p12 <- ggplot(sem.F, aes(x1.par,x2.par))+  geom_point()
p11.hist <- ggMarginal(p11, type="histogram")
p12.hist <- ggMarginal(p12, type="histogram")

cowplot::plot_grid(p11.hist,p12.hist,ncol = 1,nrow = 2)

#5. --- porownanie wspolczynnikow ekstremalnych empirycznych  i teoretycznych 
#estymacja dolnego i gornego wspolczynnika
p <- 0.01 # cut-off

(lam.C <- c(lower = fitLambda(U, p = p)[2,1],
            upper = fitLambda(U, p = p, lower.tail = FALSE)[2,1])) 

(lam.C <- c(lower = fitLambda(V, p = p)[2,1],
            upper = fitLambda(V, p = p, lower.tail = FALSE)[2,1])) 

#wspolczynniki dla wyestymowanej kopuly
BiCopPar2TailDep(cop.npar)
BiCopPar2TailDep(cop.par)

#6.--- wspolczynniki zaleznosci ekstremalnych dla pierwszych
#trzech kopul wybranych przez AIC 
nAIC3.npar <- as.numeric(rownames(AIC3.npar))
lu.coeff.AIC3 <- t(sapply(nAIC3.npar,function(i) BiCopPar2TailDep(comp.npar$models[[i]])))
rownames(lu.coeff.AIC3) <- BiCopName(AIC3.npar[,1])
lu.coeff.AIC3
