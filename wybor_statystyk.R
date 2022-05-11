#Laboratorium: BMM-automatyczny wybor r gornych statystyk
library(eva)


#1. Z każdego roku 2008-2018 wez R=10 najwiekszych temperatur
#utworz macierz/ramke danych o rozmiarze 11 x 10 - dalej nazwana ,,os10''
#temperatury w wierszach od najwiekszej do najmniejszej
#=====================================================================
md <- load(file="Wisla.Rdata")
data <- max10
head(data)

r <- 10
os10 <- matrix(NA, nrow=length(2008:2018), ncol=r)

count <- 0

for (i in 2008:2018) {
  count <- count + 1
  datar <- data[data$year == i, -c(1:4)]
  x <- sort(datar, decreasing = TRUE)
  y <- x[1:r]
  os10[count, ] <- y
}
os10

#2. Wyestymuj parametry rozkladu GEVr wykorzystujac kolejno:
#jedna, dwie,trzy,...,10 maksymalnych wartosci z roku.
#Oblicz poziomy zwrotu x20 i x50 (z przedzialami ufnosci 95%).
#(uzyj funkcji  gevrFit i gevrRl z  biblioteki 'eva')
#====================================
fit10 <- gevrFit(os10,method="mle") #wykorzystane 10 gornych statytyk
fit10$par.ests

#poziomy zwrotu
gevrRl(fit10,20)
gevrRl(fit10,50)

#parametry, poziomy zwrotu i przedzialy ufnosci
#otrzymane z wykorzystaniem roznej liczby ekstremalnych wartosci
parGEVr <- matrix(NA,nrow=10,ncol=3)
RL <- matrix(NA,nrow=10,ncol=2)
CI <- matrix(NA,nrow=10,ncol=4)

for(r in 1:10){
  
  osr <- os10[,1:r]
  fit <- gevrFit(osr,method="mle")
  
  parGEVr[r,] <- fit$par.ests
  
  x20 <- gevrRl(fit,20)
  x50 <- gevrRl(fit,50)
  
  RL[r,] <- c(x20$Estimate,x50$Estimate)
  CI[r,] <- c(x20$CI,x50$CI)
  
}

colnames(parGEVr) <- c("loc","scale","shape")
colnames(RL) <- c("rl20","rl50")

parGEVr
RL
CI

#Wyestymowane parametry na wykresach
par(mfrow=c(3,1))
plot(parGEVr[,1],type="b",xlab="r",ylab="mu")
plot(parGEVr[,2],type="b",xlab="r",ylab="beta")
plot(parGEVr[,3],type="b",xlab="r",ylab="xi")

#wyestymowane poziomy zwrotu z przedzialami ufnosci
par(mfrow=c(3,1))
plot(RL[,1],type="b",ylim=c(36.5,37.3),xlab="r",ylab="rl20 / rl50")
points(RL[,2],type="b",col=2)
grid()

plot(RL[,1],type="b",ylim=c(33.5,40),xlab="r",ylab="rl20")
points(CI[,1],type="l",col="blue")
points(CI[,2],type="l",col="blue")
grid()

plot(RL[,2],type="b",ylim=c(33,42),xlab="r",ylab="rl50",col=2)
points(CI[,3],type="l",col="blue")
points(CI[,4],type="l",col="blue")
grid()


#3. Narysuj wykresy gestosci rozkladu GEV z wyestymowanymi parametrami
#==========================================
my.col <- topo.colors(10,rev=TRUE)

par(mfrow=c(1,1))

curve(evir::dgev(x,mu=parGEVr[1,1],sigma=parGEVr[1,2],xi=parGEVr[1,3]),
      xlim=c(28,40),ylim=c(0,0.30),col=my.col[1],
      xlab=NA,ylab=NA)

for(i in 1:10){
  
  curve(evir::dgev(x,mu=parGEVr[i,1],sigma=parGEVr[i,2],xi=parGEVr[i,3]),
        col=my.col[i],add=TRUE)
}


#4. Przetestuj hipoteze zerowa, ze dla r=10 otrzymany rozklad jest
#rzeczywiscie rozkladem GEVr z wyestymowanymi parametrami.
#Wykorzystaj test Rao - bootstrap parametryczny, funkcja 'gevrPbScore' w bib. 'eva'. 
#Ustawiona jest mala liczba B probek bootsrapowych, sprawdz jak dlugo trwaja obliczenia.
#===============================================================================
#test Rao
B <- 50 #Uwaga. Jest to za mała liczba probek bootstrapowych do uzyskania wiarygodnego wyniku.

t1<-Sys.time()

testPbS <- gevrPbScore(os10, bootnum = B)

t2<-Sys.time()
print(t2-t1) 

testPbS

#test roznicy entropii (krotki czas obliczen, mniej wiarygodny dla malych probek)
testEd <- gevrEd(os10)
testEd

#5. Wyznacz optymalna liczbe r gornych statystyk pozycyjnych 
#(z maksymalnej R=10), metodami BH i ForwardStop. 
#Wykorzystaj test Rao. 
#Gotowe rozwiazanie do seryjnego testowania 
#to funkcja 'gevrSeqTests' w bibliotece eva 
#===========================================================
B <- 50 #Uwaga. Jest to za mała liczba probek  bootstrapowych do uzyskania wiarygodnego wyniku.

t1<-Sys.time()
test10 <- gevrSeqTests(os10,method="pbscore",bootnum=B)
t2<-Sys.time()
print(t2-t1) #ok. 45s

test10

#Inna mozliwosc, to test entropii (krotki czas obliczen, mniej wiarygodny dla malych probek,
#nie uwzglednia tez modelu z jedna wartoscia ekstremalna)
test10_ed <- gevrSeqTests(os10,method="ed")

#p-value dla r=10,9,...,1 
#zmieniamy kolejnosc od najwiekszej liczby statystyk do najmniejszej
p_value <- rev(test10$p.values)
length(p_value)

#metoda ForwardStop (prog odrzucen alpha=5%)
#------------------
p_FS <- rev(test10$ForwardStop)
p_FS

plot(p_FS,ylim=c(0,1),pch=19)
abline(h=0.05,col=2)
#Jakie jest optymalne r?

#methoda BH - wyliczamy q-value
#-------------------------------
R <- 10
Yi <- -log(1-p_value)
ri <- rev(1:R)
Zi0 <- Yi/ri
Zi <- cumsum(Zi0)
qi <- 1-exp(-Zi) #q-value

#Dwie metody odrzucen kolejnych qi w metodzie BH
#A. prog odrzucen kolejnych qi (klasyczny)
alpha <- 0.05
R <- 10
os <- (1:R)/R
osl <- alpha*os

#B. prog odrzucen kolejnych  qi - jako kwantyl alpha=5% rozkladu 
#i-tej statystyki pozycyjnej (z 10) rozkladu jednostajnego U(0,1)

#Metoda MC wyznaczenia kwantyli (bootstrap parametryczny)

#B1. Generujemy B=10000 probek licznosci 10 z rozkladu jednostajnego U(0,1)
#i porzadkujemy od wartosci najmniejszej do najwiekszej  
os10_unif <- matrix(NA,nrow=10000,ncol=10)

for(i in 1:10000){
  
  os10_unif[i,] <- sort(runif(10))
}

head(os10_unif)

#Dla kazdej statystyki mamy probe licznosci 10 000 z jej rozkladu - 
#wartosci w kolumnach  macierzy os10_unif.

#B2. Obliczamy kwantyle rzedu alpha=0.05 kazdej z dziesieciu statystyk pozycyjnych, 
#sa  to progi odrzucen dla kolejnych qi
osl_unif <- apply(os10_unif,2,function(x) quantile(x,alpha))

#Wyniki (q-value i progi orzucen)
plot(rev(qi),ylim=c(0,1), pch=19)
points(rev(osl),type="l",col=2 )
points(rev(osl_unif),type="l",col="pink")
#Jakie jest optymalne r?


#6.Oblicz poziomy zwrotu x20 i x50, wykorzystujac optymalna liczbe statystyk
#z punktu 5. 
#Porownaj wyniki z tymi z Analizy 1 - otrzymanymi, przy uzyciu jednej
#maksymalnej temperatury w roku. 
#=========================================================================
#Mozesz wyestymowac ponownie parametry rozkladu GEVr, przy optymalnym r
# i dalej skorzystac z funkcji 'gevrRl' (jak w  pt. 3))

#os_r <- os10[,1:r]
#fit_r <- gevrFit(os_r,method="mle")






