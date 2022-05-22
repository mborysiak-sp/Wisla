#Laboratorium: POT-automatyczny wybor progu (7.04.2021)
library(fields)
library(eva)

#1. Zaladuj maksima dzienne ze swojej stacji
#data <- ...
load(file="Wisla.Rdata")
data <- max10[,5]
data <- data[!is.na(data)]

#ustal progi na poziomie kwantyli ze zbioru A
A <- c(seq(0.75,0.97,by=0.02),seq(0.971,0.985,by=0.001))
th <- A[27]

#kolory do wykresu
d <- length(A)
my.col <- topo.colors(d,rev=TRUE)

#zrob wykres  danych i progow/ zapisz w stworzonym do tego celu
#katalogu ,,Figures''

png("Figures/wisla_plot.png", height=600,width=800)

plot(data,type="h", xlab=NA,ylab="max_day" )
for(i in 1:d){
  abline(h=th[i],col=my.col[i])
}
dev.off()

#2. Dla nadwyzek dla progow ze zbioru th wyestymuj parametry rozkladu GPD
#i wykonaj test AD - wykorzystaj funkcje ,,gpdAd'' z biblioteki 'eva'
#------------------------------------------------------------------------
p_value <- c()

for(i in 1:d){
  
  datau <- data[data>th[i]]
  fit <- gpdAd(datau)
  
  p_value[i] <- fit$p.value
  
  print(i)
}

#Narysuj otrzymane p_value - wykres zapisz od razu w katalogu Figures

png("Figures/wisla_p_value.png", height=600,width = 800)
plot(p_value, pch=19, xlab="kolejne progi")
abline(h=0.05,col=2)
dev.off()

#3. Oblicz q_value i skorygowane p_value z metody ForwardStop (p_FS).
#------------------------------------------------------------------
#q_value
Yi <- -log(1-p_value)
ri <- rev(1:d)
Zi0 <- Yi/ri
Zi <- cumsum(Zi0)
qi <- 1-exp(-Zi)

z <- pSeqStop(p_value)
p_FS <- z$ForwardStop

#Zadanie. Na jednym wykresie przedstaw p_value, q_value i p_FS.
#Wykres zapisz od razu w katalogu Figures.
png("q_value.png", height=600, width=800)
par(mfrow=c(1,1))
plot(p_value, pch=19, xlab="kolejne progi")
points(qi, col="green")
points(p_FS, pch=23, col="pink")
abline(h=0.05, col=2)
dev.off()

#4. Korzystajac z q_value i p_FS wyznacz liczbe kmax poczatkowych hipotez do odrzucenia.
#---------------------------------------------------------------------------------
#warunek odrzucanai z wykorzystaniem q-value (zobacz na wykresie)
alpha <- 0.05
os <- (1:d)/d
osl <- alpha*os

png("Figures/wisla_q_wybor.png", height=600,width = 800)
plot(qi, pch=19, col="green", xlab="kolejne progi",
     ylab="q_value")
points(osl,type="l",col=2)
dev.off()


#zrob wykres p_FS z  pozioma linia na wysokosci alpha=0.05. Odczytaj z niego kmax

plot(p_FS, ylim=c(0, 0.03))
abline(h=alpha, col=2)

#5. Wykorzystujac wyniki z pt. 4, wyznacz optymalny prog i poziomy zwrotu x20, x50.
#----------------------------------------------------------------------------------
#jesli nie odzrzucamy zadnej hipotezy przyjmij najmnizszy prog: u=th[1],
#jesli odrzucamy wszystkie przyjmij najwyzszy prog: u=th[d],
#w pozostalych przypadkach przyjmij: u=th[kmax+1]

#Dla danego progu oblicz poziom zwrotu x20 i x50 (mozesz to zrobic znanymi juz metodami
# lub wykorzystac  rozwiazanie z biblioteki eva)

#Podajemy prog i pamietamy aby w zaleznosci od wykorzystanych danych,
#ustalic odpowiednie k

fit <- gpdFit(data, threshold = th[1])
rlk <- gpdRl(fit, period = k, method = "profile",plot=FALSE)$Estimate

length(A)
qi(24)
quantile(data, A[27])