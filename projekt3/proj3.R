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
missing_wind <-c()
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
    missing_wind <- c(missing_wind, sum(is.na(station_winds)))
  }
}
median_missing <- median(missing)
median(missing_wind)
mean_missing <- mean(missing_wind)
for(i in 1:length(missing)){
  if(missing[i] > median_missing || missing_wind[i] > mean_missing){
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
load(file="D:/IT/Nauka/R/DATA/Daily_Max_August.Rdata")
library(copula)
library(VineCopula)
library(ggplot2)
library(ggExtra)
library(MASS)

output_data <- list() # komplet danych stacji i estymacji
cop_stats <-c()
it <- 1
for(station in formatted){
  X1 <- station[[2]]$maxd_temp
  X2 <- station[[2]]$maxd_wind
  X <- data.frame(X1,X2)
  X <- X[complete.cases(X),]
  V <- pobs(X)
  colnames(V) <- c("V1","V2")
  cop.npar <- BiCopSelect(V[,1],V[,2], selectioncrit="AIC", se=TRUE)
  norm.cop <- BiCop(family = cop.npar$family, par = cop.npar$par, par2=cop.npar$par2)
  cop_type <- norm.cop$familyname
  
  cop_stats <- c(cop_stats, norm.cop$familyname)
  cor_kendall <- cor(X, method="kendall")[1,2]

  
  N <- 30 * 11
  temp <- c()
  wind <- c()
  for (i in 1 : 100) {
    V <- BiCopSim(N,norm.cop)
    Z.npar <- cbind(quantile(X$X1, V[,1], na.rm=TRUE), quantile(X$X2, V[,2]))
    
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
  tryCatch({
    fit_temp <- evir::gev(temp)
    ests_t <- fit_temp$par.ests
    t_x20 <- evir::qgev(1-1/20, ests_t[[1]], ests_t[[3]], ests_t[[2]])
    t_x50 <- evir::qgev(1-1/50, ests_t[[1]], ests_t[[3]], ests_t[[2]])
  }, error = function(cond){
    t_x20 <- NA
    t_x50 <- NA
  })
  
  tryCatch({
    fit_wind <- evir::gev(wind)
    ests_w <- fit_wind$par.ests
    w_x20 <- evir::qgev(1-1/20, ests_w[[1]], ests_w[[3]], ests_w[[2]])
    w_x50 <- evir::qgev(1-1/50, ests_w[[1]], ests_w[[3]], ests_w[[2]])
  }, error = function(cond){
    w_x20 <- NA
    w_x50 <- NA
  })
  out_results <- data.frame(cop_type, cor_kendall, t_x20, t_x50, w_x20, w_x50)
  station <- append(station, out_results)
  output_data[[length(output_data) + 1]] <- station
  print(round(100 * it / length(formatted), 1))
  it <- it + 1
  
}
save(output_data, file="D:/IT/Nauka/R/DATA/Daily_Max_August_Output.Rdata")
occurence <- data.frame(table(unlist(cop_stats)))
colnames(occurence) <- c("Rodzaj kopuly", "liczba wystapien")

