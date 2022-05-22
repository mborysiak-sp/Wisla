library(fields)
library(eva)

loadValues <- function() {
  path_to_files <- "data/temp.stations-all"
  coords <- read.csv("data/coord233_alt.csv")

  L <- as.list(list.files(path = path_to_files, pattern = ".*8.csv"))
  L <- paste0(path_to_files, "\\", L)
  data0 <- lapply(L, read.csv)
  n <- length(data0)
  temps <- list()
  missing <- c()
  min_temp <- 0
  max_temp <- 40
  data_max <- 6 * 24 * 31 * 11
  for (i in 1:nrow(coords)) {
    station_temps <- c()
    for (x in 1:n) {
      if (exists(coords[i,]$station, data0[[x]])) {
        station_temps <- c(station_temps, get(coords[i,]$station, data0[[x]]))
      }
    }
    if (length(station_temps) == data_max) {
      station_temps[station_temps > max_temp] = NA
      station_temps[station_temps < min_temp] = NA
      temps[[length(temps) + 1]] <- list(station_temps, coords[i,])
      missing <- c(missing, sum(is.na(station_temps)))
    }
  }
  median_missing <- median(missing)

  for (i in 1:length(missing)) {
    if (missing[i] > median_missing) {
      temps <- temps[-i]
    }
  }

  datetime <- c()

  for (i in 1:n) {
    datetime <- c(datetime, as.character(data0[[i]]$datetime))
  }
  formatted <- list()

  library(tidyr)

  for (i in 1:length(temps)) {
    formatted[[length(formatted) + 1]] <- list(
      temps[[i]][[2]],
      data.frame(date = as.Date(datetime), max10 = temps[[i]][[1]]))
    rownames(formatted[[i]][[2]]) <- c()
    formatted[[i]][[2]] <- separate(
      formatted[[i]][[2]],
      date,
      c("year", "mth", "day"),
      convert = TRUE
    )
    formatted[[i]][[2]] <- data.frame(datetime = datetime, formatted[[i]][[2]])
    print(round(100 * i / length(temps), 1))
  }

  save(formatted, file = "data/proj2/Selected_Temp_August.Rdata")
  load(file = "data/proj2/Selected_Temp_August.Rdata")
}

# draw_map <- function(stations, thresholds) {
#   library(maps)
#   stations$th <- thresholds
#   stations_sorted <- stations[order(stations$th),]
#   pdf("data/proj2/mapazlegenda.png")
#   map('world', 'poland', fill = T, col = 'gray')
#   rbPal <- colorRampPalette(c('blue', 'red'))
#   stations_sorted$col <- rbPal(50)[as.numeric(cut(ths$th, breaks = 50))]
#   for (station in stations_sorted) {
#     lon <- station[[1]]$lon
#     lat <- station[[1]]$lat
#     points(lon, lat, pch = 19, col = station$col)
#   }
#   dev.off()
#
# }

#
# iterate_over_stations <- function() {
#   load(file = "data/proj2/Selected_Temp_August.Rdata")
#
#   for (station in formatted)
#   {
#     data <- station[2]$max10[, 5]
#     data <- data[!is.na(data)]
#
#     A <- c(seq(0.85, 0.97, by = 0.02), seq(0.971, 0.985, by = 0.001))
#     th <- quantile(data, A)
#
#     k20 <- 20 * 31 * 24 * 6
#     k50 <- 50 * 31 * 24 * 6
#     th[1]
#     fit <- gpdFit(data, threshold = th[1])
#     rl20 <- gpdRl(fit, period = k20, method = "profile", plot = FALSE)$Estimate
#     rl50 <- gpdRl(fit, period = k50, method = "profile", plot = FALSE)$Estimate
#
#     # x20pp1[[length(x20pp1) + 1]] <- rl20
#     # x50pp1[[length(x50pp1) + 1]] <- rl50
#     sprintf("Station %s", station[1]$station)
#     sprintf("Fit %f", rl20)
#     sprintf("Rlk %f", rl50)
#   }
# }
#
# calculate_hiphotesis_to_reject_count <- function(d) {
#   alpha <- 0.05
#   os <- (1:d) / d
#   osl <- alpha * os
# }
#
# calculate_q_value <- function(p_value, d) {
#   #3. Oblicz q_value i skorygowane p_value z metody ForwardStop (p_FS).
#   #------------------------------------------------------------------
#   #q_value
#   Yi <- -log(1 - p_value)
#   ri <- rev(1:d)
#   Zi0 <- Yi / ri
#   Zi <- cumsum(Zi0)
#   qi <- 1 - exp(-Zi)
#   qi
# }
#
# calculate_corrected_p_value <- function(p_value) {
#   #3. Oblicz q_value i skorygowane p_value z metody ForwardStop (p_FS).
#   #------------------------------------------------------------------
#   #q_value
#   z <- pSeqStop(p_value)
#   p_FS <- z$ForwardStop
#   p_FS
# }
#
# calculate_GPD_parameters <- function(data, th) {
#   #2. Dla nadwyzek dla progow ze zbioru th wyestymuj parametry rozkladu GPD
#   #i wykonaj test AD - wykorzystaj funkcje ,,gpdAd'' z biblioteki 'eva'
#   #------------------------------------------------------------------------
#   p_value <- c()
#
#   for (i in 1:d) {
#     datau <- data[data > th[i]]
#     fit <- gpdAd(datau)
#     p_value[i] <- fit$p.value
#     print(i)
#   }
# }

load(file = "data/proj2/Selected_Temp_August.Rdata")

x20pp1 <- NULL
x50pp1 <- NULL
thresholds <- NULL
for (station in formatted)
{
  data <- station[[2]]$max10
  data <- data[!is.na(data)]
  A <- c(seq(0.85, 0.97, by = 0.02), seq(0.971, 0.985, by = 0.001))
  th <- quantile(data, A)
  d <- length(A)

  datau <- NULL
  p_value <- c()
  for (i in 1:length(A)) {
    datau <- data[data > th[i]]
    r <- tryCatch({
      fit <- gpdAd(as.numeric(datau))
    }, error = function(cond) {
      return(NA)
    })
    p_value[i] <- fit$p.value
  }

  z <- pSeqStop(p_value)
  p_FS <- z$ForwardStop
  k <- c()
  for (i in 1:length(p_FS)) {
    if (p_FS[i] >= 0.05) {
      k <- c(k, p_FS[i])
    }
  }

  retinf <- 0
  if (is.null(k)) {
    u <- th[d]
  }else if (length(k) == d) {
    u <- th[1]
    retinf <- 1
  }else {
    u <- th[match(max(k), p_FS)]
    retinf <- 2
  }

  thresholds[[length(thresholds) + 1]] <- u

  k20 <- 20 * 31 * 24 * 6
  k50 <- 50 * 31 * 24 * 6
  fit <- gpdFit(data, threshold = u)
  r <- tryCatch({
    rl20 <- gpdRl(fit, period = k20, method = "profile", plot = FALSE)$Estimate
  }, error = function(cond) {
    return(NA)
  })
  r <- tryCatch({
    rl50 <- gpdRl(fit, period = k50, method = "profile", plot = FALSE)$Estimate
  }, error = function(cond) {
    return(NA)
  })

  x20pp1[[length(x20pp1) + 1]] <- rl20
  x50pp1[[length(x50pp1) + 1]] <- rl50

  row <- data.frame(station[[1]]$station, u, rl20, rl50)
  csv_fname <- "data/proj2/resultPOT.csv"
  write.table(row, file = csv_fname, sep = ",",
              append = TRUE, quote = FALSE,
              col.names = FALSE, row.names = FALSE)
}

color.bar <- function(lut, min, max = -min, nticks = 11, ticks = seq(min, max, len = nticks), title = '') {
  margin <- 20
  par(mar = c(margin, margin, margin, margin))
  scale <- (length(lut) - 1) / (max - min)
  plot(c(0, 10), c(min, max), type = 'n', bty = 'n', xaxt = 'n', xlab = '', yaxt = 'n', ylab = '', main = title)
  axis(2, ticks, las = 1, cex.axis = 5, font = 2,)
  for (i in 1:(length(lut) - 1)) {
    y <- (i - 1) / scale + min
    rect(0, y, 10, y + 1 / scale, col = lut[i], border = NA)
  }
}

library(maps)

stations <- formatted

for (i in 1:(min(length(stations), length(thresholds)))) {
  stations[[i]][[3]] <- thresholds[[i]]
}

thresholds_sorted <- sort(unlist(thresholds), decreasing = FALSE)
stations_sorted <- stations[order(sapply(stations, "[[", 3))]
png("Figures/proj2/mapa.png", width = 2048, height = 2048)
map('world', 'poland', fill = T, col = 'gray')
rbPal <- colorRampPalette(c('blue', 'yellow', 'red'))
step <- 7
col <- rbPal(step)[as.numeric(cut(thresholds_sorted, breaks = step))]
for (i in 1:(length(stations_sorted))) {
  stations_sorted[[i]][[4]] <- col[i]
}
for (station in stations_sorted) {
  lon <- station[[1]]$lon
  lat <- station[[1]]$lat
  points(lon, lat, pch = 20, col = station[[4]], cex = 10)
  # par(ps=5)
  # text(lon, lat, station[[3]])
}
dev.off()
png("Figures/proj2/legenda.png", width = 2048, height = 2048)
color.bar(rbPal(step), thresholds_sorted[[1]], thresholds_sorted[[length(thresholds_sorted)]], nticks = step)
dev.off()

#podpunkt 2 (BMM)
library(eva)
library(gamlss)
load(file = "data/proj2/Selected_Temp_August.Rdata")
x20pp2 <- list()
x50pp2 <- list()
parGEVr <- NULL
for(ii in 1:length(formatted)){
  t1<-Sys.time()
  temps <- formatted[[ii]][[2]][complete.cases(formatted[[ii]][[2]]), ]
  r <- 10
  os10 <- matrix(NA,nrow=length(2008:2018),ncol=r)
  count <- 1
  for(i in 2008:2018){
    datar <- temps[temps$year == i, 5]
    datar <- datar[!is.na(datar)]
    x <- sort(datar,decreasing = TRUE)
    y <- x[1:r]
    os10[count,] <- y
    count <- count + 1
  }
  B <- 5
  r <- tryCatch({
    test10 <- gevrSeqTests(os10, method="pbscore", bootnum=B)

    p_FS <- rev(test10$ForwardStop)

    length(Filter(function(x) x > 0.05, p_FS))
    }, error = function(cond) {
    return(NA)
  })
  x20 <- NA
  x50 <- NA
  if(!is.na(r)){
    if(r > 0){
      osr <- os10[,1:r]
      fit <- gevrFit(osr,method="mle")

      parGEVr[r,] <- fit$par.ests

      x20 <- gevrRl(fit,20)$Estimate
      x50 <- gevrRl(fit,50)$Estimate
    }
  }
  x20pp2[[length(x20pp2) + 1]] <- x20
  x50pp2[[length(x50pp2) + 1]] <- x50


  print(round(100 * ii / length(formatted), 1))
  t2<-Sys.time()
  print(t2-t1)
}
save(x20pp2, x50pp2, file="data/proj2/BMM_Temp_August.Rdata")

