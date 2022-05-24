library(fields)
library(eva)

loadValues <- function() {
  path_to_files <- "data/temp.stations-all"
  coords <- read.csv("data/coord233_alt.csv")

  L <- as.list(list.files(path = path_to_files, pattern = ".*7.csv"))
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

  k20 <- 20
  k50 <- 50
  fit <- gpdFit(data, threshold = u)
  r <- tryCatch({
    rl20 <- gpdRl(fit, period = k20, method = "profile", plot = FALSE)$Estimate
  }, error = function(cond) {
    rl20 <- NA
  })
  r <- tryCatch({
    rl50 <- gpdRl(fit, period = k50, method = "profile", plot = FALSE)$Estimate
  }, error = function(cond) {
    rl50 <- NA
  })

  x20pp1[[length(x20pp1) + 1]] <- rl20
  x50pp1[[length(x50pp1) + 1]] <- rl50

  row <- data.frame(station[[1]]$station, u, rl20, rl50)
  csv_fname <- "data/proj2/resultLipiecPOT.csv"
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

for (i in 1:(min(length(stations), length(x20pp1)))) {
  stations[[i]][[3]] <- x20pp1[[i]]
}

thresholds_sorted <- sort(unlist(x20pp1), decreasing = FALSE)
stations_sorted <- stations[order(sapply(stations, "[[", 3))]
png("Figures/proj2/mapa_progi_pot.png", width = 2048, height = 2048)
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
}
dev.off()
png("Figures/proj2/legenda_progi_pot.png", width = 2048, height = 2048)
color.bar(rbPal(step), thresholds_sorted[[1]], thresholds_sorted[[length(thresholds_sorted)]], nticks = step)
dev.off()

#podpunkt 2 (BMM)
library(eva)
library(gamlss)
load(file = "data/proj2/Selected_Temp_August.Rdata")
x20pp2 <- list()
x50pp2 <- list()
r_list <- list()
parGEVr <- NULL
for (ii in 1:length(formatted)) {
  t1 <- Sys.time()
  temps <- formatted[[ii]][[2]][complete.cases(formatted[[ii]][[2]]),]
  r <- 10
  os10 <- matrix(NA, nrow = length(2008:2018), ncol = r)
  count <- 1
  for (i in 2008:2018) {
    datar <- temps[temps$year == i, 5]
    datar <- datar[!is.na(datar)]
    x <- sort(datar, decreasing = TRUE)
    y <- x[1:r]
    os10[count,] <- y
    count <- count + 1
  }
  parGEVr <- matrix(NA, nrow = 10, ncol = 3)
  B <- 100
  r <- tryCatch({
    test10 <- gevrSeqTests(os10, method = "pbscore", bootnum = B)

    p_FS <- rev(test10$ForwardStop)

    r <- length(Filter(function(x) x > 0.05, p_FS))
  }, error = function(cond) {
    return(NA)
  })
  x20 <- NA
  x50 <- NA
  if (!is.na(r)) {
    if (r > 0) {
      osr <- os10[, 1:r]
      fit <- gevrFit(osr, method = "mle")

      parGEVr[r,] <- fit$par.ests

      x20 <- gevrRl(fit, 20)$Estimate
      x50 <- gevrRl(fit, 50)$Estimate
      if (x20 > 40) {
        x20 <- NA
      }
      if (x50 > 40) {
        x50 <- NA
      }
    }
  }

  x20pp2[[length(x20pp2) + 1]] <- x20
  x50pp2[[length(x50pp2) + 1]] <- x50
  r_list[[length(r_list) + 1]] <- r

  print(round(100 * ii / length(formatted), 1))
  t2 <- Sys.time()
  print(t2 - t1)
}
save(x20pp2, x50pp2, file = "data/proj2/BMM_Temp_August100.Rdata")
load(file = "data/proj2/BMM_Temp_August.Rdata")

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

for (i in 1:(min(length(stations), length(x50pp2)))) {
  stations[[i]][[3]] <- x50pp2[[i]]
}

thresholds_sorted <- sort(unlist(x50pp2), decreasing = FALSE)
stations_sorted <- stations[order(sapply(stations, "[[", 3))]
png("Figures/proj2/mapa_bmm50.png", width = 2048, height = 2048)
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
  if (is.na(station[[3]])) {
    points(lon, lat, pch = 20, cex = 10)
  } else {
    points(lon, lat, pch = 20, col = station[[4]], cex = 10)
  }
}
dev.off()
png("Figures/proj2/legenda_bmm50.png", width = 2048, height = 2048)
color.bar(rbPal(step), thresholds_sorted[[1]], thresholds_sorted[[length(thresholds_sorted)]], nticks = step)
dev.off()
