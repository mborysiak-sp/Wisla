library(maps)

coords <- read.csv("data/coord233_alt.csv")
city <- coords[coords$place=="WISLA",]
station_id = city['station']$station

pdf("data/proj2/mapa.png")
poland <- map('world', 'poland', fill=T, col='gray');
points(city[c('lon', 'lat')], pch=19, col=2);
dev.off()