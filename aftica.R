library(readr)
library(ggmap)
library(dplyr)
library(spatstat)
library(maptools)
library(rgdal)
library(sp)
library(splancs)
library(RColorBrewer)

set.seed(798487)

conflict.data <- read_csv(
  "africa.csv")
conflict.algeria = subset(conflict.data, conflict.data$COUNTRY == "Uganda")
conflict.algeria = subset(conflict.algeria, conflict.algeria$YEAR %in% c(1997, 2015))
bbox <- make_bbox(LONGITUDE,
                  LATITUDE,
                  data=conflict.algeria,
                  f=0.2)
conflict.algeria$anio = as.factor(conflict.algeria$YEAR)
anio = conflict.algeria$anio
africa <- get_map(bbox, source="stamen")
ggmap(africa) +
  geom_point(aes(x=LONGITUDE,
                 y=LATITUDE,
                 shape=anio),
             data=conflict.algeria,
             alpha=0.5) +
  xlim(29, 35.5) +
  ylim(-2, 4.5) +
  scale_color_gradient(limits=c(1997, 2015),
                       low="orangered1",
                       high="red4")

africa.marks = conflict.algeria[, "YEAR"]
spatial_df <- SpatialPointsDataFrame(coords = conflict.algeria[, c("LONGITUDE", "LATITUDE")], proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),
                                              data = africa.marks)

africa.transed <- spTransform(spatial_df, CRS("+init=epsg:3839"))
africa.transed$Easting <- coordinates(africa.transed)[, 1]
africa.transed$Northing <- coordinates(africa.transed)[, 2]

plot(x = africa.transed$Easting, y = africa.transed$Northing, ylab = "Norte", 
     xlab = "Este", main = "Ubicación (transformada) de Incidentes violentos en Ruanda")


africa.chull <- convexhull.xy(x = africa.transed$Easting, y = africa.transed$Northing)
africa.ppp <- ppp(x = africa.transed$Easting, y = africa.transed$Northing, africa.chull,
              marks = as.factor(africa.transed$YEAR))

unitname(africa.ppp) <- c("meter", "meters")
p1 <- plot.ppp(africa.ppp, main="Distribucion de Acontecimientos Violentos en Ruanda")  # By default, it will plot the first mark that you assigned
legend(max(coords(africa.ppp))[1] + 1000, mean(coords(africa.ppp))[2], pch = p1, legend = names(p1), 
       cex = 0.5)

anios = split.ppp(africa.ppp)
africa.1997 = anios$`1997`
africa.1997.sp <- as(africa.1997, "SpatialPoints")
africa.1997.spu <- elide(africa.1997.sp, scale=TRUE, unitsq=TRUE)
africa.1997.pppu <- as(africa.1997.spu, "ppp")

r = seq(0, sqrt(2)/6, by = 0.005)

env97 <- envelope(africa.1997.pppu, fun=Gest,
                   r=r, nrank=2, nsim=99)

plot(env97, main="Funcion G para conflictos durante 1997")


anios = split.ppp(africa.ppp)
africa.2015 = anios$`2015`
africa.2015.sp <- as(africa.2015, "SpatialPoints")
africa.2015.spu <- elide(africa.2015.sp, scale=TRUE, unitsq=TRUE)
africa.2015.pppu <- as(africa.2015.spu, "ppp")

env2015 <- envelope(africa.2015.pppu, fun=Gest,
                  r=r, nrank=2, nsim=99)

plot(env2015, main="Funcion G para conflictos durante 2015")



rf = seq(0, sqrt(2)/6, by = 0.0005)
Fenv97 = envelope(africa.1997.pppu, fun=Fest, r=rf,
                  nrank=2, nsim=999)
plot(Fenv97, main="Funcion F para conflictos durante 1997")
Fenv2015 = envelope(africa.2015.pppu, fun=Fest, r=rf,
                  nrank=2, nsim=999)
plot(Fenv2015, main="Funcion F para conflictos durante 2015")

Kenv97 = envelope(africa.1997.pppu, fun=Kest, r=rf,
                  nrank=2, nsim=999)
plot(Kenv97, main="Funcion K para conflictos durante 1997")
Kenv2015 = envelope(africa.2015.pppu, fun=Kest, r=rf,
                  nrank=2, nsim=999)
plot(Kenv2015, main="Funcion K para conflictos durante 2015")

# Calculo de la intensidad
mserw = bw.diggle(africa.1997.pppu)
bw = as.numeric(mserw)
bw
plot(density(redwood, bw=bw, kernel='gaussian'), main="Densidad con kernel Gaussiano para 1997")

mserw = bw.diggle(africa.2015.pppu)
bw = as.numeric(mserw)
bw
plot(density(redwood, bw=bw, kernel='gaussian'), main="Densidad con kernel Gaussiano para 2015")


