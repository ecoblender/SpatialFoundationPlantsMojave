#cactus vs shrub analyses
library(dplyr)
library(sp)
library(rgdal)
library(raster)
library(maptools)
library(spatstat)
library(spdep)


species <- read.csv("Data/species.csv")

t6.sp <- readOGR("Output/T6_J.shp")
t6.sp <- readOGR("Output/T5_J.shp")
t6.sp <- readOGR("Output/T4_J.shp")

t6.ppp <- as(t6.sp, "ppp")
marks.data <- as.data.frame(marks(t6.ppp))
marks.data <- marks.data %>% dplyr::select(Width)
marks.data <- marks.data %>% dplyr::select(sp_simp)
marks(t6.ppp) <- marks.data

plot(Gcross(t6.ppp, "cactus"))
plot(markconnect(t6.ppp, "shrub", "cactus"))
plot(markcorr(t6.ppp))

plot(density(split(t6.ppp)))
E <- envelope(t6.ppp, Kcross, nsim = 39, i = "cactus", j = "shrub")
plot(E, main = "test of marked Poisson model")

plot(smooth.ppp(t6.ppp))
plot(split(cut(t6.ppp, breaks = 5)))

plot(Emark(t6.ppp))

ppm(t6.ppp, ~marks, MultiStrauss(1))


#test if pattern is poisson
magicLand <- rpoispp(lambda = 0.5, win = owin(c(0, 10), c(0, 10)))
plot(magicLand)
# create a random map in a 10 by 10 window
ranMap <- rpoispp(lambda = 10, win = owin(c(0, 10), c(0, 10)))
# divide the window into six quadrants calculate chi-square
qTest <- quadrat.test(split(t6.ppp, nx = 4, ny = 10))

# visualize the result which shows observed count, expected count and SD.
plot(qTest)
qTest

count(marks.data, sp_simp)
is.multitype(t6.ppp)
at <- alltypes(t6.ppp, "Kcross")  #this takes awhile
plot(at)

#to do this we will need to remove ones with only 1

at <- Kdot(t6.ppp, "lar.tri", envelope = T)
plot(at)
lk <- localK(t6.ppp)
plot(lk)


lK.60 <- localK(t6.ppp, rvalue = 2, verbose = T)  #takes a long time!!!! Set verbose=T to track progress
bei.lk.60 <- unmark(t6.ppp) %mark% lK.60
lk.bei.smooth <- smooth.ppp(bei.lk.60, sigma = 2)
plot(lk.bei.smooth, col = topo.colors(128), main = "smoothed neighbourhood density")
contour(lk.bei.smooth, add = TRUE)


n <- dnearneigh(t6.sp, 0.5, 0.5)

k <- Kinhom(t6.ppp)
E <- envelope(t6.ppp, Kinhom, nsim = 99)
plot(E)
plot(pcf(E))
D <- density(t6.ppp)
Y <- rmpoispp(D, types=names(D))
plot(Y)
plot(D)
r.v <- seq(0,5, by = 0.05)
E <- envelope(t6.ppp, Kinhom, r = r.v,
         simulate=expression(rmpoispp(D, types=names(D))))
plot(E)

marks.data <- as.data.frame(marks(t6.ppp))
marks.data <- marks.data %>% dplyr::select(Type)
marks(t6.ppp) <- marks.data
marks.data <- as.data.frame(marks(t5.ppp))
marks.data <- marks.data %>% dplyr::select(Type)
marks(t5.ppp) <- marks.data
marks.data <- as.data.frame(marks(t4.ppp))
marks.data <- marks.data %>% dplyr::select(Type)
marks(t4.ppp) <- marks.data


g6 <- Gcross(t6.ppp, "cactus", "shrub")
plot(g6)
g5 <- Gcross(t5.ppp, "cactus", "shrub")
plot(g5)
g4 <- Gcross(t4.ppp, "cactus", "shrub")
plot(g4)

count(marks.data, Type)

lambda6 <- density(t6.ppp, bw.ppl)
pcf6 <- pcfcross.inhom(t6.ppp, "shrub", "cactus", lambda)
plot(pcf6)

lambda5 <- density(t5.ppp, bw.ppl)
pcf5 <- pcfcross.inhom(t5.ppp, "shrub", "cactus", lambda)
plot(pcf5)

lambda4 <- density(t4.ppp, bw.ppl)
pcf4 <- pcfcross.inhom(t4.ppp, "shrub", "cactus", lambda)
plot(pcf4)
