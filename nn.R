#nearest neighbour analyses
library(dplyr)
library(sp)
library(rgdal)
library(raster)
library(maptools)
library(spatstat)
library(spdep)
library(ggplot2)
library(ggthemes)

species <- read.csv("Data/species.csv")

t6.sp <- readOGR("Data/GIS/T6Fix.shp")
t6.sp <- readOGR("Data/GIS/T5.shp")
t6.sp <- readOGR("Data/GIS/T4.shp")

t6.ppp <- as(t6.sp, "ppp")
summary(t6.ppp)
t6.ppp$marks$Width <- as.numeric(as.character((t6.ppp$marks$Width)))
t6.ppp$marks$Heigh <- as.numeric(as.character((t6.ppp$marks$Heigh)))
marks.data <- as.data.frame(marks(t6.ppp))
marks.data <- left_join(marks.data, species, "Species")
marks.data <- marks.data %>% dplyr::select(sp_simp, Type, Species, Width, Height)
marks(t6.ppp) <- marks.data

nn <- nndist(t6.ppp)
tnn <- cbind(nn, marks.data)

nn

#allnn <- as.data.frame(matrix(data = NA, nrow = 0, ncol = 6))
allnn <- rbind(allnn, tnn)

n <- nnmark(t6.ppp, at = "points")
plot(n)
n

al <- cbind(tnn, n)

write.csv(al, "neighbours.csv")

hist(allnn$nn)
ggplot(allnn, aes(nn)) + geom_histogram(binwidth = 0.25, color = "black") + xlab("Distance to Nearest Neighbour (m)") + ylab("Frequency") + theme_Publication()


plot(nn ~ sp.simp, data = allnn)
a1 <- aov(nn ~ sp.simp, data = allnn)
a1
summary(a1)
TukeyHSD(a1)

m1 <- lm(nn ~ Width * sp.simp, data = allnn)
summary(m1)
aov(m1)

b <- bw.relrisk(t6.ppp)

split6 <- split(t6.ppp)
t6cacti <- split6$cactus 
t6cacti <- as(t6cacti, "ppp")
t6shrub <- split6$shrub
t6shrub <- as(t6shrub, "ppp")  
nncact <- nndist.default(t6cacti, t6shrub)
hist(nncact)

plot(t6cacti)
