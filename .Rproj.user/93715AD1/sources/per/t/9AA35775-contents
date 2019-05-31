#join counts!
library(dplyr)
library(sp)
library(rgdal)
library(raster)
library(maptools)
library(spatstat)
library(spdep)

species <- read.csv("Data/species.csv")

data <- readOGR("Data/Field Edits/T6_4_11_2017.shp")
data5 <- readOGR("Output/T5_J.shp")
data4 <- readOGR("Output/T4_J.shp")
data.sp$Width <- as.numeric(as.character((data$Width)))

data$sp_simp <- as.factor(data$sp_simp)
data$Type <- as.factor(data$Type)
data6 <- data
jc.results <- matrix(data = NA, nrow = 1, ncol = 4)
r <- jc(data4, 2, data4$sp_simp)
jc.results <- rbind(r, jc.results)
jc.results.T4 <- data.frame(names = row.names(jc.results), jc.results)

jc.results <- matrix(data = NA, nrow = 1, ncol = 4)
r <- jc(data5, 2, data5$sp_simp)
jc.results <- rbind(r, jc.results)
jc.results.T5 <- data.frame(names = row.names(jc.results), jc.results)

jc.results <- matrix(data = NA, nrow = 1, ncol = 4)
r <- jc(data6, 2, data6$sp_simp)
jc.results <- rbind(r, jc.results)
jc.results.T6 <- data.frame(names = row.names(jc.results), jc.results)

jc.results.T4 <- mutate(jc.results.T4, p = 2*pnorm(-abs(`z.value`)), transect = 4)

jc.results.T5 <- mutate(jc.results.T5, p = 2*pnorm(-abs(`z.value`)), transect = 5)
jc.results.T6 <- mutate(jc.results.T6, p = 2*pnorm(-abs(`z.value`)), transect = 6)

results <- rbind(jc.results.T4, jc.results.T5, jc.results.T6)

write.csv(jc.results.T6, "results.csv")
write.csv(results, "JoinCount2neigh.csv")







nn<-2
weight <- function(t6.sp, nn){
  t6nb <- knn2nb(knearneigh(data, k = nn))
  summary(t6nb)
  is.symmetric.nb(t6nb)
  t6sym <- make.sym.nb(t6nb)
  is.symmetric.nb(t6sym)
  t6listw <- nb2listw(t6sym)
}

jc <- function(t6.sp, nn, x){
t6nb <- knn2nb(knearneigh(t6.sp, k = nn))
summary(t6nb)
is.symmetric.nb(t6nb)
t6sym <- make.sym.nb(t6nb)
is.symmetric.nb(t6sym)
t6listw <- nb2listw(t6sym)
joincount.multi(x, t6listw)
}
#joincount.mc(t6.sp$Type, t6listw, 99)

#calculate p-values




data6 <- readOGR("Output/T6_J.shp")
data5 <- readOGR("Output/T5_J.shp")
data4 <- readOGR("Output/T4_J.shp")

t6.sp <- data5
nn <- 3

t6nb <- knn2nb(knearneigh(t6.sp, k = nn))
is.symmetric.nb(t6nb)
t6sym <- make.sym.nb(t6nb)
is.symmetric.nb(t6sym)
t6listw <- nb2listw(t6sym)

jctest <- joincount.test(data5$sp_simp, t6listw, alternative = "two.side")
j <- as.data.frame(jctest)









weight(data, 2)
data$Width <- as.numeric(as.character(data$Width))
m <- moran.mc(data$Width, t6listw, 9999)
moran.mc(data$Height, t6listw, 9999)
plot(m)
m

plot(localmoran(data$Width, t6listw))

geary.mc(t6.sp$Width, t6listw, 999)
geary.mc(t6.sp$Height, t6listw, 999)

globalG.test(t6.sp$Width, t6listw, zero.policy = NULL, alternative = "less")
globalG.test(t6.sp$Height, t6listw, zero.policy = NULL, alternative = "less")

#jc tests
joincount.test()