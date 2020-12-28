library(ade4)
library(vegan)

library(dplyr)
library(sp)
library(rgdal)
library(raster)
library(maptools)
library(spatstat)
library(spdep)

t6.sp <- readOGR("Output/T6_J.shp")
t6.sp <- readOGR("Output/T5_J.shp")
t6.sp <- readOGR("Output/T4_J.shp")

t6.ppp <- as(t6.sp, "ppp")
t6.ppp$marks$Width <- as.numeric(as.character((t6.ppp$marks$Width)))
t6.ppp$marks$Height <- as.numeric(as.character((t6.ppp$marks$Height)))

t6num <- t6.ppp
marks.data <- as.data.frame(marks(t6num))
marks.data <- marks.data %>% dplyr::select(Width) 
marks(t6num) <- marks.data


t6species <- t6.ppp
marks.data <- as.data.frame(marks(t6species))
marks.data <- marks.data %>% dplyr::select(sp_simp)
marks(t6species) <- marks.data
plot(pairs(density(split(t6species)[c(2,3,5)])))
panel.histogram(pairs(density(split(t6species)[c(2,3,5)])))

L <- density(split(t6.ppp))
df <- pairs(L, plot=FALSE)
co <- cor(df)
densitycor <- round(co, 2)
write.csv(densitycor, "densitycor.csv")


hist(coords(t6.ppp)$x)
hist(coords(t6.ppp)$y)

QC <- quadratcount(t6.ppp)
I <- intensity.quadratcount(QC, image = TRUE)
plot(I)
plot(split(t6.ppp)[c(4,6,7,8)])
warnings()
plot(intensity(t6species))
unitname(t6)

#pooled analysis of homogeneity
t6.sp <- readOGR("Output/T6_J.shp")
t5.sp <- readOGR("Output/T5_J.shp")
t4.sp <- readOGR("Output/T4_J.shp")
t6.ppp <- as(t6.sp, "ppp")
t5.ppp <- as(t5.sp, "ppp")
t4.ppp <- as(t4.sp, "ppp")

marks.data <- as.data.frame(marks(t6.ppp))
marks.data <- marks.data %>% dplyr::select(sp_simp, Width, Height, Type)
marks(t6.ppp) <- marks.data

marks.data <- as.data.frame(marks(t5.ppp))
marks.data <- marks.data %>% dplyr::select(sp_simp, Width, Height, Type)
marks(t5.ppp) <- marks.data

marks.data <- as.data.frame(marks(t4.ppp))
marks.data <- marks.data %>% dplyr::select(sp_simp, Width, Height, Type)
marks(t4.ppp) <- marks.data

QC1 <- quadrat.test(t6.ppp)
QC2 <- quadrat.test(t5.ppp)
QC3 <- quadrat.test(t4.ppp)
pool(QC1, QC2, QC3)

par(mfrow = c(1,3))
den <- density(t6.ppp, sigma = 1, diggle = TRUE)
plot(den)
den <- density(t5.ppp, sigma = 1)
plot(den)
den <- density(t4.ppp, sigma = 1)
plot(den)

den <- adaptive.density(t6.ppp, f=1)
plot(den)
den <- adaptive.density(t5.ppp, f = 1)
plot(den)
den <- adaptive.density(t4.ppp, f = 1)
plot(den)

plot(t6.ppp, which.marks ="Width")
plot(t5.ppp, which.marks ="Width")
plot(t4.ppp, which.marks ="Width")

smooth <- Smooth.ppp(t6.ppp, which.marks = "Width")
plot(smooth)
smooth <- smooth.ppp(t5.ppp, which.marks = "Width")
plot(smooth)
smooth <- smooth.ppp(t4.ppp, which.marks = "Width")
plot(smooth)

par(mfrow=c(1,1))
miplot(t6.ppp)
miplot(t5.ppp)
miplot(t4.ppp)
plot(frypoints(t6.ppp))

r <- seq(0, 20, by = 0.05)
par(mfrow=c(1,3))
L<- Lest(t6.ppp, r = r)
plot(L, xlim = c(0, 10))
L<- Lest(t5.ppp, r= r)
plot(L, xlim = c(0, 10))
L<- Lest(t4.ppp, r = r)
plot(L, xlim = c(0, 10))

K<- Kest(t6.ppp, r = r)
plot(K, xlim = c(0, 10))
K<- Kest(t5.ppp, r= r)
plot(K, xlim = c(0, 10))
K<- Kest(t4.ppp, r = r)
plot(K, xlim = c(0, 10))

pcf <- pcf(t6.ppp, method = "d")
plot(pcf, xlim = c(0, 10))
pcf <- pcf(t5.ppp, method = "d")
plot(pcf, xlim = c(0, 10))
pcf <- pcf(t4.ppp, method = "d")
plot(pcf, xlim = c(0, 10))

E <- envelope(t6.ppp, Kest, nsim=39, fix.n=TRUE)
plot(E)
E <- envelope(t5.ppp, Kest, nsim=39, fix.n=TRUE)
plot(E)
E <- envelope(t4.ppp, Kest, nsim=39, fix.n=TRUE)
plot(E)

mad.test(t6.ppp, Kest, rmax = 2, nsim = 99)
mad.test(t5.ppp, Kest, nsim = 99)
mad.test(t4.ppp, Kest, nsim = 99)

b <- bw.ppl(t6.ppp)
b
plot(b)
numata <- residualspaper$Fig1
lambda <- density(t6.ppp, bw.ppl)
plot(lambda)
t6K <- Kinhom(t6.ppp, lambda)
plot(t6K)
par(mfrow=c(1,1))
t6L <- Linhom(t6.ppp, lambda)
plot(t6L)

lambda <- density(t5.ppp, bw.ppl)
t5K <- Kinhom(t5.ppp, lambda)
plot(t5K)
par(mfrow=c(1,1))

t5L <- Linhom(t5.ppp, lambda)
plot(t5L)

lambda <- density(t4.ppp, bw.ppl)
t4K <- localKinhom(t4.ppp, lambda)
plot(t4K)

t4L <- Linhom(t4.ppp, lambda)
plot(t4L)


locK <- as.data.frame(t4K)[, fvnames(t4K, ".a")]
rr <- with(t4K, r)
locH <- hclust(dist(t(locK)))
plot(locH)

marks.data <- as.data.frame(marks(t6.ppp))
marks.data <- marks.data %>% dplyr::select(Type)
marks(t6.ppp) <- marks.data
Kc <- Kcross.inhom(t6.ppp, "cactus", "shrub")
plot(Kc)


plot(density(split(t6.ppp)))
plot(split(t6.ppp))

#seedlings
small <- subset(t6num, marks <= 15)
big <- subset(t6num, marks > 15)
a <- plot(t6num)
print(t6num$marks)
plot(small, symap=a)

marks.data <- as.data.frame(marks(t6.ppp))
marks.data <- marks.data %>% dplyr::select(Width) %>% transmute(size = ifelse(Width <=20, "small", "adult")) 
marks.data$size <- as.factor(marks.data$size)
marks(t6.ppp) <- marks.data

lambda <- density(t6.ppp, bw.ppl)
Kc <- Kcross.inhom(t6.ppp, "adult", "small", lambda)
plot(Kc)
pcf(Kc)

G <- Gcross(t6.ppp, "adult", "small", lambda)
plot(G)
print(marks(t6.ppp))

marks.data <- as.data.frame(marks(t5.ppp))
marks.data <- marks.data %>% dplyr::select(Width) %>% transmute(size = ifelse(Width <=20, "small", "adult")) 
marks.data$size <- as.factor(marks.data$size)
marks(t5.ppp) <- marks.data
G <- Gcross(t5.ppp, "adult", "small", lambda)
plot(G)
print(marks(t5.ppp))

marks.data <- as.data.frame(marks(t6.ppp))
marks.data <- marks.data %>% dplyr::select(Width) %>% mutate(size = ifelse(Width <=20, "small", "adult")) 
marks.data <- marks.data %>% mutate(size = ifelse(Width <=20, "small", "adult")) %>% dplyr::select(size, sp_simp)
marks.data <- marks.data %>% dplyr::select(sp_simp)
marks.data$size <- as.factor(marks.data$size)
marks(t6.ppp) <- marks.data

G <- Gcross(t4.ppp, "adult", "small")
plot(G)
print(marks(t4.ppp))

Gal <- Gest(t4.ppp)
plot(Gal)

marks.data <- as.data.frame(marks(t6.ppp))
marks.data <- marks.data %>% dplyr::select(sp_simp)
marks.data$sp_simp <- as.factor(marks.data$sp_simp)
marks.data$size <- as.factor(marks.data$size)
marks(t6.ppp) <- marks.data

b <- bw.relrisk(t6.ppp)
Prob <- relrisk(t6.ppp, sigma=b)
dominant <- im.apply(Prob, which.max)
species <- levels(marks(t6.ppp))
dominant <- eval.im(factor(dominant, levels=1:8, labels=species))
textureplot(dominant)

segregation.test(t6.ppp)

nncorr(t6.ppp)
marktable(t6.ppp, N=1, collapse=TRUE)
library(dixon)
dixon(as.data.frame(t6.ppp))$tablaC

plot(alltypes(t6.ppp, markconnect))

m <- markconnect(t6.ppp, "adult", "small")
plot(m)

m <- markcorr(t6.ppp)
plot(m)
m <- markcrosscorr(t6.ppp)

m1 <- ppm(t6.ppp ~ marks, MultiHard())
is.multitype(t6.ppp)
