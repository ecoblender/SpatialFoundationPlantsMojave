library(dplyr)
library(tidyr)


all.data <- read.csv("Data/shrub.data.csv")
all.data <- separate(all.data, Quadrat, c("QN", "QL"), sep = 1)
#extract each transect and map them 
t6 <- filter(all.data, Transect == "6")
t5 <- filter(all.data, Transect == "5")
t4 <- filter(all.data, Transect == "4")

#t6 first quadrat 2 is 11.6m long
# these are coordinates for first transect 622,565.377  3,849,067.992


mapshrubs <- function(data, x, y, weird){
  mutate(data, easting = ifelse(QL == "R", (x-10)+X/100, ifelse(QL == "L", x+X/100, NA))) %>%              
  mutate(., northing = ifelse(QN == "8", (y-10)+Y/100,+
                                  ifelse(QN == "7", y+Y/100,+
                                  ifelse(QN == "6", (y+10) + Y/100, 
                                  ifelse(QN == "5",(y+20) + Y/100, +
                                  ifelse(QN == "4", (y+30) + Y/100, +
                                  ifelse(QN == "3", (y+40) + Y/100, +
                                  ifelse(QN == "2", (y+40+weird) + Y/100, +
                                  ifelse(QN == "1", (y+50 + weird) +Y/100, NA)))))))))
}

t6 <- mapshrubs(t6, 622565.377, 3849067.992, 11.6)                
t5 <- mapshrubs(t5, 622549.553 , 3849146.702, 10.9)    
t4 <- mapshrubs(t4, 622480.394, 3849254.652, 11.05)
                                  
mapped.data <- rbind(t6, t5, t4)
mapped.data$Width <- as.numeric(mapped.data$Width)
mapped.data <- na.omit(mapped.data)
write.csv(mapped.data, "Data/testmap.csv")

write.csv(t6, "Data/t6.csv")
write.csv(t5, "Data/t5.csv")
write.csv(t4, "Data/t4.csv")

str(mapped.data)
#needed to use GIS for T6

 





