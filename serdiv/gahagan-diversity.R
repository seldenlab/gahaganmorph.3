data <- read.csv("gahagan-diversity.csv")
library(vegan)
S <- specnumber(data[2:3])
U <- specnumber(data[2:3])
Up <- U/length(N)*100

H <- diversity(data[2:3])
H
D1 <- diversity(data[2:3], index = "simpson")
D1
D2 <- diversity(data[2:3], index = "invsimpson")
D2
Hmax <- exp(H)
Hmax
J <- H/log(S)
J
