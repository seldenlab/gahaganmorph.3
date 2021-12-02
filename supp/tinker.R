# Mortuary Assemblage Diversity

data <- read.csv("gahagan-diagnostics.csv")

library(vegan)
# assemblage size (N)
N <- rowSums(data[3:14])
N ## assemblage sizes range from 8 - 184

# how many of each type were found?
T <- colSums(data[3:14])
T ## quantity of each type (T) found across contexts

# data by percentage
data.pct <- data[3:14]/N*100
# mean percent (Mp) of each type across assemblage
Mp <- colMeans(data.pct)
Mp

#diversity, richness
# richness (S) = number of types in assemblage
S <- specnumber(data[3:14])
S

# ubiquity (U) = number of assemblages that contain a particular type
U <- specnumber(data[3:14])
U

## percentage of sites that have each type
Up <- U/length(N)*100
Up

# Shannon diversity
H <- diversity(data[3:14])
## high diversity = more types and spread more evenly over types
H

# Simpson index
D1 <- diversity(data[3:14], 
                index = "simpson")
## probability that two artifacts drawn randomly represent different types
D1

# inverse Simpson index
D2 <- diversity(data[3:14], 
                index = "invsimpson")
## effective number of types
D2

# effective number of types for Shannon diversity index
Hmax <- exp(H)
Hmax

# evenness
# Pielou's J
J <- H/log(S)
## Shannon diversity index divided by natural log of richness
J

# ratio of effective species to richness
E <- Hmax/S
E

# summarize assemblage diversity to identify high & low diversity assemblages
library(maptools)
pch <- c(1, 3)[as.factor(data$region)]
plot(H, J, 
     pch = pch)
abline(h = median(J), 
       v = median(H), 
       lty = 2)
pointLabel(H, J, 
           rownames(data), 
           cex = .75)
leg.txt <- c(as.expression(bquote("Northern Behavioral Region")), 
             as.expression(bquote("Southern Behavioral Region")))
legend("bottomright", 
       leg.txt, 
       pch = c(1, 3))

# sample size and richness
plot(S~N, 
     ylim = c(0,12), 
     xlim = c(0, 184), 
     pch = 16)
abline(h = seq(0, 12, by = 2), 
       v = seq(0, 184, by = 25), 
       col = "black", 
       lty = 3)

# logarithmic function
data.log <- lm(S~log(N))
summary(data.log)
deviance(data.log)
xval <- seq(1, 184, by = 1)
lines(xval, 
      predict(data.log,
              data.frame(N = xval)),
      lty= 1)

# power function
data.pow <- lm(log(S)~log(N))
summary(data.pow)
sum((S-exp(fitted(data.pow)))^2)
lines(xval, 
      exp(predict(data.pow,
                  data.frame(N = xval))), 
      lty = 2)

# rarefaction curve
## how do individual assemblages compare to the composite?
xval <- seq(2, 200, by = 2)
data.rar <- rarefy(T, xval, se = TRUE)
Est <- data.rar[1, ]
Sd <- data.rar[2, ]
rare <- cbind(lower = Est-2*Sd, 
              expected = Est, 
              upper = Est+2*Sd)

plot(S~N, 
     ylim = range(rare), 
     xlim = range(xval), 
     pch = 16)
matlines(xval, rare, 
         type = "l", 
         lty = c(2, 1, 2), 
         col = "black")
identify(N, S, rownames(data))

## how different are the smaller samples from the largest one?
data.rare <- data[1,3:14]
data.184 <- rarefy(data.rare[1,], xval[1:8], se = TRUE)
Est <- data.184[1, ]
Sd <- data.184[2, ]
rare184 <- cbind(lower = Est-2*Sd,
                 expected = Est,
                 upper = Est+2*Sd)

plot(S~N,
     ylim = range(rare184),
     xlim = range(xval[1:8]),
     pch = 16)

