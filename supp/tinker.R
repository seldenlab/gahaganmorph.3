# Seriation

data <- read.csv("gahagan-diagnostics.csv")
data$context <- as.character(data$context)
rowSums(data[, 3:15])
colSums(data[, 3:15])

library(ca)
data.ca <- ca(data[, 3:15])
plot.ca(data.ca, cex=.75)

library(vegan)
data.dec <- decorana(data[, 3:15])
plot(data.dec, display="both")

library(plotrix)
rowlbl <- data$context
idx <- rev(order(data.ca$rowcoord[, 1]))
data.pct <- prop.table(as.matrix(data[, 3:15]), 1) * 100
battleship.plot(data.pct[idx, ], 
                cex.labels = .75, 
                maxxspan=.75,
                yaxlab=rowlbl[idx], 
                col="gray")

idx2 <- rev(order(data.dec$rproj[, 1]))
battleship.plot(data.pct[idx2, ], 
                cex.labels = .75, 
                maxxspan=.75,
                yaxlab=rowlbl[idx2], 
                col="gray")

plot(idx, idx2, pch=16)
abline(a=0, b=1)

# Mortuary Assemblage Diversity

data <- read.csv("gahagan-diagnostics.csv")

library(vegan)
# assemblage size (N)
N <- rowSums(data[3:15])
N ## assemblage sizes range from 8 - 184

# how many of each type were found?
T <- colSums(data[3:15])
T ## quantity of each type (T) found across contexts

# data by percentage
data.pct <- data[3:15]/N*100
# mean percent (Mp) of each type across assemblage
Mp <- colMeans(data.pct)
Mp

#diversity, richness
# richness (S) = number of types in assemblage
S <- specnumber(data[3:15])
S

# ubiquity (U) = number of assemblages that contain a particular type
U <- specnumber(data[3:15])
U

## percentage of sites that have each type
Up <- U/length(N)*100
Up

# Shannon diversity
H <- diversity(data[3:15])
## high diversity = more types and spread more evenly over types
H

# Simpson index
D1 <- diversity(data[3:15], 
                index = "simpson")
## probability that two artifacts drawn randomly represent different types
D1

# inverse Simpson index
D2 <- diversity(data[3:15], 
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

# nonlinear least squares
data.nls <- nls(S ~ SSasymp(N, Asym, R0, lrc))
summary(data.nls)
deviance(data.nls)
lines(xval, predict(data.nls, data.frame(N = xval)), lty = 3, lwd = 2)
legend("bottomright",
       c("Logarithmic", "Power (log-log)", "Asymptotic Curve"),
       lty = 1:3,
       lwd = c(1, 1, 2),
       bg = "white")

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
data.184 <- rarefy(data.rare[1,], xval[1:7], se = TRUE)
Est <- data.184[1, ]
Sd <- data.184[2, ]
rare184 <- cbind(lower = Est-2*Sd,
                 expected = Est,
                 upper = Est+2*Sd)

plot(S~N,
     ylim = range(rare184),
     xlim = range(xval[1:7]),
     pch = 16)
identify(N, S, rownames(data))




```{r h1a, out.width = "100%", dpi = 300, echo=TRUE, warning=FALSE, fig.cap="Stratigraphic position of Burial Pits 1 and 2 at the Mounds Plantation site; adapted from Webb (1975:Figure 7). Burial Pit 2 is the only burial found to be intrusive from the mound surface."}
knitr::include_graphics('images/h1a.jpg')
```

```{r h1a.box, out.width = "100%", dpi = 300, echo=TRUE, warning=FALSE, fig.cap="Difference in centroid size between Gahagan bifaces from earlier (BPs 1, 5, and 8) and later (BP 2) contexts at the Mounds Plantation site."}
knitr::include_graphics('images/mp-size.png')
```

#### _Hypothesis 1b: Temporal change in Caddo preference at George C. Davis_

To assess the temporal change in preference between caches of Gahagan bifaces recovered from Mound C at the George C. Davis site, those from Feature 134 are contrasted with those from Feature 119 [@RN3682;@RN8186]. The stratigraphic position of Feature 119 indicates that the burials in that feature occurred subsequent to those associated with Feature 134 [@RN3682;@RN8186].

```{r h1b, out.width = "100%", dpi = 300, echo=TRUE, warning=FALSE, fig.cap="Stratigraphic position of Features 119 and 134 at the George C. Davis site; adapted from Story (1997:Figure 13)."}
knitr::include_graphics('images/h1b.jpg')
```

##### _Hypothesis 1b Findings_:

- Gahagan bifaces from F134 and F119 do not differ in _shape_
- Gahagan bifaces from F134 and F119 do not differ in _size_

```{r h1b.box, out.width = "100%", dpi = 300, echo=TRUE, warning=FALSE, fig.cap="Difference in centroid size between Gahagan bifaces from earlier (F134) and later (F119) contexts at the George C. Davis site."}
knitr::include_graphics('images/gcd-size.png')
```

#### Hypothesis 1: Summary of results

Hypothesis 1 was tested using two samples of Gahagan bifaces from Caddo burial contexts at the Mounds Plantation and George C. Davis sites where stratigraphy dictates differing (earlier/later) temporal positions. There is an important cultural distinction between archaeological contexts at Mounds Plantation and George C. Davis, in that all of the Gahagan bifaces from Mounds Plantation were recovered _in association with an individual_ (Burials 1, 2, 5, and 8) [@RN8174], and all but two of the Gahagan bifaces from George C. Davis articulate with _caches included along the northern periphery of two group burials_ (Features 119 and 134) [@RN3682;@RN8186].

Previous studies have demonstrated a significant difference in shape between Gahagan bifaces recovered from the Mounds Plantation and George C. Davis sites [@RN11783]. In burials at Mounds Plantation, the Caddo were `selecting` for Gahagan bifaces that were significantly smaller in later contexts; a pattern that is not present at George C. Davis.

The results raise substantive questions regarding the cultural or social practices that were driving the temporal shift in size at Mounds Plantation. Might this size difference articulate with a shift in trading or exchange-based relationships with central Texas groups, and/or might the shift be related to a shift in functional use?
  
  At both sites, shape remains consistent and does not differ between contexts, indicating an established `shape preference` that may have shifted at Mounds Plantation due to cyclical differences in the variable social mechanisms associated with raw material procurement. It may also be the case that the size difference was related to signaling, as well as the intended audience. For instance, larger bifaces selected for inclusion with caches may suggest that the intended audience was the group; representative of a potentially ostentatious display. The smaller Gahagan bifaces from individual burials may have been personal items belonging to the deceased, intended to signal status among personal contacts.
