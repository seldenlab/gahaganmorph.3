# load libraries
library(here)
library(geomorph)
library(tidyr)
library(ggplot2)
library(wesanderson)

source('readmulti.csv.R')

# read .csv files
setwd("./data")
filelist <- list.files(pattern = ".csv")
coords <- readmulti.csv(filelist)
setwd("../")

# read qualitative data
qdata <- read.csv("qdata.csv", header = TRUE, row.names = 1)
qdata <- qdata[match(dimnames(coords)[[3]], rownames(qdata)),]

Y.gpa <- gpagen(coords,
                PrinAxes = TRUE,
                print.progress = FALSE)
## plot gpa
plot(Y.gpa)

# geomorph data frame
gdf <- geomorph.data.frame(shape = Y.gpa$coords,
                           size = Y.gpa$Csize,
                           context = qdata$context)  

# add centroid size to qdata
qdata$csz <- Y.gpa$Csize

# attributes for boxplot
csz <- qdata$csz
merged <- qdata$merged

# palette
pal <- wes_palette("Moonrise2")

# boxplot - centroid size by context
csz.temp <- ggplot(qdata, aes(x = merged, y = csz, color = merged)) +
  geom_boxplot() +
  geom_dotplot(binaxis = 'y', stackdir = 'center', dotsize = 0.3) +
  scale_color_manual(values = pal) +
  theme(legend.position = "none") +
  labs(x = 'Context', y = 'Centroid Size')

## render plot
csz.temp

# principal components analysis
pca <- gm.prcomp(Y.gpa$coords)
summary(pca)

# palette
pal <- wes_palette("Moonrise2")

# principal components analysis
pca <- gm.prcomp(Y.gpa$coords)
summary(pca)

# set plot parameters to plot by context
pch.gps.context <- c(15,17)[as.factor(qdata$context)]
col.gps.context <- pal[as.factor(qdata$context)]
col.hull.context <- c("#C27D38","#798E87")

## plot pca by context 2
pc.plot <- plot(pca, asp = 1,
                pch = pch.gps.context,
                col = col.gps.context)
shapeHulls(pc.plot,
           groups = qdata$context,
           group.cols = col.hull.context)

## Procrustes ANOVA
# MODEL: shape as a function of context + time
fit.shapecontemp <- procD.lm(shape ~ merged,
                             data = gdf,
                             print.progress = FALSE,
                             iter = 9999)

# ANOVA: do gahagan biface shapes differ by context + time?
anova(fit.shapecontemp)

# pairwise comparison of LS means = which differ?
pairwise.shapecontemp <- pairwise(fit.shapecontemp,
                                  groups = qdata$merged)
summary(pairwise.shapecontemp, 
        confidence = 0.95, 
        test.type = "dist")

# MODEL: size as a function of context + time
fit.sizecontemp <- procD.lm(size ~ merged,
                            data = gdf,
                            print.progress = FALSE,
                            iter = 9999)

# ANOVA: do gahagan biface sizes differ by context + time?
anova(fit.sizecontemp)

# pairwise comparison of LS means = which differ?
pairwise.sizecontemp <- pairwise(fit.sizecontemp,
                                 groups = qdata$merged)
summary(pairwise.sizecontemp, 
        confidence = 0.95, 
        test.type = "dist")

# subset landmark coordinates to produce mean shapes for contexts
new.coords <- coords.subset(A = Y.gpa$coords,
                            group = qdata$merged)

names(new.coords)

## plot shape means
mean <- lapply(new.coords, mshape)
plot(mean$initial_cache)
plot(mean$initial_individual)
plot(mean$subsequent_cache)
plot(mean$subsequent_individual)

# comparison plots
plotRefToTarget(mean$initial_individual,
                mean$initial_cache, 
                method = "points",
                mag = 1)

plotRefToTarget(mean$initial_individual,
                mean$subsequent_individual, 
                method = "points",
                mag = 1)

plotRefToTarget(mean$initial_individual,
                mean$subsequent_cache, 
                method = "points",
                mag = 1)

plotRefToTarget(mean$initial_cache,
                mean$subsequent_individual, 
                method = "points",
                mag = 1)

plotRefToTarget(mean$subsequent_individual,
                mean$initial_cache, 
                method = "points",
                mag = 1)

plotRefToTarget(mean$initial_cache,
                mean$subsequent_cache, 
                method = "points",
                mag = 1)
