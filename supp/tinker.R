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

# Assemblage Diversity
## Taxonomic composition

library(here)
library(tidyverse)

# read data
data <- read.csv("gahagan-diversity.csv")

## Alpha diversity
library(vegan)

# assemblage size (N)
N <- rowSums(data[11:25])
N ## assemblage sizes range from:

# how many of each type were found?
T <- colSums(data[11:25])
T

## richness (S) = number of types in assemblage
S <- specnumber(data[11:25])
S

## ubiquity (U) = number of assemblages that contain a particular type
U <- specnumber(data[11:25])
U

## Relative abundance

# mean number of observations for each type
colMeans(data[11:25])

# data by percentage
data.pct <- data[11:25]/N*100
## mean percent (Mp) of each type across assemblage
Mp <- colMeans(data.pct)
Mp

## percentage of sites that have each type
Up <- U/length(N)*100
Up


library(dplyr)
library(reshape)
library(ggplot2)
library(ggpubr)
library(ggExtra)
library(wesanderson)

# relative abundance of juvenile and adult burials
pal <- wes_palette("Moonrise2", 2, type = "continuous")

burials <- data %>% 
  select(context, total_adult, total_juvenile) %>% 
  mutate(context = paste(context, c("(S)", "(S)", "(S)", "(S)",
                                    "(N)", "(N)", "(N)", "(N)"))) %>% 
  melt(id.vars = "context")

# configure plot
plot0 <- burials %>%
  arrange(context) %>% 
  mutate(context = factor(context, levels = c(
    "16CD12-BP1 (N)", "41CE19-F134 (S)", "16RR1-BP3 (S)",
    "16CD12-BP5 (N)", "16RR1-BP2 (S)", "41CE19-F119 (S)", 
    "16CD12-BP2 (N)", "16CD12-BP8 (N)"))) %>% 
  ggplot(aes(x=context, y = value, fill = variable)) + 
  geom_bar(stat="identity", position = "fill") + 
  coord_flip() +
  scale_fill_manual(values = pal) +
  labs(x = "Context",
       y = "Relative abundance (%)",
       fill = "Category")


# relative abundance of burials by category
pal <- wes_palette("Moonrise2", 6, type = "continuous")

burials <- data %>% 
  select(context, adult_male:uid_adult, juvenile_male:uid_juvenile) %>% 
  mutate(context = paste(context, c("(S)", "(S)", "(S)", "(S)",
                                    "(N)", "(N)", "(N)", "(N)"))) %>% 
  melt(id.vars = "context")

# configure plot
plot00 <- burials %>%
  arrange(context) %>% 
  mutate(context = factor(context, levels = c(
    "16CD12-BP1 (N)", "41CE19-F134 (S)", "16RR1-BP3 (S)",
    "16CD12-BP5 (N)", "16RR1-BP2 (S)", "41CE19-F119 (S)", 
    "16CD12-BP2 (N)", "16CD12-BP8 (N)"))) %>% 
  ggplot(aes(x=context, y = value, fill = variable)) + 
  geom_bar(stat="identity", position = "fill") + 
  coord_flip() +
  scale_fill_manual(values = pal) +
  labs(x = "Context",
       y = "Relative abundance (%)",
       fill = "Category")

# relative abundance of lithics + ceramics
pal <- wes_palette("Moonrise2", 2, type = "continuous")

cerlith <- data %>% 
  select(context, lithics, ceramics) %>% 
  mutate(context = paste(context, c("(S)", "(S)", "(S)", "(S)",
                                    "(N)", "(N)", "(N)", "(N)"))) %>% 
  melt(id.vars = "context")

# configure plot
plot1 <- cerlith %>%
  arrange(context) %>% 
  mutate(context = factor(context, levels = c(
    "16CD12-BP1 (N)", "41CE19-F134 (S)", "16RR1-BP3 (S)",
    "16CD12-BP5 (N)", "16RR1-BP2 (S)", "41CE19-F119 (S)", 
    "16CD12-BP2 (N)", "16CD12-BP8 (N)"))) %>% 
  ggplot(aes(x=context, y = value, fill = variable)) + 
  geom_bar(stat="identity", position = "fill") + 
  coord_flip() +
  scale_fill_manual(values = pal) +
  labs(x = "Context",
       y = "Relative abundance (%)",
       fill = "Category")

# relative abundance of diagnostic types
pal <- wes_palette("Moonrise2", 17, type = "continuous")

diagnostics <- data %>% 
  select(context, alba:ceramic_bowl) %>% 
  mutate(context = paste(context, c("(S)", "(S)", "(S)", "(S)",
                                    "(N)", "(N)", "(N)", "(N)"))) %>% 
  melt(id.vars = "context")

# configure plot
plot2 <- diagnostics %>% 
  arrange(context) %>% 
  mutate(context = factor(context, levels = c(
    "16CD12-BP1 (N)", "41CE19-F134 (S)", "16RR1-BP3 (S)",
    "16CD12-BP5 (N)", "16RR1-BP2 (S)", "41CE19-F119 (S)", 
    "16CD12-BP2 (N)", "16CD12-BP8 (N)"))) %>% 
  ggplot( aes(x=context, y = value, fill = variable)) + 
  geom_bar(stat="identity", position = "fill") + 
  coord_flip() +
  scale_fill_manual(values = pal) +
  labs(x = "Context",
       y = "Relative abundance (%)",
       fill = "Types") +
  theme(legend.key.height = unit(0.05, "cm"))

# render figure
figure <- ggarrange(plot0, plot1, plot2,
                    labels = c("a","b", "c"),
                    ncol = 1, nrow = 3)

# plot figure
figure
```