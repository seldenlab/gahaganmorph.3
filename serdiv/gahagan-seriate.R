
data <- read.csv("gahagan-diagnostics.csv")
data$context <- as.character(data$context)
save(data, file="gahagan.RData")
table(data$com)
rowSums(data[, 3:14])
colSums(data[, 3:14])
gt4 <- rowSums(data[, 3:14]) > 3
data.sub <- data[gt4, ]
save(data.sub, file="gahagan.sub.RData")

library(ca)
data.sub.ca <- ca(data.sub[, 3:14])
plot.ca(data.sub.ca, cex=.75)

library(vegan)
data.sub.dec <- decorana(data.sub[, 3:14])
plot(data.sub.dec, display="contexts", type="text")

library(plotrix)
rowlbl <- data.sub$context
idx <- order(data.sub.ca$rowcoord[, 1])
data.sub.pct <- prop.table(as.matrix(data.sub[, 3:14]), 1) * 100
battleship.plot(data.sub.pct[idx, ], cex.labels = .75, maxxspan=.75,
                yaxlab=rowlbl[idx], col="lightblue")

idx2 <- rev(order(data.sub.dec$rproj[, 1]))
battleship.plot(data.sub.pct[idx2, ], cex.labels = .75, maxxspan=.75,
                yaxlab=rowlbl[idx2], col="gray")
plot(idx, idx2, pch=16)
abline(a=0, b=1)
