# Stacked ViolinPlot
## https://github.com/ycl6/StackedVlnPlot

library(Seurat)
library(ggplot2)
library(cowplot)

# Load Seurat obj
pbmc <- readRDS("data/pbmc_2k_v3_Seurat.rds")

features <- c("CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "LGALS3", "S100A8", "GNLY",
              "NKG7", "KLRB1", "FCGR3A", "FCER1A", "CST3")

a <- VlnPlot(pbmc, features, stack = TRUE, sort = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on y-axis")

b <- VlnPlot(pbmc, features, stack = TRUE, sort = TRUE, flip = TRUE) +
  theme(legend.position = "none") + ggtitle("Identity on x-axis")


## With ggplot2 and a data.frame object
# Load data.frame obj
pbmc <- readRDS("data/pbmc_2k_v3_df.rds")
identity <- readRDS("data/pbmc_2k_v3_Seurat_Idents.rds")

features <- c("CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "LGALS3", "S100A8", "GNLY",
              "NKG7", "KLRB1", "FCGR3A", "FCER1A", "CST3")

# Subset data.frame
pbmc <- pbmc[,features]

# Add cell ID and identity classes
pbmc$Cell <- rownames(pbmc)
pbmc$Idents <- identity

# Use melt to change data.frame format
pbmc <- reshape2::melt(pbmc, id.vars = c("Cell","Idents"), measure.vars = features,
                       variable.name = "Feat", value.name = "Expr")
head(pbmc, 10)

# Identity on x-axis
a <- ggplot(pbmc, aes(factor(Idents), Expr, fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0)) +
  ggtitle("Identity on x-axis") + xlab("Identity") + ylab("Expression Level")

# Identity on y-axis
b <- ggplot(pbmc, aes(Expr, factor(Idents), fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(cols = vars(Feat), scales = "free")  +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  ggtitle("Identity on y-axis") + xlab("Expression Level") + ylab("Identity")




## With ggplot2 and a data.frame object
# Load data.frame obj
pbmc <- readRDS("data/pbmc_2k_v3_df.rds")
identity <- readRDS("data/pbmc_2k_v3_Seurat_Idents.rds")

features <- c("CD79A", "MS4A1", "CD8A", "CD8B", "LYZ", "LGALS3", "S100A8", "GNLY",
              "NKG7", "KLRB1", "FCGR3A", "FCER1A", "CST3")

# Subset data.frame
pbmc <- pbmc[,features]

# Add cell ID and identity classes
pbmc$Cell <- rownames(pbmc)
pbmc$Idents <- identity

# Use melt to change data.frame format
pbmc <- reshape2::melt(pbmc, id.vars = c("Cell","Idents"), measure.vars = features,
                       variable.name = "Feat", value.name = "Expr")
head(pbmc, 10)



# Identity on x-axis
a <- ggplot(pbmc, aes(factor(Idents), Expr, fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Feat), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0)) +
  ggtitle("Identity on x-axis") + xlab("Identity") + ylab("Expression Level")

# Identity on y-axis
b <- ggplot(pbmc, aes(Expr, factor(Idents), fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(cols = vars(Feat), scales = "free")  +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5)) +
  ggtitle("Identity on y-axis") + xlab("Expression Level") + ylab("Identity")

# Use plot_grid to join plots
plot_grid(a, b, labels = c("A","B"))



## Features on x-axis (C) or y-axis (D)
# Features on x-axis
c <- ggplot(pbmc, aes(factor(Feat), Expr, fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Idents), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Feature on x-axis") + xlab("Feature") + ylab("Expression Level")

# Features on y-axis
d <- ggplot(pbmc, aes(Expr, factor(Feat), fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_x_continuous(expand = c(0, 0), labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(cols = vars(Idents), scales = "free") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold")) +
  ggtitle("Feature on y-axis") + xlab("Expression Level") + ylab("Feature")





## Sort identity classes and features
## Note: Some of the codes below are taken and modified from the Seurat package.
## Below demonstrates how to recreate the reordering of the identity classes and features seen in Seurat’s stacked violin plots.

# Calculate average expression per Idents, output as wide format
avg <- sapply(X = split(x = pbmc, f = pbmc$Idents),
              FUN = function(df) { return(tapply(X = df$Expr, INDEX = df$Feat, FUN = mean)) })

# L2Norm (Euclidean norm) function
L2Norm <- function(mat, MARGIN){
  normalized <- sweep(x = mat, MARGIN = MARGIN,
                      STATS = apply(X = mat, MARGIN = MARGIN,
                                    FUN = function(x){ sqrt(x = sum(x ^ 2)) }), FUN = "/")
  normalized[!is.finite(x = normalized)] <- 0
  return(normalized)
}

# Performs hierarchical clustering
idents.order <- hclust(d = dist(t(L2Norm(mat = avg, MARGIN = 2))))$order
avg <- avg[,idents.order]
avg <- L2Norm(mat = avg, MARGIN = 1)
mat <- hclust(d = dist(avg))$merge

# Order feature clusters by position of their "rank-1 idents"
position <- apply(X = avg, MARGIN = 1, FUN = which.max)
orderings <- list()
for (i in 1:nrow(mat)) {
  x <- if (mat[i,1] < 0) -mat[i,1] else orderings[[mat[i,1]]]
  y <- if (mat[i,2] < 0) -mat[i,2] else orderings[[mat[i,2]]]
  x.pos <- min(x = position[x])
  y.pos <- min(x = position[y])
  orderings[[i]] <- if (x.pos < y.pos) { c(x, y) } else { c(y, x) }
}
features.order <- orderings[[length(orderings)]]

# Update Feature and Identity factor orders
pbmc$Idents <- factor(pbmc$Idents, levels = levels(pbmc$Idents)[idents.order])
pbmc$Feat <- factor(pbmc$Feat, levels = levels(pbmc$Feat)[features.order])

# Plot stacked violin plot with reordered identity classes and features
e <- ggplot(pbmc, aes(factor(Feat), Expr, fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Idents), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  ggtitle("Identity & feature ordered") + xlab("Feature") + ylab("Expression Level")
e



## Add gene grouping annotation
## Below demonstrates how to add gene grouping annotation to sorted stacked violin plots.
# Create grouping info
df <- data.frame(x = levels(pbmc$Feat), group = c("A","A","B","B","B","B","B","C","C","C","D","D","D"), 
                 stringsAsFactors = FALSE)
df$x <- factor(df$x, levels = levels(pbmc$Feat))
df$group <- factor(df$group)
df


color <- c("cyan", "pink", "green", "darkorange")

# Same as plot e, but hide x-axis labels, change plot.margin to reduce spacing between plots
f <- ggplot(pbmc, aes(Feat, Expr, fill = Feat)) +
  geom_violin(scale = "width", adjust = 1, trim = TRUE) +
  scale_y_continuous(expand = c(0, 0), position="right", labels = function(x)
    c(rep(x = "", times = length(x)-2), x[length(x) - 1], "")) +
  facet_grid(rows = vars(Idents), scales = "free", switch = "y") +
  theme_cowplot(font_size = 12) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        plot.title = element_text(hjust = 0.5),
        panel.background = element_rect(fill = NA, color = "black"),
        plot.margin = margin(7, 7, 0, 7, "pt"),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold"),
        strip.text.y.left = element_text(angle = 0),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank()) +
  ggtitle("Feature on x-axis with annotation") + ylab("Expression Level")

# Use geom_tile() to add grouping colorings and geom_text() to add grouping labels
g <- ggplot(df, aes(x = x, y = 1, fill = group, label = group)) + geom_tile() +
  geom_text(fontface = "bold", size = 3) + theme_bw(base_size = 12) +
  scale_fill_manual(values = color) + scale_y_continuous(expand = c(0, 0)) +
  theme(legend.position = "none", panel.spacing = unit(0, "lines"),
        panel.background = element_blank(), 
        panel.border = element_blank(),
        plot.background = element_blank(), 
        plot.margin = margin(0, 7, 7, 7, "pt"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) + xlab("Feature")

# Use plot_grid to join plots
plot_grid(f, g, ncol = 1, rel_heights = c(0.78, 0.22), align = "v", axis = "lr")

## Legend is used to defind the grouping labels when the labels are too long to fit within the annotation bar.
# Change to long names
levels(df$group) = c("long long name A", "long long name B", "long long name C", "long long name D")

# guides() is used to specify some aesthetic parameters of legend key
h <- ggplot(df, aes(x = x, y = 1, fill = group)) + geom_tile() + theme_bw(base_size = 12) +
  scale_fill_manual(values = color) + scale_y_continuous(expand = c(0, 0)) +
  guides(fill = guide_legend(direction = "vertical", label.position = "right",
                             title.theme = element_blank(), keyheight = 0.5, nrow = 2)) +
  theme(legend.position = "bottom",
        legend.justification = "left",
        legend.margin = margin(0,0,0,0),
        legend.box.margin = margin(-10,5,0,0),
        panel.spacing = unit(0, "lines"),
        panel.background = element_blank(),
        panel.border = element_blank(),
        plot.background = element_blank(),
        plot.margin = margin(0, 7, 7, 7, "pt"),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
        axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.y = element_blank()) + xlab("Feature")

# Use plot_grid to join plots
plot_grid(f, h, ncol = 1, rel_heights = c(0.78, 0.22), align = "v", axis = "lr")


