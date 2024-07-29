## SingleCellSignalR

library(devtools) # githubからdownloadするのに必要．
library(SingleCellSignalR)

# User guide
# https://bioconductor.org/packages/devel/bioc/vignettes/SingleCellSignalR/inst/doc/UsersGuide.html


# Preparation
# The file containing the read counts should be placed in the working directory.

# suerat objectからのとり出し
matrix <- seuratobject[["RNA"]]@counts # clsssはdgCMatrix 
matrix.df <- as.data.frame(as.matrix(matrix))
#  count matrix (csv形式で)


# directoryにcount matrixのファイルをおいておかなくてはいけない
file <- "msenchyme_humangene.txt" #write.csvで保存
data <- data_prepare(file = file) #発現していないgeneは除かれる

clust <- clustering(data = data, 
                    n.cluster = 9, # clusterの数
                    n = 10, # 
                    method = "simlr",
                    write = FALSE,
                    pdf=FALSE)
    # n.cluster: a number, an estimation of the ideal number of clusters is computed if equal to 0
    # n: a number, the maximum to consider for an automatic determination of the ideal number of clusters
    # method: simlr -> caused the SIMLR() function of the SIMLR package 



# cluster analsysis function:
# differentially expressed genes in one cluster compared to the others are identified
clust.ana <- cluster_analysis(data = data,
                              genes = rownames(data),
                              cluster = clust$cluster,
                              write = FALSE)



# The SIMLR_Estimate_Number_of_Clusters() function determined the number of clusters, between 2 and n (n=10 above).


# generate cellular interaction lists using the cell_signaling() function
signal <- cell_signaling(data = data,
                         genes = rownames(data),
                         cluster = clust$cluster,
                         write = FALSE)


# An intercellular network can also be generated to map the overall ligand/receptor interactions
inter.net <- inter_network(data = data,
                           signal = signal,
                           genes = genes,
                           cluster = clust$cluster,
                           write = FALSE)



visualize_interactions(signal = signal)
visualize_interactions(signal = signal, show.in=c(1,4)) # cluster
visualize_interactions(signal = signal, write.in=c(1,4)) #write.inでsave



####  cell_classifier() ########################
##  cell clusters after the output of the cell_classifier()
class = cell_classifier(data=data,
                        genes=rownames(data),
                        markers = markers(c("immune")),
                        tsne=clust$`t-SNE`,
                        plot.details=TRUE,
                        write = FALSE)


# plot.details argument to TRUE to monitor the choice of the threshold of gene signature scores.

# Define the cluster vector and the cluster names 
cluster <- class$cluster
c.names <- class$c.names

# Remove undefined cells 
data <- data[,cluster!=(max(cluster))]
tsne <- clust$`t-SNE`[cluster!=(max(cluster)),]
c.names <- c.names[-max(cluster)]
cluster <- cluster[cluster!=(max(cluster))]





########  Install ###########
devtools::install_github(repo = "https://github.com/SCA-IRCM/SingleCellSignalR_v1", subdir = "SingleCellSignalR")
library(SingleCellSignalR)
