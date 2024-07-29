# Introduction to SCTransform, v2 regularization
# https://satijalab.org/seurat/articles/sctransform_v2_vignette.html

#####################################################################
# install Seurat from Github (automatically updates sctransform)
# invoke sctransform
object <- SCTransform(object, vst.flavor = "v2")

#####################################################################
# install glmGamPoi
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("glmGamPoi")
# install sctransform from Github
devtools::install_github("satijalab/sctransform", ref = "develop", force=TRUE)

#
#####################################################################
# install dataset
InstallData("ifnb")
InstallData('pbmc3k')
#確認
InstalledData()
# load dataset
LoadData("ifnb")

# split the dataset into a list of two seurat objects (stim and CTRL)
ifnb.list <- SplitObject(ifnb, split.by = "stim")

ctrl <- ifnb.list[["CTRL"]]
stim <- ifnb.list[["STIM"]]


# normalize and run dimensionality reduction on control dataset
ctrl <- SCTransform(ctrl, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE) %>%
  RunUMAP(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "pca", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.7, verbose = FALSE)

p1 <- DimPlot(ctrl, label = T, repel = T) + ggtitle("Unsupervised clustering")
p2 <- DimPlot(ctrl, label = T, repel = T, group.by = "seurat_annotations") + ggtitle("Annotated celltypes")

p1 | p2


stim <- SCTransform(stim, vst.flavor = "v2", verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)
ifnb.list <- list(ctrl = ctrl, stim = stim)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)

immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                         anchor.features = features)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT")


immune.combined.sct <- RunPCA(immune.combined.sct, verbose = FALSE)
immune.combined.sct <- RunUMAP(immune.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
immune.combined.sct <- FindNeighbors(immune.combined.sct, reduction = "pca", dims = 1:30)
immune.combined.sct <- FindClusters(immune.combined.sct, resolution = 0.3)



# Identify differential expressed genes across conditions