# SCTransformation


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# https://github.com/satijalab/seurat/issues/4085

LoadData("ifnb")

ifnb.list <- SplitObject(ifnb, split.by = "stim")
ifnb.list <- lapply(X = ifnb.list, FUN = SCTransform, method = "glmGamPoi")
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
ifnb.list <- lapply(X = ifnb.list, FUN = RunPCA, features = features)


immune.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT", 
                                         anchor.features = features, dims = 1:30, reduction = "rpca", k.anchor = 20)
all.feature <- rownames(ifnb[["RNA"]]@counts)
immune.combined.sct <- IntegrateData(anchorset = immune.anchors, normalization.method = "SCT", dims = 1:30, features.to.integrate = all.feature)
nrow(immune.combined.sct[["integrated"]]@scale.data)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


# >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
