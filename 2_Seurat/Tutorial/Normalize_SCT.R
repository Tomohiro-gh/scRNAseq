# Integration with SCTransform
wh.list.sct <- wh.list

wh.list.sct <- lapply(X = wh.list.sct,
                      FUN = SCTransform,
                      assay = 'RNA',
                      new.assay.name = 'SCT',
                      method = "glmGamPoi",
                      vars.to.regress = c("percent.mt"),
                      vst.flavor = "v2",
                      return.only.var.genes = TRUE) #default T

features <- SelectIntegrationFeatures(object.list = wh.list.sct, 
                                      nfeatures = 3000)
wh.list.sct <- PrepSCTIntegration(object.list = wh.list.sct, 
                                  anchor.features = features)

wh.list.sct <- lapply(X = wh.list.sct, FUN = RunPCA, features = features)

# perform integration (PCAやらないとエラーが出る)
# 引数にreduction="rpca"を入れるなら，上でRunPCAをかけなくてはならない
anchors2 <- FindIntegrationAnchors(object.list = wh.list.sct,
                                   normalization.method = "SCT", 
                                   anchor.features = features, 
                                   dims = 1:30, 
                                   reduction = "rpca")
wh.int.sct <- IntegrateData(anchorset = anchors2,
                            normalization.method = "SCT", 
                            features.to.integrate = features,
                            verbose = TRUE)


# Confirmation
wh.int.sct
# An object of class Seurat 
# 48086 features across 30629 samples within 3 assays 
# Active assay: integrated (3000 features, 3000 variable features)
# 2 other assays present: RNA, SCT


nrow(wh.int.sct[["integrated"]]@scale.data) #3000
length(features) #3000
DefaultAssay(wh.int.sct) <- "SCT"

nrow(wh.int.sct[["SCT"]]@scale.data) #3000
nrow(wh.int.sct[["SCT"]]@counts) #19103
nrow(wh.int.sct[["SCT"]]@data) #19103













