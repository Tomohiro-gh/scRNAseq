library(Seurat)
library(CellChat) # main
library(ggplot2)
library(dplyr)
library(openxlsx)

options(future.globals.maxSize = 1000 * 1024 ^ 2) 

#> Now analyzing
Args = commandArgs(trailingOnly = TRUE) ## 引数受け取り
# Args = commandArgs()
  print(Args)
## 第一引数:ファイルの名前
Args1_Filename = Args[[1]]
  print(paste0 ("Argument1: ", Args1_Filename))
## 第二引数:ファイルの名前
vertex.receiver = Args[[2]]
  print(paste0 ("Argument2: ", vertex.receiver))
  



# 
# 
# ## ------------------------------------------------------------------
# ## Load a script only to my environment
# Dir="/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/DoubletFinder/"
# scripts <- 
#   c("DoubletFinder_function_v2.R")
# 
# .myfunc.env = new.env()
# for(sc in scripts){
#   sys.source(paste0(Dir, sc), envir = .myfunc.env ) # localにload
# }
# attach(.myfunc.env)
# ## -----------------------------------------------------------------
#

## CellChat v2 ---------------------------------------
##  https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat-vignette.html#b-starting-from-a-seurat-object
  ptm = Sys.time()
## Part I: Data input & processing and initialization of CellChat object---------------------------------------
data.input <- subObject[["RNA"]]$data # normalized data matrix
  # For Seurat version >= “5.0.0”, get the normalized data via `seurat_object[["RNA"]]$data`
  labels <- Idents(subObject)
meta <- data.frame(labels = labels, row.names = names(labels)) # create a dataframe of the cell labels


## Create a CellChat object ---------------------------------------
cellchat <- createCellChat(object = subObject, group.by = "ident", assay = "RNA")
# Warning message:
#   In createCellChat(object = subObject, group.by = "ident", assay = "RNA") :
#   The 'meta' data does not have a column named `samples`. We now add this column and all cells are assumed to belong to `sample1`! 

cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
  levels(cellchat@idents) 
  
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group



## Set the ligand-receptor interaction database ---------------------------------------
CellChatDB <- CellChatDB.mouse # use CellChatDB.mouse if running on mouse data
  showDatabaseCategory(CellChatDB)
  # Show the structure of the database
  dplyr::glimpse(CellChatDB$interaction)

# use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, 
                           search = c("Secreted Signaling",
                                      #"Non-protein Signaling",
                                      "ECM-Receptor",
                                      "Cell-Cell Contact"),
                           key = "annotation") # use Secreted Signaling

# Only uses the Secreted Signaling from CellChatDB v1
#  CellChatDB.use <- subsetDB(CellChatDB, search = list(c("Secreted Signaling"), c("CellChatDB v1")), key = c("annotation", "version"))

# use all CellChatDB except for "Non-protein Signaling" for cell-cell communication analysis
# CellChatDB.use <- subsetDB(CellChatDB)

# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB. We do not suggest to use it in this way because CellChatDB v2 includes "Non-protein Signaling" (i.e., metabolic and synaptic signaling). 

# set the used database in the object
  cellchat@DB <- CellChatDB.use


## Preprocessing the expression data for cell-cell communication analysis ---------------------------------------
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
  
future::plan("multisession", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 692
  execution.time = Sys.time() - ptm
    print(as.numeric(execution.time, units = "secs"))
#> [1] 7.404624
# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

    
print("CellChat object has been created")

## ------------------------------------------------------------------------------------------
## Part II: Inference of cell-cell communication network  -----------------------------------
  ptm = Sys.time()
cellchat <- computeCommunProb(cellchat, type = "triMean")

## Users can filter out the cell-cell communication if there are only few cells in certain cell groups. By default, the minimum number of cells required in each cell group for cell-cell communication is 10.
cellchat <- filterCommunication(cellchat, min.cells = 10)

# Extract the inferred cellular communication network as a data frame
df.net <- subsetCommunication(cellchat) # extracitng all daataset
  dim(df.net) #[1] 983  11
  df.net$pathway_name %>% unique
# [1] "TGFb"     "BMP"      "VEGF"     "IGF"      "APELIN"   "CXCL"     "MIF"      "TRAIL"    "VISFATIN"
# [10] "ANGPT"    "EDN"      "KIT"      "SEMA3"    "GAS"      "GALECTIN" "KLK"      "IGFBP"    "NRG"     
# [19] "PROS"    

  write.xlsx(df.net, paste0("0_CellCommunication_", Args1_Filename, ".xlsx"))
  # choice 1
  #  df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
  # from cell groups 1 and 2 to cell groups 4 and 5.
  # choice 2
  # df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
  #  the inferred cell-cell communications mediated by signaling WNT and TGF

##  Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)


## Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

  execution.time = Sys.time() - ptm
    print(as.numeric(execution.time, units = "secs"))


groupSize <- as.numeric(table(cellchat@idents))

png(paste0("0_Number_of_interaction_", Args1_Filename, ".png"),
    width = 2000, height = 2000, res = 300)
  par(xpd = TRUE)
    netVisual_circle(cellchat@net$count, 
                      vertex.weight = groupSize, 
                      weight.scale = T, 
                      label.edge= F, 
                      title.name = "Number of interactions")
    dev.off()

png(paste0("0_Strength&Interaction_", Args1_Filename, ".png"),
    width = 2000, height = 2000, res = 300)
  par(xpd = TRUE)
    netVisual_circle(cellchat@net$weight, 
                    vertex.weight = groupSize, 
                    weight.scale = T, 
                    label.edge= F, 
                    title.name = "Interaction weights/strength")
    dev.off()

## 細胞ごとに分けてinteractionを表示
##  Here we also control the parameter edge.weight.max so that we can compare edge weights between differet networks.
mat <- cellchat@net$weight
  png(paste0("1_eachInteraction_", Args1_Filename, ".png"),
             width = 4000, height = 3000, res = 300)
  par(mfrow = c(3,4), xpd = TRUE)
  for (i in 1:nrow(mat)) {
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[i, ] <- mat[i, ]
        netVisual_circle(mat2, 
                        vertex.weight = groupSize, 
                        weight.scale = T, 
                        edge.weight.max = max(mat), 
                        title.name = rownames(mat)[i])
  }
    dev.off()


## save cellchat object
saveRDS(cellchat, file = paste0("cellchatObject_", Args1_Filename, ".rds"))


## ------------------------------------------------------------------------------------------
## Part III: Visualization of cell-cell communication network -----------------------------------

pathways.show.all <- cellchat@netP$pathways

for(i in 1:length(pathways.show.all)){
  # pathways.show <- c("APELIN") 
  # Hierarchy plot
  # Here we define `vertex.receive` so that the left portion of the hierarchy plot shows signaling to fibroblast and the right portion shows signaling to immune cells 
  
  # vertex.receiver = c(seq(1,4)) # a numeric vector. 
  # netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  # vertex.receiver = c(seq(1,6)) # a numeric vector. 
  # netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver)
  
  pathways.show <- pathways.show.all[i]

  # Circle plot
    png(paste0("PartIII_", pathways.show, "_", Args1_Filename, "_CirclePlot.png"),
        width = 2500, height = 2500, res = 300)
      par(mfrow=c(1,1), xpd = TRUE)
        netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
          dev.off()
  # Chord diagram
    png(paste0("PartIII_", pathways.show, "_", Args1_Filename,  "_ChordPlot.png"),
        width = 2500, height = 2500, res = 300)
      par(mfrow=c(1,1), xpd = TRUE)
        netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
          dev.off()
  # Heatmap
    png(paste0("PartIII_", pathways.show, "_", Args1_Filename, "_Heatmap.png"),
        width = 2000, height = 2000, res = 300)
      par(mfrow=c(1,1), xpd = TRUE)
          netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
            dev.off()
   # # Chord diagram
    # group.cellType <- c(rep("FIB", 4), rep("DC", 4), rep("TC", 4)) # grouping cell clusters into fibroblast, DC and TC cells
    # names(group.cellType) <- levels(cellchat@idents)
    # netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))
#> Plot the aggregated cell-cell communication network at the signaling pathway level
#> 
#> Compute the contribution of each ligand-receptor pair to the overall signaling pathway and visualize cell-cell communication mediated by a single ligand-receptor pair
    # par(mfrow=c(1,1))
    #   png(paste0("_pathway_", pathways.show, "_each_ligand-receptor_pair.png"),
    #         width = 3000, height = 3000, res = 300)
    #   netAnalysis_contribution(cellchat, signaling = pathways.show)
    #       dev.off()
    #   
    
    pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
  
  for(j in 1:length(pairLR)){ 
    
      LR.show <- pairLR[j,] # show one ligand-receptor pair
      # Hierarchy plot
      vertex.receiver = vertex.receiver # a numeric vector
    
        par(mfrow=c(1,1), xpd = TRUE)
          png(paste0("PartIII_", pathways.show, "_", LR.show, "_", Args1_Filename, "_netVisual_individual_circle.png"), 
            width = 3000, height = 3000, res = 300)
        netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "circle")
          dev.off()
        par(mfrow=c(1,1), xpd = TRUE)
          png(paste0("PartIII_", pathways.show, "_", LR.show, "_", Args1_Filename, "_netVisual_individual_chord.png"), 
            width = 3000, height = 3000, res = 300)
        netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = LR.show, layout = "chord")
          dev.off()
    # show all the interactions received by Inflam.DC
    ## netVisual_chord_gene(cellchat, sources.use = seq(1,11), targets.use = 8, legend.pos.x = 15)
          
  } # pairLRのloop
   
  #  netVisual_individual(cellchat, signaling = pathways.show,  pairLR.use = LR.show, vertex.receiver = vertex.receiver)
  #> [[1]]
  # Circle plot
  
    ## あえて上に持ってくる
    # Plot the signaling gene expression distribution using violin/dot plot
    par(mfrow=c(1,1), xpd = TRUE)
    png(paste0("PartIII_", pathways.show, "_", Args1_Filename, "_Vlnplot.png"), 
        width = 3000, height = 3000, res = 300)
    plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = TRUE, type = "violin")
    dev.off()
    
    plotGeneExpression(cellchat, signaling = pathways.show, enriched.only = TRUE, type = "violin")
    ggsave(paste0("PartIII_", pathways.show, "_", Args1_Filename, "_Vlnplot.png"), 
           width = 8, height = 8, dpi = 300)
    
    ## あえて上に持ってくる

  ## Automatically save the plots of the all inferred network for quick exploration
  # Access all the signaling pathways showing significant communications
  # pathways.show.all <- cellchat@netP$pathways
  # check the order of cell identity to set suitable vertex.receiver
  # levels(cellchat@idents)
  vertex.receiver = vertex.receiver
    # Visualize communication network associated with both signaling pathway and individual L-R pairs
    netVisual(cellchat, 
              signaling = pathways.show,
              vertex.receiver = vertex.receiver,
              layout = "hierarchy",
              height = 5)
    # Compute and visualize the contribution of each ligand-receptor pair to the overall signaling pathway
    gg <- netAnalysis_contribution(cellchat, signaling = pathways.show)
      ggsave(filename = paste0("PartIII_", pathways.show, "_", Args1_Filename, "_L-R_contribution.png"),
              plot = gg, 
              width = 6, height = 3, dpi = 300)

} # end pathway.show.all

## ------------
## Visualize cell-cell communication mediated by multiple ligand-receptors or signaling pathways

