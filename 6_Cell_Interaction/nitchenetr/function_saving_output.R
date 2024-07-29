library(nichenetr)
library(Seurat)
library(dplyr)
library(tidyverse)
library(openxlsx)

# Preparation:
location="/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/nitchenetr/"
ligand_target_matrix = readRDS(paste0(location,"ligand_target_matrix.rds"))
lr_network = readRDS(paste0(location,"lr_network.rds"))
weighted_networks= readRDS(paste0(location,"weighted_networks.rds"))

## Function


#### Result出力のFunction
## Fun.nichenetr.output.save(SeuratObject, nichenet_output, name, RecieverCells, SplitCondition)
## Code: /Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/nitchenetr/function_saving_output.R
Fun.nichenetr.output.save <- function(SeuratObject, nichenet_output, name, RecieverCells, SplitCondition){
  
  myTHEME <- theme(
    plot.title = element_text(hjust=0.5, size=14, face="bold"),
    plot.subtitle = element_text(size = 14, hjust = 0.5),
    axis.title.x = element_text(hjust=0.5, size= 16, face="bold"),
    axis.title.y = element_text(hjust=0.5, size= 16, face="bold"),
    axis.text.x = element_text(hjust=0.5, colour = "black", size= 14, face="bold"),
    axis.text.y = element_text(hjust=0.5, colour = "black", size= 14, face="bold"),
    #axis.line=element_line(colour = "black"),
    #axis.ticks=element_line(colour = "black"),
    text = element_text(size = 12))#,
  #legend.position="none")
  
  # 保存用directory
  directoryname <- paste0("NicheNet_output_",name)
  dir.create(directoryname)
  
  #1 
  nichenet_output$ligand_expression_dotplot + myTHEME
  ggsave(paste0(directoryname, "/1_top20ligand_DotPlot.png"),
         height=8, width=7, dpi=300)
  #2 
  nichenet_output$ligand_differential_expression_heatmap + 
    xlab("Target genes in Receiver cells") + 
    ylab("Ligand genes in Sender cells") +
    myTHEME
  ggsave(paste0(directoryname,"/2_top20ligand_HeatMap.png"),
         height=8, width=6, dpi=300)
  #3 
  nichenet_output$ligand_target_heatmap + 
    scale_fill_gradient2(low = "whitesmoke",  high = "royalblue") + 
    xlab("Target genes in Receiver cells") + 
    ylab("Ligand genes in Sender cells") +
    myTHEME
      ## グラフの横幅
      if(length(nichenet_output$top_targets)>35){
        w1 = 12 + length(nichenet_output$top_targets)*0.1
        }else{
        w1 = 12
      }
  ggsave(paste0(directoryname,"/3_Ligand_target_Heatmap.png"),
         height=8, width=w1, dpi=300)
  
  
  ##nichenet_output$ligand_target_matrix %>% .[1:10,1:6]
  ##nichenet_output$ligand_target_df # weight column = regulatory potential
  
  #4 dot plot 
  DotPlot(SeuratObject %>% subset(idents = RecieverCells),
          features = nichenet_output$top_targets,
          split.by = SplitCondition) + 
    RotatedAxis() +
    theme(axis.text.x = element_text(colour = "black", size= 14, face="bold"),
          axis.text.y = element_text(colour = "black", size= 14, face="bold")) + 
    xlab("") + 
    ylab("") + 
    labs(title = "Target genes in Receiver cells")
  ggsave(paste0(directoryname,"/4_Top_Targets_in_RecieverCells_Dotplot.png"),
         width=w1, height=4, dpi=300)
  ## 横幅は，#3で出したものを使用．横軸は同じ
  
  #5 Ligand - Receptor activity 
  nichenet_output$ligand_activity_target_heatmap + myTHEME
    ## グラフの横幅
    if(length(nichenet_output$top_targets)>30){
      w1 = 15 + length(nichenet_output$top_targets)*0.1
      }else{
      w1 = 15
      }
    ggsave(paste0(directoryname,"/5_Ligand_Target_Summary.png"),
           width=w1, height=10, dpi=300)
  
  #6 Ligand - Receptor Heatmap   
  nichenet_output$ligand_receptor_heatmap +
    xlab("Receptor genes in Receiver cells") + 
    ylab("Ligand genes in Sender cells") +
    myTHEME
    ## グラフの横幅
    if(length(nichenet_output$top_receptors)>30){
      w1 = 12 + length(nichenet_output$top_receptors)*0.15
    }else{
      w1 = 12
    }
    ggsave(paste0(directoryname,"/6_Ligand_Receptor_Heatmap.png"),
         height=8, width=w1, dpi=300) # widthあえて指定しない
  
  #7 Top receptors in Reciever cells
  DotPlot(SeuratObject %>% subset(idents = RecieverCells),
          features = nichenet_output$top_receptors,
          split.by = SplitCondition) + 
    RotatedAxis() + 
    theme(axis.text.x = element_text(colour = "black", size= 14, face="bold"),
          axis.text.y = element_text(colour = "black", size= 14, face="bold")) + 
    xlab("") + 
    ylab("") + 
    labs(title = "Receptor genes in Receiver cells")
    ggsave(paste0(directoryname,"/7_TopReceptors_by_Condition_DotPlot.png"),
           width=w1, height=4, dpi=300)
    
  ## 横幅は，#3で出したものを使用．横軸は同じ
  #8 BonaFide L-R heatmap
  nichenet_output$ligand_receptor_heatmap_bonafide + myTHEME
    ggsave(paste0(directoryname,"/8_BonaFide_L-R_heatmap.png"),
         height=8, width=8, dpi=300)
  
  
  ## 数値データの保存
  workbook <- createWorkbook()
  # sheetを作る
    addWorksheet(workbook, sheetName = "ligand_target_matrix")
    addWorksheet(workbook, sheetName = "ligand_target_df")
    addWorksheet(workbook, sheetName = "ligand_receptor_matrix")
    addWorksheet(workbook, sheetName = "ligand_receptor_df")
    addWorksheet(workbook, sheetName = "ligand_receptor_matrix_bonafide")
    addWorksheet(workbook, sheetName = "ligand_receptor_df_bonafide")
  # dataを書き込む
    writeData(workbook, sheet = 1, x=nichenet_output$ligand_target_matrix, rowNames = TRUE)
    writeData(workbook, sheet = 2, x=nichenet_output$ligand_target_df, rowNames = FALSE)
    writeData(workbook, sheet = 3, x=nichenet_output$ligand_receptor_matrix, rowNames = TRUE)
    writeData(workbook, sheet = 4, x=nichenet_output$ligand_receptor_df, rowNames = FALSE)
    writeData(workbook, sheet = 5, x=nichenet_output$ligand_receptor_matrix_bonafide, rowNames = TRUE)
    writeData(workbook, sheet = 6, x=nichenet_output$ligand_receptor_df_bonafide, rowNames = FALSE)
  #save
  saveWorkbook(workbook, 
               file = paste0(directoryname,"/9_Ligand_receptor_rawvales.xlsx"),
               overwrite = TRUE)
  
# 結果のobjectも保存
saveRDS(
  nichenet_output,
  paste0(directoryname,"/NicheNetOutput_",name,".rds")
)

}
## Function ここまで

## Output test
## 23/01/02
#Fun.nichenetr.output.save(wh.int, nichenet_PCtoEC, "ECtoPC", "Endothelial", "State")
#Fun.nichenetr.output.save(wh.int, nichenet_ECtoPC, "PCtoEC", "Endothelial", "State")