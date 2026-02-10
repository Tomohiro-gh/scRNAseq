library(clusterProfiler)
library(ReactomePA)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(org.Dr.eg.db)

#> 
#> version1 : 25/02/05
#> 
#> functionにしてみる
Fun.Comp.Clsuter.GOKGRA.v1 <- 
  
  function(NamedGenelists, filename, Wid = 18, Hei = 6, nCat = N){
  FileNameBody = paste0(filename)
  
  ## Theme -----------------------------------
  theme_cp <- theme(
    axis.text.y = element_text(size = 12, lineheight = 0.6),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.50, size = 12))
  ## Theme -----------------------------------

    ## GO:BP
    comp.BP <- compareCluster(geneClusters = NamedGenelists, 
                              fun = "enrichGO", 
                              ont="BP",
                              OrgDb = "org.Mm.eg.db")
      if(as.data.frame(comp.BP) |> nrow() >= 1){
        comp.BP <-
          comp.BP |>
          setReadable(OrgDb = org.Mm.eg.db, keyType="ENTREZID")
        
        comp.BP |>
          clusterProfiler::dotplot(showCategory = nCat) + 
            labs(title = paste0("CompareCluster ", FileNameBody, " GO:BP"), x = "") + 
            theme_cp +
            coord_flip()
              ggsave(paste0(FileNameBody, "_Compare.GOBP.png"),
                     width = Wid, height = Hei, dpi = 200)
      }
    
    ## GO:CC
    comp.CC <- compareCluster(geneClusters = NamedGenelists, 
                                       fun = "enrichGO", 
                                       ont="CC",
                                       OrgDb = "org.Mm.eg.db")
      if(as.data.frame(comp.CC) |>  nrow() >= 1){
        comp.CC  <- 
          comp.CC |>
          setReadable(OrgDb = org.Mm.eg.db, keyType="ENTREZID")
        
        comp.CC |>
          clusterProfiler::dotplot(showCategory = nCat) + 
              labs(title = paste0("CompareCluster ", FileNameBody, " GO:CC"), x = "") + 
              theme_cp +
              coord_flip()
                ggsave(paste0(FileNameBody, "_Compare.GOCC.png"),
                          width = Wid, height = Hei, dpi = 200)
      }
      
    ## GO:MF 
    comp.MF <- compareCluster(geneClusters = NamedGenelists, 
                                       fun = "enrichGO", 
                                       ont="MF",
                                       OrgDb = "org.Mm.eg.db")
      if(as.data.frame(comp.MF) |>  nrow() >= 1){
        comp.MF <- 
          comp.MF |>
          setReadable(OrgDb = org.Mm.eg.db, keyType="ENTREZID")
      
        p3 <- 
          clusterProfiler::dotplot(comp.MF, showCategory = nCat) + 
            labs(title = paste0("CompareCluster ", FileNameBody, " GO:MF"), x = "") + 
            theme_cp +
            coord_flip()
              plot(p3)
                ggsave(paste0(FileNameBody, "_Compare.GOMF.png"),
                        width = Wid, height = Hei, dpi = 200)
        }
    
    ## KEGG
    comp.KG <- compareCluster(geneClusters = NamedGenelists, 
                                         fun = "enrichKEGG", 
                                         organism = "mmu") 

    if(as.data.frame(comp.KG) |> nrow() >= 1){
      comp.KG <- 
        comp.KG |>
        setReadable(OrgDb = org.Mm.eg.db, keyType = "ENTREZID") 
      
      comp.KG@compareClusterResult <-
        comp.KG@compareClusterResult %>% 
        dplyr::filter(grepl("pathway|metabolism", 
                            Description, 
                            ignore.case = TRUE)) %>%
        dplyr::mutate(
          Description = str_replace_all(
            Description, pattern = " - Mus musculus \\(house mouse\\)",
            replacement = ""))
      
      p4 <- 
        clusterProfiler::dotplot(comp.KG, showCategory = nCat) + 
          labs(title = 
                 paste0("CompareCluster ", FileNameBody, " KEGG"), x = "") + 
          theme_cp + 
          coord_flip()
            plot(p4)
              ggsave(paste0(FileNameBody, "_Compare.KEGG.png"),
                      width = Wid, height = Hei, dpi = 200)
    }
     
    ## Reactome
    comp.RA <- compareCluster(geneClusters = NamedGenelists, 
                              fun = "enrichPathway", 
                              organism = "mouse",
                              readable = TRUE)  
      if(as.data.frame(comp.RA) |> nrow() >= 1){
        
        p5 <-
          clusterProfiler::dotplot(comp.RA, showCategory = nCat) + 
            labs(title = paste0("CompareCluster ", FileNameBody, " Reactome"), x = "") + 
            theme_cp + 
            coord_flip()
              plot(p5)
                ggsave(paste0(FileNameBody, "_Compare.Reactome.png"),
                        width = Wid, height = Hei, dpi = 200)
        }
        
  ## saving results -----------------------------------
  saveRDS(list(BP = comp.BP, 
               CC = comp.CC,
               MF = comp.MF,
               KEGG = comp.KG,
               Reactome = comp.RA),
          paste0(FileNameBody, "_CompareClusterResults.rds"))
  
  df.BP <- as.data.frame(comp.BP) |> dplyr::mutate(Analysis = "GO:BP")
  df.CC <- as.data.frame(comp.CC) |> dplyr::mutate(Analysis = "GO:CC")
  df.MF <- as.data.frame(comp.MF) |> dplyr::mutate(Analysis = "GO:MF")
  df.KG <- as.data.frame(comp.KG) |>
    dplyr::mutate(Analysis = "KEGG") |>
    dplyr::select(-c(2:3))
  df.RA <- as.data.frame(comp.RA) |> dplyr::mutate(Analysis = "Reactome")
  
  # dataframeをexcelで保存
  rbind(df.BP, df.CC, df.MF, df.KG, df.RA) |>
    write.xlsx(paste0(FileNameBody, "_CompareClusterResults.xlsx"))

}

