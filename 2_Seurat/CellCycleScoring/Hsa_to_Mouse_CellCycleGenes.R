# CellCycleGeneを抜き出す
library(Seurat)
library(homologene) # cc.genesがロードされる

cc.gene # seuratをロードすると見える
length(s.genes) #43
length(g2m.genes) #54

homologene::taxData
# tax_id                      name_txt
# 1   10090                  Mus musculus
# 2   10116             Rattus norvegicus
# 3   28985          Kluyveromyces lactis
# 4  318829            Magnaporthe oryzae
# 5   33169         Eremothecium gossypii
# 6    3702          Arabidopsis thaliana
# 7    4530                  Oryza sativa
# 8    4896     Schizosaccharomyces pombe
# 9    4932      Saccharomyces cerevisiae
# 10   5141             Neurospora crassa
# 11   6239        Caenorhabditis elegans
# 12   7165             Anopheles gambiae
# 13   7227       Drosophila melanogaster
# 14   7955                   Danio rerio
# 15   8364 Xenopus (Silurana) tropicalis
# 16   9031                 Gallus gallus
# 17   9544                Macaca mulatta
# 18   9598               Pan troglodytes
# 19   9606                  Homo sapiens
# 20   9615        Canis lupus familiaris
# 21   9913                    Bos taurus

homologene("VEGFA", inTax = 9606, outTax = 10090)

s.gene.mouse = 
  cc.genes$s.genes %>% 
  homologene(., inTax = 9606, outTax = 10090) %>% 
  pull("10090")

g2m.gene.mouse = 
  cc.genes$g2m.genes %>% 
  homologene(., inTax = 9606, outTax = 10090) %>% 
  pull("10090")

cc.genes.mouse = 
  list(
    s.gene.mouse = s.gene.mouse,
    g2m.gene.mouse = g2m.gene.mouse
  )

saveRDS(
  cc.genes.mouse, 
  "/Users/tomohiro/Dropbox/BioInfomatics/R_CommandList/Seurat/CellCycleScoring/CellCylce.genes.mouse.rds")

  
