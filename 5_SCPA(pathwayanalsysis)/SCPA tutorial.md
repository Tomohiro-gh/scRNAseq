## [SCPA: Single Cell Pathway Analysis](https://jackbibby1.github.io/SCPA/#tutorials)

Original paper: [J. A. Bibby et al., Cell Rep. 41, 111697 (2022).](https://www.cell.com/cell-reports/fulltext/S2211-1247(22)01571-6)

#### Data interpreatation
`compare_pathways`
```r
scpa_result <- compare_pathways(
     list(sample1, sample2, sample3),
     pathways = pathways)
```

 - https://jackbibby1.github.io/SCPA/articles/interpreting_scpa_output.html
 - https://rdrr.io/github/jackbibby1/SCPA/man/compare_pathways.html

#### values
Statistical results from the SCPA analysis. The qval should be the primary metric that is used to interpret pathway differences i.e. a higher qval translates to larger pathway differences between conditions.
If only two samples are provided, a fold change (FC) enrichment score will also be calculated. The FC statistic is generated from a running sum of mean changes in gene expression from all genes of the pathway. It's calculated from average pathway expression in population1 - population2, so a negative FC means the pathway is higher in population2.

