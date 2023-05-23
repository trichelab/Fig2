---
title: "Buteyn et al. Figure 2"
author: "Tim Triche and Zach DeBruine"
date: "2023-05-19"
output:
  html_document:
    fig_width: 8
    code_folding: hide
    highlight: tango
    toc: yes
    keep_md: true
    toc_float:
      collapsed: yes
---





## Objective

To map the transcriptional signatures of immune-cold transplant-refractory AML.

In particular, to find out why FUS::ERG-driven AML is consistently immune-cold,
to determine if specific biological processes are shared across transplant-
refractory AML cases, and to see if this outcome can be predicted at diagnosis.


## Approach

### Steps 

Using the previously fitted rank-90 non-negative matrix factorization model, 

1. Identify factors describing specific signatures for cytogenetic subtypes.

2. Investigate biological processes underlying these factors using GSEA.

3. Determine if these factors are predictive of transplant-refractory AML.


## Analysis

We begin by loading the previously fitted NMF model (refer to TARGET_NMF.Rmd). 


```r
library(RcppML)
library(singlet)

# load the model 
nmf <- readRDS("nmf_90D_model.rds") 
dim(nmf@h) 
```

```
## [1]   90 1630
```

```r
# [1]   90 1630
```

Next we load the original annotated RNAseq data and check that everyone's there.


```r
# original
dx0 <- readRDS("dx.rds")
colnames(dx0) <- make.unique(colnames(dx0))

# SummarizedExperiment -> SingleCellExperiment
library(SingleCellExperiment) 

# normalize exactly the same way as Seurat
logNorm <- function(x) log1p(sweep(x, 2, sizeFactors(x), `/`))

# a Seurat-like SingleCellExperiment
toSCE <- function(se, scale_factor = 1e4) { 
  assayNames(se)[1] <- "counts"
  sce <- as(se, "SingleCellExperiment")
  sce$librarySize <- colSums(counts(sce))
  sizeFactors(sce) <- sce$librarySize / scale_factor
  logcounts(sce) <- log1p(sweep(counts(sce), 2, sizeFactors(sce), `/`))
  metadata(sce)$logcounts_are_log1p_cp10k <- TRUE # reminder 
  return(sce)
}

# load the original raw counts 
dx <- toSCE(dx0)

# everything there?
ncol(dx) == ncol(nmf@h)
```

```
## [1] FALSE
```

```r
message(ncol(dx) - ncol(nmf@h), " samples in dx are missing from NMF model")
# 224 samples in dx are missing from NMF model
```

### Slight hiccup 

Unfortunately we are missing some of the rarer fusions! We can fix that, though.
In the process we will recover a few others that slipped through earlier on. 
When we just need to update a model with some new observations, we can use the
singlet::ProjectData() function with our existing W matrix.


```r
library(singlet)
# needs to be the most recent version from github
stopifnot(packageVersion("singlet") >= "0.99.26") 
# THIS WAS A LOT OF WORK TO FIX, PLEASE CLAP!!!1
colnames(dx) <- make.unique(colnames(dx))
colnames(nmf@h) <- make.unique(colnames(nmf@h))

# without reordering:
dx <- ProjectData(dx, nmf@w, reduction.name="NNLS") # reducedDim(dx, "NNLS")
NMF <- reducedDim(dx, "NNLS")
sampleFactors(NMF)[colnames(nmf@h), ] <- t(nmf@h)

# with reordering:
dx <- ProjectData(dx, nmf@w, reorder=TRUE, 
                  reduction.name="NNLSre", reduction.key="NMFre_") 
NMFre <- reducedDim(dx, "NNLSre")
sampleFactors(NMFre)[colnames(nmf@h), ] <- t(nmf@h)

# compare with and without:
library(ComplexHeatmap)

unre <- sampleFactors(reducedDim(dx, "NNLS")) - sampleFactors(NMF) 
re <- sampleFactors(reducedDim(dx, "NNLSre")) - sampleFactors(NMFre)

# I'm not sure how I feel about this:
Heatmap(t(unre), cluster_rows=F, cluster_columns=F, name="unreordered - NMF") +
Heatmap(t(re), cluster_rows=F, cluster_columns=F, name="reordered - NMF") 
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

```r
# The differences are certainly smaller for unreordered, though.

# Let's just use the unreordered projections for now.
reducedDim(dx, "NMF") <- NMF
```

Let's take a quick look at a UMAP of (projected) dx NMF versus the original NMF.


```r
library(uwot)

# model with the original data
umap_fit_orig <- uwot::umap(t(nmf@h), 
                            min_dist = 0.3,
                            n_components = 2,
                            ret_model = TRUE, 
                            metric = "cosine")
umap_coords <- umap_fit_orig$embedding

# clusters not unlike Seurat's 
library(bluster)
bluster_cluster_orig <- clusterRows(t(nmf@h), 
                                    NNGraphParam(cluster.fun="louvain"))

df <- data.frame("sample" = rownames(umap_coords),
                 "cluster" = bluster_cluster_orig, 
                 "fusion" = nmf@misc$covs$Primary.Fusion,
                 "umap1" = umap_coords[,1], 
                 "umap2" = umap_coords[,2])

p1 <- ggplot(df, aes(umap1, umap2, color = fusion)) + 
        geom_point() + 
        theme_bw() + 
        theme(aspect.ratio = 1,
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.text.x = element_blank(), 
              axis.text.y = element_blank(), 
              axis.title = element_text(hjust = 0, vjust = 0), 
              axis.line = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5)) +
        labs(x = "UMAP1", 
             y = "UMAP2", 
             color = "Fusion", 
             title = "RNA-seq by\nfusion status") + 
        theme(axis.ticks = element_blank()) + 
        NULL

df$fus_erg <- as.character(df$fusion)
df$fus_erg[df$fus_erg != "FUS-ERG"] <- "Other"

p2 <- ggplot(df, aes(umap1, umap2, color = fus_erg)) + 
        geom_point() + 
        theme_bw() + 
        theme(aspect.ratio = 1,
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.text.x = element_blank(), 
              axis.text.y = element_blank(), 
              axis.title = element_text(hjust = 0, vjust = 0), 
              axis.line = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5)) +
        labs(x = "UMAP1", 
             y = "UMAP2", 
             color = "Fusion", 
             title = "RNA-seq by\nfusion status") + 
        theme(axis.ticks = element_blank()) + 
        NULL

library(cowplot)
plot_grid(p1, p2, nrow = 1, align = "v")
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

```r
ggsave("dx_NMF_orig_umap.png", width = 10, height = 5)

metadata(dx)$umap_fit <- umap_fit_orig
# saveRDS(dx, file="dx_with_NMF_and_UMAP.rds")
```

It's not exactly identical to Zach's, but the FUS-ERG cases do split off. 
That's good enough for the time being. Let's add in the other fusions now.


```r
# load the covariates
dx$USI <- substr(colnames(dx), 1, 6) # avoid duplicate weirdness 
target_covs <- read.csv("target_covs.eligible.csv", row=1)
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
stopifnot(all(dx$Primary.Fusion == target_covs[dx$USI, "Primary.Fusion"]))
```

```
## Error in eval(expr, envir, enclos): object 'target_covs' not found
```

```r
# fusionGroup for panel 2A
dx$fusion <- dx$Primary.Fusion
fusions_list <- read.table("fusions.txt")[, 1]

# ERG-HNRNPH1 == HNRNPH1-ERG
dx$fusion[dx$fusion == "ERG-HNRNPH1"] <- "HNRNPH1-ERG"

# unspaceify
dx$fusion[dx$fusion == "KMT2A-MLLT3 "] <- "KMT2A-MLLT3"

# We have to label the ETV6-X fusions since they're a grab bag
ETV6fusions <- grep("ETV6", levels(factor(dx$fusion)), value=TRUE)
ETV6_X <- setdiff(ETV6fusions, "ETV6-MNX1")
dx$fusion[dx$fusion %in% ETV6_X] <- "ETV6-X"
mainExpName(dx) <- "RNA"

# now subset the data
length(which(dx$fusion %in% fusions_list)) 
```

```
## [1] 1682
```

```r
# [1] 1682
dx <- dx[, which(dx$fusion %in% fusions_list)]
dx$fusion <- factor(dx$fusion)
dx$fusion <- relevel(dx$fusion, which(levels(dx$fusion) == "None"))
ncol(dx)
```

```
## [1] 1682
```

```r
# Catalog them
library(knitr)
kable(sort(table(fusion=dx$fusion)))
```



|fusion        | Freq|
|:-------------|----:|
|KMT2A-MLLT11  |    7|
|NPM1-MLF1     |    8|
|HNRNPH1-ERG   |    9|
|RBM15-MKL1    |   13|
|RUNX1-CBFA2T3 |   13|
|FUS-ERG       |   14|
|KAT6A-CREBBP  |   14|
|KMT2A-SEPT6   |   14|
|ETV6-MNX1     |   16|
|KMT2A-MLLT1   |   21|
|ETV6-X        |   25|
|NUP98-KDM5A   |   38|
|CBFA2T3-GLIS2 |   42|
|KMT2A-ELL     |   49|
|KMT2A-MLLT4   |   49|
|DEK-NUP214    |   55|
|KMT2A-MLLT10  |   92|
|KMT2A-MLLT3   |  126|
|NUP98-NSD1    |  137|
|CBFB-MYH11    |  186|
|RUNX1-RUNX1T1 |  218|
|None          |  536|


Let's feed the additional fusions, projected onto the original NMF, to the 
original UMAP model. 


```r
# project the additional observations into the original UMAP model

u <- umap_transform(sampleFactors(reducedDim(dx, "NMF")), metadata(dx)$umap_fit)
reducedDim(dx, "UMAP") <- u
# saveRDS(dx, file="dx_subset_with_NMF_and_UMAP.rds")
```

We can fit factors to fusions:


```r
library(limma) # fit NMF factors
design <- model.matrix(~ fusion, data=colData(dx))
nmfdat <- t(sampleFactors(reducedDim(dx, "NMF")))[, rownames(design)]
fit <- eBayes(lmFit(nmfdat, design))

rownames(subset(topTreat(fit, coef="fusionFUS-ERG"), B > 0)) # 85
```

```
## [1] "85"
```

```r
which.max(featureLoadings(reducedDim(dx, "NMF"))["EZH2", ]) # 85
```

```
## [1] 85
```

```r
which.max(featureLoadings(reducedDim(dx, "NMF"))[, 85]) # STAB1
```

```
## STAB1 
##   261
```

```r
rownames(subset(topTreat(fit, coef="fusionCBFA2T3-GLIS2"), B > 0)) # 65 73 21 48
```

```
## [1] "65" "73" "21" "48" "27"
```

```r
which.max(featureLoadings(reducedDim(dx, "NMF"))["FOLR1", ]) # 65
```

```
## [1] 65
```

```r
which.max(featureLoadings(reducedDim(dx, "NMF"))[, 73]) # CCND2
```

```
## CCND2 
##  4860
```

```r
rownames(subset(topTreat(fit, coef="fusionHNRNPH1-ERG"), B > 0)) # none
```

```
## character(0)
```

```r
rownames(subset(topTreat(fit, coef="fusionETV6-MNX1"), B > 0)) # 76
```

```
## [1] "76"
```

```r
which.max(featureLoadings(reducedDim(dx, "NMF"))[, 76]) # TP53INP1
```

```
## TP53INP1 
##    11479
```

```r
rownames(subset(topTreat(fit, coef="fusionETV6-X"), B > 0)) # 21 2
```

```
## [1] "21" "2"
```

```r
which.max(featureLoadings(reducedDim(dx, "NMF"))[, 21]) # BAHCC1 
```

```
## BAHCC1 
##  47621
```

```r
which.max(featureLoadings(reducedDim(dx, "NMF"))[, 2]) # MALAT1
```

```
## MALAT1 
##  39004
```

We can also fit some more interesting outcomes, albeit only for exploration:


```r
dx$inductionFailure <- ifelse(dx$EFS.event.type.ID == "Induction failure", 1, 0)
designIF <- model.matrix(~ inductionFailure, data=colData(dx))
nmfIF <- t(sampleFactors(reducedDim(dx, "NMF")))[, rownames(designIF)]
fitIF <- eBayes(lmFit(nmfIF, designIF))
rownames(subset(topTreat(fitIF, n=90), B > 0))
```

```
##  [1] "16" "40" "87" "82" "45" "2"  "44" "46" "1"  "39" "48" "8"  "4"  "55" "90"
## [16] "15" "14" "5"  "68" "52" "7"  "12" "81" "9"  "78" "53" "89" "42" "11" "77"
## [31] "43" "31" "20" "3"  "10" "26" "6"  "80" "60" "29" "35" "34" "50" "18" "30"
## [46] "13" "47" "32" "86" "71" "33" "72" "17" "21" "28" "88" "36" "63" "22" "25"
## [61] "23" "61" "38" "62" "58" "19" "24" "41" "37" "69" "49" "27" "66" "54" "51"
## [76] "56" "57" "59" "83" "64" "70" "75" "67" "74" "76" "73" "79" "65" "85"
```

```r
dx$transplanted <- ifelse(dx$SCT.in.1st.CR == "Yes", TRUE, FALSE)
dxSCT <- dx[, which(dx$transplanted)]
dxSCT$died <- ifelse(dxSCT$OS.event.ID == "Dead", 1, 0)
designSCT <- model.matrix(~ died, data=colData(dxSCT))
nmfSCT <- t(sampleFactors(reducedDim(dxSCT, "NMF")))[, rownames(designSCT)]
fitSCT <- eBayes(lmFit(nmfSCT, designSCT))
rownames(subset(topTreat(fitSCT, n=90), B > 0))
```

```
##  [1] "40" "87" "16" "82" "44" "2"  "45" "39" "52" "68" "15" "8"  "46" "1"  "14"
## [16] "4"  "48" "55" "5"  "90" "21" "81" "37" "9"  "7"  "31" "77" "89" "12" "78"
## [31] "53" "3"  "42" "20" "10" "43" "26" "6"  "30" "29" "80" "35" "51" "28" "18"
## [46] "34" "50" "33" "60" "13" "58" "86" "11" "72" "71" "47" "62" "23" "36" "66"
## [61] "41" "83" "38" "88" "49" "27" "32" "70" "22" "69"
```

We can also quickly check things out with iSEE:


```r
if (FALSE) {

  library(iSEE)

  if (!exists("dx")) dx <- readRDS("dx_subset_with_NMF_and_UMAP.rds")

  UMAP_plot <- new("ReducedDimensionPlot", 
                   Type = "UMAP", 
                   ColorBy = "Column data", 
                   ColorByColumnData = "fusion",
                   TooltipColumnData = "fusion")

  NMF_plot <- new("ReducedDimensionPlot", 
                  Type = "NMF", 
                  XAxis = 85L,
                  YAxis = 2L,
                  SelectionAlpha = 0.01, 
                  ColorBy = "Column data", 
                  ColorByColumnData = "fusion",
                  TooltipColumnData = "fusion",
                  ColumnSelectionSource = "ReducedDimensionPlot1", 
                  ColumnSelectionDynamicSource = TRUE) 
      
  COVS_plot <- new("ColumnDataTable", 
                 HiddenColumns = setdiff(names(colData(dx)),
                                         c("fusion", "Primary.Fusion", "Sex")),
                 ColumnSelectionSource = "ReducedDimensionPlot2", 
                 ColumnSelectionDynamicSource = TRUE)

  GENE <- new("RowDataTable",
              Selected="EZH2")
           
  EXPRESSION <- new("FeatureAssayPlot", 
                    YAxisFeatureName = "EZH2", 
                    YAxisFeatureSource = "RowDataTable1",
                    YAxisFeatureDynamicSource = TRUE,
                    Assay = "logcounts", 
                    XAxis = "Column data", 
                    XAxisColumnData = "fusion", 
                    ColorByColumnData = "fusion", 
                    ColorBy = "Column data") 

  COVS2 <- new("ColumnDataTable", 
               HiddenColumns = setdiff(names(colData(dx)),
                                       c("Primary.Fusion",
                                         "EFS.event.type.ID",
                                         "OS.event.ID")),
               ColumnSelectionSource = "ReducedDimensionPlot2", 
               ColumnSelectionDynamicSource = TRUE)

  iSEE(dx, 
       initial = list(UMAP = UMAP_plot, 
                      NMF = NMF_plot,  
                      COVS = COVS_plot,
                      GENE = GENE,
                      EXPRESSION = EXPRESSION,
                      COVS2 = COVS2))

}
```
Well, that looks good. 


### Figure panels: 

A. UMAP of TARGET pAML (n = 1682) samples colored by cytogenetic subtype.


```r
df <- data.frame("sample" = colnames(dx), 
                 "fusion" = dx$fusion, 
                 "umap1" = reducedDim(dx, "UMAP")[,1],
                 "umap2" = reducedDim(dx, "UMAP")[,2])

p1 <- ggplot(df, aes(umap1, umap2, color = fusion)) + 
        geom_point() + 
        theme_bw() + 
        theme(aspect.ratio = 1,
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.text.x = element_blank(), 
              axis.text.y = element_blank(), 
              axis.title = element_text(hjust = 0, vjust = 0), 
              axis.line = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5)) +
        labs(x = "UMAP1", 
             y = "UMAP2", 
             color = "Fusion", 
             title = "RNA-seq by fusion status") + 
        theme(axis.ticks = element_blank()) + 
        NULL

df$fus_erg <- as.character(df$fusion)
df$fus_erg[df$fus_erg != "FUS-ERG"] <- "Other"

p2 <- ggplot(df, aes(umap1, umap2, color = fus_erg)) + 
        geom_point() + 
        theme_bw() + 
        theme(aspect.ratio = 1,
              panel.grid = element_blank(),
              panel.border = element_blank(),
              axis.text.x = element_blank(), 
              axis.text.y = element_blank(), 
              axis.title = element_text(hjust = 0, vjust = 0), 
              axis.line = element_line(colour = "black"),
              plot.title = element_text(hjust = 0.5)) +
        labs(x = "UMAP1", 
             y = "UMAP2", 
             color = "Fusion", 
             title = "RNA-seq by fusion status") + 
        theme(axis.ticks = element_blank()) + 
        NULL

library(cowplot)
plot_grid(p1, p2, nrow = 1, align = "v")
```

![plot of chunk unnamed-chunk-10](figure/unnamed-chunk-10-1.png)

```r
ggsave("dx_NMF_updated_umap.png", width = 10, height = 5)
```

B. Association of individual factors with FUS::ERG vs. all other pAML. 

C. Enrichnment of factors in cytogenetic subtypes of pediatric AML.

D. Pathway enrichment for factors associated with allo-SCT response.

E. Antigen presentation and immune response GSEA, FUS::ERG vs. other pAML.

