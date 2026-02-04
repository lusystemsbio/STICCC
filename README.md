# sticcc
State transition inference using cross-cell correlations


## Introduction

This R package implements STICCC (State Transition Inference from Cross-Cell Correlations, a method for inferring directional state transition vectors for single-cell gene expression data with respect to a user-supplied gene regulatory network (GRN). Rather than modeling transcriptional dynamics directly, STICCC avoids the need for time labels or splicing kinetics by exploiting cross-cell gene–gene correlation structure within local neighborhoods of state space. This makes it particularly useful for:
* datasets where multiple timepoints or unspliced transcripts are not available
* multistable systems where bidirectional transitions are expected
* datasets where the underlying GRN is reasonably well known

Intuitively, the rationale behind STICCC is:
* Neighboring cells are assumed to sample the same dynamical regime
* Expression of regulator genes should correlate with future (sign-corrected) expression of their targets (cross-cell correlation)
* Nearby cells with cross-correlated gene expression represent likely future gene expression states

The input is a GRN, structured as a table with three columns Source, Target, and Type which describe a set of interactions between genes/nodes (type 1 = activation, type 2 = inhibition). Users must also provide a gene expression matrix corresponding to the GRN (genes x cells), either from experiment or simulations. STICCC will then produce directional predictions about likely future gene expression states for each cell, as well as estimate the propensity for reversible transitions between states in the dataset. For more detailed discussion of the method and further example applications, see the article "Dissecting reversible and irreversible single cell state transitions from gene regulatory networks" referenced below.


## Installation
STICCC can be installed via GitHub using devtools or remotes: ``devtools::install_github("lusystemsbio/sticcc")``

Then call the package using ``library(STICCC)``


## Typical workflow
*1. Construct or obtain GRN*. Below is a toy example showing the format:
```
topo <- data.frame(Source=c("A","B","C"),Target=c("B","C","A"),Type=c(2,2,2))
```
*2. Simulate or obtain gene expression data*. Below is an example simulation using the package sRACIPE which will produce a gene expression matrix:
```
nSamples <- 1000
racipe <- simTopo(topo, numModels = nSamples)
```

*3. Preprocess gene expression data*. Log-normalize, optionally subset to GRN genes - below this is done using the built-in function in sRACIPE, but similar implementations are available in, e.g., Seurat, or can easily be manually implemented. The expression matrix should have genes as rows and cells as columns.
```
racipe_norm <- sracipeNormalize(racipe)
exprMat <- assay(racipe)
exprMat_norm <- assay(racipe_norm)
```

*4. Prepare STICCC object*. Provide the GRN and gene expression data, as well as hyperparameters `radius`, `minNeighbors`, `nPCs`, `plotDim` (see practical tips below)
```
stic <- sticSE(topo = topo, exprMat = exprMat, normData = exprMat_norm,
             topoName = topoName, expName = "repressilator_example")

# Cluster gene expression data for visualizations later
stic <- prepMetadata(sce = stic, exprMat = exprMat, cluster = T, k = 6)
```

*5. Create embeddings*. STICCC supports PCA, UMAP, and TSNE for projection. There is a built-in function to compute PCA, or you can manually add embedding coordinates corresponding to the name given in the parameter `plotDim`. Note that using UMAP or TSNE for distance calculations is strongly discouraged as these methods do not preserve global distances.
```
# This method will run PCA and store it appropriately
stic <- runPCA(stic)

# Alternatively, you can manually enter a supported embedding type as below:
reducedDim(stic, "UMAP") <- UMAP_coords
stic@metadata$params$plotDim <- "UMAP"

# You can also add PCA coordinates manually:
reducedDim(stic, "PCA") <- pca_df
stic@metadata$pca_data <- list(sdev=NA, rotation=pca_loadings, center=NA, scale=NA) # loadings just needed for PCA plots that show them
stic@metadata$pca_summary <- summary(pca) # optional
```

*6. Calculate pairwise distances*. To determine the sampling radius, we must estimate the maximum pairwise distance between cells. By default, distance between cells is calculated based on the first 10 PCs, and the max is estimated based on a random subset of the input data. You can also calculate it in full gene expression space or with more/fewer PCs by adjusting `useGenes` and `nComponents`. For smaller datasets, you can also set `est` to FALSE to calculate the full distance matrix. 
```
# compute grid based on PCA for later smoothing
stic <- computeGrid(stic)

# compute pairwise distance between points
stic <- computeDist(stic)

```

*7. Run STICCC*. The full v1 and v2 vectors are available in `stic@metadata$vectors` and `stic@metadata$vectors_in` respectively, while information about each sample is stored in `colData(stic)` and parameters in `stic@metadata$params`.
```
stic <- runSTICCC(stic)
```

*8. Visualization & downstream analysis*. You can adjust the size of the drawn vectors using `scalingFactor` and set `colorVar` to any valid column of the metadata. Adjust the arguments to `computeGridVectors` to change which outputs are shown from `plotGrid` (see below)
```
# plot individual vectors (v1)
plotVectors(stic, scalingFactor = 0.4,
            colorVar = "Cluster")

# grid-based smoothing of velocities
# To plot v1: inVectors=F, combine=F, how=NA
# To plot v2: inVectors=T, combine=F, how=NA
# To plot net flow: inVectors=F, combine=T, how="net"
# To plot reversibility: inVectors=F, combine=T, how="rev"
stic <- computeGridVectors(stic, unitVectors=F, inVectors=F, combine=F, how=NA)

# plot grid-smoothed vectors
plotGrid(stic, colorVar = "Cluster")

```


## Practical tips
* GRN construction: STICCC is moderately robust to errors in the input GRN based on simulation results, but neverthless supplying an incorrect GRN can produce nonsense results. Manual curation is always recommended, though several computational tools exist to infer a GRN from gene expression data. As long as the network can be arranged into the appropriate input format, these can be used in conjunction with STICCC. 

* Neighborhood radius (radius): This controls how smooth the inferred vector field is, and is a proportion of the maximum pairwise distance between cells/samples in the input dataset. Lowering it will increase the resolution, but will cause some cells in sparse regions not to have sufficient neighbors for a prediction. For simulated data and large datasets, 0.05 is a reasonable starting point, while for smaller experimental datasets you may need to increase it. You may also modify minNeighbors to allow/reject predictions for cells in sparser regions.

* Dimension reduction: PCA, tSNE, and UMAP are supported, but PCA is recommended. In this case, you should set nPCs to a reasonable number that captures the essential features of your dataset (>80% of variance ideally). Relatedly, vectors can be computed either in gene space or embedding space via the parameter useOriginalFeatures.


## Vignettes
For an example workflow, see the vignettes folder which shows the analysis of a small simulated gene circuit.
More detailed analysis scripts are available at: https://github.com/lusystemsbio/sticcc_analysis

## References
For more information, see the original article: Ramirez D, Lu M. Dissecting reversible and irreversible single cell state transitions from gene regulatory networks. Mol Syst Biol. 2026. doi:10.1038/s44320-026-00196-8.


