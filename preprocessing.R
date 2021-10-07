load10X <- function(filepath) {
  
  gem <- as.data.frame(as.matrix(Seurat::Read10X(data.dir = filepath)))
  
  # make names upper case and unique
  rownames(gem) <- make.names(toupper(rownames(gem)), unique=TRUE)
  
  gem
}

prepareData_10X <- function(
    filepath,  
    pct.Cells = 0.2,
    minNumGenes = 300,
    maxNumGenes = 5000,
    minNumHk = 65,
    maxPct.Mito = 2
){

    ## -----------------------------------------------------------------------------------
    ## Step 1: Read data
    message("Loading data")
    gem <- load10X(filepath)
  
    print(sprintf('size of data before filtering nGenes=%s, nCells=%s', nrow(gem), ncol(gem)))
  
    ## -----------------------------------------------------------------------------------
    ## Step 2: filter cells
  
    message("---- Cell filtering")
    hkGenes <- read.table("GSE139999-data-analysis/data/housekeeping_genes.txt", stringsAsFactors = FALSE)$V1
    hkGenes <- intersect(toupper(hkGenes), rownames(gem))
  
    # Find mitochondrial genes
    mitoGenes <- rownames(gem)[startsWith(rownames(gem), "MT")]
  
    meta.cellFilter <- data.frame(
        cellID = colnames(gem),
        keep = rep(TRUE, ncol(gem)),
        numHk = apply(gem[hkGenes, ], 2, function(x) sum(x > 0)),
        numGenes = apply(gem, 2, function(x) sum(x > 0)),
        pct.Mito = rep(NA,ncol(gem)),
        pct.Hk = rep(NA,ncol(gem)),
        reason = rep(NA, ncol(gem))
    )
  
    meta.cellFilter$pct.Mito <- round(apply(
        gem[mitoGenes, ], 
        2, 
        function(x) sum(x>0))*100/meta.cellFilter$numGenes, 
        2
    )
    meta.cellFilter$pct.Hk <- round(
        meta.cellFilter$numHk * 100 / meta.cellFilter$numGenes, 
        2
    )
  
    thresholds <- data.frame(
        nCellsPostFilter = NA, 
        nGenesPostFilter = NA, 
        pct.Cells = pct.Cells, 
        minNumGenes = minNumGenes, 
        maxNumGenes = maxNumGenes, 
        minNumHk = minNumHk, 
        maxPct.Mito = maxPct.Mito
    )
    
    IDX <- meta.cellFilter$numGenes < thresholds$minNumGenes | meta.cellFilter$numGenes > thresholds$maxNumGenes
  
    print(sprintf("Filtering %s cells ..", sum(IDX)))
    meta.cellFilter[IDX, "keep"] <- FALSE
    meta.cellFilter[IDX, "reason"] <- paste(meta.cellFilter[IDX, "reason"], "low_or_high_num_genes", sep = ";")
  
    IDX <- meta.cellFilter$numHk < thresholds$minNumHk
  
    meta.cellFilter[IDX, "keep"] <- FALSE
    meta.cellFilter[IDX, "reason"] <- paste(meta.cellFilter[IDX, "reason"], "low_num_hk", sep = ";")
  
    IDX <- meta.cellFilter$pct.Mito > thresholds$maxPct.Mito
  
    meta.cellFilter[IDX, "keep"] <- FALSE
    meta.cellFilter[IDX, "reason"] <- paste(meta.cellFilter[IDX, "reason"], "high_mito", sep = ";")
  
    keep_cells <- rownames(meta.cellFilter)[meta.cellFilter$keep]
  
    ## -----------------------------------------------------------------------------------
    ## Step 3: filter genes
    message("---- Gene filtering")
    
    gem_filt <- gem[, meta.cellFilter$cellID[meta.cellFilter$keep]]
  
    meta.geneFilter <- data.frame(
        geneSymbol = rownames(gem_filt),
        keep = rep(TRUE, nrow(gem_filt)),
        pct.Cells = apply(gem_filt, 1, function(x) round(sum(x>0)*100/ncol(gem_filt), 2))
    )
  
    meta.geneFilter$keep[meta.geneFilter$pct.Cells < thresholds$pct.Cells] <- FALSE
  
    keep_genes <- rownames(meta.geneFilter)[meta.geneFilter$keep]
  
    ## filter cells and genes and return filtered GEM
    fGem <- gem[keep_genes, keep_cells]
  
    print(sprintf('size of gem after filtering nGenes = %s, nCells = %s', nrow(fGem), ncol(fGem)))
    
    fGem
}









