#### Load required packages and set working directory.
{
  library(edgeR)
  library(limma)
  library(Glimma)
  library(gplots)
  library(org.Hs.eg.db)
  library(RColorBrewer)
  library(PCAtools)
  library(ggplot2)
  library(Biobase)
  library(Homo.sapiens)
  library(EnhancedVolcano)
  library(sva)
  library(rhdf5)
  library(gprofiler2)
  library(GSVA)
  library(qusage)
  library(RRHO2)
  setwd("~/data_release_files")
}

#### DGE Analysis of whole blood patient samples from Overmyer et al. (2021), GEO series GSE157103.
{
  # Read in counts table from GSE157103
  seqdata <- read.delim("GSE157103_genes_ec.tsv", stringsAsFactors = FALSE, row.names = 1)
  
  # Remove samples that had poor relative library size; removing them made the downstream analysis more robust
  seqdata <- seqdata[,-c(2,37,50,52)]
  
  # Create the DGEList object "dgeCounts" for limma/voom analysis.
  dgeCounts <- DGEList(seqdata)
  
  # Read in experiment design table and match sampleID orders between "dgeCounts" and "experimentDesign".
  # Experiment design table contains information about COVID-19 status and ICU status, which will be used to
  #   set up contrast matrices later.
  experimentDesign <- read.table("experiment_design_GSE157103.txt", header = T, sep = "\t")
  rownames(experimentDesign) <- experimentDesign$SampleID
  experimentDesign <- experimentDesign[colnames(dgeCounts$counts),]
  
  # Add ENTREZ gene names to dgeCounts and remove duplicate genes.
  geneid <- rownames(dgeCounts)
  genes <- select(org.Hs.eg.db, keys = geneid, columns = c("ENTREZID"), keytype = "SYMBOL")
  genes <- genes[!duplicated(genes$ENTREZID),]
  genes <- genes[match(rownames(dgeCounts$counts), genes$SYMBOL),]
  
  # Add group information to "dgeCounts" for setting up contrast matrices.
  group <- factor(experimentDesign$Status)
  dgeCounts$samples$group <- group
  
  # Add ENTREZ gene information to "dgeCounts".
  dgeCounts$genes <- genes
  
  # Filter out poorly expressed genes (<5 counts per million across the sample size of the smallest group).
  keptExpression <- filterByExpr(dgeCounts, group = group, min.count = 5)
  dgeCounts <- dgeCounts[keptExpression,, keep.lib.sizes=FALSE]
  
  # Perform TMM normalization on the counts data.
  dgeCounts <- calcNormFactors(dgeCounts, method = "TMM")
  
  # Perform batch effect correction with SVA.
  fullmatrix <- model.matrix(~0+group)
  nullmatrix <- model.matrix(~1, data = dgeCounts$samples)
  svaobj <- sva(dgeCounts$counts, fullmatrix, nullmatrix)
  
  # Create design and contrast matrices to compare COVID ICU vs. nonCOVID ICU (A) and COVID ICU vs. COVID nonICU (B).
  # Design matrix has 2 additional columns to account for 2 surrogate variables detected by SVA.
  design <- cbind(fullmatrix,svaobj$sv)
  colnames(design)<- gsub("group","",colnames(design))
  colnames(design) <- c(colnames(design)[1:4], "V5", "V6")
  contrastMatrixA <- makeContrasts(CovidICU_vs_nonCovidICU = COVID_ICU - nonCOVID_ICU,
                                   levels = colnames(design))
  contrastMatrixB <- makeContrasts(CovidICU_vs_CovidnonICU = COVID_ICU - COVID_nonICU,
                                   levels = colnames(design))
  
  # Perform voom normalization on "dgeCounts".
  y <- voom(dgeCounts, design)
  
  # Perform limma DGE analysis for all contrasts as described above (labelled A and B)
  # Absolute logFC > 1 and adjusted p-value < 0.05.
  fitA <- lmFit(y, design)
  fitA <- contrasts.fit(fitA, contrastMatrixA)
  efitA <- eBayes(fitA)
  etA <- decideTests(efitA)
  coefficientVarA = colnames(fitA$coefficients)
  limmaResA_full <- topTable(efitA, coef = coefficientVarA, n = dim(fitA)[1])
  limmaResA <- limmaResA_full[which(abs(limmaResA_full$logFC) > 1 & limmaResA_full$adj.P.Val < 0.05),]
  
  fitB <- lmFit(y, design)
  fitB <- contrasts.fit(fitB, contrastMatrixB)
  efitB <- eBayes(fitB)
  etB <- decideTests(efitB)
  coefficientVarB = colnames(fitB$coefficients)
  limmaResB_full <- topTable(efitB, coef = coefficientVarB, n = dim(fitB)[1])
  limmaResB <- limmaResB_full[which(abs(limmaResB_full$logFC) > 1 & limmaResB_full$adj.P.Val < 0.05),]
  
  # This block contains same code as above but instead using "COVID" and "nonCOVID" identifiers for the COVID vs. nonCOVID comparison (C)
  {
    dgeCounts <- DGEList(seqdata)
    geneid <- rownames(dgeCounts)
    genes <- select(org.Hs.eg.db, keys = geneid, columns = c("ENTREZID"), keytype = "SYMBOL")
    genes <- genes[!duplicated(genes$ENTREZID),]
    genes <- genes[match(rownames(dgeCounts$counts), genes$SYMBOL),]
    group <- factor(experimentDesign$Status_two)
    dgeCounts$samples$group <- group
    dgeCounts$genes <- genes
    keptExpression <- filterByExpr(dgeCounts, group = group, min.count = 5)
    dgeCounts <- dgeCounts[keptExpression,, keep.lib.sizes=FALSE]
    dgeCounts <- calcNormFactors(dgeCounts, method = "TMM")
    fullmatrix <- model.matrix(~0+group)
    nullmatrix <- model.matrix(~1, data = dgeCounts$samples)
    svaobj <- sva(dgeCounts$counts, fullmatrix, nullmatrix)
    # Design matrix has 2 additional columns to account for 2 surrogate variables detected by SVA.
    design <- cbind(fullmatrix,svaobj$sv)
    colnames(design)<- gsub("group","",colnames(design))
    colnames(design) <- c(colnames(design)[1:2], "V3", "V4")
    contrastMatrixC <- makeContrasts(Covid_vs_nonCovid = COVID - nonCOVID,
                                     levels = colnames(design))
    y <- voom(dgeCounts, design)
    fitC <- lmFit(y, design)
    fitC <- contrasts.fit(fitC, contrastMatrixC)
    efitC <- eBayes(fitC)
    etC <- decideTests(efitC)
    coefficientVarC = colnames(fitC$coefficients)
    limmaResC_full <- topTable(efitC, coef = coefficientVarC, n = dim(fitC)[1])
    limmaResC <- limmaResC_full[which(abs(limmaResC_full$logFC) > 1 & limmaResC_full$adj.P.Val < 0.05),]
  }
  
  # Output DGE analysis tables as csv files.
  filenameA <- paste("~/DEgenes_voom_", coefficientVarA, "_pval0.05_logFC1.csv", sep = "")
  filenameB <- paste("~/DEgenes_voom_", coefficientVarB, "_pval0.05_logFC1.csv", sep = "")
  filenameC <- paste("~/DEgenes_voom_", coefficientVarC, "_pval0.05_logFC1.csv", sep = "")
  write.csv(limmaResA, file = filenameA, row.names = T, quote = F)
  write.csv(limmaResB, file = filenameB, row.names = T, quote = F)
  write.csv(limmaResC, file = filenameC, row.names = T, quote = F)
  
  # Generate volcano plots to summarize DGE profiles for each comparison
  limmaResA_full <- limmaResA_full[order(rownames(limmaResA_full)),]
  limmaResB_full <- limmaResB_full[order(rownames(limmaResB_full)),]
  limmaResC_full <- limmaResC_full[order(rownames(limmaResC_full)),]
  limmaResA <- limmaResA[order(rownames(limmaResA)),]
  limmaResB <- limmaResB[order(rownames(limmaResB)),]
  limmaResC <- limmaResC[order(rownames(limmaResC)),]

#### DGE Analysis of bronchoalveolar lavage fluid patient samples from Liao et al. (2020), GEO series GSE145926.
{
  # scRNA-seq data must first be collapsed into bulk RNA-seq expression tables by aggregate mRNA counts from single cells to their respective patients.
  # The aggregated counts table will be stored in "countsdata".
  h5_file_list <- list.files(path = "GSE145926_expression_matrices_download", full.names = TRUE)
  idx <- data.frame("index" = c(1:33538))
  for(i in c(1:length(h5_file_list))){
    # Open each .h5 file
    data <- H5Fopen(h5_file_list[i])
    feature_list <- data$matrix$features$id
    
    # Subset indexed counts data and order them by index value
    data <- data.frame("index" = data$matrix$indices, "data" = data$matrix$data)
    data <- data[order(data$index, decreasing = FALSE),]
    
    # Aggregate counts data from single-cell to sample-level resolution
    data <- aggregate(data$data, list(data$index), sum)
    colnames(data) <- c("index", unlist(strsplit(h5_file_list[i], split = "_"))[5])
    
    # Add genes with zero counts
    data <- merge(idx, data, by = "index", all.x = TRUE)
    data[,2][which(is.na(data[,2]) == TRUE)] <- 0
    
    # Create and append sample-level gene expression to one expression matrix
    if(i == 1){
      countsdata <- data
    }else{
      countsdata <- merge(countsdata, data, by = "index")
    }
    
    # At the end, acquire the indexed list of features (same for every .h5 file) and append it to countsdata
    if(i == length(h5_file_list)){
      rownames(countsdata) <- feature_list[1:dim(countsdata)[1]]
      countsdata <- countsdata[,-1]
    }
  }
  
  # Save the aggregated counts data as a .csv for easy reference.
  write.csv(countsdata, "GSE145926_expression_matrix_aggregate.csv", col.names = TRUE, row.names = TRUE, quote = FALSE)
  
  
  # Perform DGE analysis for the following comparisons: Severe vs. Healthy (A) and Severe vs. Mild (B).
  # The subsequent DGE analysis follows the same steps as previously described. Lines 25-95 include explanations for each step.
  {
    countsdata <- read.csv("GSE145926_expression_matrix_aggregate.csv", header = TRUE, row.names = 1)
    dgeCounts <- DGEList(countsdata)
    experimentDesign <- read.table("experiment_design_GSE145926.txt",header=T,sep="\t")
    rownames(experimentDesign) <- experimentDesign$SampleID
    experimentDesign <- experimentDesign[colnames(dgeCounts$counts),]
    geneid <- rownames(dgeCounts)
    genes <- select(org.Hs.eg.db, keys = geneid, columns = c("ENTREZID", "SYMBOL"), keytype = "ENSEMBL")
    genes <- genes[!duplicated(genes$ENSEMBL),]
    genes <- genes[match(rownames(dgeCounts$counts), genes$ENSEMBL),]
    group <- factor(experimentDesign$PatientGroup)
    dgeCounts$samples$group <- group
    dgeCounts$genes <- genes
    keptExpression <- filterByExpr(dgeCounts, group = group, min.count = 5)
    dgeCounts <- dgeCounts[keptExpression,, keep.lib.sizes=FALSE]
    dgeCounts <- calcNormFactors(dgeCounts, method = "TMM")
    fullmatrix <- model.matrix(~0+group)
    nullmatrix <- model.matrix(~1, data = dgeCounts$samples)
    svaobj <- sva(dgeCounts$counts, fullmatrix, nullmatrix)
    # Design matrix has 2 additional columns to account for 2 surrogate variables detected by SVA.
    design <- cbind(fullmatrix,svaobj$sv)
    colnames(design)<- gsub("group","",colnames(design))
    colnames(design) <- c(colnames(design)[1:3], "V3", "V4")
    contrastMatrixA <- makeContrasts(Severe_vs_Healthy = Severe - Healthy,
                                     levels = colnames(design))
    contrastMatrixB <- makeContrasts(Severe_vs_Mild = Severe - Mild,
                                     levels = colnames(design))
    y <- voom(dgeCounts, design)
    fitA <- lmFit(y, design)
    fitA <- contrasts.fit(fitA, contrastMatrixA)
    efitA <- eBayes(fitA)
    etA <- decideTests(efitA)
    coefficientVarA = colnames(fitA$coefficients)
    limmaResA_full <- topTable(efitA, coef = coefficientVarA, n = dim(fitA)[1])
    limmaResA <- limmaResA_full[which(abs(limmaResA_full$logFC) > 1 & limmaResA_full$adj.P.Val < 0.05),]
    fitB <- lmFit(y, design)
    fitB <- contrasts.fit(fitB, contrastMatrixB)
    efitB <- eBayes(fitB)
    etB <- decideTests(efitB)
    coefficientVarB = colnames(fitB$coefficients)
    limmaResB_full <- topTable(efitB, coef = coefficientVarB, n = dim(fitB)[1])
    limmaResB <- limmaResB_full[which(abs(limmaResB_full$logFC) > 1 & limmaResB_full$adj.P.Val < 0.05),]
  }

  # Perform DGE analysis for the COVID vs. Healthy (C) comparison.
  # The subsequent DGE analysis follows the same steps as previously described. Lines 25-95 include explanations for each step.
  {
    dgeCounts <- DGEList(countsdata)
    geneid <- rownames(dgeCounts)
    genes <- select(org.Hs.eg.db, keys = geneid, columns = c("ENTREZID", "SYMBOL"), keytype = "ENSEMBL")
    genes <- genes[!duplicated(genes$ENSEMBL),]
    genes <- genes[match(rownames(dgeCounts$counts), genes$ENSEMBL),]
    group <- factor(experimentDesign$Status)
    dgeCounts$samples$group <- group
    dgeCounts$genes <- genes
    keptExpression <- filterByExpr(dgeCounts, group = group, min.count = 5)
    dgeCounts <- dgeCounts[keptExpression,, keep.lib.sizes=FALSE]
    dgeCounts <- calcNormFactors(dgeCounts, method = "TMM")
    fullmatrix <- model.matrix(~0+group)
    nullmatrix <- model.matrix(~1, data = dgeCounts$samples)
    svaobj <- sva(dgeCounts$counts, fullmatrix, nullmatrix)
    # Design matrix has 2 additional columns to account for 2 surrogate variables detected by SVA.
    design <- cbind(fullmatrix,svaobj$sv)
    colnames(design)<- gsub("group","",colnames(design))
    colnames(design) <- c(colnames(design)[1:2], "V3", "V4")
    contrastMatrixC <- makeContrasts(Covid_vs_Healthy = COVID - Healthy,
                                     levels = colnames(design))
    y <- voom(dgeCounts, design)
    fitC <- lmFit(y, design)
    fitC <- contrasts.fit(fitC, contrastMatrixC)
    efitC <- eBayes(fitC)
    etC <- decideTests(efitC)
    coefficientVarC = colnames(fitC$coefficients)
    limmaResC_full <- topTable(efitC, coef = coefficientVarC, n = dim(fitC)[1])
    limmaResC <- limmaResC_full[which(abs(limmaResC_full$logFC) > 1 & limmaResC_full$adj.P.Val < 0.05),]
  }
  
  # Output DGE analysis tables as csv files.
  filenameA <- paste("~/DEgenes_voom_", coefficientVarA, "_pval0.05_logFC1.csv", sep = "")
  filenameB <- paste("~/DEgenes_voom_", coefficientVarB, "_pval0.05_logFC1.csv", sep = "")
  filenameC <- paste("~/DEgenes_voom_", coefficientVarC, "_pval0.05_logFC1.csv", sep = "")
  write.csv(limmaResA, file = filenameA, row.names = T, quote = F)
  write.csv(limmaResB, file = filenameB, row.names = T, quote = F)
  write.csv(limmaResC, file = filenameC, row.names = T, quote = F)
}

# As an exploration (not explicitly part of the RepurNet pipeline), whole blood and bronchoalveolar lavage fluid DGE profiles were compared with RRHO2 (Figure 2B).
# The code below demonstrates how this was done, as RRHO2 requires specific formatting of input values.
{
  # Read custom functions and DGE profiles obtained earlier
  {
    remove_duplicate_gene <- function(x){
      dups <- x[,"SYMBOL"][which(duplicated(x[,"SYMBOL"]))]
      remove_ensembl <- vector("character", 0)
      for(i in c(1:length(dups))){
        dups_table <- x[which(x[,"SYMBOL"] %in% dups[i]),]
        dups_table <- dups_table[order(dups_table[,"adj.P.Val"], decreasing = FALSE),]
        remove_ensembl <- append(remove_ensembl, dups_table[,"ENSEMBL"][c(2:dim(dups_table)[1])])
      }
      return(x[which(!x[,"ENSEMBL"] %in% remove_ensembl),])
    }
    reformat_table <- function(x){
      up <- x[which(x[,"logFC"] > 0), c("SYMBOL", "adj.P.Val")]
      up[,"adj.P.Val"] <- -log10(up[,"adj.P.Val"])*1
      down <- x[which(x[,"logFC"] < 0), c("SYMBOL", "adj.P.Val")]
      down[,"adj.P.Val"] <- -log10(down[,"adj.P.Val"])*(-1)
      table <- rbind(up, down)
      colnames(table)[2] <- "SCORE"
      return(table)
    }
    
    ICU <- read.csv("DEgenes_voom_CovidICU_vs_nonCovidICU_pval0.05_logFC1.csv", header = TRUE, row.names = 1)
    COVID <- read.csv("DEgenes_voom_CovidICU_vs_CovidnonICU_pval0.05_logFC1.csv", header = TRUE, row.names = 1)
    Simple <- read.csv("DEgenes_voom_COVID_vs_nonCOVID_pval0.05_logFC1.csv", header = TRUE, row.names = 1)
    Severe <- read.csv("DEgenes_voom_Severe_Vs_Control_pval0.05_logFC1.csv", header = TRUE, row.names = 1)
    SvM <- read.csv("DEgenes_voom_Severe_Vs_Mild_pval0.05_logFC1.csv", header = TRUE, row.names = 1)
    SimpleBALF <- read.csv("DEgenes_voom_Infected_Vs_Healthy_pval0.05_logFC1.csv", header = TRUE, row.names = 1)
  }
  
  # Perform COVID ICU vs. nonCOVID ICU (blood) to Severe vs. Healthy (BALF) comparison (A)
  ICU <- ICU[which(rownames(ICU) %in% Severe$SYMBOL),]
  Severe <- Severe[which(Severe$SYMBOL %in% rownames(ICU)),]
  if(any(duplicated(Severe$SYMBOL))) Severe <- remove_duplicate_gene(Severe)
  ICU_Table <- reformat_table(ICU)
  Severe_Table <- reformat_table(Severe)
  ICU_RRHO2_pValue <- RRHO2_initialize(ICU_Table, Severe_Table, labels = c("COVID ICU vs. nonCOVID ICU", "Severe vs. Healthy"), 
                                       log10.ind = TRUE, method = "hyper", multipleTesting = "none", boundary = 0.05, stepsize = 1)
  
  # Perform COVID ICU vs. COVID nonICU (blood) to Severe vs. Mild (BALF) comparison (B)
  COVID <- COVID[which(rownames(COVID) %in% SvM$SYMBOL),]
  SvM <- SvM[which(SvM$SYMBOL %in% rownames(COVID)),]
  if(any(duplicated(SvM$SYMBOL))) SvM <- remove_duplicate_gene(SvM)
  COVID_Table <- reformat_table(COVID)
  SvM_Table <- reformat_table(SvM)
  COVID_RRHO2_pValue <- RRHO2_initialize(COVID_Table, SvM_Table, labels = c("COVID ICU vs. COVID nonICU", "Severe vs. Mild"), 
                                         log10.ind = TRUE, method = "hyper", multipleTesting = "none", boundary = 0.05, stepsize = 1)
  
  # Perform COVID vs. nonCOVID (blood) to COVID vs. Healthy (BALF) comparison (C)
  Simple <- Simple[which(rownames(Simple) %in% SimpleBALF$SYMBOL),]
  SimpleBALF <- SimpleBALF[which(SimpleBALF$SYMBOL %in% rownames(Simple)),]
  if(any(duplicated(SimpleBALF$SYMBOL))) SimpleBALF <- remove_duplicate_gene(SimpleBALF)
  Simple_Table <- reformat_table(Simple)
  SimpleBALF_Table <- reformat_table(SimpleBALF)
  Simple_RRHO2_pValue <- RRHO2_initialize(Simple_Table, SimpleBALF_Table, labels = c("COVID vs. nonCOVID", "COVID vs. Healthy"), 
                                          log10.ind = TRUE, method = "hyper", multipleTesting = "none", boundary = 0.05, stepsize = 1)
  
  # Make RRHO2 Heatmaps for each RRHO2 comparison
  png(filename = "figure2b_A.png", 
      width = 10,
      height = 10,
      units = "in",
      res = 300)
  RRHO2_heatmap(ICU_RRHO2_pValue,
                maximum = 3)
  dev.off()
  png(filename = "figure2b_B.png", 
      width = 10,
      height = 10,
      units = "in",
      res = 300)
  RRHO2_heatmap(COVID_RRHO2_pValue,
                maximum = 2)
  dev.off()
  png(filename = "figure2b_C.png", 
      width = 10,
      height = 10,
      units = "in",
      res = 300)
  RRHO2_heatmap(Simple_RRHO2_pValue,
                maximum = 4)
  
  dev.off()
}


