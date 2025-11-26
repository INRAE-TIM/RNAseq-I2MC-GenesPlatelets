# Formating sample file ----
samples <- data.frame(group = pData$group, rep_name = row.names(pData))

# is there duplicated genes names ?
if(isTRUE(sum(Counts$gene_name == unique(Counts$gene_name)) == length(Counts$gene_name))){
  cat('no duplicated gene_name')
} else{cat('duplicated gene_name exists')}

row.names(Counts) <- Counts$gene_name
# Formating counts file ----

d <- DGEList(counts = Counts, group = factor(samples$group))

# Filtering out genes with less than 1 cpm
keep <- rowSums(cpm(d) > 1) >= 2
d <- d[keep, ]
dim(d)

d$samples$lib.size <- colSums(d$counts)

# Intra-individual normalization ----
d <- calcNormFactors(d)

# Data exploration ----
plotMDS(d,
        method = "logFC",
        col = as.numeric(d$samples$group),
        main = "MDS Plot of Gene Expression Data")
legend("topleft",
       as.character(unique(d$samples$group)),
       col = 1:4,
       pch = 20)

# estimating overdispertion ----
#plotBCV(estimateTagwiseDisp(estimateCommonDisp(d, verbose = T)))

# GLM
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d, design.mat)
d2 <- estimateGLMTrendedDisp(d2, design.mat, method = "bin.spline")
d2 <- estimateGLMTagwiseDisp(d2, design.mat)
plotBCV(d2)
title(main = "BCV Plot with GLM Dispersion Estimates")
# Biological Coefficient of Variation

fit <- glmFit(d2, design.mat, prior.count	= 1)

# Comparisons table
adjResTable <- make_comparisons(fit, 
                                list(c(1, 0, 0, -1, 0, 0), 
                                     c(0, 1, 0, 0, -1, 0), 
                                     c(0, 0, 1, 0, 0, -1), 
                                     c(1, -1, 0, 0, 0, 0),
                                     c(1, 0, -1, 0, 0, 0),
                                     c(0, 0, 0, 1, -1, 0),
                                     c(0, 0, 0, 1, 0, -1),
                                     c(0, 1, -1, 0, 0, 0),
                                     c(0, 0, 0, 0, 1, -1)
                                )
)

# Formating to the MATRiX app format
colnames(adjResTable) <- c(
  gsub('\\.logFC','', 
       gsub('\\.logCPM','',
            gsub('\\.LR', '',
                 gsub('\\.PValue','',
                      paste(c('logFC', 'logCPM', 'LR', 'p.value'), 
                            colnames(adjResTable)[!grepl('adj.P.Val', colnames(adjResTable))], sep = '_'))))),
  
  gsub('\\.adj.P.Val','',
       paste(c('FDR'), 
             colnames(adjResTable)[grepl('adj.P.Val', colnames(adjResTable))], sep = '_'))
  
)
adjResTable <- cbind(
  gene_id = Counts$gene_id[Counts$gene_name %in% row.names(adjResTable)],
  gene_name = row.names(adjResTable), 
  adjResTable)

row.names(adjResTable) <- 1:dim(adjResTable)[1]

RNAseq <- list(ResTable = adjResTable,
               pData = data.frame(group = pData$group, row.names = row.names(pData)),
               WS = as.data.frame(cpm(d2, log = T)))
rm(list = c('Counts', 'd', 'd2', 'design.mat', 'fit', 'samples', 'keep'))
