# Automatically load some specific files from a directory
data_loader <- function(keyword)
{
  if(
    grepl('.txt',
          list.files('./Data')[grepl(keyword, list.files('./Data'))][[1]])
  )
  {
    dataList <- lapply(paste('./Data',
                             list.files('./Data')[grepl(keyword, list.files('./Data'))],
                             sep = '/'),
                       function(x){read.table(x, header = T)})
  } else if(
    grepl('.csv',
          list.files('./Data')[grepl(keyword, list.files('./Data'))][[1]])
  )
  {
    dataList <- lapply(paste('./Data',
                             list.files('./Data')[grepl(keyword, list.files('./Data'))],
                             sep = '/'),
                       function(x){read.csv2(x, header = T)})
  }
  
  dataList <- setNames(dataList,
           sapply(list.files('./Data')[grepl(keyword, list.files('./Data'))],
                  function(x)
                  {
                    paste(strsplit(x, '_')[[1]][7:8], collapse = '_')
                  },
                  USE.NAMES = F)
                      )
  
  dataList
}

# Filter a count table to avoid transcripts with 0 counts (default)
# use minCount option to select a minimal number of transcripts.
filter_count_table <- function(count_table,
                             minCount = NULL)
{
  Res <- count_table[!rowSums(count_table[,-1])==0,]
  if(is.null(minCount))
  {
    Res <- Res
  }else{
    Res <- Res[apply(Res[,-1],1,max)>minCount,]
  }
  
  row.names(Res) <- Res$Gene_id
  
  Res[,-1]
}

make_comparisons <- function(fit,
                            design)
# makes the comparison using EdgeR's glmLRT function.
# P-values adjustment is method = 'BH' to match EdgeR's results

{
  lrtList <- lapply(design, function(x) {
    glmLRT(fit, contrast = x)
  })
  names(lrtList) <- design |>
    lapply(function(x) {
      colnames(as.data.frame(fit$design))[x != 0]
    }) |>
    sapply(function(x) {
      paste0(x, collapse = '_vs_')
    })
  
  adj.P.Val <- lapply(lrtList, function(x){p.adjust(x$table$PValue, 'BH')}) |> 
    do.call(what = cbind.data.frame)
  
  colnames(adj.P.Val) <- paste(colnames(adj.P.Val), 'adj.P.Val', sep = '.')
  
  lrtList <- lapply(lrtList, function(x){x[['table']]}) |> 
    do.call(what = cbind.data.frame)
  
  row.names(adj.P.Val) <- row.names(lrtList)

  lrtList |> cbind.data.frame(adj.P.Val) 
  
}

do.lm <- function(gene_name)
{
  logCPM = as.numeric(RNAseq$WS[gene_name, ])
  TP = sqrt(setNames(MetaData$TP, MetaData$Diet_Treatment.Rep)[colnames(RNAseq$WS)])
  diets = c('HFC', 'ND')
  results = lapply(diets, function(x) {
    PC = TP[grepl(x, names(TP))]
    log_CPM = logCPM[grepl(x, names(TP))]
    lm(log_CPM ~ PC)
  })
  names(results) <- diets
  results
}

summarize_lm <- function(x){
  LMlist <- lapply(x, summary)
  
  mydf <- cbind.data.frame(
    data.frame(Diet = rep(names(LMlist), each = 2),
               lm.param = rep(c('intercept', 'slope'), length(LMlist))),
    do.call(rbind, lapply(LMlist, function(x) coef(x)))
  )
  
  row.names(mydf)<- 1:dim(mydf)[1]
  
  mydf
}

plot_venn <- function(data,
                     comps,
                     adj.P.Val = 0.05,
                     FC = 1)
  # a function to plot venn diagrams from a dataframe.
  # adjusted p-values shoud be in columns names as follow :
  # adj.P.Val_name_your_comparison
  # logFC columns should be named as follow :
  # logFC_name_you_comparison
{
  dataVenn = sapply(comps,
                    function(x){
                      ifelse(
                        as.vector(data[, paste(x, 'adj.P.Val', sep = '.')] < adj.P.Val) &
                          as.vector(abs(data[, paste(x, 'logFC', sep = '.')]) > log2(FC)),
                        1,
                        0
                      )
                    }) |> as.data.frame()
  
  colnames(dataVenn) = comps
  
  venn::venn(dataVenn)
}


library(pheatmap)

plotHM <- function(data,
                   pData,
                   ResTable,
                   adj.P.Val_threshold = .05,
                   HMcolors = c("blue", "white", "red"),
                   nclust,
                   scale = 'row',
                   group_colors = NULL,
                   show_rownames = FALSE,
                   show_colnames = FALSE,
                   fontsize_col = 3,
                   legend = FALSE,
                   annotation_legend = FALSE)
{
  
  dataHM <- t(data[,-1])
  colnames(dataHM) <- paste(data$group, row.names(data)) 
  
  dataHM <- dataHM[row.names(dataHM) %in%
                     row.names(ResTable)[rowSums(ResTable[,grepl('adj', colnames(ResTable))]<.05)>1],]
  
  if(is.null(group_colors)){group_colors = scales::viridis_pal()(length(unique(pData$Grp)))
  names(group_colors) <- unique(pData$Grp)}
  else{group_colors = group_colors
  names(group_colors) <- unique(pData$Grp)}
  
  # Pair correlate columns (samples)
  cols.cor <- cor(dataHM, use = "complete.obs")
  # Pair correlate  lines (genes)
  rows.cor <- cor(t(dataHM), use = "complete.obs")
  
  # first, do HM without row clustering ----
  pres <- pheatmap(dataHM,
                   scale = scale, 
                   color = colorRampPalette(HMcolors)(50),
                   clustering_distance_cols = as.dist(1 - cols.cor),
                   clustering_distance_rows =as.dist(1 - rows.cor),
                   clustering_method = "ward.D2",
                   show_rownames = FALSE,
                   silent = T)
  
  # Then, select the number of clusters ----
  
  choose_n_clust = nclust
  my_gene_col <- data.frame(cluster = as.numeric(cutree(tree = pres$tree_row,
                                                        k = choose_n_clust))) 
  my_gene_col$num <- 1:length(pres$tree_row$labels)
  row.names(my_gene_col) <- row.names(dataHM)
  my_gene_col <- my_gene_col[pres$tree_row$order,]
  my_gene_col$cluster = factor(my_gene_col$cluster, levels = unique(my_gene_col$cluster), ordered = T)
  my_gene_col$cluster = as.numeric(my_gene_col$cluster)

  my_sample_col <- data.frame(sample = pData$Grp)
  row.names(my_sample_col) <- colnames(dataHM)
  
  clustcol <- list(sample = group_colors,
                   cluster = rep(c("grey20","grey80"),
                                 (choose_n_clust/2)+1)[1:choose_n_clust])
  names(clustcol$cluster) <- paste("cluster", 1:choose_n_clust, sep = "")
  
  pHM<- pheatmap(dataHM,
                 scale = scale, 
                 color = colorRampPalette(HMcolors)(50),
                 clustering_distance_cols = as.dist(1 - cols.cor),
                 clustering_distance_rows =as.dist(1 - rows.cor),
                 clustering_method = "ward.D2",
                 annotation_colors = clustcol,
                 annotation_col = my_sample_col,
                 annotation_row = my_gene_col |>
                   dplyr::mutate(cluster = paste("cluster", cluster, sep = "")) |> 
                   dplyr::select(cluster),
                 show_rownames = show_rownames,
                 show_colnames = show_colnames,
                 legend = legend,
                 annotation_legend = annotation_legend,
                 treeheight_row = 10,
                 treeheight_col = 10,
                 fontisize = 3,
                 fontsize_col = fontsize_col,
                 silent = T)
  ResHM <<- my_gene_col
  return(ggpubr::ggarrange(pHM[[4]]))
  
  
}  

require(fgsea)
require(ViSEAGO)

GivMeGOid <- function(GO_names)
{
  sapply(GO_names, function(x){
    GOfuncR::get_ids(x)$go_id[
      GOfuncR::get_ids(x)$node_name == x]})
}


GSEA_plot <- function(
    ResTable, GO_name, comp)
{
  scoring <- setNames(-log10(ResTable[,paste(comp, 'adj.P.Val', sep = '.')])*
                        ResTable[,paste(comp, 'logFC', sep = '.')]/
                        abs(ResTable[,paste(comp, 'logFC', sep = '.')]),
                      row.names(ResTable))
  
  GenesTable = GOfuncR::get_anno_genes(
    database = "Mus.musculus",
    go_ids = GivMeGOid(GO_name))
  colnames(GenesTable) <- c('go_id','GeneName')
  GenesList <- lapply(unique(GenesTable$go_id),
                      function(x){
                        GenesTable$GeneName[GenesTable$go_id == x]
                      }) |> 
    setNames(GOfuncR::get_names(unique(GenesTable$go_id))$go_name)
  ggpubr::ggarrange(plotlist =  lapply(GenesList,
         function(x){
           fgsea::plotEnrichment(x, scoring)
         }))
}

multilevel_GSEA <- function(ResTable, 
                            comp,
                            ont,
                            score = 'signed adj.p.val',
                            nperm = 10000,
                            sampleSize = 101,
                            minSize = 1,
                            maxSize = Inf,
                            eps = 0,
                            scoreType = "std",
                            nproc = 0,
                            gseaParam = 1,
                            BPPARAM = NULL,
                            absEps= NULL)
{
  Bioconductor <- ViSEAGO::Bioconductor2GO()
  myGENE2GO <- ViSEAGO::annotate(
    id="org.Mm.eg.db",
    object=Bioconductor
  )
  
  SelGene <- as.data.frame(
    mapIds(org.Mm.eg.db,
           keys = row.names(ResTable),
           column='ENTREZID',
           keytype='SYMBOL',
           multiVals='first')
  )
  SelGene$GeneName <- row.names(SelGene)
  colnames(SelGene) <- c('Id','GeneName')
  SelGene <- SelGene[!is.na(SelGene$Id),]
  logFC = ResTable[row.names(ResTable) %in%
                        SelGene$GeneName,
                      paste(comp, 'logFC', sep = '.')]
  adj.P = ResTable[row.names(ResTable) %in%
                        SelGene$GeneName,
                      paste(comp, 'adj.P.Val', sep = '.')]
  if(score == 'signed adj.p.val'){
    SelGene$SCORE = (logFC/abs(logFC))*
      -log10(adj.P)
  } else if(score == 'logFC'){
    SelGene$SCORE = logFC
  }else if(score == 'adj.p'){
    SelGene$SCORE = adj.P
  }
  
  SelGene <- data.table::data.table(SelGene)[,c('Id','SCORE')]
  
  GSEAresults <- runfgsea(
    SelGene,
    myGENE2GO,
    ont = ont,
    method = "fgseaMultilevel",
    params = list(nperm = nperm, sampleSize = sampleSize, minSize = minSize, maxSize = maxSize, eps = eps,
                  scoreType = scoreType, nproc = nproc, gseaParam = gseaParam, BPPARAM = BPPARAM, absEps
                  = absEps)
  )
  
  return(GSEAresults)
}

show_GSEAres <- function(GSEAres)
{
  GSEAres@data[[1]] |> 
    tidyr::unnest(leadingEdge) |> 
    dplyr::mutate(GeneNames = as.character(mapIds(org.Mm.eg.db,
                                           leadingEdge,
                                           keytype = 'ENTREZID',
                                           column = 'SYMBOL'))) |> 
    dplyr::select(-leadingEdge) |>
    dplyr::group_by(pathway, pval, padj, log2err, ES, NES, size, genes_frequency) |> 
    dplyr::summarise(GeneNames = list(GeneNames)) |> 
    dplyr::ungroup() |> 
    dplyr::mutate(go_name = GOfuncR::get_names(pathway)$go_name)
  
}

GSEA_top_tiles <- function(GSEAres){
  GSEAres@data[[1]]|> 
    filter(padj<.05) |> 
    mutate(GO_name = GOfuncR::get_names(pathway)$go_name) |>
    mutate(GO_name = fct_reorder(GO_name, -log10(padj))) |> 
    ggplot()+
    geom_tile(aes(x = 1, y = GO_name, fill = -log10(padj)))+
    scale_fill_gradient(low = 'grey',high = 'yellow', limits = c(0,NA))+
    theme_classic()+
    theme(legend.position = 'top',
          axis.title = element_blank(),
          axis.line = element_blank(),
          axis.ticks = element_blank(),
          text = element_text(size = 8),
          legend.text = element_text(size = 8),
          axis.text.x = element_blank()
    )
}

leading_genes <- function(GSEAres){
  do.call(rbind,
          setNames(
            lapply(as.list((
              GSEAres@data[[1]] |> filter(padj < .05)
            )$leadingEdge),
            function(x) {
              AnnotationDbi::select(
                org.Mm.eg.db,
                keys = x,
                columns = c('ENTREZID', 'SYMBOL'),
                keytypes = 'ENTREZID'
              )
            }),
            str_replace_all(GOfuncR::get_names((
              GSEAres@data[[1]] |> filter(padj < .05)
            )$pathway)$go_name, ' ', '_')
          )) |>
    rownames_to_column('go_name') |>
    mutate(go_name = str_replace_all(str_remove_all(
      str_remove_all(go_name, '[:digit:]'), '\\.'
    ),
    '_', ' ')) |>
    dplyr::rename(GeneName = SYMBOL)
}

dual_volcano <- function(ResTable,
                         comp1,
                         comp2,
                         couleurs){
  signedPvals <-
    cbind.data.frame(-log10(ResTable[, paste(comp1, 'adj.P.Val', sep = '.')]) *
                       abs(ResTable[, paste(comp1, 'logFC', sep = '.')]) /
                       ResTable[, paste(comp1, 'logFC', sep = '.')],
                     -log10(ResTable[, paste(comp2, 'adj.P.Val', sep = '.')]) *
                       abs(ResTable[, paste(comp2, 'logFC', sep = '.')]) /
                       ResTable[, paste(comp2, 'logFC', sep = '.')]) |>
    setNames(c(comp1, comp2))
  
  signedPvals$Ledgend <- ifelse(
    abs(signedPvals[,comp1])>(-log10(.05))&abs(signedPvals[,comp2])<(-log10(.05)),comp1,
    ifelse(abs(signedPvals[,comp1])<(-log10(.05))&abs(signedPvals[,comp2])>(-log10(.05)),comp2,
           ifelse(abs(signedPvals[,comp1])>(-log10(.05))&abs(signedPvals[,comp2])>(-log10(.05)),'both',
                  'none')))
  
  signedPvals$GeneName <- row.names(ResTable)
  
  signedPvals |> 
    ggplot(aes(text = GeneName))+
    geom_point(aes(x = !!sym(comp1),
                   y = !!sym(comp2),
                   color = Ledgend))+
    geom_vline(xintercept = -log10(.05), linetype = 'dotdash')+  
    geom_vline(xintercept = log10(.05), linetype = 'dotdash')+
    geom_hline(yintercept = -log10(.05), linetype = 'dotdash')+
    geom_hline(yintercept = log10(.05), linetype = 'dotdash')+
    xlab(gsub('_',' ',comp1))+
    ylab(gsub('_',' ',comp2))+
    theme_classic()+
    scale_color_manual(values = couleurs)
}


ClusEnrich <- function(GO_ids,HMresults,GenesBackground)
{
  
  clustList = unique(HMresults$cluster) |> 
    as.list() |> 
    setNames(as.numeric(unique(HMresults$cluster)))
  
  GO_genes <- lapply(
    GO_ids,
    function(x){GOfuncR::get_anno_genes(x, 'Mus.musculus')}) |> 
    setNames(GOfuncR::get_names(GO_ids)$go_name)
  
  t(do.call(rbind.data.frame,
            lapply(clustList, function(y){
              lapply(GO_genes, function(x){
                phyper(
                  length(x$gene[x$gene %in% row.names(HMresults)[ResHM$cluster == y]]),
                  length(x$gene[x$gene %in% GenesBackground]),
                  length(GenesBackground)-
                    length(x$gene[x$gene %in% GenesBackground]),
                  length(row.names(ResHM)[ResHM$cluster == y]), lower.tail = F
                )
              })}))) |> 
    as.data.frame()|> 
    rownames_to_column('term')|> 
    pivot_longer(!term, names_to = 'cluster', values_to = 'pvalue')  |> 
    left_join(do.call(rbind.data.frame,GO_genes) |> 
                rownames_to_column('term') |> 
                separate(term, c('term', 'discard'), sep = '\\.') |> 
                dplyr::select(-discard) |> 
                group_by(term) |> 
                left_join(HMresults |> 
                            rownames_to_column('gene') |> 
                            dplyr::select(-num)) |> 
                drop_na() |> 
                group_by(term, cluster)|> 
                mutate(cluster = as.character(cluster)) |> 
                summarize(nGenes = n()) |> 
                mutate(term = str_replace_all(term, ' ', '.')))
  
  
}


my_GSEA_plot <- function(ResTable,
                         comparison,
                         GSEA_results,
                         go_id,
                         score)
{
  SelGene <- as.data.frame(
    mapIds(org.Mm.eg.db,
           keys = row.names(ResTable),
           column='ENTREZID',
           keytype='SYMBOL',
           multiVals='first')
  )
  SelGene$GeneName <- row.names(SelGene)
  colnames(SelGene) <- c('Id','GeneName')
  SelGene <- SelGene[!is.na(SelGene$Id),]
  logFC = ResTable[row.names(ResTable) %in%
                     SelGene$GeneName,
                   comparison]
  adj.P = ResTable[row.names(ResTable) %in%
                     SelGene$GeneName,
                   comparison]
  if(score == 'signed adj.p.val'){
    SelGene$SCORE = (logFC/abs(logFC))*
      -log10(adj.P)
  } else if(score == 'logFC'){
    SelGene$SCORE = logFC
  }else if(score == 'adj.p'){
    SelGene$SCORE = adj.P
  }
  
  SelGene <- data.table::data.table(SelGene)
  
  GOgenes = topGO::inverseList(slot(myGENE2GO,'BP'))[go_id][[1]]
  forEnrich <- SelGene$GeneName[SelGene$Id %in% GOgenes[GOgenes %in% SelGene$Id]]
  
  p <-
    plotEnrichment(
      forEnrich,
      setNames(ResTable[, comparison],
               row.names(ResTable))
    )
  
  NES = gsub('_vs','', gsub('.logFC','.NES', comparison))
  FDR = gsub('_vs','', gsub('.logFC','.FDR_pvalue', comparison))
  qvalue = gsub('_vs','', gsub('.logFC','.qvalue', comparison))
  
  
  p <-ggplot(data = p$data)+
    geom_hline(yintercept = 0, color = 'grey') +
    geom_line(aes(x = rank, y = ES),
              color = '#6ab533') +
    geom_rect(
      aes(
        xmin = -Inf,
        xmax = Inf,
        ymax = min(p$data$ES),
        ymin = min(p$data$ES) - (max(p$data$ES) - min(p$data$ES)) /
          4
      ),
      fill = 'white',
      color = 'grey90'
    ) +
    geom_linerange(aes(
      x = p$data$rank,
      ymin = min(p$data$ES) - (max(p$data$ES) - min(p$data$ES)) / 4,
      ymax = min(p$data$ES),
      width = .1/dim(ResTable)[1]
    )) +
    geom_rect(
      data = rbind.data.frame(
        data.frame(valeur = sort(ResTable[, comparison][ResTable[, comparison] >= 0], decreasing = T)) |>
          mutate(quartile = ifelse(
            valeur < max(valeur) / 80,
            'A',
            ifelse(
              valeur < max(valeur) / 30,
              'B',
              ifelse(valeur < max(valeur) / 10, 'C',
                     ifelse(valeur < max(valeur) /
                              5, 'D',
                            'E'))
            )
          )),
        data.frame(valeur = sort(ResTable[, comparison][ResTable[, comparison] < 0], decreasing = T)) |>
          mutate(quartile = ifelse(
            abs(valeur) < max(abs(valeur)) / 80,
            'A',
            ifelse(
              abs(valeur) < max(abs(valeur)) / 30,
              'B',
              ifelse(abs(valeur) < max(abs(valeur)) /
                       10, 'C',
                     ifelse(abs(valeur) < max(abs(
                       valeur
                     )) / 5, 'D',
                     'E'))
            )
          ))
      ) |>
        mutate(
          rank = 1:length(valeur),
          posneg = factor(
            ifelse(valeur > 0, 'pos', 'neg'),
            levels = c('pos', 'neg'),
            ordered = T
          )
        ) |>
        group_by(posneg, quartile) |>
        summarize(xmin = min(rank), xmax = max(rank)),
      aes(
        xmin = xmin,
        xmax = xmax,
        ymax = min(p$data$ES) - (max(p$data$ES) - min(p$data$ES)) / 4,
        ymin = (min(p$data$ES) - (max(p$data$ES) - min(p$data$ES)) / 4) -
          abs((min(p$data$ES) - (max(p$data$ES) - min(p$data$ES))) / 30),
        alpha = quartile,
        fill = posneg
      )
    ) +
    scale_fill_manual(values = c('red', 'blue')) +
    ggtitle(BP_sResults@data[BP_sResults@data$GO.ID == go_id, 'term']) +
    ylab('Enrichment Score') +
    xlab(comparison) +
    theme_bw() +
    theme(
      plot.background = element_rect(fill = "grey90"),
      panel.border = element_rect(colour = 'grey'),
      legend.position = 'none',
      text = element_text(size = 8)
    ) +
    geom_text(data = data.frame(
      label = paste(c('NES =', 'FDR pvalue ='),
                    round(unlist(
                      GSEA_results[GSEA_results$GO.ID == go_id,
                                        c(NES, FDR)]
                    ),
                    digits = 3),
                    collapse = '\n'),
      x = as.numeric(ifelse(
        GSEA_results[GSEA_results$GO.ID == go_id,
                          NES] > 0,
        dim(adjResTable)[1] /
          7,
        dim(adjResTable)[1] /
          1.5
      )),
      y = as.numeric(
        ifelse(
          GSEA_results[GSEA_results$GO.ID == go_id,
                            NES] > 0,
          min(p$data$ES)[1] + max(p$data$ES)[1] /
            4,
          max(p$data$ES)[1] - 0.2*(max(p$data$ES)[1] -min(p$data$ES)[1])
        )
      )
    ),
    aes(label = label, x = x, y = y),
    size = 2.1)
  
  ggpubr::ggarrange(plotlist = list(p))
}

calc_qvals <- function(GSEA_res)
{
  GSEA_name = paste(gsub('GSEA_', '', as.character(enexpr(GSEA_res))),'.',sep ='')
  
  results <- cbind.data.frame(
    GO.ID = GSEA_res@data[[1]]$pathway,
    FDR_pvalue = qvalue(GSEA_res@data[[1]]$pval)$lfdr,
    qvalue = qvalue(GSEA_res@data[[1]]$pval)$qvalues
  )
  
  names(results) <- paste(c('',rep(GSEA_name,2)),
                          names(results),
                          sep = '')
  results
  
}

# Formatting a data frame for MATRiX input
# Full regex pattern is : "^[A-Z]{3}\\.[A-Za-z0-9]{3}_vs_[A-Z]{3}\\.[A-Za-z0-9]{3}\\.logFC$" for `adjResTable`.
colnames_format.MATRiX <- function(data_frame) {
  sapply(colnames(data_frame), function(colname) {
    if (grepl("\\.logFC$", colname)) {
      new_colname <- paste("logFC", gsub('.logFC', '', colname), sep = '_')
      return(new_colname)
    }
    else if (grepl("\\.PValue$", colname)) {
      new_colname <- paste("logFC", gsub('.PValue', '', colname), sep = '_')
      return(new_colname)
    }
    else if (grepl("\\.logCPM$", colname)) {
      new_colname <- paste("logFC", gsub('.logCPM', '', colname), sep = '_')
      return(new_colname)
    }
    else if (grepl("\\.LR$", colname)) {
      new_colname <- paste("logFC", gsub('.LR', '', colname), sep = '_')
      return(new_colname)
    }
    else if (grepl("\\.adj.P.Val$", colname)) {
      parts <- unlist(strsplit(colname, "\\."))
      new_colname <- paste("FDR", gsub('.adj.P.Val', '', colname), sep = '_')
      return(new_colname)
    }
    else {
      return(colname)
    }})}

format_to_PCA <- function(RNAseq)
{
  dataPCA = cbind.data.frame(
    Grp = RNAseq$pData$group,
    t(RNAseq$WS[,!colnames(RNAseq$WS) %in% c('gene_name', 'gene_id')]))
  row.names(dataPCA) = row.names(RNAseq$pData)
  colnames(dataPCA) = c('Grp', row.names(RNAseq$WS))
  dataPCA
}

get_top_genes <- function(HMclust,
                          nTopGenes)
{
  data <- HMclust$ResTable
  sumLogFC <- rowSums(abs(data[, grepl('logFC', colnames(data))]))
  sumFDR <- rowSums(-log10(data[, grepl('FDR', colnames(data))]))
  
  scores <- data.frame(
    cluster = data$cluster,
    gene_name = data$gene_name,
    score = sumLogFC + sumFDR
  )
  
  res <- lapply(1:max(data$cluster), function(x)
  {
    A = scores[scores$cluster == x, ]
    B = A[order(-A$score),]
    B[1:nTopGenes,]
  })
  names(res) <- paste("cluster", 1:max(data$cluster), sep = '_')
  res
}

customGO_ClusPlot <- function(res_enrich, cutoff, lower_y, upper_y)
{
  res_enrich |>
    filter(p.adj < cutoff) |> 
    mutate(go_name = fct_reorder(go_name, -log10(p.adj))) |> 
    group_by(cluster) |> 
    arrange(-log10(p.adj)) |> 
    mutate(ypos = seq_along(go_name),
           cluster = gsub('_', ' ', cluster)) |> 
    ggplot()+
    geom_rect(aes(xmin = 0, 
                  xmax = -log10(p.adj),
                  ymin = ypos-lower_y,
                  ymax = ypos+upper_y,
                  fill = -log10(p.adj)))+
    geom_text(aes(x = 0,
                  y = ypos,
                  label = go_name),
              hjust = 0,
              size = 8,
              size.unit = 'pt')+
    xlab(expression(-log[10](p.adj)))+
    theme_classic()+
    theme(axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          legend.position = 'none',
          text = element_text(size = 10),
          strip.background = element_blank()
    )+
    scale_fill_viridis_c(option = 'inferno', begin = 1, end = 0.5)+
    facet_grid(rows = vars(cluster), space = 'free_y', scale = 'free_y')
}
