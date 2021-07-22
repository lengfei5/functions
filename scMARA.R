##########################################################################
##########################################################################
# Project: MS lineage embryogenesis in C elegans
# Script purpose: predicting TF regulation activity 
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Oct  7 17:31:46 2020
##########################################################################
##########################################################################
predict.TF.MARA.for.scdata = function(sub.obj, mode = c('cluster.based', 'time.bin', 'cell.based'), 
                                      id = 'manual.annot.ids', process.motif.oc = FALSE, 
                                      Y_name = 'RNA')
{
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
  library(Seurat)
  library(scater)
  library(SingleCellExperiment)
  library(scran)
  require(glmnet)
  source.my.script('scMARA_utility_functions.R')
  
  ##########################################
  # import processed tables : motif-tf mapping, motif.oc matrix,
  # known TF annotation
  ##########################################
  ll = readRDS(file = '../data/motifs_tfs/ce11_proteinCoding_genes_geneLength_transcriptLength.rds') # gene length and transript length
  #motif.oc = readRDS(file = '../data/motifs_tfs/motif_oc_all_proteinCodingGenes.rds')
  motif.oc = readRDS(file = '../data/motifs_tfs/motif_oc_scATACpeaks_all_proteinCodingGenes.rds')
  
  # modify this drosophila motif's name, because it is not for nhr-67 
  colnames(motif.oc)[which(colnames(motif.oc) == 'nhr-67.M141')] = 'M1471_1.02_M141' 
  motif.tf = readRDS(file = '../data/motifs_tfs/motif_tf_mapping.rds')
  tfs = readxl::read_xlsx('../data/motifs_tfs/Table-S2-wTF-3.0-Fuxman-Bass-Mol-Sys-Biol-2016.xlsx', sheet = 1)
  tf.mat = readRDS(file = paste0(RdataDir, 'TFs_expression_profiles_BWM.rds')) 
  
  ##########################################
  # normalized the single cell gene expression matrix with gene length 
  ##########################################
  ids = sub.obj$manual.annot.ids
  ids.uniq = unique(ids)
  ids.uniq = ids.uniq[order(ids.uniq)]
  ids.uniq = ids.uniq[order(nchar(ids.uniq))]
  
  ids.uniq = ids.uniq[grep('mixture_terminal_', ids.uniq, invert = TRUE)]
  
  # convert to SingleCellExperiment, recalculate scaling factor and normalized to fpkm
  sce = as.SingleCellExperiment(sub.obj)
  qclust <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters = qclust)
  sce <- logNormCounts(sce, log = TRUE, pseudo_count = 1)
  
  mm = match(rownames(sce), ll$gene.name)
  sce = sce[which(!is.na(mm)), ] # keep genes with correponding lengths
  transcript.length  = ll$transcript.length[mm[which(!is.na(mm) == TRUE)]]
  
  sce <- logNormCounts(sce, log = FALSE, size_factors = NULL)
  Y.fpkm <- log2(calculateFPKM(sce, lengths = transcript.length) + 1)
  
  remove(sce)
  
  ########################################################
  ########################################################
  # with manually annoted ids, there are seveval options to consider trajectories
  # 1) collapse cells in the same cell states to have pseudo-bulk
  # 2) pseudotime for trajectories, better resolution than pseud-bulk
  # 3) retain single cell 
  ########################################################
  ########################################################
  trajectory.resolution = 'pseudotime'
  
  if(trajectory.resolution == 'cluster.based'){
    ##########################################
    # average cells with the same ids (i.e. cluster-based motif activity )
    ##########################################
    Y.mat = matrix(NA, nrow = nrow(Y.fpkm), ncol = length(ids.uniq))
    colnames(Y.mat) = ids.uniq
    rownames(Y.mat) = rownames(Y.fpkm)
    for(n in 1:length(ids.uniq))
    {
      cat(ids.uniq[n], '\n')
      jj = which(ids == ids.uniq[n])
      if(length(jj) == 1) Y.mat[, n] = Y.fpkm[,jj]
      if(length(jj) > 1) Y.mat[, n] = apply(Y.fpkm[,jj], 1, mean)
    }
    
    # process.detected.tf.expression.profiles(Y.mat)
    remove(Y.fpkm) # free memory
    
    ##########################################
    # determine lineage-specific signatures 
    # i.e. selected gene sets (lineage-wide, specific, or restricted)
    # here we used markers from scRNA analysis
    ##########################################
    lineage.list = list(c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx'),
                        c('MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
    )
    
    #lineage.list = list(c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx'))
    #lineage.list = list(c('MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp'))
    lineage.list = list(c('MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx', 'MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp'))
    
    cat(length(lineage.list), ' lineages or trajectories to consider \n')
    
    trajectory.genes.with.markers = FALSE
    if(trajectory.genes.with.markers){
      # gene.sels = define.modules.for.lineags(sub.obj, Y.fpkm, lineage = lineage)
      markers = readRDS(file = paste0(RdataDir,  'AllMarkers_MST_manual.annotation.rds'))
      markers.sels = markers[which(markers$p_val<10^-3 & markers$avg_logFC > 1), ]
      cat(length(unique(markers.sels$gene)), ' marker genes selected \n')
      print(table(markers.sels$cluster))
      
      
      # prepare input matrix for tradeSeq
      ids.sels = c()
      for(m in 1:length(lineage.list)) ids.sels = unique(c(ids.sels, lineage.list[[m]]))
      
      gene.sels = markers.sels[!is.na(match(markers.sels$cluster, ids.sels)), ]
      gene.sels = gene.sels[!is.na(match(gene.sels$gene, rownames(Y.mat))), ]
      gene.sels = unique(gene.sels$gene)
      
    }else{
      gene.sels = c(Msxa.specific.genes, shared.genes.diffPattern)
      gene.sels = unique(gene.sels[!is.na(match(gene.sels, rownames(Y.mat)))])
    }
    
    cat('nb of genes to use : ', length(gene.sels), '\n')
    
    #gene.sels = unique(markers$gene[which(!is.na(match(markers$cluster, lineage)) & markers$p_val_adj<10^-5 & markers$avg_logFC > 0.7)])
    #gene.sels = gene.sels[which(!is.na(match(gene.sels, rownames(Y.mat))))]
    ##########################################
    # prepare matrix A and reponse Y and run penalized.lm
    ##########################################
    source.my.script('scMARA_utility_functions.R')
    alpha = 0;
    standardize = TRUE;
    use.lambda.min = FALSE;
    binarize.x = FALSE;
    
    for(n in 1:length(lineage.list)){
      
      # n = 1;
      lineage = lineage.list[[n]]
      
      cat(n, ' --- ')
      print(lineage)
      index.sel = match(gene.sels, rownames(Y.mat))
      Y.sel = Y.mat[index.sel, match(lineage, colnames(Y.mat))]
      y = as.matrix(Y.sel)
      
      #mm = match(rownames(Y.sel), rownames(motif.oc))
      mm = match(gene.sels, rownames(motif.oc))
      y = y[!is.na(mm), ]
      x = as.matrix(motif.oc[mm[!is.na(mm)], ])
      x[which(is.na(x) == TRUE)] = 0
      
      if(binarize.x) x = x > 0
      
      source.my.script('scMARA_utility_functions.R')
      res = run.penelized.lm(x, y, alpha = alpha, standardize = standardize, intercept = TRUE, use.lambda.min = use.lambda.min, 
                             Test = FALSE)
      
      print(res[grep('pha-4.mus|hnd-1..Tcf3|nhr-67.mus.M226|nhr-67.dm.M141|hlh-1.M175|unc-120.dm', rownames(res)),])
      
      if(n == 1) {
        keep = res[match(colnames(x), rownames(res)), ]
      }else{
        keep = cbind(keep, res[match(colnames(x), rownames(res)), ])
      }
      
    }
    
    print(keep[grep('pha-4|hnd-1..Tcf|nhr-67|hlh-1.M175|unc-120.dm', rownames(keep)),])
    
    zscore.cutoff = 2.0
    ss = apply(keep, 1, function(x) length(which(abs(x)>zscore.cutoff)))
    length(which(ss>=1))
    
    yy = keep[which(ss>=2), ] 
    
    
    yy[which(abs(yy) > 2.5)] = 2.5
    
    pdfname = paste0(resDir, "/MARA_prediction_lineages_MSxa_MSxp_try_v3.pdf")
    pdf(pdfname, width=18, height = 16)
    par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    pheatmap(yy[grep('pha-4.mus|hnd-1..Tcf3|nhr-67.mus.M226|nhr-67.dm.M141|hlh-1.M175|unc-120.dm', rownames(yy)), ], 
             cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
             scale = 'none', cluster_cols=FALSE, main = paste0("MARA prediction for controls"), 
             na_col = "white", fontsize_col = 12) 
    
    pheatmap(yy, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
             scale = 'none', cluster_cols=FALSE, main = paste0("MARA prediction"), 
             na_col = "white", fontsize_col = 12) 
    
    dev.off()
    
  }
  
  if(trajectory.resolution == 'pseudotime'){
    
    library(tradeSeq)
    library(RColorBrewer)
    library(SingleCellExperiment)
    library(slingshot)
    library(pheatmap)
    library(gridExtra)
    
    # lineage or trajectory to consider
    lineage.list = list(c('MSx', 'MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx'),
                    c('MSx', 'MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
                    )
    
    pseudotime.method = 'diffusion.map'
    cat(length(lineage.list), ' lineages or trajectories to consider \n')
    
    # prepare input matrix for tradeSeq
    ids.sels = c()
    for(m in 1:length(lineage.list)) ids.sels = unique(c(ids.sels, lineage.list[[m]]))
    cells.sels = unique(colnames(sub.obj)[!is.na(match(sub.obj$manual.annot.ids, ids.sels))])
    lineages.obj = subset(sub.obj, cells = cells.sels)
    
    counts.sel = as.matrix(lineages.obj@assays$RNA@counts)
    #crv <- newSlingshotDataSet(reducedDim = lineages.obj@reductions$umap@cell.embeddings[, c(1:2)], 
    #                           clusterLabels = lineages.obj$manual.annot.ids)
    
    pseudotime <- matrix(0, nrow = ncol(lineages.obj), ncol = length(lineage.list))
    rownames(pseudotime) = colnames(lineages.obj);
    colnames(pseudotime) = paste0('curve', c('MSxa', 'MSxp'))
    cellWeights <- pseudotime
    
    ##########################################
    # # infer pseudotime using slingshot or diffusion map
    # here the pseudotime were to be inferred and lineage-dependent genes will be identified with pseudotime
    # two options for pseudotime inferring: slingshot or diffusion map + princi_curve
    # diffusion map + princi_curve method is one trajectory one time
    ##########################################
    
    # the code were modified based on https://github.com/stevexniu/single-cell-ciona
    if(pseudotime.method == 'diffusion.map'){
      
      library(destiny)
      library(princurve)
      
      # save the pseudotime estimation 
      pdfname = paste0(resDir, "/pseudotime_estimation_v1.pdf")
      pdf(pdfname, width=12, height = 10)
      par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
      
      # loop over the trajectories
      for(n in 1:length(lineage.list)){
        
        #n = 2
        cat('lineage ', n, ': \n')
        lineage = lineage.list[[n]]
        print(lineage)
        
        ll.obj = subset(sub.obj, cells = unique(colnames(sub.obj)[!is.na(match(sub.obj$manual.annot.ids, lineage))]))
        ll.obj = FindVariableFeatures(ll.obj, selection.method = "vst", nfeatures = 2000)
        
        ll.obj = ScaleData(ll.obj, features = rownames(ll.obj))
        ll.obj <- RunPCA(object = ll.obj, features = VariableFeatures(ll.obj), verbose = FALSE, weight.by.var = FALSE)
        ElbowPlot(ll.obj, ndims = 30)
        
        nb.pcs = 10; n.neighbors = 30; min.dist = 0.3;
        ll.obj <- RunUMAP(object = ll.obj, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
        
        p1 = DimPlot(ll.obj, reduction = 'pca', group.by = 'manual.annot.ids', label = TRUE) + NoLegend()
        p2 = DimPlot(ll.obj, reduction = 'umap', group.by = 'manual.annot.ids', label = TRUE) + NoLegend()
        p1 + p2
        
        ll.pca = ll.obj@reductions$pca@cell.embeddings[, c(1:5)]
        dm <- DiffusionMap(ll.pca, sigma = 'local', n_eigs = 5, k = 100, distance = 'euclidean')
        plot(dm)
        
        #plot(dm$DC1, dm$DC2)
        dcs = as.matrix(cbind(dm$DC1, dm$DC2))
        ll.obj[["DP"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", assay = DefaultAssay(ll.obj))
        p1 = DimPlot(ll.obj, reduction = 'DP', group.by = 'manual.annot.ids')
        plot(p1)
        
        dcs = dcs[order(dcs[, 1]), ]
        princurve = principal_curve(dcs, start = dcs, smoother = 'smooth_spline', stretch = 2)
        
        plot(dcs)
        lines(princurve$s[order(princurve$lambda),], lty=1,lwd=4,col="purple",type = "l")
        whiskers(dcs, princurve$s)
        
        pseudotime.scaling = function(X) {
          return((X - min(X))/diff(range(X)))
        }
        
        pseudot = pseudotime.scaling(princurve$lambda)
        pseudot = pseudot[match(colnames(ll.obj), names(pseudot))] # match back with cell names
        
        #plot(pseudot, as.numeric(as.character(ll.obj$timingEst)),  cex = 0.5)
        Idents(ll.obj) = ll.obj$manual.annot.ids
        ll.obj$pseudotime = pseudot
        ll.obj$timingEst = as.numeric(as.character(ll.obj$timingEst))
        
        p2 = FeatureScatter(ll.obj, feature1 = "pseudotime", feature2 = "timingEst")
        plot(p2)
        # save the pseudotime in the initialized matrix
        mm = match(colnames(ll.obj), rownames(pseudotime))
        pseudotime[mm, n] = pseudot
        cellWeights[mm, n] = 1
        
      }
      
      dev.off()
      
    }
    
    # save(counts.sel, pseudotime, cellWeights, file = paste0(RdataDir, 'input_Matrix_for_tradeSeq.Rdata'))
    
    ##########################################
    # idnetify trajectory-associated genes using tradeSeq
    # the orignial code was from https://statomics.github.io/tradeSeq/articles/tradeSeq.html
    # updated version of analysis for multiple conditions were found 
    # https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html
    ##########################################
    require(tictoc)
    load(file = paste0(RdataDir, 'input_Matrix_for_tradeSeq.Rdata'))
    
    palette(brewer.pal(8, "Dark2"))
    #data(countMatrix, package = "tradeSeq")
    #counts <- as.matrix(countMatrix)
    #rm(countMatrix)
    #data(crv, package = "tradeSeq")
    #data(celltype, package = "tradeSeq")
    set.seed(5)
    tic()
    icMat <- evaluateK(counts = counts.sel, k = 3:7,
                       pseudotime = pseudotime, cellWeights = cellWeights,
                       nGenes = 200, verbose = T)
    
    toc()
    
    
    # subsetting the genes of intest rather than all genes (too slow)
    nb.knots = 4;
    
    ss = apply(counts.sel, 1, function(x) length(which(x>10)))
    genes.sel = which(ss>50)
    length(genes.sel)
    
    BPPARAM <- BiocParallel::bpparam()
    BPPARAM # lists current options
    BPPARAM$workers <- 4 # use 2 cores
    
    tic()
    set.seed(7)
    #pseudotime <- slingPseudotime(crv, na = FALSE)
    #cellWeights <- slingCurveWeights(crv)
    sce <- fitGAM(counts = counts.sel, pseudotime = pseudotime, cellWeights = cellWeights,
                  nknots = nb.knots, verbose = TRUE, parallel=TRUE, BPPARAM = BPPARAM, genes = genes.sel)
    toc()
    
    save(sce, file = paste0(RdataDir, 'fitGAM_output_tradeSeq_v3.Rdata'))
    
    load(file = paste0(RdataDir, 'fitGAM_output_tradeSeq_v3.Rdata'))
    
    mean(rowData(sce)$tradeSeq$converged)
    
    # Within-lineage comparisons : Assess if which genes are lineage- or pseudotime-dependant (basically expression is not constant) 
    # see more details in https://kstreet13.github.io/bioc2020trajectories/articles/workshopTrajectories.html
    assocRes <- associationTest(sce, lineages = TRUE, l2fc = log2(2))
    
    Msxa.genes <-  rownames(assocRes)[
      which(p.adjust(assocRes$pvalue_1, "fdr") <= 0.05)
      ]
    Msxp.genes <-  rownames(assocRes)[
      which(p.adjust(assocRes$pvalue_2, "fdr") <= 0.05)
      ]
    length(Msxa.genes)
    length(Msxp.genes)
    
    library(UpSetR)
    UpSetR::upset(fromList(list(Msxa = Msxa.genes, Msxp = Msxp.genes)))
    
    yhatSmooth <- predictSmooth(sce, gene = Msxa.genes, nPoints = 50, tidy = FALSE)
    heatSmooth <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                           cluster_cols = FALSE,
                           show_rownames = FALSE,
                           show_colnames = FALSE)
    
    shared.dynamic.genes = intersect(Msxa.genes, Msxp.genes)
    Msxa.dynamic.genes = setdiff(Msxa.genes, shared.dynamic.genes)
    Msxp.dynamic.genes = setdiff(Msxp.genes, shared.dynamic.genes)
    
    length(shared.dynamic.genes)
    length(Msxa.dynamic.genes)
    length(Msxp.dynamic.genes)
    
    gene.example = 'tbx-35'
    plotSmoothers(sce, assays(sce)$counts, gene = gene.example, alpha = 1, border = TRUE) + ggtitle(gene.example)
    
    # Between-lineage comparisons 
    # details in https://statomics.github.io/tradeSeq/articles/tradeSeq.html
    #endRes <- diffEndTest(sce) # Discovering differentiated cell type markers, used to define transiently changed genes 
    #o <- order(endRes$waldStat, decreasing = TRUE)
    #sigGene <- names(sce)[o[1]]
    #plotSmoothers(sce, assays(sce)$counts, sigGene)
    
    patternRes <- patternTest(sce, l2fc = 0)
    
    oPat <- order(patternRes$waldStat, decreasing = TRUE)
    head(rownames(patternRes)[oPat], 20)
    plotSmoothers(sce, assays(sce)$counts, gene = rownames(patternRes)[oPat][1])
    
    ##########################################
    # define lineage-shared and -specific gene groups by combining the association test and pattern test
    ##########################################
    patternRes$padj <- p.adjust(patternRes$pvalue, "fdr")
    mean(patternRes$padj <= 0.05, na.rm = TRUE)
    
    sum(patternRes$padj <= 0.01 & patternRes$fcMedian > 1.0, na.rm = TRUE)
    patternRes = patternRes[order(patternRes$waldStat, decreasing = TRUE), ]
    
    plotSmoothers(sce, assays(sce)$counts, gene = 'crt-1')
    
    genes.diffPattern = rownames(patternRes)[which(patternRes$padj <= 0.01 & patternRes$fcMedian > 1.0)]
    #genes.samePattern = setdiff(rownames(patternRes), genes.diffPattern)
    
    shared.genes.samePattern = setdiff(shared.dynamic.genes, genes.diffPattern)
    shared.genes.diffPattern = setdiff(shared.dynamic.genes, shared.genes.samePattern)
    #Msxa.specific.genes = unique(c(Msxa.dynamic.genes, shared.genes.diffPattern))
    #Msxp.specific.genes = unique(c(Msxp.dynamic.genes, shared.genes.diffPattern))
    
    Msxa.specific.genes = Msxa.dynamic.genes
    Msxp.specific.genes = Msxp.dynamic.genes
    
    for(gene.module in c('shared.genes.samePattern', 'shared.genes.diffPattern', 'Msxa.specific.genes', 'Msxp.specific.genes'))
    {
      if(length(dev.list()!=0)) { dev.off()}
      
      pdfname = paste0(resDir, "/lineage_dependent_genes_MSxa_MSxp_", gene.module, ".pdf")
      
      pdf(pdfname, width=16, height = 10)
      par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
      
      ggs = eval(parse(text= paste0(gene.module)))
      for(g in ggs){
        cat(g, '\n')
        kk = which(tfs$`Public name` == g)
        gtitle = g
        if(length(kk) >0 ) gtitle = paste0(gtitle, '-- TF') 
        p1 = plotSmoothers(sce, assays(sce)$counts, gene = g, alpha = 1, border = TRUE,
                           nPoints = 100, size = 1) + ggtitle(gtitle)
        plot(p1)
      }
      dev.off()
      
    }
    
    
    # heatmap showing different groups
    yhatSmooth <- predictSmooth(sce, gene = Msxa.specific.genes, nPoints = 50, tidy = FALSE)
    heatSmooth_MSxa <- pheatmap(t(scale(t(yhatSmooth[, 1:50]))),
                               cluster_cols = FALSE,
                               show_rownames = FALSE, show_colnames = FALSE, main = "MSxa", legend = FALSE,
                               silent = TRUE
    )
    
    matchingHeatmap_MSxp <- pheatmap(t(scale(t(yhatSmooth[heatSmooth_MSxa$tree_row$order, 51:100]))),
                                     cluster_cols = FALSE, cluster_rows = FALSE,
                                     show_rownames = FALSE, show_colnames = FALSE, main = "MSxp",
                                     legend = FALSE, silent = TRUE
    )
    
    grid.arrange(heatSmooth_MSxa[[4]], matchingHeatmap_MSxp[[4]], ncol = 2)
    
    
    # Example on combining patternTest with diffEndTest results
    Define.transient.genes = FALSE
    if(Define.transient.genes){
      library(ggplot2)
      patternRes$Gene <- rownames(patternRes)
      patternRes$pattern <- patternRes$waldStat
      patternRes <- patternRes[, c("Gene", "pattern")]
      
      endRes$Gene <- rownames(endRes)
      endRes$end <- endRes$waldStat
      endRes <- endRes[, c("Gene", "end")]
      
      compare <- merge(patternRes, endRes, by = "Gene", all = FALSE)
      compare$transientScore <- 
        rank(-compare$end, ties.method = "min")^2 + rank(compare$pattern, ties.method = "random")^2
      
      ggplot(compare, aes(x = log(pattern), y = log(end))) +
        geom_point(aes(col = transientScore)) +
        labs(x = "patternTest Wald Statistic (log scale)",
             y = "diffEndTest Wald Statistic (log scale)") +
        scale_color_continuous(low = "yellow", high = "red") +
        theme_classic()
      
      topTransient <- compare[which.max(compare$transientScore), "Gene"]
      plotSmoothers(sce, assays(sce)$counts, gene = topTransient)
      
    }
    
    
    ## plot the gene profiles in pseudotime
    plot.TFs.in.pseudotime = FALSE
    if(plot.TFs.in.pseudotime){
      gene.tfs = intersect(rownames(assocRes), tfs$`Public name`)
      
      pdfname = paste0(resDir, "/lineage_dependant_genes_MSxa_MSxp_only_TFs.pdf")
      pdf(pdfname, width=16, height = 10)
      par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
      
      for(g in gene.tfs){
        cat(g, '\n')
        kk = which(tfs$`Public name` == g)
        gtitle = g
        if(length(kk) >0 ) gtitle = paste0(gtitle, '-- TF') 
        p1 = plotSmoothers(sce, assays(sce)$counts, gene = g, alpha = 1, border = TRUE,
                           nPoints = 100, size = 1) + ggtitle(gtitle)
        plot(p1)
      }
      
      dev.off()
      
    }
    
  }
  
  
}
