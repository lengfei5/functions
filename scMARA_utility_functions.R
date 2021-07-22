##########################################################################
##########################################################################
# Project: C elegans embryogenesis for MS lineage
# Script purpose: utility functions for scMARA
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Fri Oct 23 18:27:14 2020
##########################################################################
##########################################################################
# main function 
run.penelized.lm = function(x, y, alpha = 0, intercept=TRUE, standardize=FALSE,  standardize.response=FALSE, use.lambda.min = TRUE,
                            zscore.output = TRUE, Test = FALSE, Test.zscore.cutoff = 2.0)
{
  # intercept=TRUE; standardize=TRUE;  standardize.response=FALSE; alpha = 0.5; use.lambda.min = TRUE; zscore.output = TRUE;
  
  # check the x and y
  cat(ncol(x), ' motifs \n')
  cat(nrow(x), ' gene selected \n')
  
  cat(length(which(apply(x, 2, sum) == 0)), '  motifs without target genes  \n')
  cat(length(which(apply(x, 1, sum) == 0)), ' genes without motifs \n')
  
  ### standardize matrix of motif occurrence makes more sense because the absolute number of motif occurrence is not precise.
  if(is.null(ncol(y))){
    family = 'gaussian'
    cat('responses nb : 1  \n')
  }else{
    family = 'mgaussian'
    cat('responses nb :',  ncol(y), '\n')
  }
  
  cat('-- start the penalized linear regression -- \n')
  ### use Cross-validation to select tuning paprameter
  cv.fit=cv.glmnet(x, y, family=family, grouped=FALSE, 
                   alpha=alpha, nlambda=200, standardize=standardize, 
                   standardize.response=standardize.response, intercept=intercept, relax = FALSE)
  plot(cv.fit)
  
  if(use.lambda.min){
    s.optimal = cv.fit$lambda.min
  }else{
    s.optimal = cv.fit$lambda.1se
  }
  fit=glmnet(x,y,alpha=alpha, lambda=s.optimal, family=family, 
             standardize=standardize, standardize.response=standardize.response, intercept=intercept, 
             relax = FALSE)
  #plot(fit, xvar = "lambda", label = TRUE, type.coef = "coef")
  
  # extract fitting results for either multiple response or single response; 
  # in particular, we decided to use only ridge penalized linear regression here
  if(family == 'mgaussian'){
    keep = as.data.frame(coef.glmnet(fit, s = s.optimal))
    keep = keep[-1, ] # remove intecept
    colnames(keep) = names(fit$beta)
    if(alpha > 0.0){
      rownames(keep) = rownames(fit$beta[[2]])
      res = keep;
      ## collect result from the elastic-net
      #kk = apply(keep, 1, function(x) !all(x==0))
      #keep = keep[kk, ]
      #colnames(x)[which(fit$beta[[1]][,optimal]!=0)]
      #colnames(x)[which(fit$beta[[2]][,optimal]!=0)]
      
      # rerun lm with selected features
      relax.fitting.lm = FALSE
      if(relax.fitting.lm){
        fit.lm = lm(y ~ x[, match(rownames(keep), colnames(x))])
        res = data.frame(fit.lm$coefficients)
        res = res[-1, ] # remove intercept
        rownames(res) = rownames(keep)
        
        pheatmap(res, cluster_rows=TRUE, show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
                 scale = 'column', cluster_cols=FALSE, main = paste0("motif activity"), 
                 na_col = "white", fontsize_col = 10)
      }
      
    }else{
      if(zscore.output) keep = apply(keep, 2, scale)
      rownames(keep) = rownames(fit$beta[[2]])
      #ss = apply(keep, 1, function(x) !all(x==0))
      #keep = keep[ss, ]
      #head(rownames(keep)[order(-abs(keep$MSxp))], 10)
      #head(rownames(keep)[order(-abs(keep$MSxa))], 10)
      res = keep
      
      if(Test){
        ss = apply(keep, 1, function(x) length(which(abs(x) > Test.zscore.cutoff)))
        keep = keep[which(ss>0), ]
        print(keep)
      }
    }
    
  }else{
    keep = coef.glmnet(fit, s = s.optimal)
    motif.names = keep@Dimnames[[1]]
    motif.names = motif.names[-1]
    keep = keep[-1]
    names(keep) = motif.names
    
    if(zscore.output) keep = scale(keep)
    #o1 = order(-abs(keep))
    #keep = keep[o1]
    #motif.names = motif.names[o1]
    #res = data.frame(motif = motif.names, scores = keep, stringsAsFactors = FALSE)
    res = keep
    
    if(Test) print(res[which(abs(res$scores) > Test.zscore.cutoff), ])
    
    # if(alpha > 0.0) {
    #   keep = keep[which(keep != 0)]
    #   print(names(keep))
    # }else{
     
  }
  
  #for(n in 1:ncol(keep)) keep[,n] = scale(keep[,n], center = TRUE, scale = TRUE)
  #keep[which(abs(keep)<0.05)] = 0
  #keep = data.frame(keep, stringsAsFactors = FALSE)
  return(res)
  
}

################################################################################################################
# utility functions for sctf_MARA
# part-1 : worm gene annotation convertion
# part-2 : motif occurrency matrix preparation from fimo output
# part-3 : motif-tf mapping table
# part-4 : pwm clustering based on sequence similarity (here the binding site overlapping is not considered as before) 
# and modify motif-tf mapping 
# part-5: prepare the detecte tf expression profiles
################################################################################################################

########################################################
########################################################
# Section II : process target genes (Respone Y), gene expression matrix
# 
# process worm gene tss to have unique tss for each gene. 
# to do this, we just pick the tss furthest from the gene start so that the promoter can cover as much as regulatory elements
# in addition, save the gene length, transcript length for the scRNA-seq length normalization
########################################################
########################################################
process.worm.gene.tss = function()
{
  rm(list=ls())
  setwd('/Volumes/groups/cochella/jiwang/annotations')
  
  tss = read.table('ce11_tss.bed', header = FALSE, sep = '\t')
  load('BioMart_WBcel235.Rdata')
  annot = annot[which(annot$Gene.type == 'protein_coding'), ]
  
  # filter non-protein-coding genes
  mm = match(tss$V4, annot$Gene.stable.ID)
  tss = tss[which(!is.na(mm)), ]
  
  gg.uniq = unique(tss$V4)
  keep = rep(NA, length(gg.uniq))
  names(keep) = gg.uniq
  gg.counts = table(tss$V4)
  
  gg.with.singleTss = names(gg.counts)[which(gg.counts == 1)]
  #gg.with.multiTss = names(gg.counts)[which(gg.counts > 1)]
  keep[match(gg.with.singleTss, names(keep))] = match(gg.with.singleTss, tss$V4) 
  
  nn = which(is.na(keep))
  for(n in nn)
  {
    kk = which(tss$V4 == names(keep)[n])
    if(length(kk) == 1){
      cat('Error \n')
    }else{
      cat(n, '--', as.character(names(keep)[n]), '\n')
      #jj = which(annot$)
      if(unique(tss$V6[kk]) == '+'){
        keep[n] = kk[which(tss$V2[kk] == max(tss$V2[kk]))][1]
      }else{
        keep[n] = kk[which(tss$V2[kk] == min(tss$V2[kk]))][1]
      }
    }
  }
  
  tss = tss[keep, ]
  #write.table(tss, file = 'ce11_tss_curated_singleTss_perGene_proteinCoding.bed', sep = '\t', col.names = FALSE, row.names = FALSE,
  #            quote = FALSE)
  
  
  ## save the gene length and transcript length (averge of isoform lengths)
  aa = data.frame(annot$Gene.stable.ID, annot$Gene.Start..bp., annot$Gene.End..bp., annot$Transcript.length..including.UTRs.and.CDS., 
                  annot$Gene.name, annot$Gene.type, stringsAsFactors = FALSE)
  
  aa$gene.length = abs(aa$annot.Gene.End..bp. - aa$annot.Gene.Start..bp.)
  
  gg.uniq = unique(aa$annot.Gene.stable.ID)
  keep = aa[match(gg.uniq, aa$annot.Gene.stable.ID), ]
  colnames(keep) = c('wormbase.id', 'gene.start', 'gene.end', 'transcript.length', 'gene.name', 'gene.type', 'gene.length')
  
  for(n in 1:nrow(keep))
  {
    # n = 1
    jj = which(aa$annot.Gene.stable.ID == keep$wormbase.id[n])
    if(length(jj) > 1){
      cat(n, '--', as.character(keep$wormbase.id[n]), '\n')
      keep$transcript.length[n] = as.integer(median(aa$annot.Transcript.length..including.UTRs.and.CDS.[jj]))
    }
  }
  #saveRDS(keep, file = 'ce11_proteinCoding_genes_geneLength_transcriptLength.rds')
}

display.gene.expression.MS.lineage = function(tf.mat)
{
  library("treeio")
  library("ggtree")
  
  bwm_tree = readRDS(file = paste0(dataDir, 'BWM_tree_for_visualization.rds'))
  
  source.my.script('make_lineage_ggtree_Viscello.R')
  ids.names = colnames(tf.mat)
  ids.names[which(ids.names == "MSxppapp/MSxpappp")] = 'MSxpappp'
  
  pdfname = paste0(resDir, "/BWM_lineage_expressed_TFs_profiles_v2.pdf")
  pdf(pdfname, width=10, height = 8)
  par(cex =0.7, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  options(warn=-1)
  for(n in 1:nrow(tf.mat))
  #for(n in 1:10)
  {
    cat(n, ' -- ', rownames(tf.mat)[n], '\n')
    
    tf.expr = tf.mat[n, ]
    bwm_tree$value = tf.expr[match(bwm_tree$lineage, ids.names)]
    out.tree = make_lineage_ggtree(bwm_tree, root = 'MS', color.annot = "value") + 
      ggtitle(rownames(tf.mat)[n])
    plot(out.tree)
    
  }
  
  options(warn=0)
  dev.off()
  
  #nwk <- system.file("extdata", "sample.nwk", package="treeio")
  #tree <- read.tree(nwk)
  # ggplot(tree, aes(x, y)) + geom_tree() + theme_tree()
  # ggtree(tree, color="firebrick", size=2, linetype="dotted")
  # ggtree(tree, ladderize=TRUE)
  # ggtree(tree, branch.length="none")
  # 
  # 
  # beast_file <- system.file("examples/MCC_FluA_H3.tree", 
  #                           package="ggtree")
  # beast_tree <- read.beast(beast_file)
  # ggtree(beast_tree, mrsd="2013-01-01") + theme_tree2()
  
  #ggtree(tree_tbl, )
  
}

process.detected.tf.expression.profiles = function(Y.mat)
{
  # subset Y.mat for TFs
  jj = match(tfs$`Public name`, rownames(Y.mat))
  jj = jj[!is.na(jj)]
  tf.mat = Y.mat[jj, ]
  cutoff.tf = 1;
  ss = apply(tf.mat, 1, function(x) !all(x<cutoff.tf))
  tf.mat = tf.mat[ss, ]
  
  save.tf.profiles.across.lineage = FALSE
  if(save.tf.profiles.across.lineage){
    pdfname = paste0(resDir, "/TFs_standardized_fpkm_in_BWM.pdf")
    pdf(pdfname, width=12, height = 50)
    par(cex =0.3, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    pheatmap(tf.mat, cluster_rows=TRUE, 
             show_rownames=TRUE, show_colnames = TRUE, breaks = NA,
             scale = 'row',
             cluster_cols=FALSE, 
             main = paste0("standardized fpkm of TFs "), 
             na_col = "white",
             #color = cols, 
             #annotation_col = my_sample_col,
             #gaps_row = c(1:nrow(map)-1),
             fontsize_col = 10
    )
    dev.off()
    
    
    write.csv(tf.mat, file = paste0(tabDir, 'detected_TFs_in_BWM.csv'), row.names = TRUE)
  }
  
  saveRDS(tf.mat, file = paste0(RdataDir, 'TFs_expression_profiles_BWM.rds')) 
  
}

##########################################
# important function to define genes relevant to lineage
##########################################
define.modules.for.lineags = function(sub.obj, Y.fpkm, lineage = c('MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp'))
{
  ids = sub.obj$manual.annot.ids
  ids.uniq = unique(ids)
  ids.uniq = ids.uniq[order(nchar(ids.uniq))]
  
  ids.uniq = ids.uniq[grep('mixture_terminal_', ids.uniq, invert = TRUE)]
  
  if(mode == 'cluster.mode'){
    cat('-- cluster mode : averging the gene expression in clusters -- \n')
    
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
    
    select.dyn.genes.with.FindAllMarker.MST = FALSE
    if(select.dyn.genes.with.FindAllMarker.MST){
      run.FindAllMarkers = FALSE # it takes ~ 30 minutes
      if(run.FindAllMarkers){
        Idents(sub.obj) = sub.obj$manual.annot.ids
        markers.new <- FindAllMarkers(sub.obj, only.pos = FALSE, min.pct = 0.1, logfc.threshold = 0.25)
        saveRDS(markers.new, file = paste0(RdataDir,  'AllMarkers_MST_manual.annotation.rds'))
      }else{
        markers = readRDS(file = paste0(RdataDir,  'AllMarkers_MST_manual.annotation.rds'))
      }
      
      gene.sels = markers[which(markers$p_val<10^-3 & markers$avg_logFC > 0.5), ]
      print(table(gene.sels$cluster))
      
      gene.sels = unique(gene.sels$gene)
      gene.sels = gene.sels[which(!is.na(match(gene.sels, rownames(Y.mat))))]
      cat(length(gene.sels), ' marker genes were identified \n')
      
      means = apply(Y.mat, 1, mean)
      vars = apply(Y.mat, 1, var)
      
      #top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
      #DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
      mm = match(gene.sels, rownames(Y.mat))
      
      plot(means[mm], vars[mm], cex = 0.3)
      
      pheatmap(Y.mat[mm, ], cluster_rows=TRUE, show_rownames=FALSE, show_colnames = TRUE, breaks = NA,
               gaps_col = c(6, 20),
               scale = 'row', cluster_cols=FALSE, main = paste0("dynamic genes by marker finding"), 
               na_col = "white", fontsize_col = 10
      )
      
      pheatmap(Y.mat[mm, ], cluster_rows=TRUE, show_rownames=FALSE, show_colnames = TRUE, breaks = NA,
               gaps_col = c(6, 20), cutree_rows = 5, 
               scale = 'row', cluster_cols=FALSE, main = paste0("dynamic genes by marker finding"), 
               na_col = "white", fontsize_col = 10
      )
      
      
      cluster.coexprssed.genes = FALSE
      if(cluster.coexprssed.genes){
        # cluster co-expressing modules 
        library(tidyverse)  # data manipulation
        library(cluster)    # clustering algorithms
        library(factoextra) # clustering algorithms & visualization
        
        df <- t(scale(t(Y.mat[mm, ])))
        head(df)
        #distance <- get_dist(df, method = 'euclidean')
        #fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
        
        # Elbow method to determine nb of clusters
        set.seed(123)
        fviz_nbclust(df, kmeans, method = "wss", k.max = 30)
        
        # set.seed(123)
        # gap_stat <- clusGap(df, FUN = kmeans, nstart = 25,
        #                     K.max = 20, B = 50)
        # 
        # print(gap_stat, method = "firstmax")
        # fviz_gap_stat(gap_stat)
        
        pheatmap(df, cluster_rows=TRUE, show_rownames=FALSE, show_colnames = TRUE, breaks = NA,
                 #clustering_distance_rows = "euclidean", 
                 cutree_rows = 20,
                 scale = 'none', cluster_cols=FALSE, main = paste0("dynamic genes by marker finding"), 
                 na_col = "white", fontsize_col = 10
        )
        #set.seed(123)
        #final <- kmeans(df, 4, nstart = 25)
        #print(final)
      }
      
    }
    
    ##########################################
    # select dynamic genes for reponse Y
    # 1) with fano 
    # 2) ratio betwen daughter and mother, a lot of pairs to take care
    # 3) by lineage e.g. MSx, MSxa, MSxap, MSxapp, MSxappp, MSxapppp, MSxappppx (e.g. with gam)
    ##########################################
    select.dyn.genes.with.fano = FALSE
    if(select.dyn.genes.with.fano){
      ss = apply(Y.mat, 1, mean)
      fano = apply(Y.mat, 1, var)/ss
      plot(ss, fano, cex = 0.6);
      abline(h = c(0.5,  0.7, 1.0), col = 'blue', lwd=1.2)
      length(which(fano > 1.5))
      length(which(fano > 1.0))
      length(which(fano > 0.7))
      length(which(fano > 0.5))
      #length(which(fano > 0.3))
      
      Y.sel = Y.mat[which(fano > 1.5), ]
    }
    
    select.dyn.genes.with.pair.ratios = FALSE
    if(select.dyn.genes.with.pair.ratios){
      
      Y.mat = as.data.frame(Y.mat)
      rownames(Y.sel) = rownames(Y.mat)
      colnames(Y.sel) = c('MSxa', 'MSxp')
      
      hist(Y.sel, breaks = 100);abline(v = c(-1, 1))
      cutoff = 1;
      sels = apply(Y.sel, 1, function(x) sum(abs(x)> cutoff)>1)
      cat(sum(sels), ' gene were selected \n')
      Y.sel = Y.sel[sels, ]
      
    }
  
  }
  
  ##########################################
  # the following code is modified based on original code from https://github.com/stevexniu/single-cell-ciona
  # in addition, the following code is redundant with part of code in scMARA.R (to do)
  ##########################################
  if(mode == 'time.bin'){
    library(destiny)
    library(princurve)
    
    lineage = c('MSx', 'MSxa', 'MSxap', 'MSxapp', 'MSxappp', 'MSxapppp', 'MSxappppx')
    lineage = c('MSx', 'MSxp', 'MSxpp', 'MSxppp', 'MSxpppp', 'MSxppppp')
    
    cells.sels = unique(colnames(sub.obj)[!is.na(match(sub.obj$manual.annot.ids, lineage))])
    ll.obj = subset(sub.obj, cells = cells.sels)
    ll.obj = FindVariableFeatures(ll.obj, selection.method = "vst", nfeatures = 1000)
    ll.obj = RunPCA(ll.obj, features = VariableFeatures(ll.obj), npcs = 50, verbose = FALSE, weight.by.var = TRUE)
    
    ll.pca = ll.obj@reductions$pca@cell.embeddings[, c(1:50)]
    dm <- DiffusionMap(ll.pca, sigma = 'local', n_eigs = 5)
    #plot(dm)
    #plot(dm$DC1, dm$DC2)
    dcs = as.matrix(cbind(dm$DC1, dm$DC2))
    ll.obj[["DP"]] <- CreateDimReducObject(embeddings = as.matrix(dcs), key = "DC_", assay = DefaultAssay(ll.obj))
    DimPlot(ll.obj, reduction = 'DP', group.by = 'manual.annot.ids')
    
    princurve = principal_curve(dcs, start = dcs, smoother = 'lowess', stretch = 2)
    
    plot(dcs)
    lines(princurve$s[order(princurve$lambda),], lty=1,lwd=4,col="purple",type = "l")
    
    whiskers(dcs, princurve$s)
    
    pseudot = princurve$lambda
    #plot(order(-diffmap$X[,3]), ll.obj$timingEst,  cex = 0.5)  
    plot(pseudot, ll.obj$timingEst,  cex = 0.5)  
    #mat.dist = as.matrix(dist(ll.pca, method = 'euclidean'))
    # Run diffusion map and return top 50 dimensions
    #set.seed(1)
    #diffmap = diffuse(mat.dist, maxdim=50)
    
    #diffmap.embedding = as.matrix(diffmap$X)
    #rownames(diffmap.embedding) = colnames(ll.obj)
    #colnames(diffmap.embedding) = paste0('diffumap_', c(1:ncol(diffmap.embedding)))
    
    #ll.obj[['diffmap']] = Seurat::CreateDimReducObject(embeddings=diffmap.embedding , key='diffmap_', assay='RNA')
    
    # Save first two diffusion map coordinators
    #shp@tsne.rot[1:2]=data.frame(shp.diff$X[,1:2],row.names = ll.obj@cell.names)
    #colnames(shp@tsne.rot)=c("tSNE_1","tSNE_2")
    # Visualize top two diffusion map components
    #tsne.pseudo(shp, do.label = F,label.cex.text = 1,name.y = "Diffusion Map Coordinator 2",name.x = "Diffusion Map Coordinator1",label.cols.use = c("green","yellow","orange1","orange4","orange4"),label.pt.size = 1.5,xlim=c(-0.1,0.05),ylim=c(-0.05,0.05))
    #legend("topleft",legend=c("12TVC","14STVC","16SHP","18SHP","20SHP"),col= c("green","yellow","orange1","orange4","orange4"),pch = 16,cex=0.5,pt.cex = 1)
    #plot(diffmap$X[, c(1:2)])
    # Fit the first two diffusion map components with principal curve
    
    
    df=data.frame(princurve$s[order(princurve$lambda),]);colnames(df) = c("x","y")
    
    ggplot(data=df,aes(x,y))+
      geom_line(size=1.5,colour="black")+
      geom_density2d(aes(colour=..level..),bins=6) + 
      scale_colour_gradient(low="darkgray",high="white",3) +
      xlim(-0.084,0.08) + ylim(-0.05,0.05) +
      geom_point(data=data.frame(shp@tsne.rot,color=shp@ident),aes(tSNE_1,tSNE_2),size=2,color=c(rep("green",table(shp@ident)[1]),rep("yellow",table(shp@ident)[2]),rep("orange1",table(shp@ident)[3]),rep("orange4",table(shp@ident)[4]),rep("orange4",table(shp@ident)[5]))) +
      theme_classic() +  
      theme(legend.position="none",axis.title=element_text(size = rel(1)),axis.text=element_blank(), axis.ticks = element_blank(),axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'))+ xlab(label = "Diffusion Component 1") + ylab(label = "Diffusion Component 2") +     geom_line(size=1.5,colour="black")
    
    
  }
  
  
}

##########################################
# function to infer pseudotime 
##########################################
pseudotime.inferrence = function()
{
  if(pseudotime.method == 'slingshot'){
    pca <- prcomp(t(log1p(ll.obj@assays$RNA@data)[match(VariableFeatures(ll.obj), rownames(ll.obj)), ]), scale. = FALSE)
    rd1 <- pca$x[,1:2]
    plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
    #plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)
    
    library(uwot)
    
    rd2 <- umap(t(log1p(ll.obj@assays$RNA@data)[match(VariableFeatures(ll.obj), rownames(ll.obj)), ]))
    colnames(rd2) <- c('UMAP1', 'UMAP2')
    
    plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)
    DimPlot(ll.obj, reduction = 'umap', group.by = 'manual.annot.ids')
    
  }
}


##########################################
# compare mRNA and unspliced matrix 
##########################################
compare.mRNA.and.unspliced.matrix = function(sub.obj)
{
  ll = readRDS(file = '../data/motifs_tfs/ce11_proteinCoding_genes_geneLength_transcriptLength.rds') # gene length and transript length
  
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
  
  # convert to SingleCellExperiment, recalculate scaling factor and normalized to fpkm
  pm = readRDS(file = paste0(RdataDir, 'unsplicedData_processed_by_velocyto_scranNormalized_BWM.rds'))
  
  sce = as.SingleCellExperiment(pm)
  qclust <- quickCluster(sce)
  sce <- computeSumFactors(sce, clusters = qclust)
  sce <- logNormCounts(sce, log = TRUE, pseudo_count = 1)
  
  mm = match(rownames(sce), ll$gene.name)
  sce = sce[which(!is.na(mm)), ] # keep genes with correponding lengths
  intron.length  = (ll$gene.length - ll$transcript.length)[mm[which(!is.na(mm) == TRUE)]]
  
  sce <- logNormCounts(sce, log = FALSE, size_factors = NULL)
  unspliced <- log2(calculateFPKM(sce, lengths = intron.length) + 1)
  
  remove(sce)
  
  ##########################################
  # average the gene expressoin for ids
  ##########################################
  ids = sub.obj$manual.annot.ids
  ids.uniq = unique(ids)
  ids.uniq = ids.uniq[order(ids.uniq)]
  ids.uniq = ids.uniq[order(nchar(ids.uniq))]
  
  ids.uniq = ids.uniq[grep('mixture_terminal_', ids.uniq, invert = TRUE)]
  
  Y.mat = matrix(NA, nrow = nrow(Y.fpkm), ncol = length(ids.uniq))
  colnames(Y.mat) = ids.uniq
  rownames(Y.mat) = rownames(Y.fpkm)
  
  for(n in 1:length(ids.uniq))
  {
    cat(ids.uniq[n], '\n')
    jj = which(sub.obj$manual.annot.ids == ids.uniq[n])
    if(length(jj) == 1) Y.mat[, n] = Y.fpkm[,jj]
    if(length(jj) > 1) Y.mat[, n] = apply(Y.fpkm[,jj], 1, mean)
  }
  
  P.mat = matrix(NA, nrow = nrow(unspliced), ncol = length(ids.uniq))
  colnames(P.mat) = ids.uniq
  rownames(P.mat) = rownames(unspliced)
  
  for(n in 1:length(ids.uniq))
  {
    cat(ids.uniq[n], '\n')
    jj = which(pm$annoted.ids == ids.uniq[n])
    if(length(jj) == 1) P.mat[, n] = unspliced[,jj]
    if(length(jj) > 1) P.mat[, n] = apply(unspliced[,jj], 1, mean)
  }
  
  remove(unspliced)
  remove(Y.fpkm)
  
  # select dynamic genes
  markers = readRDS(file = paste0(RdataDir,  'AllMarkers_MST_manual.annotation.rds'))
  gene.sels = markers[which(markers$p_val<10^-3 & markers$avg_logFC > 0.5), ]
  print(table(gene.sels$cluster))
  
  gene.sels = unique(gene.sels$gene)
  gene.sels = gene.sels[which(!is.na(match(gene.sels, rownames(Y.mat))))]
  cat(length(gene.sels), ' marker genes were identified \n')
  
  means = apply(Y.mat, 1, mean)
  vars = apply(Y.mat, 1, var)
  
  #top.markers <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
  #DoHeatmap(sub.obj, features = top.markers$gene, size = 5, hjust = 0, label = TRUE) + NoLegend()
  mm = match(gene.sels, rownames(Y.mat))
  
  plot(means[mm], vars[mm], cex = 0.3)
  
  yy1 = Y.mat[mm, ]
  yy1 = t(scale(t(yy1)))
  
  d <- dist(yy1, method = "euclidean") # distance matrix
  fit <- hclust(d, method="complete")
  #plot(fit) # display dendogram
  groups <- cutree(fit, k=10) 
  yy1 = yy1[order(groups), ]
  my_gene_clusters = as.data.frame(groups)
  #my_gene_col <- cutree(tree = as.dendrogram(my_hclust_gene), k = 2)
  yy1[which(abs(yy1)>2.)] = 2.
  
  pdfname = paste0(resDir, "/mRNA_preRNA_comparison.pdf")
  pdf(pdfname, width=12, height = 8)
  par(cex =0.5, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  pheatmap(yy1, cluster_rows=FALSE, show_rownames=FALSE, show_colnames = TRUE, 
           gaps_col = c(6, 20), breaks = NA,
           #annotation_row = my_gene_clusters,
           scale = 'none', cluster_cols=FALSE, main = paste0("mRNA - dynamic genes by marker finding"), 
           na_col = "white", 
           fontsize_col = 10
  )
  
  jj = match(rownames(yy1), rownames(P.mat))
  yy2 = P.mat[jj, ]
  yy2 = t(apply(yy2, 1, scale))
  yy2[which(abs(yy2)>2.)] = 2.
  colnames(yy2) = colnames(P.mat)
  
  pheatmap(yy2, cluster_rows=FALSE, show_rownames=FALSE, show_colnames = TRUE, breaks = NA,
           gaps_col = c(6, 20), 
           scale = 'none', cluster_cols=FALSE, main = paste0("pre-RNA - dynamic genes by marker finding"), 
           na_col = "white", fontsize_col = 10
  )
  
  nb.genes = apply(P.mat[jj, ], 2, function(x) return(length(which(x>10^-6))))
  
  dev.off()
  
  
}

########################################################
########################################################
# Section III : process motifs and scanning regions 
# designa matrix A in Ax = Y
# 
########################################################
########################################################

# manually add extra motifs for hnd-1, pha-4, unc-120 and nhr-67 from dm, mus and homo
convert.cisbp.format.to.meme = function()
{
  library(universalmotif)
  #pwmDIr = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/extra_pwm'
  pwm.cisbp = '../data/test/PWM.txt'
  xx = read_cisbp(file = pwm.cisbp, skip = 0)
  yy = convert_type(xx, "PWM")
  
  write_meme(yy, file = '../data/test/PWM_converted.meme', overwrite = TRUE)
  
}

generate.logos.for.motifs.pwm = function()
{
  library('universalmotif')
  pwm = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/All_PWMs_JASPAR_CORE_2016_TRANSFAC_2015_CIS_BP_2015_curated_extra.meme'
  #pwm = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/All_PWMs_JASPAR_CORE_2016_TRANSFAC_2015_CIS_BP_2015.meme'
  
  meme= read_meme(file = pwm, skip = 0, readsites = FALSE, readsites.meta = FALSE)
  motifs = convert_motifs(meme, class = "universalmotif-universalmotif")
  # comparisons <- compare_motifs(motifs, method = "PCC", min.mean.ic = 0,
  #                               score.strat = "a.mean")
  # write.table(comparisons, file = '../data/motifs_tfs/pwm_similarity_correction_PCC.txt', sep = '\t', col.names = TRUE, 
  #             row.names = TRUE, quote = FALSE)
  
  p1 = view_motifs(motifs[[1]], use.type = 'ICM')
  plot(p1)
  
  for(n in 1:length(motifs))
  {
    cat(n, '\n')
    pdfname = paste0('../data/motifs_tfs/pwm_logos/', motifs[[n]]@name, '.pdf')
    pdf(pdfname, width=8, height = 6)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    p1 = view_motifs(motifs[[n]], use.type = 'ICM')
    plot(p1)
    
    dev.off()
  }
  
}

## modify the motif names with associated TFs
process.motif.tf.mapping = function()
{
  motif.tf = read.table('/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/motifs_tfs_mapping_curated_extra.txt', 
                        header = FALSE, sep = ' ')
  motif.tf = motif.tf[, c(2:3)]
  colnames(motif.tf) = c('motifs', 'tfs')
  
  # manually modify the motif names
  motif.tf = data.frame(motif.tf, stringsAsFactors = FALSE)
  motif.tf$motifs.new = motif.tf$motifs
  motif.tf$tfs.new = motif.tf$tfs
  
  xx = motif.tf
  #xx$motifs.new = gsub('1.02', '', xx$motifs.new)
  #xx$motifs.new = gsub('1.02', '', xx$motifs.new)
  #xx$tfs.new = paste0(xx$tfs.new, '_', xx$tfs)
  #xx$tfs.new = gsub('-', '', xx$tfs.new)
  xx$tfs.new = gsub(':', '.', xx$tfs.new)
  xx$tfs.new = gsub('/', '.', xx$tfs.new)
  xx$tfs.new = gsub("\\(","", xx$tfs.new)
  xx$tfs.new = gsub("\\)","", xx$tfs.new)
  xx$tfs.new = gsub("_Homo_sapiens_DBD*.*",".homo", xx$tfs.new)
  xx$tfs.new = gsub("_Caenorhabditis_briggsae_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_Drosophila_melanogaster_DBD*.*",".dm", xx$tfs.new)
  xx$tfs.new = gsub("_Mus_musculus_DBD*.*",".mus", xx$tfs.new)
  xx$tfs.new = gsub("_Brugia_pahangi_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_Wuchereria_bancrofti_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_PBM_CONSTRUCTS_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub("_Tetraodon_nigroviridis_DBD*.*","", xx$tfs.new)
  xx$tfs.new = gsub('Mus', '.mus', xx$tfs.new)
  xx$tfs.new = gsub('Dm', '.dm', xx$tfs.new)
  xx$tfs.new = gsub('Homo', '.homo', xx$tfs.new)
  xx$tfs.new = gsub('Hand1', 'hnd-1', xx$tfs.new)
  xx$tfs.new = gsub('Foxa1', 'pha-4', xx$tfs.new)
  xx$tfs.new = gsub('SRF', 'unc-120', xx$tfs.new)
  xx$tfs.new = gsub('Srf', 'unc-120', xx$tfs.new)
  xx$tfs.new = gsub('Nr2e1', 'nhr-67', xx$tfs.new)
  xx$tfs.new = gsub('NR2E1', 'nhr-67', xx$tfs.new)
  
  xx$motifs.new = paste0(xx$motifs.new, '_', xx$tfs.new)
  
  motif.tf = xx
  
  saveRDS(motif.tf, file = '../data/motifs_tfs/motif_tf_mapping.rds')
  
}

# remove motif redundancies
remove.motifs.redundancy.by.similarity.clustering = function()
{
  library(tidyverse)  # data manipulation
  library(cluster)    # clustering algorithms
  library(factoextra) # clustering visualization
  require(graphics)
  library(universalmotif)
  
  motif.tf = readRDS(file = '../data/motifs_tfs/motif_tf_mapping.rds')
  motif.tf$motifs.new = paste0(motif.tf$tfs.new, '_', motif.tf$motifs)
  
  pwm = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/All_PWMs_JASPAR_CORE_2016_TRANSFAC_2015_CIS_BP_2015_curated_extra.meme'
  #pwm = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/PWMs_C_elegans/All_PWMs_JASPAR_CORE_2016_TRANSFAC_2015_CIS_BP_2015.meme'
  
  meme= read_meme(file = pwm, skip = 0, readsites = FALSE, readsites.meta = FALSE)
  motifs = convert_motifs(meme, class = "universalmotif-universalmotif")
  cat(length(motifs), ' motifs \n')
  
  pwm.corr <- compare_motifs(motifs, method = "PCC", min.mean.ic = 0,
                                score.strat = "a.mean")
  
  newName = motif.tf$motifs.new[match(rownames(pwm.corr), motif.tf$motifs)]
  rownames(pwm.corr) = newName
  colnames(pwm.corr) = newName
  
  comparisons <- 1 - pwm.corr
  dd <- as.dist(comparisons)
  
  # Hierarchical clustering using Complete Linkage
  hc <- hclust(dd, method = "ward.D2" )
  
  # Plot the obtained dendrogram
  #plot(hc, cex = 0.6, hang = -1)
  #sub_grp <- cutree(hc, h = 0.1)
  pdfname = paste0(resDir, "/pwm_celegans_similarity_clustering.pdf")
  pdf(pdfname, width=20, height = 30)
  par(cex =0.5, mar = c(3,0.8,2,5)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  #plot(hc, cex = 0.5, hang = -1)
  plot(as.dendrogram(hc), cex=0.5, horiz=TRUE)
  abline(v = c(0.05, 0.1, 0.15), col = c('blue', 'red', 'green'))
  #rect.hclust(hc, h = hc.cutoff, border="darkred")
  #groups <- 
  length(unique(cutree(hc, h = 0.05)))
  length(unique(cutree(hc, h = 0.1)))
  length(unique(cutree(hc, h = 0.15)))
  length(unique(cutree(hc, h = 0.2)))
  length(unique(cutree(hc, h = 0.25)))
  
  dev.off()
  
  change.pwm.logo.names = FALSE
  if(change.pwm.logo.names){
    logoDir = '../data/motifs_tfs/pwm_logos'
    logo.file = list.files(path = logoDir, pattern = '*.pdf', full.names = TRUE)
    for(n in 1:nrow(motif.tf))
    {
      # n = 1 
      cmd = paste0('mv ', logo.file[grep(motif.tf$motifs[n], logo.file)], ' ', logoDir, '/',  motif.tf$motifs.new[n], '.pdf')
      system(cmd)
    }
    
  }
  #fviz_nbclust(diss = comparisons, FUN = hcut, method = "wss")
  #fviz_nbclust(df, FUN = hcut, method = "silhouette")
  
  ##########################################
  # merge motifs using height = 0.1 and change motif names
  ##########################################
  hc.cutoff = 0.1
  
  groups <- cutree(hc, h = hc.cutoff)
  motif.tf = data.frame(motif.tf, group = groups, stringsAsFactors = FALSE)
  motif.tf$names = NA
  for(nn in unique(motif.tf$group))
  {
    # nn = 5
    kk = which(motif.tf$group == nn)
    motif.tf$names[kk] = paste0(paste0(unique(motif.tf$tfs.new[kk]), collapse = '_'), '.M', nn)
    
  }
  
  # save motif-to-tf mapping with redundancy removal information
  saveRDS(motif.tf, file = '../data/motifs_tfs/motif_tf_mapping.rds') 
  
}

##########################################
# after running FIMO, make motif occurrency matrix  
##########################################
make.motif.oc.matrix.from.fimo.output = function()
{
  library(data.table)
  motif.tf = readRDS( '../data/motifs_tfs/motif_tf_mapping.rds')
  fimo.out = '/Volumes/groups/cochella/jiwang/Databases/motifs_TFs/scATAC_earlyEmbryo_cel/fimo_out/fimo.tsv'
  fimo = fread(fimo.out, header = TRUE)
  motif.oc = table(fimo$motif_id, fimo$sequence_name, useNA = 'ifany')
  motif.oc = t(motif.oc)
  
  print(head(rownames(motif.oc)))
  
  ##########################################
  # associate the scanned regions with gene
  # here we restrict the assignment to protein-coding genes using ChIPpeakAnno 
  # https://www.bioconductor.org/packages/release/bioc/vignettes/ChIPpeakAnno/inst/doc/pipeline.html
  ##########################################
  assign.regions.to.genes = TRUE
  if(assign.regions.to.genes){
    ## loading packages
    #library(ChIPseeker)
    library(ChIPpeakAnno)
    library(GenomicFeatures)
    library(GenomicRanges)
    library(plyranges)
    
    annot = read.csv(file = '/Volumes/groups/cochella/jiwang/annotations/BioMart_WBcel235_noFilters.csv') # all annotation included
    annot = annot[grep('protein_coding', annot$Gene.type), ] # keep only protein-coding genes
    #library(TxDb.Celegans.UCSC.ce11.ensGene)
    #txdb <- TxDb.Celegans.UCSC.ce11.ensGene
    annot = annot[match(unique(annot$Gene.name), annot$Gene.name), ]
    xx = data.frame(seqnames = paste0('chr', annot$Chromosome.scaffold.name), start = annot$Gene.start..bp., end= annot$Gene.end..bp., 
                    gene_id = as.character(annot$Gene.name),
                    strand=c("."), score=0, stringsAsFactors = FALSE)
    xx = makeGRangesFromDataFrame(xx, keep.extra.columns = TRUE)
    names(xx) = as.character(annot$Gene.name)
    
    get.peak.coord = function(x){
      x = unlist(strsplit(as.character(x), '[:]'))
      chr = x[1]
      x = unlist(strsplit(as.character(x[2]), '-'))
      return(c(chr, x))
    }
    peaks = rownames(motif.oc)
    peaks = t(sapply(peaks, get.peak.coord))
    colnames(peaks) = c('chr', 'start', 'end')
    peaks = data.frame(peaks, strand=c("."), score=0, stringsAsFactors = FALSE)
    peaks = makeGRangesFromDataFrame(peaks)
    #peaks = data.frame(, gsub('^*:')))
    #peak <- readPeakFile(peakfile = peak.file, as = 'GRanges')
    peakAnno = annotatePeakInBatch(peaks, 
                        AnnotationData=xx, 
                        output='nearestLocation', 
                        bindingRegion=c(-2000, 500))
    #peakAnno <- annotatePeak(peak = peaks, tssRegion=c(-2000, 2000), level = 'gene', TxDb=txdb)
    assign = as.data.frame(peakAnno)
    
    #mm = match(assign$geneId, annot$wormbase.id)
    assign$genes = assign$feature
    rownames(assign) = assign$peak
    
    # keep scanned regions with mapped protein coding genes
    jj = which(!is.na(assign$genes))
    assign = assign[jj, ]
    motif.oc = motif.oc[jj, ]
    
    gene.uniq = unique(assign$genes)
    cat(nrow(motif.oc),  'scanned regions assigned to :', length(gene.uniq), ' genes \n')
    
    mocc = motif.oc[match(gene.uniq, assign$genes), ]
    rownames(mocc) = gene.uniq
    
    ## motif occurrency gene * motif
    for(n in 1:nrow(mocc))
    {
      kk = which(assign$genes == rownames(mocc)[n])
      if(length(kk)>1) {
        cat(n, '\n')
        mocc[n, ] = apply(motif.oc[kk, ], 2, sum) 
      }
    }
    
    motif.oc = mocc
    remove(mocc)
    
  }else{
    load(file = '../data/Hashimsholy_et_al/annotMapping_ensID_Wormbase_GeneName.Rdata')
    
    mm = match(rownames(motif.oc), geneMapping$Wormbase)
    rownames(motif.oc) = geneMapping$Gene.name[mm]
    #kk = match(rownames(motif.oc, ))
    
    ss1 = apply(motif.oc, 1, sum)
    cat(length(which(ss1 == 0)), 'genes without scanned motifs \n')
    ss2 = apply(motif.oc, 2, sum)
    
  }
  
  ##########################################
  # remove motif redundancy by merging occurrence for motifs in the same cluster
  ##########################################
  #mm = match(colnames(motif.oc), motif.tf$motifs)
  #colnames(motif.oc) = motif.tf$motifs.new[mm]
  names = unique(motif.tf$names[match(colnames(motif.oc), motif.tf$motifs)])
  xx = matrix(0, ncol = length(names), nrow = nrow(motif.oc))
  rownames(xx) = rownames(motif.oc)
  colnames(xx) = names
  
  ## to get motif occurrency gene * non-redundant-motif
  ## among the redundant motifs in the same cluster, the one with median of total occurrence is chosen.
  for(n in 1:ncol(xx))
  {
    # n = 3
    cat(n, ' -- ', colnames(xx)[n],  '\n')
    mtf = motif.tf$motifs[which(motif.tf$names == colnames(xx)[n])]
    
    kk = match(mtf, colnames(motif.oc))
    kk = kk[!is.na(kk)]
    
    if(length(kk) == 0){
      cat('Error : no motif found \n')
      
    }else{
      if(length(kk) == 1){
        
        xx[,n] = motif.oc[, kk]
        
      }else{
        cat('>>>>>>>>>',  length(kk), 'columns found \n')
        ss = apply(motif.oc[,kk], 2, sum)
        ss.o = ss[order(ss)]
        kk.sel = kk[which(ss == ss.o[ceiling(length(ss)/2)])]
        xx[,n] = motif.oc[, kk.sel[1]]
      }
    }
  }
  
  motif.oc = xx;
  remove(xx)
  saveRDS(motif.oc, file = '../data/motifs_tfs/motif_oc_scATACpeaks_all_proteinCodingGenes.rds')
  
}



########################################################
########################################################
# Section : package installation
# 
########################################################
########################################################
install.magic = FALSE
if(install.magic){
  install.packages("Rmagic")
  system('which python')
  system('python --version')
  system('pip install --user magic-impute') 
  # at the same time the magic-impute was installed in default python 
  # /usr/local/bin/python
  
}

##########################################
# test MAGIC and it works for the example
# original code from https://github.com/KrishnaswamyLab/MAGIC
##########################################
Test.Rmagic = FALSE
if(Test.Rmagic){
  library(Rmagic)
  library(ggplot2)
  data(magic_testdata)
  
  ss = apply(magic_testdata, 2, sum)
  magic_testdata = magic_testdata[, ss>0]
  MAGIC_data <- Rmagic::magic(magic_testdata, genes=c("VIM", "CDH1", "ZEB1"), verbose = 1)
  
  ggplot(MAGIC_data) +
    geom_point(aes(x=VIM, y=CDH1, color=ZEB1))
  
}

Test.phateR = FALSE
if(Test.phateR){
  ##########################################
  # test PHATE from https://github.com/KrishnaswamyLab/phateR
  # FAQ
  # 
  # Should genes (features) by rows or columns?
  #   
  #   To be consistent with common dimensionality reductions such as PCA (stats::prcomp) and t-SNE (Rtsne::Rtsne), we require that cells (observations) be rows and genes (features) be columns of your input data.
  # 
  # Can I run PHATE with Seurat?
  #   
  #   PHATE was removed from Seurat in version 3. You can install a version of Seurat with RunPHATE included by following the instructions at https://github.com/satijalab/seurat/pull/1172#issuecomment-564782167.
  # 
  # I have installed PHATE in Python, but phateR says it is not installed!
  #   
  #   Check your reticulate::py_discover_config("phate") and compare it to the version of Python in which you installed PHATE (run which python and which pip in a terminal.) Chances are reticulate can’t find the right version of Python; you can fix this by adding the following line to your ~/.Renviron:
  #   
  #   PATH=/path/to/my/python
  # 
  # You can read more about Renviron at https://CRAN.R-project.org/package=startup/vignettes/startup-intro.html.
  # Help
  # 
  # Please let us know of any issues at the GitHub repository. If you have any questions or require assistance using PHATE, please read the documentation at https://CRAN.R-project.org/package=phateR/phateR.pdf or by running help(phateR::phate) or contact us at https://krishnaswamylab.org/get-help.
  ##########################################
  library(phateR)
  #> Loading required package: Matrix
  data(tree.data)
  plot(prcomp(tree.data$data)$x, col=tree.data$branches)
  # runs phate
  tree.phate <- phate(tree.data$data)
  summary(tree.phate)
  #> PHATE embedding
  #> k = 5, alpha = 40, t = auto
  #> Data: (3000, 100)
  #> Embedding: (3000, 2)
  
  # plot embedding
  palette(rainbow(10))
  plot(tree.phate, col = tree.data$branches)
  
}
