##################################################
##################################################
## Project: collections of functions for single-cell RNAseq data analysis
## Script purpose:
## Usage example: 
## Author: Jingkui Wang (jingkui.wang@imp.ac.at)
## Date of creation: Thu Feb 22 14:51:03 2018
##################################################
##################################################
process.countTable = function(all, design, additional.filter = NULL)
{
  newall = data.frame(as.character(all[,1]), stringsAsFactors = FALSE)
  
  for(n in 1:nrow(design))
  {
    #n = 1;
    ## found the matching ID in design matrix
    if(!is.null(additional.filter)){
      jj = intersect(grep(design$SampleID[n], colnames(all)), grep(additional.filter, colnames(all)));
    }else{
      jj = grep(design$SampleID[n], colnames(all));
    }
    
    ## deal with several mapping 
    if(length(jj)==1) {
      #index = c(index,jj)
      newall = data.frame(newall, all[, jj])
    }else{
      cat(length(jj), " samples found for ID", design$SampleID[n], "\n")
      cat("start to merge those samples considering them as technical replicates...\n")
      newall = data.frame(newall, apply(as.matrix(all[, jj]), 1, sum))
    }
  }
  
  colnames(newall)[1] = "gene";
  jj = which(colnames(design)=="SampleID")
  o1 = c(setdiff(c(1:ncol(design)), jj),jj)
  #design = design
  colnames(newall)[-1] = apply(design[, o1], 1, function(x) paste0(x, collapse = "_"))
  
  #if(time.series){
  #  colnames(newall)[-1] = paste0(design$stage, "_", design$treatment, "_", design$SampleID)
  #}else{
  #}
  
  return(newall)
}

aggrate.nf.QCs = function(dirs.all, modules2cat = c("star", "featureCounts"))
{
  # dir = '../../../Ariane/R7116_R7130_scrnaseq/results_all/multiqc_data_1'
  # modules2cat = c("star", "featureCounts")
  
  res = NULL;
  for(fold in dirs.all){
    if(dir.exists(fold)){
      
      yy = NULL;
      for (mod in modules2cat)
      {
        # fold = dirs.all[1]; mod = modules2cat[1];
        files = list.files(path = fold, pattern = "*.txt", full.names = TRUE)
        files = files[grep(mod, files)]
        if(length(files)==0){
          cat("Error -- no file found for: ",  mod, "in dir --", fold, "\n")
        }else{
          for(f in files){
            if(is.null(yy)){
              yy = data.frame(read.delim(f, sep = "\t", header = TRUE), stringsAsFactors = FALSE)
            }else{
              test = data.frame(read.delim(f, sep = "\t", header = TRUE), stringsAsFactors = FALSE)
              yy = merge(yy, test, by = "Sample", suffixes = c("",".1"))
            }
          } 
        }
      } 
      
    }else{
      cat('Error -- ', fold, "not exist--\n")
    }
    
    if(is.null(res)){
      res = yy;
    }else{
      res = rbind(res, yy)
    }
  }
  
  return(res)
}

##########################################
# Convert Ensem ID to gene names
# a) keep only rows for mapped transripts
# b) map the gene names to the gene symbols
##########################################
convertGeneNames.forCountTable = function(aa, discard.gene.with.zero.reads = TRUE, annot = NULL)
{
  
  aa = aa[ grep("^__", aa$gene, invert = TRUE), ] ## remove features that are not well aligned
  
  ## load annotation and change the gene names
  if(is.null(annot)){
    annot = read.csv(file = 
                       "/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/annotations/BioMart_WBcel235_noFilters.csv", 
                     header = TRUE)
  }else{
    
  }
  
  gg.Mt = annot$Gene.name[which(annot$Chromosome.scaffold.name=="MtDNA")]
  gg.ribo = annot$Gene.name[which(annot$Gene.type=="rRNA")]
  
  ggs = as.character(aa$gene);
  mm = match(aa$gene, annot$Gene.stable.ID)
  jj = which(!is.na(mm)==TRUE & annot$Gene.name[mm] != "");
  ggs[jj] = as.character(annot$Gene.name[mm[jj]]);
  
  counts = aa
  rownames(counts) <- aa[, 1]
  counts <- as.matrix(counts[,-1])
  
  if(discard.gene.with.zero.reads){
    cat("-- discard genes with zero reads -- \n")
    ss = apply(counts, 1, sum)
    keep.genes = which(ss>0)
    counts = counts[keep.genes, ]
    ggs = ggs[keep.genes]
    ggs.unique = unique(ggs)
    
    rownames(counts) = ggs
  }
  
  return(counts)
}

find.particular.geneSet = function(geneSet = "Mt", annot = NULL)
{
  if(is.null(annot)){
    annot = read.csv(file = 
                       "/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/annotations/BioMart_WBcel235_noFilters.csv", 
                     header = TRUE)
  }
  
  if(geneSet == "Mt"){
    return(annot$Gene.name[which(annot$Chromosome.scaffold.name=="MtDNA")])
  }else{
    if(geneSet == "ribo"){
      return(annot$Gene.name[which(annot$Gene.type=="rRNA")])
    }else{
      stop("no geneSet found for :", geneSet, "\n")
    }
  }
    
}

##########################################
# function for collapse technical replicates in the lane level
# 
##########################################
compare.techinical.replicates = function(design.tr, counts.tr, filter.cells = FALSE, filter.genes = TRUE, check.correlations = TRUE)
{
  sinfos.uniq = unique(design.tr$seqInfos)
  
  ## add some new features for design for quality controls
  design.tr$log10_Total = log10(design.tr$total_reads)
  design.tr$percent_rRNA = design.tr$rRNA / design.tr$Total
  
  ## make SCE object and remove genes with zero reads detected
  sce <- SingleCellExperiment(assays = list(counts = counts.tr), 
                              colData = as.data.frame(design.tr), 
                              rowData = data.frame(gene_names = rownames(counts.tr), feature_symbol = rownames(counts.tr)))
  
  #is.spike <- grepl("^ERCC", rownames(sce))
  is.mito <- rownames(sce) %in% gg.Mt;
  is.ribo <- rownames(sce) %in% gg.ribo;
  summary(is.mito)
  summary(is.ribo)
  
  sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito, Ribo=is.ribo))
  
  # check some general statistics for each request and lane
  #cols = c(rep('gray', 2), 'red', 'darkblue', 'darkred', 'blue', 'red')
  ps1 = plotColData(sce, y = "log10_Total", x = "seqInfos") + ggtitle("total nb of reads")
  ps2 = plotColData(sce, y="uniquely_mapped_percent", x="seqInfos") + ggtitle("% of uniquely mapped ")
  ps3 = plotColData(sce, y="percent_assigned", x="seqInfos") + ggtitle("% of assigned")
  ps4 = plotColData(sce, y="pct_counts_Ribo", x="seqInfos") + ggtitle("% of rRNA contamination")
  ps5 = plotColData(sce, y="pct_counts_Mt", x="seqInfos") + ggtitle("% of Mt")
  
  ps6 = plotColData(sce, y="log10_total_counts", x="seqInfos") + ggtitle("total nb of reads mapped to transcripts")
  ps7 = plotColData(sce, y="total_features_by_counts", x="seqInfos") + ggtitle("total nb of genes")
  
  ps8 = plotColData(sce,
              x = "log10_total_counts",
              y = "log10_total_features_by_counts",
              #colour_by = "percent_mapped",
              colour_by = "seqInfos",
              size_by = "pct_counts_Mt"
  ) + scale_x_continuous(limits=c(4, 7)) +
    scale_y_continuous(limits = c(2.5, 4.1)) +
    geom_hline(yintercept=log10(c(500, 1000, 5000)) , linetype="dashed", color = "darkgray", size=0.5) +
    geom_vline(xintercept = c(4:6), linetype="dotted", color = "black", size=0.5)
  
  # filter cell and genes
  if(filter.cells){
    threshod.total.counts.per.cell = 10^4
    threshod.nb.detected.genes.per.cell = 1000;
    
    filter_by_total_counts <- (sce$total_counts > threshod.total.counts.per.cell)
    table(filter_by_total_counts)
    filter_by_expr_features <- (sce$total_features_by_counts > threshod.nb.detected.genes.per.cell)
    table(filter_by_expr_features)
    filter_by_MT = sce$pct_counts_Mt < 7.5
    table(filter_by_MT)
    
    sce$use <- (
      filter_by_expr_features & # sufficient features (genes)
        filter_by_total_counts & # sufficient molecules counted
        # filter_by_ERCC & # sufficient endogenous RNA
        filter_by_MT # remove cells with unusual number of reads in MT genes
    )
    table(sce$use)
    sce = sce[, sce$use]
  }
  if(filter.genes){
    num.cells <- nexprs(sce, byrow=TRUE)
    ave.counts <- calcAverage(sce)
    genes.to.keep <- num.cells > 5 & ave.counts >= 10  & ave.counts <10^6  # detected in >= 2 cells, ave.counts >=5 but not too high
    summary(genes.to.keep)
    # remove mt and ribo genes
    genes.to.keep = genes.to.keep & ! rownames(sce) %in% gg.Mt & ! rownames(sce) %in% gg.ribo
    summary(genes.to.keep)
    
    sce <- sce[genes.to.keep, ]
    
  }
  
  
  for(kk in 1:8)
  {
    eval(parse(text = paste0("plot(ps", kk, ")")))
  }
  
  # select cells having technical replicates and normalize them  
  sce.qc = sce
  reducedDim(sce.qc) <- NULL
  endog_genes <- !rowData(sce.qc)$is_feature_control
  
  
  Methods.Normalization = c("cpm", "DESeq2", "scran")
  for(method in Methods.Normalization)
  {
    
    set.seed(1234567)
    cat('testing normalization method -- ', method, "\n")
    
    main = paste0(method, " normalization");
        
    if(method == "cpm") { ### cpm
      assay(sce.qc, "logcounts") <- log2(calculateCPM(sce.qc, use_size_factors = FALSE) + 1)
    }
    
    if(method == "DESeq2"){
      source("scRNAseq_functions.R")
      sizeFactors(sce.qc) = calculate.sizeFactors.DESeq2(counts(sce.qc))
      sce.qc <- normalize(sce.qc, exprs_values = "counts", return_log = TRUE)
    }
    
    if(method == "scran"){
      ## scran normalization (not working here, because negative scaling factor found)
      qclust <- quickCluster(sce.qc, min.size = 50)
      sce.qc <- computeSumFactors(sce.qc, clusters = qclust)
      sce.qc <- normalize(sce.qc, exprs_values = "counts", return_log = TRUE)
    }
    
    ps.norms = scater::plotPCA(
      sce.qc[endog_genes, ],
      run_args = list(exprs_values = "logcounts"), 
      size_by = "total_counts",
      #size_by = "total_features_by_counts",
      colour_by = "seqInfos"
    ) + ggtitle(paste0("PCA --", method))
    plot(ps.norms)
    
    ## check the correction of the same cells from different technical replicates
    bcs = unique(sce.qc$barcodes)
    
    if(check.correlations){
      correlations = c()
      
      #par(mfrow=c(2,3)) 
      for(n in 1:length(bcs))
      {
        # n = 1
        kk = which(sce.qc$barcodes == bcs[n])
        xx = as.data.frame(logcounts(sce.qc[, kk[match(sinfos.uniq, sce.qc$seqInfos[kk])]]))
        ss = apply(xx, 1, sum)
        xx = xx[ss>0, ]
        colnames(xx) = paste0(sinfos.uniq, ".", method)
        
        if(n <= 10) pairs(xx, lower.panel=NULL, upper.panel=panel.fitting)
        
        if(length(sinfos.uniq) == 2) correlations = c(correlations, cor(xx[,1], xx[, 2]))
        if(length(sinfos.uniq) == 3) correlations = rbind(correlations, c(cor(xx[, 1], xx[, 2]), cor(xx[, 1], xx[, 3]), cor(xx[, 2], xx[, 3])))
      }
      
      if(length(sinfos.uniq) == 2) {
        names(correlations) = paste0("cor_", sinfos.uniq, collapse = "_")
        hist(correlations, breaks = 10)
      }
      if(length(sinfos.uniq) == 3) {
        colnames(correlations) = paste0("cor_", sinfos.uniq[c(1,2,3)], sinfos.uniq[c(2,3,1)])
        pairs(correlations, lower.panel=NULL, upper.panel=panel.fitting)
      }
      #colnames(correlations) = c('rep0.vs.hiseq.rep1', 'rep0.vs.hiseq.rep2', 'hiseq.rep1.vs.hiseq.rep_2')
    }
    
  }
  
}

merge.techinical.replicates = function(design, counts, compare.replicates = FALSE,
                                       sampleInfos.techRep = list(c("R7130_HHG5KBGX9_1", "R7130_HLWTCBGX9_1"),
                                                                  c("R7130_HHGHNBGX9_1", "R7130_CCYTEANXX_4", "R7133_CD2GTANXX_5")))
{
  
  for(n in 1:length(sampleInfos.techRep))
  {
    # check if technical replicates can be found
    techRep = sampleInfos.techRep[[n]]
    mm = match(techRep, unique(design$seqInfos))
    if(any(is.na(mm))) stop("Missed technical replicates : ", sampleInfos.techRep[which(is.na(mm))], "\n")
    
    sels = match(design$seqInfos, techRep)
    sels = which(!is.na(sels))
    
    design.keep = design[-sels, ]
    counts.keep = counts[, -sels]
    
    design.sels = design[sels, ]
    counts.sels = counts[, sels]
    
    if (compare.replicates){
    cat("-- start to compare technical replicates for lanes :", techRep, "\n")
    compare.techinical.replicates(design.tr=design.sels, counts.tr=counts.sels)
    }
    
    cat("-- start to merge technical replicates for lanes :", techRep, "\n")
    
    bcs = unique(design.sels$barcodes)
    
    design.merged = design.sels[match(bcs, design.sels$barcodes), ]
    design.merged$seqInfos = paste0(techRep[1], "_merged")
    
    # double chekc the barcode order
    for(index.bc in 1:length(bcs)){
      if(bcs[index.bc] != design.merged$barcodes[index.bc]) {
        cat("barcode order is wrong \n")
      }
    }
    
    merged = matrix(0, nrow = nrow(counts.sels), ncol = length(bcs))
    
    for(m in 1:length(techRep))
    {
      jj.trep = which(design.sels$seqInfos == techRep[m])
      counts.trep = counts.sels[, jj.trep]
      colnames(counts.trep) = design.sels$barcodes[jj.trep]
      counts.trep = counts.trep[, match(bcs, colnames(counts.trep))]
      counts.trep = as.matrix(counts.trep)
      counts.trep[which(is.na(counts.trep))] = 0
      merged = merged + counts.trep
    }
    colnames(merged) = design.merged$samples
    
    design = rbind(design.keep, design.merged)
    counts = cbind(counts.keep, merged)
  }
  
  return(list(design = design, counts = counts))
  
}


panel.cor <- function(x, y, digits=2, prefix="", cex.cor) 
{
  usr <- par("usr"); on.exit(par(usr)) 
  par(usr = c(0, 1, 0, 1)) 
  r <- abs(cor(x, y)) 
  txt <- format(c(r, 0.123456789), digits=digits)[1] 
  txt <- paste(prefix, txt, sep="") 
  if(missing(cex.cor)) cex <- 0.8/strwidth(txt) 
  
  test <- cor.test(x,y) 
  # borrowed from printCoefmat
  Signif <- symnum(test$p.value, corr = FALSE, na = FALSE, 
                   cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " ")) 
  
  text(0.5, 0.5, txt, cex = cex * r) 
  text(.8, .8, Signif, cex=cex, col=2) 
}

panel.fitting = function (x, y, bg = NA, pch = par("pch"), cex = 0.3, col='black') 
{
  #x = yy[,1];y=yy[,2];
  #kk = which(x>0 & y>0); x=x[kk];y=y[kk]
  lims = range(c(x, y), na.rm = TRUE)
  points(x, y, pch = 1, col = col, cex = cex, xlim=lims, ylim=lims)
  abline(0, 1, lwd=1.5, col='red')
  R = cor(x, y, use="na.or.complete", method='pearson')
  text(lims[2]*0.2, lims[2]*0.9, paste0('R = ', signif(R, d=2)), cex=1., col='red')
  #jj = which(!is.na(x) & !is.na(y))
  #fit = lm(y[jj] ~ x[jj])
  #slope=summary(fit)$coefficients[1]
  #slope = fit$coefficients[2]
  #intercept = fit$coefficients[1]
  #pval=summary(fit)$coefficients[4]
  #abline(intercept, slope, lwd=1.2, col='darkblue', lty=3)
  #text(lims[2]*0.1, lims[2]*0.7, paste0('slop = ', signif(slope, d=2)), cex=1., col='blue')
  #text(lims[2]*0.1, lims[2]*0.6, paste0('pval = ', signif(pval, d=2)), cex=1., col='blue')
  #ok <- is.finite(x) & is.finite(y)
  #if (any(ok)) 
  #lines(stats::lowess(x[ok], y[ok], f = span, iter = iter), 
  #        col = col.smooth, ...)
}

plot.pair.comparison.plot = function(xx, linear.scale = TRUE, main = "pairwise.comparision"){
  yy = as.matrix(xx)
  if(linear.scale){
    yy[which(yy==0)] = NA;
    yy = log2(yy)
  }
  pairs(yy, lower.panel=NULL, upper.panel=panel.fitting, main = main)
  
}

########################################################
########################################################
# Section : integrate the FACS information into the metadata
# 
########################################################
########################################################
Integrate.FACS.Information = function(sce, processing.FACS.Info = TRUE)
{
  
  find_flowcell_lane = function(x)
  {
    x = gsub(".csv", "", x)  
    x = basename(x)
    x = unlist(strsplit(as.character(x), "_"))
    x = paste0(x[c((length(x)-1), length(x))], collapse = "_")
    return(x)
  }
  
  if(processing.FACS.Info){
    ##########################################
    # collecting the facs infos from tables
    ##########################################
    path2FACS = "/Volumes/groups/cochella/git_aleks_jingkui/scRNAseq_MS_lineage/data/facs_data/FACS_indexData/"
    ff1 = paste0(path2FACS, "/barcodes_to_wellNames_96w_robot_test_plus_barcodes_CCVBPANXX_1.csv")
    facs = read.csv(ff1, header = TRUE, row.names = 1)
    facs = data.frame(facs, stringsAsFactors = FALSE)
    flane = find_flowcell_lane(ff1)
    facs$flowcell_lane = flane
    
    keep = facs
    
    ffs = list.files(path = path2FACS, pattern = "*.csv", full.names = TRUE)
    ffs = ffs[grep(flane, ffs, invert = TRUE)]
    
    jj = grep("barcodes_to_wellNames", ffs)
    wells = ffs[jj]
    facs = ffs[-jj]
    
    wells.infos = list()
    for(n in 1:length(wells))
    {
      xx = read.csv(wells[n], header = TRUE)
      xx = data.frame(xx$well_name, xx$original, stringsAsFactors = FALSE)
      colnames(xx) = c("index.well", "bc")
      wells.infos[[n]] = xx
      if(grepl("384w", wells[n])){
        names(wells.infos)[n] = "all"
      }else{
        names(wells.infos)[n] = find_flowcell_lane(wells[n])
      }  
    }
    #well2barcodes = read.csv(well2barcodes, header = TRUE)
    #wells = data.frame(well2barcodes$well_name, well2barcodes$original, stringsAsFactors = FALSE)
    keep2 = c()
    for(n in 1:length(facs))
    {
      # n = 2
      flane = find_flowcell_lane(facs[n])
      yy = read.csv(facs[n], header = TRUE)
      if(flane == "CCVTBANXX_8"){
        mapping = wells.infos$CCVTBANXX_8
        yy = data.frame(yy[, c(1:8)], rep(NA, nrow(yy)), rep(NA, nrow(yy)), yy$Index)
        colnames(yy)[c(9:11)] = colnames(keep)[9:11]
      }else{
        mapping = wells.infos$all
      }
      yy = data.frame(yy, mapping[match(yy$Index, mapping$index.well),], stringsAsFactors = FALSE)
      colnames(yy)[c(12:13)] = colnames(keep)[12:13]
      yy$flowcell_lane = flane
      keep2 = rbind(keep2, yy)
    }
    
    keep = rbind(keep, keep2)
    colnames(keep)[11:12] = c("Index.well", "Index.well.new")
    keep$flowcell_lane_bc = paste0(keep$flowcell_lane, "_", keep$barcode)
    
    #save(keep, file = paste0(RdataDir, "merged_FACS_information_all.Rdata"))
    ##########################################
    # integrate the facs info into the sce object 
    ##########################################
    nb.cells = rep(NA, ncol(sce))
    FSC = rep(NA, ncol(sce))
    BSC = rep(NA, ncol(sce))
    Index.well = rep(NA, ncol(sce))
    GFP = rep(NA, ncol(sce))
    
    for(n in 1:ncol(sce))
    {
      kk = which(keep$flowcell_lane_bc == paste0(sce$flowcell.lane[n], "_", sce$barcodes[n]))
      cat(n, " : ", length(kk), "\n" )
      if(length(kk)>0){
        nb.cells[n] = length(kk)
        Index.well[n] = unique(keep$Index.well.new[kk])
        FSC[n] = mean(keep$FSC.A[kk])
        BSC[n] = mean(keep$BSC.A[kk])
        GFP[n] = mean(keep$FITC.A.Compensated[kk])
      }
    }
    
    sce$nb.cells = nb.cells
    sce$index.well = Index.well
    sce$FSC = FSC
    sce$BSC = BSC
    sce$GFP = GFP
    # save(sce, file = paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_facsInfos.Rdata'))
    
  }
  
  return(sce)
}

##########################################
# test umap parameters for BWM cells or general seurat object
##########################################
test.umap.params.for.BWM.cells = function(sub.obj, 
                                          pdfname = 'BWM_UMAP_explore_parameters.pdf',
                                          group.by = 'manual.annot.ids', with_legend = FALSE,
                                          weight.by.var = FALSE,
                                          nfeatures.sampling = c(3000, 5000),
                                          nb.pcs.sampling = c(5, 10, 30, 50), 
                                          n.neighbors.sampling = c(5, 10, 30, 50),
                                          min.dist.sampling = c(0.01, 0.1, 0.3, 0.5), 
                                          spread.sampling = 1
)
{
  if (length(dev.list()!=0)) {dev.off()}
  pdfname = paste0(resDir,'/', pdfname);while (!is.null(dev.list()))  dev.off();
  pdf(pdfname, width=18, height = 12)
  par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  
  for(nfeatures in nfeatures.sampling)
  {
    sub.obj <- FindVariableFeatures(sub.obj, selection.method = "vst", nfeatures = nfeatures, verbose = FALSE)
    sub.obj = ScaleData(sub.obj, features = rownames(sub.obj), verbose = FALSE)
    sub.obj <- RunPCA(object = sub.obj, features = VariableFeatures(sub.obj), verbose = FALSE, npcs = max(50, nb.pcs.sampling), 
                      weight.by.var = weight.by.var)
    
    for(nb.pcs in nb.pcs.sampling)
    {
      for(n.neighbors in n.neighbors.sampling)
      {
        for(min.dist in min.dist.sampling)
        {
          for(spread in spread.sampling){
            cat('nfeatures - ', nfeatures,  ', nb.pcs - ', nb.pcs, ', n.neighbors - ', n.neighbors, 
                ', min.dist - ', min.dist, ', spread - ', spread,  '\n')
            # nfeatures = 5000;
            # nb.pcs = 50 # nb of pcs depends on the considered clusters or ids 
            # n.neighbors = 50;
            # min.dist = 0.05; spread = 1;
            sub.obj <- RunUMAP(object = sub.obj, reduction = 'pca', reduction.name = "umap", dims = 1:nb.pcs, 
                               spread = spread, n.neighbors = n.neighbors, 
                               min.dist = min.dist, verbose = FALSE)
            
            if(with_legend){
              pp = DimPlot(sub.obj, group.by = group.by, reduction = 'umap', label = TRUE, label.size = 6, 
                           pt.size = 2, repel = TRUE) + 
                ggtitle(paste0('nfeatures - ', nfeatures,  ', nb.pcs - ', nb.pcs, ', n.neighbors - ', n.neighbors, 
                               ', min.dist - ', min.dist, ', spread - ', spread))
            }else{
              pp = DimPlot(sub.obj, group.by = group.by, reduction = 'umap', label = TRUE, label.size = 6, 
                           pt.size = 2, repel = TRUE) + 
                NoLegend() + 
                ggtitle(paste0('nfeatures - ', nfeatures,  ', nb.pcs - ', nb.pcs, ', n.neighbors - ', n.neighbors, 
                               ', min.dist - ', min.dist, ', spread - ', spread))
            }
            
            plot(pp)
          }
        }
      }
      
    }
    
  }
  
  
  dev.off()
  
}

##########################################################################
##########################################################################
# BIG SECTION:
# Script purpose: functions for scRNA-seq normalization
# Usage example: 
# Author: Jingkui Wang (jingkui.wang@imp.ac.at)
# Date of creation: Wed Nov 20 12:15:18 2019
##########################################################################
##########################################################################
########################################################
########################################################
# Section : Compare normalization methods for scRNA-seqe 
# to choose the good one
########################################################
########################################################
# several common functions for normalizations 
# from Hemberg lab
# hemberg-lab.github.io/scRNA.seq.course/cleaning-the-expression-matrix.html#normalization-theory
calc_cpm <- function (expr_mat, spikes = NULL) 
{
  norm_factor <- colSums(expr_mat[-spikes, ])
  return(t(t(expr_mat)/norm_factor)) * 10^6
}

Down_Sample_Matrix <- function (expr_mat)
{
  min_lib_size <- min(colSums(expr_mat))
  down_sample <- function(x) {
    prob <- min_lib_size/sum(x)
    return(unlist(lapply(x, function(y) {
      rbinom(1, y, prob)
    })))
  }
  down_sampled_mat <- apply(expr_mat, 2, down_sample)
  return(down_sampled_mat)
}

cal_uq_Hemberg = function (expr_mat, spikes = NULL) 
{
  UQ <- function(x) {
    quantile(x[x > 0], 0.75)
  }
  if(!is.null(spikes)){
    uq <- unlist(apply(expr_mat[-spikes, ], 2, UQ))
  }else{
    uq <- unlist(apply(expr_mat, 2, UQ))
  }
  
  norm_factor <- uq/median(uq)
  return(t(t(expr_mat)/norm_factor))
}

calculate.sizeFactors.DESeq2 = function(expr_mat)
{
  # expr_mat = counts(sce.qc)
  require('DESeq2')
  condition <- factor(rep("A", ncol(expr_mat)))
  dds <- DESeqDataSetFromMatrix(expr_mat, DataFrame(condition), design = ~ 1)
  dds <- estimateSizeFactors(dds)
  
  return(sizeFactors(dds))
  
}

test.normalization = function(sce, Methods.Normalization = c("cpm", "DESeq2", "scran", "seurat", "sctransform"), use.HVGs = TRUE)
{
  #Methods.Normalization = "DESeq2"
  for(method in Methods.Normalization)
  {
    sce.qc = sce
    sce.qc$library.size = apply(counts(sce.qc), 2, sum)
    set.seed(1234567)
    
    method = 'cpm'
    
    cat('test normalization method -- ', method, "\n")
    main = paste0(method);
    
    if(method == "raw") { # raw log counts
      assay(sce.qc, "logcounts") <- log2(counts(sce.qc) + 1)
    }
    
    if(method == "cpm") { ### cpm
      #assay(sce.qc, "logcounts") <- log2(calculateCPM(sce.qc, use_size_factors = FALSE) + 1)
      sce.qc = logNormCounts(sce.qc, size_factors = NULL, log = TRUE, pseudo_count=1, center_size_factors = TRUE)
    }
    if(method == "UQ"){
      logcounts(sce.qc) <- log2(cal_uq_Hemberg(counts(sce.qc)) + 1)
    }
    
    if(method == "DESeq2"){
      sizeFactors(sce.qc) = calculate.sizeFactors.DESeq2(counts(sce.qc))
      sce.qc <- normalize(sce.qc, exprs_values = "counts", return_log = TRUE)
    }
    if(method == "downsample") {
      assay(sce.qc, "logcounts") <- log2(Down_Sample_Matrix(counts(sce.qc)) + 1)
    }
    
    if(method == "scran"){
      #min.size = 100
      ## scran normalization (not working here, because negative scaling factor found)
      qclust <- quickCluster(sce.qc,  method = 'igraph')
      sce.qc <- computeSumFactors(sce.qc, clusters = qclust)
      sce.qc <- logNormCounts(sce.qc)
    }
    
    if(method == "TMM"|method == "DESeq2"|method == "UQ"|method == "scran"){
      summary(sizeFactors(sce.qc))
      range(sizeFactors(sce.qc))
      plot(sce.qc$library.size/1e6, sizeFactors(sce.qc), log="xy", main = paste0(method), 
           xlab="Library size (millions)", ylab="Size factor",
           pch=16)
      #legend("bottomright", col=c("black"), pch=16, cex=1.2, legend = "size factor from scran vs total library size")
    }
    
    if(method != 'Seurat' & method != "sctransform"){
      
      ## to have a fair comparison, here we also need use HVGs for PCAs, UMAP and t-SNE
      HVGs = find.HVGs(sce.qc, Norm.Vars.per.batch = FALSE, method = "scran", ntop = 2000)
      
      sce.qc = scater::runPCA(sce.qc)
      param.perplexity = 10; set.seed(100)
      sce.qc = scater::runTSNE(sce.qc,dimred = 'PCA', perplexity = param.perplexity)
      sce.qc = scater::runUMAP(sce.qc, dimred = 'PCA')
      #sce.zeisel <- runTSNE(sce.zeisel, dimred="PCA", perplexity=5)
      #runUMAP(sce.zeisel, dimred="PCA")
      scater::plotReducedDim(sce.qc, dimred = "PCA", colour_by = 'request')
      scater::plotReducedDim(sce.qc, dimred = "UMAP", colour_by = 'request')
      scater::plotReducedDim(sce.qc, dimred = "TSNE", colour_by = 'request')
      
      # p1 = scater::plotPCA(
      #   sce.qc[endog_genes, ],
      #   run_args = list(exprs_values = "logcounts"), 
      #   size_by = "total_counts",
      #   #size_by = "total_features_by_counts",
      #   colour_by = "request"
      # ) + ggtitle(paste0("PCA -- ", main))
      # 
      # 
      # p2 = plotTSNE(
      #   sce.qc[endog_genes, ],
      #   run_args = list(exprs_values = "logcounts", perplexity = param.perplexity), 
      #   size_by = "total_counts",
      #   #size_by = "total_features_by_counts",
      #   colour_by = "seqInfos"  
      # ) + ggtitle(paste0("tSNE - perplexity = ", param.perplexity, "--", main))
      # 
      # p3 = plotUMAP(
      #   sce.qc[endog_genes, ],
      #   run_args = list(exprs_values = "logcounts"), 
      #   size_by = "total_counts",
      #   #size_by = "total_features_by_counts",
      #   colour_by = "seqInfos"
      # ) + ggtitle(paste0("UMAP -- ", main))
      # 
      # plot(p1); plot(p2); plot(p3)
      
    }
    
    if(method == 'seurat'| method == 'sctransform'){
      ms0 = as.Seurat(sce, counts = 'counts', data = NULL, assay = "RNA")
      
      if(method == 'sctransform'){
        ms <- SCTransform(object = ms0) # new normalization from Seurat
        
        ms <- RunPCA(object = ms, verbose = FALSE)
        #ms <- FindNeighbors(object = ms, dims = 1:20)
        #ms <- FindClusters(object = ms)
        
        ElbowPlot(ms)
        
        ms <- RunUMAP(object = ms, reduction = 'pca', dims = 1:20, n.neighbors = 30, umap.method = "uwot",
                      metric = 'correlation', min.dist = 0.25)
        DimPlot(ms, reduction = "umap", group.by = 'request')
        
        
        ms <- RunUMAP(object = ms, reduction = 'pca', dims = 1:20, n.neighbors = 30, umap.method = "uwot",
                      metric = 'cosine', min.dist = 0.25)
        DimPlot(ms, reduction = "umap", group.by = 'request')
        
        
        ms <- RunTSNE(object = ms, reduction = 'pca', dims = 1:20)
        #ms <- RunUMAP(object = ms, reduction = 'MNN', dims = 1:20, n.neighbors = 30)
        DimPlot(ms, reduction = "tsne", group.by = 'request')
        
      }
      
      
      if(method == 'seurat'){
        #scale.factor = 10^4
        #scale.factor = mean(sce.qc$library.size)
        ms.logtransform <- NormalizeData(ms0, assay = "RNA", normalization.method = 'LogNormalize')
        ms.logtransform <- FindVariableFeatures(ms.logtransform, selection.method = "vst", nfeatures = 2000)
        ms.logtransform <- ScaleData(ms.logtransform, features = rownames(ms.logtransform))
        
        ms.logtransform <- RunPCA(object = ms.logtransform, verbose = FALSE)
        ElbowPlot(ms.logtransform)
        
        ms.logtransform <- RunUMAP(object = ms.logtransform, reduction = 'pca', dims = 1:20, n.neighbors = 30, min.dist = 0.25)
        DimPlot(ms.logtransform, reduction = "umap", group.by = 'request')
        
        ms.logtransform = RunTSNE(ms.logtransform, reduction = 'pca', dims = 1:10, tsne.method = 'Rtsne')
        DimPlot(ms.logtransform, reduction = "tsne", group.by = 'request')
      }
      
      ##########################################
      # Compare sctransform and standard normalization
      ##########################################
      Compare.two.norm.in.seurat.and.scran = FALSE
      if(Compare.two.norm.in.seurat){
        library(ggplot2)
        library(cowplot)
        p1 =DimPlot(ms.logtransform, reduction = "umap", group.by = 'request')
        p2 = DimPlot(ms, reduction = "umap", group.by = 'request')
        plot_grid(p1, p2)
      }
      
      Dissect.normalization.RunPCA.RunUMAP.in.Seurat.by.comparing.scater = FALSE
      if(Dissect.normalization.RunPCA.RunUMAP.in.Seurat.by.comparing.scater)
      {
        ## compare seurat normalizatioin vs 
        sce.qc = sce
        sce.qc$library.size = apply(counts(sce.qc), 2, sum)
        sce.qc = logNormCounts(sce.qc, size_factors = NULL, log = TRUE, pseudo_count=1, center_size_factors = FALSE)
        plot(sce.qc$library.size/1e6, sizeFactors(sce.qc), log="xy", main = paste0(method), 
             xlab="Library size (millions)", ylab="Size factor")
        
        HVGs = VariableFeatures(ms.logtransform)
        mm = match(VariableFeatures(ms.logtransform), rownames(sce.qc))
        sce.qc = scater::runPCA(sce.qc,  subset_row = mm, scale = TRUE)
        
        sce.qc = scater::runUMAP(sce.qc, dimred = 'PCA', n_dimred = 1:20,  n_neighbors = 20, metric = "cosine")
        scater::plotReducedDim(sce.qc, dimred = "UMAP", colour_by = 'request')
        param.perplexity = 50; set.seed(100)
        sce.qc = scater::runTSNE(sce.qc,dimred = 'PCA', perplexity = param.perplexity)
        scater::plotReducedDim(sce.qc, dimred = "TSNE", colour_by = 'request')
        
        kk = 1
        xx = log(counts(sce)[,kk]/sce.qc$library.size[kk]*10^4 +1)
        plot(xx, ms.logtransform@assays$RNA@data[, kk], log= 'xy')
        abline(0, 1, col='red')
        
        yy = log2((counts(sce)[,kk])/sizeFactors(sce.qc)[1] +1)
        plot(yy, logcounts(sce.qc)[, kk], log= 'xy')
        abline(0, 1, col='red')
        
        plot(logcounts(sce.qc)[, kk],  ms.logtransform@assays$RNA@data[, kk]/log(2))
        abline(0, 1, col='red')
        
        plot(reducedDim(sce.qc, "PCA")[, kk], Reductions(ms.logtransform, slot = 'pca')@cell.embeddings[,kk])
        
      }
      
    }
    
  }
}

compare.scran.seurat.sctransform = function(sce, using.HVGs = TRUE)
{
  sce.qc = sce
  sce.qc$library.size = apply(counts(sce.qc), 2, sum)
  set.seed(1234567)
  
  ms0 = as.Seurat(sce.qc, counts = 'counts', data = NULL, assay = "RNA")
  ss = sce.qc$library.size
  
  nfeatures = 3000
  n.neighbors = 30
  min.dist = 0.25
  nb.pcs = 20
  ### seurat normalization (cpm actually) 
  ms.logtransform <- NormalizeData(ms0, assay = "RNA", normalization.method = 'LogNormalize', scale.factor = mean(ss))
  ms.logtransform <- FindVariableFeatures(ms.logtransform, selection.method = "vst", nfeatures = nfeatures)
  
  ms.logtransform <- ScaleData(ms.logtransform, features = rownames(ms.logtransform))
  
  ms.logtransform <- RunPCA(object = ms.logtransform, verbose = FALSE)
  ElbowPlot(ms.logtransform)
  
  ms.logtransform <- RunUMAP(object = ms.logtransform, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
  p0 = DimPlot(ms.logtransform, reduction = "umap", group.by = 'request') + ggtitle("seurat")
  
  ### cpm normalized in scater
  cpm = logNormCounts(sce.qc, size_factors = NULL, log = TRUE, pseudo_count=1, center_size_factors = FALSE)
  plot(cpm$library.size/1e6, sizeFactors(cpm), log="xy", main = 'cpm', xlab="Library size (millions)", ylab="Size factor")
  
  ms1 = as.Seurat(cpm, counts = 'counts', data = 'logcounts', assay = "RNA")
  ms1 = FindVariableFeatures(ms1, selection.method = 'vst', nfeatures = nfeatures)
  ms1 = ScaleData(ms1, features = rownames(ms1))
  ms1 = RunPCA(ms1, verbose = FALSE)
  ms1 = RunUMAP(ms1, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
  p1 = DimPlot(ms1, reduction = "umap", group.by = 'request') + ggtitle('cpm')
  
  ### scran normalization
  qclust <- quickCluster(sce.qc)
  sce.norm <- computeSumFactors(sce.qc, clusters = qclust)
  sce.norm <- logNormCounts(sce.norm, log = TRUE, pseudo_count = 1)
  plot(sce.norm$library.size/1e6, sizeFactors(sce.norm), log="xy", xlab="Library size (millions)", ylab="Size factor")
  ms2 = as.Seurat(sce.norm, counts = 'counts', data = 'logcounts', assay = "RNA")
  ms2 = FindVariableFeatures(ms2, selection.method = 'vst', nfeatures = nfeatures)
  ms2 = ScaleData(ms2, features = rownames(ms1))
  ms2 = RunPCA(ms2, verbose = FALSE)
  ms2 = RunUMAP(ms2, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
  p2 = DimPlot(ms2, reduction = "umap", group.by = 'request') + ggtitle('scran')
  
  
  ms3 <- SCTransform(object = ms0, variable.features.n =nfeatures) # new normalization from Seurat
  ms3 <- RunPCA(object = ms3, verbose = FALSE)
  #ElbowPlot(ms)
  ms3 <- RunUMAP(object = ms3, reduction = 'pca', dims = 1:nb.pcs, n.neighbors = n.neighbors, min.dist = min.dist)
  p3 = DimPlot(ms3, reduction = "umap", group.by = 'request') + ggtitle('sctransform')
  
  plot_grid(p0, p1, p2, p3, nrow = 2)
  
}

normalized.counts.using.scran.old = function(sce)
{
  set.seed(1000)
  clusters <- quickCluster(sce, min.size = 100, method="igraph")
  table(clusters)
  
  sce <- computeSumFactors(sce, clusters = clusters)
  
  ## quick check for size factors calculated by scran
  summary(sizeFactors(sce))
  plot(sce$total_counts/1e6, sizeFactors(sce), log="xy",
       xlab="Library size (millions)", ylab="Size factor", pch=16)
  
  sce <- normalize(sce, exprs_values = "counts", return_log = TRUE)
  
}

########################################################
########################################################
# Section : HVGs selection
# 
########################################################
########################################################
find.HVGs = function(sce, Norm.Vars.per.batch = TRUE, method = "scran", ntop = NULL)
{
  if(method == "scran"){
    if(!Norm.Vars.per.batch){
      ## here we use batches as blockes, i.e. calculated mean and variances separately for each batches and then fit the trend
      ## In doing so, we implicitly assume that the trend is the same between batches
      ## https://bioconductor.org/packages/3.10/workflows/vignettes/simpleSingleCell/inst/doc/var.html#41_using_the_block=_argument
      fit <- trendVar(sce, block=sce$batches, parametric=TRUE, assay.type="logcounts", use.spikes=FALSE)
      dec <- decomposeVar(sce, fit)
      
      par(mfrow=c(1,2))
      cols = c(1:length(bc.uniq)) 
      nb.genes = min(table(sce$batches))
      
      matplot(fit$means[c(1:nb.genes), ], fit$vars[c(1:nb.genes), ], col=cols,
              xlab="Mean log-expression", ylab="Variance of log-expression")
      curve(fit$trend(x), add=TRUE, col="red")
      
      plot(dec$mean, dec$total, xlab="Mean log-expression", 
           ylab="Variance of log-expression", pch=16)
      curve(fit$trend(x), col="dodgerblue", add=TRUE)
      
      tmp.sce <- sce
      tmp.sce$log_size_factor <- log(sizeFactors(sce))
      plotColData(tmp.sce, x="batches", y="log_size_factor")
      #p2 <- plotColData(tmp.416B, x="Plate", y="log_size_factor_ERCC")
      #multiplot(p1, p2, cols=2)
      
      dec.sorted <- dec[order(dec$bio, decreasing=TRUE), ]
      head(dec.sorted)
      
      if(is.null(ntop)){
        # here HVGs selected with FDR < 0.01
        gene.chosen <- rownames(dec.sorted)[which(dec.sorted$FDR < 0.05)] 
      }else{
        gene.chosen = rownames(dec.sorted)[1:ntop] 
      }
      
    }else{
      ## here we search batch-specific features
      ## https://bioconductor.org/packages/3.10/workflows/vignettes/simpleSingleCell/inst/doc/var.html#42_fitting_batch-specific_trends
      sce.bc = sce
      ## recommendation in the function help: Centering should be performed by running multiBlockNorm before calling this function. 
      sce.bc <- multiBlockNorm(sce.bc, block = sce.bc$batches)
      
      par(mfrow=c(1,1))
      plot(sizeFactors(sce), sizeFactors(sce.bc), log='xy'); abline(0, 1, lwd =2, col = 'red') # did not change anything here, weired
      
      comb.out <- multiBlockVar(sce.bc, block=sce.bc$batches, assay.type="logcounts", trend.args=list(parametric=TRUE, use.spikes=FALSE))
      
      #comb.out = multiBlockVar(sce.bc, block = sce.bc$batches, trend.args = list(use.spikes = FALSE))
      head(comb.out[,1:6])
      
      par(mfrow=c(1, length(bc.uniq)))
      #is.spike <- isSpike(sce.416B.2)
      for (plate in unique(sce.bc$batches)) {
        cur.out <- comb.out$per.block[[plate]]
        plot(cur.out$mean, cur.out$total, pch=16, cex=0.6, xlab="Mean log-expression", 
             ylab="Variance of log-expression", main=plate, ylim = c(0, 25), xlim = c(0, 15))
        curve(metadata(cur.out)$trend(x), col="dodgerblue", lwd=2, add=TRUE)
        #points(cur.out$mean[is.spike], cur.out$total[is.spike], col="red", pch=16)
      }
      
      dec.bc.sorted <- comb.out[order(comb.out$bio, decreasing=TRUE), ]
      head(dec.bc.sorted)
      
      if(is.null(ntop)){
        #length(which(dec.sorted$bio>0))
        # here HVGs selected with FDR<0.01
        #gene.chosen.bc <- which(dec.bc.sorted$p.value < 0.05)
        #gene.chosen.bc <- which(dec.bc.sorted$FDR < 0.1)
        gene.chosen.bc = rownames(dec.bc.sorted)[which(dec.bc.sorted$bio>0)]
        #length(which(dec.sorted$bio>0)) 
      }else{
        gene.chosen.bc = rownames(dec.bc.sorted)[1:ntop]
      }
      
      #length(gene.chosen.bc)
      
      ### compare the block and batch-specific HVGs selection
      #cat(length(gene.chosen), length(gene.chosen.bc), length(intersect(gene.chosen, gene.chosen.bc)), "\n")
      #library(VennDiagram)
      #venn.diagram(
      #  x = list(gene.chosen, gene.chosen.bc),
      #  category.names = c("batch-block" , "batch-specific"), filename = paste0(resDir, "/batch_block_specific.png"), output = TRUE)
      gene.chosen = gene.chosen.bc
      
    }
    cat("nb of HGVs : ", length(gene.chosen), "\n")
    
  }
  
  if(method == "Brenneck")
  {
    library(M3Drop)
    if(!Norm.Vars.per.batch){
      expr_matrix =  exp(logcounts(sce))
      Brennecke_HVG <- BrenneckeGetVariableGenes(
        expr_mat = expr_matrix,
        spikes = NA,
        fdr = 0.3,
        minBiolDisp = 0.
      )
      
      gene.chosen = Brennecke_HVG
    }
    
  }
  
  return(gene.chosen)
}

########################################################
########################################################
# Section : funcitons for cell cycle scoring and correction
# At the end, Seurat is used 
########################################################
########################################################
find.cc.markers.homologues = function()
{
  detach("package:Seurat", unload=TRUE)
  require(Seurat)
  #s.genes = c("cdk-4", "evl-18") # from GO:1901987 http://amigo1.geneontology.org/cgi-bin/amigo/term-assoc.cgi?term=GO:1902808&speciesdb=all&taxid=6239
  #g2m.genes = xx # GO:1902751
  # A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
  # segregate this list into markers of G2/M phase and markers of S phase
  homologues = read.delim("/Volumes/groups/cochella/jiwang/annotations/cellCycle_genes_worm/BioMart_worm_human_homologe.txt", sep = "\t",
                          header = TRUE)
  #homologues = homologues[which(homologues$Human.orthology.confidence..0.low..1.high.==1), ]
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  s.genes = homologues$Gene.name[match(s.genes, homologues$Human.gene.name)]
  g2m.genes = homologues$Gene.name[match(g2m.genes, homologues$Human.gene.name)]
  s.genes = s.genes[which(!is.na(s.genes))]
  g2m.genes = g2m.genes[which(!is.na(g2m.genes))]
  
  c3.genes = list(s.genes = s.genes, g2m.genes = g2m.genes)
  return(c3.genes)
}

find.cc.markers.GO = function()
{
  s.genes = read.delim("../../../../annotations/cellCycle_genes_worm/GO-0006260_worm_dnasynthesis.txt", 
                       sep = "\t", header = FALSE)
  gm1 = read.delim("../../../../annotations/cellCycle_genes_worm/GO-0006260_worm_G2Mphase.txt", 
                   sep = "\t", header = FALSE)
  gm2 = read.delim("../../../../annotations/cellCycle_genes_worm/GO-0006260_worm_Mphase.txt", 
                   sep = "\t", header = FALSE)
  gm3 = read.delim("../../../../annotations/cellCycle_genes_worm/GO-0006260_worm_G2phase.txt", 
                   sep = "\t", header = FALSE)
  s.genes = unique(s.genes[, 3])
  g2m.genes = c(unique(as.character(gm1[,3])), unique(as.character(gm2[, 3])), unique(as.character(gm3[, 3])))
  
  c3.genes = list(s.genes = s.genes, g2m.genes = g2m.genes)
  return(c3.genes)
}

find.cellcycle.markers = function(list.sel = 'homologues')
{
  ##########################################
  # different strategies to find cell cycle markers of c. elegans for Seurat
  # 1st method): using homologue between human and c. elegans
  # 2nd method): find a list of 
  ##########################################
  if(list.sel == 'homologues') c3.genes = find.cc.markers.homologues() 
  if(list.sel == "curated") c3.genes = find.cc.markers.GO()
  s.genes = c3.genes$s.genes
  g2m.genes = c3.genes$g2m.genes
  
  # manually add genes from wormbook http://www.wormbook.org/chapters/www_cellcyclereguln/cellcyclereguln.html
  #s.genes = c(as.character(s.genes), c("cye-1")
  s.genes = unique(c(as.character(s.genes), c('cye-1', 'cya-1', 'evl-18')))
  g2m.genes = unique(c(as.character(g2m.genes), c('cdk-1', 'mat-1', 'mat-2', 'mat-3', 'emb-27', 'emb-30', 'mdf-1', 'san-1')))
  
  c3.genes = list(s.genes = s.genes, g2m.genes = g2m.genes)
  return(c3.genes)
  
  # test the code from https://github.com/hbc/macrae_ghazizadeh_zebrafish_heart_sc/blob/master/seurat_cluster/seurat_cluster_adapted_WT.Rmd
  # unfornately it does not work, because several pacakges can not be properly installed
  Test.query.cellCycle.markers = FALSE
  if(Test.query.cellCycle.markers){
    require(plotly)
    require(remotes)
    annot <- basejump::annotable("Danio rerio") %>% 
      dplyr::select(c(ensgene,symbol)) %>% 
      dplyr::mutate(symbol = toupper(symbol)) 
    cell_cycle_markers <- bcbioSingleCell::cellCycleMarkers[[camel("mus musculus")]] %>% 
      dplyr::mutate(symbol = toupper(symbol)) %>% dplyr::inner_join(annot,by = "symbol") %>% 
      dplyr::select(-c(ensgene.x)) %>% dplyr::rename(ensgene = ensgene.y)
    stopifnot(is.data.frame(cell_cycle_markers))
    markdownHeader("S phase markers", level = 3)
    s_genes <- cell_cycle_markers %>%
      filter(phase == "S") %>%
      pull("ensgene")
    print(s_genes)
    markdownHeader("G2/M phase markers", level = 3)
    g2m_genes <- cell_cycle_markers %>%
      filter(phase == "G2/M") %>%
      pull("ensgene")
    print(g2m_genes)
    saveData(cell_cycle_markers, s_genes, g2m_genes, dir = data_dir)
  }
  
}

cellCycle.correction = function(ms, method = "seurat")
{
  if(method == "seurat"){
    
    pdfname = paste0(resDir, "/scRNAseq_cellCycle_regression_Seurat.pdf")
    pdf(pdfname, width=12, height = 6)
    
    #library(scater)
    #library(scran)
    # install loomR from GitHub using the remotes package remotes::install_github(repo = 'mojaveazure/loomR', ref = 'develop')
    #library(loomR)
    #library(Seurat)
    # convert sce to seurat object (see https://satijalab.org/seurat/v3.0/conversion_vignette.html)
    #seurat = as.Seurat(sce, counts = "counts", data = "logcounts")
    #Idents(seurat) <- colnames(seurat) # quite important this assignment for cell identity
    #seurat <- FindVariableFeatures(seurat, selection.method = "vst")
    
    detach("package:scater", unload=TRUE)
    detach("package:scran", unload=TRUE)
    
    seurat = ms;
    
    # Identify the 10 most highly variable genes
    top10 <- head(VariableFeatures(seurat), 25)
    plot1 <- VariableFeaturePlot(seurat)
    plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
    CombinePlots(plots = list(plot1, plot2))
    
    #seurat <- ScaleData(seurat, features = rownames(seurat), model.use = "linear") # standardize the data (x - mean(x))/sd(x)
    #seurat <- RunPCA(seurat, features = VariableFeatures(seurat), ndims.print = 6:10, nfeatures.print = 10)
    #DimPlot(seurat, reduction = "pca")
    # DimHeatmap(seurat, dims = c(1, 2))
    
    #source("scRNAseq_functions.R")
    c3.genes = find.cellcycle.markers(list.sel = "homologues")
    
    s.genes <- c3.genes$s.genes
    g2m.genes <- c3.genes$g2m.genes 
    s.genes = s.genes[which(!is.na(match(s.genes, rownames(seurat))))]
    g2m.genes = g2m.genes[which(!is.na(match(g2m.genes, rownames(seurat))))]
    
    seurat <- CellCycleScoring(seurat, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    
    # view cell cycle scores and phase assignments
    head(seurat[[]])
    
    p00 = RidgePlot(seurat, features = c("cdk-1", "cdk-4", "cyd-1", "cye-1", "cya-1", "wee-1.3"), ncol = 2)
    plot(p00)
    
    #seurat <- RunPCA(seurat, features = VariableFeatures(seurat), verbose = TRUE)
    #DimPlot(seurat, reduction = 'pca')
    
    seurat <- RunPCA(seurat, features = c(as.character(s.genes), as.character(g2m.genes)))
    DimPlot(seurat, reduction = 'pca')
    
    # regress out the cell cycle
    seurat1 <- ScaleData(seurat, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(seurat))
    
    seurat1 <- RunPCA(seurat1, features = VariableFeatures(seurat1), nfeatures.print = 10)
    DimPlot(seurat1)
    
    # When running a PCA on only cell cycle genes, cells no longer separate by cell-cycle phase
    seurat1 <- RunPCA(seurat1, features = c(s.genes, g2m.genes))
    p1 = DimPlot(seurat1)
    plot(p1)
    # regressing out the difference between the G2M and S phase scores
    seurat$CC.Difference <- seurat$S.Score - seurat$G2M.Score
    seurat2 <- ScaleData(seurat, vars.to.regress = "CC.Difference", features = rownames(seurat))
    
    # cell cycle effects strongly mitigated in PCA
    seurat2 <- RunPCA(seurat2, features = VariableFeatures(seurat2), nfeatures.print = 10)
    DimPlot(seurat2)
    
    # when running a PCA on cell cycle genes, actively proliferating cells remain distinct from G1
    # cells however, within actively proliferating cells, G2M and S phase cells group together
    seurat2 <- RunPCA(seurat2, features = c(s.genes, g2m.genes))
    p2 = DimPlot(seurat2)
    plot(p2)
    
    # save cell cycle scoring and corrected matrix
    library(scater)
    sce$S.Score = seurat$S.Score
    sce$G2M.Score = seurat$G2M.Score
    sce$Phase = seurat$Phase
    #sce$Phase.GO = seurat$old.ident
    sce$CC.Difference = seurat$CC.Difference
    
    #xx = as.data.frame(seurat@assays$RNA@scale.data); rownames(xx) = rownames(sce)
    #assay(sce, "logcounts_seurat") <- xx
    xx = as.data.frame(seurat1@assays$RNA@scale.data); rownames(xx) = rownames(sce)
    assay(sce, "logcounts_seurat_cellcycleCorrected") <- xx
    xx = as.data.frame(seurat2@assays$RNA@scale.data); rownames(xx) = rownames(sce)
    assay(sce, "logcounts_seurat_SG2MCorrected") <- xx
    
    save(sce, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2.Rdata'))
    
    dev.off()
    
    # return(sce)
    
    ##########################################
    # a peek of source code for cellCycle.scoring function from Seurat
    # https://github.com/ChristophH/in-lineage/blob/master/R/lib.R
    ##########################################
    Example.for.Seurat = FALSE
    if(Example.for.Seurat){
      library('Matrix')
      library('parallel')
      library('MASS')
      library('diffusionMap')
      library('FNN')
      library('igraph')
      library('princurve')
      library('ggplot2')
      library('inline')
      library('gplots')
      
      # for cell cycle score
      get.bg.lists <- function(goi, N, expr.bin) {
        res <- list()
        goi.bin.tab <- table(expr.bin[goi])
        for (i in 1:N) {
          res[[i]] <- unlist(lapply(names(goi.bin.tab), function(b) {
            sel <- which(expr.bin == as.numeric(b) & !(names(expr.bin) %in% goi))
            sample(names(expr.bin)[sel], goi.bin.tab[b])
          }))
        }
        return(res)
      }
      
      enr.score <- function(expr, goi, bg.lst) {
        goi.mean <- apply(expr[goi, ], 2, mean)
        bg.mean <- sapply(1:length(bg.lst), function(i) apply(expr[bg.lst[[i]], ], 2, mean))
        return((goi.mean - apply(bg.mean, 1, mean)) / apply(bg.mean, 1, sd))
      }
      
      get.cc.score <- function(cm, N=100, seed=42, 
                               s.gene.file='./annotation/s_genes.txt',
                               g2m.gene.file='./annotation/g2m_genes.txt')
      {
        set.seed(seed)
        cat('get.cc.score, ')
        cat('number of random background gene sets set to', N, '\n')
        
        min.cells <- 5
        
        cells.mols <- apply(cm, 2, sum)
        gene.cells <- apply(cm>0, 1, sum)
        cm <- cm[gene.cells >= min.cells, ]
        
        gene.mean <- apply(cm, 1, mean)
        
        breaks <- unique(quantile(log10(gene.mean), probs = seq(0,1, length.out = 50)))
        gene.bin <- cut(log10(gene.mean), breaks = breaks, labels = FALSE)
        names(gene.bin) <- rownames(cm)
        gene.bin[is.na(gene.bin)] <- 0
        
        regev.s.genes <- read.table(file=s.gene.file, header=FALSE, stringsAsFactors=FALSE)$V1
        regev.g2m.genes <- read.table(file=g2m.gene.file, header=FALSE, stringsAsFactors=FALSE)$V1
        
        goi.lst <- list('S'=rownames(cm)[!is.na(match(toupper(rownames(cm)), regev.s.genes))],
                        'G2M'=rownames(cm)[!is.na(match(toupper(rownames(cm)), regev.g2m.genes))])
        
        n <- min(40, min(sapply(goi.lst, length)))
        goi.lst <- lapply(goi.lst, function(x) x[order(gene.mean[x], decreasing = TRUE)[1:n]])
        
        bg.lst <- list('S'=get.bg.lists(goi.lst[['S']], N, gene.bin),
                       'G2M'=get.bg.lists(goi.lst[['G2M']], N, gene.bin))
        
        all.genes <- sort(unique(c(unlist(goi.lst, use.names=FALSE), unlist(bg.lst, use.names=FALSE))))
        
        expr <- log10(cm[all.genes, ]+1)
        
        s.score <- enr.score(expr, goi.lst[['S']], bg.lst[['S']])
        g2m.score <- enr.score(expr, goi.lst[['G2M']], bg.lst[['G2M']])
        
        phase <- as.numeric(g2m.score > 2 & s.score <= 2)
        phase[g2m.score <= 2 & s.score > 2] <- -1
        
        return(data.frame(score=s.score-g2m.score, s.score, g2m.score, phase))
      }
    }
    
  }
  
  ##########################################
  # scran method is based on trained classifier for mouse or human
  # so at the end it not usable
  ##########################################
  if(method == 'scran'){
    set.seed(100)
    library(scran)
    mm.pairs <- readRDS(system.file("exdata", "mouse_cycle_markers.rds", 
                                    package="scran"))
    assignments <- cyclone(sce, mm.pairs, gene.names=rowData(sce)$ENSEMBL)
    plot(assignments$score$G1, assignments$score$G2M, 
         xlab="G1 score", ylab="G2/M score", pch=16)
    sce$phases <- assignments$phases
    table(sce$phases)
  }
  
  if(method == "scLVM"){
    ##########################################
    # select the python verson to use for Rstudio
    # https://cran.r-project.org/web/packages/reticulate/vignettes/versions.html
    # still does not work at the end
    # we change the strategy: prepare the tables and run scLVM in the conda version
    ##########################################
    # system("python --version")
    # system("which python")
    # 
    # library(reticulate)
    # use_python("/Users/jiwang/anaconda3/envs/scLVM/bin/python")
    # #use_condaenv(condaenv = "scLVM", conda = "/Users/jiwang/anaconda3/condabin/conda", required = TRUE)
    # py_config()
    # 
    # system("python --version")
    # system("which python")
    # 
    # Sys.setenv(PATH = paste("/Users/jiwang/anaconda3/envs/scLVM/bin", Sys.getenv("PATH"),sep=":"))
    
    #install.packages("rPython", type = "source")
    #install.packages("/Users/jiwang/src_packages/scLVM_0.99.3.tar.gz", repos = NULL, type="source")
    
    ##########################################
    # example code from scLVM R tutorial
    # https://github.com/PMBio/scLVM/blob/master/R/tutorials/scLVM_vignette.Rmd
    ##########################################
    library(rPython)
    library(genefilter)
    library(statmod)
    require(ggplot2)
    library(gplots)
    require(DESeq2)
    library(scLVM)
    
    #limix_path = '/Users/jiwang/anaconda2/envs/scLVM/bin/python'
    #configLimix(limix_path)
    
    data(data_Tcells)
    help(data_Tcells)
    
    #dataMouse[ 1:5, 1:4 ]
    geneTypes <- factor( c( ENSM="ENSM", ERCC="ERCC" )[
      substr( rownames(dataMouse), 1, 4 ) ] )
    #2. calculate normalisation for counts
    countsMmus <- dataMouse[ which( geneTypes=="ENSM" ), ]
    countsERCC <- dataMouse[ which( geneTypes=="ERCC" ), ]
    lengthsMmus <- dataMouse[ which( geneTypes=="ENSM" ), 1 ]
    lengthsERCC <- dataMouse[ which( geneTypes=="ERCC" ), 1 ]
    sfERCC <- estimateSizeFactorsForMatrix( countsERCC )
    sfMmus <- sfERCC #also use ERCC size factor for endogenous genes
    #normalise read counts
    nCountsERCC <- t( t(countsERCC) / sfERCC )
    nCountsMmus <- t( t(countsMmus) / sfERCC )
    
    countsMmus = counts(sce)
    sfERCC = estimateSizeFactorsForMatrix(countsMmus)
    sfMmus <- sfERCC
    nCountsMmus = t( t(countsMmus) / sfERCC )
    #use spike in to find tehcnical noise. 
    # If no spike-ins are available, we can also use the endogenous read counts for fitting the mean-CV2 relation using a log-linear fit in the log-space.
    # Alternatively, we can fit the mean-variance relationship in the log-space using local 2nd order polynomial regression (loess).
    #techNoise = fitTechnicalNoise(nCountsMmus,nCountsERCC=nCountsERCC, fit_type = 'counts')  
    techNoiseLogFit = fitTechnicalNoise(nCountsMmus, fit_type = 'log', use_ERCC = FALSE, plot=TRUE) 
    #techNoiseLogVarFit = fitTechnicalNoise(nCountsMmus, fit_type = 'logvar', use_ERCC = FALSE, plot=TRUE) 
    
    #call variable genes
    #is_het = getVariableGenes(nCountsMmus, techNoiseLogFit$fit, method = "fit", 
    #                          threshold = 0.1, fit_type="log",sfEndo=sfMmus, sfERCC=sfERCC)
    #table(is_het)
    
    #we an also do this for the other fits
    is_hetLog = getVariableGenes(nCountsMmus, techNoiseLogFit$fit, plot=TRUE)
    table(is_hetLog)
    #is_hetLogVar = getVariableGenes(nCountsMmus, techNoiseLogVarFit$fit, plot=TRUE)
    #table(is_hetLogVar)
    
    #get cell cycle genes from GO
    cc.genes = read.delim("../../../../annotations/cellCycle_genes_worm/GO_0007049_genes.txt", 
                          sep = "\t", header = FALSE)
    cc.genes = unique(cc.genes[,3])
    
    #rename a few variables
    Y = t(log10(nCountsMmus+1)) #normalised trandformed read counts
    genes_het_bool = as.vector(is_hetLog) #variable genes
    #genes_het_bool[]
    geneID = rownames(nCountsMmus) #gene IDs
    tech_noise = as.vector(techNoiseLogFit$techNoiseLog) #technical noise
    ens_ids_cc <- cc.genes
    
    index.cc = match(cc.genes, geneID)
    index.cc = index.cc[which(!is.na(index.cc))]
    ##########################################
    # can not proceed anymore and save tables for python in conda
    ##########################################
    #construct and initialize new scLVM object
    sclvm = new("scLVM")
    sclvm = init(sclvm,Y=Y,tech_noise = tech_noise)
    
    # CellCycleARD = fitFactor(sclvm,geneSet = ens_ids_cc, k=20,use_ard = TRUE)
    
    write.table(Y, file = paste0(tabDir, "gene_expression_matrx_4scLVM.txt"), sep = "\t", row.names = FALSE, col.names = FALSE)
    write.table(tech_noise, file = paste0(tabDir, "tech_noise_4scLVM.txt"), sep = "\t",row.names = FALSE, col.names = TRUE, quote = FALSE )
    write.table(geneID, file =paste0(tabDir, "geneNames_4scLVM.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE )
    write.table(which(genes_het_bool==TRUE), file =paste0(tabDir, "index_hetgenes_4scLVM.txt"), sep = "\t", 
                row.names = FALSE, col.names = TRUE, quote = FALSE)
    write.table(index.cc, file =paste0(tabDir, "index_ccgenes_4scLVM.txt"), sep = "\t", 
                row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  if(method == "ccRemover"){
    # see examplel from https://cran.r-project.org/web/packages/ccRemover/vignettes/ccRemover_tutorial.html
    # this method is pretty slow and probably risky (to test), because simply PCA was done and 
    # then compare cell cycle genes were compares with control genes and significantly different PCs were removed
    require(ccRemover)
    t.cell_data = logcounts(sce)
    head(t.cell_data[,1:5])
    
    summary(apply(t.cell_data,1, mean))
    mean_gene_exp <- rowMeans(t.cell_data)
    t_cell_data_cen <- t.cell_data - mean_gene_exp
    summary(apply(t_cell_data_cen,1,mean))
    
    gene_names <- rownames(t_cell_data_cen)
    # cell_cycle_gene_indices <- gene_indexer(gene_names, species = "mouse", 
    #                                         name_type = "symbols" )
    # length(cell_cycle_gene_indices)
    # if_cc <- rep(FALSE,nrow(t_cell_data_cen)) 
    # if_cc[cell_cycle_gene_indices] <- TRUE
    # summary(if_cc)
    
    cc.genes = read.delim("../../../../annotations/cellCycle_genes_worm/GO_0007049_genes.txt", header = FALSE)
    cc.genes = unique(cc.genes$V3)
    cc.index = match(cc.genes, gene_names)
    cc.index = cc.index[which(!is.na(cc.index))]
    
    if_cc <- rep(FALSE,nrow(t_cell_data_cen))
    if_cc[cc.index] <- TRUE
    summary(if_cc)
    
    dat <- list(x=t_cell_data_cen, if_cc=if_cc)
    xhat <- ccRemover(dat, bar=TRUE, max_it = 6, nboot = 100)
    
    xhat <- xhat + mean_gene_exp
    
    pca1 = prcomp(t(t.cell_data[if_cc,]), scale. = TRUE)
    pca2 = prcomp(t(xhat[if_cc,]), scale. = TRUE)
    par(mfrow = c(1, 2))
    plot(pca1$x[, c(2:3)])
    plot(pca2$x[, c(1:2)])
    
  }
  
}



########################################################
########################################################
# Section : function for bach correction
# 
########################################################
########################################################
Check.MNN.pairs = function(mnn.out, fscs)
{
  #library(igraph)
  out = mnn.out$pairs
  set.seed(1001)
  diffs = abs(sample(fscs, size = 20000, replace = TRUE) - sample(fscs, size = 20000, replace = TRUE))
  labels = rep('random', length(diffs))
  #bc.values = as.character(mnn.out$batch)
  
  for(n in 1:length(out))
  {
    pairs = as.data.frame(out[[n]])
    diffs = c(diffs, abs(fscs[pairs[,1]] - fscs[pairs[,2]]))
    labels = c(labels, rep(paste0(n, "- merge batch #", mnn.out$order[(n+1)]), nrow(pairs))) 
    # g <- graph_from_data_frame(pairs, directed = FALSE)
    # #bipartite.mapping(g)
    # V(g)$type <- bipartite_mapping(g)$type
    # V(g)$label <- V(g)$name # set labels.
    # #plot(g)
    # plot(g, vertex.label.cex = 0.8, vertex.label.color = "black", cex = 0.1)
    # plot(g, layout=layout.bipartite, vertex.size=1, vertex.label.cex=0.6)
  }
  
  diffs.data = data.frame(diffs, labels, stringsAsFactors = TRUE)
  diffs.data$labels = as.factor(diffs.data$labels)
  p = ggplot(diffs.data, aes(x=labels, y=diffs)) + 
    geom_violin(trim=TRUE)
  
  plot(p)
  # boxplot(diffs ~ labels, data=diffs)
  
}


batchCorrection_fastMNN = function(sce, Use.fastMNN = TRUE)
{
  library(scRNA.seq.funcs)
  library(scater)
  library(SingleCellExperiment)
  library(scran)
  library(kBET)
  set.seed(1234567)
  options(stringsAsFactors = FALSE)
  
  load(file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_SCE_seuratCellCycleCorrected_v2_facsInfos_timingEstGroups.Rdata'))
  
  sce$timingEst = as.factor(sce$timingEst)
  sce$timingEst.group = as.factor(sce$timingEst.group)
  
  # choose the batches (either plates or request)
  # here we choose the request as batch
  batches = sce$request
  bc.uniq = unique(batches)
  sce$batches <- batches
  
  Use.fastMNN = TRUE
  Norm.Vars.per.batch = TRUE # HVGs for each batch or not
  Rescale.Batches = FALSE # scaling data in each batch or not
  k.mnn = 20
  cos.norm = TRUE
  nb.pcs = 50
  
  batch.sequence.to.merge = c('R7130', 'R8729', 'R8612', 'R8526', 'R7926', # 3 plates for each request
                              'R6875','R8613','R8348') # 1 plate for each request
  
  order2correct = match(batch.sequence.to.merge, bc.uniq)
  #c(12, 13, # the same
  #                10, 9, 8, 5, 6, 7,  8, 11, 1, 2, 3, 4)
  #order2correct = c(15,14,13,12,11,10,9,8, 1, 2, 3, 4, 5,6,7 )
  #order2correct = c(3, 4, 1, 2)
  
  ## double chekc  the mering order in the batch correction
  source('customizedClustering_functions.R')
  kk = match(sce$request, c('R7130', 'R8729', 'R8612', 'R8526', 'R7926'))
  kk = match(sce$request, c('R6875','R8613','R8348'))
  plotColData(sce[,which(!is.na(kk))],
              x = "FSC_log2",
              y = "BSC_log2",
              colour_by = "request",
              point_size = 1
  )
  
  #cat("merging order for batch correction :\n", paste0(bc.uniq[order2correct], collapse = "\n"), "\n")
  for(n in 1:length(order2correct)){
    
    #n = 11
    kk = order2correct[n]
    
    p = plotColData(sce[, which(sce$batches== bc.uniq[kk])],
                    x = "FSC_log2",
                    y = "BSC_log2",
                    colour_by = "timingEst",
                    point_size = 1
                    
    )
    plot(p)
    
    cat('#', n, 'batch:',  bc.uniq[kk], ': ', length(which(sce$batches == bc.uniq[kk])), 'cells \n')
    
  }
  
  
  source("scRNAseq_functions.R")
  HVGs = find.HVGs(sce, Norm.Vars.per.batch = Norm.Vars.per.batch, method = "scran", ntop = 2000)
  gene.chosen = match(HVGs, rownames(sce))
  
  cat("nb of HGV : ", length(gene.chosen), "\n")
  
  
  
  if(Use.fastMNN){
    ## rescaling for each batch is recommended by the author
    ## We adjust the size factors with multiBatchNorm() to make them more comparable across batches.
    ## This mitigates differences in scale and variance in the log-expression values between batches, especially between technologies.
    ## https://www.bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/multibatch.html#3_feature_selection_across_batches
    if(Rescale.Batches){
      ttxt = c("nout = multiBatchNorm(")
      for(n in 1:length(bc.uniq)){
        if(n != length(bc.uniq)) {
          ttxt = paste0(ttxt, "sce[, which(sce$batches == '", bc.uniq[n], "')], ")
        } else {
          ttxt = paste0(ttxt, "sce[, which(sce$batches == '", bc.uniq[n], "')])")
        }
      }
      eval(parse(text = ttxt)) ## rescaling done
      
      kk2check = 1
      par(mfrow=c(1,1))
      plot(sizeFactors(nout[[kk2check]]), sizeFactors(sce[, which(sce$batches == bc.uniq[kk2check])])); abline(0, 1, col='red')
    }
    
    original = list()
    fscs = c()
    #original0 = list()
    for(n in 1:length(bc.uniq)){
      #xx = nout[[n]];
      if(Rescale.Batches){
        original[[n]] = logcounts((nout[[n]][gene.chosen, ]))
      }else{
        original[[n]] = logcounts((sce[gene.chosen, which(sce$batches == bc.uniq[n])]))
      }
      fscs = c(fscs, sce$FSC_log2[which(sce$batches == bc.uniq[n])])
    }
    
    # Slightly convoluted call to avoid re-writing code later.
    # Equivalent to fastMNN(GSE81076, GSE85241, k=20, d=50, approximate=TRUE)
    # Comments from Aaron: https://www.bioconductor.org/packages/devel/workflows/vignettes/simpleSingleCell/inst/doc/batch.html#5_examining_the_effect_of_correction
    # The k= parameter specifies the number of nearest neighbours to consider when defining MNN pairs. This should be interpreted as the minimum frequency of each cell type or state in each batch.
    # Larger values will improve the precision of the correction by increasing the number of MNN pairs. It also provides some robustness to violations of the assumption that the batch vector is orthogonal to the biological subspace (Haghverdi et al. 2018), by allowing the neighbour search to ignore biological variation in each batch to identify the correct MNN pairs.
    # However, larger values of k can also reduce accuracy by allowing incorrect MNN pairs to form between cells of different types. Thus, we suggest starting with the default k and increasing it if one is confident that the same cell types are not adequately merged across batches. This is better than starting with a large k as incorrect merging is much harder to diagnose than insufficient merging.
    # When BSPARAM=IrlbaParam(deferred=TRUE), fastMNN() uses methods from the irlba package to perform the principal components analysis quickly. While the run-to-run differences should be minimal, it does mean that set.seed() is required to obtain fully reproducible results. The deferred= argument instructs fastMNN() to sacrifice some numerical precision for greater speed.
    set.seed(1001)
    mnn.out <- do.call(fastMNN, c(original, list(k=k.mnn, cos.norm = cos.norm, d=nb.pcs, auto.order=order2correct,
                                                 approximate=TRUE)))
    dim(mnn.out$corrected)
    mnn.out$batch
    Rle(mnn.out$batch)
    #metadata(mnn.out)$merge.info$pairs[[1]]
    reducedDim(sce, "MNN") <- mnn.out$corrected
    sce$mnn_Batch <- as.character(mnn.out$batch)
    sce
    
    ##########################################
    # check the effectiveness of batch correction with MNN
    # 1) check the MNN pairs and lost variances
    # 2) visualization wiht PCA, and UMAP before and after correction
    # 3) kBET test
    ##########################################
    pdfname = paste0(resDir, "/scRNAseq_filtered_test_MNNbatchCorrection_effectiveness_v1.pdf")
    pdf(pdfname, width=14, height = 8)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    set.seed(1001)
    nb.MNN.to.use = 20
    sce <- runUMAP(sce, use_dimred="MNN", n_dimred = nb.MNN.to.use, ncomponents = 2, scale_features = FALSE,
                   method = c("umap-learn"))
    
    p = plotUMAP(sce, ncomponents = 2, colour_by="timingEst.group", size_by = "FSC_log2", point_size= 0.01) + ggtitle("Corrected")
    plot(p)
    
    fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=12))
    p1 = plotUMAP(sce, colour_by="pha-4", size_by = "FSC_log2") +
      fontsize + ggtitle("MNN corrected")
    p2 = plotUMAP(sce, colour_by="hnd-1", size_by = "FSC_log2") +
      fontsize + ggtitle("MNN corrected")
    multiplot(p1, p2, cols = 2)
    
    sce = runTSNE(sce, use_dimred="MNN", perplexity = 30, n_dimred = nb.MNN.to.use)
    p = plotTSNE(sce, colour_by="timingEst", size_by = "FSC_log2") +
      fontsize + ggtitle("MNN corrected")
    
    p1 = plotTSNE(sce, colour_by="pha-4", size_by = "FSC_log2") +
      fontsize + ggtitle("MNN corrected")
    p2 = plotTSNE(sce, colour_by="hnd-1", size_by = "FSC_log2") +
      fontsize + ggtitle("MNN corrected")
    multiplot(p1, p2, cols = 2)
    
    # mnn.out$pairs
    source("scRNAseq_functions.R")
    Check.MNN.pairs(mnn.out, fscs)
    
    #omat <- do.call(cbind, original)
    #sce.qc <- SingleCellExperiment(list(logcounts=omat))
    set.seed(1000)
    with.var <- do.call(fastMNN, c(original,
                                   list(k=k.mnn, cos.norm = cos.norm, d=nb.pcs, auto.order=order2correct, approximate=TRUE,
                                        compute.variances=TRUE)))
    with.var$lost.var
    
    fsc.boundary = c(27.2, 28.6)
    bsc.boundary = c(24.5, 26.8)
    plotColData(sce,
                x = "FSC_log2",
                y = "BSC_log2",
                colour_by = "mnn_Batch",
                shape_by = "mnn_Batch"
                
    ) + geom_vline(xintercept = fsc.boundary, linetype="dashed", color = "blue", size=0.5) +
      geom_hline(yintercept= bsc.boundary, linetype="dashed", color = "red", size=0.5)
    
    sel.tmp = which(sce$FSC_log2 > fsc.boundary[1] & sce$FSC_log2 < fsc.boundary[2]
                    & sce$BSC_log2 > bsc.boundary[1] & sce$BSC_log2 < bsc.boundary[2])
    #sce.tmp = sce[gene.chosen, which(sce$mnn_Batch > 2)]
    sce.tmp = sce[, sel.tmp]
    #sce.tmp = sce
    dim(sce.tmp)
    
    sce.tmp <- runPCA(sce.tmp, ncomponents = 50, ntop=nrow(sce.tmp), method="irlba", exprs_values = "logcounts", scale_features = TRUE)
    
    plotPCA(sce.tmp, colour_by="batches") + ggtitle("Original")
    
    dff = as.data.frame(reducedDim(sce.tmp, "MNN"))
    colnames(dff) = paste0("PC", c(1:ncol(dff)))
    dff$batches = sce.tmp$batches
    ggp = ggplot(data=dff, aes(PC1, PC2, color=batches)) + geom_point(size=2)
    plot(ggp);
    
    # Using irlba to set up the t-SNE, for speed.
    #set.seed(100)
    #osce <- runTSNE(sce.tmp, use_dimred="PCA", perplexity = 20)
    #ot <- plotTSNE(osce, colour_by="mnn_Batch") + ggtitle("Original")
    #set.seed(100)
    #csce <- runTSNE(sce, use_dimred="MNN", perplexity = 20)
    #ct <- plotTSNE(csce, colour_by="mnn_Batch") + ggtitle("Corrected")
    #multiplot(ot, ct, cols=2)
    
    set.seed(100)
    osce <- runUMAP(sce.tmp, use_dimred="PCA", n_dimred = nb.MNN.to.use)
    ot <- plotUMAP(osce, colour_by="mnn_Batch", size_by = "FSC_log2") + ggtitle("Original")
    set.seed(100)
    csce <- runUMAP(sce.tmp, use_dimred="MNN", n_dimred = nb.MNN.to.use)
    ct <- plotUMAP(csce, colour_by="mnn_Batch", size_by = "FSC_log2") + ggtitle("Corrected")
    multiplot(ot, ct, cols=2)
    
    set.seed(100)
    osce <- runTSNE(sce.tmp, use_dimred="PCA", perplexity = 20, n_dimred = nb.MNN.to.use)
    ot <- plotTSNE(osce, colour_by="mnn_Batch", size_by = "FSC_log2") + ggtitle("Original")
    set.seed(100)
    csce <- runTSNE(sce.tmp, use_dimred="MNN", perplexity = 20, n_dimred = nb.MNN.to.use)
    ct <- plotTSNE(csce, colour_by="mnn_Batch", size_by = "FSC_log2") + ggtitle("Corrected")
    multiplot(ot, ct, cols=2)
    
    kbet.orig <- kBET(
      df = t(reducedDim(sce.tmp, "PCA")),
      batch = sce.tmp$mnn_Batch,
      heuristic = TRUE,
      do.pca = FALSE,
      verbose = TRUE,
      addTest = FALSE,
      n_repeat = 200,
      plot = TRUE)
    
    kbet.bc = kBET(
      df = t(reducedDim(sce.tmp, "MNN")),
      batch = sce.tmp$mnn_Batch,
      do.pca = FALSE,
      heuristic = FALSE,
      verbose = TRUE,
      addTest = FALSE,
      n_repeat = 200,
      plot = TRUE)
    
    dev.off()
    
  }
  
  
}

batchCorrection_mnnCorrect = function()
{
  if(Use.mnnCorrect){
    set.seed(100)
    require(BiocParallel)
    require(stats)
    bpp <- MulticoreParam(5)
    bpp
    
    mnn.out2 = mnnCorrect(R6879, R7116, R7130, k = 20, sigma = 0.1, cos.norm.in = TRUE, cos.norm.out = TRUE, order = c(3, 2, 1), 
                          svd.dim = 0, var.adj = FALSE, subset.row = gene.chosen, pc.approx = TRUE, BPPARAM=bpp)
    
    head(mnn.out2$pairs[[2]])
    
    ## check the pairs of cells used between batches
    Check.SNN.pairs = FALSE
    if(Check.SNN.pairs){
      for(n in 2:length(mnn.out$pairs))
      {
        #n = 3
        paires = data.frame(mnn.out2$pairs[[n]])
        
        cat("cell in the batch", n, ":  ",
            "total pairs --", nrow(paires), ",",  
            "paired cell in batch", n-1, "-- ", length(unique(paires[which(paires[,3]==(n-1)), 2])), ",", 
            "paired cell in batch", n+1, "-- ", length(unique(paires[which(paires[,3]==(n+1)), 2])), "\n"
        )
      }
    }
    
    res1 <- mnn.out2$corrected[[1]]
    res2 <- mnn.out2$corrected[[2]]
    res3 <- mnn.out2$corrected[[3]]
    
    dimnames(res1) <- dimnames(R6879)
    dimnames(res2) <- dimnames(R7116)
    dimnames(res3) <- dimnames(R7130)
    res = cbind(res1, res2, res3)
    
    assay(sce, "corrected") <- res
    
    #sce = sce.qc
    #save(sce, mnn.out, mnn.out2, file=paste0(RdataDir, version.DATA, '_QCed_cells_genes_filtered_normalized_batchCorrectMNN_SCE.Rdata'))
    
    ##########################################
    # double check the relatinship between fastMNN and mnnCorrect 
    ##########################################
    #osce = runPCA(sce, ncomponents = 50, ntop=Inf, method="irlba", exprs_values = "logcounts", feature_set = gene.chosen)
    csce <- runPCA(sce, ncomponents = 50, method="irlba", exprs_values = "corrected", feature_set = gene.chosen,
                   scale_features = TRUE, detect_outliers = FALSE)
    
    xx = as.data.frame(reducedDim(csce, "MNN"));
    yy = as.data.frame(reducedDim(csce, "PCA"))
    
    head(xx[, c(1:10)])
    head(yy[, c(1:10)])
    
    par(mfrow = c(1, 2))
    plot(xx[, c(1:2)], xlab = 'PC1', ylab = "PC2", main = "lowDim output from fastMNN")
    plot(yy[, c(1:2)], xlab = 'PC1', ylab = "PC2", main = "PCA from mnnCorret output")
    #plot(xx[, c(3:4)], xlab = 'PC1', ylab = "PC2", main = "lowDim output from fastMNN")
    #plot(yy[, c(3:4)], xlab = 'PC1', ylab = "PC2", main = "PCA from mnnCorret output")
    
    
    pdfname = paste0(resDir, "/scRNAseq_filtered_test_batchCorrection_fastMNNlowDim_vs_mnnCorrectOutput.pdf")
    pdf(pdfname, width=18, height = 8)
    par(cex =0.7, mar = c(3,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    
    od = scater::plotPCA(
      osce,
      run_args = list(exprs_values = "logcounts"), 
      size_by = "total_counts",
      #size_by = "total_features_by_counts",
      colour_by = "Batch"
    ) + ggtitle(paste0("PCA -- origine "))
    
    cd = scater::plotPCA(
      csce,
      run_args = list(exprs_values = "corrected"),
      size_by = "total_counts",
      #size_by = "total_features_by_counts",
      colour_by = "request"
    ) + ggtitle(paste0("PCA -- corrected by mnn"))
    
    multiplot(od, cd, cols=2)
    
    dev.off()
    
  }
}

batchCorrection_Scanorama = function()
{
  ## quick test for one batch correction method Scanorama
  ## example code from https://github.com/brianhie/scanorama
  #Sys.setenv(PATH = paste0(c("/Users/jiwang/anaconda3/envs/Scanorama/bin", Sys.getenv("PATH")), collapse = ":"))
  #Sys.setenv(RETICULATE_PYTHON = "/Users/jiwang/anaconda3/envs/Scanorama/bin/python")
  reticulate::py_config()
  #install.packages("reticulate")
  
  library(reticulate)
  #Sys.which("python")
  #reticulate::py_config()
  #reticulate::use_python("/Users/jiwang/anaconda3/envs/Scanorama/bin/python")
  reticulate::py_config()
  
  #use_virtualenv("/Users/jiwang/anaconda3/envs/Scanorama")
  #use_condaenv("Scanorama")
  #reticulate::py_config()
  
  #library(reticulate)
  #xx = import('numpy')
  #Sys.which(xx)
  scanorama <- import('scanorama')
  
  #library(reticulate)
  #scanorama <- import('scanorama')
  
  # List of data sets (matrices of cells-by-genes)
  # List of gene lists:
  datasets = list()
  genes_list = list()
  for(n in 1:length(original)){
    datasets[[n]] <- t(original[[n]])
    genes_list[[n]] <- as.list(rownames(sce)[gene.chosen])
  }
  
  # Integration.
  #integrated.data <- scanorama$integrate(datasets, genes_list)
  
  # Batch correction
  ## Note that reticulate has trouble returning sparse matrices, so you should set the return_dense flag to TRUE 
  ## (which returns the corrected data as R matrix objects) when attempting to use Scanorama's correct() method in R.
  ## This will increase memory usage, however, especially for very large data sets.
  corrected.data <- scanorama$correct(datasets, genes_list, return_dense=TRUE)
  
  # Integration and batch correction.
  #integrated.corrected.data <- scanorama$correct(datasets, genes_list, return_dimred=TRUE, return_dense=TRUE)
}

## test function for Hamony
test.harmony.function = function(ms)
{
  #library(harmony)
  #pbmcsca <- RunHarmony(ms, group.by.vars = "request", assay.use="SCT")
  #pbmcsca <- RunUMAP(pbmcsca, reduction = "harmony", dims = 1:20)
  #pbmcsca <- FindNeighbors(pbmcsca, reduction = "harmony", dims = 1:30) %>% FindClusters()
  #DimPlot(pbmcsca, reduction = 'harmony',  group.by = c("request"))
  #DimPlot(pbmcsca, reduction = 'umap',  group.by = c("request"))
  
}

########################################################
########################################################
# Section : test batch correction in Seurat
# 1) MNN calling from Seurat 
# 2) CCA from Seurat
########################################################
########################################################
compare.CCA.fastMNN.inSeurat = function()
{
  library(scater)
  library(SingleCellExperiment)
  library(scran)
  
  sce <- as.SingleCellExperiment(ms)
  batches = sce$request 
  bc.uniq = unique(batches)
  sce$batches <- batches
  
  Norm.Vars.per.batch = TRUE # HVGs for each batch or not 
  Rescale.Batches = FALSE # scaling data in each batch or not 
  k.mnn = 20
  cos.norm = TRUE
  nb.pcs = 50
  
  batch.sequence.to.merge = c('R7130', 'R8612', 'R8526', 'R7926', # 3 plates for each request
                              'R6875','R7116','R8613','R8348') # 1 plate for each request
  order2correct = match(batch.sequence.to.merge, bc.uniq) 
  
  source("scRNAseq_functions.R")
  #HVGs = find.HVGs(sce, Norm.Vars.per.batch = Norm.Vars.per.batch, method = "scran", ntop = 2000)
  HVGs = VariableFeatures(ms)
  gene.chosen = match(HVGs, rownames(sce))
  
  
  cat("nb of HGV : ", length(gene.chosen), "\n")
  
  original = list()
  fscs = c()
  #original0 = list()
  for(n in 1:length(bc.uniq)){
    #xx = nout[[n]];
    if(Rescale.Batches){
      original[[n]] = logcounts((nout[[n]][gene.chosen, ]))
    }else{
      original[[n]] = logcounts((sce[gene.chosen, which(sce$batches == bc.uniq[n])])) 
    }
    fscs = c(fscs, sce$FSC_log2[which(sce$batches == bc.uniq[n])])
  }
  
  set.seed(1001)
  mnn.out <- do.call(fastMNN, c(original, list(k=k.mnn, cos.norm = cos.norm, d=nb.pcs, auto.order=order2correct,
                                               approximate=TRUE)))
  
  reducedDim(sce, "sct_MNN") <- mnn.out$corrected
  sce$mnn_Batch <- as.character(mnn.out$batch)
  sce
  
  pbmc = as.Seurat(sce)
  
  
  ms <- RunUMAP(object = pbmc, reduction = 'PCA', dims = 1:20, n.neighbors = 30)
  DimPlot(ms, reduction = "umap", group.by = 'request')
  
  ms <- RunUMAP(object = pbmc, reduction = 'MNN', dims = 1:20, n.neighbors = 30)
  DimPlot(ms, reduction = "umap", group.by = 'request')
  
  ms <- RunUMAP(object = pbmc, reduction = 'sct_MNN', dims = 1:20, n.neighbors = 30)
  DimPlot(ms, reduction = "umap", group.by = 'request')
  
  
  pbmc.list <- SplitObject(ms, split.by = "request")
  pbmcsca <- RunFastMNN(object.list = pbmc.list, features = 2000, reduction.name = 'mnn_sct', 
                        k=20, cos.norm=TRUE, ndist=3, d=50, approximate=FALSE, auto.order = order2correct)
  
  pbmcsca <- RunUMAP(pbmcsca, reduction = "mnn", dims = 1:30)
  pbmcsca <- FindNeighbors(pbmcsca, reduction = "mnn", dims = 1:30)
  pbmcsca <- FindClusters(pbmcsca)
  DimPlot(pbmcsca, group.by = c("Method", "ident", "CellType"), ncol = 3)
  
  
  pbmc.list <- SplitObject(ms, split.by = "request")
  for (i in names(pbmc.list)) {
    pbmc.list[[i]] <- SCTransform(pbmc.list[[i]], verbose = FALSE)
  }
  
  pbmc.features <- SelectIntegrationFeatures(object.list = pbmc.list, nfeatures = 3000)
  pbmc.list <- PrepSCTIntegration(object.list = pbmc.list, anchor.features = pbmc.features)
  k.filter <- min(sapply(pbmc.list, ncol)) 
  pbmc.anchors <- FindIntegrationAnchors(object.list = pbmc.list, normalization.method = "SCT", 
                                         anchor.features = pbmc.features, k.filter = k.filter)
  
  pbmc.integrated <- IntegrateData(anchorset = pbmc.anchors, normalization.method = "SCT")
  
  pbmc.integrated <- RunPCA(object = pbmc.integrated, verbose = FALSE)
  pbmc.integrated <- RunUMAP(object = pbmc.integrated, dims = 1:30)
  
  DimPlot(pbmc.integrated, reduction = "umap", group.by = 'request')
  
  plots <- DimPlot(pbmc.integrated, group.by = c("Method", "CellType"), combine = FALSE)
  plots <- lapply(X = plots, FUN = function(x) x + theme(legend.position = "top") + guides(color = guide_legend(nrow = 4, 
                                                                                                                byrow = TRUE, override.aes = list(size = 2.5))))
  CombinePlots(plots)
  
}







