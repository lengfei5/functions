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






