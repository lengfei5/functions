#### Functions
process.countTable = function(all, design, select.counts = NULL)
{
  newall = data.frame(as.character(all[, 1]), stringsAsFactors = FALSE)
  
  for(n in 1:nrow(design))
  {
    #n = 1;
    ## found the matching ID in design matrix
    if(!is.null(select.counts)){
      jj = intersect(grep(design$SampleID[n], colnames(all)), grep(select.counts, colnames(all)));
    }else{
      jj = grep(design$SampleID[n], colnames(all));
    }
    
    ## deal with several mapping or zero mapping
    if(length(jj) == 0){
      cat("no count column found for ", design$SampleID[n], "\n")
    }else{
      if(length(jj)==1) {
        #index = c(index,jj)
        cat("*** 1 count column found for ", design$SampleID[n], ":", colnames(all)[jj], "\n")
        newall = data.frame(newall, all[, jj])
      }else{
        cat(length(jj), " samples found for ID", design$SampleID[n], ": ",  colnames(all)[jj], " -- ")
        cat("Merge those columns as technical replicates...\n")
        newall = data.frame(newall, apply(as.matrix(all[, jj]), 1, sum))
      }
    }
  }
  
  colnames(newall)[1] = "gene";
  colnames(newall)[-1] = apply(design, 1, function(x) paste0(x, collapse = "_"))
  
  names = colnames(newall)
  names = sapply(names, function(x) gsub('-', '.', x))
  names = sapply(names, function(x) gsub(' ', '', x))
  colnames(newall) = names
  # if(time.series){
  #   colnames(newall)[-1] = paste0(design$stage, "_", design$treatment, "_", design$SampleID)
  # }else{
  #   colnames(newall)[-1] = paste0(design$genotype, "_", design$tissue.cell, "_", design$treatment, "_",  design$SampleID)
  # }
  
  return(newall)
}

process.countTable.v1 = function(all, design) ## old version of process.countTable function
{
  index = c()
  for(n in 1:nrow(design))
  {
    #n = 1;
    jj = intersect(grep(design$SampleID[n], colnames(all)), grep("Total.count", colnames(all)));
    if(length(jj)==1) {
      index = c(index,jj)
    }else{print(paste0("ERROR for sample--", design$SampleID[n]))}
  }
  
  newall = data.frame(as.character(all[,1]),  as.matrix(all[, index]), stringsAsFactors = FALSE)
  colnames(newall)[1] = "gene";
  colnames(newall)[-1] = paste0(design$genotype, "_", design$tissue.cell, "_", design$SampleID)
  
  return(newall)
}

cat.countTable = function(xlist)
{
  ## input is list.files for count tables (including directories and file names)
  counts = NULL
  for(n in 1:length(xlist)){
    # n = 1
    ff = read.delim(xlist[n], sep='\t', header = TRUE, as.is = c(1));
    if(n==1){
      ggs = unique(ff[, 1]);
      counts = data.frame(ggs, ff[match(ggs, ff[, 1]) , -1], stringsAsFactors = FALSE);
    }else{
      ggs = unique(c(counts[, 1],ff[, 1]));
      counts = data.frame(ggs, counts[match(ggs, counts[, 1]), -1], ff[match(ggs, ff[, 1]) , -1], stringsAsFactors = FALSE);
    }
  };
  
  colnames(counts)[1] = 'gene'
  return(counts)
  
}

compare.readCounts.umiFr.umiNum =function(design, aa, spikes){
  
  for(n in 1:nrow(design)){
    id = design$SampleID[n]
    par(mfrow=c(1,2))
    plot(aa[, intersect(grep(id, colnames(aa)), grep("Total.count", colnames(aa)))], 
         aa[, intersect(grep(id, colnames(aa)), grep("Total.UMIfr.count", colnames(aa)))], 
         log='xy', main= paste0(design[n, ], collapse = "_"), xlab = 'read counts', ylab =' umi.frations',
         cex = 0.4
         )
    points(spikes[, c(grep(paste0("Total.spikeIn.", id), colnames(spikes)), 
                      grep(paste0("Total.UMI.spikeIn.", id), colnames(spikes)))], col = 'darkblue', cex = 1., pch=16)
    abline(0, 1, lwd=1.2, col = 'red')
    
    plot(aa[, intersect(grep(id, colnames(aa)), grep("Total.count", colnames(aa)))], 
         aa[, intersect(grep(id, colnames(aa)), grep("Total.UMInum.count", colnames(aa)))], 
         log='xy', main= paste0(design[n, ], collapse = "_"), xlab = 'read counts', ylab =' umi.counts',
         cex = 0.4
    )
    points(spikes[, c(grep(paste0("Total.spikeIn.", id), colnames(spikes)), 
                      grep(paste0("Total.UMI.spikeIn.", id), colnames(spikes)))], col = 'darkblue', cex = 1., pch=16)
    abline(0, 1, lwd=1.2, col = 'red')
    
  }
  

}


find.mirName = function(x)
{
  test = unlist(strsplit(as.character(x), '-'));
  return(paste0(test[-c(1,length(test))], collapse = '-'))
} 

find.mirName.dm = function(x)
{
  test = unlist(strsplit(as.character(x), '-'));
  return(paste0(test[-c(which(test=="RM"),length(test))], collapse = '-'))
} 



Compare.three.Controls.Ovary.dm = function(read.count, design.matrix)
{
  #read.count=read.count[, kk]; design.matrix = design.matrix[, index.qc]
  require(lattice);
  require(ggplot2)
  require('DESeq2');
  library("vsn");
  library("pheatmap");
  library("RColorBrewer");
  library("dplyr"); 
  library("ggplot2")
  #load(file=paste0('Rdata/Screen_countData_sgRNA_Gene_clean_mapping', version.data, '.Rdata'))
  # kk = grep('Ab', colnames(bb))
  if(ncol(design.matrix)>2){cc = apply(design.matrix[, -1], 1, paste0, collapse="_")
  }else{cc = design.matrix[, -1]}
  #o1 = order(cc)
  #read.count = read.count[o1,]
  #cc = cc[o1]
  raw = as.matrix(read.count)
  #xx = raw
  dim(raw)
  raw[which(is.na(raw))] = 0
  xx = raw;
  
  countData = ceiling(raw)
  conds = factor(paste0(colnames(design.matrix)[-1], collapse = " + "))
  eval(parse(text = paste0("dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ ", conds, ")")))
  #dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ condition + time)
  #dds <- dds[ rowSums(counts(dds)) > 10, ]
  dds <- estimateSizeFactors(dds)
  fpm = fpm(dds, robust = TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  pca=plotPCA(vsd, intgroup = colnames(design.matrix)[-1], returnData = FALSE)
  print(pca)
  
  ## filter lowly expressed miRNAs
  sels = find.expressed.mature.miRNA(dds)
  rownames(dds) = sels$miRNA;
  dds = dds[sels$expressed, ];
  dds <- estimateSizeFactors(dds)
  #fpm = fpm(dds, robust = TRUE)
  vsd <- varianceStabilizingTransformation(dds, blind = FALSE)
  pca=plotPCA(vsd, intgroup = colnames(design.matrix)[-1], returnData = FALSE)
  print(pca)
  cpm = fpm(dds, robust=TRUE)
  cpm = log2(cpm+0.25)
  
  pairs(cpm[, grep('_untreated', colnames(cpm))], lower.panel=NULL, upper.panel=panel.fitting)
  pairs(cpm[, grep('_treated', colnames(cpm))], lower.panel=NULL, upper.panel=panel.fitting)
  
}


mixture.gaussian.fitting = function(x, method='Normal', k=2, mu.init=NULL, sigma.init=NULL, lambda.init=c(0.3, 0.7), m.constr=NULL, sigma.constr=NULL)
{
    print("to finish")
}

Filtering.Expressed.Genes = function(countData, conds, Use.only.notreated=FALSE, posterior.cutoff=0.5, fraction.detected.samples=0.5)
{
  require('edgeR')
  library(mixtools)
  #cat(design.matrix)
  y0 = countData;
  group = conds
  #group <- apply(design.matrix[, -1], 1, paste0, collapse = "_")
  if(Use.only.notreated){
    jj = which(design.matrix$condition=="notreated");
    y0 = y0[,jj]; group = group[jj];
  }
  y <- DGEList(counts=y0, group=group)
  #y1 <- calcNormFactors(y, method=c('TMM'))
  #y2 <- calcNormFactors(y, method=c('RLE'))
  #y3 <- calcNormFactors(y, method=c('none'))
  y = calcNormFactors(y, method=c("upperquartile"), p=0.5)
  cpm = cpm(y, normalized.lib.sizes=TRUE, log=TRUE, prior.count=0.25)
  
  expressed = matrix(NA, nrow=nrow(cpm), ncol= ncol(cpm));
  par(mfrow=c(2,2))
  #unique.group= unique(group)
  for(n in 1:ncol(cpm)){
    #cc = group[n]
    cat(group[n], '\n');
    data = as.numeric(unlist(cpm[, n]));
    #jj2 = index.outliers(data);
    jj1 = which(y0[, n]>=1);
    jj2 = setdiff(c(1:nrow(expressed)), jj1);
    #samples = design.matrix[which(group==cc), 1]
    #index = c()
    #for(s in samples){index = c(index, which(colnames(cpm)==s));}
    #index = unique(index);
    #if(length(index)>1){data = as.numeric(unlist(apply(cpm[, index], 1, mean)))
    #}else{data = as.numeric(unlist(cpm[, index]));}
    fit = normalmixEM(data[jj1], lambda = c(0.3, 0.7), mu=c(0, 10), sigma = NULL, mean.constr=NULL, k=2, maxrestarts=20, maxit = 1500)
    plot.mixEM(fit, whichplots = 2, main2=paste0('fit distribution of expression for ', colnames(cpm)[n]))
    mus = fit$mu;
    expressed[jj1, n] = fit$posterior[, which(mus==max(mus))]>posterior.cutoff;
    expressed[jj2, n] = FALSE
  }
  
  par(mfrow=c(1,1))
  EEs = apply(expressed, 1, sum)
  return(EEs>=(ncol(expressed)*fraction.detected.samples))
}

index.outliers = function(data.xx)
{
  c = 1.5
  #data.xx = c(2, 3, 6, 9, 13, 18, 21, 106)
  Q1 = quantile(data.xx, 0.25, type=5)
  Q3 = quantile(data.xx, 0.75, type=5)
  IQD = Q3 - Q1
  lower = Q1 - c*IQD
  upper = Q3 + c*IQD
  index = which(data.xx<lower|data.xx>upper)
  #boxplot(data.xx);abline(h=Q1);abline(h=Q3);
}

Check.3p.5p.arm.switch = function(all)
{
  xx = all[order(all$gene), ]
  ggs = sapply(xx$gene, find.mirName);
  ggs = data.frame(ggs, xx$gene, stringsAsFactors = FALSE)
  yy = as.matrix(xx[, -1]);
  yy[which(is.na(yy))] = 0
  
  Compare.5p.3p = function(x){
    x = as.numeric(x);
    if(x[1]>x[2]) return(1)
    if(x[1]==x[2]) return(0)
    if(x[1]<x[2]) return(-1)
  }
  ggs.uniq = unique(ggs[, 1])
  counts.ps = c()
  names = c()
  for(n in 1:length(ggs.uniq)){
    #n = 1
    jj = which(ggs[, 1]==ggs.uniq[n])
    if(length(jj)>1){
      names = c(names, ggs.uniq[n])
      test = yy[jj, ]
      #means = apply(test, 1, mean)
      ps = apply(test, 2, Compare.5p.3p)
      counts.ps = rbind(counts.ps, c(length(which(ps==1)), length(which(ps==0)), length(which(ps==(-1)))));
    }
  }
  colnames(counts.ps) = c("3p", "3/5p", "5p")
  rownames(counts.ps) = names
  
  length(which(counts.ps[,1]==0 | counts.ps[,3]==0))
  length(which(counts.ps[,1]<=1 | counts.ps[,3]<=1))
  
  mm = match(rownames(counts.ps), list.expressed)
  mm = which(!is.na(mm))
  
  length(which(counts.ps[mm,1]<=1 | counts.ps[mm,3]<=1))
  length(which(counts.ps[mm,1]==0 | counts.ps[mm,3]==0))
}


Compare.total.median.normalization = function()
{
  TEST.median.total.normalization = FALSE
  if(TEST.median.total.normalization)
  {
    #cat("nb of expressed miRNAs --", nrow(dds.expressed) )
    cat("size factor is -- ", sizeFactors(dds.expressed), "\n")
    ss = apply(counts(dds.expressed.mature), 2, sum)
    cat("size factor by total expressed mature read counts --- ", ss/mean(ss), "\n")
    ss = apply(counts(dds.expressed), 2, sum)
    cat("size factor by total expressed read counts --- ", ss/mean(ss), "\n")
    #dds <- estimateSizeFactors(dds)
    
    pdfname = paste0(resDir, "Comparision_total_median_normalization_expressed_mature_", specifity, ".pdf") #save all plots during data processing
    pdf(pdfname, width = 12, height = 6)
    cols = rep("blue", nrow(dds.expressed.mature));
    cols[grep("mmu-", rownames(dds.expressed.mature))] = 'red'
    par(mfrow=c(1,2))
    
    for(n in 1:6)
    {
      jj.test = c((2*n-1), 2*n)
      cpm = data.frame(fpm(dds.expressed.mature, robust = FALSE))
      plot(cpm[, jj.test], col=cols, log='xy', main='total normalization'); abline(0, 1, lwd=2.0)
      cpm = data.frame(fpm(dds.expressed.mature, robust = TRUE))
      plot(cpm[, jj.test], col=cols, log='xy', main='median normalization'); abline(0, 1, lwd=2.0)
    }
    dev.off()
  }
  
}

my.cpm.normalization = function(countData)
{
  cpm = matrix(NA, nrow = nrow(countData), ncol = ncol(countData));
  colnames(cpm) = colnames(countData)
  rownames(cpm) = rownames(countData)
  for(n in 1:ncol(countData))
  {
    cpm[,n] = countData[,n]/sum(countData[,n])*10^6
  }
  
  return(cpm)
}

calculate.scaling.factors.using.spikeIns = function(countData, concentrations = c(0.5, 2.5, 5.0, 15, 25, 35, 50, 250),
                                                    index.spikeIn=c(1:8), read.threshold=5, 
                                                    method="Ratio", plot.spikeIns.summary = TRUE)
{
  # counData = raw;
  if(method == "Ratio"){
    ss = apply(as.matrix(countData), 2, sum);
    ss.spikeIns = apply(as.matrix(countData[index.spikeIn, ]), 2, sum)
    cpm = my.cpm.normalization(countData)
    cpm.spikeIn = cpm[index.spikeIn, ]
    res = matrix(NA, ncol = ncol(cpm), nrow = nrow(cpm))
    rownames(res) = rownames(countData);
    colnames(res) = colnames(countData);
    counts.spikeIn = countData[index.spikeIn, ]
    
    if(nrow(cpm.spikeIn) != length(concentrations)) stop("wrong number of spike in or concentrations !!!")
    
    scaling.factors = rep(NA, ncol(cpm.spikeIn))
    for(n in 1:length(scaling.factors))
    {
      #n = 1;
      jj = which(counts.spikeIn[, n]>=read.threshold)
      if(length(jj)>=1)
      {
        scaling.factors[n] = median((concentrations[jj]) / (cpm.spikeIn[jj, n])); 
        res[,n] = cpm[,n]*scaling.factors[n];
        
        if(plot.spikeIns.summary){
          
          pie(c(ss.spikeIns[n], (ss[n]-ss.spikeIns[n])),labels = c('spikeIns', 'non-spikeIns'), 
              col=c('red', 'blue'), main=colnames(cpm)[n]) 
          
          #kk = grep("spikeIn", rownames(cpm))
          plot((cpm.spikeIn[, n]+10^-6),  concentrations, log='xy', 
               xlab="cpm", ylab="molecular per mug", main=paste0(colnames(cpm)[n], '--scalingFactor : ', signif(scaling.factors[n], d=2)), cex=3.0, col='darkgreen', pch= 16);
          abline(h=10^-6, lwd=2.0, col='darkgray')
          points(range((cpm.spikeIn[, n]+10^-6)), range((cpm.spikeIn[, n]+10^-6))*scaling.factors[n], type = 'l', lwd=3.0, col='darkblue', lty=1)
        }
        ## NOT use log scale any more 
        #if( ! log.scale ){
        #  fit = lm(concentrations[jj] ~ 0 + cpm.spikeIn[jj,n])
        #  norms[n] = fit$coefficients
        #}else{
        #  norms[n] = exp(median(log(concentrations[jj]) - log(cpm.spikeIn[jj, n])));
        #norms[n]
        #}
      }else{
        cat ("NO spike ins detected or pass the threshod of read number \n")
      }
    }
    
  }
  
  if(method == "DESeq2"){ # test DESeq2
    require('DESeq2')
    dds.spike= dds[kk, ]
    dds.spike <- dds.spike[ rowSums(counts(dds.spike)) > 5, ]
    dds.spike <- estimateSizeFactors(dds.spike)
    norms = sizeFactors(dds.spike);
  }
  
  if(method == "RUVg"){ # test RUVg method
    library(RUVSeq)
    zfGenes = assay(dds);
    filter <- apply(zfGenes, 1, function(x) length(x[x>5])>=2)
    filtered <- zfGenes[filter,]
    genes <- rownames(filtered)[grep("^cel", rownames(filtered))]
    spikes <- rownames(filtered)[grep("^spike", rownames(filtered))]
    
    x <- as.factor(design$genotype)
    set <- newSeqExpressionSet(as.matrix(filtered), phenoData = data.frame(x, row.names=colnames(filtered)))
    set
    
    library(RColorBrewer)
    colors <- brewer.pal(3, "Set2")
    plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    plotPCA(set, col=colors[x], cex=1.2)
    
    set <- betweenLaneNormalization(set, which="upper")
    plotRLE(set, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    plotPCA(set, col=colors[x], cex=1.2)
    
    set1 <- RUVg(set, spikes, k=2)
    pData(set1)
    plotRLE(set1, outline=FALSE, ylim=c(-4, 4), col=colors[x])
    plotPCA(set1, col=colors[x], cex=1.2)
  }
  
  ## compute the size factor for DESeq2
  ss.geo.mean = prod(ss)^(1/length(ss))
  norms = 1/scaling.factors * ss / ss.geo.mean
  
  return(list(scaling.factors = scaling.factors, norms4DESeq2 = norms, cpm = cpm, normalization.spikeIn = res))
  
}

## this function is to merge technical replicates using all and design as inputs
## TODO improve this funciton, because it is not quite clear for the return  
Merge.techinical.replicates.using.design.countTable = function(design, all, id.list=list(c("57751", "57753"), c("57752", "57754")))
{
  for(n in 1:length(rep.technical))
  {
    index = c()
    for(id in rep.technical[[n]])
    {
      #print(id)
      index = c(index, which(design$SampleID==id))
    }
    
    design$SampleID[index[1]] = paste0(design$SampleID[index], collapse = ".")
    ss = apply(all[, (index+1)], 1, function(x) sum(x, na.rm = TRUE))
    all[, (index[1]+1)] = ss;
    colnames(all)[(index[1]+1)] = paste0(design$genotype[index[1]], "_", design$tissue.cell[index[1]], "_", design$treatment[index[1]], "_",  design$SampleID[index[1]])
    design = design[-index[-1], ]
    all = all[, -(index[-1]+1)]
  }
  
  return(list(design = design, 
              countTable = all))
}

Merge.techinical.replicates = function(stats, rep.technical = list(c("57751", "57753"), c("57752", "57754")))
{
  for(n in 1:length(rep.technical))
  {
    index = c()
    for(id in rep.technical[[n]])
    {
      #print(id)
      index = c(index, which(colnames(stats)==id))
    }
    
    #design$SampleID[index[1]] = paste0(design$SampleID[index], collapse = ".")
    ss = apply(stats[, index], 1, function(x) sum(x, na.rm = TRUE))
    stats[, (index[1])] = ss;
    colnames(stats)[index[1]] = paste0(colnames(stats)[index], collapse = ".")
    #paste0(design$genotype[index[1]], "_", design$tissue.cell[index[1]], "_", design$treatment[index[1]], "_",  design$SampleID[index[1]])
    #design = design[-index[-1], ]
    stats = stats[, -index[-1]]
  }
  
  return(stats)
  
}

######################################
######################################
## Section: functions to identify expressed miRNAs and also mature arms for those expressed miRNAs
######################################
######################################
find.expressed.mature.miRNA.using.cpm.threshold = function(countData, cpm.threshold=10, species = 'cel')
{
  #############################
  ## here we are not using DESeq2 to calculate the cpm  
  ## we are using the countData to calculate the mean and then determine if they are expressed or not 
  #############################
  # countData = all[, (kk+1)]; cpm.threshold=10; species = "cel";
  raw = as.matrix(countData)
  raw[which(is.na(raw))] = 0
  cpm = my.cpm.normalization(raw)
  colnames(cpm) = paste0(colnames(cpm), '.cpm')
  ## filter lowly expressed miRNAs using cpm=10 and select the mature arms always baesed on the untreated samples
  
  #cpm = fpm(dds, robust = FALSE)
  if(species == "cel"){
    ggs = sapply(rownames(cpm), find.mirName)
  }else{
    if(species == "dm"){
      ggs = sapply(rownames(cpm), find.mirName.dm)
    }else{
      cat("unknown species ---\n")
    }
  }
  
  if(ncol(cpm)>1) {
    mean.cpm = apply(cpm, 1, mean);
  }else{
    mean.cpm = cpm;
  }
  hist(log2(mean.cpm), breaks = 40)
  abline(v= log2(cpm.threshold), lwd=2.0, col='blue')
  gene.expressed = unique(ggs[which(mean.cpm > cpm.threshold)])
  
  #cpm.untreated = cpm[, which(dds$treatment=="untreated")]
  #cpm.treated = cpm[, which(dds$treatment=="treated")]
    
  #if(length(which(dds$treatment=="treated"))>1) {mean.treated = apply(cpm[, which(dds$treatment=="treated")], 1, mean);
  #}else {mean.treated = cpm[, which(dds$treatment=="treated")];}
  
  #if(sampletoUse=="untreated") 
  #if(sampletoUse=="treated") gene.expressed = unique(ggs[which(mean.treated>cpm.threshold)])
  #if(sampletoUse=="all"){
  #  mean.all = apply(as.matrix(cbind(mean.treated, mean.untreated)), 1, mean);
  #  gene.expressed = unique(ggs[which(mean.all>cpm.threshold)])
  #} 
  #if(sampletoUse=="untreated.or.treated") gene.expressed = unique(ggs[which(mean.treated>cpm.threshold|mean.untreated>cpm.threshold)])
  index.expressed = which(!is.na(match(ggs, gene.expressed)))
  
  mirna.all = data.frame(rownames(cpm), ggs, cpm, mean=mean.cpm)
  mirna.expressed = mirna.all[index.expressed, ]
  colnames(mirna.expressed)[c(1:2)] = c("miRNA", "gene")
  
  expressed.mature.mirna = find.mature.ones.for.prefixed.expressed.miRNAs(mirna.expressed)
  
  return(expressed.mature.mirna)
  #return(list(miRNA=sels$miRNA, expressed=sels$expressed, cpm=cpm))
}

## updated version of identify.expressed.miRNAs for time series or stages 
## this function is to identify the expressed miRNAs using the cpm threshold
identify.expressed.miRNAs.for.stages = function(countData, design.matrix, cpm.threshold=10) 
{
  cat('identify list of expressed miRNAs for each stage\n')
  require('DESeq2')
  dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ treatment + stage)
  cpm =  fpm(dds, robust = FALSE)
  
  colnames(cpm) = paste0(colnames(cpm), ".cpm")
  
  test.theshold = FALSE
  if(test.theshold){
    par(mfrow=c(4,4))
    for(n in 1:ncol(cpm))
    {
      hist(log2(cpm[which(cpm[,n]>0), n]), breaks = 30,  xlab="log2(cpm)", main=colnames(cpm)[n], cex=3.0, col='darkgray', pch= 16);
      ggs = rownames(cpm)[which(cpm[,n]>10)]
      ggs = ggs[grep('spikeIn_', ggs, invert = TRUE)]
      ggs = sapply(ggs, function(x) gsub("-3p", '', x))
      ggs = sapply(ggs, function(x) gsub("-5p", '', x))
      cat("nb of miRNAs above the threshold--", length(unique(ggs)), "\n")
      abline(v=log2(10), lwd=2.0, col='red')
      #points(range((cpm[kk, n]+10^-6)), range((cpm[kk, n]+10^-6))/norms[n], type = 'l', lwd=3.0, col='darkblue', lty=1)
    }
  }
  
  cpm.untreated = cpm[, grep("_untreated", colnames(cpm))]
  kk.spikeIn = grep("spikeIn_", rownames(cpm.untreated))
  if(length(kk.spikeIn)>0) cpm.untreated = cpm.untreated[-kk.spikeIn, ]
  
  times = unique(design.matrix$stage[which(design.matrix$treatment=="untreated")])
  compare.threshold = c()
  for(n in 1:nrow(cpm.untreated))
  {
    test = 0;
    for(t in times) 
    {
      #cat(test, "---", mean(cpm.untreated[n, grep(t, colnames(cpm.untreated))]))
      if(mean(cpm.untreated[n, grep(t, colnames(cpm.untreated))])>cpm.threshold)  test = test + 1;
    }
    compare.threshold = c(compare.threshold, test)
  }
  
  expressed = data.frame(rownames(cpm.untreated), sapply(rownames(cpm.untreated), find.mirName), compare.threshold, cpm.untreated,
                         stringsAsFactors = FALSE)
  colnames(expressed)[1:3] = c('miRNA', 'gene', 'nb.stage.above.threshold' )
  gg.expressed = unique(expressed$gene[which(expressed$nb.stage.above.threshold>0)])
  
  mm = match(expressed$gene, gg.expressed)
  expressed = expressed[which(!is.na(mm)==TRUE), ]
  
  return(expressed)
  
}

find.mature.ones.for.prefixed.expressed.miRNAs = function(list.expressed.miRNAs)
{
  list.expressed = list.expressed.miRNAs
  colnames(list.expressed)[c(1,2)] = c("miRNA", "gene")
  list.expressed = list.expressed[order(list.expressed$miRNA), ]
  ggs.uniq = unique(list.expressed$gene)
  
  mature = rep(NA, nrow(list.expressed))
  kk = grep('.cpm', colnames(list.expressed))
  
  for(n in 1:length(ggs.uniq))
  {
    jj = which(list.expressed$gene==ggs.uniq[n])
    if(length(jj)>1){
      index.max = apply(list.expressed[jj, kk], 2, which.max)
      nb.first.max = length(which(index.max==1))
      nb.sec.max = length(which(index.max==2))
      #cat(n, ": ", as.character(ggs.uniq[n]), "--",  nb.first.max, "--", nb.sec.max, "\n")
      if(nb.first.max>nb.sec.max){
        mature[jj[1]] = TRUE; mature[jj[2]] = FALSE;  
      }else{
        mature[jj[1]] = FALSE; mature[jj[2]] = TRUE; 
      }
    }else{
      #cat(n,": ", as.character(ggs.uniq[n]),  "-- no selection \t")
      mature[jj] = TRUE;
    }
  }
  expressed.miRNAs = data.frame((list.expressed[, c(1, 2)]), mature=(mature), list.expressed[, -c(1:2)], stringsAsFactors = FALSE)
  return(expressed.miRNAs)
}

######################################
######################################
## Section: other functions
######################################
######################################
find.replicates.by.removing.ID = function(x)
{
  infos = unlist(strsplit(x, "_"))
  return(paste0(infos[-length(infos)], collapse = '_'))
}

average.biological.replicates = function(cpm)
{
  # cpm = cpm.piRNA.bc;
  samples = sapply(colnames(cpm), find.replicates.by.removing.ID, USE.NAMES = FALSE)
  samples.uniq = unique(samples)
  
  if(length(samples.uniq) == length(samples)){
    cat('---no replicates exist---')
  }else{
    cpm.mean = matrix(NA, nrow = nrow(cpm), ncol=length(samples.uniq))
    rownames(cpm.mean) = rownames(cpm)
    colnames(cpm.mean) = samples.uniq;
    
    for(n in 1:ncol(cpm.mean))
    {
      kk = which(samples == colnames(cpm.mean)[n])
      if(length(kk)>1){
        cpm.mean[ ,n] = apply(as.matrix(cpm[,kk]), 1, median)
      }else{
        if(length(kk)==1) cpm.mean[, n] = cpm[,kk]
      }
    }
  }
  
  return(cpm.mean)
}

variance.biological.replicates = function(cpm)
{
  # cpm = cpm.piRNA.bc.prot;
  samples = sapply(colnames(cpm), find.replicates.by.removing.ID, USE.NAMES = FALSE)
  samples[grep("_untreated", samples)] = "whole.body"
  samples.uniq = unique(samples)
  
  if(length(samples.uniq) == length(samples)){
    cat('---no replicates exist---')
  }else{
    cpm.mean = matrix(NA, nrow = nrow(cpm), ncol=length(samples.uniq))
    rownames(cpm.mean) = rownames(cpm)
    colnames(cpm.mean) = samples.uniq;
    
    cpm.vars = cpm.mean;
    
    for(n in 1:ncol(cpm.mean))
    {
      kk = which(samples == colnames(cpm.mean)[n])
      if(length(kk)>1){
        cpm.mean[ ,n] = apply(as.matrix(cpm[,kk]), 1, median)
        cpm.vars[, n] = apply(as.matrix(cpm[,kk]), 1, var)
      }else{
        if(length(kk)==1) cpm.mean[, n] = cpm[,kk]
      }
    }
  }
  
  cpm.vars[which(cpm.vars<=0)] = NA;
  
  for(n in 1:ncol(cpm.mean)){
    x = (as.vector(cpm.mean[,n]))
    y = (as.vector(cpm.vars[,n]))
    jj = which(!is.na(x) & !is.na(y))
    x = x[jj]; y = y[jj];
    plot((x), (y), log='xy', cex=0.5)
    library(MASS);
    library(stats)
    fit = loess(y ~ x)
    #abline(fit$coefficients[1], fit$coefficients[2])
    # points(range(exp(x)), range(exp(x)^fit$coefficients[2] * exp(fit$coefficients[1])), type = "l", col = 'red')
    points(range((x)), predict(fit, range(x)), type = "l", col = 'red', lwd =2 )
    cat(fit$coefficients[1], '--',  fit$coefficients[2], "\n")
  }
  
  Check.variance.across.samples = FALSE
  if(Check.variance.across.samples){
    xlim = range(cpm.mean)
    ylim = range(cpm.vars, na.rm = TRUE)
    
    plot(1, 1, type = "n", xlim = xlim, ylim = ylim, log = "xy", xlab = "mean in liner scale", ylab = "variance in linear scale")
    for(n in 1:ncol(cpm.vars)){
      points(cpm.mean[,n], cpm.vars[,n], col = n, cex=0.4)
    }
    points(xlim, (xlim^2*0.1), type = "l", col = 'darkblue', lwd=2.0)
    
  }
  
  
  cpm.vars.hat = cpm.mean^fit$coefficients[2] * exp(fit$coefficients[1]) 
  return(cpm.vars.hat)
}


########################################################
########################################################
# Section: functions to correct batch 
########################################################
########################################################
remove.batch.by.logratios = function(cpm, design.matrix)
{
  # cpm = cpm.piRNA
  yy = cpm
  
  for(n in 1:(nrow(design.matrix)/2))
  {
    jj = c((2*n-1), 2*n) 
    jj.ref = jj[which(design.matrix$treatment[jj] == "untreated")]
    jj.treated = setdiff(jj, jj.ref)
    cat(colnames(cpm)[jj.ref], "-- vs. --", colnames(cpm)[jj.treated], "\n")
    
    yy[, jj.treated] = log2((yy[ ,jj.treated]+2^-6)/(yy[, jj.ref]+2^-6));
    
  }
  
  yy = yy[, which(design.matrix$treatment=="treated")]
  
  
}


remove.batch.using.N2.untreated = function(cpm, design.matrix, method = "linear.model")
{
  # cpm = cpm.piRNA
  logcpm = log2(cpm + 2^-6)
  
  if(method == 'linear.model'){
    ## use the linear model for data the log2 scale (sample by sample, NOT gene by gene) to remove the batch effect; use N2 untreated condtion as references
    ## (NOT used here)
    cat('Warnings -- use limma or ComBat instead of this method ')
    
    #reference = "N2_whole.body_untreated"
    #genotypes = sapply(colnames(cpm), function(x) unlist(strsplit(x, "_"))[1], USE.NAMES = FALSE)
    #ns =  sapply(colnames(cpm), function(x) unlist(strsplit(x, "_"))[2], USE.NAMES = FALSE)
    #treatment = sapply(colnames(cpm), function(x) unlist(strsplit(x, "_"))[3], USE.NAMES = FALSE)
    pheno = data.frame(design.matrix, stringsAsFactors = FALSE)
    
    jj.N2.untreated = which(pheno$genotype=="N2" & pheno$treatment == "untreated" & pheno$batch == 1)
    y = apply(logcpm[, jj.N2.untreated], 1, mean)
    batchs = unique(pheno$batch)
    batchs = batchs[which(batchs != 1)]
    for(n in 1:length(batchs))
    {
      #n = 16;
      jj.untreated = which(pheno$batch == batchs[n] & pheno$treatment == "untreated")
      if(length(jj.untreated)>0)
      {
        if(length(jj.untreated) == 1) x = logcpm[, jj.untreated]
        if(length(jj.untreated)>1) x = apply(logcpm[, jj.untreated], 1, mean)
       
        fit = lm(y ~ x)
        plot(x, y, xlab = "untreated", ylab = "N2");
        abline(fit, col='red', lwd=2.0)
        abline(0, 1, col='blue', lwd=2.0)
        jj.to.correct = which(pheno$batch == batchs[n])
        for(j in jj.to.correct) {
          logcpm[, j]  = logcpm[, j]*fit$coefficients[2] + fit$coefficients[1];
          cat(colnames(cpm)[j], "-intercept-", fit$coefficients[1], " - slop -", fit$coefficient[2],  "\n")   
        }
      }
    }
    logcpm.bc = logcpm;
  }
  
  if(method == 'limma'){
    
    cat('remove the batch effect using limma \n')
    require('limma')
    design.tokeep = design.matrix
    design.tokeep$tissue.cell[which(design.tokeep$treatment == "untreated")] = 'whole.body'
    design.tokeep$tissue.cell[which(design.tokeep$genotype=="N2" & design.tokeep$treatment=="treated")] = "background"
    design.tokeep<-model.matrix(~0 + tissue.cell,  data = design.tokeep)
    logcpm.bc = removeBatchEffect(logcpm, batch = design.matrix$batch, design = design.tokeep)
    
  }
  if(method == 'combat'){
    ## here we use the combat to remove the batch effect 
    ## the combat requires the log2cpm
    cat('remove the batch effect using ComBat \n')
    require("sva")
    # example from the ComBat function in the R package 'sva'
    TEST.example = FALSE
    if(TEST.example){
      library(bladderbatch)
      data(bladderdata)
      dat <- bladderEset[1:50,]
      pheno = pData(dat)
      edata = exprs(dat)
      batch = pheno$batch
      mod = model.matrix(~as.factor(cancer), data=pheno)
      combat_edata3 = ComBat(dat=edata, batch=batch, mod=mod, par.prior=TRUE, ref.batch=3)
    }
    
    batch = design.matrix$batch;
    design.tokeep = design.matrix
    design.tokeep$tissue.cell[which(design.tokeep$treatment == "untreated")] = 'whole.body'
    design.tokeep$tissue.cell[which(design.tokeep$genotype=="N2" & design.tokeep$treatment=="treated")] = "background"
    mod = model.matrix(~ as.factor(tissue.cell), data = design.tokeep);
    #conds = data.frame(rep(c("untreated", "treated"), ncol(cpm)/2))
    #colnames(conds) = 'treatment'
    #mod = model.matrix(~ as.factor(treatment), conds)
    logcpm.bc = ComBat(dat=logcpm, batch=batch, mod=mod, par.prior=TRUE, ref.batch = NULL)
  }
  
  return(2^logcpm.bc)
}

Test.piRNA.normalization.batch.removal = function(cpm, design.matrix)
{
  ## here is a function to test piRNA normalization and batchRemoval for the gene expression matrix in deconvolution analysis
  ## there will be one PCA plot and four samll plots for the positive controls
  # cpm = cpm.piRNA
  main.names = deparse(substitute(cpm)) 
  if(any(cpm==0)) cpm = cpm + 2^-6 
  
  require(lattice);
  require(ggplot2)
  require('DESeq2');
  library("vsn");
  library("pheatmap");
  library("RColorBrewer");
  library("dplyr");
  
  ## pca plots only for untreated samples
  sels = which(design.matrix$treatment == "untreated")
  pca = prcomp(t(log2(cpm[, sels])), scale. = TRUE, center = TRUE)
  pca2save = data.frame(pca$x, condition=design.matrix$treatment[sels], 
                        batch = design.matrix$batch[sels], 
                        name=colnames(cpm)[sels], 
                        tissue = design.matrix$tissue.cell[sels])
  
  #ggp = ggplot(data=pca2save, aes(PC1, PC2, label = batch, color = tissue)) + geom_point(size=4) +
  #  geom_text(hjust = 0.1, nudge_y = 0.2, size=5) +
  #  ggtitle(paste0("PCA - ", main.names))
  #plot(ggp);
  #pca=plotPCA(cpm, intgroup = colnames(design.matrix)[c(3, 5, 7)], returnData = FALSE)
  #print(pca)
  
  ## ratio between treated and untreated to test if piRNA normalization makes sense
  ns = unique(design.matrix$tissue.cell)
  ratios = matrix(NA, ncol = length(ns), nrow = nrow(cpm))
  rownames(ratios) = rownames(cpm);
  colnames(ratios) = ns
  total = ratios;
  tcs = ratios;
  for(n in 1:length(ns))
  {
    #jj = which(design.matrix$tissue.cell==ns[n])
    total[,n] = apply(cpm[, which(design.matrix$tissue.cell==ns[n] & design.matrix$treatment=="untreated")], 1, mean)
    tcs[,n] = apply(cpm[, which(design.matrix$tissue.cell==ns[n] & design.matrix$treatment=="treated")], 1, mean)
    ratios[,n] = tcs[,n]/total[,n] 
  }
  
  ###############################
  # check ratios between treated and untreated samples and other two samples
  ###############################
  #par(mfrow=c(1, 3))
  #hist(log10(ratios), xlab = "log10(treated/untreated)", main = paste0(main.names), breaks = 100)
  #abline(v=c(-3:0), col='darkred', lwd=2.0)
  
  #plot(tcs[,which(colnames(tcs) == "Dopaminergic.neurons")], tcs[,which(colnames(tcs) == "Ciliated.sensory.neurons")], log='xy', 
  #     xlab='Dopaminergic (in log)', ylab='Ciliated (in log)', main = main.names)
  #abline(0, 1, lwd=2.0, col='red')
  #plot(tcs[,which(colnames(tcs) == "mechanosensory.neurons" )], tcs[,which(colnames(tcs) == "unc-86.expressing.neurons")], log='xy', 
  #     xlab='mechanosensory (in log)', ylab='unc-86 (in log)', main = main.names)
  #abline(0, 1, lwd=2.0, col='red')
  
  ###############################
  # check three examples lsy-6, mir-791 and mir-790 across all samples
  # different replicates should be also displayed
  ###############################
  par(mfrow=c(1, 3))
  
  for(gg in c("lsy-6", "mir-791", "mir-790"))
  {
    #gg = "lsy-6"
    par(cex =0.7, mar = c(8,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    ## check lsy-6 in ASE, Glutatamatergic and Ciliated neurons
    kk = which(rownames(cpm)==gg)
    
    for(m in c(1)){
      if(m == 1) {log = 'y'; lims = range(cpm[kk,]);}
      if(m == 2) {log = ''; lims = range(cpm[kk,]);}
      if(m == 3) {log = ''; lims = range(cpm[kk, which(design.matrix$treatment=="treated")]);}
      plot(c(1:length(ns)), tcs[kk, ], type= 'l', col='darkblue', cex=1.0, log=log, ylim =lims, main = paste0(gg ," in ", main.names), xlab=NA, 
           ylab = 'normalizaed by piRNAs', axe = FALSE)
      #points(c(1:length(ns)), tcs[kk, ], type = "b", cex=1.0, col = 'darkblue')
      points(c(1:length(ns)), total[kk, ], type= 'l', col='black', cex=1.0)
      #points(c(1:length(ns)), total[kk, ], type= 'b', col= "black", cex= 1.0)
      
      for(n in 1:length(ns))
      {
        index.ns = which(design.matrix$tissue.cell==ns[n] & design.matrix$treatment=="treated")
        points(rep(n, length(index.ns)), cpm[kk, index.ns], type = "p", col='darkblue', cex=1.5, pch =16)
        # add sample ids 
        if(ns[n]=="Pan.neurons"){
          text(rep(n, length(index.ns)), cpm[kk, index.ns], design.matrix$SampleID[index.ns], pos = 2, offset = 0.5, cex = 0.8)
        }
        
        index.ns = which(design.matrix$tissue.cell==ns[n] & design.matrix$treatment=="untreated")
        points(rep(n, length(index.ns)), cpm[kk, index.ns], type = "p", col='black', cex=1.5, pch = 0)
      }
      #legend("topright", col=c('darkblue', "black"),  bty = "n", legend = c("treated", "untreated"), lty=1 )
      axis(2, las= 1)
      ns.short = sapply(ns, function(x) unlist(strsplit(x, "[.]"))[1], USE.NAMES = FALSE)
      ns.short[c(3:5, 8,10, 11)] = c("Seroton", "Dopamin", "Gluta", "mecha", 'pharyn', "cholin")
      axis(1, at=c(1:length(ns)), labels = ns.short, las=2,cex=0.5)
      box()
      abline(h=c(50, 100, 200, 500, 1000), lwd=0.7, col='red')
      
    }
  }
  
}


calculate.pvalues.two.groups.overlapping = function(nb.total, nb.group.A, nb.group.B, nb.overlapping)
{
  total = as.numeric(nb.total)
  q = as.numeric(nb.overlapping)
  m = as.numeric(nb.group.A)
  nn = total - m;
  k = as.numeric(nb.group.B)
  rr = q/(k*m/total)
  pvals = phyper((q-1), m, nn, k, lower.tail = FALSE, log.p = FALSE)
  
  return(c(rr, pvals))
}


####################
## function to merge the neuron classes for fraction matrix 
####################
proportions.matrix.merging.neuronClass = function(proportions, using.binary.matrix = TRUE, Test.by.plot = TRUE)
{
  cat("merge the neuron classes if they are not distinguishable ...\n")
  
  # Test.by.plot = TRUE
  xx = proportions;
  xx[which(xx>0)] = 1
  
  mydata = t(xx)
  mydata <- na.omit(mydata) # listwise deletion of missing
  #mydata <- scale(mydata) # standardize variables 
  
  d <- dist(mydata, method = "euclidean") # distance matrix
  fit <- hclust(d, method="ward.D2")
  groups <- cutree(fit, h=0) # cut tree into 5 clusters
  nb.clusters = length(unique(groups))
  
  # start to merge neuron classes if they have zero distance in the presence-absence matrix
  newdata = matrix(NA, ncol = nb.clusters, nrow = nrow(proportions))
  rownames(newdata) = rownames(proportions)
  colnames(newdata) = rep('x', ncol(newdata))
  #names.merged = c()
  for(n in 1:nb.clusters)
  {
    # n = 1
    kk = which(groups == n)
    cat("cluster", n, "--", names(groups)[kk],"\n")
    colnames(newdata)[n] = paste0(names(groups)[kk], collapse = ".")
    if(length(kk)>1){
      newdata[,n] = apply(proportions[, kk], 1, sum)
    }else{
      if(length(kk)==0){
        cat("Error-- no neuron classes found \n")
      }else{
        newdata[,n] = proportions[, kk]
      }
    }
  }
  
  if(Test.by.plot){
    library("pheatmap")
    library("RColorBrewer")
    
    pheatmap(xx, 
             cluster_rows=FALSE, 
             show_rownames=TRUE, show_colnames = TRUE,
             cluster_cols=TRUE, 
             color = c("lightgray", "blue"), legend = FALSE)
    
    
    plot(fit, cex=0.8) # display dendogram
    
    # draw dendogram with red borders around the 5 clusters
    rect.hclust(fit, k=nb.clusters, border="red")
    abline(h= 0, col='darkblue', lwd=2.0)
    #abline(h= 1, col='darkblue', lwd=2.0, lty=2.0)
    
    yy = newdata;
    yy[which(yy>0)] = 1
    
    pheatmap(yy, 
             cluster_rows=FALSE, 
             show_rownames=TRUE, show_colnames = TRUE,
             cluster_cols=TRUE, 
             color = c("lightgray", "blue"), legend = FALSE)
    
    
  }
  
  return(newdata)
  
}

expressionMatrix.grouping = function(xx, using.logscale = TRUE)
{
  expression = xx[, -c(1, 2)]
  for(n in 1:ncol(expression)) 
  {
    expression[,n] = log2(expression[,n]/xx$background)
  }
  
  iris.scaled <- scale(expression)
  
  require(cluster)
  fviz_nbclust((iris.scaled), hcut, method = "silhouette",
               hc_method = "complete", k.max = 120)
  
  library(cluster)
  #set.seed(123)
  #gap_stat <- clusGap(iris.scaled, FUN = kmeans, nstart = 25,
  #                    K.max = 10, B = 50)
  
  set.seed(123)
  gap_stat <- clusGap(iris.scaled, FUN = hcut, K.max = 120, B = 50)
  #print(gap_stat, method = "firstmax")
  
  # Plot gap statistic
  plot(gap_stat, frame = FALSE, xlab = "Number of clusters k")
  abline(v = 3, lty = 2)
  
  fviz_gap_stat(gap_stat)
  
  # Print the result
  print(gap_stat, method = "firstmax")
  
}

calculate.sizeFactors4piRNAs = function(read.count, design.matrix, lowlyExpressed.readCount.threshold = 10)
{
  require('DESeq2');
  if(ncol(design.matrix)>2){cc = apply(design.matrix[, -1], 1, paste0, collapse="_")
  }else{cc = design.matrix[, -1]}
  
  raw = as.matrix(read.count)
  raw[which(is.na(raw))] = 0
  
  countData = ceiling(raw)
  conds = factor(paste0(colnames(design.matrix)[-1], collapse = " + "))
  eval(parse(text = paste0("dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ ", conds, ")")))
  #dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ condition + time)
  if(is.null(lowlyExpressed.readCount.threshold))  lowlyExpressed.readCount.threshold = 10
  dds <- dds[ rowSums(counts(dds)) >= lowlyExpressed.readCount.threshold, ]
  
  dds <- estimateSizeFactors(dds)
  
  return(sizeFactors(dds))
  
}

Compare.piRNA.siRNA.spikeIns.for.scaling.factors = function(library.sizes, stats, countData, design.matrix, 
                                                            compareSpikeIn = TRUE, cpm.pairwise.compare=FALSE, 
                                                            piRNA.sizeFctors=NULL)
{
  ###############################
  ## check the correlation of stats
  ###############################
  if(is.null(piRNA.sizeFctors)){
    par(mfrow=c(1, 1))
    library(corrplot)
    col1 <- colorRampPalette(c("#7F0000","red","#FF7F00","yellow","white", 
                               "cyan", "#007FFF", "blue","#00007F"))
    xx = as.matrix(stats)
    xx[which(xx==0)] = NA
    M <- cor(xx, use = "na.or.complete")
    #corrplot(M, method="circle", type = 'upper', order="hclust")
    corrplot(M, method="ellipse", order="hclust", tl.cex=1.2, cl.cex=0.7, tl.col="black", 
             addrect=ceiling(ncol(xx)/2), col=col1(100), rect.col=c('green'), rect.lwd=2.0)
    
    par(cex =0.7, mar = c(5,5,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
    par(mfrow=c(2, 2))
    plot(stats$piRNA, stats$siRNA, log='xy', xlab = "piRNA total reads", ylab = "siRNA total reads");
    abline(log10(median(stats$siRNA/stats$piRNA)), 1, lwd=2.0, col='red')
    
    plot(stats$piRNA, stats$piRNA_AS, log='xy', xlab = "piRNA total reads", ylab = "piRNA_AS total reads");
    abline(log10(median(stats$piRNA_AS/stats$piRNA)), 1, lwd=2.0, col='red')
    
    #source('RNAseq_Quality_Controls.R')
    #pairs(stats, lower.panel=NULL, upper.panel=panel.fitting)
    plot(library.sizes, stats$piRNA, log = 'xy', xlab = "miRNA total reads", ylab="piRNA total reads")
    plot(library.sizes, stats$siRNA, log = 'xy', xlab = "miRNA total reads", ylab="siRNA total reads")
    
    cc = apply(design.matrix[, c(3,5)], 1, paste0, collapse="_");
    cc.uniq = unique(cc);
    cols = match(cc, cc.uniq)
    
    par(cex = 1.8, las = 1, mgp = c(1.6,0.5,0), mar = c(6,18,2,0.8)+0.1, tcl = -0.3)
    par(mfrow=c(1, 1))
    barplot(library.sizes/stats$good_reads*100, horiz = TRUE, xlim = c(0, 100),
            names.arg = colnames(countData), las=1, col = cols, 
            main='miRNA reads/good reads', xlab='percentages (%)')
    
    barplot(stats$piRNA/stats$good_reads*100, horiz = TRUE, xlim = c(0, 100),
            names.arg = colnames(countData), las=1, col = cols, 
            main='piRNA reads/good reads', xlab='percentages (%)')
    
    barplot(stats$siRNA/stats$good_reads*100, horiz = TRUE, xlim = c(0, 100),
            names.arg = colnames(countData), las=1, col = cols, 
            main='siRNA reads/good reads', xlab='percentages (%)')
  }
  ###############################
  # compare spike-in normalization and piRNA-normalization 
  ###############################
  if(compareSpikeIn){
    spike.files = list.files(path = "../data/normalized_piRNAs", pattern = "spikeIns_count_table.txt", full.names = TRUE)
    spikes = NULL;
    for(j in 1:length(spike.files)){
      if(j == 1){
        spikes = read.delim(spike.files[j], sep="\t", header = TRUE, row.names = 1)
      }else{
        spikes = rbind(spikes, read.delim(spike.files[j], sep="\t", header = TRUE, row.names = 1))
      }
    }
    spikes = t(as.matrix(spikes))
    spikes = data.frame(gene=rownames(spikes), spikes, stringsAsFactors = FALSE) 
    
    index.samples.with.spikeIns = which(design.matrix$tissue.cell=="Pan.neurons" & design.matrix$genotype=="WT")
    spikes = process.countTable(all=spikes, design = design.matrix[index.samples.with.spikeIns, ], select.Total.count = FALSE)
    rownames(spikes) = spikes$gene;
    spikes = spikes[, -1]
    
    mm = match(colnames(spikes), colnames(countData))
    
    pans = rbind(as.matrix(spikes), as.matrix(countData[,mm]));
    stats.pans = stats[mm, ]
    library.sizes.pans = library.sizes[mm]
    index.spikeIn = grep("spikeIn", rownames(pans))[c(1:8)]
    design.matrix.pans = design.matrix[mm, ]
    piRNA.sizeFctors.pans = piRNA.sizeFctors[mm]
    
    # here the concentration is amol per mug of total RNA
    concentrations = c(0.05, 0.25, 0.5, 1.5, 2.5, 3.5, 5, 25)*100
    
    ## calculate scaling factor using spike-ins
    #source("miRNAseq_functions.R")
    par(mfrow=c(2,2))
    res.spike.in = calculate.scaling.factors.using.spikeIns(pans, concentrations = concentrations, 
                                                            index.spikeIn = index.spikeIn, read.threshold = 10)
    cpm = res.spike.in$cpm;
    res = res.spike.in$normalization.spikeIn
    colnames(cpm) = paste0(colnames(cpm), ".cpm")
    colnames(res) = paste0(colnames(res), ".normSpike")
    
    ## double check total number of reads with and without spike-ins
    par(mfrow=c(2,2))
    ss = apply(pans[-c(1:8),], 2, sum)
    plot(ss, library.sizes.pans, log='', xlab='read nbs of expressed miRNAs', ylab='read nbs of all miRNAs')
    abline(0, 1, lwd=2.0, col='red')
    
    ss = apply(pans, 2, sum);
    plot(ss, library.sizes.pans, log='', xlab='read nbs of expressed miRNAs + spikeIns', ylab='read nbs of all miRNAs')
    abline(0, 1, lwd=2.0, col='red')
    
    ## double check the cpm and spike-in normalization
    ss = apply(pans, 2, sum)
    #par(mfrow=c(1,2))
    plot(pans[,1]/ss[1]*10^6, cpm[,1], log='xy', main = 'confirming cpm caclulation'); abline(0, 1, lwd=2.0, col='red')
    plot(pans[,1]/ss[1]*10^6*res.spike.in$scaling.factors[1], res[,1], log='xy', main = "confirming spike-in normalization"); 
    abline(0, 1, lwd=2.0, col='red')
    
    #plot(raw[,1]/ss[1]*10^6/norms[1], res[,1], log='xy');abline(0, 1, lwd=2.0, col='red')
    sf.p = stats.pans$piRNA/10^6
    sf.s = res.spike.in$norms4DESeq2
    #sf.s = sf.s/sf.s[1]
    if(is.null(piRNA.sizeFctors)){
      par(mfrow=c(1,1))
      plot(sf.s, sf.p, log='xy', main = "spikeIns norm vs piRNAs norm", col='darkgreen', pch=16, cex=1.5, xlim = c(4, 400), ylim = c(0.2, 50),
           xlab= 'sizeFactor.spikeIns', ylab='sizeFactor.piRNAs')
      abline(log10(median(sf.p/sf.s)), 1, lwd=2.0, col="red")
      text(sf.s, sf.p, labels = colnames(pans), offset = 0.5, pos = 1, cex = 0.7)
      
      kk = c(grep("71822", colnames(pans)), grep("71823", colnames(pans)))
      
      plot(sf.s[-kk], sf.p[-kk], log='xy', main = "spikeIns norm vs piRNAs norm", col='darkgreen', pch=16, cex=1.5, 
           xlim = c(4, 400), ylim = c(0.2, 50),
           xlab= 'sizeFactor.spikeIns', ylab='sizeFactor.piRNAs')
      abline(log10(median(sf.p/sf.s)[-kk]), 1, lwd=2.0, col="red")
      text(sf.s[-kk], sf.p[-kk], labels = colnames(pans)[-kk], offset = 0.5, pos = 1, cex = 0.7)
      cat("correction between spike-in and piRNA -- ", cor(sf.s[-kk], sf.p[-kk]), "\n")
      
      
    }else{
      par(mfrow=c(2,2))
      sf.pp = piRNA.sizeFctors.pans
      plot(sf.s, sf.p, log='xy', main = "spikeIns norm vs piRNAs library", col='darkgreen', pch=16, cex=1.5, xlim = c(4, 400), ylim = c(0.2, 100),
           xlab= 'sizeFactor.spikeIns', ylab='piRNA library.size')
      #abline(0, (median(sf.p/sf.s)), lwd=2.0, col="red")
      text(sf.s, sf.p, labels = colnames(pans), offset = 0.5, pos = 1, cex = 0.7)
      
      plot(sf.p, sf.pp, log='xy', main = "piRNA library vs sizeFactors ", col='darkgreen', pch=16, cex=1.5,
           xlab= 'piRNA library size', ylab='sizeFactor.piRNAs')
      #abline(0, (median(sf.p/sf.s)), lwd=2.0, col="red")
      text(sf.p, sf.pp, labels = colnames(pans), offset = 0.5, pos = 1, cex = 0.7)
     
      plot(sf.s, sf.pp, log='xy', main = "spikeIns norm vs piRNAs sizeFactors", col='darkgreen', pch=16, cex=1.5, xlim = c(4, 400), ylim = c(0.1, 10),
           xlab= 'sizeFactor.spikeIns', ylab='sizeFactor.piRNAs')
      #abline(0, (median(sf.p/sf.s)), lwd=2.0, col="red")
      text(sf.s, sf.pp, labels = colnames(pans), offset = 0.5, pos = 1, cex = 0.7)
      
      kk = c(grep("71822", colnames(pans)), grep("71823", colnames(pans)))
      
      plot(sf.s[-kk], sf.pp[-kk], log='xy', main = "spikeIns norm vs piRNAs sizeFactors", col='darkgreen', pch=16, cex=1.5, xlim = c(4, 400), ylim = c(0.1, 10),
           xlab= 'sizeFactor.spikeIns', ylab='sizeFactor.piRNAs')
      #abline(0, (median(sf.pp[-kk]/sf.s[-kk])), lwd=2.0, col="red")
      text(sf.s[-kk], sf.pp[-kk], labels = colnames(pans)[-kk], offset = 0.5, pos = 1, cex = 0.7)
      
    }
    
    # compare the cpm, spike-in normalization and piRNA normalization
    par(mfrow=c(1,1))
    source("RNAseq_Quality_Controls.R")
    o1 = order(design.matrix.pans$treatment, design.matrix.pans$promoter);
    plot.pair.comparison.plot(cpm[, o1], main = "rab-3 WT pairwise comparison for cpm ")
    #my.plotPCA(cpm, design.matrix = design.matrix.pans)
    plot.pair.comparison.plot(res[, o1], main = "rab-3 pairwise comparison for spike-in normalization")
    #my.plotPCA(res, design.matrix = design.matrix.pans)
    if(is.null(piRNA.sizeFctors)){
      sizefactors = as.numeric(stats.pans$piRNA)
    }else{
      sizefactors = as.numeric(piRNA.sizeFctors.pans)  
      cat("-- use piRNA size factors --\n")
    }
    
    cpm.piRNA = pans;
    for(n in 1:ncol(cpm.piRNA)) cpm.piRNA[,n] = cpm.piRNA[,n]/sizefactors[n]*10^6
    
    plot.pair.comparison.plot(cpm.piRNA[, o1], main = "rab-3 pairwise comparison for piRNA-normalization")
    
  }
  
  ###############################
  # check the cpms and spike-in normalized data
  ###############################
  cpms = my.cpm.normalization(countData = countData)
  if(is.null(piRNA.sizeFctors)){
    sizefactors = as.numeric(stats$piRNA)
  }else{
    cat("-- using piRNA size factors \n")
    sizefactors = as.numeric(piRNA.sizeFctors)
  }
  
  cpm.piRNA = countData;
  for(n in 1:ncol(cpm.piRNA)) cpm.piRNA[,n] = cpm.piRNA[,n]/sizefactors[n]*10^6
  
  if(cpm.pairwise.compare){
    my.plotPCA(cpms, design.matrix)
  }
  my.plotPCA(cpm.piRNA, design.matrix = design.matrix)
  
  tcs = unique(design.matrix$tissue.cell)
  for(tc in tcs){
    jj = which(design.matrix$tissue.cell == tc)
    jj = jj[order(design.matrix$treatment[jj])]
    if(cpm.pairwise.compare)  plot.pair.comparison.plot(cpms[,jj], main = paste0(tc, "-- pairwise comparison for cpm"))
    plot.pair.comparison.plot(cpm.piRNA[,jj], main = paste0(tc, "-- pairwise comparison for piRNA-normalization"))
  }
  
  ###############################
  # compare the rab-3 untreated with other untreated samples
  ###############################
  Compare.pans.vs.other.untreated = FALSE
  if(Compare.pans.vs.other.untreated){
    kk = grep("untreated_", colnames(cpm.piRNA))
    par(mfrow=c(3,2))
    for(n in kk){
      if(n != 55 & n != 57){
        plot(cpm.piRNA[, c(n, 55)], log='xy'); abline(0, 1, lwd=2.0, col='red')
        plot(cpm.piRNA[, c(n, 57)], log='xy'); abline(0, 1, lwd=2.0, col='red')
      }
      
    }
  }
  
}

my.plotPCA = function(xx, design.matrix){
  library("vsn");
  library("pheatmap");
  library("RColorBrewer");
  library("dplyr");  
  library("ggplot2")
  library(sva)
  library(limma)
  require('RColorBrewer')
  
  # xx = cpms;
  yy = log2(xx + 2^-6)
  
  pcs <- prcomp(yy, center = FALSE, scale = FALSE) 
  scores = as.data.frame(pcs$rotation)
  #conds = design.matrix$conds
  
  ## configure data frame for ggplots
  #pca2save = as.data.frame(plotPCA(vsd, intgroup = c("condition", "batch"), returnData = TRUE))
  pca2save = data.frame(PC1 = scores$PC1, PC2 = scores$PC2)
  pca2save = data.frame(pca2save, tissues = design.matrix$tissue.cell, conds = design.matrix$treatment)
  #pca2save$genotype = design.matrix$condition;
  #pca2save$conds = paste0(pca2save$conds, ".", design.matrix$batch)
  
  ggp = ggplot(data=pca2save, aes(PC1, PC2, shape= tissues, color=conds)) + 
    geom_point(size=2.5) +
    scale_shape_manual(values=seq(20, 1, by = -1))
  plot(ggp);
  
}

Plot.ProprotionMatrix.ExpressionMatrix = function(proportions.sel, expression.sel, fitting.space = "log2", compare.Pan.vs.otherfiveSamples = FALSE)
{
  library("pheatmap")
  library("RColorBrewer")
  
  xx = proportions.sel;
  xx[which(xx>0)] = 1
  
  yy = expression.sel;
  #yy[yy<0] = 0
  if(fitting.space == 'linear') {
    logaxis = 'xy';
    yy = t(log2(t(yy)/yy[which(rownames(yy)=='background'), ]));
    xx = xx[-1, -1];  yy = yy[-1,  ]
  }else{logaxis = ''}
  
  pheatmap(xx, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE,
           cluster_cols=TRUE, 
           color = c("lightgray", "blue"), legend = FALSE)
  
  pheatmap(yy, cluster_rows=FALSE, show_rownames=TRUE, show_colnames = TRUE, 
           cluster_cols=TRUE, 
           color = colorRampPalette(rev(brewer.pal(n = 7, name="RdYlBu")))(100))
  
  if(compare.Pan.vs.otherfiveSamples){
    ## double check the expression matrix
    par(mfrow = c(1, 1))
    mm = match(c("Cholinergic", "Glutamatergic",  "GABAergic",  "Dopaminergic", "Serotonergic"), rownames(yy))
    ranges = range(c(yy[which(rownames(yy)=="Pan.neurons"), ], apply(as.matrix(yy[mm, ]), 2, sum)))
    
    plot(yy[which(rownames(yy)=="Pan.neurons"), ], apply(as.matrix(yy[mm, ]), 2, sum), xlab = "Pan.neurons", 
         ylab = "sum of Cho, Glut, GABA, Dop and Ser", main = "compare pan vs other five samples", xlim = c(-1, 10), ylim = c(-1, 25))
    abline(0, 1, col="red", lwd=2.0)
    abline(0.5, 1, col="red", lwd=1.0, lty = 2)
    abline(-0.5, 1, col="red", lwd=1.0, lty =2)
    abline(1, 1, col="red", lwd=1.0, lty=3)
    abline(-1, 1, col="red", lwd=1.0, lty=3)
    #text(yy[which(rownames(yy)=="Pan"), ], apply(as.matrix(yy[mm, ]), 2, sum), labels = colnames(yy), cex = 0.7,
    #     pos = 1, offset = 0.4)
    #plot(yy[which(rownames(yy)=="Pan.neurons"), ], apply(as.matrix(yy[mm, ]), 2, sum), xlab = "Pan.neurons", 
    #     ylab = "sum of Cho, Glut, GABA, Dop and Ser", xlim = ranges, ylim = ranges)
    #abline(0, 1, col="red", lwd=2.0)
    #abline(h=1, col="darkgray", lwd=2.0)
    
    text(yy[which(rownames(yy)=="Pan.neurons"), ], apply(as.matrix(yy[mm, ]), 2, sum), labels = colnames(yy), cex = 0.8,
         pos = 1, offset = 0.4)
    
    par(mfrow = c(1, 2))
    plot(t(expression.sel[match(c("Dopaminergic", "Ciliatedsensory"), rownames(expression.sel)), ]), log=logaxis)
    abline(0, 1, lwd=2.0, col='red')
    plot(t(expression.sel[match(c("Mechanosensory",  "unc.86"), rownames(expression.sel)), ]), log=logaxis)
    abline(0, 1, lwd=2.0, col='red')
  }
  
}

########################################################
########################################################
# Section : functions of comparing the pan.neuron and other five samples
#  and check miRNA examples
########################################################
########################################################
plot.comparison.pan.vs.other.five.samples = function(cpm, logscale = FALSE)
{
  if(logscale){
    cpm.piRNA.bc.meanrep = average.biological.replicates(log2(cpm)); 
  }else{
    cpm.piRNA.bc.meanrep = average.biological.replicates(cpm);
  }
  
  jj = grep('_untreated', colnames(cpm.piRNA.bc.meanrep))
  total = apply(cpm.piRNA.bc.meanrep[, jj], 1, median)
  xx = data.frame(total, cpm.piRNA.bc.meanrep[, -jj])
  ncs = sapply(colnames(xx)[-c(1:2)], function(x) unlist(strsplit(x, "_"))[2], USE.NAMES = FALSE)
  ncs = sapply(ncs, function(x) gsub("*.neurons", "", x), USE.NAMES = FALSE)
  
  colnames(xx) = c('whole.body', 'background', ncs)
  
  expression = xx[, -c(1:2)]
  
  # N2 as background in both linear and log scale 
  for(n in 1:ncol(expression)) expression[,n] = expression[,n] - xx$background
  
  enriched.list = read.table(file = paste0(resDir, "/tables/Enrichment_Matrix_13samples_66genes_with_clusters_for_neuronClasses.txt"), 
                             sep = "\t", header = TRUE, row.names = 1)
  enriched.list = colnames(enriched.list)
  enriched.list = sapply(enriched.list, function(x) gsub("[.]", "-", x), USE.NAMES = FALSE)
  mm = match((enriched.list), rownames(expression))
  
  expression.sel = t(expression[mm, ])
  #expression.sel = log2(expression.sel)
  
  ##########################################
  # for both log and liear scale, the signals below background are converted to 0 (in log scale) or very samll number (10^-6)
  ##########################################
  yy = expression.sel;
  mm = match(c("cholinergic", "Glutamatergic",  "GABAergic",  "Dopaminergic", "Serotonergic"), rownames(yy))
  
  #ranges = range(c(yy[which(rownames(yy)=="Pan.neurons"), ], apply(as.matrix(yy[mm, ]), 2, sum)))
  if(logscale){
    yy[yy<0] = 0;
    sel = paste0("convert to zeros if below background (log2 scale)")
    cat(sel, "\n")
    xlim = c(-1, 10); ylim = c(-1, 25)
    plot(yy[which(rownames(yy)=="Pan"), ], apply(as.matrix(yy[mm, ]), 2, sum), xlab = "Pan.neurons", 
         ylab = "sum of Cho, Glut, GABA, Dop and Ser", main = sel, xlim = xlim, ylim = ylim)
    abline(0, 1, col="red", lwd=2.0)
    abline(0.5, 1, col="red", lwd=1.0, lty = 2)
    abline(-0.5, 1, col="red", lwd=1.0, lty =2)
    abline(1, 1, col="red", lwd=1.0, lty=3)
    abline(-1, 1, col="red", lwd=1.0, lty=3)
    text(yy[which(rownames(yy)=="Pan"), ], apply(as.matrix(yy[mm, ]), 2, sum), labels = colnames(yy), cex = 0.5,
         pos = 1, offset = 0.4)
  }else{
    #if(convert.negatives == 0){}
    yy[yy<0] = 10^-2;
    sel = paste0("convert 10-2 if below background (linear scale)")
    cat(sel, "\n")
    
    xlim = range(yy[which(rownames(yy)=="Pan"), ]);
    ylim = range(apply(as.matrix(yy[mm, ]), 2, sum))
    plot(yy[which(rownames(yy)=="Pan"), ], apply(as.matrix(yy[mm, ]), 2, sum), xlab = "Pan.neurons", 
         ylab = "sum of Cho, Glut, GABA, Dop and Ser", main = sel, xlim = xlim, ylim = ylim, log = "xy")
    ratios = median(apply(as.matrix(yy[mm, ]), 2, sum) /yy[which(rownames(yy)=="Pan"), ])
    abline(0, 1, col="red", lwd=2.0)
    points(xlim, ratios*xlim, col="darkgreen", lwd=2.0, type = "l")
    abline(0.5, 1, col="red", lwd=1.0, lty = 2)
    abline(-0.5, 1, col="red", lwd=1.0, lty =2)
    abline(1, 1, col="red", lwd=1.0, lty=3)
    abline(-1, 1, col="red", lwd=1.0, lty=3)
    text(yy[which(rownames(yy)=="Pan"), ], apply(as.matrix(yy[mm, ]), 2, sum), labels = colnames(yy), cex = 0.5,
         pos = 1, offset = 0.4)
  }
  
  #for(convert.negatives in c(0, 1)){ }
}

plot.miRNA.examples = function(cpm,  gg = c("lsy-6"), design.matrix, logscale = FALSE)
{
  ns = unique(design.matrix$tissue.cell)
  
  if(logscale){
    cpm = log2(cpm)
  }
  
  kk = which(rownames(cpm)==gg)
  index.bg = which(design.matrix$tissue.cell == "whole.body" & design.matrix$treatment=="treated");
  bg.mean = median(cpm[kk, index.bg]);
  ns = setdiff(ns, "whole.body")
  ns.mean = c()
  
  if(logscale){
    lims = range(cpm[kk, ] - bg.mean)
    main = paste0(gg, "-- log2 scale");
  }else{
    lims = range(cpm[kk, grep("_treated", colnames(cpm))] - bg.mean)
    main = paste0(gg, "-- linear scale");
  }
  
  plot(c(1:length(ns)), rep(1, length(ns)), type= 'n', col='darkblue', cex=1.0, log='', ylim =lims, main = main, xlab=NA, 
       ylab = 'normalizaed by piRNAs', axe = FALSE)
  
  for(n in 1:length(ns))
  {
    index.ns = which(design.matrix$tissue.cell==ns[n] & design.matrix$treatment=="treated")
    points(rep(n, length(index.ns)), (cpm[kk, index.ns] - bg.mean), type = "p", col='darkblue', cex=1., pch =16)
    ns.mean = c(ns.mean, median(cpm[kk, index.ns] - bg.mean))
    # add sample ids 
    if(ns[n]=="Pan.neurons"){
      text(rep(n, length(index.ns)), (cpm[kk, index.ns] - bg.mean), design.matrix$SampleID[index.ns], pos = 2, offset = 0.5, cex = 0.8)
      abline(h = median((cpm[kk, index.ns] - bg.mean)), lwd=1.5, col = 'darkblue')
    }
    index.ns = which(design.matrix$tissue.cell==ns[n] & design.matrix$treatment=="untreated")
    points(rep(n, length(index.ns)), cpm[kk, index.ns], type = "p", col='black', cex=1., pch = 0)
    
  }
  points(c(1:length(ns)), ns.mean, type = "l", cex=1.0, col = 'darkgreen')
  points(c(1:length(ns)), ns.mean, type = "p", cex=1.0, col = 'darkgreen')
  axis(2, las= 1)
  ns.short = sapply(ns, function(x) unlist(strsplit(x, "[.]"))[1], USE.NAMES = FALSE)
  axis(1, at=c(1:length(ns)), labels = ns.short, las=2,cex=0.5)
  box()
  #abline(h=c(0, 2, 5), lwd=0.7, col='red')
  abline(h=0, lwd=2.0, col='darkgray')
  if(logscale)  abline(h=c(-1, 1), lwd=1.5, col='darkgray', lty=1)
  
}

Compare.pan.neuron.vs.other.five.samples.And.check.miRNA.examples = function(cpm.piRNA.bc, design.matrix)
{
  ##########################################
  # compare the pan.neuron vs the sum of other five samples
  # which is a control showing that our data is good for the linear model
  ##########################################
  selections = c("all")
  for(sel in selections)
  {
    par(mfrow = c(1, 2))
    plot.comparison.pan.vs.other.five.samples(cpm.piRNA.bc, logscale = TRUE)
    plot.comparison.pan.vs.other.five.samples(cpm.piRNA.bc)
  }
  
  ##########################################
  # check lsy-6 and other examples
  ##########################################
  examples = c("lsy-6", "mir-791", "mir-793",  "mir-792",
               "mir-1821", "mir-83", "mir-124", "mir-1", 
               "mir-243", "lin-4", "mir-249", "mir-789-2",
               "mir-1830", "mir-392", "mir-238", "mir-254",
               "mir-245", "mir-34", "mir-794", "mir-1020")

  par(mfrow=c(1, 2))
  par(cex =0.7, mar = c(8,3,2,0.8)+0.1, mgp = c(1.6,0.5,0),las = 0, tcl = -0.3)
  for(gg in examples)
  {
    # gg = "lsy-6"
    plot.miRNA.examples(cpm = cpm.piRNA.bc, gg, design.matrix, logscale = TRUE)
    plot.miRNA.examples(cpm = cpm.piRNA.bc, gg, design.matrix, logscale = FALSE)
  }
}

calibrate.promoter.methylation.efficiency = function(cpm, design.matrix, upper = 2.0, lower.quantile = 0.15)
{
  # cpm = cpm.piRNA.bc 
  cpm = log2(cpm);
  
  enriched =  read.table(file = paste0("../results/tables_for_decomvolution", 
                                       "/tables/Enrichment_Matrix_13samples_allgenes_with_clusters_for_neuronClasses_20181203.txt"), 
                         sep = "\t", header = TRUE, row.names = 1)
  enriched.list = read.table(file = paste0("../results/tables_for_decomvolution",
                                           "/tables/Enrichment_Matrix_13samples_66genes_with_clusters_for_neuronClasses.txt"), 
                             sep = "\t", header = TRUE, row.names = 1)
  enriched.list = colnames(enriched.list)
  enriched.list = sapply(enriched.list, function(x) gsub("[.]", "-", x), USE.NAMES = FALSE)
  
  ##########################################
  # start the background calibration for the promoter-specific backgrounds 
  ##########################################
  jj = which(design.matrix$treatment == "treated")
  treated = cpm[, jj];
  untreated = cpm[, -jj];
  #treated.mean = average.biological.replicates.for.promoters(treated)
  design.treated = design.matrix[jj, ]
  
  prots = paste0(design.treated$tissue.cell, "_", design.treated$promoter)
  prots.uniq = unique(prots)
  
  design.treated$prots = prots
  yy = c()
  for(n in 1:length(prots.uniq)) {
    yy = cbind(yy, apply(treated[, which(prots == prots.uniq[n])], 1, median))
  }
  colnames(yy) = prots.uniq
  
  ##########################################
  # here a careful selection of non-enriched miRNAs were selected with the aim of having as many as possible
  # and also outlier should be removed for the sake of robustness
  ##########################################
  require(MASS)
  intercepts = c(0, rep(NA, (ncol(yy)-1)))
  
  jj = which(is.na(match(rownames(yy), enriched.list)))
  #jj = c(1:nrow(yy))
  #jj = which(apply(yy, 1, mean)<5 & apply(yy, 1, mean)>-2)
  #par(mfrow=c(3, 5))
  comps = c()
  for(n in 2:ncol(yy)){
    neurons = unlist(strsplit(as.character(colnames(yy)[n]), "_"))[1]
    rrs = yy[jj,n] - yy[jj,1]
    #upper = quantile(rrs, 0.60, names = FALSE); lower = quantile(rrs, 0.10, names = FALSE)
    #upper = 10; 
    lower = quantile(rrs, lower.quantile, names = FALSE)
    #upper = max(rrs); lower = min(rrs);
    
    jj.middle = jj[intersect(which(rrs >= lower), which(rrs <= upper))]
    
    cat("nb of non-enriched miRNAs used to calibrate the promoter efficiency -- ", length(jj.middle), "\n")
    
    rfit = rlm((yy[jj.middle, n] - yy[jj.middle, 1]) ~ 1 )
    intercepts[n] = rfit$coefficients
    
    comps = rbind(comps, cbind(rep(n, length(jj.middle)), 
                               (yy[jj.middle, n] - yy[jj.middle, 1]), 
                               (yy[jj.middle, n] - yy[jj.middle, 1]) - intercepts[n]))
    
    cat(median(rrs[which(rrs>lower & rrs< upper)]), " -- ",  intercepts[n], "\n")
        
    kk = which(design.treated$prots == colnames(yy)[n])
    for(k in kk) treated[, k] = treated[,k] - intercepts[n];
    
  }
  names(intercepts) = colnames(yy)
  
  load(file = paste0("../results/tables_for_decomvolution/Rdata/",
                     "Tables_Sample_2_Promoters_mapping_neurons_vs_neuronClasses_FractionMatrix_plus_mergedFractionMatrix", 
                     "_miRNAs_neurons_20180525", ".Rdata"))
  nbcells = apply(as.matrix(newprop), 1, sum)
  nbcells = c(2, 6, 6, 70, 56, 15, 6, 30, 20, 107, 16, 28, 32, 220, 220)
 
  intercepts = intercepts[-1]
  # cat(length(nbcells), "--", length(intercepts), "\n")
  
  par(mfrow=c(1, 1))
  plot(nbcells, intercepts, log='x', xlim = c(1, 500), ylim = c(-1, 1), 
       xlab = "nb.cells", ylab = "estimated bias due to promoter efficiency")
  text(nbcells, intercepts, names(intercepts), pos=1, cex = 0.8, offset = 0.2)
  abline(h = c(-1, 0, 1), col = 'red', lwd=2.0)
  
  par(mfrow = c(1, 2))
  boxplot(comps[, 2] ~ comps[, 1], names = names(intercepts), las = 2, col = c(1:length(intercepts)))
  abline(h=0, col='darkgray', lwd=2.0)
  boxplot(comps[, 3] ~ comps[, 1], names = names(intercepts), las = 2, col = c(1:length(intercepts)))
  abline(h=0, col='darkgray', lwd=2.0)
  #cpm.piRNA.bc.bgc = data.frame(untreated, treated, stringsAsFactors = FALSE)
  cpm.piRNA.bc.bgc = cbind(untreated, treated)
  cpm.piRNA.bc.bgc = cpm.piRNA.bc.bgc[, match(colnames(cpm), colnames(cpm.piRNA.bc.bgc))]
  
  return(cpm.piRNA.bc.bgc)
  
}

find.non.enriched.miRNAs = function(neurons, enrich.matrix, fc.cutoff = 0, pval.cutoff = 0.05)
{
  ggs = c()
  sels = grep(neurons, colnames(enrich.matrix))
  if(length(sels)==0) sels = grep(gsub("-", ".", neurons), colnames(enrich.matrix))
  for(n in 1:nrow(enrich.matrix)){
    #deleted = rep(FALSE, length(sels))
    index.pval = sels[grep("_pvalue", colnames(enrich.matrix)[sels])]
    index.fc = sels[grep("_log2FC", colnames(enrich.matrix)[sels])]
    if(all(!is.na(enrich.matrix[n, sels]))){
      if(all(enrich.matrix[n, index.fc]<fc.cutoff) && enrich.matrix[n, index.pval]> pval.cutoff) ggs = c(ggs, rownames(enrich.matrix)[n])
    }
    #for(m in 1:length(sels))
    #{
    #  if(grep("_log2FC", colnames(yy)[sels[m]])){
    #    if(enrich.matrix[n, m]<0) deleted
    #  } 
    #}
  }
  
  return(ggs)
}

merge.interactionDiff = function(res1, res2){
  ## merge both res1 and res2
  res1 = as.data.frame(res1); 
  res2 = as.data.frame(res2);
  
  ## change the sign (log2(WT.treated) - log2(WT.untreated)) - (log2(N2.treated) - log2(N2.untreated))
  # res2$log2FoldChange = -res2$log2FoldChange; 
  ii.enrich = which(res1$log2FoldChange>0) 
  res1[ii.enrich, ] = res2[ii.enrich, ] ## replaced by res2 if enriched; keep if depleted
  
  ## merge the table with comparison with N2 and the one without
  #colnames(res1) = paste0(colnames(res1), ".without.N2")
  #colnames(res2) = paste0(colnames(res2), ".with.N2")
  res = data.frame(res1, stringsAsFactors = FALSE)
  colnames(res) = paste0(colnames(res), ".with.N2")
  
  return(res)
  
}

enrichmentAnalysis.for.mutant.Norm.piRNA = function(countData, design.matrix = design, sfs=NULL)
{
  require(DESeq2)
  # countData = raw[, sels]; design.matrix = design[sels, ]; sfs = sizefactors[sels];
  dds <- DESeqDataSetFromMatrix(countData, DataFrame(design.matrix), design = ~ genotype + treatment + genotype:treatment)
  sizeFactors(dds) = sfs/median(sfs)
  dds$genotype <- relevel(dds$genotype, ref="N2") ## reference is N2 and untreated
  dds$treatment = relevel(dds$treatment, ref="untreated")
  
  
  dds = estimateDispersions(dds, fitType = "parametric")
  par(cex = 1.0, las = 1, mgp = c(2,0.2,0), mar = c(3,2,2,0.2), tcl = -0.3)
  plotDispEsts(dds, ylim=c(0.001, 10), cex=1.0)
  
  cpm = fpm(dds, robust = TRUE)
  
  dds = nbinomWaldTest(dds, betaPrior = FALSE)
  resultsNames(dds)
  
  ##########################################
  # find enriched miRNAs for WT: 
  # Here we use N2 and untreated condition as reference 
  ##########################################
  ## the condition effect for genotype WT
  #res1 = results(dds, contrast=list( c("treatment_treated_vs_untreated", "genotypeWT.treatmenttreated")))
  ## test if condition effect in WT is greater than N2
  ## i.e. test the significance log2(WT.treated/WT.untreated) - log2(N2.treated/N2.untreated)
  ## https://support.bioconductor.org/p/60543/ (one of good explanation for this)
  keep1 = results(dds, name = "genotypeWT.treatmenttreated", lfcThreshold = 0, altHypothesis = "greater" )
  
  #keep1 = merge.interactionDiff(res1, res2)
  colnames(keep1) = paste0("WT.", colnames(keep1))
  #res1 = results(dds, contrast=c("treatment","treated","untreated"))
  #summary(res)
  #plot(res$log2FoldChange, -log10(as.numeric(res$pvalue)), xlab='log2(FoldChange)', ylab='-log10(pvalue)', cex=0.8, 
  #     main=paste0("WT --", " (NOT using N2)"))
  #abline(v=0, lwd=2.0, col='black')
  #abline(h=c(5, 10), lwd=1.0, lty=2, col='blue')
  #text(res$log2FoldChange, -log10(res$pvalue), rownames(res), cex=0.7, offset = 0.3, pos = 3)
  
  #res1 = res
  #res2 = results(dds, name = "genotypeN2.treatmenttreated", lfcThreshold = 0, altHypothesis = "less")
 
  ##########################################
  # find enriched miRNAs for WT_tax4_mutant_ks28
  ##########################################
  #res1 = results(dds, contrast=list( c("treatment_treated_vs_untreated", "genotypeWT_tax4_mutant_ks28.treatmenttreated")))
  keep2 = results(dds, name = "genotypeWT_tax4_mutant_ks28.treatmenttreated", lfcThreshold = 0, altHypothesis = "greater" )
  #keep2 = merge.interactionDiff(res1, res2)
  colnames(keep2) = paste0("tax4_mutant_ks28.", colnames(keep2))
  ##########################################
  # find enriched miRNAs for WT_tax4_mutant_p678
  ##########################################
  #res1 = results(dds, contrast=list( c("treatment_treated_vs_untreated", "genotypeWT_tax4_mutant_p678.treatmenttreated")))
  keep3 = results(dds, name = "genotypeWT_tax4_mutant_p678.treatmenttreated", lfcThreshold = 0, altHypothesis = "greater" )
  #keep3 = merge.interactionDiff(res1, res2)
  colnames(keep3) = paste0("tax4_mutant_p678.", colnames(keep3))
  
  ##########################################
  # test if condition effect is different in WT_tax4_mutant_ks28 compare to WT
  # and if condition effect is different in WT_tax4_mutant_p678 compare to WT
  ##########################################
  keep4 = as.data.frame(results(dds, contrast = list("genotypeWT_tax4_mutant_ks28.treatmenttreated", "genotypeWT.treatmenttreated")))
  colnames(keep4) = paste0("tax4_mutant_ks28_vs_WT.", colnames(keep4))
  keep5 = as.data.frame(results(dds, contrast = list("genotypeWT_tax4_mutant_p678.treatmenttreated", "genotypeWT.treatmenttreated")))
  colnames(keep5) = paste0("tax4_mutant_p678_vs_WT.", colnames(keep5))
  
  keeps = data.frame(cpm, keep1[,c(2,5)], keep2[, c(2, 5)], keep3[, c(2,5)], keep4[, c(2,5)], keep5[, c(2,5)], stringsAsFactors = FALSE)
  
  keeps = keeps[order(-keeps$tax4_mutant_ks28_vs_WT.log2FoldChange), ]
  
  return(keeps)
  
}

######################################
######################################
## Section (additional) : 
# other utility plots or tables
# calculate pvalue for overlapping groups
######################################
######################################
calculate.pval.for.overlapping.groups = function(xx)
{
  if(calculate.pval.for.overlapping.groups){
    
    library(openxlsx)
    aa = read.xlsx("../results/decomvolution_results/pvalue_overlapping_groups.xlsx", sheet = 2, colNames = TRUE, detectDates = TRUE)
    aa = aa[which(!is.na(aa[, 1])==TRUE), ] 
    aa = data.frame(aa)
    aa$TOTAL.NUMBER.OF.EXPRESSED.miRNAs = 123
    
    aa = aa[, c(1:5)]
    source("miRNAseq_functions.R")
    xx = c()  
    for(n in 1:nrow(aa))
    {
      #n = 1
      xx = rbind(xx, calculate.pvalues.two.groups.overlapping(aa$TOTAL.NUMBER.OF.EXPRESSED.miRNAs[n], 
                                                              aa$groupA[n],
                                                              aa$groupB[n], 
                                                              aa$OVERLAP[n]))
    }
    
    aa = data.frame(aa, observersion.expected=xx[,1], pval=xx[, 2], stringsAsFactors = FALSE)
    write.xlsx(aa, file='..//results/decomvolution_results/pvalue_overlapping_groups_byJingkui.xlsx')
    
  }
}


