############# 
# cat tables 
#############
cat.countTable = function(xlist, countsfrom = 'featureCounts')
{
  ## input is list.files for count tables (including directories and file names)
  if(countsfrom == "featureCounts"){
    counts = NULL
    for(n in 1:length(xlist)){
      # n = 1
      cat(xlist[n], '--\n')
      
      ff = read.delim(xlist[n], sep='\t', header = TRUE, as.is = c(1), comment.char = "#");
      ff = ff[, c(1, ncol(ff))]
      
      if(n==1){
        ggs = unique(ff[, 1]);
        counts = data.frame(ggs, ff[match(ggs, ff[, 1]) , -1], stringsAsFactors = FALSE);
      }else{
        ggs = unique(c(counts[, 1], ff[, 1]));
        counts = data.frame(ggs, counts[match(ggs, counts[, 1]), -1], ff[match(ggs, ff[, 1]) , -1], stringsAsFactors = FALSE);
      }
    }
    
  }else{
    counts = NULL
    if(countsfrom == "htseq"){
      for(n in 1:length(xlist))
      {
        # n = 1
        cat(n, '\t')
        cat(xlist[n], '\n')
        
        if(n==1){
          counts = read.table(xlist[n], header = FALSE);
          colnames(counts) = c('RefseqID', basename(xlist[n]))
          counts = data.frame(counts, stringsAsFactors = FALSE)
        }else{
          test = read.table(xlist[n], header = FALSE);
          colnames(test) = c('RefseqID',  basename(xlist[n]))
          mm = match(counts$RefseqID, test$RefseqID)
          counts = data.frame(counts, test[mm, -grep('RefseqID', colnames(test))])
        }
      }
    }
    
    if(countsfrom == 'slamdunk'){
      for(n in 1:length(xlist))
      {
        # n = 1
        cat(n, '\t')
        cat(xlist[n], '\n')
        
        if(n==1){
          counts = read.table(xlist[n], header = TRUE);
          counts = data.frame(counts$Name, counts$ReadCount, stringsAsFactors = FALSE)
          colnames(counts) = c('RefseqID', basename(xlist[n]))
          #counts = data.frame(counts, stringsAsFactors = FALSE)
        }else{
          test = read.table(xlist[n], header = TRUE);
          test = data.frame(test$Name, test$ReadCount, stringsAsFactors = FALSE)
          colnames(test) = c('RefseqID',  basename(xlist[n]))
          #mm = match(counts$RefseqID, test$RefseqID)
          counts = data.frame(counts, test[, 2])
        }
      }
    }
  }
  
  colnames(counts)[1] = 'gene'
  colnames(counts)[-1] = basename(xlist)
  return(counts)
  
}

merge.countTables.htseq = function(file_list){
  
  tables <- lapply(file_list, function(f) {
    tmp <- read.delim(f)
    colnames(tmp) = c('gene', basename(f))
    return(tmp)
  })
  xx <- Reduce(function(dtf1, dtf2) merge(dtf1, dtf2, by = "gene", all = TRUE), tables)
  return(xx)
}


process.countTable = function(all, design, special.column = NULL, ensToGeneSymbol = FALSE, mergeCounts.forSameGeneSysbols = TRUE)
{
  index = c()
  index.design.with.found.sample = c()
  for(n in 1:nrow(design))
  {
    #n = 1;
    jj = grep(design$SampleID[n], colnames(all))
    if(!is.null(special.column))  jj = intersect(jj, grep(special.column, colnames(all)));
   
    if(length(jj)==1) {
      index = c(index,jj)
      index.design.with.found.sample = c(index.design.with.found.sample, n)
    }else{
      print(paste0("ERROR for sample--", design$SampleID[n]))
    }
  }
  
  newall = data.frame(as.character(all[, 1]),  as.matrix(all[, index]), stringsAsFactors = FALSE)
  colnames(newall)[1] = "gene"
  kk = which(colnames(design)=='SampleID')
  kk = c(setdiff(c(1:ncol(design)), kk), kk)
  colnames(newall)[-1] = apply(design[, kk], 1, function(x)  paste(x, collapse = "_"))[index.design.with.found.sample]
  colnames(newall) = sapply(colnames(newall), function(x) gsub(' ', '', x))
  #paste(design$gen, "_", design$tissue.cell, "_", design$SampleID)
  
  if(ensToGeneSymbol){
    annot = read.csv(file = "/Volumes/groups/cochella/jiwang//annotations/BioMart_WBcel235_noFilters.csv", header = TRUE, stringsAsFactors = FALSE)
    #annot = data.frame(annot, stringsAsFactors = FALSE)
    #######
    # the following is only for protein-coding genes
    #######
    # options(warn=2)
    # result <- try(load(file='/Volumes/groups/cochella/jiwang/annotations/BioMart_WBcel235.Rdata'), silent=TRUE)
    # if (class(result) == "try-error"){
    #   load(file='/Volumes/groups-1/cochella/jiwang/annotations/BioMart_WBcel235.Rdata')
    # }else{
    #   load(file='/Volumes/groups/cochella/jiwang/annotations/BioMart_WBcel235.Rdata')
    # }
    # options(warn = -1)
    
    mm = match(newall$gene, as.character(annot$Gene.stable.ID))
    jj = which(!is.na(mm))
    newall$gene[jj] = as.character(annot$Gene.name[mm[jj]])
    
    ss = apply(as.matrix(newall[, -1]), 1, sum)
    newall = newall[which(ss>0), ]
    
  }
  
  return(newall)
}

compare.readCount.UMI =function(design, aa, normalized = FALSE){
  if(!normalized){
    for(n in 1:nrow(design)){
      id = design$SampleID[n]
      #par(mfrow=c(1,2))
      plot(aa[, intersect(grep(id, colnames(aa)), grep(".readCount", colnames(aa)))], 
           aa[, intersect(grep(id, colnames(aa)), grep(".UMI", colnames(aa)))], 
           log='xy', main= paste0(design[n, ], collapse = "_"), xlab = 'read counts', ylab =' umi.counts',
           cex = 0.4
      )
      #points(spikes[, c(grep(paste0("Total.spikeIn.", id), colnames(spikes)), 
      #                  grep(paste0("Total.UMI.spikeIn.", id), colnames(spikes)))], col = 'darkblue', cex = 1., pch=16)
      abline(0, 1, lwd=1.2, col = 'red')
      
      #plot(aa[, intersect(grep(id, colnames(aa)), grep("Total.count", colnames(aa)))], 
      #     aa[, intersect(grep(id, colnames(aa)), grep("Total.UMInum.count", colnames(aa)))], 
      #     log='xy', main= paste0(design[n, ], collapse = "_"), xlab = 'read counts', ylab =' umi.counts',
      #     cex = 0.4
      #)
      # points(spikes[, c(grep(paste0("Total.spikeIn.", id), colnames(spikes)), 
      #                   grep(paste0("Total.UMI.spikeIn.", id), colnames(spikes)))], col = 'darkblue', cex = 1., pch=16)
      # abline(0, 1, lwd=1.2, col = 'red')
    }
  }else{
    xx = process.countTable(all=aa, design = design, special.column = ".UMI")
    gene.filtered = which(!is.na(xx$gene))
    xx = xx[gene.filtered, ]
    raw = ceiling(as.matrix(xx[, -1]))
    raw[which(is.na(raw))] = 0
    rownames(raw) = xx$gene
    dds = DESeqDataSetFromMatrix(raw, DataFrame(factor(rep('A', ncol(raw)))), design = ~ 1)
    gene.filtered2 = rowSums(counts(dds)) >= 10;
    
    dds <- dds[gene.filtered2, ]
    dds <- estimateSizeFactors(dds)
    umi = fpm(dds, robust = TRUE)
    
    yy = process.countTable(all=aa, design = design, special.column = ".readCount")
    yy = yy[gene.filtered, ]
    raw = ceiling(as.matrix(yy[, -1]))
    raw[which(is.na(raw))] = 0
    rownames(raw) = yy$gene
    dds = DESeqDataSetFromMatrix(raw, DataFrame(factor(rep('A', ncol(raw)))), design = ~ 1)
    dds = dds[match(rownames(umi), rownames(dds)), ]

    dds <- estimateSizeFactors(dds)
    reads = fpm(dds, robust = TRUE)
    
    conds = paste0(design$stage, "_", design$condition)
    
    for(cc in unique(conds)){
      cat(cc, "\n")
      jj = which(conds == cc)
      plot.pair.comparison.plot(umi[, jj], main = paste0(cc, " - umi"))
      plot.pair.comparison.plot(reads[, jj], main = paste0(cc, " - read"))
      
      vars.umi = apply(log2(umi[, jj]), 1, var)
      vars.read = apply(log2(reads[, jj]), 1, var)
      
      plot(vars.umi, vars.read, log='xy', cex= 0.4);
      abline(0, 1, lwd=1.0, col='red')
      
      # plot(apply(log2(umi[, jj]), 1, mean), vars.umi/vars.read, log='y')
    }
  }
}

Merge.techinical.replicates.using.design.countTable = function(design, all, rep.technical=list(c("57751", "57753"), c("57752", "57754")))
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



