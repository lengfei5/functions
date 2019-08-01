linear.fitting = function(y, scale='log', tt=c(rep(0, 3), rep(1, 3), rep(2, 3)))
{
  if(scale == 'linear')  y = log2(y);
  fit =lm(y ~ tt)
  slope=summary(fit)$coefficients[2, 1]; 
  intercept = summary(fit)$coefficients[1, 1];
  pval=summary(fit)$coefficients[2, 4]
  return(c(intercept, slope, pval))
}


f2min.sigmoid = function(tt, x, params=c(0.1, 0.1, 0.5, 0.5))
{
  s0 = params[1];a=params[2];
  t0 = params[3];b=params[4];
  error = sum((x-(s0+a/(1+exp(-(tt-t0)*b))))^2)
  #print(params)
  #print(error)
  return(error)
}

which.min = function(y)
{
  y = as.numeric(y);
  return(which(y==min(y)))
}

fit.all.models = function(x, tt=c(c(rep(0, 3), rep(1, 3), rep(2, 3))))
{
  #x = (yy[105, ]);tt=c(c(rep(0, 3), rep(1, 3), rep(2, 3)));
  # x = exprs[c(2), grep("B10", colnames(exprs))]; tt=rep(c(0, 6, 24, 48), each=3)
  x = as.numeric(x);
  #cat(tt, '\n');
  #cat(x, "\n")
  if(length(unique(x))==1){
    cat('EXACT same values for all time points...', x,  '\n')
    res.MS = c(0, rep(NA, 8), 1, 0, 0,  1, 1)
  }else{
    ## M1 constant (no changes)
    err1 = sum((x-mean(x))^2) 
    
    ## M2 one line
    fit =lm(x ~ tt)
    slope=summary(fit)$coefficients[2, 1]; 
    intercept = summary(fit)$coefficients[1, 1];
    err2 = sum(fit$residuals^2)
    
    ## M3 sigmoid curve
    par.init = c(min(x), (max(x)-min(x)+0.05), 12, slope);
    if(slope>0) {
      lower = c(-6, 0, 3, (slope/10));
      upper = c(15, (max(x)-min(x))*2, 36, slope*10);
    }else{
      lower = c(-6, 0, 3, slope*10);
      upper = c(15, (max(x)-min(x))*2, 36, slope/10);
    }
    
    opt = optim(par.init, f2min.sigmoid, x=x, tt=tt, method = 'L-BFGS-B', lower = lower , upper = upper, control = list(maxit=200))
    res.fit = opt$par
    err3 = opt$value
    
    ## BIC scores
    N = length(x)
    score.m1 = N*log(err1/N) + 1*log(N) 
    score.m2 = N*log(err2/N) + 2*log(N) 
    score.m3 = N*log(err3/N) + 4*log(N) 
    
    scores = c(score.m1, score.m2, score.m3)
    scores.relavtive = scores-min(scores)
    prob.model = exp(-0.5*scores.relavtive)
    prob.model = prob.model/sum(prob.model)
    
    res.MS = c(err1, err2, intercept, slope, err3, res.fit, prob.model, which.min(scores), prob.model[which.min(scores)])
  }
  
  res.MS = c(mean(x), res.MS)
  names(res.MS) = c('mean.m1','err.m1', 'err.m2', 'intercept.m2', 'slope.m2',
                    'err.m3',  's0.m3', 'a.m3', 't0.m3', 'b.m3',
                    'prob.m1', 'prob.m2', 'prob.m3', 'best.model', 'prob.best.model')
  
  ## test michelis-menten curve
  TEST.sigmoid.function = FALSE
  if(TEST.sigmoid.function)
  {
    tt = seq(0, 20, length.out = 200)
    s0 = 0.1;a = 20;t0 = 5;
    b = -10;
    plot(tt, (s0+a/(1+exp(-(tt-t0)*b))), type='l', col='green', lwd=2.0)
    abline(h=s0, col='red')
    abline(h=(s0+a), col='blue')
    abline((-b*t0+s0+a/2), b, col='orange', lwd=2.0)
  }
  
  return(res.MS)
  
}

plot.sigmoid.fitting = function(ts, mselect, tt=c(1:10), index=NULL, ylim=NULL)
{
  #ts = exprs[c(1:10), grep("B10", colnames(exprs))]
  mselect = data.frame(mselect)
  if(is.null(index)) index = c(1:nrow(mselect));
  #index = 2
  for(ii in index)
  {
    t = seq(range(tt)[1], range(tt)[2], length.out = 100)
    if(mselect$best.model[ii]==3){
      s0 = mselect$s0.m3[ii]
      a = mselect$a.m3[ii];
      t0 = mselect$t0.m3[ii];
      b = mselect$b.m3[ii];
      yy = (s0+a/(1+exp(-(t-t0)*b)))
      mains = paste0("Sigmoid fitting : s0 - ", signif(s0, d=2), "; a - ", signif(a, d=2), "; t0 - ", signif(t0, d=2), "; b(slop) - ", signif(b, d=2))
    }
    if(mselect$best.model[ii]==1){
      mains = paste0("flat line fitting  ")
      yy = rep(mean(ts[ii, ]), length(t));
    } 
    if(mselect$best.model[ii]==2) {
      mains = paste0("line fitting: intercept - ", signif(mselect$intercept.m2[ii], d=2), "; slope - ", signif(mselect$slope.m2[ii], d=2))
      yy = mselect$intercept.m2[ii] + t*mselect$slope.m2[ii];
    }
    
    mains = paste0(rownames(ts)[ii], "--- M", mselect$best.model[ii], "\n", mains)
    if(is.null(ylim)) {
      lims = range(c(ts[ii,], yy), na.rm = TRUE);
    }else{
      lims = range(c(ts[ii,], yy), na.rm = TRUE);
      lims = c(min(lims[1], ylim[1]), max(lims[2], ylim[2]))
    }
    plot(tt, ts[ii,], type='p', col='blue', xlab='time', ylab='time serie signals', ylim=lims, 
         main= mains)
    points(t, yy, type='l', col='darkgreen', lwd=2.0)
    #abline(h=s0, col='black', lty=2)
    #abline(h=(s0+a), col='black', lty=2)
    #abline((-b*t0+s0+a/2), b, col='orange', lwd=2.0)
    
    #s0 = params[1];a=params[2];
    #t0 = params[3];b=params[4];
    #error = sum((x-(s0+a/(1+exp(-(tt-t0)/b))))^2)
  }
  
}



