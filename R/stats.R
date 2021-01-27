###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch


## ROC curves
.ezROC = function(truth, data, cuts=quantile(data, probs=(1:1000)/1000), biggerIsBetter=FALSE){
  
  x1 = data[truth == 1]
  x0 = data[truth == 0]
  if (biggerIsBetter){
    cuts = sort(unique(cuts), decreasing=TRUE) ## cuts are sorted by decreasing stringency
    tpr = unlist(lapply(cuts, function(cp, val){mean(val > cp)}, val=x1))
    fpr = unlist(lapply(cuts, function(cp, val){mean(val > cp)}, val=x0))
  } else {
    cuts = sort(unique(cuts), decreasing=FALSE) ## cuts are sorted by decreasing stringency
    tpr = unlist(lapply(cuts, function(cp, val){mean(val < cp)}, val=x1))
    fpr = unlist(lapply(cuts, function(cp, val){mean(val < cp)}, val=x0))   
  }
  return(list(sens=tpr, spec=1-fpr, cuts=cuts))
}

## AUC stats
## see also caTools::trapz
## or google "R auc"
.ezAUC = function(rocList){
  x <- c(1 - rocList$spec, 1)
  y <- c(rocList$sens, rocList$sens[length(rocList$sens)])
  abs(as.double(diff(x) %*% ((y[-1] + y[-length(y)])/2))[])
}


.ezROCtest = function(){
  n = 100
  negVal = rnorm(n=n, mean=1, sd=1)
  posVal = rnorm(n=n, mean=2, sd=1)
  truth = c(rep(1, n), rep(0, n))
  data = c(negVal, posVal)
  idx = sample(x=1:(2*n), size=2*n, replace=FALSE)
  cuts = round(c(-1000:1000)/ 10) /10
  cuts = sample(cuts, size=length(cuts), replace=FALSE)
  roc = .ezROC(truth, data[idx], cuts=cuts, biggerIsBetter=FALSE)
  plot(1 - roc$spec, roc$sens, type="l")
}
