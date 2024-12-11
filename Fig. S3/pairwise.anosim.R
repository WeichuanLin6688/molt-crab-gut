pairwise.anosim <-function(x,factors,sim.method,permutations = 999, p.adjust.m)
{
  library(vegan)
  co = as.matrix(combn(unique(factors),2))
  pairs = c()
  R = c()
  p.value = c()
  for(elem in 1:ncol(co)){
    an = anosim(x[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem])),],
                
                factors[factors %in%c(as.character(co[1,elem]),as.character(co[2,elem]))],distance =sim.method,permutations = 999);
    
    pairs =c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    R = c(R,an$statistic);
    p.value = c(p.value,an$signif)
  }
  p.adjusted =p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,R,p.value,p.adjusted)
  return(pairw.res)
}