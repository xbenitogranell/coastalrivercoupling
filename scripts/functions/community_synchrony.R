#https://zenodo.org/record/818096#.YjnwterMKUk

community.sync.window <- function (d, twin=5, nrands=999) {
  years = as.numeric(rownames(d))
  starts=seq(from=1, to=length(years)-twin+1, by=1)
  ends=starts+twin-1
  mid=starts+(twin-1)/2
  
  results=matrix(nrow=length(starts), ncol=6, NA)
  for (i in 1:length(starts)) {
    locs=which(rownames(d) %in% years[starts[i]:ends[i]])
    c = community.sync(d[locs,], nrands=nrands)
    results[i,] = c(years[starts[i]], years[ends[i]],years[mid[i]], c$meancorr, c$obs, c$pval)
  }
  colnames(results)=c("start", "end","mid", "meancorr", "sync", "pval")
  results=as.data.frame(results)
  return(results)
}
