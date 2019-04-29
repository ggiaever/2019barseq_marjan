############ reformats rownames of counts data matrices; i don't know how they got messed up in the first plate!

mybarseqrows = function(df){
xcnt = as.matrix(df[,2:ncol(df)])
rownames(xcnt) = df$tag
up = grep('uptag',rownames(xcnt))
rownames(xcnt) = gsub('uptag:','',rownames(xcnt))
rownames(xcnt) = gsub('downtag:','',rownames(xcnt))
rownames(xcnt)[up] = paste0(rownames(xcnt)[up],":uptag")
rownames(xcnt)[-up] = paste0(rownames(xcnt)[-up],":downtag")
xcnt
}

###### summing up and down tags for raw count matrix  -- modified --
###### such that genes that have more than 2 tags (cause they were constructed more than once) are not overly represented; also serves as a nice crosscheck; ultimately it would be nice to add the known dubious overlaps for all of the YXX0XXCW genes
###### the "fdata" file read in below really needs updating

mysumtags = function(mat, field = "sgd_gene") {

  geneDir =  "/home/common/barseq_gg_pdata_exp_info/summary_barseq"
  file = paste(geneDir,"2019_apr21_fdata.txt",sep = '/')
  df = read.delim(paste(geneDir,"2019_apr21_fdata.txt",sep = '/'),
    stringsAsFactors = F,check.names = F)

  match = match(rownames(mat),  df$strain)

  vect =  df[match, field]

  x1 = apply(mat, 2, tapply, vect, sum)

  print(dim(mat))

  print(dim(x1))

  x1
}

############ generic function for translatin rownames into different formats
mygetgene = function(mat,match1 = 'strain',get = 'sgd_gene', df = fdata){

  geneDir =  "/home/common/barseq_gg_pdata_exp_info/summary_barseq"
  file = paste(geneDir,"2019_apr21_fdata.txt",sep = '/')
  df = read.delim(paste(geneDir,"2019_apr21_fdata.txt",sep = '/'),
    stringsAsFactors = F,check.names = F)

  wmat = which(names(df)%in% match1)
  names(df)[wmat] = match1
  wget = which(names(df)%in% get)
  match = match(rownames(mat),df[,wmat])
  rownames(mat) = df[match,wget]
  mat = mat[order(rownames(mat)),]
}

############### function for checking dataframes quickly

mylens = function(df){
  lens = sapply(df,unique)
  lens1 = sapply(lens,length)
  lens1
}


############### fitness plots

p10 <- function(matI,
  x,
  sig = 1,
  col = col,
  ylab = "log2:Ratio",
  xlab = "gene",
  las = 2,
  font = 3,
  cex = 0.9,
  cex.main = 2,
  cex.axis = 1.1,
  ...)  {
  col = ifelse (matI[, x] > sig, 1,      ifelse (matI[, x] < -sig, 3, 2))
  w <- which(matI[, x]  > sig | matI[, x]  < -sig, arr.ind = T)

  posw = ifelse (w > nrow(matI) - 0.1 * nrow(matI) ,
    2,
    ifelse (w <  nrow(matI) + 0.1 * nrow(matI), 4 ,  2))
  plot(
    matI[, x],
    col = col ,
    main = paste(colnames(matI)[x], "fitness"),
    ylab = ylab,
    xlab = xlab,
    las = las,
    cex.axis = cex.axis,
    cex.main = cex.main,
    ...
  )
  if (length(w != 0))
    text(
      w,
      matI[w, x],
      names(w),
      pos = posw,
      cex = cex,
      font = font,
      ...
    )
  abline(
    h = sig ,
    col = "red",
    lty = 2,
    lwd = 2
  )
  abline(
    h = -sig ,
    col = "red",
    lty = 2,
    lwd = 2
  )
}

#######################


###### summing up and down tags for raw count matrix  -- modified --

############ generic function for translatin rownames into different formats
mygetgene = function(mat,match1 = 'strain',get = 'sgd_gene', df = fdata){

  geneDir =  "/home/common/barseq_gg_pdata_exp_info/summary_barseq"
  file = paste(geneDir,"2019_apr21_fdata.txt",sep = '/')
  df = read.delim(paste(geneDir,"2019_apr21_fdata.txt",sep = '/'),
    stringsAsFactors = F,check.names = F)

  wmat = which(names(df)%in% match1)
  names(df)[wmat] = match1
  wget = which(names(df)%in% get)
  match = match(rownames(mat),df[,wmat])
  rownames(mat) = df[match,wget]
  mat = mat[order(rownames(mat)),]
}

############### nice tool to plot counts by group or by points

myplotcnts = function(mat,rows,group){
    group = factor(group)
    lev = levels(group)
    plot(t(mat[rows,])~group)
    stripchart(t(mat[rows,])~group,col = 1:length(lev),add=T,vertical = T,pch=19,cex = 1.2,method = "jitter")
    mat[rows,]
  }

#########################
myplotpts = function(mat,rows,group){
function(x) {
col = as.numeric(as.factor(x))
col
}
group = factor(group)
lev = levels(group)
plot(mat[rows,]~group,las=2,main = rows,ylab = "counts",xlab="group")
points(mat[rows,]~group,col=color(group))
mat[rows,]
}

######### various low count filter, none of which i'm comfortable with ##################
myall_less50 = function(mat,limit=50){
  w = which(rowSums(mat < limit) == ncol(mat))

  print(dim(mat))

  if (length(w) > 0) mat = mat[-w,]

  print(dim(mat))

  #print(c(length(w),ncol(mat)))

  mat
}
###############

mylowcnts = function(count,group, limit = 0.5){

  mytapply = function(mat,vect,FUN,margin = 1){
    mat1 = apply(mat,margin,tapply,vect,FUN)
    if(margin == 1) mat1 = t(mat1)
    mat1
  }

  cpm = cpm(count)

  xmax = mytapply(cpm(count),group,max)

  thresh = xmax > limit

  keep = rowSums(thresh) == ncol(xmax)

  count = count[keep,]

  print(dim(count))

  count

}

########convenient to check it data is swaying the right direction overall by summing all of the replicates and sweepingn out the medians
########
mysumcond = function(mat,vect,func,marg = 1){

  row = rownames(mat)

  x1 = apply(mat,marg,tapply,vect,func)

  print(dim(mat))
  print(dim(x1))
  if(marg == 1) x1 = t(x1)
  x1
}


