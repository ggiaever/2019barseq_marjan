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






########## added plotting function wrapper

mySARplots = function (obj = out.DESeq2, condRef=condRef, add = F, limit = 60)
{
  if (!I("figures" %in% dir()))
    dir.create("figures", showWarnings = FALSE)

  myMAFP(obj = out.DESeq2$results, condRef, outFile = T, add = add, limit = limit)
  myvolcanoPlot(obj = out.DESeq2$results, condRef, outFile = T, add = add, limit = limit)
  myfitplots(obj = out.DESeq2$results, condRef, outFile = T, col = mycolors)

}

mySAR_QCplots = function (counts,group)
{
  if (!I("figures" %in% dir()))
    dir.create("figures", showWarnings = FALSE)
  mycolors = c(
    "darkorange",
    "dodgerblue",
    "limegreen",
    "navy",
    "mediumpurple"  ,
    "royalblue3",
    "darkolivegreen4",
    "firebrick",
    "cyan",
    "hotpink",
    "plum4",
    "blue",
    "magenta2",
    "skyblue",
    "green",
    "red",
    "steelblue",
    "tomato",
    "purple",
    "yellow"
  )

  mynormdensity(counts,group,col=mycolors)
  myplotMDS(counts,group,col=mycolors)
  myclusterPlot(counts,group,col=mycolors)


}


mywritereport = function (target, counts, out.DESeq2, summaryResults, majSequences,
  workDir, projectName, author, targetFile, rawDir, featuresToRemove,
  varInt, condRef, batch, fitType, cooksCutoff, independentFiltering,
  alpha, pAdjustMethod, typeTrans, locfunc, colors)
{
  rmarkdown::render(input = "2019_apr27_working_deseq_notebook.Rmd", output_file = paste0(projectName,
    "_report.html"), output_dir = workDir, intermediates_dir = workDir,
    knit_root_dir = workDir, run_pandoc = TRUE, quiet = TRUE,
    clean = F)
  cat("HTML report created\n")
}
############### modified SARTools runDESeq2 to allow lfcShrink
myrunDESeq3 = function (counts,
  target,
  varInt,
  batch = NULL,
  locfunc = "median",
  fitType = "parametric",
  pAdjustMethod = "BH",
  cooksCutoff = TRUE,
  independentFiltering = TRUE,
  alpha = 0.05,
  condRef = condRef, Shrink = T, lfcThreshold = 1, type = "apeglm",
  ...)
{
  dds <- DESeqDataSetFromMatrix(countData = counts,
    colData = target,
    design = formula(paste(
      "~", ifelse(!is.null(batch), paste(batch,
        "+"), ""), varInt
    )))
  cat("Design of the statistical model:\n")
  cat(paste(as.character(design(dds)), collapse = " "), "\n")

  dds <- estimateSizeFactors(dds, locfunc = eval(as.name(locfunc)))

  cat("\nNormalization factors:\n")

  print(sizeFactors(dds))

  dds <- estimateDispersions(dds, fitType = fitType)

  dds <- nbinomWaldTest(dds, ...)

  w = which(levels(colData(dds)[, varInt]) == condRef)

  levelRef <- levels(colData(dds)[, varInt])[w]

  levelTest <- levels(colData(dds)[, varInt])[-w]

  nam = resultsNames(dds)[-1]

  results <- list()

  if (Shrink) {

    for (i in 1:length(levelTest)) {
      results[[paste0(levelTest[i], "_vs_", levelRef)]] <-
        DESeq2::lfcShrink(
          dds,
          coef = nam[i],
          lfcThreshold = lfcThreshold, type = type
        )
      cat(paste("Comparison", levelTest[i], "vs", levelRef, "done\n"))
    }

  }

  else {
    for ( i in 1:length(levelTest)) {

      results[[paste0(levelTest[i], "_vs_", levelRef)]] <- DESeq2::results(
        dds,
        contrast = c(varInt, levelTest[i], levelRef),
        pAdjustMethod = pAdjustMethod,
        cooksCutoff = cooksCutoff,
        independentFiltering = independentFiltering,
        alpha = alpha,lfcThreshold = lfcThreshold
      )
      cat(paste("Comparison", levelTest[i], "vs", levelRef, "done\n"))
    }
  }

  return(list(
    dds = dds,
    results = results,
    sf = sizeFactors(dds)
  ))
}
###############
########### SARTools modified density plots -- before AND after normalization on raw counts, not VST variance stabilized transformed (VST); just like using raw counts
######

mynormdensity = function(counts, group, col = mycolors, outFile = T)
{
  ####################################
  mylogdensity <-
    function (counts,
      group,
      col = mycolors,
      main = "density of counts distribution")
    {
      fill = col[as.numeric(as.factor(group))]
      counts <- counts[rowSums(counts) > 0,]
      plot(
        density(log2(counts[, 1] + 1)),
        las = 1,
        lwd = 2,
        main = main,
        xlab = expression(log[2] ~ (raw ~ count + 1)),
        ylim = c(0,
          max(apply(counts, 2, function(x) {
            max(density(log2(x + 1))$y)
          })) * 1.05),
        col = col[as.numeric(as.factor(group))][1]
      )
      for (i in 2:ncol(counts)) {
        lines(density(log2(counts[, i] + 1)), col = col[as.numeric(as.factor(group))][i],
          lwd = 2)
      }
      legend(
        "topleft",
        legend = unique(group),
        fill = unique(col[as.numeric(as.factor(group))]),
        bty = "n"
      )

    }

  ###############
  myseqnorm = function(mat,log = F,...){
    require(edgeR)
    dge = DGEList(mat,remove.zeros = T)
    dge = calcNormFactors(dge,method = 'TMM')
    if(log) mat = cpm(dge,log=T,prior.count = 1,norm.lib.sizes = T)
    if(log == F) mat = cpm(dge,norm.lib.sizes = T)
    mat
  }

  ###

  if (outFile)
    png(
      filename = "figures/dens_norm.png",
      width = 1800 * 2,
      height = 1800,
      res = 300
    )

  par(mfrow = c(1, 2))
  mylogdensity(counts, group, main = "density by pcr date raw cnts")
  mylogdensity(myseqnorm(counts), group, main = "density of normalized counts")

  if (outFile)
    dev.off()

}


###################
##### very slighty modified SARTools clusterPlot to include colors
################################
myclusterPlot = function (counts,
  group,
  col = mycolors,
  outFile = TRUE)
{
  suppressPackageStartupMessages(library(dendextend))

  ###############
  myseqnorm = function(mat,log = F,...){
    require(edgeR)
    dge = DGEList(mat,remove.zeros = T)
    dge = calcNormFactors(dge,method = 'TMM')
    if(log) mat = cpm(dge,log=T,prior.count = 1,norm.lib.sizes = T)
    if(log == F) mat = cpm(dge,norm.lib.sizes = T)
    mat
  }

  ###

  cnts = myseqnorm(counts, log = T)
  hc <- hclust(dist(t(cnts)), method = "ward.D")
  if (outFile)
    png(
      filename = "figures/cluster.png",
      width = 1800,
      height = 1800,
      res = 300
    )

  dend <- as.dendrogram(hc)
  colorCodes <- col[as.numeric(as.factor(group))]
  labels_colors(dend) <- colorCodes[order.dendrogram(dend)]
  par(mfrow = c(1, 1), mar = c(14, 2, 1.4, 1))
  plot(dend,
    ylab = "height",
    las = 2,
    main = "cluster: ward euclidean distance")
  if (outFile)
    dev.off()
}
#######################3
##### very slighty modified SARTools MDS plot to include colors
myplotMDS = function(counts,
  group,
  col = mycolors,
  outFile = T) {
  require(limma)
  ###############
  myseqnorm = function(mat,log = F,...){
    require(edgeR)
    dge = DGEList(mat,remove.zeros = T)
    dge = calcNormFactors(dge,method = 'TMM')
    if(log) mat = cpm(dge,log=T,prior.count = 1,norm.lib.sizes = T)
    if(log == F) mat = cpm(dge,norm.lib.sizes = T)
    mat
  }

  ###

  if (outFile)
    png(
      filename = "figures/MDS.png",
      width = 1800,
      height = 1800,
      res = 300
    )
  cnts = myseqnorm(counts, log = T)
  plotMDS(
    cnts,
    col = as.numeric(as.factor(group)),
    main = 'MDS by group',
    gene.selection = 'pairwise'
  )
  legend(
    "bottomright",
    legend = unique(group),
    lty = 1,
    col = unique(col[as.numeric(as.factor(group))]),
    lwd = 2,
    bty = "n"
  )
  if (outFile)
    dev.off()
}

##########################
# our standard fitness plots: logRatio vs gene
# would be nice to incorporate removal of big outliers for the report...

myfitplots = function(obj = out.DESeq2$results,
  condRef = condRef,
  outFile = T, col = mycolors) {
  # fitness plots function
  #

  palette(mycolors)
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

  ## clugy function to plot residuals using p10 fitness plot function as it need a matrix
  mydftomat = function(res) {
    n = sapply(res, is.numeric)
    mx = as.matrix(res[, n])
    mx = mx[order(rownames(mx)),]
    #p10(mx, 2)
    mx
  }

  lres = obj[grep(condRef, names(obj))]

  nam = names(lres)

  lres = lapply(lres, myres)

  nam = names(lres)

  lfit = NULL

  for (i in 1:length(nam))
    lfit[[i]] = mydftomat(lres[[i]])

  w = sapply(lfit, function(x)
    x = grep("log2_ratio", colnames(x)))

  names(lfit) = nam


  for (i in 1:length(lfit)) colnames(lfit[[i]]) = gsub("log2_ratio", nam[i], colnames(lfit[[i]]))

  nrow = length(lfit)

  if (outFile)
    png(
      filename = "figures/FDPlots.png",
      width = 2250,
      height = 2250 * nrow,
      res = 300
    )
  par(mfrow = c(nrow, 1), pch = 19)
  for (i in 1:length(lfit)) {
    p10(lfit[[i]], w[i])
  }

  for (i in 1:length(lfit)) p10(lfit[[i]], w[i])
  if (outFile)
    dev.off()
}

################ sets up the results from differents objects for down stream plotting
################ awkward but helpful

myres = function(res, delete = T) {
  require(dplyr)
  if (class(res) != "data.frame")
    res = data.frame(res, stringsAsFactors = F)

  edge = which(names(res) %in% c('logFC', 'logCPM', 'FDR'))
  if (length(edge) == 3) {
    names(res)[edge] = c('log2_ratio', 'logCPM', 'padj')
    res$log2_ratio = -res$log2_ratio
  }

  dseq = which(names(res) %in% c('log2FoldChange', 'padj', 'svalue'))
  if (length(dseq) == 2) {
    res$logCPM = log2(res[, 'baseMean'] + 1)
    names(res)[dseq] = c('log2_ratio', 'padj')
    res$log2_ratio = -res$log2_ratio
  }

  dimma = which(names(res) %in% c('logFC', 'AveExpr', 'adj.P.Val'))
  if (length(dimma) == 3) {
    names(res)[dimma] = c('log2_ratio', 'logCPM', 'padj')
    res$log2_ratio = -res$log2_ratio
  }

  if (all(c(length(edge), length(dseq), length(dimma)) == 0))
    stop(paste("can't identify object class"))

  top = res

  top$gene = rownames(top)

  top$sig = abs(top$log2_ratio) > 1 & top$padj < 0.05

  wna = which(is.na(top$padj))

  if (length(wna) > 0 & delete)
    top =  top[-wna,]

  cat('no. of sig genes ', length(which(top$sig == 1)), "\n")

  top = top %>% arrange(padj)

  rownames(top) = top$gene

  #print(myvolc(top))

  top
}

###### plot logRatios as a function of average logCPM for each gene
###### similar to MA plots but with more detail
## add = T means to add points if there are no significant fitness defects
## setting to F means don't add more points beyond those called significant
## in either case, the limit of labeled points is set to 40 to avoid messy plots

myggFP = function(res,
  x = 'logCPM',
  y = 'log2_ratio',
  limit = 60, add = F,
  main = 'log2 ratio vs mean counts (FDR < 0.05)') {
  library(ggplot2)
  require(ggrepel)
  require(dplyr)
  wavg = which(names(res) == x)
  names(res)[wavg] = 'logCPM'
  wlr = which(names(res) == y)
  names(res)[wlr] = 'log2_ratio'
  q = quantile(res$logCPM, probs = seq(0, 1, 0.1))
  cuts = cut(res$logCPM, q, include.lowest = T)
  lev = round(2 ^ q[-length(q)], -2)
  lev[1] = paste0('<', lev[2])
  lev[length(lev)] = paste0('>', lev[length(lev) - 1])
  levels(cuts) = lev
  res$bin = cuts
  wsig = which(res$sig == T)

  if (length(wsig) > limit ) {
    res$sig = F
    wsig = order(res$padj)[1:limit]
    res$sig[wsig] = T
  }

  if (length(wsig) < limit & add) {
    wsig = order(res$padj)[1:limit]
    res$sig[wsig] = T
  }


  print(length(wsig))
  g = ggplot(res, aes(x = bin, y = log2_ratio, fill = bin)) + geom_boxplot() + theme_bw() +
    ggtitle(main) + theme(
      plot.title = element_text(
        lineheight = 1.2,
        size = 14,
        face = "bold"
      ),
      axis.title = element_text(size = 14, face = 'bold'),
      legend.position = 'none',
      axis.text = element_text(size = 14)
    ) +  labs(x = "decile counts")
  g1 = g + geom_hline(
    aes(yintercept = 1),
    col = 'red',
    linetype = 'dashed',
    size = 1
  ) +
    geom_hline(
      aes(yintercept = -1),
      col = 'red',
      linetype = 'dashed',
      size = 1
    )
  if (length(wsig) > 0) {
    g1 = g1 + geom_text_repel(
      data = subset(res, sig == TRUE),
      aes(x = bin, y = log2_ratio, label = gene),
      point.padding = 0.25,
      segment.alpha = 0.2
    )
    g1 + geom_point(
      data = subset(res, sig == TRUE),
      aes(
        x = bin,
        y = log2_ratio,
        col = factor(sign(log2_ratio)),
        shape =  factor(sign(log2_ratio))
      ),
      size = 2.5
    )
  }
}

###### wrapper for above; myggFP

myMAFP = function(obj = out.DESeq2$results,
  condRef = condRef,
  outFile = T, add = F, limit = 60) {
  # gut check to see if logRatios are coming from lowcounts

  lres = obj[grep(condRef, names(obj))]

  nam = names(lres)

  lres = lapply(lres, myres)

  nam = names(lres)

  nrow <- length(lres)

  main = 'log2 ratio vs mean counts (FDR < 0.05)'

  main2 = paste(names(lres), main, sep = ": ")

  lggfp = NULL

  for (i in 1:length(lres))
    lggfp[[i]] = myggFP(lres[[i]], main = paste(main2[i]), limit = limit,add = add)

  if (outFile)
    png(
      filename = "figures/MAPlot.png",
      width = 1800,
      height = 1800 * nrow,
      res = 300
    )

  multiplot(plotlist = lggfp)

  if (outFile)
    dev.off()
}
## volcano plot of logRatios -- highly susceptible to outliers
## add = T means to add points if there are no significant fitness defects
## setting to F means don't add more points beyond those called significant
## in either case, the limit of labeled points is set to 40 to avoid messy plots

myvolc = function(res, main = 'log2 ratio vs FDR',add = F, limit = 40){

  wsig = which(res$sig == T)

  if (length(wsig) > limit) {
    res$sig = F
    wsig = order(res$padj)[1:limit]
    res$sig[wsig] = T
  }

  if (length(wsig) < limit & add) {
    wsig = order(res$padj)[1:limit]
    res$sig[wsig] = T
  }

  wsig = which(res$sig == T)
  print(length(wsig))
  require(ggplot2)
  require(ggrepel)
  volc = ggplot(res, aes(-log10(padj),log2_ratio)) + geom_point(aes(col = sig)) + theme_bw() +
    scale_color_manual(values=c("navy", "limegreen")) + theme(legend.position = 'none') + ggtitle(main)
  volc2 = volc + geom_text_repel(data = subset(res,sig == TRUE),aes(-log10(padj),log2_ratio,label = gene))
  print(volc2 + geom_hline(aes(yintercept = 1),col = 'red',linetype = 'dashed',size = 1) +
      geom_hline(aes(yintercept = -1),col = 'red',linetype = 'dashed',size = 1) +
      theme(plot.title = element_text(lineheight=1.2, size = 16,face="bold")) +
      theme(axis.title = element_text(size = 15,face = 'bold'), axis.text = element_text(size = 15)))
  volc3 = volc2 + geom_hline(aes(yintercept = 1),col = 'red',linetype = 'dashed',size = 1) +
    geom_hline(aes(yintercept = -1),col = 'red',linetype = 'dashed',size = 1) +
    theme(plot.title = element_text(lineheight=1.2, size = 16,face="bold")) +
    theme(axis.title = element_text(size = 15,face = 'bold'), axis.text = element_text(size = 15))
  volc3
}

###### wrapper for above; myvolc
####################
myvolcanoPlot = function(obj = out.DESeq2$results,
  condRef = condRef,
  outFile = T,add = F, limit = 60) {
  lres = obj[grep(condRef, names(obj))]

  nam = names(lres)

  lres = lapply(lres, myres)

  nrow <- length(lres)

  main = 'log2 ratio vs FDR'

  main2 = paste(names(lres), main, sep = ": ")

  lvolc = NULL

  for (i in 1:length(lres))
    lvolc[[i]] = myvolc(lres[[i]], main = paste(main2[i]),limit = limit,add = add)


  if (outFile)
    png(
      filename = "figures/volcanoPlot.png",
      width = 1800,
      height = 1800 * nrow,
      res = 300
    )

  multiplot(plotlist = lvolc)

  if (outFile)
    dev.off()
}

#######stolen ggplot multiplot function for ggplots above

multiplot <-
  function(...,
    plotlist = NULL,
    file,
    cols = 1,
    layout = NULL) {
    library(grid)
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    numPlots = length(plots)
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots / cols)),
        ncol = cols,
        nrow = ceiling(numPlots / cols))
    }
    if (numPlots == 1) {
      print(plots[[1]])
    }

    else {
      # Set up the page
      grid.newpage()
      pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <-
          as.data.frame(which(layout == i, arr.ind = TRUE))
        print(plots[[i]],
          vp = viewport(
            layout.pos.row = matchidx$row,
            layout.pos.col = matchidx$col
          ))
      }
    }
  }

mydftomat = function(res) {
  n = sapply(res, is.numeric)
  mx = as.matrix(res[, n])
  mx = mx[order(rownames(mx)),]
  #p10(mx, 2)
  mx
}