function (dataList, uniqGeneSelMethod = "dprime", calWithLimma = TRUE, 
  combinedPval = FALSE, filterData = TRUE) 
{
  args = list(dataList)
  if (filterData == TRUE) {
    dataList = inputExpCheck(dataList)
    print("Data list pre-check finished!")
  }
  library(limma)
  gse = names(dataList)
  effect.list = vector("list", length(dataList))
  names(effect.list) = gse
  entrezGenes = NULL
  geneSymbols = NULL
  expressedSamples <- function(x) {
    sum(!is.na(x))
  }
  limmaCal = function(x, design, expressedCaseNum, expressedCtrNum) {
    x = x[expressedCaseNum >= 2 & expressedCtrNum >= 2, 
      ]
    fit <- lmFit(x, design)
    contrast.matrix <- makeContrasts(PD - Control, levels = design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    df = fit2$df.total
    x <- topTable(fit2, coef = 1, adjust = "fdr", number = nrow(x), 
      sort.by = "none")
    cbind(x, df)
  }
  ttestCal = function(case.exp, ctr.exp, caseNum, controlNum) {
    logFC = rowMeans(case.exp, na.rm = TRUE)/rowMeans(ctr.exp, 
      na.rm = TRUE)
    all.exp = cbind(case.exp, ctr.exp)
    avrexpr = rowMeans(all.exp, na.rm = TRUE)
    test.stat <- t(apply(all.exp, 1, function(x) {
      case.exp = x[1:caseNum]
      ctr.exp = x[(caseNum + 1):(caseNum + controlNum)]
      condition = (sum(!is.na(case.exp)) < 2) | (sum(!is.na(ctr.exp)) < 
        2)
      if (condition) {
        return(c(NA, NA, NA, NA))
      }
      else {
        res = t.test(case.exp, ctr.exp)
        t = res$statistic["t"]
        df = res$parameter["df"]
        pval = res$p.value
        adj.P.Val = p.adjust(pval, method = "fdr")
        return(c(t, df, pval, adj.P.Val))
      }
    }))
    colnames(test.stat) = c("t", "df", "pval", "adj.P.Val")
    output = data.frame(logFC, avrexpr, test.stat)
    return(output)
  }
  effectSize <- function(x) {
    pval = x$adj.P.Val
    t = x$t
    n1 = x$expressed.case
    n2 = x$expressed.ctr
    m = x$df
    n = n1 + n2
    df = x$df
    nprime = (n1 * n2)/(n1 + n2)
    d = t/sqrt(nprime)
    cm = gamma(m/2)/gamma((m - 1)/2) * sqrt(2/m)
    dprime = cm * d
    delta = d
    terme1 = m/((m - 2) * nprime)
    vard = terme1 + d^2 * (terme1 * nprime - 1/cm^2)
    vardprime = cm^2 * (terme1 + dprime^2 * (terme1 * nprime - 
      1/cm^2))
    logFC = x$logFC
    y = cbind(logFC, t, pval, n, df, d, vard, dprime, vardprime)
    y = data.frame(y)
    rownames(y) = rownames(x)
    y
  }
  uniqGeneId <- function(x, method) {
    x = x[order(x$Entrez.Gene, abs(x[, method]), decreasing = TRUE), 
      ]
    if (method == "pval") {
      x = x[order(x$Entrez.Gene, x[, method]), ]
    }
    entrezID = unique(x$Entrez.Gene)
    id = match(entrezID, x$Entrez.Gene)
    x = x[id[!is.na(id)], ]
    x
  }
  metaEffect <- function(x) {
    if (nrow(x) < 2) {
      metaZscore = sign(x$t) * qnorm(1 - 0.5 * x$pval)
      metaPval = x$pval
    }
    else {
      x$w = 1/x$vardprime
      wsum = sum(as.vector(x$w))
      mu = sum(as.vector(x$w * x$dprime), na.rm = TRUE)/wsum
      muvar = 1/wsum
      metaZscore = mu/sqrt(muvar)
      metaPval = 2 * (1 - pnorm(abs(metaZscore)))
    }
    c(metaZscore, metaPval)
  }
  metaPval <- function(x) {
    if (nrow(x) < 2) {
      metaZscore = sign(x$t) * qnorm(1 - 0.5 * x$pval)
      metaPval = x$pval
    }
    else {
      pvalonesided = pt(x$t, df = x$df)
      nbstudies = length(pvalonesided)
      nbrepAll = sum(as.vector(x$n))
      weight = sqrt(x$n/nbrepAll)
      metaZscore = sum(weight * qnorm(1 - pvalonesided))
      metaPval = 2 * (1 - pnorm(abs(metaZscore)))
    }
    c(metaZscore, metaPval)
  }
  for (i in 1:length(dataList)) {
    exp.table = dataList[[i]][[1]]
    sampleInfo = dataList[[i]][[2]]
    designInfo = dataList[[i]][[3]]
    exp.dat = exp.table[, 1:sum(sampleInfo)]
    caseSamples = sampleInfo[1][1, 1]
    controlSamples = sampleInfo[2][1, 1]
    case.exp = exp.dat[, c(1:caseSamples)]
    ctr.exp = exp.dat[, -c(1:caseSamples)]
    expressed.case <- apply(case.exp, 1, expressedSamples)
    expressed.ctr <- apply(ctr.exp, 1, expressedSamples)
    meanCase = round(apply(case.exp, 1, mean, na.rm = TRUE), 
      6)
    meanCtr = round(apply(ctr.exp, 1, mean, na.rm = TRUE), 
      6)
    case = c(rep(1, caseSamples), rep(0, ncol(exp.dat) - 
      caseSamples))
    control = c(rep(0, caseSamples), rep(1, ncol(exp.dat) - 
      caseSamples))
    if (calWithLimma) {
      test.stat = limmaCal(exp.dat, designInfo, expressed.case, 
        expressed.ctr)
    }
    else {
      test.stat = ttestCal(case.exp, ctr.exp, caseSamples, 
        controlSamples)
    }
    test.stat = cbind(test.stat, meanCase, meanCtr, expressed.case, 
      expressed.ctr)
    test.stat$NoCase = ncol(case.exp)
    test.stat$NoCtr = ncol(ctr.exp)
    Probesets = exp.table[, (ncol(exp.table) - 1):ncol(exp.table)]
    test.stat = cbind(Probesets, test.stat)
    stat.effect = effectSize(test.stat)
    stat.effect = cbind(exp.table[, (ncol(exp.table) - 1):ncol(exp.table)], 
      stat.effect)
    stat.effect = uniqGeneId(stat.effect, uniqGeneSelMethod)
    rownames(stat.effect) = stat.effect$Entrez.Gene
    effect.list[[i]] = stat.effect[, -1]
    entrezGenes = c(entrezGenes, rownames(stat.effect))
    entrezGenes = unique(entrezGenes)
    geneSymbols = rbind(geneSymbols, stat.effect[, c(1, 
      2)])
    print(paste("Study", i, "statistics calculation completed!", 
      sep = " "))
  }
  print("In-study calculation completed!")
  id = match(entrezGenes, geneSymbols[, "Entrez.Gene"])
  geneSymbols = geneSymbols[id[!is.na(id)], ]
  rownames(geneSymbols) = entrezGenes
  print("Unique genes identified!")
  print("Meta-analysis begins!")
  meta.matrix = matrix(data = NA, nrow = length(entrezGenes), 
    ncol = 2 * length(effect.list) + 4)
  meta.matrix = data.frame(meta.matrix)
  rownames(meta.matrix) = entrezGenes
  colnames(meta.matrix) = c(paste(gse, "FC", sep = "_"), paste(gse, 
    "statPval", sep = "_"), "metaZscore", "metaPval", "effect")
  for (entrezgene in entrezGenes) {
    effect.matrix = NULL
    logFC = NULL
    statPval = NULL
    effect = NULL
    significance = NULL
    for (study in names(dataList)) {
      x = effect.list[[study]][entrezgene, ]
      thisEffect = "?"
      thisFC = NA
      thisPval = NA
      thisSignificance = "?"
      if (!is.na(x$dprime)) {
        effect.matrix = rbind(effect.matrix, x)
        thisEffect = ifelse(x$t > 0, "+", "-")
        thisSignificance = ifelse(x$pval <= 0.05, "!", 
          "#")
        thisFC = x$logFC
        thisPval = x$pval
      }
      effect = paste(effect, thisEffect, sep = "")
      significance = paste(significance, thisSignificance, 
        sep = "")
      logFC = c(logFC, as.numeric(thisFC))
      statPval = c(statPval, as.numeric(thisPval))
    }
    if (combinedPval) {
      meta.matrix[entrezgene, ] = c(logFC, statPval, metaPval(effect.matrix), 
        effect, significance)
    }
    else {
      meta.matrix[entrezgene, ] = c(logFC, statPval, metaEffect(effect.matrix), 
        effect, significance)
    }
  }
  print("Meta-analysis completed!")
  new.meta.matrix = apply(meta.matrix[, -c(ncol(meta.matrix) - 
    1, ncol(meta.matrix))], 2, function(y) {
    as.numeric(y)
  })
  new.meta.matrix = cbind(as.data.frame(new.meta.matrix), 
    as.data.frame(meta.matrix[, c(ncol(meta.matrix) - 1, 
      ncol(meta.matrix))], stringsAsFactors = FALSE))
  colnames(new.meta.matrix)[c(ncol(new.meta.matrix) - 1, ncol(new.meta.matrix))] = c("effect", 
    "significance")
  meta.matrix = new.meta.matrix
  meta.matrix$metaPvalBonf = as.numeric(meta.matrix$metaPval) * 
    nrow(meta.matrix)
  meta.matrix$metaPvalFDR <- p.adjust(as.numeric(meta.matrix$metaPval), method = "fdr")
  meta.matrix$metaPvalBonf[meta.matrix$metaPvalBonf > 1] = 1
  meta.matrix = meta.matrix[order(abs(as.numeric(meta.matrix$metaPval)), 
    decreasing = FALSE), ]
  meta.matrix = cbind(geneSymbols[rownames(meta.matrix), ], 
    meta.matrix[!is.na(id), ])
}