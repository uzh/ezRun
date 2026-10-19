###################################################################
# Functional Genomics Center Zurich
# This code is distributed under the terms of the GNU General
# Public License Version 3, June 2007.
# The terms are available here: http://www.gnu.org/licenses/gpl.html
# www.fgcz.ch

## Downstream analysis of HUMAnN tables: GSEA on MaAsLin3 statistics and
## species-contribution decomposition of stratified features.

humann_feature_ids <- function(features) {
  sub(":.*$", "", features)
}

humann_load_stratified <- function(path) {
  m <- read.table(path, header = TRUE, sep = "\t", check.names = FALSE,
                  comment.char = "", quote = "", row.names = 1)
  colnames(m) <- sub("_concat.*$", "", colnames(m))
  colnames(m) <- sub("_(Abundance|Abundance-CPM|Coverage|Abundance-RPKs).*$",
                     "", colnames(m))
  isStrat <- grepl("|", rownames(m), fixed = TRUE)
  if (!any(isStrat)) {
    return(NULL)
  }
  m <- as.matrix(m[isStrat, , drop = FALSE])
  parts <- strsplit(rownames(m), "|", fixed = TRUE)
  taxon <- vapply(parts, function(x) x[2], character(1))
  taxon <- sub("^s__", "", sub("\\.t__.*$", "", taxon))
  list(values = m,
       feature = vapply(parts, function(x) x[1], character(1)),
       taxon = taxon)
}

humann_read_metaphlan_species <- function(paths) {
  profiles <- list()
  for (s in names(paths)) {
    df <- read.delim(paths[[s]], comment.char = "#", header = FALSE,
                     stringsAsFactors = FALSE, quote = "")
    keep <- grepl("\\|s__[^|]+$", df[[1]])
    v <- as.numeric(df[[3]][keep])
    names(v) <- sub("^.*\\|s__", "", df[[1]][keep])
    profiles[[s]] <- v[!duplicated(names(v))]
  }
  species <- unique(unlist(lapply(profiles, names)))
  m <- matrix(0, nrow = length(species), ncol = length(profiles),
              dimnames = list(species, names(profiles)))
  for (s in names(profiles)) {
    m[names(profiles[[s]]), s] <- profiles[[s]]
  }
  m
}

humann_rank_statistic <- function(res, grouping) {
  sub <- res[res$metadata == grouping & res$model == "abundance", , drop = FALSE]
  nullHyp <- rep(0, nrow(sub))
  if ("null_hypothesis" %in% colnames(sub)) {
    nullHyp <- suppressWarnings(as.numeric(sub$null_hypothesis))
    nullHyp[is.na(nullHyp)] <- 0
  }
  stat <- (as.numeric(sub$coef) - nullHyp) / as.numeric(sub$stderr)
  names(stat) <- sub$feature
  stat <- stat[is.finite(stat)]
  stat[!duplicated(names(stat))]
}

humann_read_metacyc_sets <- function(path) {
  rows <- list()
  for (line in readLines(path, warn = FALSE)) {
    tokens <- strsplit(trimws(line), "\\s+")[[1]]
    if (length(tokens) < 2) {
      next
    }
    reactions <- sub("^-", "", tokens[-1])
    reactions <- unique(reactions[!reactions %in% c("(", ")", ",", "+", "")])
    if (length(reactions) == 0) {
      next
    }
    rows[[length(rows) + 1]] <- data.frame(term = tokens[1], gene = reactions,
                                           stringsAsFactors = FALSE)
  }
  do.call(rbind, rows)
}

humann_read_metacyc_names <- function(path) {
  if (!file.exists(path)) {
    return(NULL)
  }
  df <- read.delim(path, header = FALSE, stringsAsFactors = FALSE, quote = "")
  data.frame(term = df[[1]], name = df[[2]], stringsAsFactors = FALSE)
}

humann_read_gmt <- function(paths) {
  term2gene <- list()
  term2name <- list()
  for (p in paths) {
    for (line in readLines(p, warn = FALSE)) {
      tokens <- strsplit(line, "\t", fixed = TRUE)[[1]]
      if (length(tokens) < 3) {
        next
      }
      genes <- unique(tokens[-(1:2)])
      genes <- genes[nzchar(genes)]
      term2gene[[length(term2gene) + 1]] <- data.frame(term = tokens[1], gene = genes,
                                                       stringsAsFactors = FALSE)
      term2name[[length(term2name) + 1]] <- data.frame(term = tokens[1], name = tokens[2],
                                                       stringsAsFactors = FALSE)
    }
  }
  if (length(term2gene) == 0) {
    return(list(term2gene = NULL, term2name = NULL))
  }
  list(term2gene = do.call(rbind, term2gene),
       term2name = unique(do.call(rbind, term2name)))
}

humann_run_gsea <- function(stats, term2gene, term2name = NULL,
                            minSize = 5, maxSize = 500) {
  if (!requireNamespace("clusterProfiler", quietly = TRUE)) {
    stop("clusterProfiler is required for functional enrichment.")
  }
  stats <- sort(stats, decreasing = TRUE)
  term2gene <- term2gene[term2gene$gene %in% names(stats), , drop = FALSE]
  if (nrow(term2gene) == 0) {
    return(data.frame())
  }
  set.seed(42)
  res <- clusterProfiler::GSEA(geneList = stats, TERM2GENE = term2gene,
                               TERM2NAME = if (is.null(term2name)) NA else term2name,
                               minGSSize = minSize, maxGSSize = maxSize,
                               pvalueCutoff = 1, eps = 0, seed = TRUE,
                               by = "fgsea", verbose = FALSE)
  as.data.frame(res)
}

humann_gsea_plot <- function(gsea, title, nShow = 20) {
  if (is.null(gsea) || nrow(gsea) == 0) {
    return(NULL)
  }
  d <- head(gsea[order(gsea$p.adjust, gsea$pvalue), , drop = FALSE], nShow)
  label <- ifelse(d$Description == d$ID, d$ID, paste0(d$Description, " (", d$ID, ")"))
  d$label <- factor(label, levels = rev(label))
  ggplot2::ggplot(d, ggplot2::aes(x = NES, y = label, size = setSize, color = p.adjust)) +
    ggplot2::geom_point() +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    ggplot2::scale_color_gradient(low = "firebrick", high = "grey70") +
    ggplot2::labs(title = title, x = "Normalised enrichment score", y = NULL,
                  size = "Set size", color = "Adjusted p") +
    ggplot2::theme_bw(base_size = 11)
}

.humann_group_effect <- function(y, md, grouping, sampleGroup, covariates,
                                 keep = NULL) {
  na <- c(coef = NA_real_, p = NA_real_)
  df <- md[, c(grouping, covariates), drop = FALSE]
  df$y <- y
  if (!is.null(keep)) {
    df <- df[keep, , drop = FALSE]
  }
  df <- df[is.finite(df$y), , drop = FALSE]
  if (nrow(df) < 6 || length(unique(df[[grouping]])) < 2) {
    return(na)
  }
  rhs <- paste(sprintf("`%s`", c(grouping, covariates)), collapse = " + ")
  fit <- tryCatch(lm(as.formula(paste("y ~", rhs)), data = df),
                  error = function(e) NULL)
  if (is.null(fit)) {
    return(na)
  }
  cf <- summary(fit)$coefficients
  term <- c(paste0(grouping, sampleGroup), sprintf("`%s`%s", grouping, sampleGroup))
  term <- term[term %in% rownames(cf)]
  if (length(term) == 0) {
    return(na)
  }
  c(coef = cf[term[1], "Estimate"], p = cf[term[1], "Pr(>|t|)"])
}

humann_contributors <- function(strat, features, md, grouping, sampleGroup,
                                refGroup, covariates = character(0),
                                speciesAbund = NULL, nTop = 5) {
  covariates <- setdiff(intersect(covariates, colnames(md)), grouping)
  samples <- intersect(colnames(strat$values), rownames(md))
  md <- md[samples, , drop = FALSE]
  md[[grouping]] <- factor(as.character(md[[grouping]]),
                           levels = c(refGroup, sampleGroup))
  isSample <- md[[grouping]] == sampleGroup
  values <- strat$values[, samples, drop = FALSE]
  taxonTotals <- rowsum(values, strat$taxon)
  rows <- list()
  plotRows <- list()
  for (f in features) {
    idx <- strat$feature == f
    if (!any(idx)) {
      next
    }
    fm <- rowsum(values[idx, , drop = FALSE], strat$taxon[idx])
    total <- sum(fm)
    ranked <- names(sort(rowMeans(fm), decreasing = TRUE))
    top <- head(setdiff(ranked, "unclassified"), nTop)
    for (tx in top) {
      contrib <- fm[tx, ]
      speciesTotal <- taxonTotals[tx, ]
      present <- contrib > 0 & speciesTotal > 0
      ce <- .humann_group_effect(log2(contrib + 1), md, grouping, sampleGroup, covariates)
      se <- .humann_group_effect(log2(contrib / speciesTotal), md, grouping, sampleGroup,
                                 covariates, keep = present)
      ae <- c(coef = NA_real_, p = NA_real_)
      if (!is.null(speciesAbund) && tx %in% rownames(speciesAbund)) {
        a <- speciesAbund[tx, match(samples, colnames(speciesAbund))]
        ae <- .humann_group_effect(log2(a * 1e4 + 1), md, grouping, sampleGroup, covariates)
      }
      rows[[length(rows) + 1]] <- data.frame(
        feature = f, taxon = tx,
        meanRef = mean(contrib[!isSample]), meanSample = mean(contrib[isSample]),
        shareOfTotal = if (total > 0) sum(contrib) / total else NA_real_,
        presentRef = sum(contrib[!isSample] > 0),
        presentSample = sum(contrib[isSample] > 0),
        contributionCoef = ce[["coef"]], contributionP = ce[["p"]],
        withinSpeciesShareCoef = se[["coef"]], withinSpeciesShareP = se[["p"]],
        abundanceCoef = ae[["coef"]], abundanceP = ae[["p"]],
        stringsAsFactors = FALSE)
    }
    other <- setdiff(rownames(fm), c(top, "unclassified"))
    for (g in c(top, if (length(other) > 0) "Other species",
                if ("unclassified" %in% rownames(fm)) "unclassified")) {
      v <- if (g == "Other species") colSums(fm[other, , drop = FALSE]) else fm[g, ]
      plotRows[[length(plotRows) + 1]] <- data.frame(
        feature = f, taxon = g, group = c(refGroup, sampleGroup),
        mean = c(mean(v[!isSample]), mean(v[isSample])),
        stringsAsFactors = FALSE)
    }
  }
  list(table = do.call(rbind, rows), plotData = do.call(rbind, plotRows))
}

humann_contributor_plot <- function(plotData) {
  if (is.null(plotData) || nrow(plotData) == 0) {
    return(NULL)
  }
  plotData$taxon <- factor(plotData$taxon, levels = unique(plotData$taxon))
  plotData$feature <- factor(plotData$feature, levels = unique(plotData$feature))
  ggplot2::ggplot(plotData, ggplot2::aes(x = mean, y = group, fill = taxon)) +
    ggplot2::geom_col(width = 0.7) +
    ggplot2::facet_wrap(~ feature, ncol = 1, scales = "free_x") +
    ggplot2::labs(x = "Group mean (table units)", y = NULL, fill = NULL) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(legend.position = "bottom")
}
