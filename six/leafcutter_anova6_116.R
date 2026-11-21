# Differences between the six datasets (H1-H3, M1-M3) of one immune cell type, per alternative splicing event, with
# the LeafCutter Dirichlet-multinomial model (as differential_splicing_anova_HN3.R, Ner-Gaon et al. 2023: all junctions
# of the cluster together, min_samples_per_intron = 2, min_samples_per_group = 0, max_cluster_size = 1000,
# min_coverage = 20), but with the six datasets as six groups: full model ~ factor(dataset), null model ~ 1,
# likelihood-ratio test with 5 df. (differential_splicing() builds the design as cbind(1, x), i.e. a single numeric
# covariate - a linear trend over the group numbers - so it is not used directly.) Each model is fitted from the
# 'smart' and from a 'random' initialization and the fit with the larger likelihood ratio is kept, because the smart
# initialization can stop at a poor optimum for multi-group designs.
# Usage: Rscript leafcutter_anova6_116.R <variant> <cell> [cores]
#   input  ensembl116/events/<variant>_6ds/AS_clusters_value.txt, sample_design.tsv
#   output ensembl116/events/<variant>_6ds_anova/leafcutter_anova6/<cell>.tsv (cluster, status, loglr, df, p)
suppressMessages({library(leafcutter); library(foreach); library(doMC)})
args <- commandArgs(trailingOnly = TRUE)
variant <- args[1]; cell <- args[2]; cores <- if (length(args) > 2) as.integer(args[3]) else 4
registerDoMC(cores)
E <- '/gpfs0/tals/projects/Analysis/human_mouse_exons/ensembl116/events/'
datasets <- c('H1', 'H2', 'H3', 'M1', 'M2', 'M3')
MAX_CLUSTER_SIZE <- 1000; MIN_SAMPLES_PER_INTRON <- 2; MIN_SAMPLES_PER_GROUP <- 0; MIN_COVERAGE <- 20; TIMEOUT <- 60

design <- read.table(paste0(E, variant, '_6ds/sample_design.tsv'), header = TRUE, sep = '\t', stringsAsFactors = FALSE)
design <- design[design$cell == cell & design$dataset %in% datasets, ]
counts <- read.table(paste0(E, variant, '_6ds/AS_clusters_value.txt'), header = TRUE, sep = '\t', row.names = 1,
                     check.names = FALSE)
design <- design[design$sample %in% colnames(counts), ]
counts <- as.matrix(counts[, design$sample])
x <- factor(design$dataset, levels = datasets)
clu <- sub('.*:', '', rownames(counts))
clu_names <- unique(clu)
if (length(args) > 3) clu_names <- head(clu_names, as.integer(args[4]))   # test runs only

fit <- function(X, X0, y) {
  best <- NULL
  for (init in c('smart', 'random')) {
    r <- tryCatch(R.utils::withTimeout(leafcutter:::dirichlet_multinomial_anova_mc(X, X0, y, init = init),
                                       timeout = TIMEOUT, onTimeout = 'silent'), error = function(e) NULL)
    if (!is.null(r) && is.finite(r$loglr) && (is.null(best) || r$loglr > best$loglr)) best <- r
  }
  best
}

res <- foreach(cn = clu_names, .combine = rbind, .errorhandling = 'pass') %dopar% {
  out <- function(status, loglr = NA, df = NA, p = NA)
    data.frame(cluster = cn, status = status, loglr = loglr, df = df, p = p, stringsAsFactors = FALSE)
  rows <- which(clu == cn)
  if (length(rows) > MAX_CLUSTER_SIZE) return(out('Too many introns in cluster'))
  if (length(rows) <= 1) return(out('<=1 junction in cluster'))
  y <- t(counts[rows, , drop = FALSE])
  tot <- rowSums(y)
  use <- tot > 0
  if (sum(use) <= 1) return(out('<=1 sample with coverage>0'))
  if (sum(tot[use] >= MIN_COVERAGE) <= 1) return(out('<=1 sample with coverage>min_coverage'))
  xs <- droplevels(x[use]); y <- y[use, , drop = FALSE]
  keep_j <- colSums(y > 0) >= MIN_SAMPLES_PER_INTRON
  if (sum(keep_j) < 2) return(out('<2 introns used in >=min_samples_per_intron samples'))
  y <- y[, keep_j, drop = FALSE]
  ta <- table(xs[rowSums(y) >= MIN_COVERAGE])
  if (sum(ta >= MIN_SAMPLES_PER_GROUP) < 2 || nlevels(xs) < 2) return(out('Not enough valid samples'))
  X <- model.matrix(~ xs); colnames(X) <- NULL
  X0 <- matrix(1, nrow(y), 1)
  r <- fit(X, X0, y)
  if (is.null(r)) return(out('fit failed or timeout'))
  out('Success', r$loglr, r$df, pchisq(2 * r$loglr, r$df, lower.tail = FALSE))
}
od <- paste0(E, variant, '_6ds_anova/leafcutter_anova6/')
dir.create(od, showWarnings = FALSE, recursive = TRUE)
write.table(res, paste0(od, cell, '.tsv'), sep = '\t', quote = FALSE, row.names = FALSE)
cat(variant, cell, ':', nrow(res), 'clusters;', sum(res$status == 'Success'), 'tested\n')
print(table(res$status))
