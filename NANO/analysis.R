library(tidyverse)
library(GenomicRanges)
library(data.table)
library(patchwork)
library(boot)
library(multcomp)
library(limma)
library(ggrepel)
library(rjson)
library(gtools)
library(RColorBrewer)

setwd("/scratch/fmorandi/internal/John/cGAS_KO/PAPER_CODE/NANO")

paths = list()
paths$data = "./data"
paths$anno_RE = "/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/RepeatMaskerOut/GRCm39_repeats.saf"
paths$results = "./results"
dir.create(paths$results, showWarnings = F)

#### PLOTTING SETTINGS ####

w_in = 7.5 # 8.5 without margin
h_in = 10 # 11 without margin

#### FUNCTIONS ####

read_nanopore = function(nanopath, anno_re, meta) {
  # See modbam2bed docs for meaning of each column
  stats = list()
  # Make re ranges
  anno_ranges = GRanges(
    seqnames = anno_re$chr,
    names = anno_re$fragment_id,
    ranges = IRanges(
      start = anno_re$start,
      end = anno_re$end))
  # Lists to store one item for each file
  data = list()
  cnames = c("chr", "start", "end", "mod_type", "score", "strand", 
             "cov", "perc", "unmod", "mod", "filt", "nocall", "alt")
  files = dir(paste0(nanopath, "/basemods"), pattern="bed", full.names=T)
  fnames = str_extract(files, "/([[:alnum:]]+).cpg.acc.bed", group=1)
  stats$anycov = c() # Any coverage of the CpG
  stats$quant = c() # Covered and methyl quantified
  stats$noquant = c() # Covered but no quant
  stats$subdel = c() # Substitution or deletion
  stats$anycovRE = c() # Same as above but after collapsing to REs
  stats$quantRE = c()
  stats$noquantRE = c()
  stats$subdelRE = c()
  for (i in 1:length(files)) {
    # Load modbases file
    tmp = fread(files[i], data.table=F)
    tmp = tmp[,-c(7:9)]
    colnames(tmp) = cnames
    # Create a few utility variables
    tmp = tmp %>%
      mutate(quant = unmod + mod + alt) %>% # Quantified CpGs
      mutate(noquant = filt + nocall) %>% # Not quantified
      mutate(subdel = cov - quant - noquant) # Substitution or deletion
    # Save some stats
    stats$anycov[fnames[i]] = sum(tmp$cov)
    stats$quant[fnames[i]] = sum(tmp$quant)
    stats$noquant[fnames[i]] = sum(tmp$noquant)
    stats$subdel[fnames[i]] = sum(tmp$subdel)
    # Collapse to repeat instances
    cpg_ranges = GRanges(
      seqnames = tmp$chr,
      ranges = IRanges(
        start = tmp$start,
        end = tmp$end))
    overlaps = findOverlaps(cpg_ranges, anno_ranges)
    overlaps = data.frame(overlaps)
    overlaps$rep_id=anno_re$rep_id[overlaps$subjectHits]
    tmp[overlaps$queryHits, "rep_id"] = overlaps$rep_id
    tmp = tmp %>%
      dplyr::filter(!is.na(rep_id)) %>%
      dplyr::select(rep_id, cov, quant, noquant, subdel, unmod, mod, alt)
    tmp[c("superf", "fam", "rep_id")] = split_rep_id(tmp$rep_id)
    tmp = tmp %>%
      dplyr::select(-rep_id) %>%
      group_by(superf, fam) %>%
      summarize_all(sum) %>%
      dplyr::filter(cov > 0)
    # Save some stats
    stats$anycovRE[fnames[i]] = sum(tmp$cov)
    stats$quantRE[fnames[i]] = sum(tmp$quant)
    stats$noquantRE[fnames[i]] = sum(tmp$noquant)
    stats$subdelRE[fnames[i]] = sum(tmp$subdel)
    data[[fnames[i]]] = tmp
  }
  # Make homogenous lists
  data2 = list(
    met = list(),
    cov = list()
  )
  for (bc in names(data)) {
    tmp = subset(data[[bc]], quant > 0)
    data2$met[[bc]] = data2$cov[[bc]] = tmp[c("superf", "fam")]
    data2$met[[bc]][bc] = tmp$mod
    data2$cov[[bc]][bc] = tmp$quant
  }
  # Concatenate lists to tables and rename 
  name_conv = c("superf", "fam", meta$Sample)
  names(name_conv) = c("superf", "fam", meta$Barcode)
  for (n2 in names(data2)) {
    data2[[n2]] = Reduce(function(x, y) merge(x, y, by=c("superf", "fam"), all=T), data2[[n2]])
    colnames(data2[[n2]]) = name_conv[colnames(data2[[n2]])]
    data2[[n2]][is.na(data2[[n2]])] = 0
  }
  data2$re_info = data2$cov[, c("superf", "fam")]
  rownames(data2$re_info) = paste0("re", 1:nrow(data2$re_info))
  for (n2 in c("met", "cov")) {
    rownames(data2[[n2]]) = rownames(data2$re_info)
    data2[[n2]]$superf = NULL
    data2[[n2]]$fam = NULL
    data2[[n2]] = data.frame(t(data2[[n2]]))
  }
  out = list(
    data = data2,
    stats = stats
  )
  return(out)
}

split_rep_id = function(rep_id) {
  # Split rep_id into superf, fam, id properly
  tmp = rep_id
  rep_id = str_extract(tmp, "[[:digit:]]+$")
  tmp = str_replace(tmp, "/[[:digit:]]+$", "")
  fam = str_extract(tmp, "[^/]+$")
  tmp = str_replace(tmp, "/[^/]+$", "")
  superf = tmp
  out = data.frame(
    superf = superf,
    fam = fam,
    rep_id = rep_id
  )
  return(out)
}

my_volcano = function(table, logFC_col, pval_col, sig_col, ymax = NULL) {
  xmax = quantile(abs(table[, logFC_col]), 1-5/nrow(table))
  if (is.null(ymax)) {
    ymax = quantile(-log10(0.00001+table[, pval_col]), 1-1/nrow(table))
  }
  tmp = table %>%
    group_by_at(sig_col) %>%
    summarize(n=n()) %>%
    mutate(x = c(-xmax*0.9, 0, xmax*0.9)) %>%
    mutate(hjust = c(0, 0.5, 1))
  p = ggplot(table, aes(x=.data[[logFC_col]], y=-log10(.data[[pval_col]]), color=.data[[sig_col]]))+
    geom_point(size=0.1)+
    lims(x=c(-xmax, xmax), y=c(0,ymax))+
    geom_text(data=tmp, aes(x=x, y=ymax*0.95, label=n))+
    scale_color_manual(values=c("#5555cc", "#999999", "#cc5555"))+
    labs(x="logFC", y="-log10(pval)", color="Significance")
  return(p)
}

#### PREPROCESS ####

# # Load annotation
# anno_re = fread(paths$anno_RE, col.names=c("rep_id", "chr", "start", "end", "strand"),
#                 skip=1, data.table=F)
# anno_re = distinct(anno_re, chr, start, end, .keep_all=T)
# anno_re$fragment_id = paste0(anno_re$rep_id, ".", 1:nrow(anno_re))
# 
# anno_re2 = anno_re %>%
#   mutate(length = end - start)
# anno_re2[c("superf", "fam", "rep_id")] = split_rep_id(anno_re2$rep_id)
# anno_re2 = anno_re2 %>%
#   group_by(superf, fam) %>%
#   summarize(n=n(), mean_length=mean(length))
# write.table(anno_re2, paste0(paths$results, "/re_stats.tsv"), sep="\t")
# 
# # Load metadata
# meta = read.table(paste0("./meta.txt"), header=T)
# 
# # Load data and related stats
# tmp = read_nanopore(paths$data, anno_re, meta)
# data = tmp$data
# stats = tmp$stats
# 
# # Save stats
# sink(paste0(paths$results, "/sample_stats.txt"))
# cat(toJSON(stats, indent=1))
# closeAllConnections()
# 
# # Define sufficient coverage
# # Hmmm, there seem to be big coverage differences hmmmmmm
# hist(log10(colMeans(data$cov)), breaks=50)
# hist(log10(sapply(data$cov, min)), breaks=50)
# th_min = 10 # Minimum coverage requirement
# th_mean = 50 # Mean coverage requirement
# keep1 = sapply(data$cov, min) > th_min
# keep2 = colMeans(data$cov) > th_mean
# sprintf("%i, %i, %i", sum(keep1), sum(keep2), sum(keep1 & keep2))
# passing = keep1 & keep2
# data$met_hc = data$met[, passing]
# data$cov_hc = data$cov[, passing]

# save(anno_re2, data, meta, stats, file=paste0(paths$results, "/prepro.Rdata"))

#### LOAD DATA ####

load(file=paste0(paths$results, "/prepro.Rdata"))

anno_re2 = anno_re2 %>%
  mutate(class = gsub("/.*", "", superf), .before = superf)

meta$Genotype = factor(meta$Genotype)
meta$Depth = rowSums(data$cov)
plot(rowSums(data$cov), rowSums(data$cov_hc))
abline(b=1, a=0)

# Select data (All or high coverage REs only)
Met = data$met_hc
Cov = data$cov_hc
Perc = Met / Cov
re_info = data$re_info %>%
  mutate(class = gsub("/.*", "", superf), .before = superf)
re_info = re_info[colnames(Met), ]

#### PCA ####

tmp = Perc[, colSums(is.na(Perc)) == 0]
all.equal(rownames(tmp), meta$Sample)
pca = prcomp(tmp, scale=T)
pca = cbind(meta, pca$x)

ggplot(pca, aes(x=PC1, y=PC2, color=Genotype))+
  geom_point()
ggsave(paste0(paths$results, "/pca.pdf"), width = 0.5*w_in, height = 0.25*h_in)

#### DIFFERENTIAL METHYLATION ####

design = model.matrix(~0+Genotype, data=meta)
contrasts = makeContrasts(
  cGASKOvsWT = GenotypecGASKO - GenotypeWT,
  levels = design)

##### 5mC binomial #####

results5mC_binom = data.frame()
for (re in colnames(Met)) {
  fit = glm(as.matrix(cbind(Met[re], Cov[re]-Met[re]))~0+Genotype,
            data=meta, family = binomial())
  tmp = summary(glht(fit, t(contrasts)))$test
  tmp = data.frame(
    coef = tmp$coefficients,
    pval = tmp$pvalues)
  results5mC_binom[re, colnames(tmp)] = tmp
}

results5mC_binom = merge(re_info, results5mC_binom, by=0)
results5mC_binom$padj = p.adjust(results5mC_binom$pval, method="BH")
results5mC_binom$sig = "Not Sig."
results5mC_binom[results5mC_binom$padj < 0.05 & results5mC_binom$coef > 0, "sig"] = "Sig. Up"
results5mC_binom[results5mC_binom$padj < 0.05 & results5mC_binom$coef < 0, "sig"] = "Sig. Down"
results5mC_binom$sig = factor(results5mC_binom$sig, levels = c("Sig. Down", "Not Sig.", "Sig. Up"))

# Add meanCov and meanPerc
results5mC_binom$meanCov = colMeans(Cov)[results5mC_binom$Row.names]
results5mC_binom$meanPerc = colMeans(Perc, na.rm = T)[results5mC_binom$Row.names]

##### 5mC linear #####

results5mC_lin = data.frame()
for (re in colnames(Met)) {
  fit = glm(unlist(Perc[re])~0+Genotype, data=meta, family = gaussian())
  tmp = summary(glht(fit, t(contrasts)))$test
  tmp = data.frame(
    coef = tmp$coefficients,
    pval = tmp$pvalues)
  results5mC_lin[re, colnames(tmp)] = tmp
}

results5mC_lin = merge(re_info, results5mC_lin, by=0)
results5mC_lin$padj = p.adjust(results5mC_lin$pval, method="BH")
results5mC_lin$sig = "Not Sig."
results5mC_lin[results5mC_lin$padj < 0.05 & results5mC_lin$coef > 0, "sig"] = "Sig. Up"
results5mC_lin[results5mC_lin$padj < 0.05 & results5mC_lin$coef < 0, "sig"] = "Sig. Down"
results5mC_lin$sig = factor(results5mC_lin$sig, levels = c("Sig. Down", "Not Sig.", "Sig. Up"))

# Add meanCov and meanPerc
results5mC_lin$meanCov = colMeans(Cov)[results5mC_lin$Row.names]
results5mC_lin$meanPerc = colMeans(Perc, na.rm = T)[results5mC_lin$Row.names]

##### Compare methods #####

tmp = data.frame(
  "coef_lin" = results5mC_lin$coef,
  "coef_binom" = results5mC_binom$coef,
  "pval_lin" = results5mC_lin$pval,
  "pval_binom" = results5mC_binom$pval,
  "meanCov" = results5mC_binom$meanCov,
  "meanPerc" = results5mC_binom$meanPerc) 
p1=ggplot(tmp, aes(coef_lin, coef_binom))+
  geom_point()+
  ylim(-1.5, 1.5)+
  labs(x="Coef (LM)", y="Coef (Binomial GLM)")+
  geom_smooth(method="lm")+
  stat_cor()
p2=ggplot(tmp, aes(-log10(pval_lin), -log10(pval_binom)))+
  geom_point()+
  labs(x="Pval (LM)", y="Pval (Binomial GLM)")+
  geom_smooth(method="lm")+
  geom_abline(color="red", linetype="dashed")+
  stat_cor()
p1|p2
ggsave(paste0(paths$results, "/lin_vs_binom.pdf"), width = w_in, height = 0.3*h_in)

##### Volcano #####

# === Binomial ===

table(results5mC_binom$class)
tmp = results5mC_binom %>% # Restrict to TEs only for this
  dplyr::filter(class %in% c("DNA", "LTR", "RC", "Retroposon", "LINE", "SINE"))

test = wilcox.test(tmp$coef, conf.int = T) # binom or lin makes no diff
med = median(tmp$coef)
pval = ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
my_volcano(tmp, "coef", "pval", "sig", ymax=7)+
  labs(x="Coefficient")+
  ggtitle("TE methylation", subtitle = sprintf("Median coef = %.3f, pvalue %s", med, pval))
ggsave(paste0(paths$results, "/volcano_te_binom.pdf"), width=0.7*w_in, height=0.4*h_in)

# === Linear ===

tmp = results5mC_lin %>% # Restrict to TEs only for this
  dplyr::filter(class %in% c("DNA", "LTR", "RC", "Retroposon", "LINE", "SINE"))

test = wilcox.test(tmp$coef, conf.int = T) # binom or lin makes no diff
med = median(tmp$coef)
pval = ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
my_volcano(tmp, "coef", "pval", "sig", ymax=7)+
  labs(x="Coefficient")+
  ggtitle("TE methylation", subtitle = sprintf("Median coef = %.3f, pvalue %s", med, pval))
ggsave(paste0(paths$results, "/volcano_te_lin.pdf"), width=0.7*w_in, height=0.4*h_in)

#### CLASS ####

# Which classes of TEs are biased towards methyl loss? 

min_n_fam = 10

# === Binomial ===

class_stats = results5mC_binom %>%
  group_by(class) %>%
  summarize(medi_coef=median(coef), n=n())%>%
  dplyr::filter(n>10) %>%
  dplyr::filter(class != "Simple_repeat")
for (i in 1:nrow(class_stats)) {
  test = wilcox.test(results5mC_binom[results5mC_binom$class == class_stats$class[i], "coef"])
  class_stats[i, "sig"] = stars.pval(test$p.value)
}
yrange = quantile(results5mC_binom$coef, c(0.01, 0.99))
p1 = results5mC_binom %>%
  merge(., class_stats, by="class") %>%
  ggplot(., aes(x=class, y=coef, fill=medi_coef))+
  geom_boxplot(outlier.shape = NA)+
  geom_text(data=class_stats, aes(x=class, label=sig), y=yrange[2], size=5) +
  ylim(yrange)+
  theme(axis.title.x = element_blank())+
  scale_fill_gradient2(low="blue", high="red", midpoint = 0)+
  labs(y="Coefficient", fill="Median\ncoef.")+
  ggtitle("Binomial model")

# === Linear ===

class_stats = results5mC_lin %>%
  group_by(class) %>%
  summarize(medi_coef=median(coef), n=n())%>%
  dplyr::filter(n>10) %>%
  dplyr::filter(class != "Simple_repeat")
for (i in 1:nrow(class_stats)) {
  test = wilcox.test(results5mC_lin[results5mC_lin$class == class_stats$class[i], "coef"])
  class_stats[i, "sig"] = stars.pval(test$p.value)
}
yrange = quantile(results5mC_lin$coef, c(0.01, 0.99))
p2 = results5mC_lin %>%
  merge(., class_stats, by="class") %>%
  ggplot(., aes(x=class, y=coef, fill=medi_coef))+
  geom_boxplot(outlier.shape = NA)+
  geom_text(data=class_stats, aes(x=class, label=sig), y=yrange[2], size=5) +
  ylim(yrange)+
  theme(axis.title.x = element_blank())+
  scale_fill_gradient2(low="blue", high="red", midpoint = 0)+
  labs(y="Coefficient", fill="Median\ncoef.")+
  ggtitle("Linear model")

pdf(paste0(paths$results, "/coef_by_class.pdf"), width=0.5*w_in, height=0.3*h_in)
list(p1, p2)
dev.off()

#### LINE1 ####

ps = list()

# === Binomial ===

results_l1 = results5mC_binom %>%
  dplyr::filter(superf == "LINE/L1") %>%
  merge(anno_re2[c("superf", "fam", "mean_length")], ., by=c("superf", "fam")) %>%
  mutate(fam = fct_reorder(fam, mean_length)) 
p1 = ggplot(results_l1, aes(x=as.numeric(fam), y=mean_length))+
  geom_line()+
  theme(axis.title.x=element_blank())+
  scale_x_continuous(expand=expansion(mult=c(0,0)))+
  labs(y="Avg. copy\nlength (bp)")+
  ggtitle("Binomial model")
p2 = ggplot(results_l1, aes(x=fam, y=coef, fill=meanCov))+
  geom_bar(stat="identity")+
  labs(x="LINE1 families sorted by average copy length", y="Coef")+
  theme(axis.text.x = element_text(angle=90, hjust=1))
ps[["binom"]] = p1 / p2 + plot_layout(heights = c(1,4))

# === Linear ===

results_l1 = results5mC_lin %>%
  dplyr::filter(superf == "LINE/L1") %>%
  merge(anno_re2[c("superf", "fam", "mean_length")], ., by=c("superf", "fam")) %>%
  mutate(fam = fct_reorder(fam, mean_length)) 
p1 = ggplot(results_l1, aes(x=as.numeric(fam), y=mean_length))+
  geom_line()+
  theme(axis.title.x=element_blank())+
  scale_x_continuous(expand=expansion(mult=c(0,0)))+
  labs(y="Avg. copy\nlength (bp)")+
  ggtitle("Linear model")
p2 = ggplot(results_l1, aes(x=fam, y=coef, fill=meanCov))+
  geom_bar(stat="identity")+
  labs(x="LINE1 families sorted by average copy length", y="Coef")+
  theme(axis.text.x = element_text(angle=90, hjust=1))
ps[["lin"]] = p1 / p2 + plot_layout(heights = c(1,4))

pdf(paste0(paths$results, "/coef_line1.pdf"), width=2*w_in, height=0.5*h_in)
ps
dev.off()

#### L1Md ####

p1 = results5mC_binom %>%
  dplyr::filter(grepl("L1Md", fam)) %>%
  mutate(type = str_extract(fam, "(L1Md[:alnum:]+)_", group=1))%>%
  mutate(subtype = str_extract(fam, "L1Md[:alnum:]+_(.+)", group=1))%>%
  mutate(SigStar = stars.pval(padj)) %>%
  ggplot(., aes(subtype, coef, label=SigStar))+
  geom_bar(stat="identity")+
  geom_text(vjust=1, size=4)+
  theme(axis.title.x=element_blank())+
  facet_grid(~type, scales = "free_x", space="free_x")+
  labs(y="Coefficient")+
  ggtitle("Binomial model")
p2 = results5mC_lin %>%
  dplyr::filter(grepl("L1Md", fam)) %>%
  mutate(type = str_extract(fam, "(L1Md[:alnum:]+)_", group=1))%>%
  mutate(subtype = str_extract(fam, "L1Md[:alnum:]+_(.+)", group=1))%>%
  mutate(SigStar = stars.pval(padj)) %>%
  ggplot(., aes(subtype, coef, label=SigStar))+
  geom_bar(stat="identity")+
  geom_text(vjust=1, size=4)+
  theme(axis.title.x=element_blank())+
  facet_grid(~type, scales = "free_x", space="free_x")+
  labs(y="Coefficient")+
  ggtitle("Linear model")

pdf(paste0(paths$results, "/coef_l1md.pdf"), width=1*w_in, height=0.3*h_in)
list(p1, p2)
dev.off()

#### SINE ####

p1=results5mC_binom %>%
  dplyr::filter(class=="SINE") %>%
  dplyr::filter(meanCov > 500) %>%
  mutate(type = str_remove(superf, "SINE/"))%>%
  mutate(type = str_replace(type, "-", "\n"))%>%
  mutate(SigStar = stars.pval(padj)) %>%
  ggplot(., aes(fam, coef, label=SigStar))+
  geom_bar(stat="identity")+
  geom_text(vjust=1, size=4)+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1))+
  facet_grid(~type, scales = "free_x", space="free_x")+
  labs(y="Coefficient")+
  ggtitle("Binomial model")
p2=results5mC_lin %>%
  dplyr::filter(class=="SINE") %>%
  dplyr::filter(meanCov > 500) %>%
  mutate(type = str_remove(superf, "SINE/"))%>%
  mutate(type = str_replace(type, "-", "\n"))%>%
  mutate(SigStar = stars.pval(padj)) %>%
  ggplot(., aes(fam, coef, label=SigStar))+
  geom_bar(stat="identity")+
  geom_text(vjust=1, size=4)+
  theme(axis.title.x=element_blank(),
        axis.text.x = element_text(angle=45, hjust = 1))+
  facet_grid(~type, scales = "free_x", space="free_x")+
  labs(y="Coefficient")+
  ggtitle("Linear model")

pdf(paste0(paths$results, "/coef_sine.pdf"), width=1*w_in, height=0.3*h_in)
list(p1, p2)
dev.off()

#### INDIVIDUAL FEATURES ####

all.equal(rownames(Perc), meta$Sample)
which_re = "B2_Mm1a"
which_id = rownames(re_info)[re_info$fam == which_re]
meta %>%
  mutate(re = Perc[, which_id]) %>%
  ggplot(., aes(x=Genotype, y=re))+
  geom_boxplot()+
  geom_jitter(width=0.2)+
  ggtitle(which_re)+
  stat_compare_means(method = "t.test")

#### CHECKPOINT ####

# save.image(file=paste0(paths$results, "/checkpoint.Rdata"))
load(file=paste0(paths$results, "/checkpoint.Rdata"))

#### PAPER FIGURES ####

theme_set(theme_classic())

# Volcano
tmp = results5mC_binom %>% # Restrict to TEs only for this
  dplyr::filter(class %in% c("DNA", "LTR", "RC", "Retroposon", "LINE", "SINE"))
test = wilcox.test(tmp$coef, conf.int = T) # binom or lin makes no diff
med = median(tmp$coef)
pval = ifelse(test$p.value < 0.001, "< 0.001", round(test$p.value, 3))
p1 = my_volcano(tmp, "coef", "pval", "sig", ymax=7)+
  labs(x="Coefficient")+
  ggtitle("TE methylation", subtitle = sprintf("Median coef = %.3f, pvalue %s", med, pval))+
  guides(color="none")+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank())+
  labs(y=expression(-log[10]*"(pval)"))

# TE classes
class_stats = results5mC_binom %>%
  group_by(class) %>%
  summarize(medi_coef=median(coef), n=n())%>%
  dplyr::filter(n>10) %>%
  dplyr::filter(class != "Simple_repeat")
for (i in 1:nrow(class_stats)) {
  test = wilcox.test(results5mC_binom[results5mC_binom$class == class_stats$class[i], "coef"])
  class_stats[i, "sig"] = stars.pval(test$p.value)
}
yrange = quantile(results5mC_binom$coef, c(0.01, 0.99))
p2=results5mC_binom %>%
  merge(., class_stats, by="class") %>%
  ggplot(., aes(x=class, y=coef, fill=medi_coef))+
  geom_boxplot(outlier.shape = NA)+
  geom_text(data=class_stats, aes(x=class, label=sig), y=yrange[2], size=5) +
  geom_hline(yintercept=0, linetype="dashed")+
  ylim(yrange)+
  scale_fill_gradient2(low="blue", high="red", midpoint = 0)+
  labs(y="Coefficient", fill="Median\ncoef.")+
  theme(axis.title.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank())
p1|p2

# L1Md behavior
tmp = results5mC_binom %>%
  dplyr::filter(grepl("L1Md", fam)) %>%
  mutate(type = str_extract(fam, "(L1Md[:alnum:]+)_", group=1))%>%
  mutate(subtype = str_extract(fam, "L1Md[:alnum:]+_(.+)", group=1))%>%
  mutate(SigStar = stars.pval(padj)) %>%
  group_by(type) %>%
  mutate(subtype_rank = (rank(subtype)-1)/7) %>%
  ungroup() %>%
  mutate(base_col = brewer_pal(palette = "Dark2")(length(unique(type)))[match(type, unique(type))]) %>%
  mutate(fill_col = darken(base_col, amount = subtype_rank))
p3 = ggplot(tmp, aes(subtype, coef, label=SigStar))+
  geom_bar(stat="identity", aes(fill = fill_col))+
  geom_text(vjust=1, size=4)+
  facet_grid(cols = vars(type), scales = "free_x", space="free_x")+
  scale_fill_identity()+
  labs(y="Coefficient")+
  theme(
    axis.title.x = element_blank(),
    strip.text = element_text(size = 7),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank())

# SINE behavior
p4 = results5mC_binom %>%
  dplyr::filter(class=="SINE") %>%
  dplyr::filter(meanCov > 500) %>%
  mutate(type = str_remove(superf, "SINE/"))%>%
  mutate(type = str_replace(type, "-", "\n"))%>%
  mutate(type = str_replace(type, "Alu", "B1")) %>%
  mutate(SigStar = stars.pval(padj)) %>%
  ggplot(., aes(fam, coef, label=SigStar, fill=type))+
  geom_bar(stat="identity")+
  geom_text(vjust=1, size=4)+
  facet_grid(~type, scales = "free_x", space="free_x")+
  scale_fill_manual(values=desaturate(brewer.pal(5, "Set1"), 0.2))+
  guides(fill="none")+
  labs(y="Coefficient")+
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle=45, hjust = 1),
    strip.text = element_text(size = 7),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank())
p3/p4

p=(p1|p2)/p3/p4
ggsave(plot = p, paste0(paths$results, "/figure_nano.pdf"), width = w_in, height=0.8*h_in)
