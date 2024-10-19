library(tidyverse)
library(edgeR)
library(clusterProfiler)
library(ReactomePA)
library(BiocParallel)
library(rasterpdf)
library(patchwork)
library(enrichplot)
library(GOSemSim)

setwd("/scratch/fmorandi/internal/John/cGAS_KO/PAPER/RNA2")

paths = list()
paths$data = "./pipeline_out"
paths$results = "./results"
paths$tables = paste0(paths$results, "/tables")
paths$paper_figs = paste0(paths$results, "/prelim_paper_figs")

dir.create(paths$results, showWarnings = F)
dir.create(paths$tables, showWarnings = F)

#### FIGURE STORAGE (added for report) ####

figs = list()

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701

#### FUNCTIONS ####

my_volcano = function(table, logFC_col, pval_col, sig_col, clip_axes=T, clip_quantile=0.999) {
  if (clip_axes) {
    xmax = quantile(abs(table[, logFC_col]), clip_quantile) # Cutoff changed for report, was 0.999
    ymax = quantile(-log10(table[, pval_col]), clip_quantile)
  } else {
    xmax = max(abs(table[, logFC_col]))
    ymax = max(-log10(table[, pval_col]))
  }
  tmp = table %>%
    group_by_at(sig_col, .drop=F) %>%
    dplyr::summarize(n=n()) %>%
    mutate(x = c(-xmax*0.9, 0, xmax*0.9)) %>%
    mutate(hjust = c(0, 0.5, 1))
  p = ggplot(table, aes(x=.data[[logFC_col]], y=-log10(.data[[pval_col]]), color=.data[[sig_col]]))+
    geom_point(size=0.1)+
    lims(x=c(-xmax, xmax), y=c(0,ymax))+
    geom_text(data=tmp, aes(x=x, y=ymax*0.95, label=n))+
    scale_color_manual(values=c("#5555cc", "#999999", "#cc5555"), drop=F)+
    labs(x="logFC", y="-log10(pval)", color="Significance")
  return(p)
}

my_gseGO = function(table, logFC_column, orgdb, ontology="BP") {
  options(warn = 1)
  # Make gene list
  gene_list = table[[logFC_column]]
  names(gene_list) = table$Geneid
  gene_list = sort(gene_list, decreasing = TRUE)
  # Add a tiny amount to break ties
  # gene_list[duplicated(gene_list)] = gene_list[duplicated(gene_list)] - 1e-5
  # Run GSEA
  set.seed(1337)
  res = gseGO(
    geneList = gene_list,
    ont = ontology,
    keyType = "ENSEMBL",
    verbose = TRUE,
    OrgDb = orgdb,
    pvalueCutoff = 1.1,
    BPPARAM = SerialParam())
  return(res)
}

my_gseReactome = function(table, logFC_column, org) {
  options(warn = 1)
  # Make gene list
  gene_list = table[[logFC_column]]
  names(gene_list) = table$entrezgene_id
  gene_list = sort(gene_list, decreasing = TRUE)
  # Add a tiny amount to break ties
  # gene_list[duplicated(gene_list)] = gene_list[duplicated(gene_list)] - 1e-5
  # Run GSEA
  set.seed(1337)
  res = gsePathway(
    geneList=gene_list,
    verbose = TRUE,
    organism = org,
    pvalueCutoff = 1.1,
    BPPARAM = SerialParam())
  return(res)
}

convert_ids_in_string = function(df, geneid, symbol) {
  conversions = symbol
  names(conversions) = geneid
  conversions = conversions[!is.na(names(conversions))]
  for (i in 1:nrow(df)) {
    this_str = df$core_enrichment[i]
    old_ids = unlist(strsplit(this_str, "/"))
    new_ids = conversions[old_ids]
    df[i, "core_enrichment"] = paste(new_ids, collapse="/")
  }
  return(df)
}

#### LOAD DATA ####

load(paste0(paths$results, "/prepro.Rdata"))
meta$Genotype = factor(meta$Genotype, c("WT", "cGASKO"))

#### KO VALIDATION ####

all.equal(colnames(norm_ge), meta$SampleID)

meta %>%
  mutate(Cgas = unlist(norm_ge["Cgas", ])) %>%
  ggplot(., aes(Genotype, Cgas))+
  geom_boxplot()+
  theme(axis.title.x = element_blank())
ggsave(paste0(paths$results, "/ko_validation.png"), width = w*1, height=h*0.3, units = "mm")

#### PCA ####

##### Genes #####

pca = prcomp(t(log1p(norm_ge)), center=T, scale.=T)
pca = merge(meta, pca$x[, c("PC1", "PC2")], by.x="SampleID", by.y=0)
ps = list()
for (var in c("Genotype", "Sex")) { 
  ps[[var]] = ggplot(pca, aes(x=PC1, y=PC2, color=.data[[var]]))+
    geom_point()
}

ps = align_patches(ps)
pdf(paste0(paths$results, "/pca_ge.pdf"), width=1*w_in, height=0.5*h_in)
ps
dev.off()

##### <Extra for report> #####

all.equal(meta$SampleID, colnames(norm_ge))
tmp = norm_ge[rowSums(norm_ge) > 0, ] # Some columns are all-zero so remove
pca = prcomp(t(log1p(tmp)), center=T, scale.=T)
pca = merge(meta, pca$x[, c("PC1", "PC2")], by.x="SampleID", by.y=0)
pca = pca %>%
  mutate(Genotype2 = str_remove(Genotype, "cGAS")) %>%
  mutate(Genotype2 = factor(Genotype2, levels=c("WT", "KO")))
figs$pca_combined = ggplot(pca, aes(x=PC1, y=PC2, color=Genotype2, shape=Sex))+
  geom_point()+
  labs(color = "Genotype")
figs$pca_sex = ggplot(pca, aes(x=PC1, y=PC2, color=Sex))+
  geom_point()+
  labs(color = "Sex")
figs$pca_genotype = ggplot(pca, aes(x=PC1, y=PC2, color=Genotype2))+
  geom_point()+
  labs(color = "Genotype")

##### Repeats #####

pca = prcomp(t(log1p(norm_re)), center=T, scale.=T)
pca = merge(meta, pca$x[, c("PC1", "PC2")], by.x="SampleID", by.y=0)
ps = list()
for (var in c("Genotype", "Sex")) { 
  ps[[var]] = ggplot(pca, aes(x=PC1, y=PC2, color=.data[[var]]))+
    geom_point()
}

ps = align_patches(ps)
pdf(paste0(paths$results, "/pca_re.pdf"), width=1*w_in, height=0.5*h_in)
ps
dev.off()

#### DIFFERENTIAL EXPRESSION ####

all.equal(colnames(counts_ge), meta$SampleID)
all.equal(colnames(counts_re), meta$SampleID)
counts = rbind(counts_ge, counts_re) # Rejoin ge an re

##### Together #####

# Based on comparison with qPCR, FISH, IF data, this method does not appear appropriate

# dge = DGEList(counts, samples = meta)
# design = model.matrix(~0+Genotype+Sex, data=dge$samples)
# contrasts = makeContrasts(
#   cGASKOvsWT=GenotypecGASKO-GenotypeWT,
#   levels = design)
# 
# # Filter further
# keep = filterByExpr(dge, design=design)
# perc_kept = 100 * colSums(counts[keep, ]) / colSums(counts)
# hist(perc_kept, breaks=100)
# dge = dge[keep, ]
# dim(dge)
# 
# # TMM normalization
# dge = calcNormFactors(dge)
# ggplot(dge$samples, aes(x=norm.factors))+
#   geom_histogram()
# 
# dge = estimateDisp(dge, design)
# fit = glmFit(dge, design)
# 
# resultsDE = list()
# for (cont in colnames(contrasts)) {
#   print(cont)
#   resultsDE[[cont]] = as.data.frame(glmLRT(fit, contrast=contrasts[, cont])) %>%
#     dplyr::select(-logCPM, -LR) %>%
#     dplyr::rename(pval = PValue) %>%
#     mutate(padj = p.adjust(pval, method="BH")) %>%
#     rename_with(function(x) paste0(x, "_", cont)) %>%
#     rownames_to_column("geneID")
# }
# resultsDE = Reduce(function(x, y) merge(x, y, by="geneID"), resultsDE)
# 
# # Split results into genes and res
# resultsDE_ge = merge(ginfo, resultsDE, by.x=0, by.y="geneID") %>%
#   dplyr::rename("geneID"="Row.names")
# resultsDE_re = merge(rinfo, resultsDE, by.x=0, by.y="geneID") %>%
#   dplyr::rename("geneID"="Row.names")

##### Genes #####

dge = DGEList(counts_ge, samples = meta) 
design = model.matrix(~0+Genotype+Sex, data=dge$samples)
contrasts = makeContrasts(
  cGASKOvsWT=GenotypecGASKO-GenotypeWT,
  levels = design)

# Filter further
keep = filterByExpr(dge, design=design)
perc_kept = 100 * colSums(counts_ge[keep, ]) / colSums(counts_ge)
hist(perc_kept, breaks=100)
dge = dge[keep, ]
dim(dge)

# TMM normalization
dge = calcNormFactors(dge)
ggplot(dge$samples, aes(x=norm.factors))+
  geom_histogram()

dge = estimateDisp(dge, design)
fit = glmFit(dge, design)

resultsDE_ge = list()
for (cont in colnames(contrasts)) {
  print(cont)
  resultsDE_ge[[cont]] = as.data.frame(glmLRT(fit, contrast=contrasts[, cont])) %>%
    dplyr::select(-logCPM, -LR) %>%
    dplyr::rename(pval = PValue) %>%
    mutate(padj = p.adjust(pval, method="BH")) %>%
    rename_with(function(x) paste0(x, "_", cont)) %>%
    rownames_to_column("geneID")
}
resultsDE_ge = Reduce(function(x, y) merge(x, y, by="geneID"), resultsDE_ge)
resultsDE_ge = merge(ginfo, resultsDE_ge, by.x=0, by.y="geneID") %>%
  dplyr::rename("geneID"="Row.names")

##### Repeats #####

dge = DGEList(counts_re, samples = meta, lib.size = meta$annotated) # <!>

# Filter further
keep = filterByExpr(dge, design=design)
perc_kept = 100 * colSums(counts_re[keep, ]) / colSums(counts_re)
hist(perc_kept, breaks=100)
dge = dge[keep, ]
dim(dge)

# TMM normalization
dge = calcNormFactors(dge) #<!>
ggplot(dge$samples, aes(x=norm.factors))+
  geom_histogram()

dge = estimateDisp(dge, design)
fit = glmFit(dge, design)

resultsDE_re = list()
for (cont in colnames(contrasts)) {
  print(cont)
  resultsDE_re[[cont]] = as.data.frame(glmLRT(fit, contrast=contrasts[, cont])) %>%
    dplyr::select(-logCPM, -LR) %>%
    dplyr::rename(pval = PValue) %>%
    mutate(padj = p.adjust(pval, method="BH")) %>%
    rename_with(function(x) paste0(x, "_", cont)) %>%
    rownames_to_column("repID")
}
resultsDE_re = Reduce(function(x, y) merge(x, y, by="repID"), resultsDE_re)
resultsDE_re = merge(rinfo, resultsDE_re, by.x=0, by.y="repID") %>%
  dplyr::rename("repID"="Row.names")

#### DEFINE SIGNIFICANCE ####

max_padj = 0.05
min_logFC = 0.5
for (cont in colnames(contrasts)) {
  logFC = resultsDE_ge[paste0("logFC_", cont)]
  padj = resultsDE_ge[paste0("padj_", cont)]
  sig = abs(logFC) > min_logFC & padj < max_padj
  resultsDE_ge[paste0("sig_", cont)] = "Not Sig."
  resultsDE_ge[sig & logFC > 0, paste0("sig_", cont)] = "Sig. Up"
  resultsDE_ge[sig & logFC < 0, paste0("sig_", cont)] = "Sig. Down"
  resultsDE_ge[paste0("sig_", cont)] = factor(
    resultsDE_ge[, paste0("sig_", cont)],
    levels = c("Sig. Down", "Not Sig.", "Sig. Up"))
  resultsDE_ge = relocate(resultsDE_ge, paste0("sig_", cont), .after=paste0("padj_", cont))
}

for (cont in colnames(contrasts)) {
  logFC = resultsDE_re[paste0("logFC_", cont)]
  padj = resultsDE_re[paste0("padj_", cont)]
  sig = abs(logFC) > min_logFC & padj < max_padj
  resultsDE_re[paste0("sig_", cont)] = "Not Sig."
  resultsDE_re[sig & logFC > 0, paste0("sig_", cont)] = "Sig. Up"
  resultsDE_re[sig & logFC < 0, paste0("sig_", cont)] = "Sig. Down"
  resultsDE_re[paste0("sig_", cont)] = factor(
    resultsDE_re[, paste0("sig_", cont)],
    levels = c("Sig. Down", "Not Sig.", "Sig. Up"))
  resultsDE_re = relocate(resultsDE_re, paste0("sig_", cont), .after=paste0("padj_", cont))
}

#### VOLCANOS ####

##### Genes #####

ps = list()
for (cont in colnames(contrasts)) {
  ps[[cont]] = my_volcano(
    resultsDE_ge, paste0("logFC_", cont), paste0("pval_", cont), paste0("sig_", cont), clip_axes=F)+
    ggtitle(cont, subtitle=sprintf("padj < %.2f, |logFC| > %.2f", max_padj, min_logFC))
}

ps = align_patches(ps)
raster_pdf(paste0(paths$results, "/volcanos_ge.pdf"), width=1*w_in, height=0.5*h_in, res=300)
ps
dev.off()

##### <Extra for report> #####

figs$volcano = my_volcano(
  resultsDE_ge, paste0("logFC_", "cGASKOvsWT"), paste0("pval_", "cGASKOvsWT"), paste0("sig_", "cGASKOvsWT"), clip_axes=T, clip_quantile=0.999)+
  ggtitle("cGASKOvsWT", subtitle=sprintf("padj < %.2f, |logFC| > %.2f", max_padj, min_logFC))

##### Repeats #####

pdf(paste0(paths$results, "/volcanos_re.pdf"), width=2*w_in, height=0.5*h_in)
for (cont in colnames(contrasts)) {
  ps = list()
  # bps = list()
  for (re_class in c("LINE", "SINE", "LTR", "DNA")) {
    tmp = resultsDE_re %>%
      dplyr::filter(class == re_class)
    pval = round(wilcox.test(tmp[, paste0("logFC_", cont)])$p.value, 3)
    pval = ifelse(pval == 0, "<0.001", paste0("=", pval))
    median_logFC = median(tmp[, paste0("logFC_", cont)])
    ps[[re_class]] = my_volcano(
      tmp, paste0("logFC_", cont), paste0("pval_", cont), paste0("sig_", cont), clip_axes=F)+
      ggtitle(re_class, subtitle=sprintf("Median logFC=%f, p%s", median_logFC, pval))
  }
  p=(wrap_plots(ps, nrow=1)&labs(color=sprintf("Significance\npadj < %.2f, |logFC| > %.2f", max_padj, min_logFC)))+
    plot_layout(guides="collect")+
    plot_annotation(title=cont)
  print(p)
}
dev.off()

#### FRIR (added for report) ####

# For ATAC this was based on the number of reads in repeats over the number of total reads counted
# When using a GTF file containing both peaks and repeats
# In other words, how many of the reads found over features are over REs
# The equivalent of that is the following, with annotated as denominator

all.equal(colnames(counts_re), meta$SampleID)
meta$RiR = colSums(counts_re)

ggplot(meta, aes(Genotype, RiR/input_reads))+
  geom_boxplot()
ggplot(meta, aes(Genotype, RiR/uniquely_mapped))+
  geom_boxplot()
figs$frir = meta %>%
  mutate(Genotype = fct_recode(Genotype, "KO" = "cGASKO")) %>%
  mutate(Genotype = factor(Genotype, levels=c("WT", "KO"))) %>%
  ggplot(., aes(Genotype, RiR/annotated))+
  geom_boxplot()+
  labs(y="Fraction of reads in repeats")+
  stat_compare_means(label.y = 0.11)+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))

# TEs or REs?
rinfo$selfish = rinfo$class %in% c("LINE", "SINE", "LTR", "DNA", "DNA?", "Retroposon", "RC")
meta$RiT = colSums(counts_re[rinfo$selfish, ])
ggplot(meta, aes(Genotype, RiT/annotated))+
  geom_boxplot()+
  labs(y="Fraction of reads in repeats")

#### LINE LENGTH EFFECT ####

resultsL1 = resultsDE_re %>%
  dplyr::filter(superf=="LINE/L1") %>%
  mutate(young = grepl("L1Md", fam)) %>%
  mutate(fam=fct_reorder(fam, young))
hist(resultsL1$meanExpr, breaks=50) # Mean CPM (not log)
# resultsL1 = subset(resultsL1, meanExpr > 3)

ps = list()
for (cont in colnames(contrasts)) {
  ps[[cont]]=ggplot(resultsL1, aes(x=fam, y=.data[[paste0("logFC_", cont)]], fill=young))+
    geom_bar(stat="identity")+
    theme(axis.text.x = element_text(angle=90, hjust=1, size=7))
}

pdf(paste0(paths$results, "/LINE_length_effect.pdf"), width=2*w_in, height=0.5*h_in)
ps
dev.off()

#### ZOOM ON L1Md (added for report) ####

# figs$l1md_zoom = 
figs$l1md_zoom = resultsL1 %>%
  dplyr::filter(grepl("L1Md", fam)) %>%
  mutate(type = str_extract(fam, "(L1Md[:alnum:]+)_", group=1))%>%
  mutate(subtype = str_extract(fam, "L1Md[:alnum:]+_(.+)", group=1))%>%
  mutate(SigStar = stars.pval(pval_cGASKOvsWT)) %>%
  ggplot(., aes(subtype, logFC_cGASKOvsWT, label=SigStar))+
  geom_bar(stat="identity")+
  geom_text(vjust=0, size=6)+
  theme(axis.title.x=element_blank())+
  facet_grid(~type, scales = "free_x", space="free_x")
figs$l1mda_zoom = resultsL1 %>%
  mutate(L1MdA = grepl("L1MdA", fam)) %>%
  dplyr::filter(L1MdA) %>%
  ggplot(., aes(fam, logFC_cGASKOvsWT))+
  geom_bar(stat="identity")+
  theme(axis.title.x=element_blank())
figs$l1mda123_zoom = resultsL1 %>%
  dplyr::filter(fam %in% c("L1MdA_I", "L1MdA_II", "L1MdA_III")) %>%
  mutate(SigStar = stars.pval(pval_cGASKOvsWT)) %>%
  ggplot(., aes(fam, logFC_cGASKOvsWT, label=SigStar))+
  geom_bar(stat="identity")+
  geom_text(vjust=-0.5, size=8)+
  theme(axis.title.x=element_blank())+
  ylim(0, 0.30)

#### GSEA ####

# gsea = list()
# 
# for (ont in c("BP", "MF", "CC")) {
#   for (cont in colnames(contrasts)) {
#     gsea[[paste0("GO_", ont)]][[cont]] = my_gseGO(resultsDE_ge, paste0("logFC_", cont), org.Mm.eg.db, ont)
#   }
# }
# for (cont in colnames(contrasts)) {
#   gsea$Reactome[[cont]] = my_gseReactome(resultsDE_ge, paste0("logFC_", cont), "mouse")
# }

# save(gsea, file=paste0(paths$results, "/gsea.Rdata"))
load(paste0(paths$results, "/gsea.Rdata"))

#### GSEA DOTPLOTS ####

for (ont in names(gsea)) {
  for (cont in colnames(contrasts)) {
    gsea[[ont]][[cont]] = data.frame(gsea[[ont]][[cont]])
    gsea[[ont]][[cont]]$Comparison = cont
    if (ont == "Reactome") {
      gsea[[ont]][[cont]] = convert_ids_in_string(gsea[[ont]][[cont]], ginfo$entrezgene_id, ginfo$mgi_symbol)
    } else {
      gsea[[ont]][[cont]] = convert_ids_in_string(gsea[[ont]][[cont]], ginfo$Geneid, ginfo$mgi_symbol)
    }
  }
}

max_padj_gsea = 0.05

for (ont in names(gsea)) {
  ps = list()
  for (cont in colnames(contrasts)) {
    tmp = gsea[[ont]][[cont]] %>%
      dplyr::filter(p.adjust < max_padj_gsea) %>%
      mutate(Description = fct_reorder(Description, abs(NES))) %>%
      mutate(Sign = ifelse(sign(NES) > 0, "Activated", "Suppressed")) %>%
      group_by(Sign) %>%
      slice_max(abs(NES), n=20)
    ps[[cont]] = ggplot(tmp, aes(x=abs(NES), y=Description, color=p.adjust))+
      geom_point(size=5)+
      scale_color_gradient(low="red", high="blue")+
      ggtitle(paste(ont, "GSEA:", cont))+
      theme(axis.text.y = element_text(size = 10))
    if (nrow(tmp) > 0) {
      ps[[cont]] =  ps[[cont]] + facet_wrap(~Sign)
    }
  }
  ps = align_patches(ps)
  pdf(paste0(paths$results, "/gsea_", ont, ".pdf"), width=w_in*2, height=h_in)
  print(ps)
  dev.off()
}

View(gsea$Reactome$cGASKOvsWT)

#### CONDENSED GSEA (added for report) ####


##### Biological process #####

desc_len = 70
nterms = 15

# # Both signs, by NES
# figs$gsea_sortNES_both15 = gsea$GO_BP$cGASKOvsWT %>%
#   dplyr::filter(p.adjust < 0.05) %>%
#   mutate(Description = fct_reorder(Description, abs(NES))) %>%
#   mutate(Sign = ifelse(sign(NES) > 0, "Activated", "Suppressed")) %>%
#   group_by(Sign) %>%
#   slice_max(abs(NES), n=15) %>%
#   ungroup() %>%
#   mutate(Description = if_else(str_length(Description) > 60, str_c(str_sub(Description, 1, 57), "..."), as.character(Description))) %>%
#   mutate(Description = fct_reorder(Description, NES)) %>%
#   ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
#   geom_point(pch = 21)+
#   scale_fill_gradient2(low="blue", high="red")+
#   theme(axis.title = element_blank(),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank())
# Both signs, by pval
figs$gsea_sortPval_both15=gsea$GO_BP$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  mutate(Sign = ifelse(sign(NES) > 0, "Activated", "Suppressed")) %>%
  group_by(Sign) %>%
  arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  slice_head(n=nterms) %>%
  ungroup() %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, sign(NES)*-log(p.adjust))) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
# Upreg, by pval
figs$gsea_sortPval_up15 = gsea$GO_BP$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(NES > 0) %>%
  arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  slice_head(n=nterms) %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, -log(p.adjust))) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
# Downreg, by pval
figs$gsea_sortPval_down15 = gsea$GO_BP$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(NES < 0) %>%
  arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  slice_head(n=nterms) %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, log(p.adjust))) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

# === 10 terms each ===

nterms = 10
figs$gsea_sortPval_both10=gsea$GO_BP$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  mutate(Sign = ifelse(sign(NES) > 0, "Activated", "Suppressed")) %>%
  group_by(Sign) %>%
  arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  slice_head(n=nterms) %>%
  ungroup() %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, sign(NES)*-log(p.adjust))) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
# Upreg, by pval
figs$gsea_sortPval_up10 = gsea$GO_BP$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(NES > 0) %>%
  arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  slice_head(n=nterms) %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, -log(p.adjust))) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
# Downreg, by pval
figs$gsea_sortPval_down10 = gsea$GO_BP$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(NES < 0) %>%
  arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  slice_head(n=nterms) %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, log(p.adjust))) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())

# === 10 terms each ===

nterms = 8
figs$gsea_sortPval_both8=gsea$GO_BP$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  mutate(Sign = ifelse(sign(NES) > 0, "Activated", "Suppressed")) %>%
  group_by(Sign) %>%
  arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  slice_head(n=nterms) %>%
  ungroup() %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, sign(NES)*-log(p.adjust))) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
# Upreg, by pval
figs$gsea_sortPval_up8 = gsea$GO_BP$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(NES > 0) %>%
  arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  slice_head(n=nterms) %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, -log(p.adjust))) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
# Downreg, by pval
figs$gsea_sortPval_down8 = gsea$GO_BP$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  dplyr::filter(NES < 0) %>%
  arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  slice_head(n=nterms) %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, log(p.adjust))) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())


##### Cellular component #####

# Max text shorter here
desc_len = 24
nterms = 5
figs$gsea_CC_sortPval = gsea$GO_CC$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  mutate(Sign = ifelse(sign(NES) > 0, "Activated", "Suppressed")) %>%
  group_by(Sign) %>%
  arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  slice_head(n=nterms) %>%
  ungroup() %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, sign(NES)*-log(p.adjust))) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())
  

#### EXTRA GSEA (added for report) ####

my_read_gmt = function(file) {
  lines = readLines(file)
  lines = strsplit(lines, "\t")
  ids = vapply(lines, function(y) y[1], character(1))
  descs = vapply(lines, function(y) y[2], character(1))
  genes = lapply(lines, "[", -c(1:2))
  names(genes) = ids
  names(descs) = ids
  gmt = stack(genes)
  gmt$desc = descs[gmt$ind]
  colnames(gmt) = c("gene", "term", "desc")
  return(gmt[, c("term", "gene", "desc")])
}

my_gse_custom = function(table, logFC_column, sym_col, gmt) { # Have to specify symbol col because humans 
  options(warn = 1)
  # Make gene list
  gene_list = table[[logFC_column]]
  names(gene_list) = table[[sym_col]]
  gene_list = sort(gene_list, decreasing = TRUE)
  # Run GSEA
  set.seed(1337)
  res = GSEA(
    geneList=gene_list,
    TERM2GENE = gmt,
    verbose = TRUE,
    pvalueCutoff = 1.1,
    minGSSize = 5, # Important or some small custom lists were not shown
    maxGSSize = 1400,
    BPPARAM = SerialParam())
  return(res)
}

gmt = my_read_gmt("./custom_gsets/SAUL_SEN_MAYO.v2023.2.Mm.edited.gmt")
sen_mayo = my_gse_custom(resultsDE_ge, "logFC_cGASKOvsWT", "geneID", gmt) %>%
  as.data.frame()
tmp = resultsDE_ge %>% # There are inflammatory genes but not necessarily senescence related
  dplyr::filter(padj_cGASKOvsWT < 0.05)

figs$sen_mayo = sen_mayo %>%
  mutate(Description = str_remove(Description, "SenMayo_")) %>%
  mutate(Description = fct_reorder(Description, NES)) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(pvalue), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank())+
  ggtitle("SenMayo")

#### SAVE ####

write.csv(resultsDE_ge, paste0(paths$tables, "/resultsDE.csv"))
write.csv(resultsDE_re, paste0(paths$tables, "/resultsDE_re.csv"))

#### SAVE FIG OBJ (added for report) ####

# save(figs, file=paste0(paths$results, "/figs_for_report.Rdata"))
load(paste0(paths$results, "/figs_for_report.Rdata"))

#### FIGURE ASSEMBLY (added for report) ####

fformat = ".png"
fformat = ".pdf"

figs$pca_combined #+ theme(legend.position = "top")
ggsave(paste0(paths$paper_figs, "/pca_combined", fformat), width = w*0.5, height=h*0.3, units = "mm")
figs$pca_genotype
ggsave(paste0(paths$paper_figs, "/pca_genotype", fformat), width = w*0.5, height=h*0.3, units = "mm")
figs$pca_sex
ggsave(paste0(paths$paper_figs, "/pca_sex", fformat), width = w*0.5, height=h*0.3, units = "mm")

figs$volcano + 
  guides(color="none") + 
  labs(title=NULL, subtitle=NULL)
ggsave(paste0(paths$paper_figs, "/volcano", fformat), width = w*0.5, height=h*0.35, units = "mm")

figs$frir
ggsave(paste0(paths$paper_figs, "/frir", fformat), width = w*0.5, height=h*0.4, units = "mm")

figs$sen_mayo
ggsave(paste0(paths$paper_figs, "/gsea_SenMayo", fformat), width = w*0.6, height=h*0.4, units = "mm")


figs$l1md_zoom+
  scale_y_continuous(expand = expansion(mult=c(0.05, 0.2)))+
  theme(strip.text = element_text(size = 7))
ggsave(paste0(paths$paper_figs, "/zoom_L1Md", fformat), width = w*1, height=h*0.3, units = "mm")
figs$l1mda_zoom+
  theme(axis.text.x = element_text(angle=90, hjust = 1))
ggsave(paste0(paths$paper_figs, "/zoom_L1MdA", fformat), width = w*0.5, height=h*0.4, units = "mm")

figs$l1mda123_zoom
ggsave(paste0(paths$paper_figs, "/zoom_L1MdA123", fformat), width = w*0.5, height=h*0.35, units = "mm")

for (plt in names(figs)[grepl("gsea_sort",names(figs))]) {
  this_h = ifelse(grepl("both", plt), 0.68, 0.34) * h
  this_h = ifelse(grepl("10", plt), this_h * 0.66, this_h)
  this_h = ifelse(grepl("8", plt), this_h * 0.53, this_h)
  p = figs[[plt]]+
    theme(panel.grid.major.x = element_blank(), 
          panel.grid.minor.x = element_blank())
  ggsave(plot = p, paste0(paths$paper_figs, "/", plt, fformat), width = w*1, height=this_h, units = "mm")
}

figs$gsea_CC_sortPval+
  theme(panel.grid.major.x = element_blank(), 
        panel.grid.minor.x = element_blank())
ggsave(paste0(paths$paper_figs, "/gsea_CC_sortPval", fformat), width = w*0.55, height=0.3*h, units = "mm")
