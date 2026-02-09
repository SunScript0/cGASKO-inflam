library(tidyverse)
library(edgeR)
library(clusterProfiler)
library(ReactomePA)
library(BiocParallel)
library(rasterpdf)
library(patchwork)
library(enrichplot)
library(GOSemSim)
library(scales)
library(colorspace)

### Notes ###
# Gotta compare the new de results. They seem very similar but fewer are sig. Maybe filtering of genes changes slightly

setwd("/scratch/fmorandi/internal/John/cGAS_KO/PAPER_CODE/RNA")

paths = list()
paths$data = "./pipeline_out"
paths$results = "./results"
paths$tables = paste0(paths$results, "/tables")
paths$paper_figs = paste0(paths$results, "/prelim_paper_figs")

dir.create(paths$results, showWarnings = F)
dir.create(paths$tables, showWarnings = F)
dir.create(paths$paper_figs, showWarnings = F)

#### PLOTTING SETTINGS ####

w_in = 7.5 # 8.5 without margin
h_in = 10 # 11 without margin

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
ggsave(paste0(paths$results, "/ko_validation.pdf"), width = w_in*0.3, height=h_in*0.3)

#### PCA ####

pca = prcomp(t(log1p(norm_ge)), center=T, scale.=T)
pca = merge(meta, pca$x[, c("PC1", "PC2")], by.x="SampleID", by.y=0)
ps = list()
for (var in c("Genotype", "Sex")) { 
  ps[[var]] = ggplot(pca, aes(x=PC1, y=PC2, color=.data[[var]]))+
    geom_point()
}

ps = align_patches(ps)
pdf(paste0(paths$results, "/pca_ge.pdf"), width=0.6*w_in, height=0.3*h_in)
ps
dev.off()

#### DIFFERENTIAL EXPRESSION ####

all.equal(colnames(counts_ge), meta$SampleID)
all.equal(colnames(counts_re), meta$SampleID)

design = model.matrix(~0+Genotype+Sex, data=meta)
contrasts = makeContrasts(
  cGASKOvsWT=GenotypecGASKO-GenotypeWT,
  levels = design)

##### Genes #####

dge = DGEList(counts_ge, samples = meta) 

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

#### TRANSPOSONS ####

##### Frir #####

# For ATAC this was based on the number of reads in repeats over the number of total reads counted
# When using a GTF file containing both peaks and repeats
# In other words, how many of the reads found over features are over REs
# The equivalent of that is the following, with annotated as denominator

all.equal(colnames(counts_re), meta$SampleID)
meta$RiR = colSums(counts_re)

# Reads in repears in general
p1 = ggplot(meta, aes(Genotype, RiR/input_reads, fill=Genotype))+
  geom_boxplot()
p2 = ggplot(meta, aes(Genotype, RiR/uniquely_mapped, fill=Genotype))+
  geom_boxplot()
p3 = ggplot(meta, aes(Genotype, RiR/annotated, fill=Genotype))+
  geom_boxplot()

# Reads in transposons
rinfo$selfish = rinfo$class %in% c("LINE", "SINE", "LTR", "DNA", "DNA?", "Retroposon", "RC")
meta$RiT = colSums(counts_re[rinfo$selfish, ])
p4 = ggplot(meta, aes(Genotype, RiT/input_reads, fill=Genotype))+
  geom_boxplot()
p5 = ggplot(meta, aes(Genotype, RiT/uniquely_mapped, fill=Genotype))+
  geom_boxplot()
p6 = ggplot(meta, aes(Genotype, RiT/annotated, fill=Genotype))+
  geom_boxplot()

((p1|p2|p3)/(p4|p5|p6)) +
  plot_layout(guides="collect")&
  theme(axis.title.x = element_blank())
ggsave(paste0(paths$results, "/frir.pdf"), width = 0.85*w_in, height = 0.5*h_in)

##### Line length effect #####

resultsL1 = resultsDE_re %>%
  dplyr::filter(superf=="LINE/L1") %>%
  mutate(young = grepl("L1Md", fam)) %>%
  mutate(fam=fct_reorder(fam, young))
ggplot(resultsL1, aes(meanExpr))+
  geom_histogram(bins=100)+
  coord_cartesian(xlim=c(0,200))
resultsL1 = subset(resultsL1, meanExpr > 5)

ps = list()
for (cont in colnames(contrasts)) {
  ps[[cont]]=ggplot(resultsL1, aes(x=fam, y=.data[[paste0("logFC_", cont)]], fill=young))+
    geom_bar(stat="identity")+
    theme(axis.text.x = element_text(angle=90, hjust=1, size=7))
}

pdf(paste0(paths$results, "/LINE_length_effect.pdf"), width=1.6*w_in, height=0.5*h_in)
ps
dev.off()

##### L1Md #####

resultsL1 %>%
  dplyr::filter(grepl("L1Md", fam)) %>%
  mutate(type = str_extract(fam, "(L1Md[:alnum:]+)_", group=1))%>%
  mutate(subtype = str_extract(fam, "L1Md[:alnum:]+_(.+)", group=1))%>%
  mutate(SigStar = stars.pval(pval_cGASKOvsWT)) %>%
  ggplot(., aes(subtype, logFC_cGASKOvsWT, label=SigStar))+
  geom_bar(stat="identity")+
  geom_text(vjust=0, size=6)+
  theme(axis.title.x=element_blank())+
  facet_grid(~type, scales = "free_x", space="free_x")+
  scale_y_continuous(expand = expansion(mult=c(0.05, 0.2)))+
  theme(strip.text = element_text(size = 7))
ggsave(paste0(paths$results, "/zoom_L1Md.pdf"), width = w_in*1, height=h_in*0.3)

#### GSEA ####

##### Run #####

# gsea = list()
# 
# # Gene ontology
# for (ont in c("BP", "MF", "CC")) {
#   for (cont in colnames(contrasts)) {
#     gsea[[paste0("GO_", ont)]][[cont]] = my_gseGO(resultsDE_ge, paste0("logFC_", cont), org.Mm.eg.db, ont)
#   }
# }
# 
# # Reactome
# for (cont in colnames(contrasts)) {
#   gsea$Reactome[[cont]] = my_gseReactome(resultsDE_ge, paste0("logFC_", cont), "mouse")
# }
# 
# # Sen Mayo
# gmt = my_read_gmt("./custom_gsets/SAUL_SEN_MAYO.v2023.2.Mm.edited.gmt")
# for (cont in colnames(contrasts)) {
#   gsea$SenMayo[[cont]] = my_gse_custom(resultsDE_ge, paste0("logFC_", cont), "geneID", gmt)
# }
# 
# save(gsea, file=paste0(paths$results, "/gsea.Rdata"))
load(paste0(paths$results, "/gsea.Rdata"))

##### Dotplots #####

for (ont in names(gsea)) {
  for (cont in colnames(contrasts)) {
    gsea[[ont]][[cont]] = data.frame(gsea[[ont]][[cont]])
    gsea[[ont]][[cont]]$Comparison = cont
    if (ont == "Reactome") {
      gsea[[ont]][[cont]] = convert_ids_in_string(gsea[[ont]][[cont]], ginfo$entrezgene_id, ginfo$mgi_symbol)
    } else if (ont == "SenMayo") {
      next # Ids are already good
    }else {
      gsea[[ont]][[cont]] = convert_ids_in_string(gsea[[ont]][[cont]], ginfo$Geneid, ginfo$mgi_symbol)
    }
  }
}

max_padj_gsea = 0.05

# For the larger dbs
for (ont in c("GO_BP", "GO_MF", "GO_CC", "Reactome")) {
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

gsea$SenMayo$cGASKOvsWT %>%
  mutate(Description = str_remove(Description, "SenMayo_")) %>%
  mutate(Description = fct_reorder(Description, NES)) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(pvalue), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+
  ggtitle("SenMayo")
ggsave(paste0(paths$results, "/gsea_SenMayo.pdf"), width=0.55*w_in, height=0.35*h_in)

#### SAVE ####

# save.image(paste0(paths$results, "/checkpoint.Rdata"))

# write.csv(resultsDE_ge, paste0(paths$tables, "/resultsDE_ge.csv"))
# write.csv(resultsDE_re, paste0(paths$tables, "/resultsDE_re.csv"))
# 
# for (ont in names(gsea)) {
#   write.csv(gsea[[ont]]$cGASKOvsWT, paste0(paths$tables, "/gsea_", ont, ".csv"))
# }

load(paste0(paths$results, "/checkpoint.Rdata"))

#### PAPER PANELS ####

fformat = ".pdf"

set_theme(theme_classic())
ccodes = c(
  "WT" = "#3A647E",
  "cGAS KO" = "#FE7F2D"
)

geom_bar_auto = function(data, aes_map, ccodes, binwidth_rel = 1/50, dots=T, error=T) {
  yvar = as_label(aes_map$y)
  print(binwidth_rel)
  binwidth = binwidth_rel * max(data[[yvar]], na.rm = TRUE)
  p = ggplot(data, aes_map)+
    geom_bar(stat="summary", color="black")+
    scale_fill_manual(values=ccodes)+
    scale_y_continuous(expand = expansion(mult=c(0, 0.15)))+
    theme(legend.position = "none",
          axis.title.x = element_blank())
  if (dots) p = p + geom_dotplot(binaxis = "y", stackdir = "center", fill='black', binwidth = binwidth)
  if (error) p = p + geom_errorbar(stat = "summary",width = 0.2,linewidth = 0.5)
  return(p)
}

##### Volcano #####

my_volcano(resultsDE_ge, "logFC_cGASKOvsWT", "pval_cGASKOvsWT", "sig_cGASKOvsWT", clip_axes=T, clip_quantile=0.999)+
  guides(color="none")+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank())+
  labs(x=expression(log[2]*"FC"), y=expression(-log[10]*"(pval)"))
ggsave(paste0(paths$paper_figs, "/volcano", fformat), width = w_in*0.45, height=h_in*0.33)

##### GSEA: BP #####

View(gsea$GO_BP$cGASKOvsWT)
View(gsea$GO_CC$cGASKOvsWT)

desc_len = 70
nterms = 8 

# Top 8 by NES
gsea$GO_BP$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  mutate(Sign = ifelse(sign(NES) > 0, "Activated", "Suppressed")) %>%
  group_by(Sign) %>%
  arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  # arrange(p.adjust, -abs(NES)) %>% # First sort by NES, then padj
  slice_head(n=nterms) %>%
  ungroup() %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, sign(NES)*-log(p.adjust) + 0.001 * NES)) %>% # +eps*NES resolves ties
  ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank())
ggsave(paste0(paths$paper_figs, "/gsea_BP_v1", fformat), width = w_in*0.85, height=h_in*0.38)

# Top 8 by significance
gsea$GO_BP$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  mutate(Sign = ifelse(sign(NES) > 0, "Activated", "Suppressed")) %>%
  group_by(Sign) %>%
  # arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  arrange(p.adjust, -abs(NES)) %>% # First sort by NES, then padj
  slice_head(n=nterms) %>%
  ungroup() %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, sign(NES)*-log(p.adjust) + 0.001 * NES)) %>% # +eps*NES resolves ties
  ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank())
ggsave(paste0(paths$paper_figs, "/gsea_BP_v2", fformat), width = w_in*0.85, height=h_in*0.38)

##### GSEA: CC #####

desc_len = 25
nterms = 5 # Top 5 by NES
gsea$GO_CC$cGASKOvsWT %>%
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
  scale_size_continuous(breaks = scales::breaks_pretty(n = 4))+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank())
ggsave(paste0(paths$paper_figs, "/gsea_CC", fformat), width = w_in*0.5, height=h_in*0.33)

##### GSEA: SenMayo #####

gsea$SenMayo$cGASKOvsWT %>%
  mutate(Description = str_remove(Description, "SenMayo_")) %>%
  mutate(Description = str_replace_all(Description, "_", " ")) %>%
  mutate(Description = fct_reorder(Description, NES)) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(pvalue), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank())
ggsave(paste0(paths$paper_figs, "/gsea_SenMayo", fformat), width = w_in*0.55, height=h_in*0.33)

##### L1Md #####

tmp = resultsL1 %>%
  filter(grepl("L1Md", fam)) %>%
  mutate(type = str_extract(fam, "(L1Md[:alnum:]+)_", group=1))%>%
  mutate(subtype = str_extract(fam, "L1Md[:alnum:]+_(.+)", group=1))%>%
  mutate(SigStar = stars.pval(pval_cGASKOvsWT)) %>%
  group_by(type) %>%
  mutate(subtype_rank = (rank(subtype)-1)/7) %>%
  ungroup() %>%
  # mutate(base_col = hue_pal()(length(unique(type)))[match(type, unique(type))]) %>%
  mutate(base_col = brewer_pal(palette = "Dark2")(length(unique(type)))[match(type, unique(type))]) %>%
  mutate(fill_col = darken(base_col, amount = 0.5 * subtype_rank))
ggplot(tmp, aes(subtype, logFC_cGASKOvsWT, label = SigStar)) +
  geom_bar(stat = "identity", aes(fill = fill_col)) +
  geom_text(vjust = 0, size = 6) +
  facet_grid(~type, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_identity() +  # use literal hex colors from fill_col
  labs(y = expression(log[2]*"FC")) +
  guides(fill = "none") +
  theme(
    axis.title.x = element_blank(),
    strip.text = element_text(size = 7),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank())
ggsave(paste0(paths$paper_figs, "/L1Md", fformat), width = w_in*1, height=h_in*0.33)

##### FRIR, FRG #####

all.equal(meta$SampleID, colnames(counts_ge))
meta = meta %>%
  mutate(Genotype = fct_recode(Genotype, "cGAS KO" = "cGASKO")) %>%
  mutate(RiG = colSums(counts_ge))

ggplot(meta, aes(Genotype, RiR/(uniquely_mapped + multi_mapped), fill=Genotype))+
  geom_boxplot()+
  labs(y="Fraction of reads in repeats")+
  stat_compare_means(label.y = 0.11, method="t.test")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
  scale_fill_manual(values=ccodes)+
  guides(fill="none")+
  theme(axis.title.x = element_blank())
ggsave(paste0(paths$paper_figs, "/frir_vs_mapped", fformat), width = w_in*0.3, height=h_in*0.33)
ggplot(meta, aes(Genotype, RiR/annotated, fill=Genotype))+
  geom_boxplot()+
  labs(y="Fraction of reads in repeats")+
  stat_compare_means(label.y = 0.11, method="t.test")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
  scale_fill_manual(values=ccodes)+
  guides(fill="none")+
  theme(axis.title.x = element_blank())
ggsave(paste0(paths$paper_figs, "/frir_vs_annotated", fformat), width = w_in*0.3, height=h_in*0.33)

ggplot(meta, aes(Genotype, RiG/(uniquely_mapped + multi_mapped), fill=Genotype))+
  geom_boxplot()+
  labs(y="Fraction of reads in genes")+
  stat_compare_means(label.y = 0.82, method="t.test")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
  scale_fill_manual(values=ccodes)+
  guides(fill="none")+
  theme(axis.title.x = element_blank())
ggsave(paste0(paths$paper_figs, "/frig_vs_mapped", fformat), width = w_in*0.3, height=h_in*0.33)



#### PAPER ASSEMBLIES ####

fformat = ".pdf"

theme_set(theme_classic())
ccodes = c(
  "WT" = "#3A647E",
  "cGAS KO" = "#FE7F2D"
)

w_in = 7.5 # 8.5 without margin
h_in = 10 # 11 without margin


##### Figure 2 #####

p1 = my_volcano(resultsDE_ge, "logFC_cGASKOvsWT", "pval_cGASKOvsWT", "sig_cGASKOvsWT", clip_axes=T, clip_quantile=0.999)+
  guides(color="none")+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank())+
  labs(x=expression(log[2]*"FC"), y=expression(-log[10]*"(pval)"))

# NES and p.adj range so we can get a common legend for GSEA
desc_len = 30
tmp1 = gsea$GO_CC$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  mutate(Sign = ifelse(sign(NES) > 0, "Activated", "Suppressed")) %>%
  group_by(Sign) %>%
  arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  slice_head(n=5) %>%
  ungroup() %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, sign(NES)*-log(p.adjust))) # +eps*NES resolves ties
desc_len = 35
tmp2 = gsea$GO_BP$cGASKOvsWT %>%
  dplyr::filter(p.adjust < 0.05) %>%
  mutate(Sign = ifelse(sign(NES) > 0, "Activated", "Suppressed")) %>%
  group_by(Sign) %>%
  # arrange(-abs(NES), p.adjust) %>% # First sort by NES, then padj
  arrange(p.adjust, -abs(NES)) %>% # First sort by NES, then padj
  slice_head(n=8) %>%
  ungroup() %>%
  mutate(Description = str_replace(Description, "tumor necrosis factor", "TNF")) %>%
  mutate(Description = if_else(str_length(Description) > desc_len, str_c(str_sub(Description, 1, desc_len-3), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, sign(NES)*-log(p.adjust) + 0.001 * NES))
range_nes = range(c(tmp1$NES, tmp2$NES))
range_padj = range(-log10(c(tmp1$p.adjust, tmp2$p.adjust)))

p2 = ggplot(tmp1, aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red", limits=range_nes)+
  scale_size_continuous(limits=range_padj)+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank())+
  labs(size=expression(-log[10]*("p.adj")))
p3 = ggplot(tmp2, aes(x=1, y=Description, size=-log10(p.adjust), fill=NES))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red", limits=range_nes)+
  scale_size_continuous(limits=range_padj)+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank())+
  labs(size=expression(-log[10]*("p.adj")))+
  guides(fill="none", size="none") # Using the same as p2 so i'll just plot once

(p1+p2+plot_spacer()+free(p3, side="l", type = "space"))+plot_layout(widths=c(8.2,1.8), heights = c(5,8), guides="collect")
ggsave(paste0(paths$paper_figs, "/figure2_assembly", fformat), width = w_in, height=h_in*0.6)

##### Figure S4 #####

p1 = gsea$SenMayo$cGASKOvsWT %>%
  mutate(Description = str_remove(Description, "SenMayo_")) %>%
  mutate(Description = str_replace_all(Description, "_", " ")) %>%
  mutate(Description = fct_reorder(Description, NES)) %>%
  mutate(sig = stars.pval(pvalue)) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(pvalue), fill=NES, label=sig))+
  geom_point(pch = 21)+
  geom_text(color="white")+
  scale_fill_gradient2(low="blue", high="red")+
  labs(size=expression(-log[10]*("p.val")))+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        legend.position = "bottom")
p2 = ggplot(meta, aes(Genotype, RiG/(uniquely_mapped + multi_mapped), fill=Genotype))+
  geom_boxplot()+
  labs(y="Fraction of reads in genes")+
  stat_compare_means(label.y = 0.82, method="t.test")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
  scale_fill_manual(values=ccodes)+
  guides(fill="none")+
  theme(axis.title.x = element_blank())
p3 = ggplot(meta, aes(Genotype, RiR/(uniquely_mapped + multi_mapped), fill=Genotype))+
  geom_boxplot()+
  labs(y="Fraction of reads in repeats")+
  stat_compare_means(label.y = 0.11, method="t.test")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
  scale_fill_manual(values=ccodes)+
  guides(fill="none")+
  theme(axis.title.x = element_blank())
p4 = resultsL1 %>%
  filter(grepl("L1Md", fam)) %>%
  mutate(type = str_extract(fam, "(L1Md[:alnum:]+)_", group=1))%>%
  mutate(subtype = str_extract(fam, "L1Md[:alnum:]+_(.+)", group=1))%>%
  mutate(SigStar = stars.pval(pval_cGASKOvsWT)) %>%
  group_by(type) %>%
  mutate(subtype_rank = (rank(subtype)-1)/7) %>%
  ungroup() %>%
  # mutate(base_col = hue_pal()(length(unique(type)))[match(type, unique(type))]) %>%
  mutate(base_col = brewer_pal(palette = "Dark2")(length(unique(type)))[match(type, unique(type))]) %>%
  mutate(fill_col = darken(base_col, amount = 0.5 * subtype_rank)) %>%
  ggplot(., aes(subtype, logFC_cGASKOvsWT, label = SigStar)) +
  geom_bar(stat = "identity", aes(fill = fill_col)) +
  geom_text(vjust = 0, size = 6) +
  facet_grid(~type, scales = "free_x", space = "free_x") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.2))) +
  scale_fill_identity() +  # use literal hex colors from fill_col
  labs(y = expression(log[2]*"FC")) +
  guides(fill = "none") +
  theme(
    axis.title.x = element_blank(),
    strip.text = element_text(size = 7),
    strip.background = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.line = element_blank())
p=(p1+plot_spacer()+p2+p3+plot_layout(nrow=1, widths=c(1.8,1,3,3)))
ggsave(plot = p, paste0(paths$paper_figs, "/figureS4_assembly1", fformat), width = w_in, height=h_in*0.35)

ggsave(plot = p4, paste0(paths$paper_figs, "/figureS4_assembly2", fformat), width = w_in, height=h_in*0.3)
  