library(tidyverse)
library(data.table)
library(GenomicRanges)
library(GenomicFeatures)
library(patchwork)
library(ggpubr)
library(ChIPseeker)
library(edgeR)
library(AnnotationHub)
library(clusterProfiler)
library(org.Mm.eg.db)
library(BSgenome.Mmusculus.UCSC.mm39)
library(TFBSTools)
library(seqinr)
library(zoo)
library(ggrastr)

setwd("/scratch/fmorandi/internal/John/cGAS_KO/PAPER_CODE/ATAC/")

paths = list()
paths$data = "./pipeline_out"
paths$results = "./results"
paths$gtf = "/scratch/fmorandi/external/references/GRCm39-mm39-Ensembl/Mus_musculus.GRCm39.108.gtf"
paths$tables = paste0(paths$results, "/tables")
paths$seqs = paste0(paths$results, "/seqs")
paths$paper_figs = paste0(paths$results, "/prelim_paper_figs")

dir.create(paths$results, showWarnings = F)
dir.create(paths$tables, showWarnings = F)
dir.create(paths$seqs, showWarnings = F)
dir.create(paths$paper_figs, showWarnings = F)

#### PLOTTING SETTINGS ####

w_in = 7.5 # 8.5 without margin
h_in = 10 # 11 without margin

#### FUNCTIONS ####

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

tpm = function(counts, lengths) {
  tmp = 1e3 * sweep(counts, 2, lengths, "/")
  tmp = 1e6 * sweep(tmp, 1, rowSums(tmp), "/")
  return(tmp)
}

my_volcano = function(table, col, v1, v2, title=NULL, xmax=4, ymax=10, raster=F) {
  tmp = table %>%
    group_by_at(col, .drop = F) %>%
    summarize(n=n()) %>%
    mutate(x = c(-xmax*0.9, 0, xmax*0.9))
  if (raster) {
    p = ggplot(table, aes(x=.data[[v1]], y=-log10(.data[[v2]]), color=.data[[col]]))+
      geom_point_rast(size=0.1)+
      lims(x=c(-xmax, xmax), y=c(0,ymax))+
      ggtitle(title)+
      geom_text(data=tmp, aes(x=x, y=ymax, label=n))
  } else {
    p = ggplot(table, aes(x=.data[[v1]], y=-log10(.data[[v2]]), color=.data[[col]]))+
      geom_point(size=0.1)+
      lims(x=c(-xmax, xmax), y=c(0,ymax))+
      ggtitle(title)+
      geom_text(data=tmp, aes(x=x, y=ymax, label=n))
  }
  if (!is.null(title)) p = p + ggtitle(title)
  return(p)
}

convert_ids_in_string = function(df, id_col, id, symbol) {
  conversions = symbol
  names(conversions) = id
  conversions = conversions[!is.na(names(conversions))]
  for (i in 1:nrow(df)) {
    this_str = df[i, id_col]
    old_ids = unlist(strsplit(this_str, "/"))
    new_ids = conversions[old_ids]
    df[i, id_col] = paste(new_ids, collapse="/")
  }
  return(df)
}

#### PREPROCESS ####

# meta = read.table("meta.txt", header = T)
# 
# data_ocr = fread(paste0(paths$data, "/07_peaksets_and_tables/counts.tsv"),
#                       header=T, data.table=F)
# data_ocr_re = fread(paste0(paths$data, "/07_peaksets_and_tables/counts_ocrs_and_rtes.tsv"),
#                     header=T, data.table=F)
# 
# # Separate feature info
# all.equal(data_ocr$Geneid, data_ocr_re$Geneid[1:nrow(data_ocr)])
# rinfo = data_ocr_re[, 1:6]
# data_ocr = data_ocr[, -c(1:6)]
# data_ocr_re = data_ocr_re[, -c(1:6)]
# 
# # Tidy sample names
# all.equal(colnames(data_ocr), colnames(data_ocr_re))
# colnames(data_ocr) = str_extract(colnames(data_ocr), "/([[:alnum:]]+)_clean.bam", group=1)
# colnames(data_ocr_re) = colnames(data_ocr)
# 
# # Collapse repeat instances
# rinfo$is_rep = grepl("/", rinfo$Geneid)
# data_ocr_re1 = data_ocr_re[!rinfo$is_rep, ]
# data_ocr_re2 = data_ocr_re[rinfo$is_rep, ]
# rinfo1 = rinfo[!rinfo$is_rep, ]
# rinfo2 = rinfo[rinfo$is_rep, ]
# rinfo2[c("superf", "fam", "rep_id")] = split_rep_id(rinfo2$Geneid)
# data_ocr_re2[c("superf", "fam", "Length")] = rinfo2[c("superf", "fam", "Length")]
# data_ocr_re2 = data_ocr_re2 %>%
#   mutate(count = 1) %>%
#   group_by(superf, fam) %>%
#   summarise_all(sum)
# data_ocr_re2 = data_ocr_re2[rowSums(data_ocr_re2[,-c(1,2, 11, 12)]) != 0, ]
# rinfo2 = data_ocr_re2[c("superf", "fam", "Length", "count")] %>%
#   mutate(is_rep = T) %>%
#   mutate(Geneid = paste(superf, fam, sep="/")) %>%
#   mutate(Chr = NA, Start = NA, End = NA, Strand = NA)
# rinfo1 = rinfo1 %>%
#   mutate(superf = NA, fam = NA, count=1)
# rinfo = rbind(
#   rinfo1[c("Geneid", "Chr", "Start", "End", "Strand", "Length", "count", "is_rep", "superf", "fam")],
#   rinfo2[c("Geneid", "Chr", "Start", "End", "Strand", "Length", "count", "is_rep", "superf", "fam")])
# data_ocr_re = rbind(
#   data_ocr_re1,
#   data_ocr_re2[,-c(1,2, 11, 12)])
# rm(rinfo1, rinfo2, data_ocr_re1, data_ocr_re2)
# 
# # Tidy up names
# colnames(rinfo) = c("name", "chr", "start", "end", "strand", "length", "count", "is_rep", "superf", "fam")
# rownames(rinfo) = paste0("r",  1:nrow(rinfo), ifelse(rinfo$is_rep, "_rep", "_ocr"))
# rownames(data_ocr) = rownames(rinfo)[!rinfo$is_rep]
# rownames(data_ocr_re) = rownames(rinfo)
# 
# # Total reads, fraction within repeats
# data_ocr = data.frame(t(data_ocr))
# data_ocr_re = data.frame(t(data_ocr_re))
# all.equal(meta$ID, rownames(data_ocr))
# meta$RIP = rowSums(data_ocr)
# meta$RIPnR = rowSums(data_ocr_re)
# meta$FRIR = rowSums(data_ocr_re[, rinfo$is_rep]) / meta$RIPnR
# 
# # TPM normalization
# hist(unlist(data_ocr), breaks=c(seq(0,1000,10), 1e6), xlim = c(0,1000))
# norm_ocr = tpm(data_ocr, rinfo[colnames(data_ocr), "length"])
# hist(unlist(norm_ocr), breaks=c(seq(0,40,0.5), 1e4), xlim = c(0,40))
# hist(unlist(data_ocr_re), breaks=c(seq(0,1000,10), 1e6), xlim = c(0,1000))
# norm_ocr_re = tpm(data_ocr_re, rinfo$length)
# hist(unlist(norm_ocr_re), breaks=c(seq(0,40,0.5), 1e4), xlim = c(0,40))
# 
# # ChipSeeker
# rinfo$start = as.integer(rinfo$start)
# rinfo$end = as.integer(rinfo$end)
# tmp = subset(rinfo, !is_rep)
# ranges =  GRanges(
#   seqnames = tmp$chr,
#   ranges = IRanges(tmp$start, tmp$end)
# )
# txdb = makeTxDbFromGFF(paths$gtf)
# anno = annotatePeak(ranges, TxDb = txdb, tssRegion = c(-1000, 1000))
# annoDF = as.data.frame(anno) %>%
#   mutate(annotation = str_replace_all(annotation, " \\(.*\\)", "")) %>%
#   mutate(annotation = str_replace_all(annotation, "Distal ", "")) %>%
#   mutate(annotation2 = fct_recode(annotation,
#                                   "Other Genic" = "Downstream",
#                                   "Other Genic" = "3' UTR",
#                                   "Other Genic" = "5' UTR"), .after=annotation) %>%
#   dplyr::select(-(width:strand), -(geneChr:geneStrand), -transcriptId)
# rinfo = rinfo %>%
#   rownames_to_column("id") %>%
#   merge(., annoDF, by.x=c("chr", "start", "end"), by.y=c("seqnames", "start", "end"), all=T) %>%
#   column_to_rownames("id")
# rinfo = rinfo[colnames(data_ocr_re), ]
# 
# save(data_ocr, data_ocr_re, norm_ocr, norm_ocr_re, rinfo, meta, file=paste0(paths$results, "/prepro.Rdata"))

#### QC SUMMARY ####

# fc_report = read.table(paste0(paths$data, "/07_peaksets_and_tables/counts.tsv.summary"),
#                        header=T, row.names=1)
# colnames(fc_report) = str_extract(colnames(fc_report), "([A-Z]+[0-9]+)_clean\\.bam", group=1)
# 
# logs = dir(paste0(paths$data, "/A_mapping_logs"), full.names=T)
# 
# qc = data.frame()
# for (log in logs) {
#   f = read_file(log)
#   sname = str_extract(log, "([A-Z]+[0-9]+)\\.txt", group=1)
#   # TrimGalore! prints the total number of pairs processed after validation
#   #   This is the number of raw pairs
#   #   Sequences that were trimmed too much are removed
#   s = str_extract(f, "Total number of sequences analysed: ([0-9]+)", group=1)
#   qc[sname, "raw"] = as.numeric(s)
#   # Bowtie2 summary includes the number of pairs processed
#   #   This corresponds to the number of pairs after trimming
#   s = str_extract(f, "([0-9]+) reads; of these:", group=1)
#   qc[sname, "trimmed"] = as.numeric(s)
#   # Flagstat prints the number of properly paired reads
#   #   I can confirm that this is the number of reads seen from now on
#   #   By comparing to the remove_duplicates log
#   s = str_extract(f, "([0-9]+) \\+ [0-9]+ properly paired", group=1)
#   qc[sname, "proper_pair"] = as.numeric(s) / 2
#   # I had the pipeline print the number of mito reads from idxstats
#   s = str_extract(f, "Found ([0-9]+) mitochondrial pairs", group=1)
#   qc[sname, "mitochondrial"] = as.numeric(s)
#   # Picard MarkDuplicates prints the number of dupes in the log (divide by 2 to get pairs)
#   s = str_extract(f, "Marking ([0-9]+) records as duplicates", group=1)
#   qc[sname, "duplicates"] = as.numeric(s) / 2
#   # I had the pipeline print the number of clean reads
#   s = str_extract(f, "Clean BAM contains ([0-9]+) pairs", group=1)
#   qc[sname, "clean"] = as.numeric(s)
# }
# 
# # Verify that qc and fc_report match up
# all.equal(colnames(fc_report), rownames(qc))
# all.equal(unname(colSums(fc_report)/2), qc$clean)
# # Get frip
# fc_report = data.frame(t(fc_report))
# qc$cut_sites = fc_report$Assigned + fc_report$Unassigned_NoFeatures
# qc$cut_sites_in_peak = fc_report$Assigned
# qc$frip = qc$cut_sites_in_peak / qc$cut_sites
# 
# write.table(qc, paste0(paths$results, "/qc_summary.tsv"), sep="\t", quote=F)

#### LOAD DATA ####

load(paste0(paths$results, "/prepro.Rdata"))
meta$Genotype = factor(meta$Genotype, levels=c("WT", "KO"))

# Square root transformed tables, for nicer distributions
# summary(unlist(norm_ocr))
# summary(unlist(norm_ocr_re))
# sort(unique(unlist(norm_ocr)))[1:5]
# sort(unique(unlist(norm_ocr_re)))[1:5]
# norm_ocr_sqrt = sqrt(norm_ocr)
# norm_ocr_re_sqrt = sqrt(norm_ocr_re)
# summary(unlist(norm_ocr_sqrt))
# summary(unlist(norm_ocr_re_sqrt))
# save(norm_ocr_sqrt, norm_ocr_re_sqrt, file=paste0(paths$results, "/data_sqrt.Rdata"))
load(paste0(paths$results, "/data_sqrt.Rdata"))

#### EXPLORATION ####

table(meta$Sex, meta$Genotype)

# Total reads and fraction within repeats
meta %>%
  ggplot(., aes(x=Genotype, y=FRIR, fill=Genotype))+
  geom_boxplot()+
  stat_compare_means(label.y=0.37)+
  labs(y="Fraction of reads in repeats")+
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2)))+
  guides(fill="none")
ggsave(paste0(paths$results, "/frir.pdf"), width=0.3*w_in, height=0.25*h_in)
summary(lm(FRIR~Genotype+Sex, data=meta))

# PCA for OCR data
pca = prcomp(norm_ocr_sqrt, center=T, scale.=T)
pca = merge(meta, pca$x, by.x="ID", by.y=0)
ggplot(pca, aes(x=PC1, y=PC2, color=Genotype, shape=Sex))+
  geom_point()
ggsave(paste0(paths$results, "/pca_ocrs.pdf"), width=0.5*w_in, height=0.25*h_in)

# PCA for OCR+REP data
pca = prcomp(norm_ocr_re_sqrt, center=T, scale.=T)
pca = merge(meta, pca$x, by.x="ID", by.y=0)
ggplot(pca, aes(x=PC1, y=PC2, color=Genotype, shape=Sex))+
  geom_point()
ggsave(paste0(paths$results, "/pca_ocrs_and_res.pdf"), width=0.5*w_in, height=0.25*h_in)

#### REMOVE UNANNOTATED OCRS ####

# Not interested in a few unannotated alternative chromosomes
alt_chr = rownames(rinfo[is.na(rinfo$annotation) & !rinfo$is_rep, ])
rinfo = rinfo[!rownames(rinfo) %in% alt_chr, ]
data_ocr = data_ocr[, !colnames(data_ocr) %in% alt_chr]
norm_ocr = norm_ocr[, !colnames(norm_ocr) %in% alt_chr]
norm_ocr_sqrt = norm_ocr_sqrt[, !colnames(norm_ocr_sqrt) %in% alt_chr]
data_ocr_re = data_ocr_re[, !colnames(data_ocr_re) %in% alt_chr]
norm_ocr_re = norm_ocr_re[, !colnames(norm_ocr_re) %in% alt_chr]
norm_ocr_re_sqrt = norm_ocr_re_sqrt[, !colnames(norm_ocr_re_sqrt) %in% alt_chr]

# In the remaining rinfo, if anno == NA then its a repeat
rinfo[is.na(rinfo$annotation), "annotation"] = "Repeat"
rinfo$annotation2 = as.character(rinfo$annotation2)
rinfo[is.na(rinfo$annotation2), "annotation2"] = "Repeat"

#### DIFFERENTIAL ACCESSIBILITY ####

all.equal(rownames(data_ocr_re), meta$ID)
dge = DGEList(t(data_ocr_re), samples = meta)
dge = calcNormFactors(dge)

meta$Genotype = factor(meta$Genotype, levels = c("WT", "KO"))
design = model.matrix(~Genotype+Sex, data = meta)

# # Filter low accessibility
# keep = filterByExpr(dge, design=design)
# sum(keep) # Almost all are kept, so I wont filter

# Fit model
dge = estimateDisp(dge, design)
fit = glmFit(dge, design)
res_da = as.data.frame(glmLRT(fit, coef="GenotypeKO"))

# Save DA decisions
all.equal(rownames(res_da), rownames(rinfo))
rinfo[c("logFC", "pval")] = res_da[c("logFC", "PValue")]
rinfo = rinfo %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  mutate(sig = padj < 0.05) %>%
  mutate(sigs = interaction(sig, sign(logFC))) %>%
  mutate(sigs = fct_recode(sigs, "Sig+" = "TRUE.1", "Sig-" = "TRUE.-1",
                                 "NotSig" = "FALSE.1", "NotSig" = "FALSE.-1")) %>%
  mutate(sigs = factor(sigs, levels=c("Sig-", "NotSig", "Sig+")))

#### VOLCANO + ENRICHMENT ####

p1=rinfo %>%
  dplyr::filter(annotation != "Repeat") %>%
  my_volcano(., "sigs", "logFC", "pval", ymax=7.5)+
  scale_color_manual(values=c("#5555cc", "#999999", "#cc5555")) +
  guides(color="none")
ggsave(paste0(paths$results, "/volcano_ocrs.pdf"), plot = p1, width=0.6*w_in, height=0.4*h_in)

p2=rinfo %>%
  dplyr::filter(annotation == "Repeat") %>%
  dplyr::filter(superf != "Simple_repeat") %>%
  my_volcano(., "sigs", "logFC", "pval", xmax=2.7, ymax=2.1)+
  scale_color_manual(values=c("#999999", "#5555cc","#cc5555")) +
  guides(color="none")
ggsave(paste0(paths$results, "/volcano_res.pdf"), plot = p2, width=0.6*w_in, height=0.4*h_in)

rinfo %>%
  group_by(sigs, annotation2) %>%
  summarize(count = n()) %>%
  ggplot(., aes(x=sigs, y=count, fill=annotation2))+
  geom_bar(position="fill", stat="identity")+
  scale_y_continuous(expand = c(0,0))+
  labs(fill = "Annotation")
ggsave(paste0(paths$results, "/da_annotation.pdf"), width=0.6*w_in, height=0.3*h_in)

#### REPEAT BEHAVIOUR ####

rinfo_l1 = rinfo %>%
  dplyr::filter(is_rep, superf == "LINE/L1") %>%
  mutate(fam = fct_reorder(fam, length/count)) %>%
  mutate(L1Md = grepl("L1Md", fam))

# --- Not filtering by min number of genomic copies ---
p1=ggplot(rinfo_l1, aes(x=fam, y=length/count, group=1))+
  geom_line()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())+
  labs(y="Average length")
p2=rinfo_l1 %>%
  mutate(type = ifelse(L1Md, "L1Md", "Other")) %>%
  ggplot(., aes(x=fam, y=logFC, fill=type))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle=90, hjust=1, size=6))+
  labs(x="LINE1 families sorted by average copy length")
p1/p2 + plot_layout(heights = c(1,2))
ggsave(paste0(paths$results, "/line1s_unfilt.pdf"), width=2*w_in, height=0.5*h_in)

# --- Removing very low copy number line1s ---
min_copies=250
ggplot(rinfo_l1, aes(count))+
  geom_histogram(bins=100)+
  geom_vline(xintercept = min_copies)
table(rinfo_l1$count < min_copies)
rinfo_l1 %>%
  dplyr::filter(count <= min_copies) %>%
  pull(fam) %>%
  sort()

p1 = rinfo_l1 %>%
  dplyr::filter(count > min_copies) %>%
  ggplot(., aes(x=fam, y=length/count, group=1))+
  geom_line()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank())+
  labs(y="Avg. copy\nlength (bp)")
p2 = rinfo_l1 %>%
  dplyr::filter(count > min_copies) %>%
  mutate(type = ifelse(L1Md, "L1Md", "Other")) %>%
  ggplot(., aes(x=fam, y=logFC, fill=type))+
  geom_bar(stat="identity")+
  theme(axis.text.x = element_text(angle=90, size=6, hjust=1))+
  labs(x="LINE1 families sorted by average copy length")
p1/p2 + plot_layout(heights = c(1,2))
ggsave(paste0(paths$results, "/line1s_filt.pdf"), width=2*w_in, height=0.5*h_in)

#### GENE ONTOLOGY ####

rinfo2 = subset(rinfo, !is_rep)

go = list()

# All peaks
universe = unique(rinfo2$geneId)
genes_up = rinfo2 %>%
  dplyr::filter(sigs == "Sig+") %>%
  distinct(geneId, .keep_all = T)
genes_down = rinfo2 %>%
  dplyr::filter(sigs == "Sig-") %>%
  distinct(geneId, .keep_all = T)
intersect(genes_up$geneId, genes_down$geneId) # Negligible
go$up = enrichGO(genes_up$geneId, OrgDb = org.Mm.eg.db, ont = "BP", 
                 universe = universe, keyType = "ENSEMBL", pvalueCutoff = 0.05)
go$down = enrichGO(genes_down$geneId, OrgDb = org.Mm.eg.db, ont = "BP", 
                   universe = universe, keyType = "ENSEMBL", pvalueCutoff = 0.05)
go$up = as.data.frame(go$up)
go$down = as.data.frame(go$down)

# Promoters
universe = unique(rinfo2[rinfo2$annotation == "Promoter", "geneId"])
genes_up = rinfo2 %>%
  dplyr::filter(annotation == "Promoter", sigs == "Sig+") %>%
  distinct(geneId, .keep_all = T)
genes_down = rinfo2 %>%
  dplyr::filter(annotation == "Promoter", sigs == "Sig-") %>%
  distinct(geneId, .keep_all = T)
intersect(genes_up$geneId, genes_down$geneId) # None
go$up_pro = enrichGO(genes_up$geneId, OrgDb = org.Mm.eg.db, ont = "BP", 
                 universe = universe, keyType = "ENSEMBL", pvalueCutoff = 0.05)
go$down_pro = enrichGO(genes_down$geneId, OrgDb = org.Mm.eg.db, ont = "BP", 
                   universe = universe, keyType = "ENSEMBL", pvalueCutoff = 0.05)
go$up_pro = as.data.frame(go$up_pro)
go$down_pro = as.data.frame(go$down_pro)

sapply(go, nrow)

# Convert ensembl ids to symbols and save
ensids = unique(rinfo2$geneId)
id_conv = data.frame(
  ens = ensids,
  symbol = mapIds(org.Mm.eg.db, keys = ensids, keytype = "ENSEMBL", 
                 column = "SYMBOL", multiVals = "first")) %>%
  drop_na()
for (res in names(go)) {
  go[[res]] = convert_ids_in_string(go[[res]], "geneID", id_conv$ens, id_conv$symbol)
  write.csv(go[[res]], file=paste0(paths$tables, "/go_results_", res, ".csv"))
}

rinfo$symbol = id_conv[rinfo$geneId, "symbol"]
write.csv(rinfo, file=paste0(paths$tables, "/resultsDA_all.csv"))
write.csv(subset(rinfo, annotation2 == "Promoter"), file=paste0(paths$tables, "/resultsDA_pro.tsv"))

#### CHECKPOINT ####

# save.image(paste0(paths$results, "/checkpoint.Rdata"))
load(paste0(paths$results, "/checkpoint.Rdata"))

#### QUICK GO PLOT ####

# Only upreg, because not enough promoters are closed for sig enrichment
go$up_pro %>%
  dplyr::filter(p.adjust < 0.05) %>% # Always ensure sig
  mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
  mutate(BgRatio = sapply(BgRatio, function(x) eval(parse(text=x)))) %>%
  mutate(Enrichment = GeneRatio / BgRatio) %>%
  slice_max(Enrichment, n=20) %>% # Most enriched
  mutate(Description = fct_reorder(Description, Enrichment)) %>%
  ggplot(., aes(x=Description, y=Enrichment, fill = -log10(pvalue)))+
  geom_bar(stat="identity")+
  coord_flip()+
  scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
  scale_fill_continuous(low="blue", high="red")+
  ggtitle("Most enriched GO terms associated with promoters that open")
ggsave(paste0(paths$results, "/go_top20_by_enrich.pdf"), width=2*w_in, height=0.6*h_in)

go$up_pro %>%
  dplyr::filter(p.adjust < 0.05) %>% # Always ensure sig
  slice_min(p.adjust, n=20) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
  mutate(BgRatio = sapply(BgRatio, function(x) eval(parse(text=x)))) %>%
  mutate(Enrichment = GeneRatio / BgRatio) %>%
  slice_min(p.adjust, n=20) %>% # Most significant
  mutate(Description = fct_reorder(Description, -pvalue)) %>%
  ggplot(., aes(x=Description, y=Enrichment, fill = -log10(pvalue)))+
  geom_bar(stat="identity")+
  coord_flip()+
  scale_y_continuous(expand = expansion(mult=c(0,0.1)))+
  scale_fill_continuous(low="blue", high="red")+
  ggtitle("Most significant GO terms associated with promoters that open")
ggsave(paste0(paths$results, "/go_top20_by_pval.pdf"), width=2*w_in, height=0.6*h_in)

#### SENMAYO ####

# Nothing sig and in case of down not even enough for the test to run so skip

#### MOTIF ANALYSIS ####

##### Write seqs #####

rinfo3 = rinfo2 %>%
  dplyr::filter(!grepl("\\.", chr)) %>% # Remove weird chromosomes (not many)
  dplyr::filter(length < 2000) # Remove abnormal length (not many)

tmp = list(
  up = subset(rinfo3, sigs == "Sig+"),
  up_pro = subset(rinfo3, sigs == "Sig+" & annotation2 == "Promoter"),
  down = subset(rinfo3, sigs == "Sig-"),
  down_pro = subset(rinfo3, sigs == "Sig-" & annotation2 == "Promoter"),
  bg = subset(rinfo3, sigs == "NotSig"),
  bg_pro = subset(rinfo3, sigs == "NotSig" & annotation2 == "Promoter"))
set.seed(1337)
tmp$bg_subsamp = tmp$bg[sample(nrow(tmp$bg), round(0.1*nrow(tmp$bg))), ]
sapply(tmp, nrow)
peak_seqs = list()
# Sequence in a 500bp window. Works better for enrichment
for (nm in names(tmp)) {
  mid = (tmp[[nm]]$start+tmp[[nm]]$end)/2
  peak_range_fixed = GRanges(
    seqnames = paste0("chr", tmp[[nm]]$chr),
    ranges = IRanges(
      start = mid-250,
      end = mid+250
    )
  )
  peak_seq_fixed = getSeq(BSgenome.Mmusculus.UCSC.mm39, peak_range_fixed)
  peak_seqs[[paste0(nm, "_fixed")]] = peak_seq_fixed
  write.fasta(sequences = as.list(peak_seq_fixed),
              names = rownames(tmp[[nm]]),
              file = paste0(paths$seqs, "/", nm, "_fixed.fa"))
}
# # Exact sequence of the ocr. Doesn't work well with enrichment.
# for (nm in names(tmp)) {
#   peak_range = GRanges(
#     seqnames = paste0("chr", tmp[[nm]]$chr),
#     ranges = IRanges(
#       start = tmp[[nm]]$start,
#       end = tmp[[nm]]$end
#     )
#   )
#   peak_seq = getSeq(BSgenome.Mmusculus.UCSC.mm39, peak_range)
#   peak_seqs[[nm]] = peak_seq
#   write.fasta(sequences = as.list(peak_seq),
#               names = rownames(tmp[[nm]]),
#               file = paste0(paths$seqs, "/", nm, ".fa"))
# }

##### Make figures #####

# Motif pos up
tmp = fread("./motif/site_count_up.txt", data.table = F, header = F)
rownames(tmp) = NULL
cuti = which(grepl("DB", tmp$V1))
motif_pos_up = list()
for (i in 1:(length(cuti)-1)) {
  df = tmp[(cuti[i]+1):(cuti[i+1]-1), ]
  colnames(df) = tmp[cuti[i], ]
  df = data.frame(sapply(df, as.numeric))
  motif_pos_up[[colnames(df)[3]]] = df
}
# Motif pos down
tmp = fread("./motif/site_count_down.txt", data.table = F, header = F)
rownames(tmp) = NULL
cuti = which(grepl("DB", tmp$V1))
motif_pos_down = list()
for (i in 1:(length(cuti)-1)) {
  df = tmp[(cuti[i]+1):(cuti[i+1]-1), ]
  colnames(df) = tmp[cuti[i], ]
  df = data.frame(sapply(df, as.numeric))
  motif_pos_down[[colnames(df)[3]]] = df
}

motif_pos_up$MEME.1 %>%
  mutate(movingMean = rollmean(HACTTCCTCTTTYN, k=25, fill=0) / sum(HACTTCCTCTTTYN)) %>%
  dplyr::filter(abs(DB.0.MOTIF) < 230) %>%
  ggplot(., aes(DB.0.MOTIF, movingMean))+
  geom_line()+
  labs(x="Position", y="Probability")+
  scale_y_continuous(breaks=c(0, 0.002, 0.004, 0.006), minor_breaks = F)+
  theme(axis.title = element_text(size=8))
ggsave(paste0(paths$results, "/motif_pos_spi1.pdf"), width = w_in*0.3, height=h_in*0.1)

motif_pos_down$MEME.1 %>%
  mutate(movingMean = rollmean(TGTTTHYTTTGGC, k=25, fill=0) / sum(TGTTTHYTTTGGC)) %>%
  dplyr::filter(abs(DB.0.MOTIF) < 230) %>%
  ggplot(., aes(DB.0.MOTIF, movingMean))+
  geom_line()+
  labs(x="Position", y="Probability")+
  scale_y_continuous(breaks=c(0, 0.002, 0.004, 0.006), minor_breaks = F)+
  theme(axis.title = element_text(size=8))
ggsave(paste0(paths$results, "/motif_pos_fox.pdf"), width = w_in*0.3, height=h_in*0.1)


# Needs to be space separated
motifs = readJASPARMatrix("./motif/pfms.txt", matrixClass = "PFM")

pdf(paste0(paths$results, "/motif_spi1.pdf"), width = w_in*0.6, height=h_in*0.25)
seqLogo(toICM(motifs[[1]], pseudocounts=0.8, schneider=TRUE), xfontsize = 10, yaxis=F)
dev.off()

pdf(paste0(paths$results, "/motif_fox.pdf"), width = w_in*0.6, height=h_in*0.25)
seqLogo(toICM(motifs[[2]], pseudocounts=0.8, schneider=TRUE), xfontsize = 10, yaxis=F)
dev.off()

#### SEX CHECK ####

all.equal(rownames(norm_ocr), meta$ID)
rinfo[rinfo$symbol %in% c("Xist", "Uty"), ]
meta %>%
  mutate(Xist = norm_ocr[, "r183518_ocr"]) %>%
  ggplot(., aes(Sex, Xist))+
  geom_boxplot()+
  geom_jitter(width=0.2)+
  facet_wrap(~Genotype)
meta %>%
  mutate(Uty = norm_ocr[, "r185191_ocr"]) %>%
  ggplot(., aes(Sex, Uty))+
  geom_boxplot()+
  geom_jitter(width=0.2)+
  facet_wrap(~Genotype)

#### PAPER ASSEMBLIES ####

theme_set(theme_classic())

p1 = rinfo %>%
  dplyr::filter(annotation != "Repeat") %>%
  my_volcano(., "sigs", "logFC", "pval", ymax=7.5, raster=T)+
  scale_color_manual(values=c("#5555cc", "#999999", "#cc5555")) +
  guides(color="none")+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank())+
  labs(x=expression(log[2]*"FC"), y=expression(-log[10]*"(pval)"))
p2 = go$up_pro %>%
  dplyr::filter(p.adjust < 0.05) %>%
  slice_min(p.adjust, n=8) %>%
  mutate(GeneRatio = sapply(GeneRatio, function(x) eval(parse(text=x)))) %>%
  mutate(BgRatio = sapply(BgRatio, function(x) eval(parse(text=x)))) %>%
  mutate(Enrichment = GeneRatio / BgRatio) %>%
  mutate(Description = str_replace(Description, "tumor necrosis factor", "TNF")) %>%
  mutate(Description = if_else(str_length(Description) > 50, str_c(str_sub(Description, 1, 47), "..."), as.character(Description))) %>%
  mutate(Description = fct_reorder(Description, -log10(p.adjust))) %>%
  ggplot(., aes(x=1, y=Description, size=-log10(p.adjust), fill=Enrichment))+
  geom_point(pch = 21)+
  scale_fill_gradient2(low="blue", high="red")+
  scale_size(limits=c(4.5, 6.5), breaks = c(4.5, 5.5, 6.5))+
  theme(axis.title = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank())+
  labs(size=expression(-log[10]*("p.adj")))
p = (p1|p2)+plot_layout(widths=c(8,2))
ggsave(plot = p, paste0(paths$paper_figs, "/top_row.pdf"), width = w_in, height=0.3*h_in)

pdf(paste0(paths$paper_figs, "/motif_spi1.pdf"), width = w_in*0.4, height=h_in*0.2)
seqLogo(toICM(motifs[[1]], pseudocounts=0.8, schneider=TRUE), xfontsize = 10, yaxis=F)
dev.off()
pdf(paste0(paths$paper_figs, "/motif_fox.pdf"), width = w_in*0.4, height=h_in*0.2)
seqLogo(toICM(motifs[[2]], pseudocounts=0.8, schneider=TRUE), xfontsize = 10, yaxis=F)
dev.off()

motif_pos_up$MEME.1 %>%
  mutate(movingMean = rollmean(HACTTCCTCTTTYN, k=25, fill=0) / sum(HACTTCCTCTTTYN)) %>%
  dplyr::filter(abs(DB.0.MOTIF) < 230) %>%
  ggplot(., aes(DB.0.MOTIF, movingMean))+
  geom_line()+
  labs(x="Position", y="Probability")+
  scale_x_continuous(expand = expansion(mult=0))+
  scale_y_continuous(breaks=c(0, 0.002, 0.004, 0.006), minor_breaks = F)+
  theme(axis.title = element_text(size=8))
ggsave(paste0(paths$paper_figs, "/motif_pos_spi1.pdf"), width = 0.32*w_in, height=0.115*h_in)
motif_pos_down$MEME.1 %>%
  mutate(movingMean = rollmean(TGTTTHYTTTGGC, k=25, fill=0) / sum(TGTTTHYTTTGGC)) %>%
  dplyr::filter(abs(DB.0.MOTIF) < 230) %>%
  ggplot(., aes(DB.0.MOTIF, movingMean))+
  geom_line()+
  labs(x="Position", y="Probability")+
  scale_x_continuous(expand = expansion(mult=0))+
  scale_y_continuous(breaks=c(0, 0.002, 0.004, 0.006), minor_breaks = F)+
  theme(axis.title = element_text(size=8))
ggsave(paste0(paths$paper_figs, "/motif_pos_fox.pdf"), width = 0.32*w_in, height=0.115*h_in)


p1 = rinfo_l1 %>%
  dplyr::filter(count > min_copies) %>%
  ggplot(., aes(x=fam, y=length/count, group=1))+
  geom_line()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank())+
  labs(y="Avg. copy\nlength (bp)")
p2 = rinfo_l1 %>%
  dplyr::filter(count > min_copies) %>%
  mutate(type = ifelse(L1Md, "L1Md*", "Other")) %>%
  ggplot(., aes(x=fam, y=logFC, fill=type))+
  geom_col(width=1)+
  labs(x="LINE1 families sorted by average copy length")+
  scale_fill_manual(values=c("#bb33bb", "#999999"))+
  labs(y=expression(log[2]*"FC"), fill="Subfamily")+
  theme(axis.text.x = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
        axis.line = element_blank(),
        legend.position = "bottom")
p1/p2 + plot_layout(heights = c(1,3))
ggsave(paste0(paths$paper_figs, "/line1s_filt.pdf"), width=w_in, height=0.3*h_in)
