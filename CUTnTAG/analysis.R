library(data.table)
library(tidyverse)
library(patchwork)
library(edgeR)
library(GenomicRanges)
library(ggpubr)
library(clusterProfiler)
library(BiocParallel)
library(gtools)

setwd("/scratch/fmorandi/internal/John/cGAS_KO/PAPER_CODE/CUTnTAG")

#### PLOTTING SETTINGS #####

w_in = 7.5 # 8.5 without margin
h_in = 10 # 11 without margin

#### LOAD ####

##### Coverage #####

files = dir("./pipeline_out/08_bin_counts_dups_kept/", pattern = "*.bedgraph", full.names = T)
# files = dir("./pipeline_out/08_bin_counts_dups_rem/", pattern = "*.bedgraph", full.names = T)

# Combine count files
counts = list()
for (f in files) {
  sname = str_remove(f, ".*/")
  sname = str_remove(sname, "_counts_10k.bedgraph")
  counts[[sname]] = fread(f, data.table = F)[, 1:4]
  colnames(counts[[sname]]) = c("chr", "start", "end", sname)
}
counts = Reduce(function(x, y) merge(x, y, by=c("chr", "start", "end")), counts)
counts = counts %>%
  arrange(chr, start, end)

# Remove non canonical chromosomes
counts = counts[!grepl("\\.", counts$chr), ]

# Separate bin info and counts
rownames(counts) = paste0("chr", counts$chr, ":", counts$start, "-", counts$end)
binfo = counts[, 1:3]
counts = counts[, -c(1:3)]

# Remove zero coverage bins
keep = rowSums(counts) > 0
mean(keep)
binfo = binfo[keep, ]
counts = counts[keep, ]
all.equal(rownames(binfo), rownames(counts))

meta = read.table("./meta.tsv", sep="\t", header = T)%>%
  mutate(Genotype = factor(Genotype, levels=c("WT", "cGAS_KO")))
counts = counts[, meta$FileName]

norm = 1e6 * sweep(counts, 2, colSums(counts), FUN="/")
colSums(norm)

metak9 = meta %>%
  dplyr::filter(Target == "H3K9me3")
countsk9 = counts[, metak9$FileName]
normk9 = norm[, metak9$FileName]

metak27 = meta %>%
  dplyr::filter(Target == "H3K27me3")
countsk27 = counts[, metak27$FileName]
normk27 = norm[, metak27$FileName]

##### Annotations #####

cytobands = read.table("./annotation/cytobands2.bed", sep="\t")
colnames(cytobands) = c("chr", "start", "end", "name", "type")
cytobands = cytobands %>%
  dplyr::filter(!grepl("_", chr)) %>%
  mutate(chr = str_replace(chr, "chrM", "MT"))%>%
  mutate(chr = str_remove(chr, "chr")) %>%
  mutate(gvalue = as.numeric(str_extract(type, "\\d+"))) %>%
  mutate(gvalue = replace_na(gvalue, 0))

gaps = read.table("./annotation/gaps.bed", sep="\t")
colnames(gaps) = c("bin", "chr", "start", "end", "ix", "n", "size", "type", "bridge")
gaps = gaps %>%
  dplyr::filter(!grepl("_", chr)) %>%
  mutate(chr = str_replace(chr, "chrM", "MT"))%>%
  mutate(chr = str_remove(chr, "chr"))

blacklist = read.table("./annotation//mm39.excluderanges.bed", sep="\t")
colnames(blacklist) = c("chr", "start", "end", "length", "strand", "type")
blacklist = blacklist %>%
  dplyr::filter(!grepl("_", chr)) %>%
  mutate(chr = str_replace(chr, "chrM", "MT"))%>%
  mutate(chr = str_remove(chr, "chr"))

table(cytobands$chr)
table(gaps$chr)
table(blacklist$chr)
table(binfo$chr)

table(cytobands$type)
table(gaps$type)
table(blacklist$type)

sum(gaps$size) / 3e9 # 2% of genome
sum(blacklist$length) / 3e9 # 5% of genome

ranges_bins = GRanges(
  seqnames = binfo$chr,
  ranges = IRanges(start = binfo$start, end=binfo$end))
ranges_cytobands = GRanges(
  seqnames = cytobands$chr,
  ranges = IRanges(start = cytobands$start, end=cytobands$end))
ranges_gaps = GRanges(
  seqnames = gaps$chr,
  ranges = IRanges(start = gaps$start, end=gaps$end))
ranges_blacklist = GRanges(
  seqnames = blacklist$chr,
  ranges = IRanges(start = blacklist$start, end=blacklist$end))

hits = findOverlaps(ranges_bins, ranges_cytobands)
overlaps = data.frame(
  bin = queryHits(hits),
  cytoband = subjectHits(hits),
  width = width(pintersect(ranges_bins[queryHits(hits)], ranges_cytobands[subjectHits(hits)]))
  ) %>%
  arrange(bin, -width) %>%
  distinct(bin, .keep_all=T)
binfo[overlaps$bin, "cytoband"] = cytobands[overlaps$cytoband, "type"]
binfo[overlaps$bin, "gvalue"] = cytobands[overlaps$cytoband, "gvalue"]
binfo %>%
  dplyr::filter(chr == 10) %>%
  ggplot(., aes(start, y=0, fill=gvalue))+
  geom_tile()+
  scale_fill_gradient(low="white", high="black")

hits = findOverlaps(ranges_bins, ranges_gaps)
overlaps = data.frame(
  bin = queryHits(hits),
  gap = subjectHits(hits),
  width = width(pintersect(ranges_bins[queryHits(hits)], ranges_gaps[subjectHits(hits)]))
) %>%
  arrange(bin, -width) %>%
  distinct(bin, .keep_all=T)
binfo[overlaps$bin, "gap"] = gaps[overlaps$gap, "type"]

hits = findOverlaps(ranges_bins, ranges_blacklist)
overlaps = data.frame(
  bin = queryHits(hits),
  blacklist = subjectHits(hits),
  width = width(pintersect(ranges_bins[queryHits(hits)], ranges_blacklist[subjectHits(hits)]))
) %>%
  arrange(bin, -width) %>%
  distinct(bin, .keep_all=T)
binfo[overlaps$bin, "blacklist"] = blacklist[overlaps$blacklist, "type"]

##### ChromHMM #####

# chromHMM_lung = fread("./annotation/encode3RenChromHmmLungP0_mm39.bed", data.table = F)
# colnames(chromHMM_lung) = c("chr", "start", "end", "state")
# chromHMM_lung$chr = str_remove(chromHMM_lung$chr, "chr")
# unique(chromHMM_lung$chr)
# 
# table(chromHMM_lung$state)
# 
# 
# ranges_bins = GRanges(
#   seqnames = binfo$chr,
#   ranges = IRanges(start = binfo$start, end=binfo$end))
# 
# ranges_hmm = GRanges(
#   seqnames = chromHMM_lung$chr,
#   ranges = IRanges(start = chromHMM_lung$start, end = chromHMM_lung$end),
#   state = chromHMM_lung$state)
# 
# hits = findOverlaps(ranges_bins, ranges_hmm)
# overlap_width = width(pintersect(ranges_bins[queryHits(hits)], ranges_hmm[subjectHits(hits)]))
# overlaps = data.frame(
#   bin_id = queryHits(hits),
#   state = ranges_hmm$state[subjectHits(hits)],
#   overlap_width = overlap_width
# )
# dominant_state = overlaps %>%
#   group_by(bin_id, state) %>%
#   summarise(overlap_width = sum(overlap_width), .groups = "drop") %>%
#   group_by(bin_id) %>%
#   slice_max(overlap_width, n = 1, with_ties = FALSE)
# 
# binfo$chromHMM_state = NA
# binfo$chromHMM_state[dominant_state$bin_id] = dominant_state$state
# binfo$chromHMM_simple = str_extract(binfo$chromHMM_state, "^([^-]+)")
# 
# table(binfo$chromHMM_state)
# table(binfo$chromHMM_simple)

##### ChromHMM 2 #####

chromHMM_universal = fread("./annotation/mm39_chrom_hmm_100_segments.bed", data.table = F)
colnames(chromHMM_universal) = c("chr", "start", "end", "state")

unique(chromHMM_universal$chr)
chromHMM_universal = chromHMM_universal %>%
  mutate(chr = str_remove(chr, "chr")) %>%
  dplyr::filter(!grepl("_", chr)) %>%
  mutate(chr = ifelse(chr == "M", "MT", chr)) %>%
  mutate(simple = str_extract(state, "\\d+_([^\\d+]+)\\d+", group=1))

table(chromHMM_universal$state)
table(chromHMM_universal$simple)
head(chromHMM_universal, n=10)

chromHMM_universal = data.table(chromHMM_universal)
setDT(chromHMM_universal)
setorder(chromHMM_universal, chr, start, end)

chromHMM_universal[, grp := rleid(chr, simple, start == data.table::shift(end) & chr == data.table::shift(chr))]

chromHMM_universal = chromHMM_universal[, .(start = min(start), end = max(end)), by = .(chr, simple, grp)][, grp := NULL]
summary(chromHMM_universal$end - chromHMM_universal$start) # Still small intervals, gotta be efficient in the next steps

ranges_bins = GRanges(
  seqnames = binfo$chr,
  ranges = IRanges(start = binfo$start, end=binfo$end))

ranges_hmm = GRanges(
  seqnames = chromHMM_universal$chr,
  ranges = IRanges(start = chromHMM_universal$start, end = chromHMM_universal$end),
  state = chromHMM_universal$simple)

hits = findOverlaps(ranges_bins, ranges_hmm)
overlap_width = width(pintersect(ranges_bins[queryHits(hits)], ranges_hmm[subjectHits(hits)]))
overlaps = data.table(
  bin_id = queryHits(hits),
  state = ranges_hmm$state[subjectHits(hits)],
  overlap_width = overlap_width
)

# Pick dominant state per bin
# dominant_state = overlaps %>% # Not fast enough
#   group_by(bin_id, state) %>%
#   summarise(overlap_width = sum(overlap_width), .groups = "drop") %>%
#   group_by(bin_id) %>%
#   slice_max(overlap_width, n = 1, with_ties = FALSE)
dominant_state = overlaps[, .(overlap_width = sum(overlap_width)), by = .(bin_id, state)]
dominant_state = dominant_state[dominant_state[, .I[which.max(overlap_width)], by = bin_id]$V1]

binfo$chromHMM_state = NA
binfo$chromHMM_state[dominant_state$bin_id] = dominant_state$state
binfo$chromHMM_state = str_remove(binfo$chromHMM_state, "^m")
binfo$chromHMM_simple = str_extract(binfo$chromHMM_state, "(TSS|ZNF|Prom|Tx|ReprPC|OpenC|Enh|GapArtf|HET|Quies)")

sort(table(binfo$chromHMM_state))
sort(table(binfo$chromHMM_simple))

#### EXPLORE ####

tmp = t(norm[rowSums(norm) > 0, ])
pca = prcomp(tmp, center=T, scale=T)
pca = merge(meta, pca$x[,1:2], by.x="FileName", by.y=0)
p1=ggplot(pca, aes(PC1, PC2, color=Target))+
  geom_point()+
  theme(legend.position = "bottom")

tmp = t(normk9[rowSums(normk9) > 0, ])
pca = prcomp(tmp, center=T, scale=T)
pca = merge(meta, pca$x[,1:2], by.x="FileName", by.y=0)
p2=ggplot(pca, aes(PC1, PC2, color=Genotype))+
  geom_point()+
  ggtitle("H3K9me3")

tmp = t(normk27[rowSums(normk27) > 0, ])
pca = prcomp(tmp, center=T, scale=T)
pca = merge(meta, pca$x[,1:2], by.x="FileName", by.y=0)
p3=ggplot(pca, aes(PC1, PC2, color=Genotype))+
  geom_point()+
  ggtitle("H3K27me3")

p1+(p2+p3+plot_layout(nrow=2, guides="collect")) + plot_layout(widths=c(2,1))
ggsave("./results/pca.pdf", width=1*w_in, height=0.55*h_in)

#### INTERCORS ####

intercors = data.frame(cor(norm)) %>%
  rownames_to_column("S1") %>%
  pivot_longer(cols=-S1, names_to = "S2") %>%
  mutate(Target1 = str_extract(S1, "_(.*)", group=1)) %>%
  mutate(Target2 = str_extract(S2, "_(.*)", group=1)) %>%
  mutate(Target1 = factor(Target1), Target2 = factor(Target2)) %>%
  mutate(S1 = fct_reorder(S1, as.numeric(Target1)))%>%
  mutate(S2 = fct_reorder(S2, as.numeric(Target2)))
ggplot(intercors, aes(S1, S2, fill=value, label=round(value, 2)))+
  geom_tile() +
  geom_text()+
  scale_fill_gradient2(low="blue", high="red")+
  theme(axis.title=element_blank(),
        axis.text.x=element_text(angle=45, hjust=1))
ggsave("./results/coverage_intercors.pdf", width=1*w_in, height=0.6*h_in)

#### NORM TO IGG ####

# DOES NOT WORK! WT IgG too sparse

#### SMOOTHED PLOTS ####

##### Version 1 #####

plot_chromosome = function(binfo, norm, which_chr, win_flank=40) {
  chr_all_smooth = chr_all = cbind(binfo, norm) %>%
    dplyr::filter(chr == which_chr)
  chr_all_smooth[, colnames(norm)] = NA
  for(i in (win_flank+1):(nrow(chr_all_smooth)-win_flank)){
    tmp = chr_all[(i-win_flank):(i+win_flank), colnames(norm)]
    chr_all_smooth[i, colnames(norm)] = colMeans(tmp)
  }
  chr_all_smooth = chr_all_smooth %>%
    dplyr::filter(!is.na(chromHMM_simple))
  
  xrange=c(min(chr_all_smooth$start), max(chr_all_smooth$end))
  p1=chr_all_smooth %>%
    pivot_longer(cols=-c(chr:chromHMM_simple), names_to = "sample") %>%
    dplyr::filter(!grepl("IgG", sample)) %>%
    mutate(Genotype = ifelse(grepl("cGas", sample), "cGAS_KO", "WT")) %>%
    mutate(Target = ifelse(grepl("K9", sample), "H3K9me3", "H3K9me27")) %>%
    group_by(Genotype, Target, chr, start, end) %>%
    summarize(mn = mean(value), sd = sd(value)) %>%
    ggplot(., aes(start, mn, color=Genotype))+
    geom_line()+
    facet_grid(rows=vars(Target), scale="free_y")+
    coord_cartesian(xlim=xrange)+
    theme(axis.title.x = element_blank())
  p2=chr_all_smooth %>%
    dplyr::filter(gap %in% c("centromere", "telomere")) %>%
    ggplot(., aes(xmin=start, ymin=0, xmax=end, ymax = 1, fill=gap))+
    geom_rect()+
    coord_cartesian(xlim=xrange)+
    theme(axis.text.x=element_blank())
  p3 = ggplot(chr_all_smooth, aes(xmin=start, ymin=0, xmax=end, ymax = 1, fill=gvalue))+
    geom_rect()+
    scale_fill_gradient(low="white", high="black")+
    coord_cartesian(xlim=xrange)+
    theme(axis.text.x=element_blank())
  p4 = ggplot(chr_all_smooth, aes(xmin=start, ymin=0, xmax=end, ymax = 1, fill=chromHMM_simple))+
    geom_rect()+
    coord_cartesian(xlim=xrange)
  
  anno_bar = (p2+p3+p4 & 
                theme(axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                      panel.grid = element_blank(),
                      axis.ticks.y = element_blank()))+
    plot_layout(ncol = 1)
  p=p1+anno_bar+plot_layout(heights = c(4,1), guides="collect")
  return(p)
}

p=plot_chromosome(binfo, norm, "19")
ggsave(plot = p, "./results/chr19_coverage.pdf", width = w_in, height=h_in*0.7)

p=plot_chromosome(binfo, norm, "18")
ggsave(plot = p, "./results/chr18_coverage.pdf", width = w_in, height=h_in*0.7)


##### Version 2 #####

plot_chromosome2 = function(binfo, norm, which_chr, win_flank=40*10000) {
  
  chr_all_smooth = chr_all = cbind(binfo, norm) %>%
    dplyr::filter(chr == which_chr) %>%
    dplyr::filter(!is.na(chromHMM_simple))
  chr_all_smooth[, colnames(norm)] = NA
  for(i in 1:nrow(chr_all_smooth)){
    this_pos = mean(chr_all[i, "end"], chr_all[i, "start"])
    tmp = chr_all %>%
      dplyr::filter(start > this_pos - win_flank, end < this_pos + win_flank)
    chr_all_smooth[i, colnames(norm)] = colMeans(tmp[, colnames(norm)])
    majority_state = names(sort(table(tmp[, "chromHMM_simple"]), decreasing = T))[1]
    chr_all_smooth[i, "majority_state"] = majority_state
  }
  
  xrange=c(min(chr_all_smooth$start), max(chr_all_smooth$end))
  p1=chr_all_smooth %>%
    pivot_longer(cols=-c(chr:chromHMM_simple, majority_state), names_to = "sample") %>%
    dplyr::filter(!grepl("IgG", sample)) %>%
    mutate(Genotype = ifelse(grepl("cGas", sample), "cGAS_KO", "WT")) %>%
    mutate(Target = ifelse(grepl("K9", sample), "H3K9me3", "H3K9me27")) %>%
    group_by(Genotype, Target, chr, start, end) %>%
    summarize(mn = mean(value), sd = sd(value)) %>%
    ggplot(., aes(start, mn, color=Genotype))+
    geom_line()+
    facet_grid(rows=vars(Target), scale="free_y")+
    coord_cartesian(xlim=xrange)+
    theme(axis.title.x = element_blank())
  p2 = ggplot(chr_all_smooth, aes(xmin=start, ymin=0, xmax=end, ymax = 1, fill=gvalue))+
    geom_rect()+
    scale_y_continuous(expand = expansion(mult=0))+
    scale_fill_gradient(low="white", high="black")+
    coord_cartesian(xlim=xrange)+
    theme(axis.text.x=element_blank(),
          panel.border = element_rect())
  p3 = ggplot(chr_all_smooth, aes(xmin=start, ymin=0, xmax=end, ymax = 1, fill=majority_state))+
    geom_rect()+
    scale_y_continuous(expand = expansion(mult=0))+
    coord_cartesian(xlim=xrange)+
    theme(panel.border = element_rect())
  
  anno_bar = (p2+p3 & 
                theme(axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                      panel.grid = element_blank(),
                      axis.ticks.y = element_blank()))+
    plot_layout(ncol = 1)
  p=p1+anno_bar+plot_layout(heights = c(4,1), guides="collect")&scale_x_continuous(expand = expansion(mult=0))
  return(p)
}

p=plot_chromosome2(binfo, norm, "19")
ggsave(plot = p, "./results/chr19_coverage2.pdf", width = w_in, height=h_in*0.7)

##### For paper #####

which_chr = "19"
win_flank=40*10000

chr_all_smooth = chr_all = cbind(binfo, norm) %>%
  dplyr::filter(chr == which_chr) %>%
  dplyr::filter(!is.na(chromHMM_simple))
chr_all_smooth[, colnames(norm)] = NA
for(i in 1:nrow(chr_all_smooth)){
  this_pos = mean(chr_all[i, "end"], chr_all[i, "start"])
  tmp = chr_all %>%
    dplyr::filter(start > this_pos - win_flank, end < this_pos + win_flank)
  chr_all_smooth[i, colnames(norm)] = colMeans(tmp[, colnames(norm)])
  majority_state = names(sort(table(tmp[, "chromHMM_simple"]), decreasing = T))[1]
  chr_all_smooth[i, "majority_state"] = majority_state
}

xrange=c(min(chr_all_smooth$start), max(chr_all_smooth$end))
tmp = chr_all_smooth %>%
  pivot_longer(cols=-c(chr:chromHMM_simple, majority_state), names_to = "sample") %>%
  dplyr::filter(!grepl("IgG", sample)) %>%
  mutate(Genotype = ifelse(grepl("cGas", sample), "cGAS_KO", "WT")) %>%
  mutate(Genotype = factor(Genotype, levels=c("WT", "cGAS_KO"))) %>%
  mutate(Target = ifelse(grepl("K9", sample), "H3K9me3", "H3K27me3")) %>%
  group_by(Genotype, Target, chr, start, end) %>%
  summarize(mn = mean(value), sd = sd(value))

p1= tmp %>%
  dplyr::filter(Target == "H3K9me3") %>%
  ggplot(., aes(start, mn, color=Genotype))+
  geom_line()+
  coord_cartesian(xlim=xrange)+
  scale_color_manual(values=c("#3A647E", "#FE7F2D"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(),
        axis.line = element_blank())+
  labs(y = "Avg. H3K9me3 CPM")
p2 = ggplot(chr_all_smooth, aes(xmin=start, ymin=0, xmax=end, ymax = 1, fill=majority_state))+
  geom_rect()+
  scale_y_continuous(expand = expansion(mult=0))+
  # scale_fill_brewer(palette="Set2")+
  scale_fill_manual(values = scales::col_darker(
    scales::brewer_pal(type="qual", "Set2")(7), amount = 8
  ))+
  coord_cartesian(xlim=xrange)+
  theme(panel.border = element_rect(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        panel.grid = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line = element_blank())+
  labs(fill="Chromatin\nState", x="Chr 19 coordinate")
p=p1+p2+plot_layout(heights = c(6,1), guides="collect") &
  scale_x_continuous(expand = expansion(mult=0))
ggsave(plot = p, paste0("./results/prelim_paper_figs/chr19_h3k9me3.pdf"), width = w_in, height=0.25*h_in)

p1= tmp %>%
  dplyr::filter(Target == "H3K27me3") %>%
  ggplot(., aes(start, mn, color=Genotype))+
  geom_line()+
  coord_cartesian(xlim=xrange)+
  scale_color_manual(values=c("#3A647E", "#FE7F2D"))+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(),
        axis.line = element_blank())+
  labs(y = "Avg. H3K27me3 CPM")
p=p1+p2+plot_layout(heights = c(6,1), guides="collect") &
  scale_x_continuous(expand = expansion(mult=0))
ggsave(plot = p, paste0("./results/prelim_paper_figs/chr19_h3k27me3.pdf"), width = w_in, height=0.25*h_in)

#### DIFFERENTIAL TESTING ####

all.equal(colnames(countsk9), metak9$FileName)
all.equal(colnames(countsk27), metak27$FileName)

##### H3K27me3 #####

dge = DGEList(countsk27, samples = metak27)
design = model.matrix(~Genotype, data=dge$samples)

dge = dge[rowSums(dge$counts) > 0, ]

dge = calcNormFactors(dge)
dge$samples$norm.factors
dge = estimateDisp(dge, design)

fit = glmFit(dge, design)
dres_k27 = data.frame(glmLRT(fit, coef="GenotypecGAS_KO")) %>%
  dplyr::select(logFC, PValue) %>%
  dplyr::rename("pval" = "PValue") %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  merge(binfo, ., by=0) %>%
  arrange(chr, start, end)

dres_k27$sig = "Not Sig."
dres_k27[dres_k27$padj < 0.05 & dres_k27$logFC > 0, "sig"] = "Sig. Up"
dres_k27[dres_k27$padj < 0.05 & dres_k27$logFC < 0, "sig"] = "Sig. Down"
dres_k27$sig = factor(dres_k27$sig, levels = c("Sig. Down", "Not Sig.", "Sig. Up"))
table(dres_k27$sig)

##### H3K9me3 #####

dge = DGEList(countsk9, samples = metak9)
design = model.matrix(~Genotype, data=dge$samples)

dge = dge[rowSums(dge$counts) > 0, ]

dge = calcNormFactors(dge)
dge$samples$norm.factors
dge = estimateDisp(dge, design)

fit = glmFit(dge, design)
dres_k9 = data.frame(glmLRT(fit, coef="GenotypecGAS_KO")) %>%
  dplyr::select(logFC, PValue) %>%
  dplyr::rename("pval" = "PValue") %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  merge(binfo, ., by=0) %>%
  arrange(chr, start, end)

dres_k9$sig = "Not Sig."
dres_k9[dres_k9$padj < 0.05 & dres_k9$logFC > 0, "sig"] = "Sig. Up"
dres_k9[dres_k9$padj < 0.05 & dres_k9$logFC < 0, "sig"] = "Sig. Down"
dres_k9$sig = factor(dres_k9$sig, levels = c("Sig. Down", "Not Sig.", "Sig. Up"))
table(dres_k9$sig)

##### H3K27me3 by chromatin state #####

all.equal(rownames(binfo), rownames(countsk27))
counts_k27_by_state = countsk27 %>%
  mutate(state = binfo$chromHMM_state) %>%
  group_by(state) %>%
  summarize_all(sum) %>%
  drop_na() %>%
  column_to_rownames("state")

dge = DGEList(counts_k27_by_state, samples = metak27)
design = model.matrix(~Genotype, data=dge$samples)

dge = calcNormFactors(dge)
dge$samples$norm.factors
dge = estimateDisp(dge, design)

fit = glmFit(dge, design)
dres_k27_by_state = data.frame(glmLRT(fit, coef="GenotypecGAS_KO")) %>%
  dplyr::select(logFC, PValue) %>%
  dplyr::rename("pval" = "PValue") %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  rownames_to_column("state")

dres_k27_by_state$sig = "Not Sig."
dres_k27_by_state[dres_k27_by_state$padj < 0.05 & dres_k27_by_state$logFC > 0, "sig"] = "Sig. Up"
dres_k27_by_state[dres_k27_by_state$padj < 0.05 & dres_k27_by_state$logFC < 0, "sig"] = "Sig. Down"
dres_k27_by_state$sig = factor(dres_k27_by_state$sig, levels = c("Sig. Down", "Not Sig.", "Sig. Up"))
table(dres_k27_by_state$sig)

##### H3K9me3 by chromatin state #####

all.equal(rownames(binfo), rownames(countsk9))
counts_k9_by_state = countsk9 %>%
  mutate(state = binfo$chromHMM_state) %>%
  group_by(state) %>%
  summarize_all(sum) %>%
  drop_na() %>%
  column_to_rownames("state")

dge = DGEList(counts_k9_by_state, samples = metak9)
design = model.matrix(~Genotype, data=dge$samples)

dge = calcNormFactors(dge)
dge$samples$norm.factors
dge = estimateDisp(dge, design)

fit = glmFit(dge, design)
dres_k9_by_state = data.frame(glmLRT(fit, coef="GenotypecGAS_KO")) %>%
  dplyr::select(logFC, PValue) %>%
  dplyr::rename("pval" = "PValue") %>%
  mutate(padj = p.adjust(pval, method="BH")) %>%
  rownames_to_column("state")

dres_k9_by_state$sig = "Not Sig."
dres_k9_by_state[dres_k9_by_state$padj < 0.05 & dres_k9_by_state$logFC > 0, "sig"] = "Sig. Up"
dres_k9_by_state[dres_k9_by_state$padj < 0.05 & dres_k9_by_state$logFC < 0, "sig"] = "Sig. Down"
dres_k9_by_state$sig = factor(dres_k9_by_state$sig, levels = c("Sig. Down", "Not Sig.", "Sig. Up"))

#### CHECKPOINT ####

# save.image("./results/checkpoint.Rdata")
load("./results/checkpoint.Rdata")
table(dres_k9_by_state$sig)

#### LOCAL ENRICHMENT ####

dres_OoE = function(dres, win_flank) {
  win_size = win_flank * 2 + 1
  print(paste("Win size:", round(win_size * 1e4 / 1e6, 2), "Mb"))
  sig_up = as.integer(dres$sig == "Sig. Up")
  sig_down = as.integer(dres$sig == "Sig. Down")
  chrs = dres$chr
  # cumulative sums per chromosome
  split_idx = split(seq_along(chrs), chrs)
  for (chr in names(split_idx)) {
    idx = split_idx[[chr]]
    if (length(idx) < win_size) next   # skip short chromosomes
    
    su = sig_up[idx]
    sd = sig_down[idx]
    
    # cumulative sums (pad with 0 at start)
    csu = c(0, cumsum(su))
    csd = c(0, cumsum(sd))
    
    # window sums: difference of cumsums
    nup = csu[(win_size+1):length(csu)] - csu[1:(length(csu)-win_size)]
    ndown = csd[(win_size+1):length(csd)] - csd[1:(length(csd)-win_size)]
    
    # insert back into result
    valid_idx = idx[(win_flank+1):(length(idx)-win_flank)]
    dres[valid_idx, "nup_win"] = nup
    dres[valid_idx, "ndown_win"] = ndown
  }
  p_up = mean(dres$sig == "Sig. Up")
  p_down = mean(dres$sig == "Sig. Down")
  dres$z_up = (dres$nup_win - win_size*p_up) / sqrt(win_size*p_up*(1-p_up))
  dres$z_down = (dres$ndown_win - win_size*p_down) / sqrt(win_size*p_down*(1-p_down))
  dres$p_up= pbinom(dres$nup_win - 1, size = win_size, prob = p_up, lower.tail = FALSE)
  dres$p_down = pbinom(dres$ndown_win - 1, size = win_size, prob = p_down, lower.tail = FALSE)
  return(dres)
}

dres_k9 = dres_OoE(dres_k9, 100)
dres_k27 = dres_OoE(dres_k27, 100)

get_banks = function(dres, indicator) {
  out = list(
    "chr" = c(),
    "start" = c(),
    "end" = c()
  )
  in_region = F
  for (i in 1:length(indicator)) {
    if (indicator[i] & !in_region) {
      in_region = T
      out$chr = c(out$chr, dres$chr[i])
      out$start = c(out$start, dres$start[i])
    } else if(!indicator[i] & in_region) {
      in_region = F
      out$end = c(out$end, dres$end[i])
    }
  }
  if (in_region) {
    out$end = c(out$end, dres$end[i])
  }
  return(data.frame(out))
}
banks_k27_down = get_banks(dres_k27, !is.na(dres_k27$p_down) & dres_k27$p_down < 0.05)
banks_k27_up = get_banks(dres_k27, !is.na(dres_k27$p_up) & dres_k27$p_up < 0.05)
banks_k9_down = get_banks(dres_k9, !is.na(dres_k9$p_down) & dres_k9$p_down < 0.05)
banks_k9_up = get_banks(dres_k9, !is.na(dres_k9$p_up) & dres_k9$p_up < 0.05)
write.table(banks_k27_down, "./results/banks/banks_k27_down.bed", col.names=F, sep="\t", quote=F, row.names = F)
write.table(banks_k27_up, "./results/banks/banks_k27_up.bed", col.names=F, sep="\t", quote=F, row.names = F)
write.table(banks_k9_down, "./results/banks/banks_k9_down.bed", col.names=F, sep="\t", quote=F, row.names = F)
write.table(banks_k9_up, "./results/banks/banks_k9_up.bed", col.names=F, sep="\t", quote=F, row.names = F)

p=dres_k9 %>%
  mutate(chr = factor(chr, levels=c(as.character(1:19), c("X", "Y")))) %>%
  ggplot(.)+
  geom_line(aes(start, -log10(p_up)), color="red")+
  geom_line(aes(start, -log10(p_down)), color="blue")+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  facet_wrap(~chr, scales = "free", ncol=3)
ggsave(plot=p, "./results/local_enrichment_H3K9me3.pdf", width=1*w_in, height=1*h_in)

p=dres_k27 %>%
  mutate(chr = factor(chr, levels=c(as.character(1:19), c("X", "Y")))) %>%
  ggplot(.)+
  geom_line(aes(start, -log10(p_up)), color="red")+
  geom_line(aes(start, -log10(p_down)), color="blue")+
  geom_hline(yintercept = -log10(0.05), linetype="dashed")+
  facet_wrap(~chr, scales = "free", ncol=3)
ggsave(plot=p, "./results/local_enrichment_H3K27me3.pdf", width=1*w_in, height=1*h_in)

#### ENTROPY ####

all.equal(colnames(norm), meta$FileName)

p=meta %>%
  mutate(entropy = apply(norm/1e6, 2, TFBSTools::shannon.entropy)) %>%
  dplyr::filter(Target != "IgG") %>%
  ggplot(., aes(Genotype, entropy))+
  geom_boxplot()+
  facet_wrap(~Target)+
  stat_compare_means(method="t.test")
p+scale_y_continuous(expand = expansion(mult=c(0.1, 0.2)))
ggsave("./results/entropy.pdf", width=0.5*w_in, height=0.3*h_in)

#### CHROM STATE ENRICHMENT ####

# GSEA TOO SLOW
# my_gseGOcustom = function(table, logFC_column, state_column) {
#   options(warn = 1)
#   # Make gmt
#   gmt = data.frame(
#     "term" = table[[state_column]],
#     "gene" = table$Row.names
#   )
#   # Make gene list
#   gene_list = table[[logFC_column]]
#   names(gene_list) = table$Row.names
#   gene_list = sort(gene_list, decreasing = TRUE)
#   # Run GSEA
#   set.seed(1337)
#   res = GSEA(
#     geneList=gene_list,
#     TERM2GENE = gmt,
#     verbose = TRUE,
#     pvalueCutoff = 1.1,
#     maxGSSize = 200000,
#     BPPARAM = SerialParam())
#   return(res)
# }
# 
# gsea_state = my_gseGOcustom(dres_k9, "logFC", "chromHMM_state")

p1=dres_k27_by_state %>%
  mutate(simple = str_extract(state, "(TSS|ZNF|Prom|Tx|ReprPC|OpenC|Enh|GapArtf|HET|Quies)")) %>%
  mutate(star = stars.pval(padj)) %>%
  ggplot(., aes(state, logFC, label=star, fill=simple))+
  geom_col()+
  geom_text(aes(y=1.1*logFC), size=5)+
  facet_wrap(~simple, scales = "free_x", space="free_x")+
  ggtitle("H3K27me3")
p2=dres_k9_by_state %>%
  mutate(simple = str_extract(state, "(TSS|ZNF|Prom|Tx|ReprPC|OpenC|Enh|GapArtf|HET|Quies)")) %>%
  mutate(star = stars.pval(padj)) %>%
  ggplot(., aes(state, logFC, label=star, fill=simple))+
  geom_col()+
  geom_text(aes(y=1.1*logFC), size=5)+
  facet_wrap(~simple, scales = "free_x", space="free_x")+
  ggtitle("H3K9me3")
p1/p2&theme(axis.title=element_blank())
ggsave("./results/chr_states_logFC.pdf", width=1*w_in, height=0.5*h_in)


dres_k9_by_state %>%
  mutate(simple = str_extract(state, "(TSS|ZNF|Prom|Tx|ReprPC|OpenC|Enh|GapArtf|HET|Quies)")) %>%
  mutate(star = stars.pval(padj)) %>%
  ggplot(., aes(state, logFC, label=star, fill=simple))+
  geom_col()+
  geom_text(aes(y=1.1*logFC), size=5)+
  facet_grid(cols=vars(simple), scales = "free_x", space="free_x")+
  guides(fill="none")+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=1, angle=90),
        strip.background = element_blank(),
        panel.border = element_rect(),
        axis.line = element_blank())
ggsave("./results/prelim_paper_figs/chr_states_h3k9me3.pdf", width=1*w_in, height=0.25*h_in)

dres_k27_by_state %>%
  mutate(simple = str_extract(state, "(TSS|ZNF|Prom|Tx|ReprPC|OpenC|Enh|GapArtf|HET|Quies)")) %>%
  mutate(star = stars.pval(padj)) %>%
  ggplot(., aes(state, logFC, label=star, fill=simple))+
  geom_col()+
  geom_text(aes(y=1.1*logFC), size=5)+
  facet_grid(cols=vars(simple), scales = "free_x", space="free_x")+
  guides(fill="none")+
  theme_classic()+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(hjust=1, angle=90),
        strip.background = element_blank(),
        panel.border = element_rect(),
        axis.line = element_blank())
ggsave("./results/prelim_paper_figs/chr_states_h3k27me3.pdf", width=1*w_in, height=0.25*h_in)


all.equal(dres_k27_by_state$state, dres_k9_by_state$state)
data.frame(
  "state" = dres_k27_by_state$state,
  "logFC_k27" = dres_k27_by_state$logFC,
  "logFC_k9" = dres_k9_by_state$logFC) %>%
  ggplot(., aes(logFC_k27, logFC_k9, label=state))+
  geom_text()+
  stat_cor(label.y = -0.5)+
  geom_smooth(method="lm")
ggsave("./results/chr_states_logFC_cor.pdf", width=0.7*w_in, height=0.5*h_in)

#### LOGFC BY QUARTILE ####

avg_k9_ctrl = rowMeans(normk9[dres_k9$Row.names, metak9[metak9$Genotype == "WT", "FileName"]])
dres_k92 = dres_k9 %>%
  mutate(cltr_quartile = quantcut(avg_k9_ctrl, q = 4))%>%
  mutate(cltr_quartile_id = as.factor(as.numeric(cltr_quartile)))

avg_k27_ctrl = rowMeans(normk27[dres_k27$Row.names, metak27[metak27$Genotype == "WT", "FileName"]])
dres_k272 = dres_k27 %>%
  mutate(cltr_quartile = quantcut(avg_k27_ctrl, q = 4))%>%
  mutate(cltr_quartile_id = as.factor(as.numeric(cltr_quartile)))

p1 = ggplot(dres_k92, aes(cltr_quartile_id, logFC, fill=cltr_quartile_id))+
  geom_boxplot(outliers = F)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ggtitle("H3K27me3")
p2 = ggplot(dres_k272, aes(cltr_quartile_id, logFC, fill=cltr_quartile_id))+
  geom_boxplot(outliers = F)+
  geom_hline(yintercept = 0, linetype="dashed")+
  ggtitle("H3K9me3")
p1+p2+plot_layout(guides="collect")
ggsave("./results/logFC_by_quartile.pdf", width=1*w_in, height=0.4*h_in)

#### CHECK PROPER NORM ####

# Do the norm factors from EdgeR corrctly fix the skw due to some higly covered bins at telos etc...?
chr19_all = list()
win_flank = 40
for(samp in c("WT1", "WT2", "WT3", "cGAS1", "cGAS2", "cGAS3")) {
  tmp = norm %>%
    dplyr::filter(chr == 19) %>%
    pivot_longer(cols=-c("chr", "start", "end"), names_to = "sample") %>%
    dplyr::filter(sample == samp) %>%
    arrange(chr, start, end) %>%
    mutate(value = value)
  # mutate(value = log1p(value))
  for(i in (win_flank+1):(nrow(tmp)-win_flank)){
    tmp[i, "mean"] = mean(tmp$value[(i-win_flank):(i+win_flank)])
    tmp[i, "sd"] = sd(tmp$value[(i-win_flank):(i+win_flank)])
  }
  chr19_all[[samp]] = tmp
}
chr19_all = do.call(rbind, chr19_all)
chr19_all[, "norm_f"] = dge$samples[chr19_all$sample, "norm.factors"]

p_before=chr19_all %>%
  mutate(Genotype = ifelse(grepl("cGAS", sample), "cGAS", "WT")) %>%
  group_by(Genotype, chr, start, end) %>%
  summarize(mn = mean(mean), sd = sd(mean)) %>%
  ggplot(., aes(start, mn, color=Genotype))+
  geom_line()
p_after=chr19_all %>%
  mutate(Genotype = ifelse(grepl("cGAS", sample), "cGAS", "WT")) %>%
  mutate(mean = mean / norm_f) %>%
  group_by(Genotype, chr, start, end) %>%
  summarize(mn = mean(mean), sd = sd(mean)) %>%
  ggplot(., aes(start, mn, color=Genotype))+
  geom_line()
p_before / p_after
# Yes looks good

#### SAVE FILES ####

write_csv(dres_k9_by_state, "./results/tables/dres_h3k9me3.csv")
write_csv(dres_k27_by_state, "./results/tables/dres_h3k27me3.csv")
