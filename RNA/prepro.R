library(tidyverse)
library(data.table)
library(biomaRt)

setwd("/scratch/fmorandi/internal/John/cGAS_KO/PAPER_CODE/RNA")

paths = list()
paths$data = "./pipeline_out"
paths$results = "./results"
dir.create("./results", showWarnings = F)

#### PLOTTING SETTINGS ####

w = 174 # mm
h = 230
w_in = w*0.0393701
h_in = h*0.0393701

#### LOAD DATA ####

# Read counts
files = dir(paste0(paths$data, "/05_counts"), pattern="*cntTable", full.names = T)
counts = list()
for (f in files) {
  sname = str_extract(f, "/([^/]+).cntTable", group=1)
  counts[[sname]] = fread(f, col.names = c("Feature", sname))
}
counts = Reduce(function(x,y) merge(x, y, by="Feature"), counts)
counts = column_to_rownames(counts, "Feature")

# Take out gene info
ginfo = data.frame(
  "Geneid" = rownames(counts), 
  "is_re" = grepl(":", rownames(counts)),
  row.names = rownames(counts))

# Keep genes exressed in at least 3 samples
counts = counts[rowSums(counts > 0) > 3, ]
ginfo = ginfo[ginfo$Geneid %in% rownames(counts), ]

# Read meta
meta = read.table("meta.txt", sep="\t", header=T)
if (!"FileName" %in% colnames(meta)) {
  meta$FileName = meta$SampleID
}

# Read qc info from aligment
qc = list()
files = dir(paste0(paths$data, "/04_mapped"), pattern="*Log.final.out", full.names = T)
for (f in files) {
  sname = str_extract(f, "/([^/]+)_Log.final.out", group=1)
  tmp = read_file(f)
  input_reads = str_extract(tmp, "Number of input reads \\|\\t([[:digit:]]+)", group=1)
  uniquely_mapped = str_extract(tmp, "Uniquely mapped reads number \\|\\t([[:digit:]]+)", group=1)
  multi_mapped = str_extract(tmp, "Number of reads mapped to multiple loci \\|\\t([[:digit:]]+)", group=1)
  unmapped = str_extract_all(tmp, "Number of reads unmapped.*([[:digit:]]+)", simplify = T)
  unmapped = str_extract_all(unmapped, "[[:digit:]]+", simplify = T)
  qc[[sname]] = list(
    input_reads = as.numeric(input_reads),
    uniquely_mapped = as.numeric(uniquely_mapped),
    multi_mapped = as.numeric(multi_mapped),
    unmapped = sum(as.numeric(unmapped))
  )
}

# Read qc info from counting
files = dir(paste0(paths$data, "/B_TEcounts_logs"), pattern="*txt", full.names = T)
for (f in files) {
  sname = str_extract(f, "/\\d+-\\d+-\\d+_([^/]+).txt", group=1)
  tmp = read_file(f)
  annotated = str_extract(tmp, "Total annotated reads = (\\d+) ", group=1)
  unannotated = str_extract(tmp, "Total unannotated reads = (\\d+) ", group=1)
  qc[[sname]][["annotated"]] = as.numeric(annotated)
  qc[[sname]][["unannotated"]] = as.numeric(unannotated)
}

# Merge qc and meta
qc = as.data.frame(do.call(rbind, qc)) %>%
  mutate_all(as.numeric)
meta = merge(meta, qc, by.x="FileName", by.y=0)
counts = counts[, meta$FileName]
write.table(meta, paste0(paths$results, "/qc_summary.tsv"), sep="\t", quote=F, row.names=F)
rm(unmapped, annotated, input_reads, multi_mapped, unannotated, uniquely_mapped, qc)

#### RENAME SAMPLES ####

conv = meta$SampleID
names(conv) = meta$FileName

colnames(counts) = conv[colnames(counts)]

#### SPLIT GINFO AND RINFO ####

counts_re = counts[ginfo$is_re, ]
counts_ge = counts[!ginfo$is_re, ]

rinfo = ginfo[ginfo$is_re, ]
ginfo = ginfo[!ginfo$is_re, ]

rm(counts)

#### CONVERT GENE IDS ####

# Fetch gene symbols and entrez ids
mart = useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
conv = getBM(attributes = c("ensembl_gene_id", "mgi_symbol", "entrezgene_id"),
             filter = "ensembl_gene_id",
             values = ginfo$Geneid, mart = mart)

for (id_type in c("mgi_symbol", "entrezgene_id")) {
  conv2 = conv %>%
    dplyr::filter(!duplicated(ensembl_gene_id)) %>%
    dplyr::filter(!duplicated(id_type))
  rownames(conv2) = conv2$ensembl_gene_id
  ginfo[[id_type]] = conv2[ginfo$Geneid, id_type]
}

# Discard genes with no mgi_symbol
ginfo = drop_na(ginfo)
ginfo = ginfo[!ginfo$mgi_symbol == "", ]
ginfo = ginfo[!duplicated(ginfo$mgi_symbol), ]
counts_ge = counts_ge[ginfo$Geneid, ]
rownames(ginfo) = ginfo$Geneid
rownames(counts_ge) = ginfo[rownames(counts_ge), "mgi_symbol"]
rownames(ginfo) = ginfo$mgi_symbol

#### SPLIT REPEAT NAMES ####

tmp = str_split_fixed(rinfo$Geneid, ":", n = 3)
rinfo$class = tmp[, 3]
rinfo$superf = paste(tmp[, 3], tmp[, 2], sep="/")
rinfo$fam = tmp[, 1]

#### NORMALIZE ####

lib_sizes = colSums(counts_ge)
norm_ge = 1e6 * (counts_ge / lib_sizes)
norm_re = 1e6 * (counts_re / lib_sizes)
ginfo$meanExpr = rowMeans(norm_ge)
rinfo$meanExpr = rowMeans(norm_re)

#### SAVE ####

save(counts_ge, counts_re, norm_ge, norm_re, ginfo, rinfo, meta,
     file=paste0(paths$results, "/prepro.Rdata"))

