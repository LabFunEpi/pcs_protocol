library(BPCells)
library(readr)
library(BSgenome) 
library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(TFBSTools)
library(BiocParallel)
library(Seurat)
library(Signac)
library(chromVAR)
library(motifmatchr)
library(data.table)
library(dplyr)
library(tidyr)
library(Matrix)
library(magrittr)
library(stringr)
library(tidyverse)
register(SerialParam())

# part A: read data and prep Seurat object 
mat = open_matrix_10x_hdf5("counts_10x.h5", feature_type = "Peaks")
md = read_csv("metadata.csv")
chrom_assay = CreateChromatinAssay(counts = mat, sep = c("-", "-"))
obj = CreateSeuratObject(counts = chrom_assay, assay = "ATAC", genome = "hg38")
obj$celltype = md$celltype
obj$grouping1 = md$grouping1
obj$grouping2 = md$grouping2
obj = obj %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 20) %>% RunSVD()
saveRDS(obj, "seurat_obj.rds")

# part B: Prepare ReMap data
# download remap BED file from https://remap.univ-amu.fr/download_page and extract it (outside R)
# gunzip remap2022_nr_macs2_hg38_v1_0.bed.gz

# additions to remap BED file (peaks in addons.bed)
# cut -d":" -f1 remap2022_nr_macs2_hg38_v1_0.bed > remap2022_nr_macs2_hg38_v1_0_MOD.bed
# grep -E '^chr([0-9]+|X|Y|M)[[:space:]]' remap2022_nr_macs2_hg38_v1_0_MOD.bed > remap2022_nr_macs2_hg38_v1_0_MOD_clean.bed
# cat remap2022_nr_macs2_hg38_v1_0_MOD_clean.bed addons.bed | sort -k1,1 -k2,2n > myremap.bed

# part C: create remap assay (OPTIONAL - if not doing this, continue with filtered <- readRDS("seurat_obj.rds") in part D)
# first set up genome stuff
genome <- BSgenome.Hsapiens.UCSC.hg38
keepBSgenomeSequences <- function(genome, seqnames) {
  stopifnot(all(seqnames %in% seqnames(seqinfo(genome))))
  genome@user_seqnames <- setNames(seqnames, seqnames)
  genome@seqinfo <- genome@seqinfo[seqnames]
  genome
}                     
sequences_to_keep <- paste0("chr", c(1:22, "X", "Y"))
genome <- keepBSgenomeSequences(genome, sequences_to_keep)

# load seurat object
filtered <- readRDS("seurat_obj.rds")

# export peaks BED from atac assay
data.frame(Peak = rownames(filtered)) %>%
  separate(col = "Peak", sep = "-", into = c("chr", "start", "end")) %>%
  write.table("multiome_peaks.bed", sep="\t", quote=FALSE, 
              row.names=FALSE, col.names=FALSE)

# run the following command outside R
# bedtools intersect -a multiome_peaks.bed -b myremap.bed -wa -wb > peaks_x_remap.txt

hits <- fread("peaks_x_remap.txt", col.names = c("peak", "TF"))
hits <- hits[grepl("^chr([0-9]+|X|Y|M)$", sub("-.*", "", hits$peak))]

# build the named list and drop any TFs with no overlapping peaks
TF_features <- split(hits$peak, paste0("remap_", hits$TF))
TF_features <- TF_features[lengths(TF_features) > 0]

# add TF_features as an assay in the seurat object
filtered <- AddChromatinModule(filtered, features = TF_features, 
                               genome = genome, assay = "ATAC")
remapscores <- filtered@meta.data %>% select(starts_with("remap-")) %>% t()
rownames(remapscores) <- str_sub(rownames(remapscores), 7L, -1L)
filtered[["remap"]] <- CreateAssayObject(data = remapscores)
filtered@meta.data <- filtered@meta.data %>% select(!starts_with("remap-"))
saveRDS(filtered, "seurat_obj_remap.rds")


# part D: defining the PCS from the atac assay
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(tidyverse)
# workaround to avoid seqinfo inherited method error - ignore if you don't have the error!
#setMethod("seqinfo", signature(x = "ChromatinAssay"), function(x) {
#  GenomeInfoDb::seqinfo(x@ranges)
#})

# load seurat object
filtered <- readRDS("seurat_obj_remap.rds")
DefaultAssay(filtered) <- "ATAC"
# subset to celltype of interest
Tumor_Epi <- subset(filtered, subset = celltype == "Epithelial_cells")

# get differentially accessible peaks for both comparisons
daps_NaiveNACT <- FindMarkers(Tumor_Epi, ident.1 = "NACT", ident.2 = "Naive",min.pct = 0 , group.by = "grouping1") %>% rownames_to_column("peak")
write.table(daps_NaiveNACT, "daps_NACT_vs_Naive.tsv",  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

daps_SensRes <- FindMarkers(Tumor_Epi, ident.1 = "Resistant", ident.2 = "Sensitive",min.pct = 0 , group.by = "grouping2") %>% rownames_to_column("peak")
write.table(daps_SensRes, "daps_Resistant_vs_Sensitive.tsv",  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# get intersection of peaks from both comparisons - this is the persister signature peaks
NACTNaiv_P <- read.table("daps_NACT_vs_Naive.tsv", sep = "\t", header = TRUE) %>% drop_na() %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% pull(peak)
ResiSens_P <- read.table("daps_Resistant_vs_Sensitive.tsv", sep = "\t", header = TRUE) %>% drop_na() %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% pull(peak)
length(setdiff(NACTNaiv_P, ResiSens_P))
length(base::intersect(NACTNaiv_P, ResiSens_P))
length(setdiff(ResiSens_P, NACTNaiv_P))
patient_PCS_peaks <- base::intersect(NACTNaiv_P, ResiSens_P)
write.table(data.frame(peak = patient_PCS_peaks) %>% separate(peak, c("chr", "start", "end"), sep = "-") %>% arrange(chr, as.numeric(start)), "patient_PCS_peaks.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# part E: run ReMapEnrich on the EPS peaks to identify enriched DNA binding factors
library(ReMapEnrich)
library(GenomicRanges)
library(Signac)
# we are using a custom remap catalog file, but the catalog from the ReMap2022 database will work fine: https://remap.univ-amu.fr/download_page
remapCatalog <- bedToGranges("myremap.bed")
set.seed(1)
length(patient_PCS_peaks)
en <- enrichment(patient_PCS_peaks %>% StringToGRanges(), remapCatalog, chromSizes = loadChromSizes("hg38"), byChrom = FALSE)
plotdata <- en %>% mutate(category = factor(category, levels = en$category)) %>% mutate(rn = row_number())
patient_EPS_factors <- plotdata %>% pull(category) %>% head(n = 100) %>% as.character() %>% sort()
TFs_unique <- plotdata %>% 
  pull(category) %>% 
  as.character() %>%
  sub(":.*", "", .) %>%   # strip cell line
  unique() %>%            # remove duplicates (preserves rank order)
  head(120)               # take top 120 by significance
write.table(TFs_unique, "PCS_ReMapEnrich.csv", row.names = FALSE, col.names = FALSE)

# validation with published data
pcs_paper_df <- read.csv("PCS_paper_DBF_remap_results_100TFs.csv", header = FALSE)
pcs_paper_factors <- pcs_paper_df %>% pull(V1)
length(setdiff(TFs_unique, pcs_paper_factors))
length(base::intersect(TFs_unique, pcs_paper_factors))


####part F: visualization####
library(ggplot2)
library(cowplot)
temp <- data.frame(en) %>% drop_na() %>% arrange(`q.significance`)
plotdata2 <- temp %>%
  mutate(category = factor(category, levels = temp$category)) %>%
  mutate(TF_name = sub(":.*", "", as.character(category))) %>%
  mutate(issig = case_when(
    TF_name %in% TFs_unique ~ TRUE,  
    TRUE ~ NA                          # everything else → black
  ))

ggplot(data = plotdata2 %>% arrange(issig), 
       mapping = aes(x = category, y = `q.significance`, color = issig)) +
  ggrastr::rasterise(geom_point(size = 2, stroke = 0.1, shape = 16), dpi = 400) +
  scale_color_manual(
    values = c("red"), 
    na.value = "black",
    labels = c("TRUE" = "Top 100 DBFs", "NA" = "Other DBFs"),  
    name = NULL  # no legend title
  ) +
  scale_x_discrete(expand = c(0.1, 0.1)) +
  geom_hline(yintercept = 450, linetype = "dashed", color = "red", linewidth = 0.5) +
  labs(
    x = "DNA binding factors", 
    y = expression(-log[10](adj.~p-value)),
    title = "Persister cell signature (PCS)"
  ) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  theme_cowplot() + 
  theme(
    axis.text.x = element_blank(), 
    axis.ticks.x = element_blank(), 
    legend.position = c(0.05, 0.85),
    plot.title = element_text(hjust = 0.5)  
  )

