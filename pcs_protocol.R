library(BPCells)
library(readr)
library(BSgenome)  
library(GenomicRanges)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
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
library(org.Hs.eg.db)
library(BiocParallel)
library(ggpubr)
library(ggplot2)
library(cowplot)
library(patchwork)

register(SerialParam()) 
register(SnowParam(workers = 4))
all_genes <- keys(org.Hs.eg.db, keytype = "SYMBOL")

#### part A: read data and prep Seurat object ####
mat = open_matrix_10x_hdf5("data/counts_10x.h5", feature_type = "Peaks")
md = read_csv("data/metadata.csv")
chrom_assay = CreateChromatinAssay(counts = mat, sep = c("-", "-"))
obj = CreateSeuratObject(counts = chrom_assay, assay = "ATAC", genome = "hg38")
obj$celltype = md$celltype
obj$grouping1 = md$grouping1
obj$grouping2 = md$grouping2
obj = obj %>% RunTFIDF() %>% FindTopFeatures(min.cutoff = 20) %>% RunSVD()
saveRDS(obj, "data/seurat_obj.rds")

#### part B: Prepare ReMap data ####
# download remap BED file from https://remap.univ-amu.fr/download_page and extract it (outside R) and place in data folder
# cd data
# gunzip remap2022_nr_macs2_hg38_v1_0.bed.gz
# cut -d":" -f1 remap2022_nr_macs2_hg38_v1_0.bed > remap2022_nr_macs2_hg38_v1_0_MOD.bed
# grep -E '^chr([0-9]+|X|Y|M)[[:space:]]' remap2022_nr_macs2_hg38_v1_0_MOD.bed > remap2022_nr_macs2_hg38_v1_0_MOD_clean.bed

# OPTIONAL: additions to remap BED file (peaks in addons.bed)
# cat remap2022_nr_macs2_hg38_v1_0_MOD_clean.bed addons.bed | sort -k1,1 -k2,2n > myremap.bed

#### part C: create remap assay (OPTIONAL - if not doing this, continue with filtered <- readRDS("seurat_obj.rds") in part D) ####
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
obj <- readRDS("data/seurat_obj.rds")

# export peaks BED from atac assay
data.frame(Peak = rownames(obj)) %>%
  separate(col = "Peak", sep = "-", into = c("chr", "start", "end")) %>%
  write.table("multiome_peaks.bed", sep="\t", quote=FALSE, 
              row.names=FALSE, col.names=FALSE)

# run the following command outside R
# bedtools intersect -a multiome_peaks.bed -b data/myremap.bed -wa -wb > peaks_x_remap.txt

hits <- fread("peaks_x_remap.txt", header = FALSE)
# Cols 1-3 are your peak coords, col 7 is the TF name from remap
hits <- hits[, .(peak = paste(V1, V2, V3, sep = "-"), TF = V7)]

hits <- hits[grepl("^chr([0-9]+|X|Y|M)$", sub("-.*", "", hits$peak))]

# build the named list and drop any TFs with no overlapping peaks
TF_features <- split(hits$peak, paste0("remap_", hits$TF))
TF_features <- TF_features[lengths(TF_features) > 0]

# add TF_features as an assay in the seurat object
obj <- AddChromatinModule(obj, features = TF_features, 
                               genome = genome, assay = "ATAC")
remapscores <- obj@meta.data %>% select(starts_with("remap-")) %>% t()
rownames(remapscores) <- str_sub(rownames(remapscores), 7L, -1L)
obj[["remap"]] <- CreateAssayObject(data = remapscores)
obj@meta.data <- obj@meta.data %>% select(!starts_with("remap-"))
saveRDS(obj, "data/seurat_obj_remap.rds")


#### part D: defining the PCS from the atac assay####
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(GenomicRanges)
library(tidyverse)
# workaround to avoid seqinfo inherited method error - ignore if you don't have the error!
setMethod("seqinfo", signature(x = "ChromatinAssay"), function(x) {
  GenomeInfoDb::seqinfo(x@ranges)
})

# load seurat object
obj <- readRDS("data/seurat_obj_remap.rds")
DefaultAssay(obj) <- "ATAC"
# subset to celltype of interest
Tumor_Epi <- subset(obj, subset = celltype == "Epithelial_cells")

# get differentially accessible peaks for both comparisons
daps_NaiveNACT <- FindMarkers(Tumor_Epi, ident.1 = "NACT", ident.2 = "Naive",min.pct = 0 , group.by = "grouping1") %>% rownames_to_column("peak")
write.table(daps_NaiveNACT, "results/daps_NACT_vs_Naive.tsv",  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

daps_SensRes <- FindMarkers(Tumor_Epi, ident.1 = "Resistant", ident.2 = "Sensitive",min.pct = 0 , group.by = "grouping2") %>% rownames_to_column("peak")
write.table(daps_SensRes, "results/daps_Resistant_vs_Sensitive.tsv",  sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)

# get intersection of peaks from both comparisons - this is the persister signature peaks
NACTNaiv_P <- read.table("results/daps_NACT_vs_Naive.tsv", sep = "\t", header = TRUE) %>% drop_na() %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% pull(peak)
ResiSens_P <- read.table("results/daps_Resistant_vs_Sensitive.tsv", sep = "\t", header = TRUE) %>% drop_na() %>% dplyr::filter(p_val_adj < 0.05 & avg_log2FC > 0) %>% pull(peak)
length(setdiff(NACTNaiv_P, ResiSens_P))
length(base::intersect(NACTNaiv_P, ResiSens_P))
length(setdiff(ResiSens_P, NACTNaiv_P))
patient_PCS_peaks <- base::intersect(NACTNaiv_P, ResiSens_P)
write.table(data.frame(peak = patient_PCS_peaks) %>% separate(peak, c("chr", "start", "end"), sep = "-") %>% arrange(chr, as.numeric(start)), "results/patient_PCS_peaks.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

#### part E: run ReMapEnrich on the EPS peaks to identify enriched DNA binding factors ####
library(ReMapEnrich)
library(GenomicRanges)
library(Signac)

remapCatalog <- bedToGranges("data/myremap.bed")
set.seed(1)
en <- enrichment(patient_PCS_peaks %>% StringToGRanges(), remapCatalog, chromSizes = loadChromSizes("hg38"), byChrom = FALSE)
temp <- data.frame(en) %>% drop_na() %>% arrange(-`q.significance`) %>% filter(category %in% all_genes)
plotdata <- temp %>% mutate(category = factor(category, levels = temp$category)) %>% mutate(rn = row_number())
patient_PCS_factors <- plotdata %>% pull(category)  %>% as.character()
patient_PCS_factors_120 <- plotdata %>% pull(category) %>% head(n = 120) %>% as.character()
write.table(patient_PCS_factors_120, "results/PCS_ReMapEnrich.csv", row.names = FALSE, col.names = FALSE)
# validation with published data
pcs_paper_df <- read.csv("data/PCS_paper_results.csv", header = TRUE)
pcs_paper_factors <- pcs_paper_df %>% pull(PCS)
length(setdiff(patient_PCS_factors_120, pcs_paper_factors))
length(base::intersect(patient_PCS_factors_120, pcs_paper_factors))


####part F: optional RNA section####

mat = open_matrix_10x_hdf5("data/counts_10x_RNA.h5")
md = read_csv("data/metadata.csv")
rna_obj <- CreateSeuratObject(counts = mat, assay = "RNA")
rna_obj$celltype = md$celltype
rna_obj$grouping1 = md$grouping1
rna_obj$grouping2 = md$grouping2
Tumor_Epi_rna <- subset(rna_obj, subset = celltype == "Epithelial_cells")


# Calculate the percentage of NACT tumor cells that express each gene
NACT_perc_expr <- enframe(rowSums(Tumor_Epi_rna$RNA$counts[,which(Tumor_Epi_rna$grouping1 == "NACT")] > 0) * 100 / sum(Tumor_Epi_rna$grouping1 == "NACT")) # WMI
# Get top 100 expressed PCS factors from the top 120 chromatin-enriched PCS factors
patient_PCS_factors_rna <- NACT_perc_expr %>% filter(name %in% patient_PCS_factors_120) %>% arrange(-value) %>% head(n = 100) %>% pull(name)

# validation with published data
pcs_paper_df <- read.csv("data/PCS_paper_results.csv", header = TRUE)
pcs_paper_factors <- pcs_paper_df %>% pull(PCS)
length(setdiff(patient_PCS_factors_rna, pcs_paper_factors))
length(base::intersect(patient_PCS_factors_rna, pcs_paper_factors))



####part F: visualization####
library(ggplot2)
library(cowplot)
library(patchwork)

#plot 1
temp <- data.frame(en) %>% drop_na() %>% arrange(`q.significance`) %>% filter(category %in% all_genes)
plotdata <- temp %>%
  mutate(category = factor(category, levels = temp$category)) %>%
  mutate(issig = case_when(
    category %in% patient_PCS_factors_rna ~ "PCS",
    category %in% patient_PCS_factors_120 ~ "Low expressed",
    TRUE ~ "Other"
  )) %>%
  mutate(issig = factor(issig, levels = c("Other", "Low expressed", "PCS")))  # PCS last = drawn on top

p1 <- ggplot(data = plotdata %>% arrange(issig),  # Other first, PCS last (on top)
             mapping = aes(x = category, y = `q.significance`, color = issig)) +
  ggrastr::rasterise(geom_point(size = 2, stroke = 0.1, shape = 16), dpi = 400) +
  scale_color_manual(
    values = c("Other" = "black", "Low expressed" = "grey60", "PCS" = "red"),
    breaks = c("PCS", "Low expressed"),  # only these two appear in legend
    name = NULL
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
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "white", color = NA),  
    panel.background = element_rect(fill = "white", color = NA) 
  )

#plot 2

tbl <- read.table(file = "daps_NACT_vs_Naive.tsv", sep = "\t", skip = 1) %>%
  set_colnames(c("peak", "p_val", "effect_size", "pct.1", "pct.2", "p_val_adj"))
temp1 <- tbl %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
  mutate(neg_log10_adj_pval = -log10(p_val_adj))

p2 <- ggplot(temp1, aes(x = effect_size, y = neg_log10_adj_pval)) +
  ggrastr::rasterise(geom_point(data = temp1, color = "#e78ac3", size = 1.5, stroke=0.1, shape = 16), dpi = 400) +
  annotate("rect", xmin = -2.8, xmax = -0.1, ymin = 420, ymax = 440,
           fill = "#2E8B8B", color = NA) +
  annotate("text", x = -1.45, y = 430, label = "Naïve",
           color = "white", size = 3.5, fontface = "bold") +
  annotate("rect", xmin =  0.1, xmax =  5.5, ymin = 420, ymax = 440,
           fill = "#1C5F6B", color = NA) +
  annotate("text", x =  2.8,  y = 430, label = "NACT",
           color = "white", size = 3.5, fontface = "bold") +
  scale_x_continuous(
    limits = c(-2.8, 5.5),
    breaks = c(-2.5, 0, 2.5, 5.0)
  ) +
  theme(axis.title = element_blank()) +
  labs(
    title = "Differentially accessible peaks",
    x     = expression("Avg" ~-log[2] ~ "(adj. p-value)"),
    y     = expression(-log[10] ~ "(adj. p-value)")
  )+
  theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.line    = element_line(color = "black"),
    panel.grid   = element_blank()
  )

## plot 3 

tbl <- read.table(file = "daps_Resistant_vs_Sensitive.tsv", sep = "\t", skip = 1) %>%
  set_colnames(c("peak", "p_val", "effect_size", "pct.1", "pct.2", "p_val_adj"))
temp1 <- tbl %>% 
  filter(p_val_adj < 0.05) %>% 
  mutate(p_val_adj = ifelse(p_val_adj == 0, 1e-305, p_val_adj)) %>% 
  mutate(neg_log10_adj_pval = -log10(p_val_adj))

p3 <- ggplot(temp1, aes(x = effect_size, y = neg_log10_adj_pval)) +
  ggrastr::rasterise(geom_point(data = temp1, color = "#e78ac3", size = 1.5, stroke=0.1, shape = 16), dpi = 400) +
  annotate("rect", xmin = -2.8, xmax = -0.1, ymin = 420, ymax = 440,
           fill = NA, color ="#128ecc" ) +
  annotate("text", x = -1.45, y = 430, label = "Sensitive",
           color = "#128ecc", size = 3.5, fontface = "bold") +
  annotate("rect", xmin =  0.1, xmax =  5.5, ymin = 420, ymax = 440,
           fill = NA, color = "#a62323") +
  annotate("text", x =  2.8,  y = 430, label = "Resistant",
           color = "#a62323", size = 3.5, fontface = "bold") +
  scale_x_continuous(
    limits = c(-2.8, 5.5),
    breaks = c(-2.5, 0, 2.5, 5.0)
  ) +
  theme(axis.title = element_blank()) +
  labs(
    title = "Differentially accessible peaks",
    x     = expression("Avg" ~-log[2] ~ "(adj. p-value)"),
    y     = expression(-log[10] ~ "(adj. p-value)")
  )+
  theme_classic(base_size = 12) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    axis.line    = element_line(color = "black"),
    panel.grid   = element_blank()
  )


panel <- plot_grid(
  p2, p3, p1,
  ncol = 3,
  align = "h",
  axis = "b",
  labels = c("A", "B", "C"),
  label_size = 12
)
ggsave("plots/figure_panel.tiff", plot = panel, width = 14, height = 5, units = "in")

#### part H: Add PCS as a module score for each cell####
Tumor_Epi <- AddModuleScore(
  object = Tumor_Epi,
  features = list(patient_PCS_factors_rna),
  assay = "remap",
  name = 'PCS_factors',
  ctrl = 10,
  seed = 123
)
Tumor_Epi$persist <- case_when(Tumor_Epi$PCS_factors1 > quantile(Tumor_Epi$PCS_factors1)[4] ~ "High", Tumor_Epi$PCS_factors1 < quantile(Tumor_Epi$PCS_factors1)[2] ~ "Low", TRUE ~ "Mid")

saveRDS(Tumor_Epi, "data/TumorEpi_PCSscore.rds")
Tumor_Epi <- readRDS("data/TumorEpi_PCSscore.rds")
##
plotdata <- Tumor_Epi@meta.data %>% 
  dplyr::select(persist, orig.ident, grouping2) %>%
  group_by(persist, orig.ident, grouping2) %>%
  mutate(persist = factor(persist, levels = c("Low", "Mid", "High"))) %>%
  summarize(n = n())


table1 <- Tumor_Epi@meta.data %>% dplyr::select(orig.ident, grouping2, PCS_factors1) %>%
  filter(PCS_factors1 != Inf & PCS_factors1 != -Inf)

#violin plot
p1 <- ggplot(table1, aes(x=grouping2, y=PCS_factors1)) + 
  geom_violin(colour = "black", fill = "grey80") +
  stat_compare_means(aes(label = ..p.signif..), comparisons = list(c("Sensitive", "Resistant"), c("Resistant", "NACT"))) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  stat_summary(fun=mean, geom="point", size=1, color="black") +
  labs(y = "PCS score", title = "PCS Score") +
  theme_cowplot() + 
  theme(legend.position = "none", axis.title.x = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        plot.title = element_text(hjust = 0.5))
#proportion plot
p2 <- ggplot(plotdata, aes(fill=persist, y=n, x=orig.ident)) +
  geom_bar(position="fill", stat="identity", color = "black") + 
  facet_grid(cols = vars(grouping2), scales = "free", space = "free") +
  scale_fill_manual(values = c(Low = "#fdb515", Mid = "#cccccb", High = "#3953a4"), name = "PCS Score") +
  theme_cowplot() + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1), axis.title = element_blank())


combined <- plot_grid(p1, p2, ncol = 2, rel_widths = c(0.6, 1))

pdf("plots.pdf", width = 8, height = 6)
print(combined)
ggsave("plots/figure_panel2.tiff", plot = combined, width = 10, height = 5, units = "in")
print(combined)
dev.off()