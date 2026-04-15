# Figure 6-17 Code Reference

This document pulls together the recovered and reconstructed R code used for thesis Figures 6-17.


## Common Setup
open R
```r
library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(readxl)

ax <- readRDS("Axolotl.rds")
ps <- readRDS("Polypterus.rds")
zf <- readRDS("Zebrafish.rds")
zf_map <- read_excel("zebrafish_genenames_geneIDs.xlsx")

if (!dir.exists("Plots")) dir.create("Plots")
if (!dir.exists("Thesis Stuff")) dir.create("Thesis Stuff")
```

## Figure 6. Axolotl Heatmap

This is the Axolotl heatmap with the verified gene ordering and cluster ordering used in the remade figure.

```r
source("remake_species_heatmaps.R")

# Re-save Axolotl only if needed
save_heatmap(ax_result, "Axolotl_CellTypeGeneHeatmap_scaledSpectrum_VERIFIED")
```

## Figure 7. Polypterus Heatmap

This is the Polypterus heatmap with `IGF1R = LOC120541537`, `TSC1 = LOC120536016`, and the verified ordering.

```r
source("remake_species_heatmaps.R")

# Re-save Polypterus only if needed
save_heatmap(ps_result, "Polypterus_CellTypeGeneHeatmap_scaledSpectrum_VERIFIED")
```

## Figure 8. Zebrafish Heatmap

This is the Zebrafish heatmap with the verified reordered y-axis and zebrafish-specific gene IDs.

```r
source("remake_species_heatmaps.R")

# Re-save Zebrafish only if needed
save_heatmap(zf_result, "Zebrafish_CellTypeGeneHeatmap_scaledSpectrum_reordered_VERIFIED")
```

## Helper for Figures 9-11. Top Enriched Cell Type Per Gene

Recovered from `.Rhistory` and cleaned into one reusable helper.

```r
top_enriched_celltype <- function(obj, genes, species,
                                  group_col = "seurat_clusters",
                                  assay = "RNA", slot = "data") {
  genes_present <- intersect(genes, rownames(obj))
  genes_missing <- setdiff(genes, genes_present)

  if (length(genes_present) == 0) {
    stop(paste0(species, ": none of the requested genes were found in rownames(obj)."))
  }

  avg_mat <- AverageExpression(
    obj,
    assays = assay,
    features = genes_present,
    group.by = group_col,
    slot = slot
  )[[assay]]

  df <- as.data.frame(avg_mat) |>
    rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "celltype", values_to = "avg_expr") |>
    group_by(gene) |>
    slice_max(order_by = avg_expr, n = 1, with_ties = FALSE) |>
    ungroup() |>
    mutate(species = species)

  list(
    top_table = df,
    missing_genes = genes_missing,
    avg_matrix = avg_mat
  )
}

plot_top_table <- function(df, species_title) {
  df <- df %>% arrange(avg_expr)

  ggplot(df, aes(x = avg_expr, y = gene, fill = celltype)) +
    geom_col() +
    theme_classic() +
    labs(
      title = paste0(species_title, " - Top enriched cell type per gene"),
      x = "Average expression (slot=data)",
      y = "Gene",
      fill = "Top cell type"
    )
}
```

## Figure 9. Axolotl Top Enriched Cell Type Per Gene

```r
genes_ax <- c(
  "IGF1R", "IGF2R", "AKT1", "AKT2", "TSC1", "TSC2",
  "RHEB", "MTOR", "RPTOR", "LAMTOR1", "LAMTOR2", "LAMTOR3", "LAMTOR4", "LAMTOR5"
)

res_ax <- top_enriched_celltype(ax, genes_ax, "Axolotl")
p_ax_top <- plot_top_table(res_ax$top_table, "Axolotl")

ggsave("Thesis Stuff/Axolotl_TopEnrichedCellType_perGene.png", p_ax_top, width = 9, height = 5, dpi = 300)
ggsave("Thesis Stuff/Axolotl_TopEnrichedCellType_perGene.pdf", p_ax_top, width = 9, height = 5)
```

## Figure 10. Polypterus Top Enriched Cell Type Per Gene

```r
genes_ps <- c(
  "LOC120541537", "igf2r", "akt1s1", "akt2", "LOC120536016", "tsc2",
  "rheb", "mtor", "rptor", "lamtor1", "lamtor2", "lamtor3", "lamtor4", "lamtor5"
)

res_ps <- top_enriched_celltype(ps, genes_ps, "Polypterus")
res_ps$top_table <- res_ps$top_table %>%
  mutate(
    gene = recode(
      gene,
      "LOC120541537" = "IGF1R",
      "LOC120536016" = "TSC1",
      .default = toupper(gene)
    )
  )

p_ps_top <- plot_top_table(res_ps$top_table, "Polypterus")

ggsave("Thesis Stuff/Polypterus_TopEnrichedCellType_perGene.png", p_ps_top, width = 9, height = 5, dpi = 300)
ggsave("Thesis Stuff/Polypterus_TopEnrichedCellType_perGene.pdf", p_ps_top, width = 9, height = 5)
```

## Figure 11. Zebrafish Top Enriched Cell Type Per Gene

```r
genes_zf_symbols <- c(
  "igf1ra", "igf1rb", "igf2r", "akt1", "akt2",
  "tsc1a", "tsc1b", "tsc2", "rheb", "mtor", "rptor",
  "lamtor1", "lamtor2", "lamtor3", "lamtor4", "lamtor5"
)

symbol_to_ens <- setNames(zf_map$Gene, toupper(zf_map$gene_names))
genes_zf_ens <- unname(symbol_to_ens[toupper(genes_zf_symbols)])
genes_zf_ens <- genes_zf_ens[!is.na(genes_zf_ens)]

res_zf_raw <- top_enriched_celltype(zf, genes_zf_ens, "Zebrafish")
ens_to_symbol <- setNames(zf_map$gene_names, zf_map$Gene)

res_zf <- res_zf_raw
res_zf$top_table <- res_zf$top_table %>%
  mutate(gene = toupper(ifelse(gene %in% names(ens_to_symbol), ens_to_symbol[gene], gene)))

p_zf_top <- plot_top_table(res_zf$top_table, "Zebrafish")

ggsave("Thesis Stuff/Zebrafish_TopEnrichedCellType_perGene.png", p_zf_top, width = 9, height = 5, dpi = 300)
ggsave("Thesis Stuff/Zebrafish_TopEnrichedCellType_perGene.pdf", p_zf_top, width = 9, height = 5)
```

## Figure 12. Cross-Species Comparison of mTOR Pathway Gene Localization

This cross-species comparison uses the top-enriched cell type tables from Figures 9-11.

```r
ax_summary <- res_ax$top_table |>
  count(celltype, name = "n_genes") |>
  mutate(species = "Axolotl", proportion = n_genes / sum(n_genes))

ps_summary <- res_ps$top_table |>
  count(celltype, name = "n_genes") |>
  mutate(species = "Polypterus", proportion = n_genes / sum(n_genes))

zf_summary <- res_zf$top_table |>
  count(celltype, name = "n_genes") |>
  mutate(species = "Zebrafish", proportion = n_genes / sum(n_genes))

cross_species_df <- bind_rows(ax_summary, ps_summary, zf_summary)

p_cross <- ggplot(cross_species_df, aes(x = species, y = proportion, fill = celltype)) +
  geom_col(position = "fill") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_classic() +
  labs(
    title = "Cross-Species Comparison of mTOR Pathway Gene Localization",
    x = NULL,
    y = "Proportion of pathway genes",
    fill = "Top-enriched cell type"
  )

ggsave("CrossSpecies_mTOR_Localization_BarPlot.png", p_cross, width = 10, height = 6, dpi = 300)
ggsave("CrossSpecies_mTOR_Localization_BarPlot.pdf", p_cross, width = 10, height = 6)
```

## Helper for Figures 13-17. Split UMAPs by Days Post-Amputation

This helper produces one split UMAP across the day/condition column for a given species and gene. It is based on the `FeaturePlot(..., split.by = "Condition")` pattern already present in your scripts.

```r
plot_split_umap <- function(obj, feature_id, output_stem,
                            split_col = "Condition",
                            low_color = "lightgrey",
                            high_color = "magenta") {
  p <- FeaturePlot(
    obj,
    features = feature_id,
    split.by = split_col,
    cols = c(low_color, high_color)
  )

  ggsave(paste0(output_stem, ".png"), p, width = 14, height = 5, dpi = 300)
  ggsave(paste0(output_stem, ".pdf"), p, width = 14, height = 5)
}
```

## Figure 13. AKT2 UMAPs Across Days Post-Amputation

```r
# Axolotl
plot_split_umap(ax, "AKT2", "Thesis Stuff/Axolotl_AKT2_UMAP")

# Polypterus
plot_split_umap(ps, "akt2", "Thesis Stuff/Polypterus_akt2_UMAP")

# Zebrafish
plot_split_umap(zf, "ENSDARG00000011219", "Thesis Stuff/Zebrafish_ENSDARG00000011219_akt2_UMAP")
```

## Figure 14. RPTOR UMAPs Across Days Post-Amputation

```r
# Axolotl
plot_split_umap(ax, "RPTOR", "Thesis Stuff/Axolotl_RPTOR_UMAP")

# Polypterus
plot_split_umap(ps, "rptor", "Thesis Stuff/Polypterus_rptor_UMAP")

# Zebrafish
plot_split_umap(zf, "ENSDARG00000098726", "Thesis Stuff/Zebrafish_ENSDARG00000098726_rptor_UMAP")
```

## Figure 15. LAMTOR3 UMAPs Across Days Post-Amputation

```r
# Axolotl
plot_split_umap(ax, "LAMTOR3", "Thesis Stuff/Axolotl_LAMTOR3_UMAP")

# Polypterus
plot_split_umap(ps, "lamtor3", "Thesis Stuff/Polypterus_lamtor3_UMAP")

# Zebrafish
plot_split_umap(zf, "ENSDARG00000057075", "Thesis Stuff/Zebrafish_ENSDARG00000057075_lamtor3_UMAP")
```

## Figure 16. LAMTOR4 UMAPs Across Days Post-Amputation

```r
# Axolotl
plot_split_umap(ax, "LAMTOR4", "Thesis Stuff/Axolotl_LAMTOR4_UMAP")

# Polypterus
plot_split_umap(ps, "lamtor4", "Thesis Stuff/Polypterus_lamtor4_UMAP")

# Zebrafish
plot_split_umap(zf, "ENSDARG00000045542", "Thesis Stuff/Zebrafish_ENSDARG00000045542_lamtor4_UMAP")
```

## Figure 17. IGF1R UMAPs Across Days Post-Amputation

```r
# Axolotl
plot_split_umap(ax, "IGF1R", "Thesis Stuff/Axolotl_IGF1R_UMAP")

# Polypterus
plot_split_umap(ps, "LOC120541537", "Thesis Stuff/Polypterus_IGF1R_UMAP")

# Zebrafish
plot_split_umap(zf, "ENSDARG00000027423", "Thesis Stuff/Zebrafish_ENSDARG00000027423_igf1r_UMAP")
```

## Notes

```r
# Axolotl feature names are direct gene symbols for these pathway genes.
# Polypterus uses mixed lowercase symbols plus LOC IDs for some orthologs.
# Zebrafish uses ENSDARG IDs in the object and needs mapping via zebrafish_genenames_geneIDs.xlsx.
#
# Confirm these special mappings when reproducing figures:
# Polypterus IGF1R = LOC120541537
# Polypterus TSC1  = LOC120536016
# Polypterus RICTOR = LOC120527345
```
