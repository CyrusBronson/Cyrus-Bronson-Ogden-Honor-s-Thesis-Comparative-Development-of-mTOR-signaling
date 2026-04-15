# Figure 6-17 Code Reference

This document pulls together the recovered and reconstructed R code used for thesis Figures 6-17.

## Required Files and Setup Rules

To run this code successfully, the following files must be present in the working directory:

- `Axolotl.rds`
- `Polypterus.rds`
- `Zebrafish.rds`

Required R packages:

- `Seurat`
- `dplyr`
- `tidyr`
- `ggplot2`
- `tibble`

Important:

- Run the `Common Setup` block first.
- The setup block creates these object names:
  - Axolotl = `ax`
  - Polypterus = `ps`
  - Zebrafish = `zf`
- The split UMAP code assumes the metadata column `Condition` exists in the object.

## Zebrafish Translation Notes

The Zebrafish object uses `ENSDARG` identifiers rather than plain gene symbols. To make this code fully portable, the translations used in these figures are hardcoded directly below, so the Excel mapping file is not required.

Common translations used in these figures:

- `IGF1RA` = `ENSDARG00000027423`
- `IGF1RB` = `ENSDARG00000034434`
- `IGF2R` = `ENSDARG00000006094`
- `AKT1` = `ENSDARG00000099657`
- `AKT2` = `ENSDARG00000011219`
- `TSC1A` = `ENSDARG00000026048`
- `TSC1B` = `ENSDARG00000057918`
- `TSC2` = `ENSDARG00000103125`
- `RHEB` = `ENSDARG00000090213`
- `MTOR` = `ENSDARG00000053196`
- `RPTOR` = `ENSDARG00000098726`
- `LAMTOR1` = `ENSDARG00000076464`
- `LAMTOR2` = `ENSDARG00000039872`
- `LAMTOR3` = `ENSDARG00000057075`
- `LAMTOR4` = `ENSDARG00000045542`
- `LAMTOR5` = `ENSDARG00000090194`


## Common Setup

```r
install_if_missing <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing) > 0) install.packages(missing)
}

install_if_missing(c("Seurat", "dplyr", "tidyr", "ggplot2", "tibble"))

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)

ax <- readRDS("Axolotl.rds")
ps <- readRDS("Polypterus.rds")
zf <- readRDS("Zebrafish.rds")

if (!dir.exists("Plots")) dir.create("Plots")
if (!dir.exists("Thesis Stuff")) dir.create("Thesis Stuff")

zebrafish_id_map <- c(
  IGF1RA = "ENSDARG00000027423",
  IGF1RB = "ENSDARG00000034434",
  IGF2R = "ENSDARG00000006094",
  AKT1 = "ENSDARG00000099657",
  AKT2 = "ENSDARG00000011219",
  TSC1A = "ENSDARG00000026048",
  TSC1B = "ENSDARG00000057918",
  TSC2 = "ENSDARG00000103125",
  RHEB = "ENSDARG00000090213",
  MTOR = "ENSDARG00000053196",
  RPTOR = "ENSDARG00000098726",
  LAMTOR1 = "ENSDARG00000076464",
  LAMTOR2 = "ENSDARG00000039872",
  LAMTOR3 = "ENSDARG00000057075",
  LAMTOR4 = "ENSDARG00000045542",
  LAMTOR5 = "ENSDARG00000090194"
)
```

## Helper for Figures 6-8. Full Heatmap Code

This helper is the full runnable heatmap workflow used for the verified Axolotl, Polypterus, and Zebrafish heatmaps.

```r
zebrafish_palette <- c("#F3E6CF", "#F6D365", "#6FCF97", "#1F5AA6")

make_scaled_heatmap <- function(obj,
                                gene_map,
                                species,
                                cluster_order,
                                group_col = "seurat_clusters",
                                assay = "RNA",
                                slot = "data") {
  genes_present <- gene_map[gene_map %in% rownames(obj)]

  if (length(genes_present) == 0) {
    stop(paste0(species, ": none of the requested genes were found."))
  }

  avg_mat <- AverageExpression(
    obj,
    assays = assay,
    features = unname(genes_present),
    group.by = group_col,
    slot = slot
  )[[assay]]

  rownames(avg_mat) <- names(genes_present)[match(rownames(avg_mat), genes_present)]

  df_long <- as.data.frame(avg_mat) |>
    rownames_to_column("gene") |>
    pivot_longer(-gene, names_to = "celltype", values_to = "avg_expr") |>
    group_by(gene) |>
    mutate(
      avg_expr_scaled = {
        s <- sd(avg_expr)
        if (is.na(s) || s == 0) rep(0, dplyr::n()) else (avg_expr - mean(avg_expr)) / s
      }
    ) |>
    ungroup()

  df_long$gene <- factor(df_long$gene, levels = names(gene_map))

  cluster_levels <- cluster_order[cluster_order %in% unique(df_long$celltype)]
  remaining_levels <- setdiff(unique(df_long$celltype), cluster_levels)
  display_levels <- c(cluster_levels, remaining_levels)
  df_long$celltype <- factor(df_long$celltype, levels = rev(display_levels))

  ggplot(df_long, aes(x = celltype, y = gene, fill = avg_expr_scaled)) +
    geom_tile() +
    scale_fill_gradientn(
      colours = zebrafish_palette,
      name = "Scaled\navg expr"
    ) +
    labs(
      title = paste0(species, ": Mean expression by cell type (scaled per gene)"),
      x = "Cell type / cluster",
      y = "Gene"
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}
```

## Figure 6. Axolotl Heatmap

```r
axolotl_cluster_order <- c(
  "Superficial epidermis 1", "Superficial epidermis 2",
  "Intermediate epidermis 1", "Intermediate epidermis 2", "Intermediate epidermis 3",
  "Intermediate epidermis 4", "Intermediate epidermis 5",
  "Basal epidermis 1", "Basal epidermis 2",
  "Cluster 7", "Cluster 6",
  "Endothelial cell",
  "Connective tissue 1", "Connective tissue 2", "Connective tissue 3",
  "Platelet",
  "Erythrocyte 1", "Erythrocyte 2", "Erythrocyte 3",
  "Immune cell 1", "Immune cell 2", "Immune cell 3",
  "Cluster 22",
  "Myeloid cell 1", "Myeloid cell 2", "Myeloid cell 3"
)

axolotl_gene_map <- c(
  IGF1R = "IGF1R",
  IGF2R = "IGF2R",
  AKT1 = "AKT1",
  AKT2 = "AKT2",
  TSC1 = "TSC1",
  TSC2 = "TSC2",
  RHEB = "RHEB",
  MTOR = "MTOR",
  RPTOR = "RPTOR",
  LAMTOR1 = "LAMTOR1",
  LAMTOR2 = "LAMTOR2",
  LAMTOR3 = "LAMTOR3",
  LAMTOR4 = "LAMTOR4",
  LAMTOR5 = "LAMTOR5"
)

p_ax_heatmap <- make_scaled_heatmap(
  obj = ax,
  gene_map = axolotl_gene_map,
  species = "Axolotl",
  cluster_order = axolotl_cluster_order
)

ggsave("Thesis Stuff/Axolotl_CellTypeGene_HEATMAP.png", p_ax_heatmap + coord_flip(), width = 14, height = 9, dpi = 300)
ggsave("Thesis Stuff/Axolotl_CellTypeGene_HEATMAP.pdf", p_ax_heatmap + coord_flip(), width = 14, height = 9)
```

## Figure 7. Polypterus Heatmap

```r
polypterus_cluster_order <- c(
  "Proliferating cell 1", "Proliferating cell 2", "Proliferating cell 3",
  "Superficial epidermis 1", "Superficial epidermis 2", "Superficial epidermis 3", "Superficial epidermis 4",
  "Intermediate epidermis 1", "Intermediate epidermis 2", "Intermediate epidermis 3", "Intermediate epidermis 4",
  "Basal epidermis 1", "Basal epidermis 2",
  "Specialized epidermal cell", "Goblet cell 1", "Goblet cell 2",
  "Ionocyte 1", "Ionocyte 2", "Ionocyte 3",
  "Endothelial cell",
  "Connective tissue 1", "Connective tissue 2",
  "Platelet", "Erythrocyte",
  "Muscle 1", "Muscle 2",
  "Myeloid cell", "Granulocyte 1", "Granulocyte 2", "Immune cell", "Lymphoid cell", "Glial cell"
)

polypterus_gene_map <- c(
  IGF1R = "LOC120541537",
  IGF2R = "igf2r",
  AKT1 = if ("akt1" %in% rownames(ps)) "akt1" else "akt1s1",
  AKT2 = "akt2",
  TSC1 = "LOC120536016",
  TSC2 = "tsc2",
  RHEB = "rheb",
  MTOR = "mtor",
  RPTOR = "rptor",
  LAMTOR1 = "lamtor1",
  LAMTOR2 = "lamtor2",
  LAMTOR3 = "lamtor3",
  LAMTOR4 = "lamtor4",
  LAMTOR5 = "lamtor5"
)

p_ps_heatmap <- make_scaled_heatmap(
  obj = ps,
  gene_map = polypterus_gene_map,
  species = "Polypterus",
  cluster_order = polypterus_cluster_order
)

ggsave("Thesis Stuff/Polypterus_CellTypeGene_HEATMAP.png", p_ps_heatmap + coord_flip(), width = 14, height = 9, dpi = 300)
ggsave("Thesis Stuff/Polypterus_CellTypeGene_HEATMAP.pdf", p_ps_heatmap + coord_flip(), width = 14, height = 9)
```

## Figure 8. Zebrafish Heatmap

```r
zebrafish_cluster_order <- c(
  "Proliferating cell 1", "Proliferating cell 2", "Proliferating cell 3",
  "Superficial epidermis 1", "Superficial epidermis 2",
  "Intermediate epidermis 1", "Intermediate epidermis 2", "Intermediate epidermis 3",
  "Intermediate epidermis 4", "Intermediate epidermis 5", "Intermediate epidermis 6", "Intermediate epidermis 7",
  "Basal epidermis 1", "Basal epidermis 2", "Basal epidermis 3",
  "Specialized epidermal cell 1", "Specialized epidermal cell 2", "Epidermal mucous cell",
  "Connective tissue 1", "Connective tissue 2", "Connective tissue 3",
  "Hematopoietic cell 1 1", "Hematopoietic cell 2",
  "Myeloid cell", "Metaphocyte", "Pigment cell"
)

zebrafish_gene_map <- c(
  IGF1RA = zebrafish_id_map["IGF1RA"],
  IGF1RB = zebrafish_id_map["IGF1RB"],
  IGF2R = zebrafish_id_map["IGF2R"],
  AKT1 = zebrafish_id_map["AKT1"],
  AKT2 = zebrafish_id_map["AKT2"],
  TSC1A = zebrafish_id_map["TSC1A"],
  TSC1B = zebrafish_id_map["TSC1B"],
  TSC2 = zebrafish_id_map["TSC2"],
  RHEB = zebrafish_id_map["RHEB"],
  MTOR = zebrafish_id_map["MTOR"],
  RPTOR = zebrafish_id_map["RPTOR"],
  LAMTOR1 = zebrafish_id_map["LAMTOR1"],
  LAMTOR2 = zebrafish_id_map["LAMTOR2"],
  LAMTOR3 = zebrafish_id_map["LAMTOR3"],
  LAMTOR4 = zebrafish_id_map["LAMTOR4"],
  LAMTOR5 = zebrafish_id_map["LAMTOR5"]
)

zebrafish_gene_map <- zebrafish_gene_map[!is.na(zebrafish_gene_map)]

p_zf_heatmap <- make_scaled_heatmap(
  obj = zf,
  gene_map = zebrafish_gene_map,
  species = "Zebrafish",
  cluster_order = zebrafish_cluster_order
)

ggsave("Thesis Stuff/Zebrafish_CellTypeGene_HEATMAP.png", p_zf_heatmap + coord_flip(), width = 14, height = 9, dpi = 300)
ggsave("Thesis Stuff/Zebrafish_CellTypeGene_HEATMAP.pdf", p_zf_heatmap + coord_flip(), width = 14, height = 9)
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
genes_zf_symbols <- names(zebrafish_id_map)
genes_zf_ens <- unname(zebrafish_id_map[genes_zf_symbols])

res_zf_raw <- top_enriched_celltype(zf, genes_zf_ens, "Zebrafish")
ens_to_symbol <- setNames(names(zebrafish_id_map), zebrafish_id_map)

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

## Helper for Figures 13-17. Full Split UMAP Code

This helper produces one split UMAP across the day/condition column for a given species and gene. It is based on the `FeaturePlot(..., split.by = "Condition")` pattern already present in your scripts.

```r
plot_split_umap <- function(obj, feature_id, output_stem,
                            split_col = "Condition",
                            low_color = "lightgrey",
                            high_color = "magenta",
                            width = 14,
                            height = 5) {
  if (!(split_col %in% colnames(obj@meta.data))) {
    stop(paste0("Column '", split_col, "' was not found in metadata."))
  }

  p <- FeaturePlot(
    obj,
    features = feature_id,
    split.by = split_col,
    cols = c(low_color, high_color)
  )

  ggsave(paste0(output_stem, ".png"), p, width = width, height = height, dpi = 300)
  ggsave(paste0(output_stem, ".pdf"), p, width = width, height = height)
}
```

## Figure 13. AKT2 UMAPs Across Days Post-Amputation

```r
# Run this helper first:
# plot_split_umap(...)

# Axolotl
plot_split_umap(ax, "AKT2", "Thesis Stuff/Axolotl_AKT2_UMAP")

# Polypterus
plot_split_umap(ps, "akt2", "Thesis Stuff/Polypterus_akt2_UMAP")

# Zebrafish
plot_split_umap(zf, "ENSDARG00000011219", "Thesis Stuff/Zebrafish_ENSDARG00000011219_akt2_UMAP")
```

## Figure 14. RPTOR UMAPs Across Days Post-Amputation

```r
# Run this helper first:
# plot_split_umap(...)

# Axolotl
plot_split_umap(ax, "RPTOR", "Thesis Stuff/Axolotl_RPTOR_UMAP")

# Polypterus
plot_split_umap(ps, "rptor", "Thesis Stuff/Polypterus_rptor_UMAP")

# Zebrafish
plot_split_umap(zf, "ENSDARG00000098726", "Thesis Stuff/Zebrafish_ENSDARG00000098726_rptor_UMAP")
```

## Figure 15. LAMTOR3 UMAPs Across Days Post-Amputation

```r
# Run this helper first:
# plot_split_umap(...)

# Axolotl
plot_split_umap(ax, "LAMTOR3", "Thesis Stuff/Axolotl_LAMTOR3_UMAP")

# Polypterus
plot_split_umap(ps, "lamtor3", "Thesis Stuff/Polypterus_lamtor3_UMAP")

# Zebrafish
plot_split_umap(zf, "ENSDARG00000057075", "Thesis Stuff/Zebrafish_ENSDARG00000057075_lamtor3_UMAP")
```

## Figure 16. LAMTOR4 UMAPs Across Days Post-Amputation

```r
# Run this helper first:
# plot_split_umap(...)

# Axolotl
plot_split_umap(ax, "LAMTOR4", "Thesis Stuff/Axolotl_LAMTOR4_UMAP")

# Polypterus
plot_split_umap(ps, "lamtor4", "Thesis Stuff/Polypterus_lamtor4_UMAP")

# Zebrafish
plot_split_umap(zf, "ENSDARG00000045542", "Thesis Stuff/Zebrafish_ENSDARG00000045542_lamtor4_UMAP")
```

## Figure 17. IGF1R UMAPs Across Days Post-Amputation

```r
# Run this helper first:
# plot_split_umap(...)

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
