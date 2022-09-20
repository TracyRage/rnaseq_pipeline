# Transcriptome analysis utility functions
# Load necessary packages
library(QuasR)
library(rtracklayer)
library(EDASeq)
library(pheatmap)
library(RUVSeq)
library(DESeq2)
library(tidyverse)

# Column metadata cleaning function
clean_column_metadata <- function(filename) {
  col_data <- readr::read_tsv(filename) %>% 
    column_to_rownames("sample") %>% 
    janitor::clean_names()
}

# Clean gff input file
clean_gff <- function(annotation_file, format="gff3") {
  # Import gff data
  gtf_regions <- import.gff(annotation_file,
                            format=format)
  # Clean annotation data
  gene_df <- tibble(data.frame(gtf_regions)) %>% 
    select(gene) %>% drop_na() %>% pull()
  names(gtf_regions) <- gene_df
  # Return clean annotation file
  gtfRegionList <- split(gtf_regions, names(gtf_regions))
}

# Create gene product and counts tibble 
make_tibble <- function(metadata_matrix, counts_result) {
  # Render product tibble
  gene_df  <- tibble(data.frame(metadata_matrix)) %>% 
    select(gene, product) %>%
    drop_na() %>%
    distinct()
  # Render counts tibble
  counts_df <-  tibble(rownames_to_column(data.frame(counts_result), 
                                          var="gene"))
  return(list(gene_df=gene_df, counts_df=counts_df))
}

# Mapping function
func_align <- function(genome, seqs_file, aux) {
  proj <- qAlign(seqs_file, 
                 genome,
                 paired="fr",
                 auxiliaryFile = aux,
                 clObj = cl)
}

# Prepare data from EDA and normalization
prep_eda <- function(column_data, counts_matrix) {
  # remove 'width' column from result (counts) matrix
  count_data <- as.matrix(subset(result, select=c(-width)))
  
  # create a seqExpressionSet object using EDASeq package
  set <- newSeqExpressionSet(counts=count_data,
                             phenoData = column_data)
  return(set)
}

# Choose index for unwanted variance (k, RUVs)
pick_k <- function(column_data, eda_set) {
  # make a table of samples groups from column_data
  differences <- makeGroups(column_data$group)
  # pick appropriate k index
  par(mfrow=c(2,2))
  for (k in 1:4) {
    set_s <- RUVs(eda_set, unique(rownames(set)),
                  k=k, differences) # all genes
    DESeq2::plotPCA(set_s, col=as.numeric(column_data$group),
            cex=0.9, adj=0.5,
            main=paste0("with RUVs, k = ", k),
            ylim=c(-1,1), xlim=c(-0.6, 0.6))
  }
}

# Compare RLE plots
compare_RLE <- function(column_data, initial_set, rectified_set) {
  # Plot initial set
  plotRLE(initial_set, outline=FALSE, ylim=c(-4,4),
          col=as.numeric(column_data$group),
          main="without RUVs")
  # Plot rectified set
  plotRLE(rectified_set, outline=FALSE, ylim=c(-4,4),
          col=as.numeric(column_data$group),
          main="with RUVs")
}

# Compare PCA plots
compare_PCA <- function(column_data, initial_set, rectified_set) {
  # Plot initial set
  DESeq2::plotPCA(initial_set, 
          ylim=c(-0.75, 0.75),
          xlim=c(-0.75, 0.75),
          adj=0.5,
          col=as.numeric(column_data$group),
          main="without RUVs")
  # Plot rectified set
  DESeq2::plotPCA(rectified_set, 
          ylim=c(-0.75, 0.75),
          xlim=c(-0.75, 0.75),
          adj=0.5,
          col=as.numeric(column_data$group),
          main="without RUVs")
}

# Render heatmap
render_heatmap <- function(column_data, rectified_set) {
  # Extract normalized counts that are cleared from unwanted variation
  norm_count <- normCounts(rectified_set)
  selected_genes <- names(sort(apply(norm_count, 1, var),
                               decreasing = TRUE))[1:413]
  # Render heatmap
  pheatmap(norm_count[selected_genes,],
           annotation_col = subset(column_data, select=c(-group)),
           show_rownames = FALSE,
           annotation_names_col = FALSE,
           cutree_cols = 3,
           cutree_rows = 5,
           gaps_row = 10,
           angle_col = 45,
           scale = "row")
}

# Re-run DESeq2 with the computed covariates
re_run_deseq <- function(count_data, 
                         column_data, 
                         design,
                         rectified_set) {
  # set up DESeqDataSet object
  dds <- DESeqDataSetFromMatrix(countData = count_data,
                                colData = column_data,
                                design = design)
  
  # filter low count genes
  dds <- dds[rowSums(DESeq2::counts(dds)) > 10]
  
  # insert the covariates computed by RUVs
  colData(dds) <- cbind(colData(dds),
                        pData(set_s)[rownames(colData(dds)),
                                     grep("W_[0-9]",
                                          colnames(pData(rectified_set)))])
  
  # update the design formula (save interest variable for the last)
  design(dds) <- ~ W_1 + W_2 + group 
  
  # repeat the analysis
  dds <- DESeq(dds)
  
  # extract deseq results (PHE)
  res_PHE <- results(dds, contrast = c("group", "phe", "ctrl"))
  res_PHE <- res_PHE[order(res_PHE$padj < 0.05), ]
  PHE_df <- tibble(rownames_to_column(as.data.frame(res_PHE), var="gene"))
  
  res_PY <- results(dds, contrast = c("group", "py", "ctrl"))
  res_PY <- res_PY[order(res_PY$padj < 0.05), ]
  PY_df <- tibble(rownames_to_column(as.data.frame(res_PY), var="gene"))
  
  return(list(PHE=PHE_df, PY=PY_df, dds=dds, res_PHE=res_PHE, res_PY=res_PY))
  
}

# Render log tables w/ GO and KEGG annotations
bio_function <- function(go_file, pathway_file, log_file) {
  go_edited <- read_csv(go_file) %>% 
    clean_names() %>% 
    dplyr::select(gene_name, go_term, product_description, locus_tag) %>% 
    distinct(gene_name, .keep_all = TRUE) %>% 
    rename(gene = gene_name)
  
  go_annot <- left_join(log_file, go_edited, by="gene") %>% 
    drop_na() %>% 
    dplyr::select(-c("lfcSE", "stat", "pvalue", "baseMean"))
  
  pathways <- read_csv(pathway_file) %>% 
    clean_names() %>% 
    dplyr::select(locus_tag, pathway_name) %>% 
    distinct(locus_tag, .keep_all = T)
  
  go_pathway <- left_join(go_annot, pathways, by="locus_tag") %>% 
    drop_na() %>% 
    dplyr::select(-locus_tag) %>% 
    relocate(go_term, .after = product_description)
  
  go_path_up <- go_pathway %>% 
    filter(log2FoldChange > 0.5 & padj < 0.05)
  
  go_path_down <- go_pathway %>% 
    filter(log2FoldChange < -0.5 & padj < 0.05)
  
  go_path_both <- bind_rows(go_path_down, go_path_up)
  
  
  return(list(UP=go_path_up, DOWN=go_path_down, BOTH=go_path_both))
}

# Render final heatmap
render_final_hp <- function(rectified_set, column_data, row_data) {
  # Normalize gene counts
  norm_count <- normCounts(rectified_set)
  # Select genes w/ variance != 0
  selected_genes <- names(sort(discard(apply(norm_count, 1, var), ~.x == 0),
                               decreasing = TRUE))
  norm <- tibble(rownames_to_column(as.data.frame(norm_count[selected_genes, ]), "gene"))
  
  row_annot <- left_join(row_data, norm, by="gene") %>% 
    drop_na()
  
  # Normalized dataframe
  norm_2 <- column_to_rownames(row_annot %>% select(-pathway_name), "gene")
  # Rectified rownames annotation
  row_annot_2 <- column_to_rownames(row_annot %>% select(gene, pathway_name) %>% 
                                      rename(Pathway=pathway_name), "gene")
  
  colnames(column_data)[colnames(column_data) == "treatment"] <- "Treatment"
  
  ph <- pheatmap::pheatmap(norm_2,
                           annotation_col = subset(column_data, select=c(-group)),
                           annotation_row = row_annot_2,
                           annotation_names_col = FALSE,
                           scale="row",
                           cutree_rows = 15,
                           cutree_cols = 3,
                           gaps_row = 10,
                           gaps_col = 30,
                           angle_col = 45,
                           cellheight = 8,
                           cellwidth = 40,
                           fontsize_row = 8,
                           display_numbers = TRUE)
  return(ph)
}

# Render final heatmap
render_final_hp_RO <- function(rectified_set, column_data, row_data) {
  # Normalize gene counts
  norm_count <- normCounts(rectified_set)
  # Select genes w/ variance != 0
  selected_genes <- names(sort(discard(apply(norm_count, 1, var), ~.x == 0),
                               decreasing = TRUE))
  norm <- tibble(rownames_to_column(as.data.frame(norm_count[selected_genes, ]), "Gene"))
  
  print(norm)
  
  row_annot <- left_join(row_data, norm, by="Gene") %>% 
    drop_na()
  
  print(row_annot)
  
  # Normalized dataframe
  norm_2 <- column_to_rownames(row_annot %>% select(-`Căi metabolice`), "Gene")
  
  print(norm_2)
  # Rectified rownames annotation
  row_annot_2 <- column_to_rownames(row_annot %>% select(Gene, `Căi metabolice`) %>% 
                                      rename(`Căi metabolice`=`Căi metabolice`), "Gene")
  
  print(row_annot_2)
  
  colnames(column_data)[colnames(column_data) == "treatment"] <- "Tratament"
  
  ph <- pheatmap::pheatmap(norm_2,
                           annotation_col = subset(column_data, select=c(-group)),
                           annotation_row = row_annot_2,
                           annotation_names_col = FALSE,
                           scale="row",
                           cutree_rows = 15,
                           cutree_cols = 3,
                           gaps_row = 10,
                           gaps_col = 30,
                           angle_col = 45,
                           cellheight = 8,
                           cellwidth = 40,
                           fontsize_row = 8,
                           display_numbers = TRUE,
                           filename = "analysis_results/heatmap_ro.png")
  return(ph)
}

### Get waffle chart (in RO)
get_waffle <- function(pathway_counts) {
  pathway_counts %>% 
    rename(`Căi metabolice` = pathway_name) %>% 
    dplyr::count(`Căi metabolice`, sort = T, name="counts") %>% 
    arrange(desc(counts)) %>% 
    head(16) %>% 
    ggplot(aes(fill = `Căi metabolice`, values = counts)) +
    geom_waffle(color = "white", size=1.125, n_rows = 6) +
    coord_equal() +
    theme_enhance_waffle()
}
