DEseq2_function <- function(abundance_table, metadata, col_of_interest, id_col = "pathway"){
  
  # Check if the column exists in metadata
  if(!col_of_interest %in% colnames(metadata)){
    stop(paste("Column", col_of_interest, "not found in metadata."))
  }
  
  # Create copies of data
  DESeq2_metadata <- metadata
  DESeq2_abundance_mat <- abundance_table %>% column_to_rownames(id_col)
  
  # Rename the column of interest to "Group"
  colnames(DESeq2_metadata)[colnames(DESeq2_metadata) == col_of_interest] <- "Group"
  DESeq2_metadata$Group <- as.factor(DESeq2_metadata$Group)
  
  # Generate all pairwise combinations of groups (as a list)
  group_levels <- unique(DESeq2_metadata$Group)
  if(length(group_levels) < 2){
    stop("There are fewer than 2 groups in the Group column.")
  }
  group_combinations <- utils::combn(group_levels, 2, simplify = FALSE)
  
  results_list <- list()
  
  message("Performing pairwise comparisons with DESeq2...")
  for(comb in group_combinations){
    # Subset metadata and counts for samples in this combination
    keep <- DESeq2_metadata$Group %in% comb
    sub_metadata <- DESeq2_metadata[keep, ]
    sub_counts <- DESeq2_abundance_mat[, keep, drop = FALSE]
    sub_counts <- round(sub_counts)
    
    dds <- DESeq2::DESeqDataSetFromMatrix(
      countData = sub_counts,
      colData = sub_metadata,
      design = ~ Group
    )
    dds <- BiocGenerics::estimateSizeFactors(dds, type = "poscounts")
    dds <- DESeq2::DESeq(dds)
    res <- DESeq2::results(dds)
    
    # Convert results to data frame and add identifier columns
    res_df <- as.data.frame(res)
    res_df$Comparison <- paste(comb, collapse = "_vs_")
    res_df$feature <- rownames(res_df)
    
    results_list[[paste(comb, collapse = "_vs_")]] <- res_df
  }
  
  # Combine results from all comparisons
  combined_results <- dplyr::bind_rows(results_list)
  
  return(combined_results)
}