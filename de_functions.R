# de_functions.R
#
# A library of functions to analyze and visualize differential gene expression data
# Companion file to deseq_analysis.R
# BY Sami Friedrich
# CREATED 06/01/2021
# UPDATED 04/21/2022

library(tidyverse)
library(ggplot2)
library(DESeq2)
library(xlsx)
library(readr)
library(DEGreport)
library(grid)
library(gridExtra)
library(annotables)
library(org.Hs.eg.db)
library(clusterProfiler)

# set working directory in RStudio to wherever this R script file lives
working_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(working_dir)

# human genome
grch38 <- annotables::grch38

### Iitialize library of GO terms by gene
goannot <- AnnotationDbi::select(org.Hs.eg.db, keys=keys(org.Hs.eg.db), columns="GO")
genesbygo <- split(goannot$ENTREZID, goannot$GO)

# custom color palette for plotting
paired_palette <- c("#E69F00", "#D55E00", "#009E73", "#0072B2")

# group labels
group_labels <- c("Female 20 DPH", "Female 50 DPH", "Male 20 DPH", "Male 50 DPH")

# result labels
result_labels <- c("Female development", 
                   "Male development", 
                   "20 DPH sex contrast",
                   "50 DPH sex contrast",
                   "Sex+age interaction")
names(result_labels) <- c("res_F", "res_M", "res_20", "res_50", "res_interaction")

# labels for DEG by chromosome plots
effect_titles <- c("Female development", 
                   "Male development", 
                   "20 DPH", 
                   "50 DPH", 
                   "Interaction Effect of Sex and Age")
res_names <- c("res_F", "res_M","res_20", "res_50", "res_interaction")
sig_names <- c("sig_F", "sig_M", "sig_20", "sig_50", "sig_interaction")
table_of_names <- data.frame(res_name=res_names,
                             sig_name=sig_names,
                             effect_title=effect_titles)

# Save an .rds file
save_rds <- function(vbl){
    vbl_char <- deparse(substitute(vbl))
    filename <- file.path(working_dir, "rds_objects", paste0(vbl_char,".rds"))
    saveRDS(vbl, file = filename)
}

# Split string by character, useful for gene lists from GO results
split_string <- function(string, sep="/"){
    return(unlist(strsplit(string, sep)))
}
# Save an object to a file
save_rds <- function(vbl){
    vbl_char <- deparse(substitute(vbl))
    filename <- file.path(working_dir, paste0(vbl_char,".rds"))
    saveRDS(vbl, file = filename)
}

# Subset a dataframe or DESeqResults object by gene name
gene_subset <- function(df, genelist){
    if (!"gene" %in% colnames(df)){
        df <- rownames_to_column(data.frame(df), var="gene")
    }
    data <- df[df$gene %in% genelist, ]
    return(data)
}

### Calculate results for 4 contrasts and interaction effect
#
# Note that DESeq2::results() performs independent filtering by default 
# using the mean of normalized counts as a filter statistic. Genes
# that do not pass the filter threshold will have 'NA' as their padj.
#
# The default false discovery rate correction used here is the 
# Benjamini-Hochberg adjustment.
#
# Reference levels are sex=F, age=20
#
# Calculate results of DESeqDataSet for paired contrasts and sex+age interaction
get_results <- function(padj_thresh=0.01, lfc_thresh=NULL) {
    print("Effect of AGE in FEMALES")
    if (is.numeric(lfc_thresh)) {
        res_F <- results(ddsi, name = "age_50_vs_20", alpha=padj_thresh, lfcThreshold=lfc_thresh)
    } else {
        res_F <- results(ddsi, name = "age_50_vs_20", alpha=padj_thresh)
    }
    summary(res_F)
    sig_F <- filter_results(res_F, padj_thresh)
    
    print("Effect of SEX at 20 DPH")
    if (is.numeric(lfc_thresh)) {
        res_20 <- results(ddsi, name ="sex_M_vs_F", alpha=padj_thresh, lfcThreshold=lfc_thresh)
    } else {
        res_20 <- results(ddsi, name ="sex_M_vs_F", alpha=padj_thresh)
    }
    summary(res_20)
    sig_20 <- filter_results(res_20, padj_thresh)
    
    print("Effect of AGE in MALES")
    if (is.numeric(lfc_thresh)){
        res_M <- results(ddsi, list( c("age_50_vs_20", "sexM.age50") ), alpha=padj_thresh, lfcThreshold=lfc_thresh)
    } else {
        res_M <- results(ddsi, list( c("age_50_vs_20", "sexM.age50") ), alpha=padj_thresh)
    }
    summary(res_M)
    sig_M <- filter_results(res_M, padj_thresh)
    
    print("Effect of SEX at 50 DPH")
    if (is.numeric(lfc_thresh)) {
        res_50 <- results(ddsi, list( c("sex_M_vs_F", "sexM.age50")), alpha=padj_thresh, lfcThreshold=lfc_thresh)
    } else {
        res_50 <- results(ddsi, list( c("sex_M_vs_F", "sexM.age50")), alpha=padj_thresh )
    }
    summary(res_50)
    sig_50 <- filter_results(res_50, padj_thresh)
    
    print("Interaction effect of SEX and AGE")
    if (is.numeric(lfc_thresh)) {
        res_interaction <- results(ddsi, name="sexM.age50", alpha = padj_thresh, lfcThreshold=lfc_thresh)
    } else {
        res_interaction <- results(ddsi, name="sexM.age50", alpha = padj_thresh)
    }
    summary(res_interaction)
    sig_interaction <- filter_results(res_interaction, padj_thresh)
    
    res <- list("res_F"=res_F, "res_M"=res_M, "res_20"=res_20, "res_50"=res_50, "res_interaction"=res_interaction)
    sig <- list("sig_F"=sig_F, "sig_M"=sig_M, "sig_20"=sig_20, "sig_50"=sig_50, "sig_interaction"=sig_interaction)
    return(list("res_desobj"=res, "sig_df"=sig))
}

# Subset the results to return genes with padj < padj_thresh (drops NA's) & add gene info
# Options to sort by a specificed column, or filter for chromosome
filter_results <- function(results, padj_thresh=0.01, sort_by="padj", chrom=NULL){
    results %>%
        data.frame() %>%
        rownames_to_column(var="gene") %>% 
        as_tibble() %>% 
        dplyr::filter(padj < padj_thresh) %>%
        plyr::join(gene_info, by="gene", type="left", match="first") %>%
        {if (is.character(chrom)) filter(., chromosome == chrom) else .} %>%
        arrange((get(sort_by))) 
}

# Return df of log2FoldChange by gene
# (Helper function for lfc_summary())
lfc_df <- function(sig_df, suffix="_X"){
    cols <- c("gene", "l2fc")
    cols[2] <- paste0(cols[2], suffix)
    sig_df %>%
        arrange(dplyr::desc(log2FoldChange)) %>%
        dplyr::select(gene, log2FoldChange) %>%
        `colnames<-`(cols)
}

# Combine l2FC values from all contrasts into a single dataframe
# Add info on ESTIMA clones and RA marker status
lfc_summary <- function(deg, output_filename=NULL){
    sigs <- deg$sig_df
    lfcM <- lfc_df(sigs$sig_F, "_F")
    lfcF <- lfc_df(sigs$sig_M, "_M")
    lfc20 <- lfc_df(sigs$sig_20, "_20")
    lfc50 <- lfc_df(sigs$sig_50, "_50")
    lfcint <- lfc_df(sigs$sig_interaction, "_int")
    lfc_all <- Reduce(function(x, y) merge(x, y, by="gene", all=TRUE), list(lfcM, lfcF, lfc20, lfc50, lfcint))
    lfc_all <- plyr::join(lfc_all, gene_info, by="gene")
    lfc_all$description <- sapply(lfc_all$description, function(x) paste(unlist(x), collapse=", "))
    dups <- lfc_all %>% dplyr::filter(gene %in% subset(lfc_all, duplicated(gene))$gene)  # get duplicate indices
    lfc_all <- lfc_all[!duplicated(lfc_all$gene), ]  # remove duplicates 
    row.names(lfc_all) <- lfc_all$gene  # set rownames to genes
    if (!is.null(output_filename)){
        write.xlsx(lfc_all, output_filename, row.names = FALSE, showNA = FALSE)
    }
    return(lfc_all)
}

### Subset a dataframe by a term-based search of the gene description column
term_subset <- function(search_string, df=lfc){
    if ("description" %in% names(df)){
        return(df[grep(search_string, df$description), ])
    }
    else {
        print("Provided dataframe is missing a 'description' column.")
    }
}

# Calculate the percentage of DEGs per chromosome
deg_by_chromosome <- function(sig_df){
    total_genes_per_chromosome <- plyr::count(gene_info, "chromosome")
    deg_counts <- plyr::count(sig_df, "chromosome")
    chr_counts <- plyr::join(deg_counts, total_genes_per_chromosome, by="chromosome", type="right")
    colnames(chr_counts) <- c("chromosome", "n_deg", "n_chr")
    chr_counts$n_deg[is.na(chr_counts$n_deg)] <- 0
    chr_counts$chromosome <- factor(chr_counts$chromosome, levels = chr_counts$chromosome)
    chr_counts$pct_chr_deg <- (chr_counts$n_deg / chr_counts$n_chr) * 100
    chr_counts <- chr_counts[complete.cases(chr_counts),]  # remove NAs
    chr_counts <- chr_counts[order(match(chr_counts$chromosome, chr_order)),]  # custom sort
    chr_counts$chromosome <- factor(chr_counts$chromosome, levels = chr_order)
    return(chr_counts)
}

# Make tidy dataframe of rlog (log2) normalized expression values
tidy_rld <- function(deseq_transform_obj){
    rld_mat <- assay(deseq_transform_obj)
    rld_df <- rownames_to_column(data.frame(rld_mat), "gene")
    rld_tidy <- gather(data.frame(rld_df), id, log2, -gene)
    rld_tidy$id <- sub("^X", "", rld_tidy$id)
    rld_tidy$group <- coldata[rld_tidy$id, "group"]
    return(rld_tidy)
}

########## PLOTTING FUNCTIONS ##########

# Plot % of genes on each chromosome that are DEGs
plot_deg_by_chr <- function(plotdata){
    p <- ggplot(plotdata[1:33,], aes(x=chromosome, 
                                     y=pct_chr_deg,
                                     fill=color)) +
        geom_col() +
        scale_fill_manual(values=c("#000000", "#FF0000")) +
        theme_bw() +
        xlab("Chromosome") +
        ylab("") +
        ggtitle("DEGs by chromosome") +
        theme(legend.position="none")
    return (p)
}

# Custom MA plotting function based on geneplotter::plotMA
# Uses data return option from DESeq2's plotMA() function
plot_MA <- function(shrunken_res, plot_title='', y_lim=c(-4,4), alpha=ALPHA, 
                    tri_size=1){
    ma_data <- plotMA(shrunken_res,
                      alpha = alpha,
                      returnData=TRUE)
    ma_data <- subset(ma_data, mean != 0)
    ma_data["lfc_bounded"] <- pmax(y_lim[1], pmin(y_lim[2], ma_data$lfc))
    ggplot(data=ma_data, aes(x=mean, y=lfc_bounded)) +
        geom_point(size=ifelse(ma_data$lfc<y_lim[1], tri_size, ifelse(ma_data$lfc>y_lim[2], tri_size, 0.5)),
                   col=ifelse(ma_data$isDE, "mediumblue", "gray65"), 
                   pch=ifelse(ma_data$lfc<y_lim[1], 6, ifelse(ma_data$lfc>y_lim[2], 2, 16))) +
        scale_x_log10() +
        coord_cartesian(ylim=y_lim) +
        ylab(expression('log'[2]*' fold change')) +
        xlab("Mean of normalized counts") +
        geom_hline(yintercept=0, color="black") +
        ggtitle(plot_title) +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              axis.title=element_text(size=12),
              axis.text=element_text(size=10),
              plot.title=element_text(size=12, hjust=0.5)
        )
}

# Plot log2 abundance of gene(s)
plot_abundance <- function(genelist, log2_tidy=rld_tidy, num_col=5, no_legend=TRUE){
    group_labels <- c("Female 20 DPH", "Female 50 DPH", "Male 20 DPH", "Male 50 DPH")
    data <- gene_subset(log2_tidy, genelist)
    data$gene <- factor(data$gene, levels = genelist)  # force ggplot to match plot order to genelist
    p <- ggplot(data, aes(x=group, y=log2)) +
        geom_point(aes(color=group), position=position_jitter(w = 0.1,h = 0)) +  
        labs(color="") +
        facet_wrap(~ gene, ncol=num_col, scales="free_y") + 
        xlab("") +
        ylab(expression('log'[2]*' abundance')) +
        scale_x_discrete(labels=group_labels) +
        scale_color_manual(values = paired_palette, labels=group_labels) +
        theme_bw() +
        theme(strip.background =element_rect(fill="white")) +
        theme(axis.text.x = element_text(angle=45, vjust=1, hjust=1))
    if (no_legend){
        p <- p + theme(legend.position="none")
    }
    return(p)
}

# Plot histogram of log2FC M:F for all autosomal versus Z genes
plot_AvsZ <- function(result, 
                      title = element_blank(), 
                      legend = TRUE,
                      leg_pos = c(0.8, 0.9),
                      xlimits = c(-5,5),
                      ylimits=c(0,60))
  {
  plotdata <- data.frame(result)
  plotdata <- rownames_to_column(plotdata, var = "gene")
  plotdata <- plyr::join(plotdata, gene_info, by="gene")
  plotdata <- plotdata %>% dplyr::select(gene, log2FoldChange, chromosome)
  index <- plotdata$chromosome != "Z"
  plotdata$chromosome[index] <- "A"
  plotdata <- drop_na(plotdata)
  print(tapply(plotdata$log2FoldChange, plotdata$chromosome, summary)) #summary stats
  p <- ggplot(plotdata, aes(x=log2FoldChange, fill=chromosome)) +
    geom_histogram(aes(y=c(..count..[..group..==1]/sum(..count..[..group..==1]),
                           ..count..[..group..==2]/sum(..count..[..group..==2]))*100),
                   position='dodge', binwidth=0.5) +
    scale_fill_manual(#values = c("#f27d29", "#2b7bc2"),
                      values = c("#000000", "#FB0505"),
                      labels = c("Autosome", "Z")) +
    theme_bw() +
    ylab("Percentage of genes") +
    ggtitle(title) +
    coord_cartesian(xlim = xlimits, ylim=ylimits)
  if (legend == FALSE){
   p <- p + theme(legend.title = element_blank(),
          legend.position = "none")
  }
  else {
    p <- p + theme(legend.text=element_text(size=8),
                   legend.title = element_blank(),
                   legend.margin=margin(c(2,3,2,2)),
                   legend.position = leg_pos,
                   legend.background = element_blank(),
                   legend.spacing.y = unit(0, "mm"),
                   legend.spacing.x = unit(0.5, 'mm'),
                   legend.box.background = element_rect(colour = "black")
                   )
  }
  return(p)
}

# Custom clustering plots from degPatterns() results
plot_clusters <- function(data, num_rows=1, leg_pos="bottom"){
  new_titles <- gsub("Group: ", "Cluster", as.character(data$title))
  new_titles <- gsub("- genes: ", "(n=", new_titles)
  new_titles <- paste0(new_titles, ")")
  data$title <- new_titles
  p <- ggplot(data,
              aes(age, value, color = sex, fill = sex)) +
    geom_boxplot(alpha = 0,
                 outlier.size = 0,
                 outlier.shape = 0) +
    geom_point(aes(shape=sex),
               position = position_jitterdodge(dodge.width = 0.9),
               alpha = 0.4, 
               size = 1) +
    # change the method to make it smoother
    geom_smooth(aes(group=sex), 
                method = "lm",
                se = FALSE) +
    scale_x_discrete(labels=c("20 DPH", "50 DPH")) +
    guides(fill="none") + # clear legend to then apply custom legend below
    scale_color_manual(values = c("#f27d29", "#2b7bc2"),
                       labels = c("Female", "Male")) +
    scale_shape_manual(values = c(19, 17),
                       labels = c("Female", "Male")) +
    guides(color = guide_legend(override.aes = list(linetype = 0, 
                                                    size=4, 
                                                    alpha = 0.7))) +
    facet_wrap(~title,
               nrow = num_rows) +
    ylab("Z-score of gene abundance") +
    xlab(NULL) +
    theme_bw() +
    theme(strip.background=element_rect(fill="white"),
          legend.title = element_blank(),
          legend.text=element_text(size=8),
          legend.margin=margin(c(2,3,2,2)),
          legend.position = leg_pos,
          legend.background = element_blank(),
          legend.spacing.y = unit(0, "mm"),
          legend.spacing.x = unit(0, 'mm'),
          legend.box.background = element_rect(colour = "black")
    )
  return(p)
}

########## GO ANALYSIS FUNCTIONS ########## 

# Map zebra finch genes to NCBI human orthologs
# Returns df with additional columns containing human ortholog data:
# ensgene, entrez, chr, start, end, strand, biotype, description
map_to_human_orthologs <- function(res){
    if (class(res) == "character"){
        res_tb <- as_tibble(data.frame(gene=res))
    } else if (!"gene" %in% colnames(res)){
        res_tb <- res %>%
            data.frame() %>%
            rownames_to_column(var="gene") %>% 
            as_tibble()
    } else {
        res_tb <- as_tibble(res)
    }
    idx <- grch38$symbol %in% res_tb$gene
    ids <- grch38[idx, ]
    non_duplicates <- which(duplicated(ids$symbol) == FALSE)
    ids <- ids[non_duplicates, ] 
    res_orthos <- inner_join(res_tb, ids, by=c("gene"="symbol"))
    return(res_orthos)
}

# Overrepresentation analysis
go_enrich <- function(res_orthos=NULL, 
                      padj_thresh = 0.01, 
                      genelist=NULL, 
                      ontology_type="ALL",
                      q_value_cutoff=0.05){
    universe_genes <- map_to_human_orthologs(row.names(cts))$ensgene
    # define significant genes
    if (is.character(genelist)){
        sig_genes <- map_to_human_orthologs(genelist)$ensgene
    } else {
        sig_genes <- dplyr::filter(res_orthos, padj < padj_thresh)
        sig_genes <- as.character(sig_genes$ensgene)
    }
    ego <- enrichGO(gene = sig_genes, 
                    universe = universe_genes,
                    keyType = "ENSEMBL",
                    OrgDb = org.Hs.eg.db, 
                    ont = ontology_type, 
                    pAdjustMethod = "BH", 
                    qvalueCutoff = q_value_cutoff, 
                    readable = TRUE)
    return(ego)
}

# GSEA analysis
gsea_go <- function(res_orthos, pval_cutoff=0.05, ontology="BP"){
    res_entrez <- dplyr::filter(res_orthos, entrez != "NA")
    res_entrez <- res_entrez[which(duplicated(res_entrez$entrez) == F), ]
    foldchanges <- res_entrez$log2FoldChange
    names(foldchanges) <- res_entrez$entrez
    foldchanges <- sort(foldchanges, decreasing=TRUE)
    gseaGO <- gseGO(geneList = foldchanges,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = ontology,
                     minGSSize = 10,
                     pvalueCutoff = pval_cutoff,
                     pAdjustMethod = "BH", 
                     verbose = FALSE)
    return(gseaGO)
}

# Convert EntrezIDs to gene symmbols
core_enrich_to_symbol <- function(gse_result){
    for (i in seq_along(gse_result$core_enrichment)){
        entrez_list <- unlist(strsplit(gse_result$core_enrichment[[i]], "/"))
        idx <- grch38$entrez %in% entrez_list
        genes <- unique(grch38$symbol[idx])
        genes_string <- paste(genes, collapse="/")
        gse_result$core_enrichment[i] <- genes_string
    }
    return(gse_result)
}

# Perform ORA on each contrast DEG set for a given ontology category
ora_contrasts <- function(ontology_cat){
    ora_results <- list()
    for(res_name in names(deg$res_desobj)){
        res <- deg$res_desobj[[res_name]]
        res_orthos <- map_to_human_orthologs(res)
        go_result <- go_enrich(res_orthos, 
                               padj_thresh = ALPHA,
                               q_value_cutoff = 0.05,
                               ontology_type = ontology_cat
        )
        go_result <- data.frame(simplify(go_result))
        ora_results[[res_name]] <- go_result
    }
    return(ora_results)
}

# Perform GSEA on each contrast DEG set for a given ontology category
gsea_contrasts <- function(ontology_cat){
    gsea_results <- list()
    for (res_name in names(deg$res_desobj)){
        print(paste("...Testing GO term Gene Set Enrichment for", result_labels[res_name], "DEGs"))
        res <- deg$res_desobj[[res_name]]
        res_orthos <- map_to_human_orthologs(res)
        gseaGO <- gsea_go(res_orthos, 
                          pval_cutoff = 0.05, 
                          ontology=ontology_cat)
        gsea <- simplify(gseaGO)@result
        if (dim(gsea)[1]!=0){
            gsea <- core_enrich_to_symbol(gsea)
        }
        gsea_results[[res_name]] <- gsea
    }
    return(gsea_results)
}

# Perform ORA on each cluster DEG set for a given ontology category
ora_clusters <- function(cluster_df, ontology_cat){
    ora_result <- list()
    for (clust_num in unique(cluster_df$cluster)){
        print(paste("...Testing GO term over-representation for Cluster", clust_num))
        cluster <- cluster_df %>% dplyr::filter(cluster == clust_num)
        go_result <- go_enrich(genelist = cluster$gene,
                               padj_thresh = ALPHA,
                               q_value_cutoff = 0.05,
                               ontology_type = ontology_cat)
        go_result <- data.frame(simplify(go_result))
        ora_result[[paste("Cluster",clust_num)]] = go_result
    }
    return(ora_result)
}

# Combine two sets (e.g. two ontology categories) of GO results
combine_go_results <- function(x, y){
    ora <- list()
    for (res_name in names(x)) {
        df1 <- data.frame(x[[res_name]])
        df2 <- data.frame(y[[res_name]])
        combined <- rbind(df1, df2)
        sorted <- combined[order(combined$qvalue),]
        ora[[res_name]] <- sorted
    }
    return( ora)
}

