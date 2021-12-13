# deseq_analysis.R
#
# Differential gene analysis of zebra finch RA using DESeq2 for the publication titled:
# Emergence of sex-specific transriptomes in a sexually dimorphic brain nucleus
# 
# BY Sami Friedrich
# CREATED 06/01/2021
# UPDATED 12/13/2021

options(java.parameters = "-Xmx8000m")
working_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(working_dir)
input_path <- "input_files"
figures_path <- "output_figures"
dir.create(figures_path)
tables_path <- "output_tables"
dir.create(tables_path)

##################################################
# Load data and build DESeq2 model
##################################################

### Build count matrix and feature tables
source("make_count_and_feature_tables.R")
### Load libraries and custom functions 
source("de_functions.R")
### Create DESeq2 object
ddsi <- DESeqDataSetFromMatrix(countData = cts,
                               colData = coldata,
                               design = ~ sex + age + sex:age)
ddsi <- DESeq(ddsi)
### Make log2 abundance tables for plotting
rldi <- rlog(ddsi, blind=TRUE)  # transform count data to log2 scale and normalize to library size
rld_mat <- assay(rldi)  # matrix of log2 abundances after rlog transformation
rld_tidy <- tidy_rld(rldi)  # tidy df in long format of rlog transformed log2 abundance

# Check distribution of gene dispersions
plotDispEsts(ddsi)

##################################################
# Sample clustering using PCA
##################################################
# PCA based on 500 most variant genes
plotPCA(rldi, intgroup="group", ntop=500)
data <- plotPCA(rldi, intgroup="group", ntop=500, returnData=TRUE)
p <- ggplot(data, aes(x=PC1, y=PC2)) +
  geom_point(aes(fill=group), colour="white",pch=21, size=4) +
  scale_fill_manual(values=paired_palette, labels=group_labels) +
  theme_bw() +
  xlab("PC1: 35% variance") +
  ylab("PC2: 13% variance") +
  labs(fill="Group") +
  theme(strip.background =element_rect(fill="white"),
        
        axis.title=element_text(size=12),
        axis.text=element_text(size=10),
        plot.title=element_text(size=12, hjust=0.5)
  )
p
ggsave(file.path(figures_path, "PCA_n500.pdf"),
       plot = p,
       device = "pdf",
       width = 5.5, height = 4, units = "in")

# PCA based on all genes
plotPCA(rldi, intgroup="group", ntop=nrow(rld_mat))  # plot to get % variance for PC1 and PC2
data <- plotPCA(rldi, intgroup="group", ntop=nrow(rld_mat), returnData=TRUE)
p <- ggplot(data, aes(x=PC1, y=PC2)) +
  geom_point(aes(fill=group), colour="white",pch=21, size=4) +
  scale_fill_manual(values=paired_palette, labels=group_labels) +
  theme_bw() +
  xlab("PC1: 26% variance") +
  ylab("PC2: 14% variance") +
  labs(fill="Group") +
  theme(strip.background =element_rect(fill="white"),
        
        axis.title=element_text(size=12),
        axis.text=element_text(size=10),
        plot.title=element_text(size=12, hjust=0.5)
  )
p
ggsave(file.path(figures_path, "PCA_all-genes.pdf"),
       plot = p,
       device = "pdf",
       width = 5.5, height = 4, units = "in")

##################################################
# Differential gene expression contrasts
##################################################
ALPHA = 0.01
LFC_THRESH = NULL
# Test for differential expression to specified alpha level, option to also test to lfc threshold
deg <- get_results(padj_thresh=ALPHA, lfc_thresh = LFC_THRESH)
# Save differential gene expression results to .xlsx
for (res_name in names(deg$res_desobj)){
  result <- data.frame(deg$res_desobj[[res_name]])
  result <- rownames_to_column(result, var = "gene")
  result <- plyr::join(result, gene_info, by="gene")
  result <- result[order(result$padj),]
  result$description <- lapply(result$description, 
                               FUN = function(x) paste(unlist(x), collapse=" | "))
  result$GeneID <- as.character(result$GeneID)
  write.xlsx2(result,
              file=file.path(tables_path, "Differential_expression_stats.xlsx"),
              append=TRUE,
              sheetName=result_labels[res_name],
              row.names=FALSE,
              col.names=TRUE,
              showNA=FALSE
  )
}
# Make a df of l2fc values for all genes with at least one significant contrast
lfc <- lfc_summary(deg)

##################################################
# RNA-Seq validation with RA ZEBrA markers
##################################################
zebra <- read_csv(file.path(input_path, "ZEBrA_MarkerTable_070721.csv"))
zebra_ra <- zebra %>% dplyr::filter(StructureID == "RA")  # filter for RA rows
zebra_ra$GeneID[zebra_ra$GeneID == "MST4"] <- "STK26"  # rename MST4 to STK26
zebra_ra <- dplyr::arrange(zebra_ra, GeneID)
marker_genes_in_dataset <- zebra_ra[zebra_ra$GeneID %in% rownames(deg$res_desobj$res_F),]
marker_genes_in_dataset <- marker_genes_in_dataset[,c(1,3)]
names(marker_genes_in_dataset) <- c("gene", "marker_type")
marker_genes_table <- plyr::join(marker_genes_in_dataset, lfc, by="gene")
marker_genes_table <- marker_genes_table[, c(1:2,4)]
marker_genes_table <- drop_na(marker_genes_table)
marker_genes_table$consistent <- rep(0, nrow(marker_genes_table))
for (idx in 1:nrow(marker_genes_table)){
  row <- marker_genes_table[idx,]
  if (row$marker_type == "+" & row$l2fc_M>0){
    marker_genes_table$consistent[idx] = 1
  } else if (row$marker_type == "-" & row$l2fc_M<0){
    marker_genes_table$consistent[idx] = 1
  }
}
write.xlsx2(marker_genes_table,
            file=file.path(tables_path, "ZEBrA_marker_validation.xlsx"),
            row.names=FALSE,
            col.names=TRUE,
            showNA=FALSE
)
# number positive markers consistent with increase during male development
sum(marker_genes_table$marker_type == "+" & marker_genes_table$l2fc_M>0)
# number negative markers consistent with decrease during male development
sum(marker_genes_table$marker_type == "-" & marker_genes_table$l2fc_M<0)

##################################################
# MA plots
##################################################
### Shrink LFC values for better data viz (removes noise from low count genes)
shrunk_20 <- lfcShrink(ddsi, contrast=c("group", "F_20", "M_20"), res=deg$res_desobj$res_20, type="ashr")
shrunk_50 <- lfcShrink(ddsi, contrast=c("group", "F_50", "M_50"), res=deg$res_desobj$res_50, type="ashr")
shrunk_M <- lfcShrink(ddsi, contrast=c("group", "M_20", "M_50"), res=deg$res_desobj$res_M, type="ashr")
shrunk_F <- lfcShrink(ddsi, contrast=c("group", "F_20", "F_50"), res=deg$res_desobj$res_F, type="ashr")

# Create combined figure of MA plots for all 4 contrasts with shared axes
p1 <- plot_MA(shrunk_20, "Male vs Female at 20 DPH")
p2 <- plot_MA(shrunk_50, "Male vs Female at 50 DPH")
p3 <- plot_MA(shrunk_F, "20 vs 50 DPH in Females")
p4 <- plot_MA(shrunk_M, "20 vs 50 DPH in Males")
# Remove axis titles from all plots
p = list(p1,p2,p3,p4) %>% map(~.x + labs(x=NULL, y=NULL))
ylab <- textGrob(expression('log'[2]*' fold change'),
                  rot = 90, gp = gpar(fontsize=12))
xlab <- textGrob("Mean of normalized counts",
                 gp = gpar(fontsize=12))
maplots <- grid.arrange(grobs=p, left = ylab, bottom = xlab)
ggsave(file.path(figures_path, "MAplots.pdf"),
       plot = maplots,
       device = "pdf",
       width = 8, height = 7, units = "in")

##################################################
# DEG summary stats and subsets
##################################################
### DEG LFC summary stats
mean(abs(lfc$l2fc_F), na.rm=TRUE)  # mean lfc in female development
median(abs(lfc$l2fc_F), na.rm=TRUE)  # median lfc in female development
mean(abs(lfc$l2fc_M), na.rm=TRUE)  # mean lfc in male development
median(abs(lfc$l2fc_M), na.rm=TRUE)  # median lfc in male development

mean(abs(lfc$l2fc_20), na.rm=TRUE)  # mean lfc for 20 DPH sex contrast
median(abs(lfc$l2fc_20), na.rm=TRUE)  # median lfc for 20 DPH sex contrast
mean(abs(lfc$l2fc_50), na.rm=TRUE)  # mean lfc for 50 DPH sex contrast
median(abs(lfc$l2fc_50), na.rm=TRUE)  # median lfc for 50 DPH sex contrast

shared_M_F_genes <- intersect(deg$sig_df$sig_M$gene, deg$sig_df$sig_F$gene)  # genes that change over age in M and F
shared_M_F <- lfc[shared_M_F_genes,]
shared_pos_F <- shared_M_F$l2fc_F > 0  # developmentally regulated genes in female with positive LFC 
shared_pos_M <- shared_M_F$l2fc_M > 0  # developmentally regulated genes in male with positive LFC

shared_M_F_same_direction <- shared_M_F[shared_pos_F == shared_pos_M, ]
shared_MF_same_direction_interaction <- shared_M_F_same_direction %>% dplyr::filter(!is.na(l2fc_int))
shared_MF_same_direction_no_interaction <- shared_M_F_same_direction %>% dplyr::filter(is.na(l2fc_int))
shared_M_F_opposite_direction <- shared_M_F[shared_pos_F != shared_pos_M, ]

M_only <- setdiff(deg$sig_df$sig_M$gene, deg$sig_df$sig_F$gene)  # genes that change over age only in M
F_only <- setdiff(deg$sig_df$sig_F$gene, deg$sig_df$sig_M$gene)  # genes that change over age only in F

##################################################
# Plotting DEGs by chromosome
##################################################
### Plot % DEGs by chromosome
plots <- list()
for (name in names(deg$sig_df)[3:4]){
  data <- deg$sig_df[[name]]
  p <- plot_deg_by_chr(deg_by_chromosome(data))
  p <- p + ggtitle(table_of_names[table_of_names$sig_name == name, "effect_title"])
  p <- p + coord_cartesian(ylim = c(0, 50))  # set common y axis
  plots[[name]] <- p
}
# Remove axis titles from all plots
p <- plots %>% map(~.x + labs(x=NULL, y=NULL))
# Plot together with shared axes labels
deg_by_chr <- grid.arrange(grobs=p, ncol = 1, nrow = 2, 
                    left = "% chromosome genes that are DEGs", 
                    bottom = "Chromosome")
ggsave(file.path(figures_path, "deg_by_chromosome.pdf"),
       plot = deg_by_chr,
       device = "pdf",
       width = 6.5, height = 4, units = "in")

###################################################################
# GO TERM ENRICHMENT TESTS FOR EACH CONTRAST DEG SET
###################################################################

### ORA for DEG sets
ora_bp_contrasts <- ora_contrasts(ontology_cat = "BP")
ora_mf_contrasts <- ora_contrasts(ontology_cat = "MF")
ora_results <- combine_go_results(ora_bp_contrasts, ora_mf_contrasts)
# Save as .xlsx
for (res_name in names(ora_results)){
  write.xlsx2(ora_results[[res_name]],
              file=file.path(tables_path, "ORA_contrasts.xlsx"),
              append=TRUE,
              sheetName=result_labels[res_name],
              row.names=FALSE,
              col.names=TRUE,
              showNA=FALSE
    )
}
### GSEA for DEG sets
gsea_mf_results <- gsea_contrasts(ontology_cat = "MF")
gsea_bp_results <- gsea_contrasts(ontology_cat = "BP")
gsea_results <- combine_go_results(gsea_bp_results, gsea_mf_results)
# Save as .xlsx
for (res_name in names(gsea_results)){
  write.xlsx2(gsea_results[[res_name]],
              file=file.path(tables_path, "GSEA_contrasts.xlsx"),
              append=TRUE,
              sheetName=result_labels[res_name],
              row.names=FALSE,
              col.names=TRUE,
              showNA=FALSE
  )
}

##################################################
# Cluster analysis of sex+age interaction genes
##################################################
### Hierarchical clustering
res <- deg$sig_df$sig_interaction
cluster_rlog <- rld_mat[res$gene, ] # subset rlog matrix for significant genes of particular effect
cl <- degPatterns(cluster_rlog, metadata = colData(rldi), time="age", col="sex", minc=9)  # perform clustering 
cl$plot <- degPlotCluster(cl$normalized, "age", "sex", lines = FALSE) +  # customize plot
  scale_x_discrete(labels=c("20 DPH", "50 DPH")) +
  guides(fill="none") +
  scale_color_discrete(name="Sex", labels=c("Female", "Male")) +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))
cl$plot
ggsave(file.path(figures_path, "sex+age_clusters.pdf"),
       plot = cl$plot,
       device = "pdf",
       width = 6, height = 7.5, units = "in")
### ORA of enriched GO terms in clusters
clusters <- cl$df
colnames(clusters) <- c("gene", "cluster")
ora_mf_clusters <- ora_clusters(clusters, ontology_cat = "MF")
ora_bp_clusters <- ora_clusters(clusters, ontology_cat = "BP")
cluster_ora_results <- combine_go_results(ora_bp_clusters, ora_mf_clusters)
# Save as .xlsx
for (res_name in names(cluster_ora_results)){
  write.xlsx2(cluster_ora_results[[res_name]],
              file=file.path(tables_path, "ORA_sex+age_clusters.xlsx"),
              append=TRUE,
              sheetName=res_name,
              row.names=FALSE,
              col.names=TRUE,
              showNA=FALSE
  )
}

############ Z GENE ANALYSES ############
### ZW pairs
lfc_z <- dplyr::filter(lfc, chromosome=="Z")
wz_pairs <- read_csv(file.path(input_path,"wzpairs.csv"), col_names = FALSE)
wz_pairs <- wz_pairs$X1
wz_lfc <- lfc_z[lfc_z$gene %in% wz_pairs, ]

### Summary stats of Z genes
# 20 DPH sex contrast
sig_20 <- deg$sig_df$sig_20
sum(sig_20$chromosome == "Z", na.rm = T)  # number of Z genes in 20 DPH sex contrast
sum(sig_20$chromosome == "Z" & sig_20$log2FoldChange > 0, na.rm=T) # number of male biased Z genes in 20 DPH sex contrast
sum(sig_20$chromosome == "Z", na.rm = T)/nrow(sig_20)  # % of 20 DPH DEGs that are Z genes
# % of Z genes that are male-biased at 20 DPH
sum(sig_20$chromosome == "Z" & sig_20$log2FoldChange > 0, na.rm=T) / sum(sig_20$chromosome == "Z", na.rm = T)
# 50 DPH sex contrast
sig_50 <- deg$sig_df$sig_50
sum(sig_50$chromosome == "Z", na.rm = T)  # number of Z genes in 50 DPH sex contrast
sum(sig_50$chromosome == "Z" & sig_50$log2FoldChange > 0, na.rm=T) # number of male biased Z genes in 50 DPH sex contrast
sum(sig_50$chromosome == "Z", na.rm = T)/nrow(sig_50)  # % of 50 DPH DEGs that are Z genes
# % of Z genes that are male-biased at 50 DPH
sum(sig_50$chromosome == "Z" & sig_50$log2FoldChange > 0, na.rm=T) / sum(sig_50$chromosome == "Z", na.rm = T)

### Hierarchical clustering of all Z chr DEGs
lfc_z <- dplyr::filter(lfc, chromosome=="Z")
z_genes <- lfc_z$gene[lfc_z$gene %in% rownames(rld_mat)]
cluster_rlog <- rld_mat[z_genes, ]
cl_z <- degPatterns(cluster_rlog, metadata = colData(rldi), time="age", col="sex", minc=9)
cl_z$plot <- degPlotCluster(cl_z$normalized, "age", "sex", lines = FALSE) +
  scale_x_discrete(labels=c("20 DPH", "50 DPH")) +
  guides(fill="none") + # clear legend to then apply custom legend below
  scale_color_discrete(name="Sex", labels=c("Female", "Male")) +
  theme_bw() +
  theme(strip.background =element_rect(fill="white"))
cl_z$plot
ggsave(file.path(figures_path, "z_clusters.pdf"),
       plot = cl_z$plot,
       device = "pdf",
       width = 6, height = 6, units = "in")

### Perform over-representation analysis on each cluster of Z chromosome DEGs
clusters_z <- cl_z$df
colnames(clusters_z) <- c("gene", "cluster")
z_cluster_mf_ora <- ora_clusters(clusters_z, ontology_cat = "MF")
z_cluster_bp_ora <- ora_clusters(clusters_z, ontology_cat = "BP")
z_cluster_ora_results <- combine_go_results(z_cluster_bp_ora, z_cluster_mf_ora)
# Save as .xlsx
for (res_name in names(z_cluster_ora_results)){
  write.xlsx2(z_cluster_ora_results[[res_name]],
              file=file.path(tables_path, "ORA_z_clusters.xlsx"),
              append=TRUE,
              sheetName=res_name,
              row.names=FALSE,
              col.names=TRUE,
              showNA=FALSE
  )
}

##################################################
# Gene families of interest analysis
##################################################
### Transcription factors
# GO term
go_id <- "GO:0003700" # DNA-binding transcription factor activity
go_genes <- genesbygo[go_id]
go_genes <- as.numeric(unlist(go_genes))
go_genes_symbol <- grch38$symbol[grch38$entrez %in% go_genes]
txf_go0003700 <- lfc %>% dplyr::filter(gene %in% go_genes_symbol)
# Term  based
txnfactors_lfc <- term_subset("transcription factor")
# Combine
txnf_lfc <- unique(rbind(txnfactors_lfc, txf_go0003700))
txnf_lfc$GeneID <- as.character(txnf_lfc$GeneID)
# Save
write.xlsx(txnf_lfc[,1:9], 
           file=file.path(tables_path, "functions_of_interest.xlsx"),
           append=TRUE,
           sheetName="Transcription factors",
           row.names=FALSE,
           col.names=TRUE,
           showNA=FALSE
)
### Apoptosis and cell death
# Term-based
apop_lfc <- term_subset("apoptosis")
death_lfc <- term_subset("death")
# Combine
apop_death_lfc <- unique(rbind(apop_lfc, death_lfc))
apop_death_lfc$GeneID <- as.character(apop_death_lfc$GeneID)
# Save
write.xlsx2(apop_death_lfc[,1:9],
            file=file.path(tables_path, "Functions_of_interest.xlsx"),
            append=TRUE,
            sheetName="Apoptosis regulators",
            row.names=FALSE,
            col.names=TRUE,
            showNA=FALSE
)
### Growth factor and related genes
growth_factor_lfc <- term_subset("growth factor")
growth_lfc <- term_subset("growth")
neurotrophic_lfc <- term_subset("neurotrophic")
# Combine
cellgrowth <- unique(rbind(growth_factor_lfc, growth_lfc, neurotrophic_lfc))
cellgrowth$GeneID <- as.character(cellgrowth$GeneID)
# Save
write.xlsx2(cellgrowth[,1:9],
            file=file.path(tables_path, "functions_of_interest.xlsx"),
            append=TRUE,
            sheetName="Growth factors",
            row.names=FALSE,
            col.names=TRUE,
            showNA=FALSE
)
### Steroid enzyme genes
steroid_term_lfc <- term_subset("steroid")
estrogen_subset_lfc <- term_subset("estrogen")
androgen_subset_lfc <- term_subset("androgen")
# Combine
steroid_lfc <- unique(rbind(steroid_term_lfc, estrogen_subset_lfc, androgen_subset_lfc))
steroid_lfc$GeneID <- as.character(steroid_lfc$GeneID)
# Save
write.xlsx2(steroid_lfc[,1:9],
            file=file.path(tables_path, "functions_of_interest.xlsx"),
            append=TRUE,
            sheetName="Sex steroid related_new",
            row.names=FALSE,
            col.names=TRUE,
            showNA=FALSE
)
### Voltage-gated channel genes
vg_lfc <- term_subset("voltage-gated")
vg_lfc $GeneID <- as.character(vg_lfc$GeneID)
# Save
write.xlsx2(vg_lfc[,1:9],
            file=file.path(tables_path, "functions_of_interest.xlsx"),
            append=TRUE,
            sheetName="Voltage-gated channels",
            row.names=FALSE,
            col.names=TRUE,
            showNA=FALSE
)
