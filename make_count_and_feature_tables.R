# make_count_and_feature_tables.R
#
# Builds count matrix from RNA-seq reads, and supporting data objects
# for the publication:
# Emergence of sex-specific transriptomes in a sexually dimorphic brain nucleus
#
# BY Sami Friedrich
# CREATED 06/01/2021
# UPDATED 04/21/2022

library(plyr)
library(dplyr)
library(readr)
library(tibble)

working_dir <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(working_dir)
counts_dir <- file.path(working_dir, "input_files","STAR_counts")

# Custom function to save an .rds
save_rds <- function(vbl){
    vbl_char <- deparse(substitute(vbl))
    filename <- file.path(working_dir, "rds_objects", paste0(vbl_char,".rds"))
    saveRDS(vbl, file = filename)
}

############ BUILD A COUNT MATRIX ############
counts_files <- list.files(path = counts_dir, pattern = "*ReadsPerGene.out.tab$", full.names = TRUE)
counts_list <- lapply(counts_files, read.table, skip = 4)  # skip first 4 lines (not gene counts)
cts <- as.data.frame(sapply (counts_list, function(x) x[ , 2 ]) )  # extract counts from second row
counts_files <- gsub( "_ReadsPerGene[.]out[.]tab", "", counts_files )
basefile <- sapply(strsplit(counts_files, "/"), tail, 1)
sample_ids <- sapply(strsplit(basefile, "_"), getElement, 2)
colnames(cts) <- sample_ids  # assign sample IDs to columns
row.names(cts) <- counts_list[[1]]$V1  # assign gene symbols to rownames
rm(counts_list, basefile, sample_ids)
rows_to_remove <- c("transcript_id")
cts <- cts[!(row.names(cts) %in% rows_to_remove), ] # remove "transcript_id" row
### Manually annotate select LOC genes to their gene symbols
manual_annot <- read_delim(file.path("input_files", "manual_annotations.tsv"), "\t")
row.names(cts) <- plyr::mapvalues(row.names(cts), 
                               from = manual_annot$loc, 
                               to = manual_annot$gene)

############ BUILD COLDATA (SAMPLE INFO) ############
coldata <- read.csv(file.path("input_files", "sample_info.csv"), row.names = 2)
coldata <- coldata[,-1]
coldata$group <- paste0(coldata$sex, "_", coldata$age) 
# sanity check that all samples are represented in both cts and coldata
all(rownames(coldata) %in% colnames(cts))
all(colnames(cts) %in% rownames(coldata))
# check to see if samples are in same order between coldata rows and cts columns
all(rownames(coldata) == colnames(cts))
# if not, rearrange and confirm
cts <- cts[, rownames(coldata)]
all(rownames(coldata) == colnames(cts))
#save_rds(cts)  # save to file
coldata$lane <- sub(".*_L00", "", counts_files)  # add lane number column
coldata$sex <- factor(coldata$sex)
coldata$age <- factor(coldata$age)
coldata$group <- factor(coldata$group)
coldata$lane <- factor(coldata$lane)
#save_rds(coldata)

############ BUILD GENE INFO TABLE ############
feature_file <- file.path(working_dir, "input_files", "GCF_003957565.1_bTaeGut1_v1.p_feature_table.csv")
features <- read.csv(feature_file, sep=",", fill=TRUE, header=TRUE)
# rename columns
features_colnames <- colnames(features)
features_colnames[1] <- "feature_type"
features_colnames[15] <- "gene"
colnames(features) <- features_colnames
#save_rds(features)
### Create table of gene descriptions, collapsing across gene
gene_desc <- select(features, gene, name)
gene_desc <- aggregate(gene_desc$name, list(gene_desc$gene), paste, sep=",")
colnames(gene_desc) <- c("gene", "description")
gene_desc$description <- lapply(gene_desc$description, function(x) x[x != "" & x != "\n"])  # remove empty strings
gene_info <- features %>%
    select(gene, 
           GeneID, 
           chromosome) %>%
    distinct()  %>%
    join(gene_desc, by="gene", type="left") %>%
    arrange(gene)
rm(gene_desc)
gene_info <- gene_info[-1,]  # remove first row of tRNA names from MT chr
gene_info$chromosome[gene_info$chromosome == ""] <- "unplaced_scaffold"  # annotate unplaced
### Manually annotate select LOC genes to their gene symbols
gene_info$gene <- plyr::mapvalues(gene_info$gene, 
                                  from = manual_annot$loc, 
                                  to = manual_annot$gene)
#save_rds(gene_info)

### Write gene count matrix to file with Entrez Gene IDs for GEO submission
gene_count_matrix <- rownames_to_column(cts, var = "gene_symbol")
# use named list to map gene name to Entrez GeneID
gene_dict <- gene_info$GeneID
names(gene_dict) <- gene_info$gene
EntrezGeneID <- gene_dict[gene_count_matrix$gene]
gene_count_matrix <- cbind(EntrezGeneID, gene_count_matrix[,-1])
write.table(
    gene_count_matrix,
    file.path(tables_path, "raw_gene_counts_matrix.tsv"),
    col.names = TRUE,
    row.names = FALSE,
    sep = "\t",
    quote = FALSE
)
rm(gene_count_matrix, gene_dict, EntrezGeneID)

############ BUILD VECTOR SPECIFYING CHROMOSOME ORDER ############
chr_order <- c("1", "1A", "2", "3", "4", "4A", "5", "6", "7", "8", "9", "10", 
               "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", 
               "21", "22", "23", "24", "25", "26", "27", "28", "29", "30",
               "Z", "unplaced_scaffold", "MT")
#save_rds(chr_order)


