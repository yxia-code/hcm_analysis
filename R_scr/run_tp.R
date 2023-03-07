library(nlme)
library(tidyverse)
library(ggplot2)
library(compositions)
library(ggrepel)

args <- commandArgs(TRUE)

## import data
table_path <- as.character(args[[1]])
meta_path <- as.character(args[[2]])
main_var <- as.character(args[[3]])
p_cutoff <- as.numeric(args[[4]])
out_path <- as.character(args[[5]])
adj_formula <- as.character(args[[6]])

source("R_scr/ancom_v2.1.R")

otu_data = read_tsv(table_path)
otu_id = otu_data$`#OTU ID`
otu_data = data.frame(otu_data[, -1], check.names = FALSE)
rownames(otu_data) = otu_id

feature_table = otu_data

meta_data = read_tsv(meta_path)
meta_data = meta_data %>% rename(Sample.ID = SampleID)

### Run ANCOM
res = ANCOM(feature_table = feature_table, meta_data = meta_data, struc_zero = NULL, main_var = main_var, p_adj_method = 'BH', alpha=p_cutoff, adj_formula=adj_formula, rand_formula=NULL)

### output
write_csv(res$out, paste(out_path, 'ancom_results.csv', sep='/'))
pdf(paste(out_path, 'volcano_plot.pdf', sep='/'))

### Setting volcano plot
n_taxa = nrow(feature_table)
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data
dat_ann = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.6"], label = "W[0.6]")
dat_ann2 = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.7"], label = "W[0.7]")
dat_ann3 = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.8"], label = "W[0.8]")
dat_ann4 = data.frame(x = min(res$fig$data$x), y = cut_off["detected_0.9"], label = "W[0.9]")
fig = res$fig + geom_text_repel(data=res$fig$data, aes(x=x, y=y, label=taxa_id)) +  
  geom_hline(yintercept = cut_off["detected_0.6"], linetype = "dashed") +  
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) + 
  geom_hline(yintercept = cut_off["detected_0.7"], linetype = "dashed") +  
  geom_text(data = dat_ann2, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) +
  geom_hline(yintercept = cut_off["detected_0.8"], linetype = "dashed") +  
  geom_text(data = dat_ann3, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) +
  geom_hline(yintercept = cut_off["detected_0.9"], linetype = "dashed") +  
  geom_text(data = dat_ann4, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) 
fig  
dev.off()
