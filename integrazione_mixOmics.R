# Impostazione cartella di lavoro
if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
  cat("Cartella di lavoro impostata su:", getwd(), "\n")
} else {
  cat("Impostare la cartella di lavoro manualmente.\n")
}

library(readr)
library(readxl)
library(tidyr)
library(dplyr)
library(DESeq2)
library(pcaMethods) 
library(MOFA2) 
library(ggplot2)
library(sva) 
library(org.Hs.eg.db)
library(msigdbr)
library(mixOmics)
library(clusterProfiler)
library(ggrepel)
library(mixOmics)
library(reticulate)
library(ggVennDiagram)
library(plotly)
library(patchwork)
library(grid)
library(UpSetR)
library(stringr)
library(enrichplot)
library(writexl)
library(edgeR)   # serve per LOG2-CPM (pre-processing trascrittomica)
library(patchwork)
library(fgsea) # per analisi GSEA
library(openxlsx)
library(RColorBrewer) # per plot grafico log2-cpm
library(GGally)
library(WGCNA)
allowWGCNAThreads()
use_condaenv("/Users/ceciliabergamini/miniforge3/envs/mofa_env", required = TRUE)

# Risoluzione conflitti esplicita
first  <- dplyr::first
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename

# Caricamento dei file

## Trascrittomica ##
meta2 <- read_excel("/Users/ceciliabergamini/Desktop/TESI/COL3A1_Data_&_Scripts/COL3A1_sampleAttributes_3D.xlsx") 

trascr2D <- read_table("/Users/ceciliabergamini/Desktop/TESI/COL3A1_Data_&_Scripts/COL3A1_rawCounts_mod_2D.txt")
trascr3D <- read_table("/Users/ceciliabergamini/Desktop/TESI/COL3A1_Data_&_Scripts/COL3A1_rawCounts_mod_3D.tsv")

# Eliminazione delle righe con valori mancanti per le colonne 'ensembl_gene_id', 'gene_id'
trascr2D <- trascr2D %>% drop_na(ensembl_gene_id) 
trascr3D <- trascr3D %>% drop_na(gene_id)

neworder2D <- c(10,2,1,7,14,12,13,9,11,8,6,4,3,5)
neworder3D <- c(1,7,6,4,5,8,3,12,2,13,9,10,11,14)
neworderMeta <- c(1,6,5,3,4,7,2,9,8,10,11,12,13,14)

# Seleziono i pazienti rispetto al nuovo ordine e rendo la struttura data frame
matrix2D <- as.data.frame(trascr2D[, neworder2D+1])
matrix3D <- as.data.frame(trascr3D[, neworder3D+2])
meta <- meta2[neworderMeta, ]

# Inserisco come nomi delle righe i geni
rownames(matrix2D) <- trascr2D$ensembl_gene_id
rownames(matrix3D) <- trascr3D$gene_id

# Siccome legge i valori come double, li trasformo in integer
trascr2D_init <- matrix2D %>% mutate(across(everything(), as.integer))
trascr3D_init <- matrix3D %>% mutate(across(everything(), as.integer))

colnames(trascr2D_init) <- c("19-1597_2D", "20-3154_2D", "20-3156_2D", "20-4456_2D", "20-4222_2D", "22-0547_2D", 
                        "22-0548_2D", "GAU_2D", "SF_2D", "VL_2D", "C1A_2D", "C2B_2D", "C3B_2D", "RR_2D")
colnames(trascr3D_init) <- c("19_1597_3D", "20_3154_3D", "20_3156_3D", "20_4456_3D", "20_4222_3D", "22_0547_3D",
                        "22_0548_3D", "GAU_3D", "SF_3D", "VL_3D", "C1A_3D", "C2B_3D", "C3B_3D", "RR_3D" )

meta$Condition <- as.factor(meta$Condition)


#########################################################################
##             GRAFICO OVERDISPERSION DATI TRASCR 2D E 3D              ##
#########################################################################
par(mfrow = c(1, 2))

for (nm in c("2D", "3D")) {
  counts <- as.matrix(if (nm == "2D") trascr2D_init else trascr3D_init)
  counts <- sweep(counts, 2, estimateSizeFactorsForMatrix(counts), "/")
  
  m <- rowMeans(counts)
  v <- apply(counts, 1, var)
  ok <- m > 0 & v > 0
  
  cat("\nTrascrittomica", nm, "\n")
  cat("Geni analizzati:", sum(ok), "\n")
  cat("% geni con varianza > media:", round(100 * mean(v[ok] > m[ok]), 1), "\n")
  cat("Mediana varianza/media:", round(median(v[ok] / m[ok]), 2), "\n")
  
  plot(m[ok], v[ok], log = "xy", pch = 16, cex = 0.4, col = rgb(0, 0, 0, 0.2),
       xlab = "Media", ylab = "Varianza", main = paste("Trascrittomica", nm))
  abline(0, 1, col = "red", lwd = 2)
  legend("topleft", legend = "Poisson (var = media)", col = "red", lwd = 2, bty = "n")
}

#########################################################################
###                   PRE-PROCESSING TRASCRITTOMICA                   ###
#########################################################################
cat("Numero iniziale di geni 2D:", nrow(trascr2D_init), "\n") # 33,694

dge2D <- DGEList(counts = trascr2D_init, group = meta$Condition)
L2D <- mean(dge2D$samples$lib.size) * 1e-6
M2D <- median(dge2D$samples$lib.size) * 1e-6

lcpm_raw_2D <- cpm(dge2D, log = TRUE)
keep_2D <- filterByExpr(dge2D, group = meta$Condition)

dge2D_filt <- dge2D[keep_2D, , keep.lib.sizes = FALSE]
trascr2D <- trascr2D_init[keep_2D, ] 
lcpm_filt_2D <- cpm(dge2D_filt, log = TRUE)

cat("Numero di geni 2D dopo filtraggio:", nrow(trascr2D), "\n") # 10,973

# ==========================================
# 2. ELABORAZIONE MATRICE 3D
# ==========================================
cat("Numero iniziale di geni 3D:", nrow(trascr3D_init), "\n") # 33,694

dge3D <- DGEList(counts = trascr3D_init, group = meta$Condition)
L3D <- mean(dge3D$samples$lib.size) * 1e-6
M3D <- median(dge3D$samples$lib.size) * 1e-6

lcpm_raw_3D <- cpm(dge3D, log = TRUE)
keep_3D <- filterByExpr(dge3D, group = meta$Condition)

dge3D_filt <- dge3D[keep_3D, , keep.lib.sizes = FALSE]
trascr3D <- trascr3D_init[keep_3D, ]
lcpm_filt_3D <- cpm(dge3D_filt, log = TRUE)

cat("Numero di geni 3D dopo filtraggio:", nrow(trascr3D), "\n") # 11,376


# ==========================================
# 3. PLOT 2D (Prima vs Dopo)
# ==========================================
nsamples_2D <- ncol(dge2D)
samplenames_2D <- colnames(trascr2D)
col_2D <- colorRampPalette(brewer.pal(12, "Paired"))(nsamples_2D)
cutoff_2D <- log2(10 / M2D + 2 / L2D)

par(mfrow = c(1, 2)) 

# 2D - RAW
plot(density(lcpm_raw_2D[, 1]), col = col_2D[1], lwd = 2, ylim = c(0, 0.85), las = 2, main = "", xlab = "")
title(main = "A. Dati grezzi (2D)", xlab = "Log-cpm")
abline(v = cutoff_2D, lty = 3)
for (i in 2:nsamples_2D) {
  lines(density(lcpm_raw_2D[, i]), col = col_2D[i], lwd = 2)
}
legend("topright", samplenames_2D, text.col = col_2D, bty = "n", cex = 1.2)

# 2D - FILTERED
plot(density(lcpm_filt_2D[, 1]), col = col_2D[1], lwd = 2, ylim = c(0, 0.25), las = 2, main = "", xlab = "")
title(main = "B. Dati filtrati (2D)", xlab = "Log-cpm")
abline(v = cutoff_2D, lty = 3)
for (i in 2:nsamples_2D) {
  lines(density(lcpm_filt_2D[, i]), col = col_2D[i], lwd = 2)
}
legend("topright", samplenames_2D, text.col = col_2D, bty = "n", cex = 1.2)


# ==========================================
# 4. PLOT 3D (Prima vs Dopo)
# ==========================================
nsamples_3D <- ncol(dge3D)
samplenames_3D <- colnames(trascr3D)
col_3D <- colorRampPalette(brewer.pal(12, "Paired"))(nsamples_3D)
cutoff_3D <- log2(10 / M3D + 2 / L3D)

par(mfrow = c(1, 2)) # Finestra a 2 grafici per il 3D

# 3D - RAW
plot(density(lcpm_raw_3D[, 1]), col = col_3D[1], lwd = 2, ylim = c(0, 0.85), las = 2, main = "", xlab = "")
title(main = "A. Dati grezzi (3D)", xlab = "Log-cpm")
abline(v = cutoff_3D, lty = 3)
for (i in 2:nsamples_3D) {
  lines(density(lcpm_raw_3D[, i]), col = col_3D[i], lwd = 2)
}
legend("topright", samplenames_3D, text.col = col_3D, bty = "n", cex = 1.2)

# 3D - FILTERED
plot(density(lcpm_filt_3D[, 1]), col = col_3D[1], lwd = 2, ylim = c(0, 0.25), las = 2, main = "", xlab = "")
title(main = "B. Dati filtrati (3D)", xlab = "Log-cpm")
abline(v = cutoff_3D, lty = 3)
for (i in 2:nsamples_3D) {
  lines(density(lcpm_filt_3D[, i]), col = col_3D[i], lwd = 2)
}
legend("topright", samplenames_3D, text.col = col_3D, bty = "n", cex = 1.2)


# # ==============================================================================
# # 5. Filtraggio geni poco espressi (euristica personale)
# # ==============================================================================
# prob <- 2/3
# # --- Filtraggio Trascrittomica 2D ---
# k_2D <- 50
# index_alti <- rowSums(trascr2D_init) >= k_2D
# index <- rowSums(trascr2D_init) < k_2D     # Seleziono i geni con num di conteggi inferiore a 50
# bassi_malato <- trascr2D_init[index, 1:10]  # estrazione dei conteggi per i malati
# bassi_sano <- trascr2D_init[index, 11:14]   # estrazione dei conteggi per i sani
# logfoldchanges_bassi <- log2(apply(bassi_malato, 1, mean)/apply(bassi_sano, 1, mean)) # calcolo di logfoldchange
# inf_pos  <- is.infinite(logfoldchanges_bassi) & logfoldchanges_bassi > 0
# inf_neg  <- is.infinite(logfoldchanges_bassi) & logfoldchanges_bassi < 0
# soglia_malati <- ceiling(prob * ncol(bassi_malato))
# soglia_sani <- ceiling(prob * ncol(bassi_sano))
# keep_pos <- inf_pos & rowSums(bassi_malato > 0) >= soglia_malati
# keep_neg <- inf_neg & rowSums(bassi_sano > 0) >= soglia_sani
# keep_inf <- keep_pos | keep_neg
# index_ok <- rep(FALSE, length(logfoldchanges_bassi))
# finite_idx <- is.finite(logfoldchanges_bassi)
# index_ok[finite_idx] <- abs(logfoldchanges_bassi[finite_idx]) >= 1
# inf_idx <- is.infinite(logfoldchanges_bassi)
# index_ok[inf_idx] <- keep_inf[inf_idx]
# index_final <- index_alti
# index_final[index] <- index_ok
# trascr2D_euristica <- trascr2D_init[index_final, ]
# 
# # --- Filtraggio Trascrittomica 3D ---
# k_3D <- 50
# index_alti <- rowSums(trascr3D_init) >= k_3D
# index <- rowSums(trascr3D_init) < k_3D      # Seleziono i geni con num di conteggi inferiore a 50
# bassi_malato <- trascr3D_init[index, 1:10]  # estrazione dei conteggi per i malati
# bassi_sano <- trascr3D_init[index, 11:14]   # estrazione dei conteggi per i sani
# logfoldchanges_bassi <- log2(apply(bassi_malato, 1, mean)/apply(bassi_sano, 1, mean)) # calcolo di logfoldchange
# inf_pos  <- is.infinite(logfoldchanges_bassi) & logfoldchanges_bassi > 0
# inf_neg  <- is.infinite(logfoldchanges_bassi) & logfoldchanges_bassi < 0
# soglia_malati <- ceiling(prob * ncol(bassi_malato))
# soglia_sani <- ceiling(prob * ncol(bassi_sano))
# keep_pos <- inf_pos & rowSums(bassi_malato > 0) >= soglia_malati
# keep_neg <- inf_neg & rowSums(bassi_sano > 0) >= soglia_sani
# keep_inf <- keep_pos | keep_neg
# index_ok <- rep(FALSE, length(logfoldchanges_bassi))
# finite_idx <- is.finite(logfoldchanges_bassi)
# index_ok[finite_idx] <- abs(logfoldchanges_bassi[finite_idx]) >= 1
# inf_idx <- is.infinite(logfoldchanges_bassi)
# index_ok[inf_idx] <- keep_inf[inf_idx]
# index_final <- index_alti
# index_final[index] <- index_ok
# trascr3D_euristica <- trascr3D_init[index_final, ]
# 
# # ==============================================================================
# # 6. PLOT COMPARATIVO FINALE (OVERLAY A 3 VIE)
# # ==============================================================================
# get_log_cpm_df <- function(mat, label) {
#   mean_expr <- rowMeans(cpm(as.matrix(mat), log = TRUE, prior.count = 1))
#   data.frame(LogCPM = mean_expr, Filtro = label)
# }
# 
# 
# # Preparazione dataframe per 2D
# df_2D_comp <- rbind(
#   get_log_cpm_df(trascr2D_init, "1. Prima del filtro"),
#   get_log_cpm_df(trascr2D, "2. Filtro CPM"),
#   get_log_cpm_df(trascr2D_euristica, "3. Filtro Euristico")
# )
# 
# # Preparazione dataframe per 3D
# df_3D_comp <- rbind(
#   get_log_cpm_df(trascr3D_init, "1. Prima del filtro"),
#   get_log_cpm_df(trascr3D, "2. Filtro CPM"),
#   get_log_cpm_df(trascr3D_euristica, "3. Filtro Euristico")
# )
# 
# colori_filtro <- c(
#   "1. Prima del filtro"        = "grey", 
#   "2. Filtro CPM"  = "blue", 
#   "3. Filtro Euristico"   = "orange" 
# )
# 
# plot_3way_comparison <- function(df_tot, title_label, n_orig, n_valle, n_eur) {
#   ggplot(df_tot, aes(x = LogCPM, color = Filtro, fill = Filtro)) +
#     geom_density(alpha = 0.25, linewidth = 1) +
#     scale_color_manual(values = colori_filtro) +
#     scale_fill_manual(values = colori_filtro) +
#     theme_minimal() +
#     labs(
#       title = title_label,
#       subtitle = paste0("Geni iniziali: ", n_orig, 
#                         "\nFiltro CPM: ", n_valle, 
#                         " | Filtro Euristico: ", n_eur),
#       x = "Espressione Media (log2 CPM)",
#       y = "Densità",
#       color = "Strategia", fill = "Strategia"
#     ) +
#     theme(legend.position = "bottom")
# }
# 
# p_2D_final <- plot_3way_comparison(df_2D_comp, "Trascrittomica 2D - Confronto Filtri", 
#                                    nrow(trascr2D_init), nrow(trascr2D), nrow(trascr2D_euristica))
# 
# p_3D_final <- plot_3way_comparison(df_3D_comp, "Trascrittomica 3D - Confronto Filtri", 
#                                    nrow(trascr3D_init), nrow(trascr3D), nrow(trascr3D_euristica))
# 
# # Mostra i due plot comparativi affiancati
# p_2D_final + p_3D_final


dds2D <- DESeqDataSetFromMatrix(trascr2D, colData = meta, design = ~ Condition)
dds3D <- DESeqDataSetFromMatrix(trascr3D, colData = meta, design = ~ Condition)


rld2D <- rlog(dds2D, blind = TRUE) # blind = TRUE -> TUTTI i campioni trattati come NON correlati
trascr2D <- assay(rld2D)
rld3D <- rlog(dds3D, blind=TRUE)
trascr3D <- assay(rld3D)

###############################################################################################################
#### Grafico confronto distribuzione iniziale dei dati post filtraggio vs applicazione trasformazione rlog ####
###############################################################################################################
raw2D <- as.matrix(trascr2D)
raw3D <- as.matrix(trascr3D)
rld2D <- rlog(dds2D, blind = TRUE)
rld3D <- rlog(dds3D, blind = TRUE)

rlog2D <- assay(rld2D)
rlog3D <- assay(rld3D)

# Sui dati grezzi uso log2(x+1) solo per la visualizzazione

# 4) Grafico affiancato
par(mfrow = c(2, 2))
hist(log2(raw2D + 1), breaks = 60, main = "2D - grezzi", xlab = "log2(conteggio+1)", col = "blue")
hist(rlog2D,           breaks = 60, main = "2D - rlog",   xlab = "rlog",               col = "blue")
hist(log2(raw3D + 1), breaks = 60, main = "3D - grezzi", xlab = "log2(conteggio+1)", col = "orange")
hist(rlog3D,           breaks = 60, main = "3D - rlog",   xlab = "rlog",               col = "orange")
par(mfrow = c(1, 1))

################
## Proteomica ##
################
data <- read_excel("/Users/ceciliabergamini/Desktop/TESI/COL3A1_Data_&_Scripts/proteomics_Original Matrix_Col3A1_2D_3D.xlsx")
data <- data %>% drop_na('PG.Genes')
data <- data %>% 
  mutate(across(everything(), ~ gsub('"', '', .))) %>%   
  mutate(across(everything(), ~ gsub(",", ".", .))) %>% 
  mutate(across(everything(), ~ trimws(.)))             

prot_names <- data$PG.Genes
prot_names <- sub(";.*", "", prot_names)
prot2D <- as.data.frame(data[,1:14])
prot2D <- prot2D %>% mutate(across(everything(), as.numeric))

prot3D <- as.data.frame(data[,15:28])
prot3D <- prot3D %>% mutate(across(everything(), as.numeric))

new_order2D <- c(1, 7, 6, 5, 4, 8, 3, 9, 2, 10:14)
new_order3D <- c(1, 7, 6, 4, 5, 8, 3, 9, 2, 10:14)
prot2D <- prot2D[, new_order2D]
prot3D <- prot3D[, new_order3D]
colnames(prot2D) <- c("19-1597_2D", "20-3154_2D", "20-3156_2D", "20-4456_2D", "20-4222_2D", "22-0547_2D", 
                      "22-0548_2D", "GAU_2D", "SF_2D", "VL_2D", "C1A_2D", "C2B_2D", "C3B_2D", "RR_2D")
colnames(prot3D) <- c("19_1597_3D", "20_3154_3D", "20_3156_3D", "20_4456_3D", "20_4222_3D", "22_0547_3D",
                      "22_0548_3D", "GAU_3D", "SF_3D", "VL_3D", "C1A_3D", "C2B_3D", "C3B_3D", "RR_3D" )
rownames(prot2D) <- prot_names
rownames(prot3D) <- prot_names

# Trasformazione log2
prot2D <- log2(prot2D + 1)
prot3D <- log2(prot3D + 1)

prot2D <- prot2D[(rowSums(!is.na(prot2D[, 1:10])) > 6) | (rowSums(!is.na(prot2D[, 11:14])) > 2), ]
prot3D <- prot3D[(rowSums(!is.na(prot3D[, 1:10])) > 6) | (rowSums(!is.na(prot3D[, 11:14])) > 2), ] 

valori_osservati_2D <- as.vector(as.matrix(prot2D))
prot2D[is.na(prot2D)] <- NA # faccio imputation
bpca_model <- pcaMethods::pca(prot2D, method = "bpca", nPcs = 5)
prot2D <- completeObs(bpca_model)
valori_imputati_2D <- as.vector(prot2D)

prot3D[is.na(prot3D)] <- NA # faccio imputation
bpca_model <- pcaMethods::pca(prot3D, method = "bpca", nPcs = 5)
prot3D <- completeObs(bpca_model)

## Grafico 1: Trasformazione Log2 (prima vs dopo) ##
df_grezzi <- data.frame(Intensita = as.vector(as.matrix(prot2D)), Stato = "1. Grezzi")
prot2D_log <- log2(prot2D + 1)
prot3D_log <- log2(prot3D + 1)
df_log2   <- data.frame(Intensita = as.vector(as.matrix(prot2D_log)), Stato = "2. Log2(x+1)")

p1 <- ggplot(na.omit(df_grezzi), aes(x = Intensita)) +
  geom_density(fill = "steelblue", alpha = 0.6) +
  theme_bw() +
  labs(title = "A. Intensità Grezze (Proteomica 2D)", x = "Intensità del segnale", y = "Densità")

p2 <- ggplot(na.omit(df_log2), aes(x = Intensita)) +
  geom_density(fill = "darkgreen", alpha = 0.6) +
  theme_bw() +
  labs(title = "B. Trasformazione Log2(x + 1)", x = "log2(Intensità + 1)", y = "Densità")

p_log_comp <- p1 + p2
print(p_log_comp)

## Grafico 2: Confronto imputazione BPCA (valori osservati vs completi) ##
df_obs <- data.frame(Valore = valori_osservati_2D, Tipo = "Valori Osservati (con NA)")
df_imp <- data.frame(Valore = valori_imputati_2D, Tipo = "Dopo Imputazione BPCA")

df_tot_imp <- rbind(na.omit(df_obs), df_imp)

p_imp <- ggplot(df_tot_imp, aes(x = Valore, fill = Tipo)) +
  geom_density(alpha = 0.4) +
  scale_fill_manual(values = c("tomato", "steelblue")) +
  theme_bw() +
  labs(title = "Distribuzione Valori Proteici Prima e Dopo Imputazione BPCA (2D)",
       x = "Espressione Log2", 
       y = "Densità", 
       fill = "Stato") +
theme(
  plot.title   = element_text(size = 16, face = "bold"), # Titolo del grafico
  axis.title   = element_text(size = 14, face = "bold"), # Titoli degli assi X e Y
  axis.text    = element_text(size = 12),               # Numeri sugli assi X e Y
  legend.title = element_text(size = 14, face = "bold"), # Titolo della legenda
  legend.text  = element_text(size = 12),               # Testo delle voci in legenda
  legend.key.size = unit(0.8, "cm")                     # Dimensione dei quadratini colorati
)
plot(p_imp)

##############################################
####    Intersezione tra omiche 2D e 3D    ###
##############################################
common_genes <- intersect(rownames(trascr2D), rownames(trascr3D)) # per trascrittomica: seleziono le righe comuni a 2D e 3D
trascr2D_common <- trascr2D[common_genes, ]
trascr3D_common <- trascr3D[common_genes, ]

common_proteins <- intersect(rownames(prot2D), rownames(prot3D)) # per proteomica: seleziono le righe comuni a 2D e 3D
prot2D_common <- as.matrix(prot2D[common_proteins, ])
prot3D_common <- as.matrix(prot3D[common_proteins, ])

trascr <- cbind(trascr2D_common, trascr3D_common)
prot <- cbind(prot2D_common, prot3D_common)
groups <- c(rep("2D", 14), rep("3D", 14)) # Vettore dei gruppi: uno per ogni colonna (sample)

# 1^ scenario
dataList <- list(
  "transcriptomics" = trascr,
  "proteomics" = prot
)


##############################################
####                 MOFA                  ###
##############################################
### 1. APPLICAZIONE MOFA (Dataset completo con tutti i 14 campioni)

MOFAobj <- create_mofa(dataList, groups = groups)
plot_data_overview(MOFAobj)

dataOptions <- get_default_data_options(MOFAobj)
modelOptions <- get_default_model_options(MOFAobj)
modelOptions$num_factors <- 4
trainOptions <- get_default_training_options(MOFAobj)

MOFAobj <- prepare_mofa(
  object = MOFAobj,
  data_options = dataOptions,
  model_options = modelOptions,
  training_options = trainOptions
)


MOFAmodel_no_combat <- run_mofa(
  MOFAobj, 
  outfile = "MOFA_model_COL3A1.hdf5", 
  use_basilisk = TRUE
)
save(MOFAmodel_no_combat, file = "MOFA_model_COL3A1.RData")


### 2. APPLICAZIONE MOFA (Escludendo il paziente "RR" - colonna 14 in 2D e 3D)

trascr_meno_RR <- trascr[, -c(14, 28)] 
prot_meno_RR   <- prot[, -c(14, 28)] 
groups_meno_RR <- groups[-c(14, 28)]

# NOTA: ComBat è stato rimosso per essere coerenti con la nuova pipeline senza batch correction
dataList_meno_RR <- list(
  "transcriptomics" = trascr_meno_RR,
  "proteomics"      = prot_meno_RR
)

MOFAobj_meno_RR <- create_mofa(dataList_meno_RR, groups = groups_meno_RR)

dataOptions_meno_RR  <- get_default_data_options(MOFAobj_meno_RR)
modelOptions_meno_RR <- get_default_model_options(MOFAobj_meno_RR)
modelOptions_meno_RR$num_factors <- 4
trainOptions_meno_RR <- get_default_training_options(MOFAobj_meno_RR)

MOFAobj_meno_RR <- prepare_mofa(
  object = MOFAobj_meno_RR,
  data_options = dataOptions_meno_RR,
  model_options = modelOptions_meno_RR,
  training_options = trainOptions_meno_RR
)

MOFAmodel_meno_RR_no_combat <- run_mofa(
  MOFAobj_meno_RR, 
  outfile = "MOFA_meno_RR_model_COL3A1.hdf5", 
  use_basilisk = TRUE
)

save(MOFAmodel_meno_RR_no_combat, file = "MOFA_meno_RR_model_COL3A1.RData")


####################################################################
### 3. AGGIUNTA METADATI E VISUALIZZAZIONI

# Assegnazione corretta di metadati per il modello completo (duplicati per 2D e 3D)
MOFAmodel_no_combat@samples_metadata$condition <- rep(meta$Condition, 2)
MOFAmodel_no_combat@samples_metadata$group     <- groups

# Assegnazione metadati per il modello senza RR
MOFAmodel_meno_RR_no_combat@samples_metadata$condition <- rep(meta$Condition[-14], 2)
MOFAmodel_meno_RR_no_combat@samples_metadata$group     <- groups_meno_RR


### --- A. Matrice di Correlazione tra Fattori (Ortogonalità) ---
# Dimostra l'indipendenza e la non-ridondanza dei 4 fattori latenti estratti
plot_factor_cor(MOFAmodel_meno_RR_no_combat)


### --- B. Analisi di Arricchimento Funzionale dei Top Feature Weights ---
# 1. Estrazione dei nomi dei geni (ID Ensembl) presenti nel tuo modello MOFA
ensembl_ids <- features_names(MOFAmodel_meno_RR_no_combat)$transcriptomics

# 2. Conversione da Ensembl ID a Gene Symbol
gene_map <- mapIds(
  org.Hs.eg.db,
  keys      = ensembl_ids,
  column    = "SYMBOL",
  keytype   = "ENSEMBL",
  multiVals = "first"
)

# Rimuovi eventuali NA dal mapping
gene_map_clean <- na.omit(gene_map)
valid_ensembl  <- names(gene_map_clean)
valid_symbols  <- as.character(gene_map_clean)

# 3. Scarica i Gene Sets (Reactome) con i Gene Symbol
msig_reactome <- msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:REACTOME")
pathway_list  <- split(msig_reactome$gene_symbol, msig_reactome$gs_name)

# 4. Costruzione della matrice binaria usando gli ID ENSEMBL come colonne
# (Così le colonne della matrice coincidono esattamente con i nomi delle feature in MOFA)
feature_matrix <- t(sapply(pathway_list, function(pathway_genes) {
  as.numeric(valid_symbols %in% pathway_genes)
}))

colnames(feature_matrix) <- valid_ensembl

# 5. Esecuzione corretta di run_enrichment su MOFA2
enrichment_res <- run_enrichment(
  MOFAmodel_meno_RR_no_combat,
  view         = "transcriptomics",
  feature.sets = feature_matrix,
  sign         = "all"
)

# 6. Plot del grafico di enrichment (es. per il Factor 1 o Factor 2)
# Visualizzazione dei pathway arricchiti con il parametro corretto (max.pathways)
plot_enrichment(
  enrichment_res, 
  factor       = 1, 
  max.pathways = 15
) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 15, face = "bold"),
    axis.text.y  = element_text(size = 11, face = "bold"),
    axis.text.x  = element_text(size = 13),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12)
  )

### --- A. Plot Varianza Spiegata ---
p_var_full <- plot_variance_explained(MOFAmodel_no_combat) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 15, face = "bold"),
    axis.text    = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold")
  )
print(p_var_full)

p_var_no_rr <- plot_variance_explained(MOFAmodel_meno_RR_no_combat) +
  theme(
    plot.title   = element_text(size = 18, face = "bold"),
    axis.title   = element_text(size = 15, face = "bold"),
    axis.text    = element_text(size = 13, face = "bold"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text  = element_text(size = 12),
    strip.text   = element_text(size = 14, face = "bold")
  )
print(p_var_no_rr)


colori_mixomics <- c("CTR1" = "#3897C5", "EXP1" = "orange")

### --- B. Visualizzazione Singoli Fattori (Violin Plots) ---

# 1. Modello Completo (Tutti i pazienti)
p_full <- plot_factor(MOFAmodel_no_combat, 
                      factors = 1:4, 
                      color_by = "condition",
                      dot_size = 3.5, 
                      dodge = TRUE, 
                      legend = TRUE,
                      add_violin = TRUE, 
                      violin_alpha = 0.25) +
  scale_color_manual(values = colori_mixomics) +
  scale_fill_manual(values = colori_mixomics) +
  labs(title = "MOFA: Singoli Fattori (Tutti i pazienti)", color = "Condition", fill = "Condition") +
  theme(
    plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title      = element_text(size = 15, face = "bold"),
    axis.text.x     = element_text(size = 13, face = "bold"),
    axis.text.y     = element_text(size = 13),
    legend.title    = element_text(size = 14, face = "bold"),
    legend.text     = element_text(size = 13),
    legend.key.size = unit(0.9, "cm"),
    strip.text      = element_text(size = 14, face = "bold") # Intestazione Factor1, Factor2...
  )
print(p_full)

# 2. Modello Escluso Outlier RR
p_no_rr <- plot_factor(MOFAmodel_meno_RR_no_combat, 
                       factors = 1:4, 
                       color_by = "condition",
                       dot_size = 3.5, 
                       dodge = TRUE, 
                       legend = TRUE,
                       add_violin = TRUE, 
                       violin_alpha = 0.25) +
  scale_color_manual(values = colori_mixomics) +
  scale_fill_manual(values = colori_mixomics) +
  labs(title = "MOFA: Singoli Fattori (Escluso paziente RR)", color = "Condition", fill = "Condition") +
  theme(
    plot.title      = element_text(size = 18, face = "bold", hjust = 0.5),
    axis.title      = element_text(size = 15, face = "bold"),
    axis.text.x     = element_text(size = 13, face = "bold"),
    axis.text.y     = element_text(size = 13),
    legend.title    = element_text(size = 14, face = "bold"),
    legend.text     = element_text(size = 13),
    legend.key.size = unit(0.9, "cm"),
    strip.text      = element_text(size = 14, face = "bold")
  )
print(p_no_rr)


### --- C. Plot Fattori Combinati (Pair Plots Multi-Pannello) ---
make_pairs <- function(model, palette = colori_mixomics) {
  Z  <- get_factors(model, factors = 1:4, as.data.frame = TRUE)
  Zw <- pivot_wider(Z[, c("sample", "factor", "value")],
                    names_from = "factor", values_from = "value")
  md <- samples_metadata(model)[, c("sample", "condition")]
  df <- merge(Zw, md, by = "sample")
  df$condition <- factor(trimws(as.character(df$condition)), levels = names(palette))
  df$col <- palette[as.character(df$condition)]
  
  lower_fun <- function(data, mapping, ...) {
    ggplot(data, aes(x = !!mapping$x, y = !!mapping$y)) +
      geom_point(size = 3.5, color = data$col)
  }
  
  diag_fun <- function(data, mapping, ...) {
    ggplot(data, aes(x = !!mapping$x)) +
      geom_density(data = subset(data, condition == "CTR1"),
                   fill = palette[["CTR1"]], color = palette[["CTR1"]], alpha = 0.3) +
      geom_density(data = subset(data, condition == "EXP1"),
                   fill = palette[["EXP1"]], color = palette[["EXP1"]], alpha = 0.3)
  }
  
  upper_fun <- function(data, mapping, ...) {
    ggally_cor(data = data, mapping = mapping, size = 4, stars = FALSE,
               display_grid = FALSE)
  }
  
  ggpairs(
    df, columns = 2:5,
    mapping = aes(color = condition),
    lower = list(continuous = lower_fun),
    upper = list(continuous = upper_fun),
    diag  = list(continuous = diag_fun)
  ) +
    theme(
      strip.text = element_text(size = 16, face = "bold"),
      axis.text  = element_text(size = 11, face = "bold")
    )
}

p_comb_full  <- make_pairs(MOFAmodel_no_combat)
p_comb_no_rr <- make_pairs(MOFAmodel_meno_RR_no_combat)

print(p_comb_full)
print(p_comb_no_rr)


### --- D. Top Feature Weights (Loadings Bar Plots) ---
theme_weights <- theme(
  plot.title   = element_text(size = 16, face = "bold"),
  axis.title   = element_text(size = 14, face = "bold"),
  axis.text.y  = element_text(size = 12, face = "bold"), # Nomi dei geni/proteine ingranditi
  axis.text.x  = element_text(size = 12),
  strip.text   = element_text(size = 14, face = "bold")
)

# Modello Completo
plot_top_weights(MOFAmodel_no_combat, view = "transcriptomics", factor = c(1,4), nfeatures = 20) + theme_weights
plot_top_weights(MOFAmodel_no_combat, view = "proteomics",      factor = c(2,3), nfeatures = 20) + theme_weights

# Modello senza RR
plot_top_weights(MOFAmodel_meno_RR_no_combat, view = "transcriptomics", factor = c(1,4), nfeatures = 20) + theme_weights
plot_top_weights(MOFAmodel_meno_RR_no_combat, view = "proteomics",      factor = c(2,3), nfeatures = 20) + theme_weights


##############################################
####               mixOmics                ###
##############################################

### Calcolo di MCC per valutare le performance (singola omica) ###
compute_mcc_direct_verbose <- function(X, Y, ncomp_curr, keepX_vec, dist, kx_val, stability_threshold = 0.8) {
  n_samples <- length(Y)
  preds_loo <- character(n_samples)
  
  sample_names <- if (!is.null(rownames(X))) rownames(X) else paste0("S", seq_len(n_samples))
  
  cls_pos <- levels(Y)[1] # CTR1
  cls_neg <- levels(Y)[2] # EXP1
  
  # Lista per memorizzare le feature selezionate in ciascun fold LOO
  selected_features_per_fold <- vector("list", n_samples)
  
  for (s in seq_len(n_samples)) {
    X_tr <- X[-s, , drop = FALSE]
    Y_tr <- Y[-s]
    X_te <- X[s, , drop = FALSE]
    
    mod_sub <- suppressWarnings(
      mixOmics::splsda(X = X_tr, Y = Y_tr, ncomp = ncomp_curr, keepX = keepX_vec)
    )
    
    # 1. SALVATAGGIO FEATURE SELEZIONATE IN QUESTO FOLD
    # Estraiamo le variabili selezionate sulla componente corrente
    vars_fold <- mixOmics::selectVar(mod_sub, comp = ncomp_curr)$name
    selected_features_per_fold[[s]] <- vars_fold
    
    pred_sub <- suppressWarnings(stats::predict(mod_sub, newdata = X_te, dist = dist))
    preds_loo[s] <- as.character(pred_sub$class[[dist]][1, ncomp_curr])
    
    pred_class <- NULL
    
    # 1. Estraggo la matrice delle classi per la distanza richiesta
    if (!is.null(pred_sub$class) && !is.null(pred_sub$class[[dist]])) {
      mat_cls <- pred_sub$class[[dist]]
      
      # Gestione sia con intestazioni "comp1" che "1"
      col_target <- paste0("comp", ncomp_curr)
      if (!col_target %in% colnames(mat_cls)) {
        col_target <- as.character(ncomp_curr)
      }
      
      if (col_target %in% colnames(mat_cls)) {
        pred_class <- as.character(mat_cls[1, col_target])
      }
    }
    
    # Fallback su WeightedVote se presente (es. per modelli multi-blocco)
    if ((is.null(pred_class) || length(pred_class) == 0) && !is.null(pred_sub$WeightedVote)) {
      wv <- pred_sub$WeightedVote
      if (is.list(wv) && dist %in% names(wv)) {
        mat_wv <- wv[[dist]]
        col_target <- paste0("comp", ncomp_curr)
        if (!col_target %in% colnames(mat_wv)) col_target <- as.character(ncomp_curr)
        if (col_target %in% colnames(mat_wv)) pred_class <- as.character(mat_wv[1, col_target])
      }
    }
    
    # Assegnazione con controllo di sicurezza per evitare 'length zero'
    if (!is.null(pred_class) && length(pred_class) > 0) {
      preds_loo[s] <- pred_class
    } else {
      preds_loo[s] <- NA
    }
  }
  
  # 2. CALCOLO FREQUENZA DI STABILITÀ SULLE FEATURE
  all_selected_vec <- unlist(selected_features_per_fold)
  feature_counts   <- table(all_selected_vec)
  feature_stab     <- feature_counts / n_samples # Frequenza relativa (0.0 - 1.0)
  
  # Stabilità media del modello e conto delle feature stabili
  mean_stability  <- mean(feature_stab)
  stable_features <- names(feature_stab[feature_stab >= stability_threshold])
  num_stable      <- length(stable_features)
  
  # Matrice di confusione e Metriche
  preds_factor <- factor(preds_loo, levels = levels(Y))
  is_tp <- (Y == cls_pos & preds_factor == cls_pos)
  is_tn <- (Y == cls_neg & preds_factor == cls_neg)
  is_fp <- (Y == cls_neg & preds_factor == cls_pos)
  is_fn <- (Y == cls_pos & preds_factor == cls_neg)
  
  TP <- sum(is_tp); TN <- sum(is_tn)
  FP <- sum(is_fp); FN <- sum(is_fn)
  
  fp_pazienti <- sample_names[is_fp]
  fn_pazienti <- sample_names[is_fn]
  
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else 0
  ber  <- ((1 - sens) + (1 - spec)) / 2
  b_acc <- 1 - ber
  
  num <- (TP * TN) - (FP * FN)
  den <- sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- if (is.na(den) || den == 0) 0 else num / den
  
  return(list(
    kx = kx_val, TP = TP, TN = TN, FP = FP, FN = FN, 
    SENS = round(sens, 4), SPEC = round(spec, 4), B_ACC = round(b_acc, 4),
    BER = round(ber, 4), MCC = round(mcc, 4),
    MEAN_STAB = round(mean_stability, 4),       # Stabilità media globale (0-1)
    NUM_STABLE = num_stable,                    # Quante feature superano la soglia
    STABLE_FEATS = stable_features,             # Vettore dei nomi delle feature stabili
    FEATURE_STAB_TABLE = sort(feature_stab, decreasing = TRUE), # Tabella di stabilità di ciascuna
    FP_pazienti = if(length(fp_pazienti) > 0) paste(fp_pazienti, collapse=", ") else "-",
    FN_pazienti = if(length(fn_pazienti) > 0) paste(fn_pazienti, collapse=", ") else "-"
  ))
}

compute_block_mcc_loo <- function(X, Y, design, keepX, ncomp, dist = "centroids.dist") {
  n_samples <- length(Y)
  Y <- factor(Y)
  cls_pos <- levels(Y)[1] # CTR1
  cls_neg <- levels(Y)[2] # EXP1
  nomi_blocchi <- names(X)
  
  preds_loo <- character(n_samples)
  selected_vars <- lapply(nomi_blocchi, function(b) vector("list", n_samples))
  names(selected_vars) <- nomi_blocchi
  
  for (s in seq_len(n_samples)) {
    # Leave-One-Out split
    X_tr <- lapply(X, function(b) b[-s, , drop = FALSE])
    Y_tr <- Y[-s]
    X_te <- lapply(X, function(b) b[s, , drop = FALSE])
    
    if (length(unique(Y_tr)) < 2) {
      preds_loo[s] <- NA
      next
    }
    
    # Addestramento DIABLO
    mod_sub <- suppressMessages(suppressWarnings(
      tryCatch(mixOmics::block.splsda(X = X_tr, Y = Y_tr, ncomp = ncomp, keepX = keepX, design = design),
               error = function(e) NULL)
    ))
    
    if (is.null(mod_sub)) {
      preds_loo[s] <- NA
      next
    }
    
    # Salva feature selezionate nel fold per ciascun blocco
    for (b in nomi_blocchi) {
      vars_fold <- unlist(lapply(1:ncomp, function(c) {
        sel <- tryCatch(mixOmics::selectVar(mod_sub, block = b, comp = c), error = function(e) NULL)
        if (is.null(sel)) return(NULL)
        if (b %in% names(sel)) sel[[b]]$name else sel$name
      }))
      selected_vars[[b]][[s]] <- unique(vars_fold)
    }
    
    # Predizione LOO
    pred_sub <- suppressMessages(suppressWarnings(
      tryCatch(predict(mod_sub, newdata = X_te, dist = dist), error = function(e) NULL)
    ))
    
    if (!is.null(pred_sub)) {
      pred_class <- NULL
      
      # Strategia 1: Estrarre da WeightedVote (voto pesato di DIABLO)
      if (!is.null(pred_sub$WeightedVote)) {
        wv <- pred_sub$WeightedVote
        if (is.list(wv) && dist %in% names(wv)) {
          mat_dist <- wv[[dist]]
          pred_class <- as.character(mat_dist[1, ncomp])
        } else if (is.array(wv)) {
          pred_class <- as.character(wv[1, ncomp, dist])
        }
      }
      
      # Strategia 2: Fallback su MajorityVote
      if ((is.null(pred_class) || is.na(pred_class)) && !is.null(pred_sub$MajorityVote)) {
        mv <- pred_sub$MajorityVote
        if (is.list(mv) && dist %in% names(mv)) {
          pred_class <- as.character(mv[[dist]][1, ncomp])
        }
      }
      
      # Strategia 3: Fallback sul primo blocco se il voto d'insieme non è accessibile
      if (is.null(pred_class) || is.na(pred_class)) {
        if (!is.null(pred_sub$class) && !is.null(pred_sub$class[[dist]])) {
          mat_cls <- pred_sub$class[[dist]]
          pred_class <- as.character(mat_cls[1, ncomp])
        }
      }
      
      preds_loo[s] <- ifelse(!is.null(pred_class), pred_class, NA)
    } else {
      preds_loo[s] <- NA
    }
  }
  
  # Calcolo Confusion Matrix e Metriche
  preds_factor <- factor(preds_loo, levels = levels(Y))
  
  TP <- sum(!is.na(preds_factor) & Y == cls_pos & preds_factor == cls_pos)
  TN <- sum(!is.na(preds_factor) & Y == cls_neg & preds_factor == cls_neg)
  FP <- sum(!is.na(preds_factor) & Y == cls_neg & preds_factor == cls_pos)
  FN <- sum(!is.na(preds_factor) & Y == cls_pos & preds_factor == cls_neg)
  
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else 0
  ber  <- ((1 - sens) + (1 - spec)) / 2
  b_acc <- 1 - ber
  
  num <- (TP * TN) - (FP * FN)
  den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc <- if (is.na(den) || den == 0) 0 else num / den
  
  # Stabilità media delle feature
  stab_per_block <- sapply(nomi_blocchi, function(b) {
    all_feats <- unlist(selected_vars[[b]])
    if (length(all_feats) == 0) return(0)
    mean(table(all_feats) / n_samples)
  })
  
  return(list(
    TP = TP, TN = TN, FP = FP, FN = FN,
    SENS = round(sens, 4), SPEC = round(spec, 4), 
    BER = round(ber, 4), B_ACC = round(b_acc, 4), MCC = round(mcc, 4),
    STAB_GLOBAL = round(mean(stab_per_block), 4)
  ))
}

### Valutazione LOO single-omica: le predizioni sono per indice di campione (nessun problema di allineamento) ###
evaluate_splsda_loo <- function(X, Y, keepX_vec, ncomp, dist = "centroids.dist",
                                collect_features = FALSE) {
  Y <- factor(Y)
  n <- length(Y)
  cls_pos <- levels(Y)[1]
  cls_neg <- levels(Y)[2]
  sample_names <- if (!is.null(rownames(X))) rownames(X) else paste0("S", seq_len(n))
  
  preds <- rep(NA_character_, n)
  feats <- vector("list", n)
  first_error <- NA_character_
  
  for (s in seq_len(n)) {
    mod <- tryCatch(
      suppressMessages(suppressWarnings(
        mixOmics::splsda(X = X[-s, , drop = FALSE], Y = Y[-s], ncomp = ncomp, keepX = keepX_vec))),
      error = function(e) { if (is.na(first_error)) first_error <<- paste("splsda():", conditionMessage(e)); NULL })
    if (is.null(mod)) next
    
    if (collect_features) {
      feats[[s]] <- unique(unlist(lapply(seq_len(ncomp), function(k) mixOmics::selectVar(mod, comp = k)$name)))
    }
    
    pr <- tryCatch(
      suppressWarnings(stats::predict(mod, newdata = X[s, , drop = FALSE], dist = dist)),
      error = function(e) { if (is.na(first_error)) first_error <<- paste("predict():", conditionMessage(e)); NULL })
    if (is.null(pr)) next
    
    mat <- pr$class[[dist]]
    if (is.null(mat)) {
      if (is.na(first_error)) first_error <- paste0("predict() non restituisce class$", dist)
      next
    }
    preds[s] <- as.character(mat[1, ncomp])
  }
  
  pf    <- factor(preds, levels = levels(Y))
  is_tp <- !is.na(pf) & Y == cls_pos & pf == cls_pos
  is_tn <- !is.na(pf) & Y == cls_neg & pf == cls_neg
  is_fp <- !is.na(pf) & Y == cls_neg & pf == cls_pos
  is_fn <- !is.na(pf) & Y == cls_pos & pf == cls_neg
  TP <- sum(is_tp); TN <- sum(is_tn); FP <- sum(is_fp); FN <- sum(is_fn)
  
  sens <- if ((TP + FN) > 0) TP / (TP + FN) else 0
  spec <- if ((TN + FP) > 0) TN / (TN + FP) else 0
  ber  <- ((1 - sens) + (1 - spec)) / 2
  den  <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
  mcc  <- if (is.na(den) || den == 0) 0 else ((TP * TN) - (FP * FN)) / den
  
  feat_freq <- NULL
  if (collect_features) {
    n_ok <- sum(!vapply(feats, is.null, logical(1)))
    if (n_ok > 0) feat_freq <- sort(table(unlist(feats)) / n_ok, decreasing = TRUE)
  }
  
  list(TP = TP, TN = TN, FP = FP, FN = FN, SENS = sens, SPEC = spec,
       BER = ber, B_ACC = 1 - ber, MCC = mcc,
       FP_pazienti = if (any(is_fp)) paste(sample_names[is_fp], collapse = ", ") else "-",
       FN_pazienti = if (any(is_fn)) paste(sample_names[is_fn], collapse = ", ") else "-",
       n_na = sum(is.na(pf)), first_error = first_error, feat_freq = feat_freq)
}

### Selezione keepX basata su CONTEGGI di errori invece che su tolleranze percentuali ###
seleziona_keepX_conteggi <- function(tab, criterio = c("BER", "MCC"),
                                     max_extra_FN = 1, max_FN = NULL, max_extra_FP = NULL) {
  criterio <- match.arg(criterio)
  ord <- if (criterio == "BER") order(tab$BER, -tab$MCC, tab$keepX) else order(-tab$MCC, tab$BER, tab$keepX)
  rif <- tab[ord[1], ]
  
  ok <- rep(TRUE, nrow(tab))
  if (!is.null(max_extra_FN)) ok <- ok & tab$FN <= rif$FN + max_extra_FN
  if (!is.null(max_FN))       ok <- ok & tab$FN <= max_FN
  if (!is.null(max_extra_FP)) ok <- ok & tab$FP <= rif$FP + max_extra_FP
  
  if (!any(ok)) {
    cat("   [!] Nessuna configurazione rispetta la regola (il riferimento ha già FN > max_FN): uso il riferimento.\n")
    ok[ord[1]] <- TRUE
  }
  cand   <- tab[ok, , drop = FALSE]
  scelto <- cand[which.max(cand$keepX), ]   # a parità di regola: più feature (per la GO)
  list(keepX = scelto$keepX, rif = rif, scelto = scelto, n_candidati = nrow(cand))
}

tune_and_run_dual_criterion <- function(X, Y, 
                                        list.keepX = c(5, 10, 15, 20, 25, 30, 40, 50, 60, 70, 90, 120), 
                                        ncomp = 2, 
                                        dist = "mahalanobis.dist", 
                                        max_extra_FN = 1,       # FN in più ammessi rispetto al riferimento
                                        max_FN = NULL,          # num max di di FN (es. 2 = "inferiore a 3")
                                        max_extra_FP = NULL,    # FP in più ammessi (0 = non devono aumentare)
                                        tolerance_ber = 0.01,   # usate solo se le tre regole sopra sono tutte NULL
                                        tolerance_mcc = 0.02,
                                        stability_threshold = 0.8,
                                        nome_omica = "Omica") {
  
  Y <- factor(Y)
  usa_conteggi <- !(is.null(max_extra_FN) && is.null(max_FN) && is.null(max_extra_FP))
  descr_regola <- if (usa_conteggi) paste(c(
    if (!is.null(max_extra_FN)) sprintf("FN <= FN_rif + %d", max_extra_FN),
    if (!is.null(max_FN))       sprintf("FN <= %d", max_FN),
    if (!is.null(max_extra_FP)) sprintf("FP <= FP_rif + %d", max_extra_FP)), collapse = "  E  ")
  else sprintf("tolleranze: BER +%.3f | MCC -%.3f", tolerance_ber, tolerance_mcc)
  
  cat("\n========================================================================================================================\n")
  cat("   VALUTAZIONE MULTI-CRITERIO CON STABILITÀ (LOO-CV):", nome_omica, "\n")
  cat("========================================================================================================================\n")
  cat("Classi in Y:               ", paste(levels(Y), collapse = " [Positiva] vs "), "[Negativa]\n")
  cat("Distribuzione campioni Y:  ", paste(table(Y), collapse = " / "), "\n")
  cat("Soglia di stabilità applicata: >=", stability_threshold * 100, "%\n")
  cat("Regola di selezione:       ", descr_regola, "\n")
  
  set.seed(1234)
  
  ber_matrix  <- matrix(NA, nrow = length(list.keepX), ncol = ncomp, dimnames = list(list.keepX, paste0("comp", 1:ncomp)))
  mcc_matrix  <- matrix(NA, nrow = length(list.keepX), ncol = ncomp, dimnames = list(list.keepX, paste0("comp", 1:ncomp)))
  stab_matrix <- matrix(NA, nrow = length(list.keepX), ncol = ncomp, dimnames = list(list.keepX, paste0("comp", 1:ncomp)))
  
  keepX_ber_selected <- numeric(ncomp)
  keepX_mcc_selected <- numeric(ncomp)
  tabelle <- vector("list", ncomp)
  
  stampa_sel <- function(nome, s) {
    cat(sprintf("[%s] Riferimento: keepX = %s (BER %.4f | MCC %.4f | FN %d | FP %d) --> Scelto keepX = %s (BER %.4f | MCC %.4f | FN %d | FP %d | STAB %.4f) | candidati validi: %d\n",
                nome, s$rif$keepX, s$rif$BER, s$rif$MCC, s$rif$FN, s$rif$FP,
                s$scelto$keepX, s$scelto$BER, s$scelto$MCC, s$scelto$FN, s$scelto$FP, s$scelto$STAB, s$n_candidati))
  }
  
  for (c in 1:ncomp) {
    cat("\n------------------------------------------------------------------------------------------------------------------------\n")
    cat(sprintf("   >>> COMPONENTE %d - TABELLA VALUTAZIONE LOO CON STABILITÀ <<<\n", c))
    cat("------------------------------------------------------------------------------------------------------------------------\n")
    
    # 1. Prima esecuzione trasparente per calcolare i risultati di tutta la griglia della componente c
    tab_c <- data.frame()
    
    for (i in seq_along(list.keepX)) {
      kx <- list.keepX[i]
      keepX_test <- if (c == 1) kx else c(keepX_ber_selected[1:(c-1)], kx)
      
      res <- compute_mcc_direct_verbose(X = X, Y = Y, ncomp_curr = c, 
                                        keepX_vec = keepX_test, dist = dist, kx_val = kx,
                                        stability_threshold = stability_threshold)
      
      ber_matrix[i, c]  <- res$BER
      mcc_matrix[i, c]  <- res$MCC
      stab_matrix[i, c] <- res$MEAN_STAB
      
      tab_c <- rbind(tab_c, data.frame(
        keepX = kx, TP = res$TP, TN = res$TN, FP = res$FP, FN = res$FN,
        SENS = res$SENS, SPEC = res$SPEC, BER = res$BER, MCC = res$MCC, B_ACC = res$B_ACC,
        STAB = res$MEAN_STAB, N_STAB = res$NUM_STABLE,
        FP_pazienti = res$FP_pazienti, FN_pazienti = res$FN_pazienti, stringsAsFactors = FALSE))
    }
    tabelle[[c]] <- tab_c
    
    # 2. Calcolo dinamico delle larghezze delle colonne per ospitare I NOMI COMPLETI dei pazienti
    w_fp <- max(nchar(tab_c$FP_pazienti), nchar("Pazienti FP"))
    w_fn <- max(nchar(tab_c$FN_pazienti), nchar("Pazienti FN"))
    
    fmt_header <- sprintf("%%-6s | %%-2s | %%-2s | %%-2s | %%-2s | %%-6s | %%-6s | %%-6s | %%-6s | %%-6s | %%-6s | %%-5s | %%-%ds | %%-%ds\n", w_fp, w_fn)
    fmt_row    <- sprintf("%%-6d | %%-2d | %%-2d | %%-2d | %%-2d | %%-6.4f | %%-6.4f | %%-6.4f | %%-6.4f | %%-6.4f | %%-6.4f | %%-5d | %%-%ds | %%-%ds\n", w_fp, w_fn)
    
    # 3. Stampa della tabella con larghezza adattiva e senza troncamento
    cat(sprintf(fmt_header, "keepX", "TP", "TN", "FP", "FN", "SENS", "SPEC", "BER", "MCC", "B_ACC", "STAB", "#STAB", "Pazienti FP", "Pazienti FN"))
    cat(paste(rep("-", 82 + w_fp + w_fn), collapse = ""), "\n")
    
    for (i in seq_len(nrow(tab_c))) {
      r <- tab_c[i, ]
      cat(sprintf(fmt_row, r$keepX, r$TP, r$TN, r$FP, r$FN, 
                  r$SENS, r$SPEC, r$BER, r$MCC, r$B_ACC, 
                  r$STAB, r$N_STAB, r$FP_pazienti, r$FN_pazienti))
    }
    
    cat("\n--- RIEPILOGO COMPONENTE", c, "---\n")
    
    if (usa_conteggi) {
      sel_ber <- seleziona_keepX_conteggi(tab_c, "BER", max_extra_FN, max_FN, max_extra_FP)
      sel_mcc <- seleziona_keepX_conteggi(tab_c, "MCC", max_extra_FN, max_FN, max_extra_FP)
      keepX_ber_selected[c] <- sel_ber$keepX
      keepX_mcc_selected[c] <- sel_mcc$keepX
      stampa_sel("BER", sel_ber)
      stampa_sel("MCC", sel_mcc)
    } else {
      ber_vec <- ber_matrix[, c]
      kx_ber  <- max(as.numeric(names(ber_vec)[ber_vec <= min(ber_vec, na.rm = TRUE) + tolerance_ber]))
      mcc_vec <- mcc_matrix[, c]
      kx_mcc  <- max(as.numeric(names(mcc_vec)[mcc_vec >= max(mcc_vec, na.rm = TRUE) - tolerance_mcc]))
      keepX_ber_selected[c] <- kx_ber
      keepX_mcc_selected[c] <- kx_mcc
      cat(sprintf("[BER] Scelto keepX = %d (BER: %.4f | MCC: %.4f | Stabilità: %.4f)\n",
                  kx_ber, ber_vec[as.character(kx_ber)], mcc_matrix[as.character(kx_ber), c], stab_matrix[as.character(kx_ber), c]))
      cat(sprintf("[MCC] Scelto keepX = %d (MCC: %.4f | BER: %.4f | Stabilità: %.4f)\n",
                  kx_mcc, mcc_vec[as.character(kx_mcc)], ber_matrix[as.character(kx_mcc), c], stab_matrix[as.character(kx_mcc), c]))
    }
  }
  
  # Modelli Finali
  model_ber <- mixOmics::splsda(X = X, Y = Y, ncomp = ncomp, keepX = keepX_ber_selected)
  model_mcc <- mixOmics::splsda(X = X, Y = Y, ncomp = ncomp, keepX = keepX_mcc_selected)
  
  biomarkers_ber <- unique(unlist(lapply(1:ncomp, function(c) mixOmics::selectVar(model_ber, comp = c)$name)))
  biomarkers_mcc <- unique(unlist(lapply(1:ncomp, function(c) mixOmics::selectVar(model_mcc, comp = c)$name)))
  
  cat("\n========================================================================================================================\n")
  cat("keepX scelto (BER):", paste(keepX_ber_selected, collapse = " / "), "| feature uniche:", length(biomarkers_ber), "\n")
  cat("keepX scelto (MCC):", paste(keepX_mcc_selected, collapse = " / "), "| feature uniche:", length(biomarkers_mcc), "\n")
  cat("Feature in comune tra i due criteri:", length(intersect(biomarkers_ber, biomarkers_mcc)), "\n")
  cat("========================================================================================================================\n")
  
  return(list(
    model_ber = model_ber, keepX_ber = keepX_ber_selected, biomarkers_ber = biomarkers_ber,
    model_mcc = model_mcc, keepX_mcc = keepX_mcc_selected, biomarkers_mcc = biomarkers_mcc,
    grid_ber = ber_matrix, grid_mcc = mcc_matrix, grid_stab = stab_matrix,
    tabelle = tabelle, regola = descr_regola
  ))
}


estrai_riepilogo_singola_omica <- function(res_finali, nome_omica) {
  cat("\n====================================================\n")
  cat("OMICA:", nome_omica, "\n")
  cat("====================================================\n")
  for (comp in seq_along(res_finali$tabelle)) {
    riga <- res_finali$tabelle[[comp]][1, ]
    cat(sprintf("Comp %d (keepX=%d) -> TP=%d | TN=%d | FP=%d (%s) | FN=%d (%s) | BER=%.4f | MCC=%.4f\n",
                comp, riga$keepX, riga$TP, riga$TN,
                riga$FP, riga$FP_pazienti,
                riga$FN, riga$FN_pazienti,
                riga$BER, riga$MCC))
  }
}


### Funzione per valutare che p-value ha la soluzione ottimale trovata caso SINGOLA OMICA ###
check_overfitting_grid_permutations <- function(X, Y, 
                                              list.keepX = c(25, 300, 10000), 
                                              ncomp = 2, 
                                              dist = "mahalanobis.dist", 
                                              n_perm = 1001, 
                                              seed = 123,
                                              nome_omica = "Omica") {
  
  set.seed(seed)
  Y <- factor(Y)
  n_samples <- length(Y)
  cls_pos <- levels(Y)[1]
  cls_neg <- levels(Y)[2]
  
  cat("\n================================================================================================\n")
  cat("  PERMUTATION TEST SEQUENZIALE (1001 RUNS):", nome_omica, "\n")
  cat("================================================================================================\n")
  cat("Valori keepX testati per componente:", paste(list.keepX, collapse = ", "), "\n")
  cat("Numero permutazioni:", n_perm, "| Modalità: Sequenziale (Singolo Core)\n\n")
  
  # ------------------------------------------------------------------------------------------------
  # 1. CALCOLO MCC REALE PER CIASCUN KEEPX NELLA GRIGLIA (LOO-CV sui dati veri)
  # ------------------------------------------------------------------------------------------------
  cat("--> [1/2] Calcolo performance sui dati REALI...\n")
  mcc_real <- numeric(length(list.keepX))
  names(mcc_real) <- as.character(list.keepX)
  
  for (k in seq_along(list.keepX)) {
    kx <- list.keepX[k]
    keepX_vec <- rep(kx, ncomp)
    
    preds_loo <- character(n_samples)
    for (s in seq_len(n_samples)) {
      X_tr <- X[-s, , drop = FALSE]; Y_tr <- Y[-s]
      X_te <- X[s, , drop = FALSE]
      
      mod_sub <- suppressWarnings(mixOmics::splsda(X = X_tr, Y = Y_tr, ncomp = ncomp, keepX = keepX_vec))
      pred_sub <- tryCatch(suppressWarnings(predict(mod_sub, newdata = X_te, dist = dist)), error = function(e) NULL)
      
      if (!is.null(pred_sub) && !is.null(pred_sub$class[[dist]])) {
        mat_cls <- pred_sub$class[[dist]]
        col_t <- if (paste0("comp", ncomp) %in% colnames(mat_cls)) paste0("comp", ncomp) else as.character(ncomp)
        preds_loo[s] <- as.character(mat_cls[1, col_t])
      } else { preds_loo[s] <- NA }
    }
    
    p_fac <- factor(preds_loo, levels = levels(Y))
    tp <- sum(!is.na(p_fac) & Y == cls_pos & p_fac == cls_pos)
    tn <- sum(!is.na(p_fac) & Y == cls_neg & p_fac == cls_neg)
    fp <- sum(!is.na(p_fac) & Y == cls_neg & p_fac == cls_pos)
    fn <- sum(!is.na(p_fac) & Y == cls_pos & p_fac == cls_neg)
    
    num <- (tp * tn) - (fp * fn)
    den <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc_real[k] <- if (is.na(den) || den == 0) 0 else num / den
  }
  
  # ------------------------------------------------------------------------------------------------
  # 2. CALCOLO PERMUTAZIONI IN SEQUENZIALE
  # ------------------------------------------------------------------------------------------------
  cat("--> [2/2] Avvio permutazioni casuali...\n")
  
  # Matrice per salvare le MCC simulate (righe = permutazioni, colonne = keepX)
  mcc_null_matrix <- matrix(NA, nrow = n_perm, ncol = length(list.keepX))
  colnames(mcc_null_matrix) <- as.character(list.keepX)
  
  # Barra di avanzamento
  pb <- txtProgressBar(min = 0, max = n_perm, style = 3)
  
  for (p in 1:n_perm) {
    Y_perm <- sample(Y) # Permutazione casuale delle etichette
    
    for (k in seq_along(list.keepX)) {
      kx <- list.keepX[k]
      keepX_vec <- rep(kx, ncomp)
      preds_perm <- character(n_samples)
      
      for (s in seq_len(n_samples)) {
        X_tr <- X[-s, , drop = FALSE]; Y_tr <- Y_perm[-s]
        X_te <- X[s, , drop = FALSE]
        
        mod_perm <- suppressWarnings(mixOmics::splsda(X = X_tr, Y = Y_tr, ncomp = ncomp, keepX = keepX_vec))
        pred_p <- tryCatch(suppressWarnings(predict(mod_perm, newdata = X_te, dist = dist)), error = function(e) NULL)
        
        if (!is.null(pred_p) && !is.null(pred_p$class[[dist]])) {
          mat_p <- pred_p$class[[dist]]
          col_p <- if (paste0("comp", ncomp) %in% colnames(mat_p)) paste0("comp", ncomp) else as.character(ncomp)
          preds_perm[s] <- as.character(mat_p[1, col_p])
        } else { preds_perm[s] <- NA }
      }
      
      p_fac <- factor(preds_perm, levels = levels(Y_perm))
      tp_p <- sum(!is.na(p_fac) & Y_perm == cls_pos & p_fac == cls_pos)
      tn_p <- sum(!is.na(p_fac) & Y_perm == cls_neg & p_fac == cls_neg)
      fp_p <- sum(!is.na(p_fac) & Y_perm == cls_neg & p_fac == cls_pos)
      fn_p <- sum(!is.na(p_fac) & Y_perm == cls_pos & p_fac == cls_neg)
      
      num_p <- (tp_p * tn_p) - (fp_p * fn_p)
      den_p <- sqrt(as.numeric(tp_p + fp_p) * (tp_p + fn_p) * (tn_p + fp_p) * (tn_p + fn_p))
      mcc_null_matrix[p, k] <- if (is.na(den_p) || den_p == 0) 0 else num_p / den_p
    }
    
    # Aggiorna la barra di avanzamento
    setTxtProgressBar(pb, p)
    
    # Pulisce la memoria ogni 100 permutazioni per evitare accumuli di RAM
    if (p %% 100 == 0) gc(verbose = FALSE)
  }
  close(pb)
  
  # ------------------------------------------------------------------------------------------------
  # 3. CALCOLO P-VALUE EMPIRICO E TABELLA FINALE
  # ------------------------------------------------------------------------------------------------
  p_values <- numeric(length(list.keepX))
  names(p_values) <- as.character(list.keepX)
  
  for (k in seq_along(list.keepX)) {
    kx_str <- as.character(list.keepX[k])
    # p-value empirico con correzione +1 per evitare p-value pari a zero
    p_values[k] <- (sum(mcc_null_matrix[, kx_str] >= mcc_real[kx_str]) + 1) / (n_perm + 1)
  }
  
  cat("\n------------------------------------------------------------------------------------------------\n")
  cat("  RIEPILOGO SIGNIFICATIVITÀ E OVERFITTING (1001 Permutazioni)\n")
  cat("------------------------------------------------------------------------------------------------\n")
  res_table <- data.frame(
    keepX = list.keepX,
    MCC_Reale = round(mcc_real, 4),
    MCC_Medio_Casuale = round(colMeans(mcc_null_matrix), 4),
    p_value = round(p_values, 4),
    Significativo = ifelse(p_values <= 0.05, "SÌ (Segnale Reale)", "NO (Overfitting / Rumore)")
  )
  print(res_table)
  cat("================================================================================================\n\n")
  
  return(list(
    summary_table = res_table,
    mcc_real = mcc_real,
    p_values = p_values,
    null_distributions = mcc_null_matrix
  ))
}

### Grid Search INIZIALE: usa solo BER come misura di performance ###
grid_search_splsda <- function(X, Y, ncomp = 2,
                               tuning_ncomp = FALSE, 
                               design_values = c(0.05, 0.1, 0.2, 0.3),
                               valori_trascr = c(5, 10, 20, 30, 50, 70, 90, 120),
                               valori_prot   = c(10, 15, 20, 30, 50, 70),
                               dist          = "mahalanobis.dist",
                               n_probe       = 3,
                               tolerance_ber = 0.01) { # <--- Tolleranza BER (es. 0.01 = +1% di errore max)
  
  if (length(valori_trascr) != length(valori_prot)) {
    stop("[!] Errore: 'valori_trascr' e 'valori_prot' devono avere la stessa lunghezza!")
  }
  
  nomi_blocchi   <- names(X)
  risultati      <- data.frame()
  tempi          <- data.frame()
  n_combinazioni <- length(design_values) * length(valori_trascr)
  n_completate   <- 0
  
  cat("======================================================\n")
  cat("   GRID SEARCH OTTIMIZZATA (GO ORIENTED)\n")
  cat("======================================================\n")
  cat("Design values:", paste(design_values, collapse = ", "), "\n")
  cat("Tolleranza BER per Gene Ontology:", tolerance_ber, "\n")
  cat("Combinazioni totali:", n_combinazioni, "\n\n")
  
  coppie_keepX <- data.frame(keepX_trascr = valori_trascr, keepX_prot = valori_prot)
  combinazioni <- merge(data.frame(design_val = design_values), coppie_keepX, by = NULL)
  
  # --------------------------------------------------
  # STIMA TEMPO
  # --------------------------------------------------
  cat("--- Stima tempo su", n_probe, "combinazioni ---\n")
  set.seed(42)
  probe_idx <- sample(nrow(combinazioni), min(n_probe, nrow(combinazioni)))
  tempi_probe <- c()
  
  for (idx in probe_idx) {
    t_start <- proc.time()
    design_test <- matrix(combinazioni$design_val[idx], ncol = length(nomi_blocchi), 
                          nrow = length(nomi_blocchi), dimnames = list(nomi_blocchi, nomi_blocchi))
    diag(design_test) <- 0
    keepX_test <- list(transcriptomics = rep(combinazioni$keepX_trascr[idx], ncomp),
                       proteomics = rep(combinazioni$keepX_prot[idx], ncomp))
    
    mod <- tryCatch({ block.splsda(X, Y, ncomp = ncomp, keepX = keepX_test, design = design_test) }, error = function(e) NULL)
    if (!is.null(mod)) tryCatch({ perf(mod, validation = "loo", dist = dist, progressBar = FALSE) }, error = function(e) NULL)
    
    t_elapsed <- (proc.time() - t_start)["elapsed"]
    tempi_probe <- c(tempi_probe, t_elapsed)
  }
  
  cat("Tempo stimato totale:", round((mean(tempi_probe) * n_combinazioni) / 60, 1), "min. Premi INVIO per iniziare.\n")
  readline()
  
  # --------------------------------------------------
  # GRID SEARCH COMPLETA
  # --------------------------------------------------
  for (i in 1:nrow(combinazioni)) {
    val <- combinazioni$design_val[i]
    kT  <- combinazioni$keepX_trascr[i]
    kP  <- combinazioni$keepX_prot[i]
    n_completate <- n_completate + 1
    
    cat("\n[", n_completate, "/", n_combinazioni, "] design:", val, "| kT:", kT, "| kP:", kP)
    
    design_test <- matrix(val, ncol = length(nomi_blocchi), nrow = length(nomi_blocchi), dimnames = list(nomi_blocchi, nomi_blocchi))
    diag(design_test) <- 0
    keepX_test <- list(transcriptomics = rep(kT, ncomp), proteomics = rep(kP, ncomp))
    
    set.seed(1234)
    mod <- tryCatch({ block.splsda(X, Y, ncomp = ncomp, keepX = keepX_test, design = design_test) }, error = function(e) NULL)
    if (is.null(mod)) next
    
    perf_res <- tryCatch({ perf(mod, validation = "loo", dist = dist, progressBar = FALSE) }, error = function(e) NULL)
    if (is.null(perf_res)) next
    
    ber_vals <- perf_res$WeightedVote.error.rate[[dist]]["Overall.BER", ]
    
    for (c in 1:ncomp) {
      comp_name <- paste0("comp", c)
      if (comp_name %in% names(ber_vals)) {
        risultati <- rbind(risultati, data.frame(design_val = val, keepX_trascr = kT, keepX_prot = kP, 
                                                 n_variabili = kT + kP, comp = c, BER = round(ber_vals[comp_name], 4)))
      }
    }
  }
  
  # --------------------------------------------------
  # ANALISI E SELEZIONE PER GENE ONTOLOGY
  # --------------------------------------------------
  ncomp_ottimale <- ncomp
  if(tuning_ncomp && ncomp > 1) {
    analisi_comp <- risultati %>% group_by(comp) %>% summarise(BER_m = mean(BER), .groups = "drop")
    delta_ber <- analisi_comp$BER_m[1] - analisi_comp$BER_m[2]
    if(delta_ber < 0.02) {
      cat("\n[!] Delta BER basso (", round(delta_ber, 4), "). Suggerito ncomp = 1.\n")
      ncomp_ottimale <- 1
    }
  }
  
  risultati <- risultati %>% distinct(design_val, keepX_trascr, keepX_prot, comp, .keep_all = TRUE)
  
  sommario <- risultati %>%
    group_by(design_val, keepX_trascr, keepX_prot) %>%
    summarise(BER_totale = sum(BER), 
              BER_comp1 = BER[comp == 1], 
              BER_ultima = BER[comp == ncomp_ottimale], 
              n_variabili = dplyr::first(n_variabili), .groups = "drop")
  
  # 1. Ricerca del BER minimo assoluto
  min_ber <- min(sommario$BER_ultima, na.rm = TRUE)
  soglia_max_ber <- min_ber + tolerance_ber
  
  # 2. Isolamento dei candidati per GO (Tutte le configurazioni entro la tolleranza)
  candidati_go <- sommario %>%
    filter(BER_ultima <= soglia_max_ber) %>%
    arrange(desc(n_variabili), BER_ultima) # Priorità al numero MASSIMO di variabili, poi al BER
  
  # Si tiene quella con più variabili entro la tolleranza
  best_row <- candidati_go[1, ]
  
  cat("\n======================================================\n")
  cat("   SELEZIONE MODELLO OTTIMIZZATO PER GENE ONTOLOGY\n")
  cat("======================================================\n")
  cat("BER Minimo trovato:                 ", min_ber, "\n")
  cat("BER Massimo tollerato (+", tolerance_ber, "):     ", soglia_max_ber, "\n")
  cat("Configurazioni valide nel range:    ", nrow(candidati_go), "\n")
  cat("--> Scelta Finale: kT =", best_row$keepX_trascr, "| kP =", best_row$keepX_prot, 
      "| Totale feature =", best_row$n_variabili, "| BER =", best_row$BER_ultima, "\n")
  cat("======================================================\n")
  
  best_design_mat <- matrix(best_row$design_val, ncol = length(nomi_blocchi), nrow = length(nomi_blocchi), dimnames = list(nomi_blocchi, nomi_blocchi))
  diag(best_design_mat) <- 0
  
  best_keepX_list <- list(transcriptomics = rep(best_row$keepX_trascr, ncomp_ottimale),
                          proteomics = rep(best_row$keepX_prot, ncomp_ottimale))
  
  # Grafici
  p1 <- ggplot(sommario, aes(x = factor(keepX_trascr), y = factor(design_val), fill = BER_ultima)) +
    geom_tile() + geom_text(aes(label = round(BER_ultima, 3))) +
    scale_fill_gradient(low = "steelblue", high = "tomato") + theme_bw() + labs(title = "BER Heatmap - Transcriptomics")
  
  p2 <- ggplot(sommario, aes(x = factor(keepX_prot), y = factor(design_val), fill = BER_ultima)) +
    geom_tile() + geom_text(aes(label = round(BER_ultima, 3))) +
    scale_fill_gradient(low = "steelblue", high = "tomato") + theme_bw() + labs(title = "BER Heatmap - Proteomics")
  
  p3 <- plot_ly(data = sommario, x = ~keepX_trascr, y = ~keepX_prot, z = ~BER_ultima, 
                color = ~factor(design_val), type = "scatter3d", mode = "lines+markers") %>%
    layout(scene = list(xaxis = list(title = "Selected genes (kT)"), yaxis = list(title = "Selected proteins (kP)"), zaxis = list(title = "BER")))
  
  return(list(risultati = risultati, 
              sommario = sommario, 
              candidati_go = candidati_go, # restituzione di tutti i candidati per GO
              best_design = best_design_mat,
              best_keepX = best_keepX_list, 
              best_ncomp = ncomp_ottimale, 
              best_dist = dist,
              plot_trascr = p1, 
              plot_prot = p2, 
              plot_3d = p3))
}


### Funzione per valutare le predizioni e l'MCC in LOO per block.splsda ###
evaluate_block_multi_comp_loo <- function(X, Y, design, keepX, ncomp, dist = "centroids.dist") {
  Y <- factor(Y)
  livelli <- levels(Y)
  cls_pos <- livelli[1]; cls_neg <- livelli[2]
  
  vuoto_comp <- function() list(TP=0, TN=0, FP=0, FN=0, BER=1, MCC=0, n_na=length(Y))
  
  mod <- suppressMessages(suppressWarnings(
    tryCatch(mixOmics::block.splsda(X = X, Y = Y, ncomp = ncomp, keepX = keepX, design = design),
             error = function(e) NULL)))
  if (is.null(mod)) return(list(per_comp = replicate(ncomp, vuoto_comp(), simplify = FALSE),
                                first_error = "block.splsda fallito"))
  
  perf_res <- suppressMessages(suppressWarnings(
    tryCatch(mixOmics::perf(mod, validation = "loo", dist = dist, progressBar = FALSE),
             error = function(e) NULL)))
  if (is.null(perf_res)) return(list(per_comp = replicate(ncomp, vuoto_comp(), simplify = FALSE),
                                     first_error = "perf() fallito"))
  
  pred_nrep <- tryCatch(perf_res$class$nrep1, error = function(e) NULL)
  if (is.null(pred_nrep) || length(pred_nrep) == 0) {
    return(list(per_comp = replicate(ncomp, vuoto_comp(), simplify = FALSE),
                first_error = "class$nrep1 assente"))
  }
  nomi_blocchi <- names(pred_nrep)
  
  per_comp <- vector("list", ncomp)
  for (c in seq_len(ncomp)) {
    comp_col <- paste0("comp", c)
    pm <- tryCatch(sapply(nomi_blocchi, function(b) {
      mat_b <- pred_nrep[[b]][[dist]]
      if (is.null(mat_b)) return(rep(NA_character_, length(Y)))
      if (comp_col %in% colnames(mat_b)) mat_b[, comp_col] else mat_b[, c]
    }), error = function(e) NULL)
    
    if (is.null(pm) || (is.matrix(pm) && nrow(pm) == 0)) { per_comp[[c]] <- vuoto_comp(); next }
    if (is.null(dim(pm))) pm <- matrix(pm, ncol = 1, dimnames = list(names(pm), nomi_blocchi[1]))
    rownames(pm) <- rownames(pred_nrep[[nomi_blocchi[1]]][[dist]])
    
    ord_camp <- rownames(pm)
    majority <- apply(pm, 1, function(riga) {
      riga <- riga[!is.na(riga)]
      if (length(riga) == 0) return(NA_character_)
      tb <- sort(table(riga), decreasing = TRUE)
      if (length(tb) > 1 && tb[1] == tb[2]) riga[1] else names(tb)[1]
    })
    pred_class <- factor(majority, levels = livelli)
    Y_ord <- if (!is.null(names(Y))) Y[match(ord_camp, names(Y))] else Y[match(ord_camp, rownames(X[[1]]))]
    
    is_tp <- !is.na(pred_class) & Y_ord == cls_pos & pred_class == cls_pos
    is_tn <- !is.na(pred_class) & Y_ord == cls_neg & pred_class == cls_neg
    is_fp <- !is.na(pred_class) & Y_ord == cls_neg & pred_class == cls_pos
    is_fn <- !is.na(pred_class) & Y_ord == cls_pos & pred_class == cls_neg
    TP <- sum(is_tp); TN <- sum(is_tn); FP <- sum(is_fp); FN <- sum(is_fn)
    
    sens <- if ((TP+FN) > 0) TP/(TP+FN) else 0
    spec <- if ((TN+FP) > 0) TN/(TN+FP) else 0
    ber  <- ((1-sens)+(1-spec))/2
    den  <- sqrt(as.numeric(TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
    mcc  <- if (is.na(den) || den == 0) 0 else ((TP*TN)-(FP*FN))/den
    
    per_comp[[c]] <- list(TP=TP, TN=TN, FP=FP, FN=FN, BER=ber, MCC=mcc, n_na = sum(is.na(pred_class)))
  }
  list(per_comp = per_comp, first_error = NA_character_)
}

### Grid Search con altre metrice di performance oltre BER ###
grid_search_splsda_advanced <- function(X, Y, ncomp = 2,
                                        design_values = c(0.1, 0.2),
                                        valori_trascr = c(5, 7, 10, 30, 40, 50),
                                        valori_prot   = c(5, 7, 12, 20, 30, 40),
                                        dist          = "centroids.dist",
                                        tolerance_ber = 0.01,
                                        tolerance_mcc = 0.02,
                                        max_errori_positivi = NULL,
                                        extra_FN = NULL,
                                        fallback_se_vuoto = TRUE,
                                        ask_before_run = FALSE) {
  
  Y <- factor(Y)
  nomi_blocchi <- names(X)
  if (length(nomi_blocchi) < 2) {
    stop("[!] La lista X deve contenere almeno 2 blocchi omici con nomi validi.")
  }
  
  cls_pos <- levels(Y)[1]
  cls_neg <- levels(Y)[2]
  
  design_values <- sort(unique(design_values))
  valori_trascr <- sort(unique(valori_trascr))
  valori_prot   <- sort(unique(valori_prot))
  
  combinazioni <- expand.grid(design_val    = design_values,
                              keepX_trascr  = valori_trascr,
                              keepX_prot    = valori_prot,
                              KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  n_comb <- nrow(combinazioni)
  
  risultati_comp1 <- data.frame()
  risultati_comp2 <- data.frame()
  
  cat("\n========================================================================================================================\n")
  cat("   GRID SEARCH MULTI-OMICA (LOO Nativo mixOmics) - VALUTAZIONE PER COMPONENTE\n")
  cat("========================================================================================================================\n")
  cat("Classi in Y:             ", paste(levels(Y), collapse = " [Positiva] vs "), "[Negativa]\n")
  cat("Distribuzione campioni Y:", paste(table(Y), collapse = " / "), "\n")
  cat("Distanza utilizzata:     ", dist, "\n")
  cat("Numero Componenti:       ", ncomp, "\n")
  cat(sprintf("Combinazioni: %d design x %d kT x %d kP = %d modelli\n\n",
              length(design_values), length(valori_trascr), length(valori_prot), n_comb))
  
  if (ask_before_run) {
    cat("Premi INVIO per procedere (Esc per annullare)...\n")
    readline()
  }
  
  t_start <- proc.time()
  
  # Funzione interna per estrarre le metriche da una matrice di predizione
  calcola_metriche <- function(pm, Y, val, kT, kP, cls_pos, cls_neg) {
    ord_camp <- rownames(pm)
    majority_pred <- apply(pm, 1, function(riga) {
      riga <- riga[!is.na(riga)]
      if (length(riga) == 0) return(NA_character_)
      tb <- sort(table(riga), decreasing = TRUE)
      if (length(tb) > 1 && tb[1] == tb[2]) riga[1] else names(tb)[1]
    })
    pred_class <- factor(majority_pred, levels = levels(Y))
    Y_ord <- if (!is.null(names(Y))) Y[match(ord_camp, names(Y))] else Y[match(ord_camp, rownames(X[[1]]))]
    
    is_tp <- (Y_ord == cls_pos & pred_class == cls_pos)
    is_tn <- (Y_ord == cls_neg & pred_class == cls_neg)
    is_fp <- (Y_ord == cls_neg & pred_class == cls_pos)
    is_fn <- (Y_ord == cls_pos & pred_class == cls_neg)
    
    TP <- sum(is_tp); TN <- sum(is_tn); FP <- sum(is_fp); FN <- sum(is_fn)
    fp_pazienti <- ord_camp[is_fp]
    fn_pazienti <- ord_camp[is_fn]
    
    sens <- if ((TP + FN) > 0) TP / (TP + FN) else 0
    spec <- if ((TN + FP) > 0) TN / (TN + FP) else 0
    ber  <- ((1 - sens) + (1 - spec)) / 2
    b_acc <- 1 - ber
    
    num <- (TP * TN) - (FP * FN)
    den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
    mcc <- if (is.na(den) || den == 0) 0 else num / den
    
    fp_str <- if (length(fp_pazienti) > 0) paste(fp_pazienti, collapse = ", ") else "-"
    fn_str <- if (length(fn_pazienti) > 0) paste(fn_pazienti, collapse = ", ") else "-"
    
    data.frame(
      design_val = val, keepX_trascr = kT, keepX_prot = kP, n_variabili = kT + kP,
      TP = TP, TN = TN, FP = FP, FN = FN,
      BER = round(ber, 4), MCC = round(mcc, 4), SENS = round(sens, 4), SPEC = round(spec, 4), B_ACC = round(b_acc, 4),
      FP_pazienti = fp_str, FN_pazienti = fn_str, stringsAsFactors = FALSE
    )
  }
  
  for (i in seq_len(n_comb)) {
    val <- combinazioni$design_val[i]
    kT  <- combinazioni$keepX_trascr[i]
    kP  <- combinazioni$keepX_prot[i]
    
    el <- (proc.time() - t_start)["elapsed"]
    resid <- if (i > 1) el / (i - 1) * (n_comb - i + 1) / 60 else NA
    cat(sprintf("\r  [%d/%d] design=%.2f kT=%d kP=%d | residuo stimato: %s   ",
                i, n_comb, val, kT, kP, if (is.na(resid)) "..." else sprintf("%.1f min", resid)))
    
    design_test <- matrix(val, ncol = length(nomi_blocchi), nrow = length(nomi_blocchi), 
                          dimnames = list(nomi_blocchi, nomi_blocchi))
    diag(design_test) <- 0
    
    keepX_test <- list()
    keepX_test[[nomi_blocchi[1]]] <- rep(kT, ncomp)
    keepX_test[[nomi_blocchi[2]]] <- rep(kP, ncomp)
    
    mod <- suppressMessages(suppressWarnings(
      tryCatch(mixOmics::block.splsda(X = X, Y = Y, ncomp = ncomp, keepX = keepX_test, design = design_test),
               error = function(e) NULL)
    ))
    if (is.null(mod)) next
    
    perf_res <- suppressMessages(suppressWarnings(
      tryCatch(mixOmics::perf(mod, validation = "loo", dist = dist, progressBar = FALSE),
               error = function(e) NULL)
    ))
    if (is.null(perf_res)) next
    
    pred_nrep <- tryCatch(perf_res$class$nrep1, error = function(e) NULL)
    if (is.null(pred_nrep) || length(pred_nrep) == 0) next
    
    # Estrazione isolata per COMPONENTE 1
    pm_c1 <- tryCatch({
      pm <- sapply(names(pred_nrep), function(b) pred_nrep[[b]][[dist]][, "comp1"])
      rownames(pm) <- rownames(pred_nrep[[names(pred_nrep)[1]]][[dist]])
      pm
    }, error = function(e) NULL)
    
    # Estrazione isolata per COMPONENTE 2 (Cumulativa)
    pm_c2 <- tryCatch({
      pm <- sapply(names(pred_nrep), function(b) pred_nrep[[b]][[dist]][, "comp2"])
      rownames(pm) <- rownames(pred_nrep[[names(pred_nrep)[1]]][[dist]])
      pm
    }, error = function(e) NULL)
    
    if (!is.null(pm_c1)) risultati_comp1 <- rbind(risultati_comp1, calcola_metriche(pm_c1, Y, val, kT, kP, cls_pos, cls_neg))
    if (!is.null(pm_c2)) risultati_comp2 <- rbind(risultati_comp2, calcola_metriche(pm_c2, Y, val, kT, kP, cls_pos, cls_neg))
  }
  cat("\n\n")
  
  # Funzione interna per formattare e stampare ciascuna tabella
  stampa_tabella <- function(df, titolo) {
    cat("========================================================================================================================\n")
    cat("   ", titolo, "\n")
    cat("========================================================================================================================\n")
    w_fp <- max(nchar(df$FP_pazienti), nchar("Pazienti FP"))
    w_fn <- max(nchar(df$FN_pazienti), nchar("Pazienti FN"))
    
    header_fmt <- sprintf("%%-6s | %%-5s | %%-5s | %%-5s | %%-2s | %%-2s | %%-2s | %%-2s | %%-6s | %%-6s | %%-6s | %%-6s | %%-6s | %%-%ds | %%-%ds\n", w_fp, w_fn)
    row_fmt    <- sprintf("%%-6.2f | %%-5d | %%-5d | %%-5d | %%-2d | %%-2d | %%-2d | %%-2d | %%-6.4f | %%-6.4f | %%-6.4f | %%-6.4f | %%-6.4f | %%-%ds | %%-%ds\n", w_fp, w_fn)
    larg_tot   <- 100 + w_fp + w_fn
    
    cat(sprintf(header_fmt, "Design", "kT", "kP", "Tot", "TP", "TN", "FP", "FN", "SENS", "SPEC", "BER", "MCC", "B_ACC", "Pazienti FP", "Pazienti FN"))
    cat(paste(rep("-", larg_tot), collapse = ""), "\n")
    for (i in seq_len(nrow(df))) {
      r <- df[i, ]
      cat(sprintf(row_fmt, r$design_val, r$keepX_trascr, r$keepX_prot, r$n_variabili,
                  r$TP, r$TN, r$FP, r$FN, r$SENS, r$SPEC, r$BER, r$MCC, r$B_ACC,
                  r$FP_pazienti, r$FN_pazienti))
    }
    cat(paste(rep("-", larg_tot), collapse = ""), "\n\n")
  }
  
  # Stampa formattata delle due tabelle
  if (nrow(risultati_comp1) > 0) stampa_tabella(risultati_comp1, "TABELLA PERFORMANCE: COMPONENTE 1")
  if (nrow(risultati_comp2) > 0) stampa_tabella(risultati_comp2, "TABELLA PERFORMANCE: COMPONENTE 2 (CUMULATIVA)")
  
  # Selezione Modello Finale basata sulla Componente 2 (Cumulativa)
  risultati <- risultati_comp2
  min_FN_osservato <- min(risultati$FN, na.rm = TRUE)
  
  soglia_FN <- if (!is.null(extra_FN)) min_FN_osservato + extra_FN else if (!is.null(max_errori_positivi)) max_errori_positivi else min_FN_osservato
  
  candidati_go <- risultati %>% filter(FN <= soglia_FN) %>% arrange(BER, desc(MCC), desc(n_variabili))
  if (nrow(candidati_go) == 0 && fallback_se_vuoto) {
    candidati_go <- risultati %>% filter(FN <= min_FN_osservato) %>% arrange(BER, desc(MCC), desc(n_variabili))
  }
  
  best_row <- candidati_go[1, ]
  
  cat("========================================================================================================================\n")
  cat("   SELEZIONE MODELLO FINALE (Basato su Comp 1 + Comp 2)\n")
  cat("========================================================================================================================\n")
  cat(sprintf("--> Scelta Finale: design = %.2f | kT = %d | kP = %d | Tot Feats = %d | BER = %.4f | MCC = %.4f | FN = %d | FP = %d\n",
              best_row$design_val, best_row$keepX_trascr, best_row$keepX_prot, best_row$n_variabili, 
              best_row$BER, best_row$MCC, best_row$FN, best_row$FP))
  cat("FP (", cls_neg, "-->", cls_pos, "):", best_row$FP_pazienti, "\n")
  cat("FN (", cls_pos, "-->", cls_neg, "):", best_row$FN_pazienti, "\n")
  cat("========================================================================================================================\n")
  
  best_design_mat <- matrix(best_row$design_val, ncol = length(nomi_blocchi), nrow = length(nomi_blocchi), 
                            dimnames = list(nomi_blocchi, nomi_blocchi))
  diag(best_design_mat) <- 0
  
  best_keepX_list <- list()
  best_keepX_list[[nomi_blocchi[1]]] <- rep(best_row$keepX_trascr, ncomp)
  best_keepX_list[[nomi_blocchi[2]]] <- rep(best_row$keepX_prot, ncomp)
  
  return(list(
    risultati_comp1 = risultati_comp1,
    risultati_comp2 = risultati_comp2,
    candidati_go    = candidati_go,
    best_design     = best_design_mat,
    best_keepX      = best_keepX_list,
    best_ncomp      = ncomp,
    best_dist       = dist
  ))
}

### Singola Omica: Test di permutazione per valutare il p-value statistico del modello (MCC e BER) ###
permutation_test_singola_omica <- function(X, Y, keepX_vec, ncomp = 2, 
                                            dist = "mahalanobis.dist",
                                            n_perm = 10001, max_esaustivo = 5000, seed = 42,
                                            ask_before_run = TRUE,
                                            output_dir = NULL, nome_dataset = "modello") {
  Y <- factor(Y)
  n <- length(Y)
  if (nrow(X) != n) stop("[!] nrow(X) diverso da length(Y): servono i campioni sulle righe.")
  
  config_keepX <- rep_len(keepX_vec, ncomp)
  etichetta_keepX <- paste(config_keepX, collapse = "/")
  
  livelli <- levels(Y)
  cls_pos <- livelli[1]
  cls_neg <- livelli[2]
  n_pos   <- sum(Y == cls_pos)
  
  esaustivo <- choose(n, n_pos) <= max_esaustivo
  
  if (esaustivo) {
    comb  <- combn(n, n_pos)
    n_tot <- ncol(comb)
    get_Y <- function(i) { 
      lab <- rep(cls_neg, n)
      lab[comb[, i]] <- cls_pos
      factor(lab, levels = livelli) 
    }
  } else {
    set.seed(seed)
    n_tot <- n_perm
    get_Y <- function(i) sample(Y)
  }
  
  cat("\n========================================================================================================================\n")
  cat("   PERMUTATION TEST SEPARATO PER COMPONENTE:", nome_dataset, "\n")
  cat("========================================================================================================================\n")
  cat("Configurazione keepX:", etichetta_keepX, "| ncomp:", ncomp, "| dist:", dist, "\n")
  cat("Modalità:", if (esaustivo) sprintf("ESAUSTIVA (%d assegnazioni)", n_tot) else sprintf("Monte Carlo (%d run)", n_tot), "\n\n")
  
  # Funzione interna per estrarre predizioni LOO PER OGNI COMPONENTE
  evaluate_comp_loo <- function(X_mat, Y_vec) {
    preds_comp <- matrix(NA_character_, nrow = n, ncol = ncomp)
    colnames(preds_comp) <- paste0("comp", 1:ncomp)
    
    for (s in seq_len(n)) {
      mod_sub <- suppressWarnings(mixOmics::splsda(X = X_mat[-s, , drop = FALSE], 
                                                   Y = Y_vec[-s], 
                                                   ncomp = ncomp, 
                                                   keepX = config_keepX))
      
      pred_sub <- tryCatch(suppressWarnings(stats::predict(mod_sub, newdata = X_mat[s, , drop = FALSE], dist = dist)), 
                           error = function(e) NULL)
      
      if (!is.null(pred_sub) && !is.null(pred_sub$class[[dist]])) {
        mat_cls <- pred_sub$class[[dist]]
        for (c in 1:ncomp) {
          col_t <- if (paste0("comp", c) %in% colnames(mat_cls)) paste0("comp", c) else as.character(c)
          preds_comp[s, c] <- as.character(mat_cls[1, col_t])
        }
      }
    }
    
    # Calcolo metriche separatamente per ogni componente
    res_list <- vector("list", ncomp)
    for (c in 1:ncomp) {
      p_fac <- factor(preds_comp[, c], levels = levels(Y_vec))
      tp <- sum(!is.na(p_fac) & Y_vec == cls_pos & p_fac == cls_pos)
      tn <- sum(!is.na(p_fac) & Y_vec == cls_neg & p_fac == cls_neg)
      fp <- sum(!is.na(p_fac) & Y_vec == cls_neg & p_fac == cls_pos)
      fn <- sum(!is.na(p_fac) & Y_vec == cls_pos & p_fac == cls_neg)
      
      sens <- if ((tp + fn) > 0) tp / (tp + fn) else 0
      spec <- if ((tn + fp) > 0) tn / (tn + fp) else 0
      ber  <- ((1 - sens) + (1 - spec)) / 2
      
      num <- (tp * tn) - (fp * fn)
      den <- sqrt(as.numeric(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
      mcc <- if (is.na(den) || den == 0) 0 else num / den
      
      res_list[[c]] <- list(TP = tp, TN = tn, FP = fp, FN = fn, BER = ber, MCC = mcc)
    }
    return(res_list)
  }
  
  # 1. VALUTAZIONE REALE PER COMPONENTE
  cat("--> [1/2] Calcolo performance reali per singola componente...\n")
  t0 <- proc.time()
  reali_comp <- evaluate_comp_loo(X, Y)
  t_real <- (proc.time() - t0)["elapsed"]
  
  for (c in 1:ncomp) {
    r <- reali_comp[[c]]
    cat(sprintf("   [Comp %d] REALE | TP=%d TN=%d FP=%d FN=%d | BER=%.4f MCC=%.4f\n", 
                c, r$TP, r$TN, r$FP, r$FN, r$BER, r$MCC))
  }
  
  cat(sprintf("\nStima tempo permutazioni: %.1f min (%d run)\n", (t_real * n_tot) / 60, n_tot))
  if (ask_before_run) { cat("Premi INVIO per procedere (Esc per annullare)...\n"); readline() }
  
  # 2. RUN PERMUTAZIONI
  mcc_null <- matrix(NA_real_, nrow = n_tot, ncol = ncomp)
  ber_null <- matrix(NA_real_, nrow = n_tot, ncol = ncomp)
  colnames(mcc_null) <- paste0("Comp", 1:ncomp)
  colnames(ber_null) <- paste0("Comp", 1:ncomp)
  
  t_start <- proc.time()
  for (i in seq_len(n_tot)) {
    Yp <- get_Y(i)
    r_perm <- evaluate_comp_loo(X, Yp)
    
    for (c in 1:ncomp) {
      mcc_null[i, c] <- r_perm[[c]]$MCC
      ber_null[i, c] <- r_perm[[c]]$BER
    }
    
    if (i %% 50 == 0 || i == n_tot) {
      el <- (proc.time() - t_start)["elapsed"]
      cat(sprintf("   [%d/%d] Trascorso: %.1f min | Residuo: %.1f min\n", 
                  i, n_tot, el / 60, (el / i) * (n_tot - i) / 60))
    }
  }
  
  # 3. P-VALUE SEPARATI
  pval_fnc <- function(ge) if (esaustivo) ge / n_tot else (ge + 1) / (n_tot + 1)
  
  tabella_res <- data.frame()
  for (c in 1:ncomp) {
    r_real <- reali_comp[[c]]
    p_mcc  <- pval_fnc(sum(mcc_null[, c] >= r_real$MCC - 1e-12))
    p_ber  <- pval_fnc(sum(ber_null[, c] <= r_real$BER + 1e-12))
    
    tabella_res <- rbind(tabella_res, data.frame(
      Componente = paste("Comp", c),
      keepX = config_keepX[c],
      TP = r_real$TP, FN = r_real$FN, TN = r_real$TN, FP = r_real$FP,
      BER_Reale = round(r_real$BER, 4),
      p_BER = round(p_ber, 4),
      MCC_Reale = round(r_real$MCC, 4),
      p_MCC = round(p_mcc, 4),
      Significativo = ifelse(p_mcc <= 0.05, "SÌ (p <= 0.05)", "NO")
    ))
  }
  
  cat("\n========================================================================================================================\n")
  cat("   RIEPILOGO SIGNIFICATIVITÀ PER SINGOLA COMPONENTE (", n_tot, "assegnazioni )\n")
  cat("========================================================================================================================\n")
  print(tabella_res, row.names = FALSE)
  cat("========================================================================================================================\n")
  
  invisible(list(tabella = tabella_res, mcc_null = mcc_null, ber_null = ber_null))
}

### Integrazione Omica: Test di permutazione per valutare il p-value statistico del modello (MCC e BER) ###
permutation_test_multi_omica <- function(X, Y, design, keepX, ncomp, dist = "centroids.dist",
                                         n_perm = 1000, max_esaustivo = 5000, seed = 42,
                                         ask_before_run = TRUE,
                                         output_dir = NULL, nome_dataset = "modello") {
  Y <- factor(Y)
  n <- length(Y)
  livelli <- levels(Y)
  cls_pos <- livelli[1]; cls_neg <- livelli[2]
  n_pos <- sum(Y == cls_pos)
  
  esaustivo <- choose(n, n_pos) <= max_esaustivo
  if (esaustivo) {
    comb  <- combn(n, n_pos)
    n_tot <- ncol(comb)
    get_Y <- function(i) { lab <- rep(cls_neg, n); lab[comb[, i]] <- cls_pos; factor(lab, levels = livelli) }
  } else {
    set.seed(seed)
    n_tot <- n_perm
    get_Y <- function(i) sample(Y)
  }
  
  cat("\n========================================================================================================================\n")
  cat("   PERMUTATION TEST MULTI-OMICA PER COMPONENTE:", nome_dataset, "\n")
  cat("========================================================================================================================\n")
  cat("ncomp:", ncomp, "| dist:", dist, "\n")
  cat("Modalità:", if (esaustivo) sprintf("ESAUSTIVA (%d assegnazioni, p-value esatto)", n_tot)
      else sprintf("Monte Carlo (%d permutazioni, seed %d)", n_tot, seed), "\n\n")
  
  cat("--> [1/2] Calcolo performance reali per componente...\n")
  t0 <- proc.time()
  reale <- evaluate_block_multi_comp_loo(X, Y, design, keepX, ncomp, dist)
  t_real <- (proc.time() - t0)["elapsed"]
  if (!is.na(reale$first_error)) stop("[!] Modello reale fallito: ", reale$first_error)
  
  for (c in seq_len(ncomp)) {
    r <- reale$per_comp[[c]]
    if (r$n_na > 0) stop("[!] Comp ", c, ": ", r$n_na, " predizioni LOO mancanti sul modello reale.")
    cat(sprintf("   [Comp %d] REALE | TP=%d TN=%d FP=%d FN=%d | BER=%.4f MCC=%.4f\n",
                c, r$TP, r$TN, r$FP, r$FN, r$BER, r$MCC))
  }
  cat(sprintf("\nStima tempo totale: %.1f min (%d assegnazioni)\n", t_real * n_tot / 60, n_tot))
  if (ask_before_run) { cat("Premi INVIO per procedere (Esc per annullare)...\n"); readline() }
  
  cat("\n--> [2/2] Avvio permutazioni...\n")
  mcc_null <- matrix(NA_real_, n_tot, ncomp, dimnames = list(NULL, paste0("Comp", seq_len(ncomp))))
  ber_null <- mcc_null
  n_fail <- 0
  t_start <- proc.time()
  
  for (i in seq_len(n_tot)) {
    r_perm <- evaluate_block_multi_comp_loo(X, get_Y(i), design, keepX, ncomp, dist)
    for (c in seq_len(ncomp)) {
      pc <- r_perm$per_comp[[c]]
      mcc_null[i, c] <- pc$MCC
      ber_null[i, c] <- pc$BER
      if (pc$n_na > 0) n_fail <- n_fail + 1
    }
    if (i %% 25 == 0 || i == n_tot) {
      el <- (proc.time() - t_start)["elapsed"]
      cat(sprintf("   [%d/%d] Trascorso: %.1f min | Residuo stimato: %.1f min\n", i, n_tot, el/60, el/i*(n_tot-i)/60))
    }
  }
  if (n_fail > 0) cat("   [!]", n_fail, "valutazioni componente con predizioni mancanti\n")
  
  pval <- function(ge) if (esaustivo) ge / n_tot else (ge + 1) / (n_tot + 1)
  
  tabella <- data.frame()
  for (c in seq_len(ncomp)) {
    r <- reale$per_comp[[c]]
    p_mcc <- pval(sum(mcc_null[, c] >= r$MCC - 1e-12))
    p_ber <- pval(sum(ber_null[, c] <= r$BER + 1e-12))
    tabella <- rbind(tabella, data.frame(
      Componente = paste0("Comp", c),
      TP = r$TP, FN = r$FN, TN = r$TN, FP = r$FP,
      BER = round(r$BER, 4), p_BER = round(p_ber, 4),
      MCC = round(r$MCC, 4), p_MCC = round(p_mcc, 4),
      Significativo = ifelse(p_mcc <= 0.05, "SÌ (p<=0.05)", "NO"),
      stringsAsFactors = FALSE))
  }
  
  cat("\n========================================================================================================================\n")
  cat("   RIEPILOGO PER COMPONENTE (", n_tot, "assegnazioni", if (esaustivo) "- p-value esatti" else "", ")\n")
  cat("========================================================================================================================\n")
  print(tabella, row.names = FALSE)
  cat("========================================================================================================================\n")
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    lab_txt <- sprintf("Comp %d\np_MCC=%.3f | p_BER=%.3f", seq_len(ncomp), tabella$p_MCC, tabella$p_BER)
    lab <- factor(lab_txt, levels = lab_txt)
    m1 <- "MCC (più alto = meglio)"; m2 <- "BER (più basso = meglio)"
    df  <- rbind(data.frame(comp = rep(lab, each = n_tot), metrica = m1, valore = as.vector(mcc_null)),
                 data.frame(comp = rep(lab, each = n_tot), metrica = m2, valore = as.vector(ber_null)))
    dfr <- rbind(data.frame(comp = lab, metrica = m1, valore = tabella$MCC),
                 data.frame(comp = lab, metrica = m2, valore = tabella$BER))
    safe_plot(file.path(output_dir, paste0("Permutation_", nome_dataset, ".png")),
              width = 1500, height = 400 + 350 * ncomp, expr = {
                print(ggplot(df, aes(x = valore)) +
                        geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "white") +
                        geom_vline(data = dfr, aes(xintercept = valore), color = "tomato", linetype = "dashed", linewidth = 1.1) +
                        facet_grid(comp ~ metrica, scales = "free") + theme_bw() +
                        labs(title = paste("Distribuzione nulla per componente -", nome_dataset), x = "", y = "Frequenza"))
              })
  }
  
  invisible(list(tabella = tabella, mcc_null = mcc_null, ber_null = ber_null,
                 esaustivo = esaustivo, n_tot = n_tot))
}


prepare_param_second_phase <- function(risultati_fase1){
  
  best_design <- risultati_fase1$sommario$design_val[1]
  best_keep_trascr <- risultati_fase1$sommario$keepX_trascr[1]
  best_keep_prot <- risultati_fase1$sommario$keepX_prot[1]
  
  # Calcolo delle perturbazioni
  moltiplicatore <- c(-2, -1, 0, 1, 2)
  nuovi_valori_trascr <- pmax(best_keep_trascr + (moltiplicatore*3), 5)
  nuovi_valori_prot <- pmax(best_keep_prot + (moltiplicatore*2), 5)
  nuovi_design <- pmax(pmin(best_design + c(-0.05, 0, 0.05), 1), 0)
  
  return(list(
    nuovi_valori_trascr = nuovi_valori_trascr, 
    nuovi_valori_prot = nuovi_valori_prot,
    nuovi_design = nuovi_design))
}

### Sceglie la configurazione col miglior MCC (a parità di MCC: più feature per la GO, oppure meno) ###
scegli_config_migliore <- function(perm_grid, tie = c("max_feat", "min_feat")) {
  tie  <- match.arg(tie)
  tab  <- perm_grid$tabella
  cand <- tab[tab$MCC >= max(tab$MCC) - 1e-9, ]
  cand <- cand[order(if (tie == "max_feat") -cand$n_feat else cand$n_feat), ]
  scelta <- cand$keepX[1]
  cat("Configurazione scelta (MCC max =", max(tab$MCC), "):", scelta, "| candidati a pari MCC:", nrow(cand), "\n")
  perm_grid$configs[[match(scelta, tab$keepX)]]
}

run_splsda_pipeline <- function(X, Y, design, keepX, ncomp, dist) {
  set.seed(1234)
  Y <- factor(Y)
  livelli_originari <- levels(Y)
  cls_pos <- livelli_originari[1]
  cls_neg <- livelli_originari[2]
  
  nomi_blocchi <- names(X)
  
  cat("\n======================================================\n")
  cat("         INIZIO PIPELINE splsda MULTI-OMICA           \n")
  cat("======================================================\n")
  
  for (nome in nomi_blocchi) {
    cat("\nDimensioni", nome, ":", dim(X[[nome]]))
  }
  cat("\nNumero campioni in Y:", length(Y), "\n")
  cat("Classi:", paste(table(Y), collapse = "/"), "\n")
  cat("Componenti impostate:", ncomp, "\n") 
  cat("Distanza:", dist, "\n")
  
  # Sanitizzazione keepX 
  keepX_clean <- list()
  for (nome in nomi_blocchi) {
    if (nome %in% names(keepX)) {
      vals <- keepX[[nome]]
      if (length(vals) < ncomp) {
        keepX_clean[[nome]] <- rep(vals[1], ncomp)
      } else {
        keepX_clean[[nome]] <- vals[1:ncomp]
      }
    } else {
      stop(paste("[!] Errore: Blocco", nome, "non trovato nella lista keepX!"))
    }
  }
  keepX <- keepX_clean
  
  # 1. Modello finale
  cat("\n1. Addestramento modello finale.\n")
  final.model <- tryCatch({
    block.splsda(X = X, Y = Y, ncomp = ncomp, keepX = keepX, design = design)
  }, error = function(e) {
    stop("Modello finale fallito: ", conditionMessage(e))
  })
  
  # 2. Performance LOO sul modello finale
  cat("\n2. Valutazione performance (LOO).\n")
  perf.splsda <- tryCatch({
    perf(final.model, validation = "loo", dist = dist, progressBar = FALSE)
  }, error = function(e) {
    cat("   [!] perf() fallito:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  # 2b. TABELLA METRICHE COMPLETA (TP/TN/FP/FN/SENS/SPEC/BER/MCC/B_ACC)
  metrics <- list(TP = NA, TN = NA, FP = NA, FN = NA, SENS = NA, SPEC = NA, 
                  BER = NA, B_ACC = NA, MCC = NA)
  
  if (!is.null(perf.splsda)) {
    pred_nrep <- tryCatch(perf.splsda$class$nrep1, error = function(e) NULL)
    
    if (!is.null(pred_nrep) && length(pred_nrep) > 0) {
      nomi_blocchi_pred <- names(pred_nrep)
      comp_col <- paste0("comp", ncomp)
      
      pred_matrix <- tryCatch({
        pm <- sapply(nomi_blocchi_pred, function(b) {
          mat_b <- pred_nrep[[b]][[dist]]
          if (is.null(mat_b)) return(rep(NA_character_, length(Y)))
          if (comp_col %in% colnames(mat_b)) mat_b[, comp_col] else mat_b[, ncomp]
        })
        rownames(pm) <- rownames(pred_nrep[[nomi_blocchi_pred[1]]][[dist]])
        pm
      }, error = function(e) NULL)
      
      if (!is.null(pred_matrix) && nrow(pred_matrix) > 0) {
        sample_order <- rownames(pred_matrix)
        
        majority_pred <- apply(pred_matrix, 1, function(riga) {
          riga <- riga[!is.na(riga)]
          if (length(riga) == 0) return(NA_character_)
          tb <- sort(table(riga), decreasing = TRUE)
          if (length(tb) > 1 && tb[1] == tb[2]) riga[1] else names(tb)[1]
        })
        
        pred_class <- factor(majority_pred, levels = livelli_originari)
        
        if (!is.null(names(Y))) {
          Y_ord <- Y[match(sample_order, names(Y))]
        } else {
          Y_ord <- Y[match(sample_order, rownames(X[[1]]))]
        }
        
        conf_mat <- table(True = Y_ord, Pred = pred_class)
        
        TP <- ifelse(cls_pos %in% rownames(conf_mat) && cls_pos %in% colnames(conf_mat), conf_mat[cls_pos, cls_pos], 0)
        TN <- ifelse(cls_neg %in% rownames(conf_mat) && cls_neg %in% colnames(conf_mat), conf_mat[cls_neg, cls_neg], 0)
        FP <- ifelse(cls_neg %in% rownames(conf_mat) && cls_pos %in% colnames(conf_mat), conf_mat[cls_neg, cls_pos], 0)
        FN <- ifelse(cls_pos %in% rownames(conf_mat) && cls_neg %in% colnames(conf_mat), conf_mat[cls_pos, cls_neg], 0)
        
        sens <- if ((TP + FN) > 0) TP / (TP + FN) else 0
        spec <- if ((TN + FP) > 0) TN / (TN + FP) else 0
        ber  <- ((1 - sens) + (1 - spec)) / 2
        b_acc <- 1 - ber
        
        num <- (TP * TN) - (FP * FN)
        den <- sqrt(as.numeric(TP + FP) * (TP + FN) * (TN + FP) * (TN + FN))
        mcc <- if (is.na(den) || den == 0) 0 else num / den
        
        metrics <- list(
          TP = TP, TN = TN, FP = FP, FN = FN,
          SENS = round(sens, 4), SPEC = round(spec, 4),
          BER = round(ber, 4), B_ACC = round(b_acc, 4), MCC = round(mcc, 4)
        )
      }
    }
  }
  
  cat("\n   -> Metriche finali (voto di maggioranza multi-blocco, LOO):\n")
  cat(sprintf("      TP=%s TN=%s FP=%s FN=%s | SENS=%s SPEC=%s BER=%s B_ACC=%s MCC=%s\n",
              metrics$TP, metrics$TN, metrics$FP, metrics$FN,
              metrics$SENS, metrics$SPEC, metrics$BER, metrics$B_ACC, metrics$MCC))
  
  # 3. Estrazione biomarcatori 
  cat("\n3. Estrazione dei biomarcatori.\n")
  biomarkers <- list()
  for (nome in nomi_blocchi) {
    var_blocco <- c()
    for (i in 1:ncomp) {
      res_select <- tryCatch({
        sel <- selectVar(final.model, block = nome, comp = i)
        if (nome %in% names(sel)) sel[[nome]]$name
        else if (!is.null(sel$name)) sel$name
        else NULL
      }, error = function(e) NULL)
      if (!is.null(res_select)) var_blocco <- c(var_blocco, res_select)
    }
    biomarkers[[nome]] <- unique(var_blocco)
    cat("   ->", nome, ":", length(biomarkers[[nome]]), "variabili selezionate\n")
  }
  
  # 4. Stabilità 
  cat("\n4. Analisi di stabilità (LOO).\n")
  stability_results <- tryCatch({
    stability_splsda(X = X, Y = Y, design = design, ncomp = ncomp, keepX = keepX)
  }, error = function(e) {
    cat("   [!] stability_splsda fallito:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  return(list(
    model      = final.model,
    perf_data  = perf.splsda,
    metrics    = metrics,          
    biomarkers = biomarkers,
    stability  = stability_results
  ))
}

ispeziona_errori <- function(res_pipeline, Y, dist = "centroids.dist", ncomp = 1) {
  perf_obj <- res_pipeline$perf_data
  Y <- factor(Y)
  livelli <- levels(Y)
  cls_pos <- livelli[1]  # Classe positiva (es. CTR1)
  cls_neg <- livelli[2]  # Classe negativa (es. EXP1)
  
  pred_nrep <- perf_obj$class$nrep1
  nomi_blocchi <- names(pred_nrep)
  comp_col <- paste0("comp", ncomp)
  
  # Estrazione matrice predizioni per blocco
  pred_matrix <- sapply(nomi_blocchi, function(b) {
    mat_b <- pred_nrep[[b]][[dist]]
    if (comp_col %in% colnames(mat_b)) mat_b[, comp_col] else mat_b[, ncomp]
  })
  rownames(pred_matrix) <- rownames(pred_nrep[[nomi_blocchi[1]]][[dist]])
  
  # Voto di maggioranza 
  majority_pred <- apply(pred_matrix, 1, function(riga) {
    riga <- riga[!is.na(riga)]
    tb <- sort(table(riga), decreasing = TRUE)
    if (length(tb) > 1 && tb[1] == tb[2]) riga[1] else names(tb)[1]
  })
  
  sample_order <- rownames(pred_matrix)
  Y_ord <- if (!is.null(names(Y))) Y[match(sample_order, names(Y))] else Y
  
  df_check <- data.frame(
    Campione   = sample_order,
    Classe_Vera = as.character(Y_ord),
    Predetto   = as.character(majority_pred),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      Esito = dplyr::case_when(
        Classe_Vera == cls_pos & Predetto == cls_pos ~ "TP",
        Classe_Vera == cls_neg & Predetto == cls_neg ~ "TN",
        Classe_Vera == cls_neg & Predetto == cls_pos ~ "FP",
        Classe_Vera == cls_pos & Predetto == cls_neg ~ "FN"
      )
    )
  
  fp_names <- df_check$Campione[df_check$Esito == "FP"]
  fn_names <- df_check$Campione[df_check$Esito == "FN"]
  
  cat("\n--- DETTAGLIO ERRORI DI CLASSIFICAZIONE ---\n")
  cat("Classe Positiva (1° livello):", cls_pos, "| Classe Negativa (2° livello):", cls_neg, "\n")
  cat("FP (", length(fp_names), "):", if (length(fp_names) > 0) paste(fp_names, collapse = ", ") else "Nessuno", "\n")
  cat("FN (", length(fn_names), "):", if (length(fn_names) > 0) paste(fn_names, collapse = ", ") else "Nessuno", "\n")
  
  return(df_check)
}

stability_splsda <- function(X, Y, design, ncomp, keepX) {
  set.seed(1234)
  
  n <- length(Y)
  livelli      <- levels(droplevels(Y))
  
  if (length(livelli) != 2) stop("Y deve avere esattamente 2 livelli.")
  cat("\nClassi rilevate ->", livelli[1], "vs", livelli[2])
  cat("\nCampioni totali:", n, "| Iterazioni LOO:", n, "\n")
  
  stability_genes <- vector("list", ncomp)
  stability_prots <- vector("list", ncomp)
  loading_genes   <- vector("list", ncomp)
  loading_prots   <- vector("list", ncomp)
  
  for (comp in 1:ncomp) {
    stability_genes[[comp]] <- c()
    stability_prots[[comp]] <- c()
    loading_genes[[comp]]   <- data.frame()
    loading_prots[[comp]]   <- data.frame()
  }
  
  n_ok <- 0
  
  # --------------------------------------------------
  # LOOP LOO: ogni iterazione lascia fuori 1 campione
  # --------------------------------------------------
  for (i in 1:n) {
    cat("\nLOO iterazione", i, "- campione escluso:", i)
    
    idx   <- (1:n)[-i]           
    X_sub <- lapply(X, function(block) block[idx, , drop = FALSE])
    Y_sub <- droplevels(Y[idx])
    
    # Controllo che entrambe le classi siano presenti nel subset
    if (length(unique(Y_sub)) < 2) {
      cat("\n   [!] Subset con una sola classe - salto iterazione", i, "\n")
      next
    }
    
    # Modello sul subset
    model_sub <- tryCatch({
      block.splsda(X      = X_sub,
                   Y      = Y_sub,
                   ncomp  = ncomp,
                   keepX  = keepX,
                   design = design)
    }, error = function(e) {
      cat("\n   [!] Modello fallito - Errore:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (is.null(model_sub)) next
    n_ok <- n_ok + 1
    
    # Componenti effettivamente disponibili nel modello
    max_comp <- min(ncomp, ncol(model_sub$variates$transcriptomics))
    
    for (comp in 1:max_comp) {
      
      # GENI
      gene_res <- tryCatch(
        selectVar(model_sub, block = "transcriptomics", comp = comp),
        error = function(e) NULL
      )
      
      if (!is.null(gene_res) && length(gene_res$name) > 0) {
        stability_genes[[comp]] <- c(stability_genes[[comp]], gene_res$name)
        
        loading_genes[[comp]] <- rbind(
          loading_genes[[comp]],
          data.frame(
            variabile = gene_res$name,
            loading   = gene_res$value[, 1],
            row.names = NULL
          )
        )
      }
      
      # PROTEINE
      prot_res <- tryCatch(
        selectVar(model_sub, block = "proteomics", comp = comp),
        error = function(e) NULL
      )
      
      if (!is.null(prot_res) && length(prot_res$name) > 0) {
        stability_prots[[comp]] <- c(stability_prots[[comp]], prot_res$name)
        
        loading_prots[[comp]] <- rbind(
          loading_prots[[comp]],
          data.frame(
            variabile = prot_res$name,
            loading   = prot_res$value[, 1],
            row.names = NULL
          )
        )
      }
    }
  }
  
  # --------------------------------------------------
  # FREQUENZE (quante volte selezionato / iterazioni ok)
  # --------------------------------------------------
  cat("\n\nIterazioni LOO riuscite:", n_ok, "/", n, "\n")
  
  if (n_ok == 0) {
    warning("Nessuna iterazione LOO completata. Controlla i dati e i parametri.")
    return(NULL)
  }
  
  freq_genes <- lapply(stability_genes, function(x) {
    if (length(x) == 0)
      return(data.frame(variabile = character(), frequenza = numeric()))
    tb <- sort(table(x) / n_ok, decreasing = TRUE)
    data.frame(variabile = names(tb), frequenza = as.numeric(tb), row.names = NULL)
  })
  
  freq_prots <- lapply(stability_prots, function(x) {
    if (length(x) == 0)
      return(data.frame(variabile = character(), frequenza = numeric()))
    tb <- sort(table(x) / n_ok, decreasing = TRUE)
    data.frame(variabile = names(tb), frequenza = as.numeric(tb), row.names = NULL)
  })
  
  # Loading medi per variabile
  loading_medio_genes <- lapply(loading_genes, function(df) {
    if (nrow(df) == 0) return(data.frame())
    df %>%
      group_by(variabile) %>%
      summarise(loading_medio = mean(loading), .groups = "drop") %>%
      arrange(desc(abs(loading_medio)))
  })
  
  loading_medio_prots <- lapply(loading_prots, function(df) {
    if (nrow(df) == 0) return(data.frame())
    df %>%
      group_by(variabile) %>%
      summarise(loading_medio = mean(loading), .groups = "drop") %>%
      arrange(desc(abs(loading_medio)))
  })
  
  return(list(
    freq_genes          = freq_genes,
    freq_prots          = freq_prots,
    loading_medio_genes = loading_medio_genes,
    loading_medio_prots = loading_medio_prots,
    n_iterazioni_ok     = n_ok
  ))
}

### Preparazione parametri dopo aver applicato grid_search_splsda_advanced ###
prepare_param_second_phase_advanced <- function(risultati_fase1, step_trascr = 3, step_prot = 2, step_design = 0.05) {
  
  best_design      <- risultati_fase1$best_design[1, 2]  # valore off-diagonal della matrice
  best_keep_trascr <- risultati_fase1$best_keepX$transcriptomics[1]
  best_keep_prot   <- risultati_fase1$best_keepX$proteomics[1]
  
  moltiplicatore <- c(-2, -1, 0, 1, 2)
  
  nuovi_valori_trascr <- pmax(best_keep_trascr + (moltiplicatore * step_trascr), 5)
  nuovi_valori_prot   <- pmax(best_keep_prot   + (moltiplicatore * step_prot), 5)
  nuovi_design        <- pmax(pmin(best_design + c(-step_design, 0, step_design), 1), 0.01)
  
  cat("\n--- PARAMETRI RAFFINATI PER FASE 2 ---\n")
  cat("Punto di partenza (Fase 1): design =", best_design, "| kT =", best_keep_trascr, "| kP =", best_keep_prot, "\n")
  cat("Nuovi valori kT: ", paste(unique(nuovi_valori_trascr), collapse = ", "), "\n")
  cat("Nuovi valori kP: ", paste(unique(nuovi_valori_prot), collapse = ", "), "\n")
  cat("Nuovi valori design: ", paste(nuovi_design, collapse = ", "), "\n")
  
  return(list(
    nuovi_valori_trascr = nuovi_valori_trascr,
    nuovi_valori_prot   = nuovi_valori_prot,
    nuovi_design        = nuovi_design
  ))
}

## Plot ##
safe_plot <- function(filepath, width = 1000, height = 1000, res = 150, expr) {
  png(filepath, width = width, height = height, res = res)
  tryCatch({
    invisible(capture.output(force(expr)))
  }, error = function(e) {
    cat("   [!] Plot fallito:", conditionMessage(e), "\n")
  }, finally = {
    if (dev.cur() > 1) dev.off()
  })
}

generate_splsda_report_optimized <- function(risultati_pipeline, 
                                             nome_dataset, 
                                             Y, 
                                             colori_gruppi = c("CTR1" = "#377EB8", "EXP1" = "#E41A1C"),
                                             colori_blocchi = c('darkorchid', 'lightgreen', 'steelblue', 'tomato'),
                                             labels_blocchi = NULL) {
  
  output_dir <- paste0("Report_", nome_dataset)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  model <- risultati_pipeline$model
  ncomp_used <- ncol(model$variates[[1]])
  nomi_blocchi_raw <- names(model$X)
  
  # Gestione etichette custom dei blocchi omici
  if (is.null(labels_blocchi)) {
    labels_blocchi <- setNames(tools::toTitleCase(nomi_blocchi_raw), nomi_blocchi_raw)
  } else if (is.null(names(labels_blocchi))) {
    labels_blocchi <- setNames(labels_blocchi, nomi_blocchi_raw)
  }
  
  # Assegnazione dei colori ai gruppi Y
  livelli_Y <- levels(as.factor(Y))
  if (is.null(names(colori_gruppi))) {
    names(colori_gruppi) <- livelli_Y
  }
  
  cat("\n=======================================================\n")
  cat("Generazione report ottimizzato per:", nome_dataset, "\n")
  cat("=======================================================\n\n")
  
  # -------------------------------------------------------------------------
  # 0. Correlation Plot (DIABLO)
  # -------------------------------------------------------------------------
  safe_plot(file.path(output_dir, "0_Correlation_Plot.png"), width = 1000, height = 1000, expr = {
    plotDiablo(model, ncomp = 1, legend = TRUE)
  })
  
  # -------------------------------------------------------------------------
  # 1. Sample Separation GGPLOT2 (Estensione limiti e ggrepel)
  # -------------------------------------------------------------------------
  dist_used <- ifelse(!is.null(risultati_pipeline$metrics$dist), risultati_pipeline$metrics$dist, "mahalanobis.dist")
  pred_res <- tryCatch(predict(model, newdata = model$X, dist = dist_used), error = function(e) NULL)
  
  misclassified_samples <- character()
  if (!is.null(pred_res)) {
    pred_class <- pred_res$WeightedPredict[, , ncomp_used]
    if (is.matrix(pred_class)) pred_class <- pred_class[, 1]
    misclassified_samples <- names(Y)[which(as.character(pred_class) != as.character(Y))]
  }
  
  plots_list <- list()
  
  for (nome in nomi_blocchi_raw) {
    label_attuale <- labels_blocchi[nome]
    
    mat_var <- model$variates[[nome]]
    raw_names <- rownames(mat_var)
    clean_names <- gsub("_2D|_3D", "", raw_names)
    
    df_plot <- data.frame(
      Sample_Full = raw_names,
      Sample_Clean = clean_names,
      X1 = mat_var[, 1],
      X2 = if (ncomp_used >= 2) mat_var[, 2] else rep(0, nrow(mat_var)),
      Group = as.factor(Y),
      Is_Misclassified = raw_names %in% misclassified_samples | clean_names %in% gsub("_2D|_3D", "", misclassified_samples)
    )
    
    prop_var <- model$prop_expl_var[[nome]]
    var_comp1 <- round(prop_var[1] * 100, 1)
    var_comp2 <- if (ncomp_used >= 2) round(prop_var[2] * 100, 1) else 0
    
    p <- ggplot(df_plot, aes(x = X1, y = X2, color = Group, shape = Group)) +
      geom_point(size = 4, alpha = 0.85) +
      scale_shape_manual(values = c(16, 17)) + 
      scale_color_manual(values = colori_gruppi) +
      
      geom_text_repel(
        aes(
          label = Sample_Clean,
          color = ifelse(Is_Misclassified, "Misclassified", as.character(Group)),
          fontface = ifelse(Is_Misclassified, "bold", "bold"),
          size = ifelse(Is_Misclassified, 4.5, 3.8)
        ),
        box.padding = 0.35,
        point.padding = 0.3,
        max.overlaps = Inf,
        show.legend = FALSE
      ) +
      scale_size_identity() +
      
      scale_x_continuous(expand = expansion(mult = c(0.15, 0.15))) +
      scale_y_continuous(expand = expansion(mult = c(0.15, 0.15))) +
      
      theme_bw(base_size = 13) +
      theme(
        panel.grid.major = element_line(color = "grey92", size = 0.4),
        panel.grid.minor = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(face = "bold", size = 11, color = "black"),
        legend.position = "bottom",
        plot.margin = margin(15, 15, 15, 15)
      ) +
      labs(
        title = label_attuale,
        x = paste0("X-variate 1: ", var_comp1, "% expl. var"),
        y = paste0("X-variate 2: ", var_comp2, "% expl. var"),
        color = "Group", shape = "Group"
      )
    
    plots_list[[nome]] <- p
  }
  
  safe_plot(file.path(output_dir, "1_Sample_Separation.png"), width = 1400, height = 700, expr = {
    combined_plot <- wrap_plots(plots_list, ncol = length(nomi_blocchi_raw)) + 
      plot_annotation(
        title = paste("Separazione sPLS-DA con Evidenza dei Campioni Misclassificati -", nome_dataset),
        theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
      )
    print(combined_plot)
  })
  
  # -------------------------------------------------------------------------
  # 2. Correlation Circle Plot 
  # -------------------------------------------------------------------------
  if (ncomp_used >= 2) {
    safe_plot(file.path(output_dir, "2_Correlation_Circle.png"), width = 1100, height = 1100, expr = {
      
      # Margini e parametri grafici in grassetto
      op <- par(font = 2, font.axis = 2, font.lab = 2, cex.axis = 1.3, cex.lab = 1.4, lwd = 2.2, xpd = TRUE)
      
      pch_blocchi <- rep(19, length(nomi_blocchi_raw))
      
      plotVar(
        model, 
        comp = c(1, 2), 
        cutoff = 0.6, 
        var.names = FALSE, 
        style = "graphics", 
        pch = pch_blocchi, 
        cex = rep(1.8, length(nomi_blocchi_raw)), 
        col = colori_blocchi[1:length(nomi_blocchi_raw)],
        title = "Correlation Circle (Comp 1 & 2)",
        rad.in = 0.5
      )
      
      # Aggiunta manuale della legenda per i blocchi omici
      legend("bottomleft", 
             legend = unname(labels_blocchi[nomi_blocchi_raw]), 
             col = colori_blocchi[1:length(nomi_blocchi_raw)], 
             pch = 19, 
             pt.cex = 2, 
             cex = 1.2, 
             bg = "white", 
             box.col = "gray30",
             title = "Blocchi Omici")
      
      par(op)
    })
  }
  
  # -------------------------------------------------------------------------
  # 3. Boxplot Scores
  # -------------------------------------------------------------------------
  for (comp in 1:ncomp_used) {
    for (nome in nomi_blocchi_raw) {
      label_attuale <- labels_blocchi[nome]
      titolo_box <- paste0("Score Distribution - ", label_attuale, "\n(Component ", comp, ")")
      
      safe_plot(file.path(output_dir, paste0("3_Boxplot_", nome, "_Comp", comp, ".png")), width = 850, height = 650, expr = {
        scores_vec <- as.numeric(model$variates[[nome]][, comp])
        df <- data.frame(Score = scores_vec, Gruppo = as.factor(Y))
        
        p <- ggplot(df, aes(x = Gruppo, y = Score, fill = Gruppo)) +
          geom_boxplot(alpha = 0.7, outlier.shape = NA) + 
          geom_jitter(width = 0.15, size = 3, alpha = 0.8) +
          theme_bw(base_size = 13) + 
          ggtitle(titolo_box) +
          scale_fill_manual(values = colori_gruppi) +
          theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12, lineheight = 1.1),
            axis.title = element_text(face = "bold"),
            axis.text = element_text(face = "bold", color = "black")
          ) +
          labs(x = "Gruppo", y = paste("Variate Score (Comp", comp, ")"))
        print(p)
      })
    }
  }
  
  # -------------------------------------------------------------------------
  # 4. Loadings Plot (Etichette CTR1 a sinistra, EXP1 a destra)
  # -------------------------------------------------------------------------
  grp_left  <- if ("EXP1" %in% livelli_Y) "EXP1" else livelli_Y[2]
  grp_right <- if ("CTR1" %in% livelli_Y) "CTR1" else livelli_Y[1]
  
  for (nome in nomi_blocchi_raw) {
    label_attuale <- labels_blocchi[nome]
    
    for (comp in 1:ncomp_used) {
      titolo_loadings <- paste0("Loadings: ", label_attuale, " (Comp ", comp, ")")
      
      # 1. Estrazione sicura del data.frame dalla lista restituita da block.splsda
      df_load_raw <- plotLoadings(
        model,
        block = nome,
        comp = comp,
        method = 'median',
        contrib = 'max',
        plot = FALSE
      )
      
      if (is.list(df_load_raw) && !is.data.frame(df_load_raw)) {
        df_load <- if (nome %in% names(df_load_raw)) df_load_raw[[nome]] else df_load_raw[[1]]
      } else {
        df_load <- df_load_raw
      }
      df_load <- as.data.frame(df_load)
      
      max_features <- 50
      
      df_load <- df_load %>%
        dplyr::mutate(
          Feature = rownames(df_load),
          Abs_Importance = abs(importance)
        ) %>%
        dplyr::filter(importance != 0) %>%
        dplyr::slice_max(order_by = Abs_Importance, n = max_features, with_ties = FALSE) %>%
        dplyr::arrange(dplyr::desc(Abs_Importance)) %>%   # più grande = Row_ID 1 = in basso
        dplyr::mutate(Row_ID = dplyr::row_number())
      
      num_features <- nrow(df_load)
      if (num_features == 0) next
      
      # 2. Separazione dei gruppi: EXP1 a sinistra, CTR1 a destra
      df_left  <- df_load %>% dplyr::filter(GroupContrib == grp_left)
      df_right <- df_load %>% dplyr::filter(GroupContrib == grp_right)
      
      max_val <- max(df_load$Abs_Importance) * 1.08
      
      # Dimensioni dinamiche ad alta risoluzione
      altezza_dinamica <- max(8.5, num_features * 0.28)
      max_char         <- max(nchar(df_load$Feature))
      larghezza        <- max(12.5, 8.5 + max_char * 0.20)
      
      # Posizione x delle etichette di destra (bordo destro del pannello)
      x_right_lab <- max_val * 1.02
      
      p_load <- ggplot(df_load, aes(x = importance, y = Row_ID, fill = GroupContrib)) +
        geom_vline(xintercept = 0, color = "gray40", linewidth = 0.6) +
        # Nessun bordo nero sulle barre (rimossi color e linewidth)
        geom_col(width = 0.74, orientation = "y") +
        
        # Etichette del gruppo di DESTRA (CTR1 - blu) sul bordo destro del pannello
        geom_text(
          data = df_right,
          aes(x = x_right_lab, y = Row_ID, label = Feature),
          inherit.aes = FALSE,
          hjust = 0,
          size = 5,
          fontface = "bold",
          color = colori_gruppi[grp_right]
        ) +
        
        scale_fill_manual(values = colori_gruppi, name = "Gruppo Arricchito") +
        
        # Asse y: etichette del gruppo di SINISTRA (EXP1 - arancione)
        scale_y_continuous(
          breaks = df_left$Row_ID,
          labels = df_left$Feature,
          expand = c(0, 0)
        ) +
        
        coord_cartesian(
          xlim = c(-max_val, max_val),
          ylim = c(0.5, num_features + 0.5),
          expand = FALSE,
          clip = "off"
        ) +
        
        labs(
          title = titolo_loadings,
          subtitle = paste0("A sinistra: biomarcatori ", grp_left, " | A destra: biomarcatori ", grp_right),
          x = "Peso sul Vettore di Carico (Loading)",
          y = NULL
        ) +
        theme_bw(base_size = 13) +
        theme(
          plot.title         = element_text(size = 18, face = "bold", hjust = 0.5),
          plot.subtitle      = element_text(size = 15, hjust = 0.5, color = "gray30"),
          axis.text.y.left   = element_text(size = 13, face = "bold", color = colori_gruppi[grp_left]),
          axis.ticks.y       = element_blank(),
          axis.text.x        = element_text(size = 13, color = "black"),
          axis.title.x       = element_text(size = 13, face = "bold", margin = margin(t = 8)),
          panel.grid.major.y = element_line(color = "gray92", linewidth = 0.3),
          panel.grid.minor.y = element_blank(),
          legend.position    = "bottom",
          legend.title       = element_text(size = 15, face = "bold"),
          legend.text        = element_text(size = 13),
          plot.margin        = margin(10, max(160, max_char * 10), 10, 10)
        )
      
      ggsave(
        filename = file.path(output_dir, paste0("4_Loadings_", nome, "_Comp", comp, ".png")),
        plot = p_load,
        width = larghezza,
        height = altezza_dinamica,
        dpi = 600,
        limitsize = FALSE
      )
    }
  }
  
  # -------------------------------------------------------------------------
  # Preparazione sotto-modello filtrato:
  # - Network: Top 20 per blocco 
  # - Circos:  Top 40 per blocco 
  # -------------------------------------------------------------------------
  crea_modello_filtrato <- function(mod_orig, max_n) {
    mod_tmp <- mod_orig
    filtrato <- FALSE
    for (nome in nomi_blocchi_raw) {
      load_mat <- mod_tmp$loadings[[nome]]
      imp_max <- apply(abs(load_mat[, 1:ncomp_used, drop = FALSE]), 1, max)
      if (sum(imp_max > 0) > max_n) {
        filtrato <- TRUE
        top_idx <- order(imp_max, decreasing = TRUE)[1:max_n]
        keep_mask <- seq_len(nrow(load_mat)) %in% top_idx
        mod_tmp$loadings[[nome]][!keep_mask, ] <- 0
      }
    }
    list(model = mod_tmp, is_filtered = filtrato)
  }
  
  net_prep    <- crea_modello_filtrato(model, max_n = 20)
  circos_prep <- crea_modello_filtrato(model, max_n = 40)
  
  # -------------------------------------------------------------------------
  # 5. Network (Senza titolo interno che sovrascrive i nodi)
  # -------------------------------------------------------------------------
  safe_plot(file.path(output_dir, "5_Network.png"), width = 1450, height = 1400, expr = {
    
    op <- par(font = 2, xpd = NA, mar = c(4, 6, 2, 2))
    
    network(
      net_prep$model, 
      blocks = seq_along(nomi_blocchi_raw), 
      color.node = colori_blocchi[1:length(nomi_blocchi_raw)], 
      shape.node = c("circle", "rectangle")[1:length(nomi_blocchi_raw)],
      cutoff = 0.80,
      size.node = 0.20,
      cex.node = 0.90,
      lwd.edge = 2.4,
      show.color.key = TRUE
    )
    
    legend("bottomleft", 
           inset = c(-0.22, -0.03),
           legend = unname(labels_blocchi[nomi_blocchi_raw]), 
           fill = colori_blocchi[1:length(nomi_blocchi_raw)],
           border = "black",
           cex = 1.30, 
           bg = "white", 
           box.col = "gray30", 
           title = "Blocchi Omici")
    
    par(op)
  })
  
  # -------------------------------------------------------------------------
  # 6. CIM Heatmap
  # -------------------------------------------------------------------------
  safe_plot(file.path(output_dir, "6_CIM_Heatmap.png"), width = 1600, height = 1200, expr = {
    
    mod_cim <- model
    nuovi_nomi_cim <- unname(labels_blocchi[nomi_blocchi_raw])
    
    names(mod_cim$X)        <- nuovi_nomi_cim
    names(mod_cim$variates) <- nuovi_nomi_cim
    names(mod_cim$loadings) <- nuovi_nomi_cim
    if (!is.null(mod_cim$names$blocks)) {
      mod_cim$names$blocks  <- nuovi_nomi_cim
    }
    if (!is.null(mod_cim$prop_expl_var)) {
      names(mod_cim$prop_expl_var)[seq_along(nuovi_nomi_cim)] <- nuovi_nomi_cim
    }
    
    cimDiablo(mod_cim, 
              color.blocks = colori_blocchi[1:length(nomi_blocchi_raw)], 
              margins = c(10, 22))
  })
  
  # -------------------------------------------------------------------------
  # 7. Barplot di Stabilità Feature
  # -------------------------------------------------------------------------
  if (!is.null(risultati_pipeline$stability)) {
    for (comp in 1:ncomp_used) {
      for (nome in nomi_blocchi_raw) {
        label_attuale <- labels_blocchi[nome]
        
        freq_df <- NULL
        chiavi_possibili <- c(paste0("freq_", nome), nome, paste0("freq_", tolower(label_attuale)))
        for (k in chiavi_possibili) {
          if (!is.null(risultati_pipeline$stability[[k]][[comp]])) {
            freq_df <- risultati_pipeline$stability[[k]][[comp]]
            break
          }
        }
        
        if (is.null(freq_df) || nrow(freq_df) == 0) next
        
        freq_top <- head(freq_df[order(-freq_df$frequenza), ], 15)
        
        safe_plot(file.path(output_dir, paste0("7_Stability_", nome, "_Comp", comp, ".png")), 
                  width = 1000, height = 800, expr = {
                    p <- ggplot(freq_top, aes(x = reorder(variabile, frequenza), y = frequenza)) +
                      geom_col(fill = colori_blocchi[1], alpha = 0.85) +
                      geom_hline(yintercept = 0.8, linetype = "dashed", color = "red", size = 1) +
                      coord_flip() + theme_bw(base_size = 13) +
                      labs(title = paste0("Stabilità Feature - ", label_attuale, "\n(Comp ", comp, ")"),
                           x = "Molecola / Feature", 
                           y = "Frequenza di Selezione (CV / LOO)") +
                      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12))
                    print(p)
                  })
      }
    }
  }
  
  # -------------------------------------------------------------------------
  # 8. Circos Plot (Fino a 38 feature/blocco, senza sovrapposizione con le label)
  # -------------------------------------------------------------------------
  safe_plot(file.path(output_dir, "8_Circos.png"), width = 1750, height = 1750, expr = {
    
    op <- par(mar = c(2, 2, 2, 2), font = 2)
    
    # Copia del modello con i nomi dei blocchi aggiornati secondo labels_blocchi
    mod_circos <- circos_prep$model
    nuovi_nomi <- unname(labels_blocchi[nomi_blocchi_raw])
    
    names(mod_circos$X)            <- nuovi_nomi
    names(mod_circos$variates)     <- nuovi_nomi
    names(mod_circos$loadings)     <- nuovi_nomi
    if (!is.null(mod_circos$names$blocks)) {
      mod_circos$names$blocks      <- nuovi_nomi
    }
    if (!is.null(mod_circos$prop_expl_var)) {
      names(mod_circos$prop_expl_var)[seq_along(nuovi_nomi)] <- nuovi_nomi
    }
    
    # Elenco dei nomi dei biomarcatori da mantenere in font normale (font = 1)
    nomi_biomarcatori <- unlist(mod_circos$names$colnames, use.names = FALSE)
    
    # Wrapper locali per applicare font = 1 ai biomarcatori e font = 2 al resto
    text <- function(x, y = NULL, labels = seq_along(x), font = NULL, ...) {
      if (is.null(font)) {
        font <- ifelse(as.character(labels) %in% nomi_biomarcatori, 1, 2)
      }
      graphics::text(x = x, y = y, labels = labels, font = font, ...)
    }
    
    legend <- function(..., text.font = NULL) {
      if (is.null(text.font)) text.font <- 2
      graphics::legend(..., text.font = text.font)
    }
    
    circosPlot(
      mod_circos, 
      cutoff = 0.75, 
      line = TRUE,
      color.blocks = colori_blocchi[1:length(nomi_blocchi_raw)],
      color.Y = colori_gruppi,
      color.cor = c("#E41A1C", "#377EB8"),
      size.labels = 1.25,
      size.variables = 0.92,
      size.legend = 1.4
    )
    
    par(op)
  })
  
  cat("✔ Report completato con successo. Salvato in:", output_dir, "\n\n")
}

## Plot curva BER per sPLS-DA singola omica ##
plot_tuning_curves <- function(df_tuning, titolo_omica) {
  
  # Assicuriamoci che keepX sia trattato correttamente o come fattore o su scala log
  df_tuning$keepX_factor <- factor(df_tuning$keepX, levels = unique(df_tuning$keepX))
  
  # Tema grafico unificato per la tesi (Cornice chiusa + Testi grandi)
  tema_tuning <- theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.title = element_text(size = 13, face = "bold"),
      axis.text.x = element_text(size = 11, face = "bold", color = "black", angle = 45, hjust = 1),
      axis.text.y = element_text(size = 12, face = "bold", color = "black")
    )
  
  # 1. Grafico BER (Blu)
  p_ber <- ggplot(df_tuning, aes(x = keepX_factor, y = BER, group = 1)) +
    geom_line(color = "#377EB8", linewidth = 1.2) +
    geom_point(color = "#377EB8", size = 3) +
    labs(
      title = paste("BER -", titolo_omica),
      x = "Numero di Feature (keepX)",
      y = "Balanced Error Rate (BER)"
    ) +
    tema_tuning
  
  # 2. Grafico MCC (Verde scuro)
  p_mcc <- ggplot(df_tuning, aes(x = keepX_factor, y = MCC, group = 1)) +
    geom_line(color = "#006D2C", linewidth = 1.2) +
    geom_point(color = "#006D2C", size = 3) +
    labs(
      title = paste("MCC -", titolo_omica),
      x = "Numero di Feature (keepX)",
      y = "Matthews Correlation Coefficient (MCC)"
    ) +
    tema_tuning
  
  # Affiancamento dei due grafici con patchwork
  p_finale <- p_ber | p_mcc
  return(p_finale)
}

## Estrazione geni per fare GOEA ##
get_go_gene_lists <- function(risultati_pipeline, save_files = TRUE) {
  
  genes <- risultati_pipeline$biomarkers$transcriptomics
  proteins <- risultati_pipeline$biomarkers$proteomics
  
  cat("\n--- LISTE PRONTE PER GENE ONTOLOGY ---\n")
  cat("Geni estratti (Transcriptomics): ", length(genes), "\n")
  cat("Proteine estratte (Proteomics):  ", length(proteins), "\n")
  
  if (save_files) {
    write.lines(genes, "genes_for_GO.txt")
    write.lines(proteins, "proteins_for_GO.txt")
    cat("--> File salvati: 'genes_for_GO.txt' e 'proteins_for_GO.txt'\n")
  }
  
  return(list(
    genes = genes,
    proteins = proteins
  ))
}

# ----------------------------------------------------
# 1) Applicazione del metodo separatamente per 2D e 3D
# ----------------------------------------------------

## Preparazione e Normalizzazione dei dati ##

##################################################################################################
##################################################################################################
# Riduzione del numero dei geni sotto 10,000 (Median Absolute Deviation per RNA-seq data)
mad_trascr2D <- apply(trascr2D_common, 1, mad) 
mad_trascr3D <- apply(trascr3D_common, 1, mad) 

n_top <- 10000 
top_genes_2D <- names(sort(mad_trascr2D, decreasing = TRUE))[1:min(n_top, length(mad_trascr2D))]
top_genes_3D <- names(sort(mad_trascr3D, decreasing = TRUE))[1:min(n_top, length(mad_trascr3D))]

trascr_filtered2D <- trascr2D_common[top_genes_2D, ]
trascr_filtered3D <- trascr3D_common[top_genes_3D, ]

# Sfondi (Universe) per Gene Ontology
background_universe_2D <- top_genes_2D
background_universe_3D <- top_genes_3D

## Trascrittomica: scala per feature (per gene)
trascr_filtered2D_scaled <- t(scale(t(trascr_filtered2D),
                                    center = TRUE, scale = TRUE))
trascr_filtered3D_scaled <- t(scale(t(trascr_filtered3D),
                                    center = TRUE, scale = TRUE))

# # Proteomica: scala per feature (per proteina)
prot2D_scaled <- t(scale(t(prot2D_common),
                         center = TRUE, scale = TRUE))
prot3D_scaled <- t(scale(t(prot3D_common),
                         center = TRUE, scale = TRUE))

# Controllo: le medie per gene devono essere ~0 e le deviazioni standard ~1
cat("Media trascr2D (deve essere ~0):", 
    round(mean(rowMeans(trascr_filtered2D_scaled)), 4), "\n")
cat("SD trascr2D (deve essere ~1):", 
    round(mean(apply(trascr_filtered2D_scaled, 1, sd)), 4), "\n")

cat("Media prot2D (deve essere ~0):", 
    round(mean(rowMeans(prot2D_scaled)), 4), "\n")
cat("SD prot2D (deve essere ~1):", 
    round(mean(apply(prot2D_scaled, 1, sd)), 4), "\n")

X2D <- list(transcriptomics = t(trascr_filtered2D_scaled),
            proteomics     = t(prot2D_scaled))

X3D <- list(transcriptomics = t(trascr_filtered3D_scaled),
            proteomics     = t(prot3D_scaled))

Y <- as.factor(meta$Condition)

##########################################
### ----- Metodologia di analisi ----- ###
##########################################

#####  PER OGNI OMICA SINGOLARMENTE  #####
## OBIETTIVO: identificare il potenziale discriminatiorio di ogni singolo blocco
# 0) Analisi preliminare con PCA
# 1) modello PLS iniziale classico (si tengono TUTTI i biomarcatori, per stimare ncomp)
# 2) tuning sPLS: selezione del numero di componenti (BER)
# 3) keepX parameter: selezione del numero di variabili da tenere
# 4) modello finale
# 5) valutazione delle performance 
# P.S. le fasi 2), 3), 4) vengono svolte sia con loo che con k-fold cross validation


####  (0) Analisi Preliminare con PCA  

# Requisiti: dataset centrati e scalati
# 1^ STEP: runnarla con ncomp elevato (es 10), poi elbow rule
# 2^ STEP: valutare il clustering dei campioni nello spazio delle componenti principali

### Considero le 4 omiche separate ###
cat("\n\n>>> ESECUZIONE PCA <<<\n")

pca_trascr2D <- pca(X2D$transcriptomics, ncomp = 10, center = FALSE, scale = FALSE)
pca_prot2D <- pca(X2D$proteomics, ncomp = 10, center = FALSE, scale = FALSE)
pca_trascr3D <- pca(X3D$transcriptomics, ncomp = 10, center = FALSE, scale = FALSE)
pca_prot3D <- pca(X3D$proteomics, ncomp = 10, center = FALSE, scale = FALSE)

par(mfrow = c(2,2))
plot(pca_trascr2D, main = "Trascr2D")
plot(pca_prot2D, main = "Prot2D")
plot(pca_trascr3D, main = "Trascr3D")
plot(pca_prot3D, main = "Prot3D")
par(mfrow = c(1,1))
# Commento: si puà selezionare un numero di componenti pari a 2 per ogni omica,

nomi_campioni_2D <- rownames(X2D$transcriptomics)
nomi_campioni_3D <- rownames(X3D$transcriptomics)

etichette_sani_2D <- ifelse(meta$Condition == "CTR1", nomi_campioni_2D, "")
etichette_sani_3D <- ifelse(meta$Condition == "CTR1", nomi_campioni_3D, "")

tema_assi_grandi <- theme(
  axis.title = element_text(size = 14, face = "bold"), # Titoli degli assi (PC1, PC2, ecc.)
  axis.text  = element_text(size = 12)                # Valori numerici riportati sugli assi
)

# Applicazione ai grafici con geom_text_repel
p1_mod <- p1$graph + 
  geom_text_repel(aes(label = etichette_sani_2D), size = 5, color = "black", max.overlaps = Inf) + 
  tema_assi_grandi

p2_mod <- p2$graph + 
  geom_text_repel(aes(label = etichette_sani_2D), size = 5, color = "black", max.overlaps = Inf) + 
  tema_assi_grandi

p3_mod <- p3$graph + 
  geom_text_repel(aes(label = etichette_sani_3D), size = 5, color = "black", max.overlaps = Inf) + 
  tema_assi_grandi

p4_mod <- p4$graph + 
  geom_text_repel(aes(label = etichette_sani_3D), size = 5, color = "black", max.overlaps = Inf) + 
  tema_assi_grandi

# Composizione del plot complessivo con patchwork
plot_tot <- (p1_mod | p2_mod) /
  (p3_mod | p4_mod) +
  plot_layout(guides = "collect") & 
  theme(legend.title = element_text(size = 14, face = "bold"), # Ingrandisce anche la legenda condivisa
        legend.text  = element_text(size = 12))

print(plot_tot)

perc_pca_trascr2D <- round(pca_trascr2D$prop_expl_var$X * 100, 2)
perc_pca_prot2D <- round(pca_prot2D$prop_expl_var$X * 100, 2)
perc_pca_trascr3D <- round(pca_trascr3D$prop_expl_var$X * 100, 2)
perc_pca_prot3D <- round(pca_prot3D$prop_expl_var$X * 100, 2)

tab_perc_pca <- rbind(
  "Trascr2D" = perc_pca_trascr2D,
  "Prot2D" = perc_pca_prot2D,
  "Trascr3D" = perc_pca_trascr3D,
  "Prot3D" = perc_pca_prot3D
)

colnames(tab_perc_pca) <- paste0("Comp", 1:10)
print(tab_perc_pca)

# OSSERVAZIONE: i campioni NON si raggruppano in base alla condizione, si può
# dedurre che questa caratteristica non rappresenti il fattore principale di 
# distinzione fra i diversi gruppi di campioni.

##########################################
### Singola Omica: analisi esplorativa ###
##########################################

## splsda: mode = 'regression' -> per usare i dati di espressione genica per spiegare i dati clinici
## perf: "loo" (Leave-One-Out) -> perchè pochi campioni
# BER = (errore_sani + errori_malati) / 2   errore per la classificazione
# overall = conta solo gli errori senza considerare lo sbilanciamento delle classi

##############
## Trascr2D ##
##############
res_trascr2D <- tune_and_run_dual_criterion(
  X = t(trascr_filtered2D_scaled), Y = Y,
  list.keepX = c(3, 5, 10, 12, 14, 15, 60, 90, 100, 350, 500, 700, 1000, 5000),
  ncomp = 2, dist = "mahalanobis.dist",
  max_extra_FN = 1, max_FN = 2, max_extra_FP = NULL,
  nome_omica = "Trascrittomica 2D"
)
# COMMENTO:
# TP = 0, TN = 10, FP = 0, FN = 4 fino a keepX = 3
# TP = 1, TN = 10, FP = 0, FN = 3 da keepX = 5 fino a keepX = 14 => keepX = 14
# TP = 0, TN = 10, FP = 0, FN = 4 da keepX = 15 fino a keepX = 69
# TP = 0, TN = 7, FP = 3, FN = 4 da keepX = 70 fino a keepX = 99
# TP = 0, TN = 9, FP = 1, FN = 4 da keepX = 100 

### Valutazioni performance modello con calcolo del p-value ###
perm_trascr2D <- permutation_test_singola_omica(
  X = t(trascr_filtered2D_scaled), Y = Y,
  keepX_vec = c(14,14),
  ncomp = 2, dist = "mahalanobis.dist",
  output_dir = "Report_Trascr2D", nome_dataset = "Trascr2D"
)
# Comp1: p-BER = 0.1788, p-MCC = 0.1538 => NON significativo
perm_trascr2D$tabella

# Recupero dei risultati del modello e delle performance
res_trascr2D_finali <- tune_and_run_dual_criterion(
  X = t(trascr_filtered2D_scaled), Y = Y,
  list.keepX = 14, 
  ncomp = 2, 
  dist = "mahalanobis.dist",
  nome_omica = "Trascrittomica 2D (keepX = 14)"
)

# Visualizza i geni scelti
geni_per_GO_2D <- res_trascr2D_finali$biomarkers_mcc
head(geni_per_GO_2D)
cat("Totale geni estratti per la GOEA:", length(geni_per_GO_2D), "\n") # 28
writeLines(geni_per_GO_2D, "geni_trascr2D_per_GO.txt")
## comp1: BER: 0.375  MCC: 0.4385  ##


##############
##  Prot2D  ##
##############
res_prot2D <- tune_and_run_dual_criterion(
  X = t(prot2D_scaled),
  Y = Y,
  list.keepX = c(8, 30, 70, 90, 97, 100, 120, 350),  # 97 (3TP, 10TN, 0FP, 1FN)
  ncomp = 2, dist = "mahalanobis.dist",
  max_extra_FN = 1, max_FN = 2, max_extra_FP = NULL,
  nome_omica = "Proteomica 2D"
)
# COMMENTO:
# TP = 2, TN = 10, FP = 0, FN = 2 fino a keepX = 30
# TP = 3, TN = 10, FP = 0, FN = 1 da keepX = 31 fino a keepX = 97 => keepX = 97
# TP = 3, TN = 9, FP = 1, FN = 1 da keepX = 98 

### Valutazioni performance modello con calcolo del p-value ###
perm_prot2D <- permutation_test_singola_omica(
  X = t(prot2D_scaled), Y = Y,
  keepX_vec = 97,
  ncomp = 2, dist = "mahalanobis.dist",
  output_dir = "Report_Prot2D", nome_dataset = "Prot2D"
) 
# Comp1: p-BER = 0.011  p-MCC = 0.011 => test significativo
perm_prot2D$tabella

# Recupero dei risultati del modello e delle performance
res_prot2D_finali <- tune_and_run_dual_criterion(
  X = t(prot2D_scaled), Y = Y,
  list.keepX = 97, 
  ncomp = 2, 
  dist = "mahalanobis.dist",
  nome_omica = "Proteomica 2D (keepX = 97)"
)

# Visualizza i geni scelti
prot_per_GO_2D <- res_prot2D_finali$biomarkers_mcc
head(prot_per_GO_2D)
cat("Totale proteine estratte per la GOEA:", length(prot_per_GO_2D), "\n") # 194
writeLines(prot_per_GO_2D, "prot_prot2D_per_GO.txt")
## comp1: BER: 0.175 MCC: 0.65  ##

##############
## Trascr3D ##
##############
res_trascr3D <- tune_and_run_dual_criterion(
  X = t(trascr_filtered3D_scaled),
  Y = Y,
  list.keepX = c(15, 25, 60, 100, 300, 400, 700, 1000, 2000, 3000, 5000, 10000), 
  ncomp = 2, dist = "mahalanobis.dist",
  max_extra_FN = 1, max_FN = 2, max_extra_FP = NULL,
  nome_omica = "Trascrittomica 3D"
)
# COMMENTO:
# TP = 1, TN = 8, FP = 2, FN = 3 fino a keepX = 120
# TP = 1, TN = 9, FP = 1, FN = 3 da keepX = 200 fino a keepX = 10,000 => seleziono keepX = 200 perchè sintomo di overfitting

### Valutazioni performance modello con calcolo del p-value ###
perm_trascr3D <- permutation_test_singola_omica(
  X = t(trascr_filtered3D_scaled), Y = Y,
  keepX_vec = 200,
  ncomp = 2, dist = "mahalanobis.dist",
  output_dir = "Report_Trascr3D", nome_dataset = "Trascr3D"
)
# Comp1: p-BER = 0.2228  p-MCC = 0.2068 => test NON significativo
perm_trascr3D$tabella

# Recupero dei risultati del modello e delle performance
res_trascr3D_finali <- tune_and_run_dual_criterion(
  X = t(trascr_filtered3D_scaled), Y = Y,
  list.keepX = 200, 
  ncomp = 2, 
  dist = "mahalanobis.dist",
  nome_omica = "Trascrittomica 3D (keepX = 200)"
)

# Visualizza i geni scelti
geni_per_GO_3D <- res_trascr3D_finali$biomarkers_mcc
head(geni_per_GO_3D)
cat("Totale geni estratti per la GOEA:", length(geni_per_GO_3D), "\n") # 400
writeLines(geni_per_GO_3D, "geni_per_GO_3D.txt")
## comp1: BER: 0.4250  MCC: 0.1936  ##


##############
##  Prot3D  ##
##############
res_prot3D <- tune_and_run_dual_criterion(
  X = t(prot3D_scaled),
  Y = Y,
  list.keepX = c(5, 15, 25, 30, 50, 80, 140, 146, 200),
  ncomp = 2, dist = "mahalanobis.dist",
  max_extra_FN = 1, max_FN = 2, max_extra_FP = NULL,
  nome_omica = "Proteomica 3D"
)
# COMMENTO:
# TP = 2, TN = 9, FP = 1, FN = 2 fino a keepX = 15
# TP = 3, TN = 9, FP = 1, FN = 1 da keepX = 25 fino a keepX = 146 => keepX = 146
# TP = 2, TN = 9, FP = 1, FN = 2 da keepX = 147  


### Valutazioni performance modello con calcolo del p-value ###
perm_prot3D <- permutation_test_singola_omica(
  X = t(prot3D_scaled), Y = Y,
  keepX_vec = 146,
  ncomp = 2, dist = "mahalanobis.dist",
  output_dir = "Report_Prot3D", nome_dataset = "Prot3D"
)
# Comp1: p-BER = 0.013  p-MCC = 0.013 => test NON significativo
perm_prot3D$tabella

# Recupero dei risultati del modello e delle performance
res_prot3D_finali <- tune_and_run_dual_criterion(
  X = t(prot3D_scaled), Y = Y,
  list.keepX = 146, 
  ncomp = 2, 
  dist = "mahalanobis.dist",
  nome_omica = "Proteomica 3D (keepX = 146)"
)

# Visualizza i geni scelti
prot_per_GO_3D <- res_prot3D_finali$biomarkers_mcc
head(prot_per_GO_3D)
cat("Totale proteine estratte per la GOEA:", length(prot_per_GO_3D), "\n") # 291
writeLines(prot_per_GO_3D, "prot_per_GO_3D.txt")
## comp1: BER: 0.175, MCC: 0.65  ##

# Estrazione dei valori registrati per ogni componente in termini di metriche, valori della matrice di
# confusione FP, FN e nomi dei pazienti mal classificati
estrai_riepilogo_singola_omica(res_trascr2D_finali, "Trascr2D")
estrai_riepilogo_singola_omica(res_prot2D_finali,   "Prot2D")
estrai_riepilogo_singola_omica(res_trascr3D_finali, "Trascr3D")
estrai_riepilogo_singola_omica(res_prot3D_finali,   "Prot3D")


#################################################
##### 2^ STUDIO: esclusione dell'outlier RR #####
#################################################
## Questa parte è ben visibile consultando i risultati che derivano dalla composizione della tabella di confusione  ##
## TP, FP, TN, FN. Quindi inserire i dovuti commenti sul fatto che RR risulti outlier per solo un blocco omico      ##  
## mentre per gli altri NO quando si fa analisi supervisionata, comportamento diverso invece per NON supervisionata ##
## dove sembra sempre essere outlier.                                                                               ##

# Tema personalizzato: bordo chiuso su tutti i lati, titoli in grassetto e più grandi
tema_tesi_box <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    strip.text = element_text(size = 14, face = "bold", color = "black"), # Titoli interni dei riquadri
    strip.background = element_rect(fill = "gray92", color = "black", linewidth = 1),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1), # Bordo completo chiuso (inclusa barra destra)
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 12, color = "black", face = "bold"),
    legend.title = element_text(size = 13, face = "bold"),
    legend.text = element_text(size = 12, face = "bold")
  )

# Elenco dei campioni misclassificati (FN + FP) per ciascun blocco
misclass_trascr2D <- c("C1A_2D", "C2B_2D", "C3B_2D", "RR_2D")
misclass_prot2D   <- c("C1A_2D")
misclass_trascr3D <- c("C1A_3D", "C3B_3D", "RR_3D", "SF_3D")
misclass_prot3D   <- c("C1A_3D", "GAU_3D")

# Funzione per aggiungere le etichette del testo solo sui campioni misclassificati
aggiungi_etichette_selezionate <- function(p_graph, etichette_vettore, lista_misclass) {
  df_text <- data.frame(
    x = p_graph$data$x,
    y = p_graph$data$y,
    label = etichette_vettore
  )
  df_text_filt <- df_text[df_text$label %in% lista_misclass, ]
  
  p_graph + 
    geom_text_repel(
      data = df_text_filt,
      aes(x = x, y = y, label = label),
      size = 5,
      fontface = "bold",
      color = "darkred",
      bg.color = "white",
      bg.r = 0.15,
      max.overlaps = Inf,
      box.padding = 0.5
    )
}

# Estrazione plotIndiv con forme piene (pch = c(16, 17): cerchio pieno e triangolo pieno)
p1_base <- plotIndiv(res_trascr2D_finali$model_mcc, ind.names = FALSE, legend = TRUE, style = "ggplot2", 
                     title = "Trascrittomica 2D", point.size = 4.5, col.per.group = colori_mixomics, pch = c(16, 17))$graph + tema_tesi_box

p2_base <- plotIndiv(res_prot2D_finali$model_mcc,   ind.names = FALSE, legend = TRUE, style = "ggplot2", 
                     title = "Proteomica 2D",     point.size = 4.5, col.per.group = colori_mixomics, pch = c(16, 17))$graph + tema_tesi_box

p3_base <- plotIndiv(res_trascr3D_finali$model_mcc, ind.names = FALSE, legend = TRUE, style = "ggplot2", 
                     title = "Trascrittomica 3D", point.size = 4.5, col.per.group = colori_mixomics, pch = c(16, 17))$graph + tema_tesi_box

p4_base <- plotIndiv(res_prot3D_finali$model_mcc,   ind.names = FALSE, legend = TRUE, style = "ggplot2", 
                     title = "Proteomica 3D",     point.size = 4.5, col.per.group = colori_mixomics, pch = c(16, 17))$graph + tema_tesi_box

# Applicazione delle etichette del testo
p1_mod <- aggiungi_etichette_selezionate(p1_base, rownames(X2D$transcriptomics), misclass_trascr2D)
p2_mod <- aggiungi_etichette_selezionate(p2_base, rownames(X2D$transcriptomics), misclass_prot2D)
p3_mod <- aggiungi_etichette_selezionate(p3_base, rownames(X3D$transcriptomics), misclass_trascr3D)
p4_mod <- aggiungi_etichette_selezionate(p4_base, rownames(X3D$transcriptomics), misclass_prot3D)

# Composizione del grafico 2x2 con Patchwork
plot_splsda_tot <- (p1_mod | p2_mod) / (p3_mod | p4_mod) +
  plot_annotation(
    title = "Separazione sPLS-DA con Evidenza dei Campioni Misclassificati",
    tag_levels = 'A'
  ) & theme(
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    legend.position = "bottom",
    legend.key.size = unit(1.2, "cm")
  )

# Visualizzazione e salvataggio
print(plot_splsda_tot)
ggsave("Figure_sPLSDA_All_Omics_Misclassified_Filled.png", plot = plot_splsda_tot, width = 13, height = 11, dpi = 300)

####################################################################################
## 2) CERCHI DI CORRELAZIONE (plotVar Aggregato 2x2)
####################################################################################
# Tema standardizzato per i cerchi di correlazione
tema_plotvar <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(size = 13, face = "bold"),
    axis.text = element_text(size = 12, color = "black", face = "bold")
  )

# Funzione corretta per l'estrazione di plotVar
estrai_plotvar_sicuro <- function(risultato_finali, titolo) {
  modello <- risultato_finali$model_mcc
  
  # Estraiamo il numero di componenti in modo sicuro
  n_comp <- if (!is.null(modello$ncomp)) max(modello$ncomp) else 0
  
  if (!is.null(modello) && n_comp >= 2) {
    res_pv <- tryCatch({
      mixOmics::plotVar(modello, style = "ggplot2", title = titolo, var.names = FALSE)
    }, error = function(e) NULL)
    
    p <- NULL
    if (inherits(res_pv, "ggplot")) {
      p <- res_pv
    } else if (is.list(res_pv) && !is.null(res_pv$graph) && inherits(res_pv$graph, "ggplot")) {
      p <- res_pv$graph
    }
    
    if (!is.null(p)) {
      return(p + tema_plotvar)
    }
  }
  
  # Segnaposto se il modello ha davvero solo 1 componente
  p_dummy <- ggplot() + 
    annotate("text", x = 0.5, y = 0.5, 
             label = "Modello 1D\n(Cerchio 2D non applicabile)", 
             size = 5, fontface = "bold", color = "gray30") +
    theme_bw(base_size = 14) + 
    ggtitle(titolo) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 15),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
  
  return(p_dummy)
}


####################################################################################
## 3) DISTRIBUZIONE DEGLI SCORE SULLA COMPONENTE 1 (Boxplot Puliti)
####################################################################################
crea_boxplot_pulito <- function(modello, Y_fac, titolo) {
  df <- data.frame(
    Score = modello$variates$X[, 1],
    Condizione = Y_fac
  )
  
  ggplot(df, aes(x = Condizione, y = Score, fill = Condizione)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA, size = 0.9, width = 0.5) + 
    geom_jitter(width = 0.12, size = 3, alpha = 0.9, color = "black") +
    theme_classic(base_size = 13) + 
    ggtitle(titolo) +
    scale_fill_manual(values = colori_mixomics) + 
    labs(x = "Gruppo Clinico", y = "Score Comp 1") +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
      axis.text = element_text(size = 12, face = "bold", color = "black"),
      axis.title = element_text(size = 13, face = "bold")
    )
}

bp1 <- crea_boxplot_pulito(res_trascr2D_finali$model_mcc, Y, "Trascrittomica 2D")
bp2 <- crea_boxplot_pulito(res_prot2D_finali$model_mcc,   Y, "Proteomica 2D")
bp3 <- crea_boxplot_pulito(res_trascr3D_finali$model_mcc, Y, "Trascrittomica 3D")
bp4 <- crea_boxplot_pulito(res_prot3D_finali$model_mcc,   Y, "Proteomica 3D")

boxplot_grid <- (bp1 | bp2) / (bp3 | bp4) + 
  plot_annotation(
    title = "Distribuzione degli Score sPLS-DA sulla Componente 1 per tutte le Omiche", 
    tag_levels = 'A'
  ) & theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

print(boxplot_grid)
ggsave("Figure_Boxplot_Scores_Comp1.png", plot = boxplot_grid, width = 11, height = 9, dpi = 300)


####################################################################################
## 4) FEATURE LOADINGS SULLA COMPONENTE 1 (PlotLoadings con colori_mixomics)
####################################################################################
crea_plot_loadings_gg <- function(modello, titolo, top_n = 20) {
  loadings_df <- data.frame(
    Feature = rownames(modello$loadings$X),
    Weight = modello$loadings$X[, 1]
  )
  loadings_df <- loadings_df[loadings_df$Weight != 0, ]
  loadings_df <- loadings_df[order(abs(loadings_df$Weight), decreasing = TRUE), ]
  
  if (nrow(loadings_df) > top_n) {
    loadings_df <- loadings_df[1:top_n, ]
  }
  
  loadings_df$Feature <- factor(loadings_df$Feature, levels = loadings_df$Feature[order(loadings_df$Weight)])
  loadings_df$Gruppo <- ifelse(loadings_df$Weight > 0, levels(Y)[1], levels(Y)[2])
  
  ggplot(loadings_df, aes(x = Weight, y = Feature, fill = Gruppo)) +
    geom_col(width = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = colori_mixomics) + 
    theme_minimal(base_size = 12) +
    labs(title = titolo, x = "Peso della Feature (Loading Comp 1)", y = "") +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      axis.text.y = element_text(size = 10, face = "bold"),
      axis.text.x = element_text(size = 11),
      legend.position = "none"
    )
}

pl1 <- crea_plot_loadings_gg(res_trascr2D_finali$model_mcc, "Trascrittomica 2D")
pl2 <- crea_plot_loadings_gg(res_prot2D_finali$model_mcc,   "Proteomica 2D")
pl3 <- crea_plot_loadings_gg(res_trascr3D_finali$model_mcc, "Trascrittomica 3D")
pl4 <- crea_plot_loadings_gg(res_prot3D_finali$model_mcc,   "Proteomica 3D")

loadings_grid <- (pl1 | pl2) / (pl3 | pl4) +
  plot_annotation(
    title = "Top Feature Discriminanti sulla Componente 1 (Loading Plots)",
    tag_levels = 'A'
  ) & theme(
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "bottom"
  )

print(loadings_grid)
ggsave("Figure_PlotLoadings_Comp1_ggplot.png", plot = loadings_grid, width = 14, height = 11, dpi = 300)

####################################################################################
## 5) PLOT CURVE ERRORE BER e MCC
####################################################################################
g1 <- plot_tuning_curves(res_trascr2D$tabelle[[1]], "Trascrittomica 2D")
g2 <- plot_tuning_curves(res_prot2D$tabelle[[1]],   "Proteomica 2D")
g3 <- plot_tuning_curves(res_trascr3D$tabelle[[1]], "Trascrittomica 3D")
g4 <- plot_tuning_curves(res_prot3D$tabelle[[1]],   "Proteomica 3D")

# Stampa a schermo
print(g1)
print(g2)
print(g3)
print(g4)

# Salvataggio in alta definizione
ggsave("metriche_trascr2D.png", plot = g1, width = 12, height = 5, dpi = 300)
ggsave("metriche_prot2D.png",   plot = g2, width = 12, height = 5, dpi = 300)
ggsave("metriche_trascr3D.png", plot = g3, width = 12, height = 5, dpi = 300)
ggsave("metriche_prot3D.png",   plot = g4, width = 12, height = 5, dpi = 300)


# Ripristino impostazioni grafiche standard
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2) + 0.1, oma = c(0, 0, 0, 0))
dev.off()



##########################################
###      Integrazione multi-omica      ###
##########################################
## OBIETTIVO: migliorare la classificazione e identificare le correlazioni fra omiche
# (0) Confronti PLS a coppie
# (1) Design
# (2) Tuning del numero di componenti
# (3) Tuning del numero di features
# (4) Modello finale
# (5) Plot

#### (0) Confronti pls a coppie -> per capire SE le omiche sono correlate!

# Selezione arbitraria del numero di features da tenere
list.keep_trasc = c(25, 25) 
list.keep_prot = c(25, 25)
pls_2D <- spls(X2D$transcriptomics, X2D$proteomics, keepX = list.keep_trasc, keepY = list.keep_prot, ncomp = 2)
pls_3D <- spls(X3D$transcriptomics, X3D$proteomics, keepX = list.keep_trasc, keepY = list.keep_prot, ncomp = 2)
pls_trascr <- spls(X2D$transcriptomics, X3D$transcriptomics, keepX = list.keep_trasc, keepY = list.keep_trasc, ncomp = 2)
pls_prot <- spls(X2D$proteomics, X3D$proteomics, keepX = list.keep_prot, keepY = list.keep_prot, ncomp = 2)

par(mfrow = c(2, 2), mar = c(2.5, 2.5, 2.5, 1))
legenda_pers <- function(model, title, leg_text, col_vec) {
  plotVar(model, title = title, var.names = FALSE, style = 'graphics', 
          pch = c(16, 17), col = col_vec, cex = c(1.5, 1.5))
  
  # Aggiunta legenda personalizzata
  legend("topright", legend = leg_text, col = col_vec, pch = c(16, 17), 
         cex = 0.6,          
         pt.cex = 0.8,       
         bty = "n",          
         inset = c(0.02, 0.02))
}

legenda_pers(pls_2D, "Correlazione Trascr2D e Prot2D", c("Trascr2D", "Prot2D"), c('orange', 'violet'))
legenda_pers(pls_3D, "Correlazione Trascr3D e Prot3D", c("Trascr3D", "Prot3D"), c('red', 'blue'))
legenda_pers(pls_trascr, "Correlazione Trasc2D e Trasc3D", c("Trasc2D", "Trasc3D"), c('orange', 'red'))
legenda_pers(pls_prot, "Correlazione Prot2D e Prot3D", c("Prot2D", "Prot3D"), c('violet', 'blue'))
par(mfrow = c(1,1), mar = c(5, 4, 4, 2))
## OSS: l'integrazione biologica è giustificata perchè le variabili si sovrappongono
# nello spazio delle componenti. 

# Estrazione della correlazione fra le omiche/dimensioni, risultati tra le comp1
cor(pls_2D$variates$X, pls_2D$variates$Y)         # rho = 0.98344299
cor(pls_3D$variates$X, pls_3D$variates$Y)         # rho = 0.99151946
cor(pls_trascr$variates$X, pls_trascr$variates$Y) # rho = 0.9890834 
cor(pls_prot$variates$X, pls_prot$variates$Y)     # rho = 0.98891765 

#### (1) Design 
# a causa dell'elevata correlazione si deve imporre un valore alto nella matrice 
# di design, MA poichè valori ELEVATI riducono le capacità di previsione
# => si impone un valore di 0.1/0.2 per dare priorità alla DISCRIMINAZIONE

######################################
### Trascrittomica e proteomica 2D ###
######################################
cat("\n\n>>> AVVIO ANALISI MODELLO 2D <<<\n")
# ============================================================
# FASE 1 - Esplorazione iniziale
# ============================================================
risultati_fase1_2D <- grid_search_splsda_advanced(
  X             = X2D,
  Y             = Y,
  ncomp         = 2,
  design_values = c(0.01, 0.02),
  valori_trascr = c(14, 69, 99, 150),
  valori_prot   = c(30, 97, 150),
  extra_FN = 0    # accetta al massimo 1 FN in più del meglio trovato
)

# ============================================================
# FASE 2 - Raffinamento attorno al best di Fase 1
# ============================================================
parametri_fase2 <- prepare_param_second_phase_advanced(risultati_fase1_2D)

risultati_fase2_2D <- grid_search_splsda_advanced(
  X             = X2D,
  Y             = Y,
  ncomp       = 2,
  design_values = parametri_fase2$nuovi_design,
  valori_trascr = parametri_fase2$nuovi_valori_trascr,
  valori_prot   = parametri_fase2$nuovi_valori_prot,
  max_errori_positivi = 2 
)

# Parametri finali selezionati 
best_design2D    <- risultati_fase2_2D$best_design
best_keepX2D     <- risultati_fase2_2D$best_keepX
best_dist2D      <- risultati_fase2_2D$best_dist
ncomp_ottimale2D <- risultati_fase2_2D$best_ncomp

# ============================================================
# TEST DI PERMUTAZIONE - null distribution più stabile
# ============================================================
perm_test_res <- permutation_test_multi_omica(
  X = X2D, Y = Y,
  design = best_design2D, keepX = best_keepX2D,
  ncomp = ncomp_ottimale2D, dist = best_dist2D,
  output_dir = "Report_2D", nome_dataset = "2D"
)

# ============================================================
# MODELLO FINALE + REPORT COMPLETO
# ============================================================
risultati_2D <- run_splsda_pipeline(
  X      = X2D,
  Y      = Y,
  design = best_design2D,
  keepX  = best_keepX2D,
  ncomp  = ncomp_ottimale2D,
  dist   = best_dist2D
)

print(risultati_2D$metrics) # BER=0.375, MCC=0.4385

generate_splsda_report_optimized(
  risultati_pipeline = risultati_2D, 
  nome_dataset       = "2D", 
  Y                  = Y,
  colori_gruppi      = colori_mixomics,
  labels_blocchi     = c("transcriptomics" = "Trascrittomica 2D", "proteomics" = "Proteomica 2D")
)

ispeziona_errori(risultati_2D, Y[1:14], dist = best_dist2D, ncomp = 1)
ispeziona_errori(risultati_2D, Y[1:14], dist = best_dist2D, ncomp = 2)
#generate_poster_plots(risultati_pipeline = risultati_2D, nome_dataset = "2D")

######################################
### Trascrittomica e proteomica 3D ###
######################################
cat("\n\n>>> AVVIO ANALISI MODELLO 3D <<<\n")
# ============================================================
# FASE 1 - Esplorazione iniziale
# ============================================================
risultati_fase1_3D <- grid_search_splsda_advanced(
  X             = X3D,
  Y             = Y,
  ncomp         = 2,
  #design_values = c(0.01, 0.05),
  design_values = c(0.3, 0.4),    # questi valori alti permettono di far emergere il segnale FORTE di Prot3D
  valori_trascr = c(15, 120, 200),
  valori_prot   = c(15, 25, 146, 160),
  max_errori_positivi = 2   # al massimo 2 sani su 4 classificati male
)

# ============================================================
# FASE 2 - Raffinamento attorno al best di Fase 1
# ============================================================
parametri_fase2 <- prepare_param_second_phase_advanced(risultati_fase1_3D)

risultati_fase2_3D <- grid_search_splsda_advanced(
  X             = X3D,
  Y             = Y,
  ncomp         = 2,
  design_values = parametri_fase2$nuovi_design,
  valori_trascr = parametri_fase2$nuovi_valori_trascr,
  valori_prot   = parametri_fase2$nuovi_valori_prot,
  max_errori_positivi = 2 
)

# Parametri finali selezionati 
best_design3D    <- risultati_fase2_3D$best_design
best_keepX3D     <- risultati_fase2_3D$best_keepX
best_dist3D      <- risultati_fase2_3D$best_dist
ncomp_ottimale3D <- risultati_fase2_3D$best_ncomp

# ============================================================
# TEST DI PERMUTAZIONE - null distribution più stabile
# ============================================================
perm_test_res <- permutation_test_multi_omica(
  X = X3D, Y = Y,
  design = best_design3D, keepX = best_keepX3D,
  ncomp = ncomp_ottimale3D, dist = best_dist3D,
  output_dir = "Report_3D", nome_dataset = "3D"
)


# ============================================================
# MODELLO FINALE + REPORT COMPLETO
# ============================================================
risultati_3D <- run_splsda_pipeline(
  X      = X3D,
  Y      = Y,
  design = best_design3D,
  keepX  = best_keepX3D,
  ncomp  = ncomp_ottimale3D,
  dist   = best_dist3D
)

print(risultati_3D$metrics) # BER=0.3500  MCC=0.3000

ispeziona_errori(risultati_3D, Y[1:14], dist = best_dist3D, ncomp = 1)
ispeziona_errori(risultati_3D, Y[1:14], dist = best_dist3D, ncomp = 2)

generate_splsda_report_optimized(
  risultati_pipeline = risultati_3D, 
  nome_dataset       = "3D", 
  Y                  = Y,
  colori_gruppi      = colori_mixomics,
  labels_blocchi     = c("transcriptomics" = "Trascrittomica 3D", "proteomics" = "Proteomica 3D")
)


##############################
### Trascrittomica 2D e 3D
##############################
cat("\n\n>>> AVVIO ANALISI MODELLO con Trascrittomica <<<\n")

# Funzione che rimuove i suffissi 2D e 3D dai nomi dei pazienti
clean_rownames <- function(mat) {
  rownames(mat) <- gsub("_2D|_3D", "", rownames(mat))
  return(mat)
}

# Funzione per uniformare i nomi per trasformare tutti i '-' in '_'
uniform_rownames <- function(mat) {
  rownames(mat) <- gsub("-", "_", rownames(mat))
  return(mat)
}

# Pulizia dei nomi
trascr_filtered2D_scaled_trans <- clean_rownames(t(trascr_filtered2D_scaled))
trascr_filtered3D_scaled_trans <- clean_rownames(t(trascr_filtered3D_scaled))
trascr_filtered2D_scaled_trans <- uniform_rownames(trascr_filtered2D_scaled_trans)
trascr_filtered3D_scaled_trans <- uniform_rownames(trascr_filtered3D_scaled_trans)

X_Trascr <- list(transcriptomics = trascr_filtered2D_scaled_trans,
            proteomics     = trascr_filtered3D_scaled_trans)

# ============================================================
# FASE 1 - Esplorazione iniziale
# ============================================================
risultati_fase1_Trascr <- grid_search_splsda_advanced(
  X             = X_Trascr,
  Y             = Y,
  ncomp         = 2,
  design_values = c(0.3, 0.4),
  valori_trascr = c(14, 69, 99, 150),
  valori_prot   = c(30, 97, 150),
  extra_FN      = 0    # accetta al massimo 1 FN in più del meglio trovato
)

# ============================================================
# FASE 2 - Raffinamento attorno al best di Fase 1
# ============================================================
parametri_fase2 <- prepare_param_second_phase_advanced(risultati_fase1_Trascr)

risultati_fase2_Trascr <- grid_search_splsda_advanced(
  X             = X_Trascr,
  Y             = Y,
  ncomp         = 2,
  design_values = parametri_fase2$nuovi_design,
  valori_trascr = parametri_fase2$nuovi_valori_trascr,
  valori_prot   = parametri_fase2$nuovi_valori_prot,
  max_errori_positivi = 2 
)

# Parametri finali selezionati 
best_design_Trascr    <- risultati_fase2_Trascr$best_design
best_keepX_Trascr     <- risultati_fase2_Trascr$best_keepX
best_dist_Trascr      <- risultati_fase2_Trascr$best_dist
ncomp_ottimale_Trascr <- risultati_fase2_Trascr$best_ncomp

# ============================================================
# TEST DI PERMUTAZIONE - null distribution più stabile
# ============================================================
perm_test_res <- permutation_test_multi_omica(
  X = X_Trascr, Y = Y,
  design = best_design_Trascr, keepX = best_keepX_Trascr,
  ncomp = ncomp_ottimale_Trascr, dist = best_dist_Trascr,
  output_dir = "Report_Trascr", nome_dataset = "Trascr"
)

# ============================================================
# MODELLO FINALE + REPORT COMPLETO
# ============================================================
risultati_Trascr <- run_splsda_pipeline(
  X      = X_Trascr,
  Y      = Y,
  design = best_design_Trascr,
  keepX  = best_keepX_Trascr,
  ncomp  = ncomp_ottimale_Trascr,
  dist   = best_dist_Trascr
)

print(risultati_Trascr$metrics) # BER=0.2500  MCC=0.6455

generate_splsda_report_optimized(
  risultati_pipeline = risultati_Trascr, 
  nome_dataset       = "Trascr", 
  Y                  = Y,
  colori_gruppi      = colori_mixomics,
  labels_blocchi     = c("transcriptomics" = "Trascrittomica 2D", "proteomics" = "Trascrittomica 3D")
)

ispeziona_errori(risultati_Trascr, Y[1:14], dist = best_dist_Trascr, ncomp = 1)
ispeziona_errori(risultati_Trascr, Y[1:14], dist = best_dist_Trascr, ncomp = 2)
#generate_poster_plots(risultati_pipeline = risultati_Trascr, nome_dataset = "Trascr")


##############################
### Proteomica 2D e 3D
##############################

cat("\n\n>>> AVVIO ANALISI MODELLO con Proteomica <<<\n")

# Pulizia dei nomi
prot_scaled2D_trans <- clean_rownames(t(prot2D_scaled))
prot_scaled3D_trans <- clean_rownames(t(prot3D_scaled))
prot_scaled2D_trans <- uniform_rownames(prot_scaled2D_trans)
prot_scaled3D_trans <- uniform_rownames(prot_scaled3D_trans)

X_Prot <- list(transcriptomics = prot_scaled2D_trans,
                 proteomics     = prot_scaled3D_trans)

# ============================================================
# FASE 1 - Esplorazione iniziale
# ============================================================
risultati_fase1_Prot <- grid_search_splsda_advanced(
  X             = X_Prot,
  Y             = Y,
  ncomp         = 2,
  design_values = c(0.1, 0.2),
  valori_trascr = c(30, 97, 150, 200, 300),
  valori_prot   = c(15, 25, 146, 300),
  extra_FN      = 0    # accetta al massimo 1 FN in più del meglio trovato
)

# ============================================================
# FASE 2 - Raffinamento attorno al best di Fase 1
# ============================================================
parametri_fase2 <- prepare_param_second_phase_advanced(risultati_fase1_Prot)

risultati_fase2_Prot <- grid_search_splsda_advanced(
  X             = X_Prot,
  Y             = Y,
  ncomp         = 2,
  design_values = parametri_fase2$nuovi_design,
  valori_trascr = parametri_fase2$nuovi_valori_trascr,
  valori_prot   = parametri_fase2$nuovi_valori_prot,
  max_errori_positivi = 2 
)

# Parametri finali selezionati 
best_design_Prot    <- risultati_fase2_Prot$best_design
best_keepX_Prot     <- risultati_fase2_Prot$best_keepX
best_dist_Prot      <- risultati_fase2_Prot$best_dist
ncomp_ottimale_Prot <- risultati_fase2_Prot$best_ncomp

# ============================================================
# TEST DI PERMUTAZIONE - null distribution più stabile
# ============================================================
perm_test_res <- permutation_test_multi_omica(
  X = X_Prot, Y = Y,
  design = best_design_Prot, keepX = best_keepX_Prot,
  ncomp = ncomp_ottimale_Prot, dist = best_dist_Prot,
  output_dir = "Report_Prot", nome_dataset = "Prot"
)

# ============================================================
# MODELLO FINALE + REPORT COMPLETO
# ============================================================
risultati_Prot <- run_splsda_pipeline(
  X      = X_Prot,
  Y      = Y,
  design = best_design_Prot,
  keepX  = best_keepX_Prot,
  ncomp  = ncomp_ottimale_Prot,
  dist   = best_dist_Prot
)

print(risultati_Prot$metrics) # BER=0.125  MCC=0.8257

generate_splsda_report_optimized(
  risultati_pipeline = risultati_Prot, 
  nome_dataset       = "Prot", 
  Y                  = Y,
  colori_gruppi      = colori_mixomics,
  labels_blocchi     = c("transcriptomics" = "Proteomica 2D", "proteomics" = "Proteomica 3D")
)

ispeziona_errori(risultati_Prot, Y[1:14], dist = best_dist_Prot, ncomp = 1)
ispeziona_errori(risultati_Prot, Y[1:14], dist = best_dist_Prot, ncomp = 2)
#generate_poster_plots(risultati_pipeline = risultati_Prot, nome_dataset = "Prot")


######################################
### Confronto biomarcatori 2D e 3D ###
######################################

### Estrazioni dei biomarcatori dall'analisi delle singole omiche ###
geni_per_GO_2D # 28
geni_per_GO_3D # 400
prot_per_GO_2D # 194
prot_per_GO_3D # 291
cat("Trascrittomica -> 2D:", length(geni_per_GO_2D), "| 3D:", length(geni_per_GO_3D), "| Comuni:", length(intersect(geni_per_GO_2D, geni_per_GO_3D)), "\n")
cat("Proteomica -> 2D:", length(prot_per_GO_2D), "| 3D:", length(prot_per_GO_3D), "| Comuni:", length(intersect(prot_per_GO_2D, prot_per_GO_3D)), "\n")


### Estrazione dei veri biomarcatori dall'integrazione multi-omica ###

## Studio di integrazione per DIMENSIONE (2D e 3D separato)
geni_2D <- risultati_2D$biomarkers$transcriptomics # 28 (vs prima 5)
geni_3D <- risultati_3D$biomarkers$transcriptomics # 412 (vs 44 prima)

prot_2D <- risultati_2D$biomarkers$proteomics      # 68 (vs 5 prima)
prot_3D <- risultati_3D$biomarkers$proteomics      # 328 (vs 46 prima)
 
## Studio di integrazione per OMICA (trascrittomica e proteomica separate)
geni_da_Trascr_2D <- risultati_Trascr$biomarkers$transcriptomics # 34 (vs 5 prima)
geni_da_Trascr_3D <- risultati_Trascr$biomarkers$proteomics      # 64 (vs 6 prima)

# Estrazione mirata sulla sola Comp 2 per Trascrittomica (poichè p-value significativo solo lì)
get_biomarkers_comp <- function(model, block_idx, comp_num) {
  # Estrae la matrice dei loadings per il blocco specificato
  loadings_mat <- model$loadings[[block_idx]]
  
  # Isola la colonna della componente desiderata
  comp_loadings <- loadings_mat[, comp_num]
  
  # Seleziona solo i nomi delle feature con valore diverso da zero
  features_selected <- names(comp_loadings[comp_loadings != 0])
  return(features_selected)
}


# Estrazione per il Blocco 1 (es. 2D) e Blocco 2 (es. 3D) sulla sola Componente 2
geni_Trascr_comp2_2D <- get_biomarkers_comp(risultati_Trascr$model, block_idx = 1, comp_num = 2)
geni_Trascr_comp2_3D <- get_biomarkers_comp(risultati_Trascr$model, block_idx = 2, comp_num = 2)

# Controllo lunghezze
cat("Geni Trascr 2D (Comp 2):", length(geni_Trascr_comp2_2D), "\n")
cat("Geni Trascr 3D (Comp 2):", length(geni_Trascr_comp2_3D), "\n")

prot_da_Prot_2D <- risultati_Prot$biomarkers$transcriptomics     # 312 (vs 48 prima)
prot_da_Prot_3D <- risultati_Prot$biomarkers$proteomics          # 58 (vs 92 prima)

# Salva le liste pulite per clusterProfiler / DAVID / Panther
writeLines(geni_da_Trascr_2D, "GO_geni_Trascr_2D.txt")
writeLines(prot_da_Prot_2D, "GO_prot_Prot_2D.txt")


# Calcolo dei geni comuni
geni_comuni <- intersect(geni_2D, geni_3D) 
prot_comuni <- intersect(prot_2D, prot_3D) 
geni_da_Trascr_comuni <- intersect(geni_da_Trascr_2D, geni_da_Trascr_3D) 
prot_da_Prot_comuni <- intersect(prot_da_Prot_2D, prot_da_Prot_3D) 

geni_2D_da_studi_diversi <- intersect(geni_2D, geni_da_Trascr_2D) 
geni_3D_da_studi_diversi <- intersect(geni_3D, geni_da_Trascr_3D)
prot_2D_da_studi_diversi <- intersect(prot_2D, prot_da_Prot_2D)
prot_3D_da_studi_diversi <- intersect(prot_3D, prot_da_Prot_3D) 

# Confronti con studio su singole omiche
geni_2D_intersec_singolaOmica <- intersect(geni_2D, geni_per_GO_2D) 
geni_3D_intersec_singolaOmica <- intersect(geni_per_GO_3D, geni_da_Trascr_3D)
prot_2D_intersec_singolaOmica <- intersect(prot_2D, prot_per_GO_2D) 
prot_3D_intersec_singolaOmica <- intersect(prot_3D, prot_per_GO_3D) 

cat("\n--- REPORT INTERSEZIONE TRASCRITTOMICA ---\n")
cat("Geni selezionati nel 2D:", length(geni_2D), "\n") 
cat("Geni selezionati nel 3D:", length(geni_3D), "\n")  
cat("Geni in COMUNE (Core Signature):", length(geni_comuni), "\n") # 
cat("Geni in 2D selezionati sempre:", length(geni_2D_da_studi_diversi))  
cat("Geni in 3D selezionati sempre:", length(geni_3D_da_studi_diversi)) 
cat("Geni comuni usando come due blocchi Trasrc2D e Trascr3D:", length(geni_da_Trascr_comuni)) #

cat("\n--- REPORT INTERSEZIONE PROTEOMICA ---\n")
cat("Proteine selezionate nel 2D:", length(prot_2D), "\n") 
cat("Proteine selezionate nel 3D:", length(prot_3D), "\n") 
cat("Proteine in COMUNE (Core Signature):", length(prot_comuni), "\n\n") # 0
cat("Proteine in 2D selezionate sempre:", length(prot_2D_da_studi_diversi))
cat("Proteine in 3D selezionate sempre:", length(prot_3D_da_studi_diversi)) 
cat("Proteine comuni usando come due blocchi Prot2D e Prot3D:", length(prot_da_Prot_comuni)) # 


print("Anteprima Geni in comune:")
print(head(geni_comuni))
print("Anteprima Geni 2D in comune:")
print(head(geni_2D_da_studi_diversi))
print("Anteprima Geni 3D in comune:")
print(head(geni_3D_da_studi_diversi))
print("Anteprima Geni comuni usando come blocchi Trascr2D e Trascr3D:")
print(head(geni_da_Trascr_comuni))

print("Anteprima Proteine in comune:")
print(head(prot_comuni))
print("Anteprima Proteine 2D in comune:")
print(head(prot_2D_da_studi_diversi))
print("Anteprima Proteine 3D in comune:")
print(head(prot_3D_da_studi_diversi))
print("Anteprima Proteine comuni usando come blocchi Prot2D e Prot3D:")
print(head(prot_da_Prot_comuni))

# 3. Visualizzazione tramite Diagrammi di Venn
lista_venn_2D <- list(
  "Prot 2D (vs 1)" = prot_da_Prot_2D, 
  "Trascr 2D (vs 4)" = geni_da_Trascr_2D
)

lista_venn_3D <- list(
  "Prot 3D (vs 2)" = prot_da_Prot_3D, 
  "Trascr 3D (vs 3)" = geni_da_Trascr_3D
)


# Plot 2D
plot_2D <- ggVennDiagram(lista_venn_2D, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "skyblue") +
  labs(title = "Stabilità Molecolare: Modello 2D") +
  theme(legend.position = "none")

# Plot 3D
plot_3D <- ggVennDiagram(lista_venn_3D, label_alpha = 0) +
  scale_fill_gradient(low = "white", high = "orange") +
  labs(title = "Stabilità Molecolare: Modello 3D") +
  theme(legend.position = "none")

dev.off() 
dev.new()
plot_2D + plot_3D

all_items <- unique(c(geni_2D, geni_3D, prot_2D, prot_3D, geni_da_Trascr_2D, geni_da_Trascr_3D, prot_da_Prot_2D, prot_da_Prot_3D))

df_upset <- data.frame(
  ID = all_items,
  Geni_2D = as.integer(all_items %in% geni_2D),
  Geni_3D = as.integer(all_items %in% geni_3D),
  Prot_2D = as.integer(all_items %in% prot_2D),
  Prot_3D = as.integer(all_items %in% prot_3D),
  Geni_da_Trascr_2D = as.integer(all_items %in% geni_da_Trascr_2D),
  Geni_da_Trascr_3D = as.integer(all_items %in% geni_da_Trascr_3D),
  Prot_da_Prot_2D = as.integer(all_items %in% prot_da_Prot_2D),
  Prot_da_Prot_3D = as.integer(all_items %in% prot_da_Prot_3D)
)

# Creazione del grafico UpSet
upset(df_upset, 
      sets = c("Geni_2D", "Geni_3D", "Prot_2D", "Prot_3D", 
               "Geni_da_Trascr_2D", "Geni_da_Trascr_3D", "Prot_da_Prot_2D", "Prot_da_Prot_3D"), 
      text.scale = c(1.5, 1.5, 1.5, 1.5, 1.5, 1.5), 
      
      number.angles = 0, 
      point.size = 3.5, 
      line.size = 1.5,
      mainbar.y.label = "Numero di biomarcatori condivisi", 
      sets.x.label = "Totale biomarcatori per studio",
      order.by = "freq")

#####################################################################################
####          GENE ONTOLOGY - PROTEOMICA: SINGOLA OMICA + INTEGRAZIONE           #### 
####          Lista completa dei termini testati, ordinata per p-value,          ####
####          con fasce alpha = 0.05, 0.1, 0.25, plot ed export Excel            ####
#####################################################################################

## ------------------------------------------------------------------------- ##
## 0) PARAMETRI                                                              ##
## ------------------------------------------------------------------------- ##

alpha_values <- c(0.05, 0.1, 0.25)
p_max <- 0.5
ontologie = c("BP", "CC", "MF")
file_excel <- "GO_Prot_singola_e_integrazione.xlsx"

colonne_out <- c("Cluster", "ID", "Description", "GeneRatio", "BgRatio",
                 "pvalue", "p.adjust", "soglie", "geneID")


## ------------------------------------------------------------------------- ##
## 1) DEFINIZIONE ED ELABORAZIONE BACKGROUND (UNIVERSO)                      ##
## ------------------------------------------------------------------------- ##

clean_ensembl <- function(ids) {
  ids <- sub("\\..*", "", ids)
  return(unique(na.omit(ids)))
}

# 1a) Conversione SYMBOL -> ENTREZID (Proteomica) 
converti_symbol_to_entrez <- function(id_list, nome_lista = "Dataset") {
  id_clean <- unique(na.omit(id_list))
  mapped <- bitr(
    geneID   = id_clean,
    fromType = "SYMBOL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  return(unique(mapped$ENTREZID))
}

# 1b) Conversione ENSEMBL -> ENTREZID (Trascrittomica) 
converti_ensembl_to_entrez <- function(id_list, nome_lista = "Dataset") {
  ids_clean <- clean_ensembl(id_list)
  mapped <- bitr(
    geneID   = ids_clean,
    fromType = "ENSEMBL",
    toType   = "ENTREZID",
    OrgDb    = org.Hs.eg.db
  )
  return(unique(mapped$ENTREZID))
}

## ------------------------------------------------------------------------- ##
## 2) CONVERSIONE BACKGROUND E BIOMARCATORI UNIFORMATI                       ##
## ------------------------------------------------------------------------- ##

# Background
bg_prot   <- rownames(prot2D_scaled)
bg_prot_entrez   <- converti_symbol_to_entrez(bg_prot, "BG Prot")

# Proteomica (Uniformata per tutte le componenti)
prot_per_GO_2D_entrez <- converti_symbol_to_entrez(prot_per_GO_2D, "Prot 2D Singola")
prot_per_GO_3D_entrez <- converti_symbol_to_entrez(prot_per_GO_3D, "Prot 3D Singola")
prot_da_Prot_2D_entrez <- converti_symbol_to_entrez(prot_da_Prot_2D, "Prot 2D Int Prot")
prot_da_Prot_3D_entrez <- converti_symbol_to_entrez(prot_da_Prot_3D, "Prot 3D Int Prot")

## ------------------------------------------------------------------------- ##
## 3) LISTE DI GENI DEI DUE STUDI                                            ##
## ------------------------------------------------------------------------- ##

## Studio 1: singola proteomica
lista_prot_singola <- list(
  "Proteomica 2D" = prot_per_GO_2D_entrez,
  "Proteomica 3D" = prot_per_GO_3D_entrez
)

## Studio 2: integrazione, con lista unificata 2D + 3D 
lista_prot_int <- list(
  "Prot 2D+3D unificata" = unique(c(prot_da_Prot_2D_entrez, prot_da_Prot_3D_entrez))
)

studi <- list(
  "Singola Omica" = list(lista = lista_prot_singola, universo = bg_prot_entrez),
  "Integrazione" = list(lista = lista_prot_int, universo = bg_prot_entrez)
)

## ------------------------------------------------------------------------- ##
## 4) FUNZIONI                                                               ##
## ------------------------------------------------------------------------- ##

## 4a) compareCluster completo: tiene tutti i termini testati (NON applica filtri).
run_full <- function(lista, universo, ont, cutoff = p_max) {
  compareCluster(
    geneClusters  = lista,
    fun           = "enrichGO",
    universe      = universo,
    OrgDb         = org.Hs.eg.db,
    ont           = ont,
    pAdjustMethod = "BH",
    pvalueCutoff  = 1,  
    qvalueCutoff  = 1,
    readable      = FALSE     # simboli convertiti dopo, in prepara_tabella()
  )
}

## 4b) Tabella ordinata per p-value, con colonna "range"
prepara_tabella <- function(res, soglie = alpha_values, p_max = 0.5) {
 
   if (is.null(res)) {
    return(setNames(
      data.frame(matrix(ncol = length(colonne_out), nrow = 0)),
      colonne_out
    ))
  }
  
  tab <- as.data.frame(res)
  
  if (nrow(tab) == 0) {
    return(setNames(
      data.frame(matrix(ncol = length(colonne_out), nrow = 0)),
      colonne_out
    ))
  }
  tab <- tab[tab$pvalue <= p_max, , drop = FALSE]
  tab <- tab[order(tab$Cluster, tab$pvalue), , drop = FALSE]
  
  if (nrow(tab) > 0) {
    sym <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys     = unique(unlist(strsplit(tab$geneID, "/"))),
      column   = "SYMBOL", keytype = "ENTREZID", multiVals = "first"
    )
    tab$geneID <- vapply(strsplit(tab$geneID, "/"),
                         function(g) paste(sym[g], collapse = "/"), character(1))
  }
  
  tab$soglie <- as.character(cut(
    tab$p.adjust,
    breaks = c(-Inf, soglie, Inf),
    labels = c(paste0("<= ", soglie[1]),
               paste0(soglie[-length(soglie)], " - ", soglie[-1]),
               paste0("> ", max(soglie)))
  ))
  tab$Cluster <- as.character(tab$Cluster)
  tab[, colonne_out, drop = FALSE]
}

## 4c) Conteggio termini per cluster e alpha 
riepilogo_alpha <- function(res, nomi_cluster, range = alpha_values) {
  
  df <- as.data.frame(res)
  
  out <- matrix(
    0,
    nrow = length(nomi_cluster),
    ncol = length(range),
    dimnames = list(
      nomi_cluster,
      paste0("alpha_", range)
    )
  )
  
  for (j in seq_along(range)) {
    for (i in seq_along(nomi_cluster)) {
      
      out[i, j] <- sum(
        df$Cluster == nomi_cluster[i] &
          !is.na(df$p.adjust) &
          df$p.adjust <= range[j]
      )
    }
  }
  
  return(out)
}

## 4d) Stampa a terminale con una linea a fine di ogni soglie
stampa_con_soglie <- function(tab, nomi_cluster, p_adjust_max = 0.25) {
  
  for (cl in nomi_cluster) {
    
    d <- tab[
      tab$Cluster == cl & tab$p.adjust <= p_adjust_max,
      setdiff(colonne_out, "geneID"),
      drop = FALSE
    ]
    
    cat(
      "\n--- ", cl,
      " (", nrow(d),
      " termini con pvalue <= ", p_max,
      " e p.adjust <= ", p_adjust_max, ") ---\n",
      sep = ""
    )
    
    if (nrow(d) == 0) next
    
    righe <- capture.output(print(d, row.names = FALSE))
    
    cat(righe[1], "\n", sep = "")
    
    for (i in seq_len(nrow(d))) {
      
      cat(righe[i + 1], "\n", sep = "")
      
      if (i < nrow(d) && d$soglie[i] != d$soglie[i + 1]) {
        etich <- paste0(
          " fine soglie p.adjust ",
          d$soglie[i],
          " ")
        
        cat(
          etich,
          strrep("-", max(0, nchar(righe[1]) - nchar(etich))),
          "\n",
          sep = ""
        )}}}
}

## 4e) Dotplot a una data soglia alpha (con testi ingranditi per tesi)
dotplot_alpha <- function(res, alpha, titolo) {
  df <- as.data.frame(res)
  df <- df[!is.na(df$p.adjust) & df$p.adjust <= alpha, , drop = FALSE]
  
  if (nrow(df) == 0) {
    return(
      ggplot() + 
        annotate("text", x = 0.5, y = 0.5, label = "Nessun termine significativo", size = 6) +
        ggtitle(titolo) + 
        theme_bw() + 
        theme(
          axis.text   = element_blank(), 
          axis.ticks  = element_blank(), 
          axis.title  = element_blank(),
          plot.title  = element_text(size = 16, face = "bold", hjust = 0.5)
        )
    )
  }
  
  df$Description <- factor(df$Description, levels = unique(df$Description))
  
  r <- res
  r@compareClusterResult <- df
  
  p <- dotplot(r, showCategory = 10, title = titolo, font.size = 14) +
    theme_bw(base_size = 15) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      # Titolo del singolo pannello
      plot.title   = element_text(size = 16, face = "bold", hjust = 0.5, margin = margin(b = 10)),
      # Etichette e titoli degli assi X e Y
      axis.text.y  = element_text(size = 14, color = "black", face = "bold"),
      axis.text.x  = element_text(angle = 35, hjust = 1, size = 13, color = "black", face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold", margin = margin(t = 8)),
      axis.title.y = element_blank(),
      # Testi e titoli della legenda (p.adjust e GeneRatio)
      legend.title = element_text(size = 14, face = "bold"),
      legend.text  = element_text(size = 13, color = "black"),
      legend.key.size = unit(0.8, "cm")
    )
  
  return(p)
}

## 2f) Scrive un foglio con blocchi BP / CC / MF e un sotto-blocco per cluster
scrivi_foglio <- function(wb, sheet, risultati, nomi_cluster,
                          soglie = alpha_values, p_max = 0.5) {
  addWorksheet(wb, sheet)
  st_ont  <- createStyle(textDecoration = "bold", fontSize = 14,
                         fgFill = "#1F4E78", fontColour = "#FFFFFF")
  st_cl   <- createStyle(textDecoration = "bold", fgFill = "#DDEBF7")
  st_head <- createStyle(textDecoration = "bold", border = "bottom")
  st_sig  <- createStyle(fgFill = "#E2F0D9")
  st_line <- createStyle(border = "bottom", borderStyle = "medium")
  nc   <- length(colonne_out)
  riga <- 1
  
  for (ont in names(risultati)) {
    tab <- prepara_tabella(risultati[[ont]], soglie, p_max)
    
    writeData(wb, sheet, paste("Gene Ontology -", ont), startRow = riga, colNames = FALSE)
    addStyle(wb, sheet, st_ont, rows = riga, cols = 1:nc, stack = TRUE)
    riga <- riga + 2
    
    for (cl in nomi_cluster) {
      d <- tab[tab$Cluster == cl, , drop = FALSE]
      
      writeData(wb, sheet,
                paste0(cl, "  (", nrow(d), " termini con pvalue <= ", p_max, ")"),
                startRow = riga, colNames = FALSE)
      addStyle(wb, sheet, st_cl, rows = riga, cols = 1:nc, stack = TRUE)
      riga <- riga + 1
      
      if (nrow(d) == 0) { riga <- riga + 2; next }
      
      writeData(wb, sheet, d, startRow = riga)
      addStyle(wb, sheet, st_head, rows = riga, cols = 1:nc, stack = TRUE)
      
      sig <- which(d$p.adjust <= max(soglie))              # passa almeno la soglia piu' larga
      if (length(sig))
        addStyle(wb, sheet, st_sig, rows = riga + sig, cols = 1:nc,
                 gridExpand = TRUE, stack = TRUE)
      
      fine <- which(d$soglie[-nrow(d)] != d$soglie[-1])    # linea a fine di ogni soglia
      if (length(fine))
        addStyle(wb, sheet, st_line, rows = riga + fine, cols = 1:nc,
                 gridExpand = TRUE, stack = TRUE)
      
      riga <- riga + nrow(d) + 3
    }
    riga <- riga + 1
  }
  setColWidths(wb, sheet, cols = 1:nc, widths = c(22, 13, 50, 11, 11, 11, 11, 12, 60))
}

## ------------------------------------------------------------------------- ##
## 5) Calcolo funzione                                                       ##
## ------------------------------------------------------------------------- ##
lavori <- expand.grid(studio = names(studi), ont = ontologie, stringsAsFactors = FALSE)

risultati_completi <- setNames(
  lapply(names(studi), function(s) setNames(vector("list", length(ontologie)), ontologie)),
  names(studi)
)

t0 <- Sys.time()
for (k in seq_len(nrow(lavori))) {
  s <- lavori$studio[k]; o <- lavori$ont[k]
  t1 <- Sys.time()
  risultati_completi[[s]][o] <- list(run_full(studi[[s]]$lista, studi[[s]]$universo, o))
  
  passo <- as.numeric(difftime(Sys.time(), t1, units = "secs"))
  tot   <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  cat(sprintf("[%d/%d] %s - %s: %.0f s | trascorsi %.1f min | stima rimanente ~%.1f min\n",
              k, nrow(lavori), s, o, passo, tot / 60, tot / k * (nrow(lavori) - k) / 60))
}


## ------------------------------------------------------------------------- ##
## 6) Output terminale e plot                                                ##
## ------------------------------------------------------------------------- ##

for (studio in names(studi)) {
  nomi_cluster <- names(studi[[studio]]$lista)
  
  for (ont in ontologie) {
    res <- risultati_completi[[studio]][[ont]]
    
    cat("\n\n########################################################\n")
    cat("#", studio, "-", ont, "\n")
    cat("########################################################\n")
    
    cat("\nNumero di termini per soglia alpha (p.adjust):\n")
    print(riepilogo_alpha(res, nomi_cluster))
    
    tab <- prepara_tabella(res, alpha_values, p_max)
    stampa_con_soglie(tab, nomi_cluster)
    
    pp <- lapply(alpha_values, function(a)
      dotplot_alpha(res, a, paste0(studio, " ", ont, " - alpha = ", a)))
    print(pp[[1]] | pp[[2]] | pp[[3]])
  }
}


## ------------------------------------------------------------------------- ##
## 6b) GRAFICI AGGIUNTIVI GOEA: ENRICHMENT MAP E CNETPLOT                    ##
## ------------------------------------------------------------------------- ##

# --- A) CELLULAR COMPONENT (CC) - Integrazione Proteomica ---
res_int_cc <- risultati_completi[["Integrazione"]][["CC"]]

# Conversione ENTREZID -> SYMBOL e filtro alpha = 0.05
res_int_cc_clean <- setReadable(res_int_cc, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
df_cc <- as.data.frame(res_int_cc_clean)
res_int_cc_clean@compareClusterResult <- df_cc[!is.na(df_cc$p.adjust) & df_cc$p.adjust <= 0.05, , drop = FALSE]

# Calcolo matrice di somiglianza tra termini GO
res_int_cc_sim <- pairwise_termsim(res_int_cc_clean)

# 1) Enrichment Map (CC)
p_emap_cc <- emapplot(res_int_cc_sim, 
                      showCategory = 15) +
  ggtitle("Enrichment Map: Cellular Component (Proteomica Integrata, alpha = 0.05)") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

print(p_emap_cc)
ggsave("GO_Emapplot_Integrazione_CC.png", plot = p_emap_cc, width = 12, height = 9, dpi = 300)

# 2) Cnetplot (CC) - top 5 termini GO e proteine collegate
p_cnet_cc <- cnetplot(res_int_cc_clean, 
                      showCategory = 5) +
  ggtitle("Rete Funzionale: Cellular Component e Proteine Selezionate") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

print(p_cnet_cc)
ggsave("GO_Cnetplot_Integrazione_CC.png", plot = p_cnet_cc, width = 13, height = 10, dpi = 300)


# --- B) BIOLOGICAL PROCESS (BP) - Integrazione Proteomica ---
res_int_bp <- risultati_completi[["Integrazione"]][["BP"]]

# Conversione ENTREZID -> SYMBOL e filtro alpha = 0.05
res_int_bp_clean <- setReadable(res_int_bp, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
df_bp <- as.data.frame(res_int_bp_clean)
res_int_bp_clean@compareClusterResult <- df_bp[!is.na(df_bp$p.adjust) & df_bp$p.adjust <= 0.05, , drop = FALSE]

# 3) Cnetplot (BP) - sui 3 termini significativi di Biological Process
p_cnet_bp <- cnetplot(res_int_bp_clean, 
                      showCategory = 3) +
  ggtitle("Rete Funzionale: Biological Process e Proteine Selezionate") +
  theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))

print(p_cnet_bp)
ggsave("GO_Cnetplot_Integrazione_BP.png", plot = p_cnet_bp, width = 13, height = 10, dpi = 300)


## ------------------------------------------------------------------------- ##
## 7) EXPORT EXCEL: due fogli (Singola omica / Integrazione)                 ##
## ------------------------------------------------------------------------- ##
wb <- createWorkbook()
for (studio in names(studi)) {
  scrivi_foglio(wb, studio, risultati_completi[[studio]],
                names(studi[[studio]]$lista), alpha_values, p_max)
}
saveWorkbook(wb, file_excel, overwrite = TRUE)

#######################################
#### Singola omica Prot2D e Prot3D ####
#######################################
### ------------------------ ###
### "BP": Biological Process ###
### ------------------------ ###
## Risultati Prot2D singola omica "BP" ##
# alpha = 0.05 -> 1 pathway: GO:0002181
# alpha = 0.1  -> 2 pathway: GO:0002181, GO:0006412
# alpha = 0.25 -> 2 pathway: GO:0002181, GO:0006412

## Risultati Prot3D singola omica "BP" ##
# alpha = 0.05 -> 0 pathway
# alpha = 0.1  -> 0 pathway
# alpha = 0.25 -> 0 pathway
# OSS: bisogna spostare alpha = 0.37 per avere 2 pathway ##

### ------------------------ ###
### "CC": Cellular Component ###
### ------------------------ ###
## Risultati Prot2D singola omica "CC" ##
# alpha = 0.05 -> 14 pathway: GO:0044391  GO:0005840  GO:0015934  GO:0022626  GO:0015935  GO:0000313  GO:0005761  GO:0098798  GO:0022627
#                            GO:1990904  GO:0000315  GO:0005762  GO:0022625  GO:0032040

# alpha = 0.1  -> 15 pathway: 14 sopra + GO:0072562
# alpha = 0.25 -> 16 pathway: 15 sopra + GO:0030684

## Risultati Prot3D singola omica "CC" ##
# alpha = 0.05 -> 0 pathway
# alpha = 0.1  -> 0 pathway
# alpha = 0.25 -> 0 pathway
## OSS: Nemmeno spostando alpha = 0.5 isi trovano dei pathway ##

### ------------------------ ###
### "MF": Molecular Function ###
### ------------------------ ###
## Risultati Prot2D singola omica "MF" ##
# alpha = 0.05 -> 2 pathway: GO:0003735  GO:0005198

# alpha = 0.1  -> 2 pathway:  stessi di sopra 
# alpha = 0.25 -> 2 pathway:  stessi di sopra 

## Risultati Prot3D singola omica "MF" ##
# alpha = 0.05 -> 0 pathway
# alpha = 0.1  -> 0 pathway
# alpha = 0.25 -> 0 pathway
## OSS: NON trova nessun significativo nemmeno mettendo alpha = 0.5 ##


#######################################
####    Integrazione Proteomica    ####
#######################################
### ------------------------ ###
### "BP": Biological Process ###
### ------------------------ ###
## Risultati  "BP" ##
# alpha = 0.05 ->  pathway: 3
# alpha = 0.1  ->  pathway: soliti 3
# alpha = 0.25 ->  pathway: soliti 3


### ------------------------ ###
### "CC": Cellular Component ###
### ------------------------ ###
## Risultati Prot2D singola omica "CC" ##
# alpha = 0.05 ->  pathway: 19
# alpha = 0.1  ->  pathway: 19 + 2 = 21
# alpha = 0.25 ->  pathway: 21 + 1 = 22


### ------------------------ ###
### "MF": Molecular Function ###
### ------------------------ ###
## Risultati "MF" ##
# alpha = 0.05 ->  pathway: 2
# alpha = 0.1  ->  pathway: soliti 2
# alpha = 0.25 ->  pathway: 2 + 1 = 3



#############################################################################
####                  GSEA (Gene Set Enrichment Analysis)                ####
#############################################################################

# ==============================================================================
# 0. PREPARAZIONE LIBRERIE PATHWAY (GO: Biological Process)
# ==============================================================================
# per uniformarsi allo studio già pubblicato e confrontare i risultat
msig_c5 <- msigdbr(species = "Homo sapiens", collection = "C5")

# Estrazione C2 (KEGG e REACTOME) usando la colonna corretta gs_subcollection
msig_c2 <- msigdbr(species = "Homo sapiens", collection = "C2") %>% 
  filter(grepl("KEGG|REACTOME", gs_subcollection, ignore.case = TRUE))

# Unione delle collezioni
msig_all <- rbind(msig_c5, msig_c2)

# Liste Gene Set per ENSEMBL e SYMBOL
pathways_list_ensembl <- split(x = msig_all$ensembl_gene, f = msig_all$gs_name)
pathways_list_symbol  <- split(x = msig_all$gene_symbol,  f = msig_all$gs_name)
# ==============================================================================
# UTILITY FUNCTIONS
# ==============================================================================
clean_ensembl <- function(ids) {
  ids <- sub("\\..*", "", ids)
  return(unique(na.omit(ids)))
}

crea_mappa_ensembl_to_symbol <- function(ensembl_ids) {
  clean_ids <- clean_ensembl(ensembl_ids)
  mapped <- suppressMessages(
    clusterProfiler::bitr(
      geneID   = clean_ids,
      fromType = "ENSEMBL",
      toType   = "SYMBOL",
      OrgDb    = org.Hs.eg.db
    )
  )
  setNames(mapped$SYMBOL, mapped$ENSEMBL)
}

get_plsda_rank <- function(mat_X, vec_Y) {
  mod <- plsda(X = mat_X, Y = vec_Y, ncomp = 1)
  sel <- selectVar(mod, comp = 1)$value
  vec_loadings <- setNames(sel$value, rownames(sel))
  vec_loadings <- vec_loadings + runif(length(vec_loadings), -1e-12, 1e-12)
  return(sort(vec_loadings, decreasing = TRUE))
}

run_gsea_on_block_plsda <- function(X_list, Y_vec, pathways, 
                                    minSize = NULL, maxSize = 500, 
                                    fuse_blocks = FALSE, 
                                    map_prot_to_gene = NULL) {
  
  if (is.null(names(X_list))) {
    names(X_list) <- paste0("Block_", seq_along(X_list))
  }
  
  design_mat <- matrix(0.1, ncol = length(X_list), nrow = length(X_list), 
                       dimnames = list(names(X_list), names(X_list)))
  diag(design_mat) <- 0
  
  mod_full <- block.plsda(X = X_list, Y = as.factor(Y_vec), ncomp = 1, design = design_mat)
  
  # ============================================================================
  # OPZIONE A: FUSIONE DEI LOADINGS (fuse_blocks = TRUE)
  # ============================================================================
  if (fuse_blocks) {
    cat("  -> Esecuzione GSEA su Loadings Z-standardizzati e Fusi...\n")
    
    df_list <- list()
    
    for (nome_blocco in names(X_list)) {
      loadings_raw <- mod_full$loadings[[nome_blocco]][, 1]
      
      # Z-standardizzazione intra-blocco
      z_scores <- as.numeric(scale(loadings_raw))
      names(z_scores) <- names(loadings_raw)
      
      # CORREZIONE BUG 1: Rimappatura del blocco Trascrittomico (Ensembl -> Symbol)
      is_tx_block <- grepl("trascr|transcriptomics", nome_blocco, ignore.case = TRUE)
      
      if (is_tx_block && !is.null(map_prot_to_gene)) {
        cat("     * Rimappatura ENSEMBL -> SYMBOL per:", nome_blocco, "\n")
        
        # Pulisce gli ID Ensembl rimuovendo la versione decimale (.1, .2)
        clean_ids <- clean_ensembl(names(z_scores))
        mapped_names <- map_prot_to_gene[clean_ids]
        
        valid_idx <- !is.na(mapped_names) & mapped_names != ""
        z_scores <- z_scores[valid_idx]
        names(z_scores) <- mapped_names[valid_idx]
      }
      
      df_list[[nome_blocco]] <- data.frame(
        Gene = names(z_scores),
        Z    = z_scores,
        stringsAsFactors = FALSE
      ) %>%
        group_by(Gene) %>%
        summarise(Z = mean(Z), .groups = "drop")
    }
    
    # Merging e media degli Z-score (Usa Reduce di R Base per evitare conflitti con IRanges)
    df_fused <- Reduce(function(x, y) full_join(x, y, by = "Gene"), df_list)
    z_cols   <- setdiff(colnames(df_fused), "Gene")
    
    df_fused <- df_fused %>%
      rowwise() %>%
      mutate(Z_integrated = mean(c_across(all_of(z_cols)), na.rm = TRUE)) %>%
      ungroup() %>%
      arrange(desc(Z_integrated))
    
    ranking_vec <- df_fused$Z_integrated
    names(ranking_vec) <- df_fused$Gene
    ranking_vec <- ranking_vec + runif(length(ranking_vec), -1e-12, 1e-12)
    ranking_vec <- sort(ranking_vec, decreasing = TRUE)
    
    curr_minSize <- if (is.null(minSize)) 15 else minSize
    
    gsea_res <- fgsea(
      pathways  = pathways,
      stats     = ranking_vec,
      minSize   = curr_minSize,
      maxSize   = maxSize,
      scoreType = "std"
    )
    
    return(list(
      gsea     = gsea_res,
      rankings = ranking_vec
    ))
  }
  
  # ============================================================================
  # OPZIONE B: SEPARATA PER SINGOLO BLOCCO (fuse_blocks = FALSE)
  # ============================================================================
  cat("  -> Esecuzione GSEA SEPARATA per ciascun blocco...\n")
  
  gsea_results  <- list()
  rankings_list <- list()
  
  for (i in seq_along(X_list)) {
    nome_blocco <- names(X_list)[i]
    cat("     * Elaborazione blocco:", nome_blocco, "\n")
    
    if (is.list(pathways) && nome_blocco %in% names(pathways)) {
      curr_pathways <- pathways[[nome_blocco]]
    } else if (is.list(pathways) && is.list(pathways[[1]]) && i <= length(pathways)) {
      curr_pathways <- pathways[[i]]
    } else {
      curr_pathways <- pathways
    }
    
    if (is.null(minSize)) {
      is_prot <- grepl("prot", nome_blocco, ignore.case = TRUE)
      curr_minSize <- if (is_prot) 5 else 15
    } else {
      curr_minSize <- minSize
    }
    
    loadings_raw <- mod_full$loadings[[nome_blocco]][, 1]
    
    ranking_vec <- as.numeric(scale(loadings_raw))
    names(ranking_vec) <- names(loadings_raw)
    ranking_vec <- ranking_vec + runif(length(ranking_vec), -1e-12, 1e-12)
    ranking_vec <- sort(ranking_vec, decreasing = TRUE)
    
    rankings_list[[nome_blocco]] <- ranking_vec
    
    gsea_results[[nome_blocco]] <- fgsea(
      pathways  = curr_pathways,
      stats     = ranking_vec,
      minSize   = curr_minSize,
      maxSize   = maxSize,
      scoreType = "std"
    )
  }
  
  return(list(
    gsea     = gsea_results,
    rankings = rankings_list
  ))
}

clean_gsea_results_aligned <- function(gsea_res, fdr_cutoff = 0.2, nes_cutoff = 1.0) {
  if (is.null(gsea_res) || nrow(gsea_res) == 0) {
    return(data.frame(pathway = character(), pval = numeric(), padj = numeric(), 
                      NES = numeric(), size = integer(), leadingEdge_str = character()))
  }
  
  df_clean <- gsea_res %>%
    as.data.frame() %>%
    filter(!is.na(padj), !is.na(NES)) %>%
    filter(padj < fdr_cutoff, abs(NES) >= nes_cutoff) %>%
    arrange(padj, desc(abs(NES)))
  
  if (nrow(df_clean) == 0) {
    return(data.frame(pathway = character(), pval = numeric(), padj = numeric(), 
                      NES = numeric(), size = integer(), leadingEdge_str = character()))
  }
  
  df_clean <- df_clean %>%
    mutate(leadingEdge_str = sapply(leadingEdge, function(x) paste(x, collapse = ", "))) %>%
    select(pathway, pval, padj, NES, size, leadingEdge_str)
  
  return(df_clean)
}

# ==============================================================================
# (1) LIVELLO SINGOLO BLOCCO
# ==============================================================================
rank_trascr2D <- get_plsda_rank(t(trascr_filtered2D_scaled), Y[1:14])
rank_trascr3D <- get_plsda_rank(t(trascr_filtered3D_scaled), Y[1:14])
rank_prot2D   <- get_plsda_rank(t(prot2D_scaled),   Y[1:14])
rank_prot3D   <- get_plsda_rank(t(prot3D_scaled),   Y[1:14])

gsea_trascr2D <- fgsea(pathways = pathways_list_ensembl, stats = rank_trascr2D, minSize = 15, maxSize = 500)
gsea_trascr3D <- fgsea(pathways = pathways_list_ensembl, stats = rank_trascr3D, minSize = 15, maxSize = 500)
gsea_prot2D   <- fgsea(pathways = pathways_list_symbol,  stats = rank_prot2D,   minSize = 5,  maxSize = 500)
gsea_prot3D   <- fgsea(pathways = pathways_list_symbol,  stats = rank_prot3D,   minSize = 5,  maxSize = 500)


# ==============================================================================
# (2) LIVELLO DI INTEGRAZIONE PER OMICA (2D vs 3D)
# ==============================================================================
# A. Trascrittomica (trascr2D + trascr3D)
gsea_integ_omica_trascr <- run_gsea_on_block_plsda(
  X_list        = X_Trascr, 
  Y_vec         = Y, 
  pathways      = pathways_list_ensembl, 
  fuse_blocks = FALSE
)

# B. Proteomica (prot2D + prot3D)
gsea_integ_omica_prot <- run_gsea_on_block_plsda(
  X_list        = X_Prot, 
  Y_vec         = Y, 
  pathways      = pathways_list_symbol, 
  fuse_blocks = FALSE
)


# ==============================================================================
# (3) LIVELLO DI INTEGRAZIONE PER DIMENSIONE (Trascrittomica vs Proteomica)
# ==============================================================================
# Mappatura specifica delle librerie per ciascun blocco all'interno del modello integrato
pathways_mix <- list(
  transcriptomics = pathways_list_ensembl,
  proteomics      = pathways_list_symbol
)

# A. Stato 2D (Trascrittomica + Proteomica)
gsea_integ_2D <- run_gsea_on_block_plsda(
  X_list   = X2D, 
  Y_vec    = Y, 
  pathways = pathways_mix,
  fuse_blocks = FALSE
)

# B. Stato 3D (Trascrittomica + Proteomica)
gsea_integ_3D <- run_gsea_on_block_plsda(
  X_list   = X3D, 
  Y_vec    = Y, 
  pathways = pathways_mix,
  fuse_blocks = FALSE
)


# Studio Singole omiche
top_pathways_trascr2D_from_single <- clean_gsea_results_aligned(gsea_trascr2D)
top_pathways_trascr3D_from_single <- clean_gsea_results_aligned(gsea_trascr3D)
top_pathways_prot2D_from_single <- clean_gsea_results_aligned(gsea_prot2D)
top_pathways_prot3D_from_single <- clean_gsea_results_aligned(gsea_prot3D)

nrow(top_pathways_trascr2D_from_single) # 96
nrow(top_pathways_trascr3D_from_single) # 1216
nrow(top_pathways_prot2D_from_single)   # 1120
nrow(top_pathways_prot3D_from_single)   # 877

# Studio Integrazione Trascr
top_pathways_2D_from_trascr <- clean_gsea_results_aligned(gsea_integ_omica_trascr$gsea[[1]])
top_pathways_3D_from_trascr <- clean_gsea_results_aligned(gsea_integ_omica_trascr$gsea[[2]])

nrow(top_pathways_2D_from_trascr) # 109
nrow(top_pathways_3D_from_trascr) # 1926

# Studio Integrazione Prot
top_pathways_2D_from_prot <- clean_gsea_results_aligned(gsea_integ_omica_prot$gsea[[1]])
top_pathways_3D_from_prot   <- clean_gsea_results_aligned(gsea_integ_omica_prot$gsea[[2]])

nrow(top_pathways_2D_from_prot) # 1728
nrow(top_pathways_3D_from_prot) # 1697

# Studio Integrazione 2D
top_pathways_trascr_from_2D <- clean_gsea_results_aligned(gsea_integ_2D$gsea[[1]])
top_pathways_prot_from_2D  <- clean_gsea_results_aligned(gsea_integ_2D$gsea[[2]])

nrow(top_pathways_trascr_from_2D) # 432
nrow(top_pathways_prot_from_2D)   # 2309

# Studio Integrazione 3D
top_pathways_trascr_from_3D <- clean_gsea_results_aligned(gsea_integ_3D$gsea[[1]])
top_pathways_prot_from_3D   <- clean_gsea_results_aligned(gsea_integ_3D$gsea[[2]])

nrow(top_pathways_trascr_from_3D) # 1392
nrow(top_pathways_prot_from_3D)   # 1014

####################################
#### Intersezione di Pathway GO ####
####################################

## Fissata omica e dimensione ##
pathways_comuni_trascr2D <- Reduce(intersect, list(top_pathways_trascr2D_from_single$pathway, 
                                                   top_pathways_2D_from_trascr$pathway, 
                                                   top_pathways_trascr_from_2D$pathway))
cat("Numero di pathway condivisi da TUTTI e 3 i livelli per Trascrittomica 2D:", length(pathways_comuni_trascr2D), "\n") # 40
print(pathways_comuni_trascr2D)

pathways_comuni_trascr3D <- Reduce(intersect, list(top_pathways_trascr3D_from_single$pathway, 
                                                   top_pathways_3D_from_trascr$pathway, 
                                                   top_pathways_trascr_from_3D$pathway))
cat("Numero di pathway condivisi da TUTTI e 3 i livelli per Trascrittomica 3D:", length(pathways_comuni_trascr3D), "\n") # 928
print(pathways_comuni_trascr3D)


pathways_comuni_prot2D <- Reduce(intersect, list(top_pathways_prot2D_from_single$pathway, 
                                                   top_pathways_2D_from_prot$pathway, 
                                                   top_pathways_prot_from_2D$pathway))
cat("Numero di pathway condivisi da TUTTI e 3 i livelli per Proteomica 2D:", length(pathways_comuni_prot2D), "\n") # 693
print(pathways_comuni_prot2D)


pathways_comuni_prot3D <- Reduce(intersect, list(top_pathways_prot3D_from_single$pathway, 
                                                   top_pathways_3D_from_prot$pathway, 
                                                   top_pathways_prot_from_3D$pathway))
cat("Numero di pathway condivisi da TUTTI e 3 i livelli per Proteomica 3D:", length(pathways_comuni_prot3D), "\n") # 606
print(pathways_comuni_prot3D)

## Fissata solo il tipo di omica ##
pathways_comuni_trascr <- Reduce(intersect, list(top_pathways_trascr2D_from_single$pathway, 
                                                 top_pathways_2D_from_trascr$pathway, 
                                                 top_pathways_trascr_from_2D$pathway,
                                                 top_pathways_trascr3D_from_single$pathway, 
                                                 top_pathways_3D_from_trascr$pathway, 
                                                 top_pathways_trascr_from_3D$pathway))
cat("Numero di pathway condivisi da TUTTI e 3 i livelli per Trascrittomica (indipendentemente dalla dimensione):", length(pathways_comuni_trascr), "\n") # 7
print(pathways_comuni_trascr)

pathways_comuni_prot <- Reduce(intersect, list(top_pathways_prot2D_from_single$pathway, 
                                               top_pathways_2D_from_prot$pathway, 
                                               top_pathways_prot_from_2D$pathway,
                                               top_pathways_prot3D_from_single$pathway, 
                                               top_pathways_3D_from_prot$pathway, 
                                               top_pathways_prot_from_3D$pathway))
cat("Numero di pathway condivisi da TUTTI e 3 i livelli per Proteomica (indipendentemente dalla dimensione):", length(pathways_comuni_prot), "\n") # 175
print(pathways_comuni_prot)


####################################################
### Confronto risultati con pubblicazione Micali ###
####################################################
file_correl <- "/Users/ceciliabergamini/Desktop/TESI/script /file_confronto_risultati_Micale/Supplementary_Table_GSEA_2D_Correl.xlsx"
file_gse    <- "/Users/ceciliabergamini/Desktop/TESI/script /file_confronto_risultati_Micale/Supplementary_Table_GSEA_2D_GSE27019.xlsx"

get_pathways_single_sheet <- function(file_path) {
  
  df_sheet <- read_excel(file_path, sheet = 1, col_names = FALSE)
  flat_vector <- as.character(unlist(df_sheet))
  pathways_found <- flat_vector[grepl("^(GOBP|GOCC|GOMF|KEGG|REACTOME)_", flat_vector)]
  
  return(unique(pathways_found))
}


pathways_paper_correl <- get_pathways_single_sheet(file_correl)
pathways_paper_gse    <- get_pathways_single_sheet(file_gse)
pathways_paper_totali <- unique(c(pathways_paper_correl, pathways_paper_gse))
pathways_from_single_trascr2D <- top_pathways_trascr2D_from_single$pathway
pathways2D_from_integrazione_trascr <-  top_pathways_2D_from_trascr$pathway
pathways_trascr2D_from_integration <- top_pathways_trascr_from_2D$pathway

# 1. Confronto SEPARATO con File 1 (COL3A1-correlati)
comuni_correl <- intersect(pathways_paper_correl, pathways_from_single_trascr2D)
#comuni_correl <- intersect(pathways_paper_correl, pathways2D_from_integrazione_trascr)
comuni_correl <- intersect(pathways_paper_correl, pathways_trascr2D_from_integration)
pct_correl    <- (length(comuni_correl) / length(pathways_paper_correl)) * 100

# 2. Confronto SEPARATO con File 2 (GSE270199 globali)
#comuni_gse <- intersect(pathways_paper_gse, pathways_from_single_trascr2D)
#comuni_gse <- intersect(pathways_paper_gse, pathways2D_from_integrazione_trascr)
comuni_gse <- intersect(pathways_paper_gse, pathways_trascr2D_from_integration)
pct_gse    <- (length(comuni_gse) / length(pathways_paper_gse)) * 100

# 3. Confronto AGGREGATO (Unione dei due file del paper)
pathways_totali_paper <- unique(c(pathways_paper_correl, pathways_paper_gse))
#comuni_totali         <- intersect(pathways_totali_paper, pathways_from_single_trascr2D)
#comuni_totali         <- intersect(pathways_totali_paper, pathways2D_from_integrazione_trascr)
comuni_totali         <- intersect(pathways_totali_paper, pathways_trascr2D_from_integration)
pct_totale            <- (length(comuni_totali) / length(pathways_totali_paper)) * 100

# Stampa del Summary
cat("\n======================================================================\n")
cat("                SUMMARY CONFRONTO PATHWAY GSEA                        \n")
cat("======================================================================\n")
cat("1. COL3A1-Correlati (File 1): ", length(comuni_correl), "/", length(pathways_paper_correl), 
    " (", round(pct_correl, 1), "%)\n")
cat("2. Globali 2D (File 2):      ", length(comuni_gse), "/", length(pathways_paper_gse), 
    " (", round(pct_gse, 1), "%)\n")
cat("3. UNIONE TOTALE (Paper):    ", length(comuni_totali), "/", length(pathways_totali_paper), 
    " (", round(pct_totale, 1), "%)\n")
cat("======================================================================\n")


# ==============================================================================
# 1. COSTRUZIONE DATA FRAME
# ==============================================================================

df_foglio1 <- bind_rows(
  top_pathways_trascr2D_from_single %>% mutate(Omica = "Trascrittomica 2D"),
  top_pathways_trascr3D_from_single %>% mutate(Omica = "Trascrittomica 3D"),
  top_pathways_prot2D_from_single   %>% mutate(Omica = "Proteomica 2D"),
  top_pathways_prot3D_from_single   %>% mutate(Omica = "Proteomica 3D")
) %>% select(Omica, everything())

df_foglio2 <- bind_rows(
  top_pathways_2D_from_trascr %>% mutate(Integrazione = "Trascrittomica (2D)"),
  top_pathways_3D_from_trascr %>% mutate(Integrazione = "Trascrittomica (3D)"),
  top_pathways_2D_from_prot   %>% mutate(Integrazione = "Proteomica (2D)"),
  top_pathways_3D_from_prot   %>% mutate(Integrazione = "Proteomica (3D)")
) %>% select(Integrazione, everything())

df_foglio3 <- bind_rows(
  top_pathways_trascr_from_2D %>% mutate(Integrazione = "Stato 2D (Trascrittomica)"),
  top_pathways_prot_from_2D   %>% mutate(Integrazione = "Stato 2D (Proteomica)"),
  top_pathways_trascr_from_3D %>% mutate(Integrazione = "Stato 3D (Trascrittomica)"),
  top_pathways_prot_from_3D   %>% mutate(Integrazione = "Stato 3D (Proteomica)")
) %>% select(Integrazione, everything())

## Foglio 4: contiene le diverse intersezioni (per identificare i core signature) ##
# A. Fissata Omica e Dimensione
df_intersez_omica_dim <- bind_rows(
  data.frame(Categoria = "Fissata Omica e Dimensione", Condizione = "Trascrittomica 2D", pathway = pathways_comuni_trascr2D),
  data.frame(Categoria = "Fissata Omica e Dimensione", Condizione = "Trascrittomica 3D", pathway = pathways_comuni_trascr3D),
  data.frame(Categoria = "Fissata Omica e Dimensione", Condizione = "Proteomica 2D",     pathway = pathways_comuni_prot2D),
  data.frame(Categoria = "Fissata Omica e Dimensione", Condizione = "Proteomica 3D",     pathway = pathways_comuni_prot3D)
)

# B. Fissata solo l'Omica - Indipendente dalla Dimensione (8, 188)
df_intersez_omica_only <- bind_rows(
  data.frame(Categoria = "Fissata solo Omica", Condizione = "Trascrittomica Global (2D+3D)", pathway = pathways_comuni_trascr),
  data.frame(Categoria = "Fissata solo Omica", Condizione = "Proteomica Global (2D+3D)",     pathway = pathways_comuni_prot)
)

# Unione per il Foglio 4
df_foglio4 <- bind_rows(df_intersez_omica_dim, df_intersez_omica_only)

# ==============================================================================
# 2. CREAZIONE E SALVATAGGIO WORKBOOK 
# ==============================================================================

wb <- createWorkbook()

addWorksheet(wb, "1_Singola_Omica")
addWorksheet(wb, "2_Integrazione_Omica")
addWorksheet(wb, "3_Integrazione_Dimensione")
addWorksheet(wb, "4_Intersezioni_Core")

writeData(wb, sheet = "1_Singola_Omica",           x = df_foglio1)
writeData(wb, sheet = "2_Integrazione_Omica",      x = df_foglio2)
writeData(wb, sheet = "3_Integrazione_Dimensione", x = df_foglio3)
writeData(wb, sheet = "4_Intersezioni_Core",       x = df_foglio4)

# Salvataggio sicuro
percorso_output <- "/Users/ceciliabergamini/Desktop/TESI/script /Risultati_GSEA_MultiOmica.xlsx"

saveWorkbook(wb, file = percorso_output, overwrite = TRUE)

cat("\n[✓] File Excel salvato con successo in:\n", percorso_output, "\n")


######################################################################
## Adattamento del medesimo approccio per gli altri 3 blocchi omici ##
######################################################################

# ==============================================================================
# FUNZIONE PER GSEA BASATA SU CORRELAZIONE DI PEARSON CON COL3A1 (Stile Paper)
# ==============================================================================
col3a1_ensembl <- mapIds(
  org.Hs.eg.db,
  keys      = "COL3A1",
  column    = "ENSEMBL",
  keytype   = "SYMBOL",
  multiVals = "first"
)

cat("L'ID Ensembl di COL3A1 è:", col3a1_ensembl, "\n")

# Estrazione del profilo usando l'ID Ensembl
if (!is.na(col3a1_ensembl) && col3a1_ensembl %in% rownames(trascr2D_common)) {
  col3a1_profile <- trascr2D_common[col3a1_ensembl, ]
  cat("Profilo di COL3A1 estratto con successo\n")
} else {
  cat("L'ID Ensembl di COL3A1 non è stato trovato tra le righe di trascr2D_common.\n")
}

get_pearson_col3a1_rank <- function(mat_data, col3a1_expr) {
  # mat_data ha le feature sulle righe e i campioni sulle colonne
  # Calcolo della correlazione di Pearson tra ogni gene/proteina e COL3A1
  cor_vec <- apply(mat_data, 1, function(x) cor(x, col3a1_expr, method = "pearson", use = "complete.obs"))
  
  # Rimuove eventuali NA e aggiunge jitter per evitare ties
  cor_vec <- cor_vec[!is.na(cor_vec)]
  cor_vec <- cor_vec + runif(length(cor_vec), -1e-12, 1e-12)
  
  return(sort(cor_vec, decreasing = TRUE))
}

# 1. Estrazione del profilo di espressione di COL3A1 (dalla Trascrittomica 2D)
col3a1_id_clean <- unname(col3a1_ensembl)
col3a1_profile <- trascr2D_common[col3a1_id_clean, ]
cat("[✓] Estratti con successo", length(col3a1_profile), "campioni per COL3A1 (", col3a1_id_clean, ")\n")

# 2. Generazione dei Ranking Univariati per Pearson su tutti e 4 i blocchi
rank_correl_trascr2D <- get_pearson_col3a1_rank(trascr2D_common, col3a1_profile)
rank_correl_trascr3D <- get_pearson_col3a1_rank(trascr3D_common, col3a1_profile)
rank_correl_prot2D   <- get_pearson_col3a1_rank(prot2D_common,   col3a1_profile)
rank_correl_prot3D   <- get_pearson_col3a1_rank(prot3D_common,   col3a1_profile)

# 3. Esecuzione della GSEA Univariata
gsea_correl_trascr2D <- fgsea(pathways = pathways_list_ensembl, stats = rank_correl_trascr2D, minSize = 15, maxSize = 500)
gsea_correl_trascr3D <- fgsea(pathways = pathways_list_ensembl, stats = rank_correl_trascr3D, minSize = 15, maxSize = 500)
gsea_correl_prot2D   <- fgsea(pathways = pathways_list_symbol,  stats = rank_correl_prot2D,   minSize = 5,  maxSize = 500)
gsea_correl_prot3D   <- fgsea(pathways = pathways_list_symbol,  stats = rank_correl_prot3D,   minSize = 5,  maxSize = 500)

# 4. Pulizia Risultati (FDR < 0.2)
top_correl_trascr2D <- clean_gsea_results_aligned(gsea_correl_trascr2D)
top_correl_trascr3D <- clean_gsea_results_aligned(gsea_correl_trascr3D)
top_correl_prot2D   <- clean_gsea_results_aligned(gsea_correl_prot2D)


# Conteggio pathway significativi per blocco
cat("Pathways correlati a COL3A1 trovati:\n")
cat(" - Trascrittomica 2D:", nrow(top_correl_trascr2D), "\n") # 851
cat(" - Trascrittomica 3D:", nrow(top_correl_trascr3D), "\n") # 15
cat(" - Proteomica 2D:    ", nrow(top_correl_prot2D), "\n")   # 534

# Top 5 pathway Co-espressi (NES > 0) nella Trascrittomica 2D
head(top_correl_trascr2D %>% filter(NES > 0) %>% select(pathway, NES, padj), 5)

# Top 5 pathway Inversamente correlati (NES < 0) nella Trascrittomica 2D
head(top_correl_trascr2D %>% filter(NES < 0) %>% select(pathway, NES, padj), 5)


######################################################################
##                     Export risultati GSEA                        ##
######################################################################

# ==============================================================================
# 0. CONFIGURAZIONE PERCORSI E PARAMETRI
# ==============================================================================
cartella_output <- "/Users/ceciliabergamini/Desktop/TESI/script"    
fdr_label       <- "padj < 0.2 e |NES| >= 1"

# ==============================================================================
# 1. FUNZIONE PER LA SCRITTURA FORMATTATA DEI FOGLI EXCEL
# ==============================================================================
scrivi_foglio_gsea <- function(wb, sheet, gruppi, titolo_foglio = NULL) {
  addWorksheet(wb, sheet)
  nc <- 6   # pathway, pval, padj, NES, size, leadingEdge_str
  
  # Stili di formattazione
  st_titolo  <- createStyle(textDecoration = "bold", fontSize = 16)
  st_memo    <- createStyle(textDecoration = "italic", fontSize = 10, fontColour = "#595959")
  st_summary <- createStyle(textDecoration = "bold", fontSize = 11, fgFill = "#DDEBF7")
  st_gruppo  <- createStyle(textDecoration = "bold", fontSize = 14,
                            fgFill = "#1F4E78", fontColour = "#FFFFFF")
  st_blocco  <- createStyle(textDecoration = "bold", fgFill = "#DDEBF7")
  st_head    <- createStyle(textDecoration = "bold", border = "bottom")
  
  riga <- 1
  
  # 1. INTESTAZIONE: Titolo, Memo sulle Soglie e Conteggio Totale
  if (!is.null(titolo_foglio)) {
    # Calcolo del totale complessivo dei pathway significativi nel foglio
    totale_pathways_foglio <- sum(sapply(gruppi, function(g) {
      sum(sapply(g, nrow))
    }))
    
    # Titolo
    writeData(wb, sheet, titolo_foglio, startRow = riga, colNames = FALSE)
    addStyle(wb, sheet, st_titolo, rows = riga, cols = 1, stack = TRUE)
    riga <- riga + 1
    
    # Memo sui parametri di filtraggio
    writeData(wb, sheet, paste("Parametri di filtraggio:", fdr_label), startRow = riga, colNames = FALSE)
    addStyle(wb, sheet, st_memo, rows = riga, cols = 1, stack = TRUE)
    riga <- riga + 1
    
    # Casella con il Conteggio Totale
    summary_txt <- paste0("Totale pathway significativi trovati nel foglio: ", totale_pathways_foglio)
    writeData(wb, sheet, summary_txt, startRow = riga, colNames = FALSE)
    addStyle(wb, sheet, st_summary, rows = riga, cols = 1:nc, stack = TRUE)
    
    riga <- riga + 3
  }
  
  # 2. TABELLE DEI GRUPPI E BLOCCHI
  for (g in names(gruppi)) {
    writeData(wb, sheet, g, startRow = riga, colNames = FALSE)
    addStyle(wb, sheet, st_gruppo, rows = riga, cols = 1:nc, stack = TRUE)
    riga <- riga + 2
    
    for (b in names(gruppi[[g]])) {
      d <- gruppi[[g]][[b]]
      
      # Sub-header per ciascun blocco con il rispettivo conteggio
      writeData(wb, sheet, paste0(b, "   (", nrow(d), " pathway)"),
                startRow = riga, colNames = FALSE)
      addStyle(wb, sheet, st_blocco, rows = riga, cols = 1:nc, stack = TRUE)
      riga <- riga + 1
      
      if (nrow(d) == 0) {
        writeData(wb, sheet, "nessun pathway significativo", startRow = riga, colNames = FALSE)
        riga <- riga + 3
        next
      }
      
      writeData(wb, sheet, d, startRow = riga)
      addStyle(wb, sheet, st_head, rows = riga, cols = 1:nc, stack = TRUE)
      riga <- riga + nrow(d) + 3     # due righe vuote tra un blocco e il successivo
    }
    riga <- riga + 2                 # spazio extra tra un gruppo e il successivo
  }
  
  setColWidths(wb, sheet, cols = 1:nc, widths = c(60, 11, 11, 9, 7, 80))
  freezePane(wb, sheet, firstRow = FALSE)
}

# ==============================================================================
# 2. DEFINIZIONE DELLE STRUTTURE DEI FOGLI
# ==============================================================================

# Foglio 1: le 4 omiche separate
gruppi_f1 <- list(
  "TRASCRITTOMICA - singola omica" = list(
    "Trascrittomica 2D" = top_pathways_trascr2D_from_single,
    "Trascrittomica 3D" = top_pathways_trascr3D_from_single
  ),
  "PROTEOMICA - singola omica" = list(
    "Proteomica 2D" = top_pathways_prot2D_from_single,
    "Proteomica 3D" = top_pathways_prot3D_from_single
  )
)

# Foglio 2: integrazione fissando la DIMENSIONE (trascrittomica + proteomica)
gruppi_f2 <- list(
  "STATO 2D - integrazione Trascrittomica + Proteomica" = list(
    "Trascrittomica 2D (da integrazione)" = top_pathways_trascr_from_2D,
    "Proteomica 2D (da integrazione)"     = top_pathways_prot_from_2D
  ),
  "STATO 3D - integrazione Trascrittomica + Proteomica" = list(
    "Trascrittomica 3D (da integrazione)" = top_pathways_trascr_from_3D,
    "Proteomica 3D (da integrazione)"     = top_pathways_prot_from_3D
  )
)

# Foglio 3: integrazione fissando il TIPO DI OMICA (2D + 3D)
gruppi_f3 <- list(
  "TRASCRITTOMICA - integrazione 2D + 3D" = list(
    "Trascrittomica 2D (da integrazione)" = top_pathways_2D_from_trascr,
    "Trascrittomica 3D (da integrazione)" = top_pathways_3D_from_trascr
  ),
  "PROTEOMICA - integrazione 2D + 3D" = list(
    "Proteomica 2D (da integrazione)" = top_pathways_2D_from_prot,
    "Proteomica 3D (da integrazione)" = top_pathways_3D_from_prot
  )
)

# ==============================================================================
# 3. CREAZIONE E SALVATAGGIO SICURO DEL WORKBOOK
# ==============================================================================

wb1 <- createWorkbook()

scrivi_foglio_gsea(wb1, "1_Singole_Omiche",         gruppi_f1,
                   "GSEA su ranking dei loadings PLS-DA - singole omiche")
scrivi_foglio_gsea(wb1, "2_Integr_Fissa_Dimensione", gruppi_f2,
                   "GSEA su ranking dei loadings block PLS-DA - integrazione a dimensione fissata")
scrivi_foglio_gsea(wb1, "3_Integr_Fissa_Omica",      gruppi_f3,
                   "GSEA su ranking dei loadings block PLS-DA - integrazione a omica fissata")

# Verifica ed eventuale creazione della cartella di destinazione
if (!dir.exists(cartella_output)) {
  dir.create(cartella_output, recursive = TRUE)
  cat("[+] Creata la cartella di output:", cartella_output, "\n")
}

# Costruzione del percorso pulito per il salvataggio
percorso_file <- file.path(cartella_output, "Risultati_GSEA_MultiOmica.xlsx")

# Salvataggio con sovrascrittura abilitata
saveWorkbook(wb1, file = percorso_file, overwrite = TRUE)

# Controllo e conferma a schermo dell'esito
if (file.exists(percorso_file)) {
  cat("\n[✓] FILE SALVATO CON SUCCESSO IN:\n", normalizePath(percorso_file), "\n")
} else {
  cat("\n[!] ATTENZIONE: Impossibile scrivere nella cartella specificata. Salvo una copia di backup sulla Scrivania...\n")
  backup_path <- "/Users/ceciliabergamini/Desktop/Risultati_GSEA_MultiOmica.xlsx"
  saveWorkbook(wb1, file = backup_path, overwrite = TRUE)
  cat("[✓] File salvato direttamente sul Desktop:", backup_path, "\n")
}

## ------------------------------------------------------------------------ ##
## GSEA: ordinamento dei geni/proteine in base alla correlazione con COL3A1 ##
## ------------------------------------------------------------------------ ##
profilo_stessa_dimensione <- TRUE

col3a1_ensembl <- unname(mapIds(org.Hs.eg.db, keys = "COL3A1", column = "ENSEMBL",
                                keytype = "SYMBOL", multiVals = "first"))
stopifnot(!is.na(col3a1_ensembl),
          col3a1_ensembl %in% rownames(trascr2D_common),
          col3a1_ensembl %in% rownames(trascr3D_common))

profili_col3a1 <- list(
  "2D" = trascr2D_common[col3a1_ensembl, ],
  "3D" = trascr3D_common[col3a1_ensembl, ]
)
cat("COL3A1:", col3a1_ensembl, "| campioni 2D:", length(profili_col3a1[["2D"]]),
    "| campioni 3D:", length(profili_col3a1[["3D"]]), "\n")

get_pearson_col3a1_rank <- function(mat_data, col3a1_expr, method = "pearson") {
  # mat_data: feature sulle righe, campioni sulle colonne, STESSO ordine di col3a1_expr
  stopifnot(ncol(mat_data) == length(col3a1_expr))
  if (!is.null(colnames(mat_data)) && !is.null(names(col3a1_expr)) &&
      !identical(colnames(mat_data), names(col3a1_expr)))
    message("  ! nomi dei campioni diversi tra matrice e profilo COL3A1: ",
            "verifica che l'ordine dei campioni sia lo stesso")
  
  cor_vec <- suppressWarnings(
    apply(mat_data, 1, function(x) cor(x, col3a1_expr, method = method, use = "complete.obs"))
  )
  cor_vec <- cor_vec[!is.na(cor_vec)]
  cor_vec <- cor_vec + runif(length(cor_vec), -1e-12, 1e-12)   # evita pareggi
  sort(cor_vec, decreasing = TRUE)
}

blocchi_corr <- list(
  "Trascrittomica 3D" = list(mat = trascr3D_common, dim = "3D",
                             paths = pathways_list_ensembl, minSize = 15),
  "Proteomica 2D"     = list(mat = prot2D_common,   dim = "2D",
                             paths = pathways_list_symbol,  minSize = 5),
  "Proteomica 3D"     = list(mat = prot3D_common,   dim = "3D",
                             paths = pathways_list_symbol,  minSize = 5)
)

gsea_correl <- list()
for (nm in names(blocchi_corr)) {
  b <- blocchi_corr[[nm]]
  dim_profilo <- if (profilo_stessa_dimensione) b$dim else "2D"
  cat("  -> GSEA correlazione COL3A1:", nm, "(profilo COL3A1", dim_profilo, ")\n")
  
  set.seed(123)                                   # riproducibilita'
  rk <- get_pearson_col3a1_rank(b$mat, profili_col3a1[[dim_profilo]])
  gsea_correl[[nm]] <- fgsea(pathways = b$paths, stats = rk,
                             minSize = b$minSize, maxSize = 500)
}

top_correl <- lapply(gsea_correl, clean_gsea_results_aligned)

cat("\nPathway correlati a COL3A1 (", fdr_label, "):\n", sep = "")
for (nm in names(top_correl)) cat(" -", nm, ":", nrow(top_correl[[nm]]), "\n")

## Excel 2: titolo che riporta il calcolo
titolo_corr <- paste0("GSEA su ranking = correlazione di Pearson con COL3A1",
                      " (profilo COL3A1: ",
                      if (profilo_stessa_dimensione) "stessa dimensione del blocco" else "2D",
                      "); Trascrittomica 2D esclusa (studio gia' effettuato)")

gruppi_corr <- list(
  "RANKING: correlazione di Pearson con COL3A1" = top_correl
)

wb2 <- createWorkbook()
scrivi_foglio_gsea(wb2, "GSEA_Correl_COL3A1", gruppi_corr, titolo_corr)
saveWorkbook(wb2, file.path(cartella_output, "Risultati_GSEA_Correlazione_COL3A1.xlsx"),
             overwrite = TRUE)
cat("[+] Salvato: Risultati_GSEA_Correlazione_COL3A1.xlsx\n")



# ################################################################################
# #### GSEA con integrazione dei geni che provengono da studi di integrazione ####
# ################################################################################
# 
# mappa_tx_symbol <- crea_mappa_ensembl_to_symbol(colnames(X2D$transcriptomics))
# head(mappa_tx_symbol)
# 
# cat("Mappati con successo:", sum(!is.na(mappa_tx_symbol) & mappa_tx_symbol != ""), 
#     "su", length(mappa_tx_symbol), "\n")
# 
# head(colnames(X2D$transcriptomics))
# 
# 
# gsea_integ_2D <- run_gsea_on_block_plsda(
#   X_list           = X2D, 
#   Y_vec            = Y, 
#   pathways         = pathways_list_symbol, 
#   fuse_blocks      = TRUE,
#   map_prot_to_gene = mappa_tx_symbol
# )
# 
# gsea_integ_3D <- run_gsea_on_block_plsda(
#   X_list           = X3D, 
#   Y_vec            = Y, 
#   pathways         = pathways_list_symbol, 
#   fuse_blocks      = TRUE,
#   map_prot_to_gene = mappa_tx_symbol
# )
# 
# gsea_integ_omica_trascr <- run_gsea_on_block_plsda(
#   X_list      = X_Trascr, 
#   Y_vec       = Y, 
#   pathways    = pathways_list_ensembl, 
#   fuse_blocks = TRUE
# )
# 
# gsea_integ_omica_prot <- run_gsea_on_block_plsda(
#   X_list      = X_Prot, 
#   Y_vec       = Y, 
#   pathways    = pathways_list_symbol, 
#   fuse_blocks = TRUE
# )
# 
# top_fuso_2D     <- clean_gsea_results_aligned(gsea_integ_2D$gsea)
# top_fuso_3D     <- clean_gsea_results_aligned(gsea_integ_3D$gsea)
# top_fuso_trascr <- clean_gsea_results_aligned(gsea_integ_omica_trascr$gsea)
# top_fuso_prot   <- clean_gsea_results_aligned(gsea_integ_omica_prot$gsea)
# 
# # Foglio 1: Cross-Omica (Fissata Dimensione)
# df_foglio_fuso_dim <- bind_rows(
#   top_fuso_2D %>% mutate(Integrazione = "Integrazione 2D (Trascr + Prot)"),
#   top_fuso_3D %>% mutate(Integrazione = "Integrazione 3D (Trascr + Prot)")
# ) %>% select(Integrazione, everything())
# 
# # Foglio 2: Intra-Omica (Fissata Omica)
# df_foglio_fuso_omica <- bind_rows(
#   top_fuso_trascr %>% mutate(Integrazione = "Integrazione Trascrittomica (2D + 3D)"),
#   top_fuso_prot   %>% mutate(Integrazione = "Integrazione Proteomica (2D + 3D)")
# ) %>% select(Integrazione, everything())
# 
# # Helper sicuro per le intersezioni
# crea_df_intersezione <- function(eti_confronto, vett_pathways) {
#   n_elem <- length(vett_pathways)
#   if (n_elem == 0) {
#     return(data.frame(Confronto = eti_confronto, Conteggio = 0, 
#                       pathway = "Nessun pathway in comune", stringsAsFactors = FALSE))
#   } else {
#     return(data.frame(Confronto = eti_confronto, Conteggio = n_elem, 
#                       pathway = vett_pathways, stringsAsFactors = FALSE))
#   }
# }
# 
# pathways_comuni_2D_3D  <- intersect(top_fuso_2D$pathway, top_fuso_3D$pathway)
# pathways_comuni_Omiche <- intersect(top_fuso_trascr$pathway, top_fuso_prot$pathway)
# pathways_core_global   <- Reduce(intersect, list(
#   top_fuso_2D$pathway, top_fuso_3D$pathway,
#   top_fuso_trascr$pathway, top_fuso_prot$pathway
# ))
# 
# df_foglio_intersezioni <- bind_rows(
#   crea_df_intersezione("Cross-Dimensione (2D Fuso vs 3D Fuso)", pathways_comuni_2D_3D),
#   crea_df_intersezione("Cross-Omica (Trascr Fusa vs Prot Fusa)", pathways_comuni_Omiche),
#   crea_df_intersezione("GLOBAL CORE (Integrazione Totale)", pathways_core_global)
# )
# 
# wb <- createWorkbook()
# 
# # Stili per le intestazioni del foglio
# st_titolo  <- createStyle(textDecoration = "bold", fontSize = 14, fontColour = "#1F4E78")
# st_memo    <- createStyle(textDecoration = "italic", fontSize = 10, fontColour = "#595959")
# st_summary <- createStyle(textDecoration = "bold", fontSize = 11, fgFill = "#DDEBF7", border = "TopBottomLeftRight")
# st_head    <- createStyle(textDecoration = "bold", fgFill = "#1F4E78", fontColour = "#FFFFFF", halign = "center")
# 
# # Helper per scrivere Intestazione + Memo + Conteggio
# scrivi_intestazione_foglio <- function(wb, sheet, titolo_foglio, total_pathways_count) {
#   addWorksheet(wb, sheet)
#   
#   # 1. Titolo del Foglio
#   writeData(wb, sheet, titolo_foglio, startRow = 1, startCol = 1)
#   addStyle(wb, sheet, st_titolo, rows = 1, cols = 1)
#   
#   # 2. Memo sui Parametri di Filtraggio Usati
#   memo_text <- "Parametri di Filtraggio Applicati: p.adjust (FDR) < 0.20  |  |NES| >= 1.0  |  scoreType = 'std'"
#   writeData(wb, sheet, memo_text, startRow = 2, startCol = 1)
#   addStyle(wb, sheet, st_memo, rows = 2, cols = 1)
#   
#   # 3. Casella con Conteggio Totale Pathway Significativi
#   summary_text <- paste0("  [✓] Numero totale di pathway che hanno superato il test: ", total_pathways_count, "  ")
#   writeData(wb, sheet, summary_text, startRow = 4, startCol = 1)
#   addStyle(wb, sheet, st_summary, rows = 4, cols = 1)
# }
# 
# # --- FOGLIO 1: Integrazione Cross-Omica (Fissata la Dimensione 2D vs 3D) ---
# n_tot_fuso_dim <- nrow(df_foglio_fuso_dim)
# scrivi_intestazione_foglio(wb, "1_Integrazione_Dimensione", 
#                            "GSEA Multi-Omica Integrata (Fissata Dimensione 2D / 3D)", 
#                            n_tot_fuso_dim)
# writeData(wb, sheet = "1_Integrazione_Dimensione", x = df_foglio_fuso_dim, startRow = 6, headerStyle = st_head)
# 
# 
# # --- FOGLIO 2: Integrazione Intra-Omica (Fissata l'Omica: Trascrittomica vs Proteomica) ---
# n_tot_fuso_omica <- nrow(df_foglio_fuso_omica)
# scrivi_intestazione_foglio(wb, "2_Integrazione_Omica", 
#                            "GSEA Multi-Omica Integrata (Fissata Omica Trascrittomica / Proteomica)", 
#                            n_tot_fuso_omica)
# writeData(wb, sheet = "2_Integrazione_Omica", x = df_foglio_fuso_omica, startRow = 6, headerStyle = st_head)
# 
# 
# # --- FOGLIO 3: Intersezioni e Pathway Core Integrati ---
# n_tot_intersezioni <- nrow(df_foglio_intersezioni)
# scrivi_intestazione_foglio(wb, "3_Intersezioni_Core_Fusi", 
#                            "Intersezioni tra Modelli Integrati e Global Core Signature", 
#                            n_tot_intersezioni)
# writeData(wb, sheet = "3_Intersezioni_Core_Fusi", x = df_foglio_intersezioni, startRow = 6, headerStyle = st_head)
# 
# 
# # Adattamento automatico della larghezza delle colonne per tutti i fogli
# for (s in names(wb)) {
#   setColWidths(wb, sheet = s, cols = 1:10, widths = "auto")
# }
# 
# # Salvataggio Workbook
# percorso_output <- "/Users/ceciliabergamini/Desktop/TESI/script /Risultati_GSEA_MultiOmica_FUSI.xlsx"
# saveWorkbook(wb, file = percorso_output, overwrite = TRUE)
# 
# cat("\n[✓] Excel con Memo e Conteggi salvato con successo in:\n", percorso_output, "\n")

################################################################
####                   GSEA con WGCNA                       ####
################################################################

# ==============================================================================
# 1. FUNZIONE PER IL RANKING DI RETE WGCNA (Network Screening orientato)
# ==============================================================================
# Implementa fedelmente le equazioni della tesi:
# - GS_i = |cor(x_i, Y)|^gamma
# - a_ij = |cor(x_i, x_j)|^alpha  (con a_ii = 0 per j != i)
# - q_i  = GS_i * sum_{j != i} (a_ij * GS_j)
# - r(g_i) = sign(cor(x_i, Y)) * q_i

get_wgcna_network_rank <- function(mat_data, Y_vec, alpha = 6, gamma = 1) {
  # mat_data: feature sulle righe (G), campioni sulle colonne (N)
  # Y_vec: vettore numerico del tratto clinico (es. 0/1 o 1/2)
  y_num <- as.numeric(as.factor(Y_vec))
  stopifnot(ncol(mat_data) == length(y_num))
  
  # Rimuove eventuali feature a varianza zero per evitare NA nelle correlazioni
  vars <- apply(mat_data, 1, var, na.rm = TRUE)
  mat_clean <- mat_data[!is.na(vars) & vars > 0, , drop = FALSE]
  
  # 1. Correlazione clinica e Gene Significance (GS_i)
  cor_y <- suppressWarnings(
    apply(mat_clean, 1, function(x) cor(x, y_num, method = "pearson", use = "complete.obs"))
  )
  cor_y[is.na(cor_y)] <- 0
  GS <- abs(cor_y)^gamma
  
  # 2. Matrice di adiacenza pesata A = [a_ij] con soft-thresholding alpha
  # Usa adjacency() di WGCNA su t(mat_clean) (campioni sulle righe, geni sulle colonne)
  A <- adjacency(t(mat_clean), type = "unsigned", power = alpha, corFnc = "cor")
  diag(A) <- 0  # Impone j != i nella sommatoria del network screening
  
  # 3. Network Screening Score: q_i = GS_i * sum_{j != i} (a_ij * GS_j)
  # Il prodotto matrice-vettore (A %*% GS) calcola esattamente sum_{j != i} a_ij * GS_j
  neighbor_gs <- as.numeric(A %*% GS)
  q_score <- GS * neighbor_gs
  
  # 4. Ranking orientato: r(g_i) = sign(cor(x_i, Y)) * q_i
  r_vec <- sign(cor_y) * q_score
  names(r_vec) <- rownames(mat_clean)
  
  # Jitter infinitesimo per evitare pareggi esatti in fgsea
  set.seed(123)
  r_vec <- r_vec + runif(length(r_vec), -1e-12, 1e-12)
  
  return(sort(r_vec, decreasing = TRUE))
}

# ==============================================================================
# 2. STIMA AUTOMATICA DEL SOFT-THRESHOLDING ALPHA E CALCOLO RANKING SINGOLI
# ==============================================================================
# Nota: se preferisci fissare un alpha standard (es. alpha = 6 per reti unsigned),
# puoi impostare use_auto_power <- FALSE
use_auto_power <- TRUE
powers <- c(1:10, seq(12, 20, 2))

stima_alpha <- function(mat_data, default_alpha = 6) {
  if (!use_auto_power) return(default_alpha)
  sft <- suppressMessages(
    pickSoftThreshold(t(mat_data), powerVector = powers, networkType = "unsigned", verbose = 0)
  )
  if (!is.na(sft$powerEstimate)) return(sft$powerEstimate) else return(default_alpha)
}

alpha_t2D <- stima_alpha(trascr_filtered2D_scaled)
alpha_t3D <- stima_alpha(trascr_filtered3D_scaled)
alpha_p2D <- stima_alpha(prot2D_scaled)
alpha_p3D <- stima_alpha(prot3D_scaled)

cat("Potenze soft-thresholding (alpha) selezionate ->",
    "Trascr2D:", alpha_t2D, "| Trascr3D:", alpha_t3D,
    "| Prot2D:", alpha_p2D, "| Prot3D:", alpha_p3D, "\n")

# Calcolo dei ranking WGCNA sui 4 blocchi singoli
rank_wgcna_trascr2D <- get_wgcna_network_rank(trascr_filtered2D_scaled, Y[1:14], alpha = alpha_t2D, gamma = 1)
rank_wgcna_trascr3D <- get_wgcna_network_rank(trascr_filtered3D_scaled, Y[1:14], alpha = alpha_t3D, gamma = 1)
rank_wgcna_prot2D   <- get_wgcna_network_rank(prot2D_scaled,            Y[1:14], alpha = alpha_p2D, gamma = 1)
rank_wgcna_prot3D   <- get_wgcna_network_rank(prot3D_scaled,            Y[1:14], alpha = alpha_p3D, gamma = 1)

# Esecuzione GSEA su ciascun blocco singolo
set.seed(123)
gsea_wgcna_trascr2D <- fgsea(pathways = pathways_list_ensembl, stats = rank_wgcna_trascr2D, minSize = 15, maxSize = 500)
gsea_wgcna_trascr3D <- fgsea(pathways = pathways_list_ensembl, stats = rank_wgcna_trascr3D, minSize = 15, maxSize = 500)
gsea_wgcna_prot2D   <- fgsea(pathways = pathways_list_symbol,  stats = rank_wgcna_prot2D,   minSize = 5,  maxSize = 500)
gsea_wgcna_prot3D   <- fgsea(pathways = pathways_list_symbol,  stats = rank_wgcna_prot3D,   minSize = 5,  maxSize = 500)

top_wgcna_trascr2D <- clean_gsea_results_aligned(gsea_wgcna_trascr2D)
top_wgcna_trascr3D <- clean_gsea_results_aligned(gsea_wgcna_trascr3D)
top_wgcna_prot2D   <- clean_gsea_results_aligned(gsea_wgcna_prot2D)
top_wgcna_prot3D   <- clean_gsea_results_aligned(gsea_wgcna_prot3D)

# ==============================================================================
# 3. LIVELLI DI INTEGRAZIONE WGCNA (Fissata Omica e Fissata Dimensione)
# ==============================================================================
# Per l'integrazione di rete WGCNA, i punteggi di ranking r(g_i) di ciascun blocco
# vengono Z-standardizzati intra-blocco (per renderli confrontabili su scale diverse)
# e mediati per gene/omica (approccio consensus network screening).

fondi_ranking_wgcna <- function(lista_rank, map_ensembl_to_sym = NULL) {
  df_list <- lapply(names(lista_rank), function(nm) {
    rk <- lista_rank[[nm]]
    z_rk <- as.numeric(scale(rk))
    names(z_rk) <- names(rk)
    
    # Se un blocco ha ID Ensembl e viene fornita la mappa, converte in SYMBOL
    if (!is.null(map_ensembl_to_sym) && any(grepl("^ENSG", names(z_rk)))) {
      clean_ids <- clean_ensembl(names(z_rk))
      mapped <- map_ensembl_to_sym[clean_ids]
      valid <- !is.na(mapped) & mapped != ""
      z_rk <- z_rk[valid]
      names(z_rk) <- mapped[valid]
    }
    
    data.frame(Gene = names(z_rk), Z = z_rk, stringsAsFactors = FALSE) %>%
      group_by(Gene) %>%
      summarise(Z = mean(Z), .groups = "drop")
  })
  
  df_merged <- Reduce(function(x, y) full_join(x, y, by = "Gene"), df_list)
  z_cols <- setdiff(colnames(df_merged), "Gene")
  
  df_merged <- df_merged %>%
    rowwise() %>%
    mutate(Z_consensus = mean(c_across(all_of(z_cols)), na.rm = TRUE)) %>%
    ungroup() %>%
    arrange(desc(Z_consensus))
  
  out_rk <- df_merged$Z_consensus
  names(out_rk) <- df_merged$Gene
  set.seed(123)
  out_rk <- out_rk + runif(length(out_rk), -1e-12, 1e-12)
  return(sort(out_rk, decreasing = TRUE))
}

# Mappa Ensembl -> Symbol per l'integrazione cross-omica
mappa_tx_symbol <- crea_mappa_ensembl_to_symbol(rownames(trascr_filtered2D_scaled))

# A) Integrazione Fissata Omica (2D + 3D)
rank_wgcna_integ_trascr <- fondi_ranking_wgcna(
  list(T2D = rank_wgcna_trascr2D, T3D = rank_wgcna_trascr3D)
)
rank_wgcna_integ_prot <- fondi_ranking_wgcna(
  list(P2D = rank_wgcna_prot2D, P3D = rank_wgcna_prot3D)
)

set.seed(123)
gsea_wgcna_integ_trascr <- fgsea(pathways = pathways_list_ensembl, stats = rank_wgcna_integ_trascr, minSize = 15, maxSize = 500)
gsea_wgcna_integ_prot   <- fgsea(pathways = pathways_list_symbol,  stats = rank_wgcna_integ_prot,   minSize = 5,  maxSize = 500)

top_wgcna_integ_trascr <- clean_gsea_results_aligned(gsea_wgcna_integ_trascr)
top_wgcna_integ_prot   <- clean_gsea_results_aligned(gsea_wgcna_integ_prot)

# B) Integrazione Fissata Dimensione (Trascrittomica + Proteomica)
rank_wgcna_integ_2D <- fondi_ranking_wgcna(
  list(T2D = rank_wgcna_trascr2D, P2D = rank_wgcna_prot2D),
  map_ensembl_to_sym = mappa_tx_symbol
)
rank_wgcna_integ_3D <- fondi_ranking_wgcna(
  list(T3D = rank_wgcna_trascr3D, P3D = rank_wgcna_prot3D),
  map_ensembl_to_sym = mappa_tx_symbol
)

set.seed(123)
gsea_wgcna_integ_2D <- fgsea(pathways = pathways_list_symbol, stats = rank_wgcna_integ_2D, minSize = 15, maxSize = 500)
gsea_wgcna_integ_3D <- fgsea(pathways = pathways_list_symbol, stats = rank_wgcna_integ_3D, minSize = 15, maxSize = 500)

top_wgcna_integ_2D <- clean_gsea_results_aligned(gsea_wgcna_integ_2D)
top_wgcna_integ_3D <- clean_gsea_results_aligned(gsea_wgcna_integ_3D)

# ==============================================================================
# 4. CONFRONTO TRA RANKING MULTIVARIATO (PLS-DA) E DI RETE (WGCNA)
# ==============================================================================
cat("\n======================================================================\n")
cat("          RIEPILOGO PATHWAY SIGNIFICATIVI: WGCNA vs PLS-DA             \n")
cat("======================================================================\n")
cat("1. Trascrittomica 2D -> WGCNA:", nrow(top_wgcna_trascr2D), "| PLS-DA:", nrow(top_pathways_trascr2D_from_single), "\n")
cat("2. Trascrittomica 3D -> WGCNA:", nrow(top_wgcna_trascr3D), "| PLS-DA:", nrow(top_pathways_trascr3D_from_single), "\n")
cat("3. Proteomica 2D     -> WGCNA:", nrow(top_wgcna_prot2D),   "| PLS-DA:", nrow(top_pathways_prot2D_from_single), "\n")
cat("4. Proteomica 3D     -> WGCNA:", nrow(top_wgcna_prot3D),   "| PLS-DA:", nrow(top_pathways_prot3D_from_single), "\n")
cat("----------------------------------------------------------------------\n")
cat("5. Integ. Trascr (2D+3D) -> WGCNA:", nrow(top_wgcna_integ_trascr), "\n")
cat("6. Integ. Prot   (2D+3D) -> WGCNA:", nrow(top_wgcna_integ_prot), "\n")
cat("7. Integ. Stato 2D (T+P) -> WGCNA:", nrow(top_wgcna_integ_2D), "\n")
cat("8. Integ. Stato 3D (T+P) -> WGCNA:", nrow(top_wgcna_integ_3D), "\n")
cat("======================================================================\n")

# ==============================================================================
# 5. EXPORT EXCEL RISULTATI GSEA WGCNA (usando la tua funzione scrivi_foglio_gsea)
# ==============================================================================
gruppi_wgcna_singoli <- list(
  "TRASCRITTOMICA - Network Screening WGCNA" = list(
    "Trascrittomica 2D" = top_wgcna_trascr2D,
    "Trascrittomica 3D" = top_wgcna_trascr3D
  ),
  "PROTEOMICA - Network Screening WGCNA" = list(
    "Proteomica 2D" = top_wgcna_prot2D,
    "Proteomica 3D" = top_wgcna_prot3D
  )
)

gruppi_wgcna_integrati <- list(
  "INTEGRAZIONE INTRA-OMICA (2D + 3D)" = list(
    "Trascrittomica (2D+3D)" = top_wgcna_integ_trascr,
    "Proteomica (2D+3D)"     = top_wgcna_integ_prot
  ),
  "INTEGRAZIONE CROSS-OMICA (Trascr + Prot)" = list(
    "Stato 2D (Trascr + Prot)" = top_wgcna_integ_2D,
    "Stato 3D (Trascr + Prot)" = top_wgcna_integ_3D
  )
)

wb_wgcna <- createWorkbook()
scrivi_foglio_gsea(wb_wgcna, "1_WGCNA_Singole_Omiche", gruppi_wgcna_singoli,
                   "GSEA su ranking WGCNA (Network Screening orientato) - Singole Omiche")
scrivi_foglio_gsea(wb_wgcna, "2_WGCNA_Integrazione",   gruppi_wgcna_integrati,
                   "GSEA su ranking WGCNA (Consensus Network Screening) - Studi di Integrazione")

percorso_wgcna <- file.path(cartella_output, "Risultati_GSEA_WGCNA.xlsx")
saveWorkbook(wb_wgcna, file = percorso_wgcna, overwrite = TRUE)
cat("[✓] File Excel WGCNA salvato in:", percorso_wgcna, "\n")


# Cartella dove salvare i grafici GSEA (assicurati che esista o creala)
dir_plot_gsea <- "/Users/ceciliabergamini/Desktop/TESI/script/plot_GSEA"
if (!dir.exists(dir_plot_gsea)) dir.create(dir_plot_gsea, recursive = TRUE)

# ==============================================================================
# A) SALVATAGGIO FOGLIO INTERSEZIONI CORE NEL FILE EXCEL PLS-DA (wb1)
# ==============================================================================
# Aggiunge il 4° foglio con i 175 pathway core proteomici e i 7 trascrittomici
pathways_core_plsda_wgcna <- intersect(pathways_comuni_prot, top_wgcna_integ_prot$pathway)

cat("\n[✓] Pathway del Core Proteomico PLS-DA (175) confermati anche da WGCNA:",
    length(pathways_core_plsda_wgcna), "\n")

df_intersezioni_export <- bind_rows(
  data.frame(Categoria = "Core Globale Proteomica (6 studi PLS-DA)", 
             Conteggio = length(pathways_comuni_prot), 
             pathway = pathways_comuni_prot, stringsAsFactors = FALSE),
  data.frame(Categoria = "Core Globale Trascrittomica (6 studi PLS-DA)", 
             Conteggio = length(pathways_comuni_trascr), 
             pathway = pathways_comuni_trascr, stringsAsFactors = FALSE),
  data.frame(Categoria = "Super-Core Proteomico Validato (PLS-DA + WGCNA)", 
             Conteggio = length(pathways_core_plsda_wgcna), 
             pathway = pathways_core_plsda_wgcna, stringsAsFactors = FALSE)
)

addWorksheet(wb1, "4_Intersezioni_Core")
writeData(wb1, sheet = "4_Intersezioni_Core", x = df_intersezioni_export)
setColWidths(wb1, sheet = "4_Intersezioni_Core", cols = 1:3, widths = c(45, 12, 65))
saveWorkbook(wb1, file = percorso_file, overwrite = TRUE)


# ==============================================================================
# B) PLOT 1: CONFRONTO QUANTITATIVO TRA METODI (PLS-DA vs WGCNA)
# ==============================================================================
df_confronto_metodi <- data.frame(
  Studio = rep(c(
    "Trascrittomica 2D (Singola)", "Trascrittomica 3D (Singola)",
    "Proteomica 2D (Singola)",     "Proteomica 3D (Singola)",
    "Integ. Trascrittomica (2D+3D)", "Integ. Proteomica (2D+3D)",
    "Integ. Stato 2D (T+P)",       "Integ. Stato 3D (T+P)"
  ), times = 2),
  Metodo = rep(c("Loadings PLS-DA", "Network Screening WGCNA"), each = 8),
  Pathway_Significativi = c(
    # Valori PLS-DA
    nrow(top_pathways_trascr2D_from_single), nrow(top_pathways_trascr3D_from_single),
    nrow(top_pathways_prot2D_from_single),   nrow(top_pathways_prot3D_from_single),
    nrow(top_pathways_2D_from_trascr) + nrow(top_pathways_3D_from_trascr),
    nrow(top_pathways_2D_from_prot)   + nrow(top_pathways_3D_from_prot),
    nrow(top_pathways_trascr_from_2D) + nrow(top_pathways_prot_from_2D),
    nrow(top_pathways_trascr_from_3D) + nrow(top_pathways_prot_from_3D),
    # Valori WGCNA
    nrow(top_wgcna_trascr2D), nrow(top_wgcna_trascr3D),
    nrow(top_wgcna_prot2D),   nrow(top_wgcna_prot3D),
    nrow(top_wgcna_integ_trascr), nrow(top_wgcna_integ_prot),
    nrow(top_wgcna_integ_2D),     nrow(top_wgcna_integ_3D)
  )
)

df_confronto_metodi$Studio <- factor(df_confronto_metodi$Studio, levels = rev(unique(df_confronto_metodi$Studio)))

p_bar_metodi <- ggplot(df_confronto_metodi, aes(x = Studio, y = Pathway_Significativi, fill = Metodo)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68, color = "black") +
  geom_text(aes(label = Pathway_Significativi), 
            position = position_dodge(width = 0.75), 
            hjust = -0.15, size = 4.2, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("Loadings PLS-DA" = "orangered", "Network Screening WGCNA" = "darkblue")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(
    title = "Confronto Numero di Pathway Significativi (padj < 0.20, |NES| >= 1)",
    subtitle = "Ranking multivariato (PLS-DA) vs Ranking topologico di rete (WGCNA)",
    x = NULL,
    y = "Numero di Pathway Arricchiti",
    fill = "Metrica di Ranking"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title    = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    axis.text.y   = element_text(size = 12, face = "bold", color = "black"),
    axis.text.x   = element_text(size = 12, color = "black"),
    axis.title.x  = element_text(size = 13, face = "bold", margin = margin(t = 8)),
    legend.position = "bottom",
    legend.title  = element_text(size = 13, face = "bold"),
    legend.text   = element_text(size = 12)
  )

print(p_bar_metodi)
ggsave(file.path(dir_plot_gsea, "GSEA_Confronto_PLSDA_vs_WGCNA.png"), 
       plot = p_bar_metodi, width = 12, height = 7.5, dpi = 300)


# ==============================================================================
# C) PLOT 2: DOTPLOT TOP 15 PATHWAY DEL CORE PROTEOMICO VALIDATO (PLS-DA + WGCNA)
# ==============================================================================
# Prende i top 15 pathway presenti nell'intersezione globale proteomica (175)
# e ordinati per significatività nell'integrazione proteomica
df_core_top <- top_pathways_2D_from_prot %>%
  filter(pathway %in% pathways_core_plsda_wgcna) %>%
  arrange(padj, desc(abs(NES))) %>%
  slice_head(n = 15) %>%
  mutate(
    # Pulisce e accorcia leggermente i nomi lunghi per leggibilità grafica
    pathway_clean = gsub("^(GOBP|GOCC|GOMF|KEGG|REACTOME)_", "", pathway),
    pathway_clean = gsub("_", " ", pathway_clean),
    pathway_clean = stringr::str_wrap(pathway_clean, width = 42)
  )

p_core_dot <- ggplot(df_core_top, aes(x = NES, y = reorder(pathway_clean, abs(NES)))) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
  geom_point(aes(size = size, color = padj)) +
  scale_color_gradient(low = "#D73027", high = "#4575B4", name = "p.adjust (FDR)") +
  scale_size_continuous(range = c(4, 10), name = "Dimensione\nGene Set") +
  labs(
    title = "Top 15 Pathway del Core Proteomico Conservato (PLS-DA + WGCNA)",
    subtitle = "Pathway stabili tra 2D, 3D e tutti i modelli integrativi",
    x = "Normalized Enrichment Score (NES)",
    y = NULL
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title    = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    axis.text.y   = element_text(size = 11, face = "bold", color = "black"),
    axis.text.x   = element_text(size = 12, color = "black"),
    axis.title.x  = element_text(size = 13, face = "bold", margin = margin(t = 8)),
    legend.title  = element_text(size = 12, face = "bold"),
    legend.text   = element_text(size = 11)
  )

print(p_core_dot)
ggsave(file.path(dir_plot_gsea, "GSEA_Dotplot_Core_Proteomico.png"), 
       plot = p_core_dot, width = 13, height = 8.5, dpi = 300)


# ==============================================================================
# D) PLOT 3: SCATTERPLOT CONCORDANZA DIREZIONALE NES (PLS-DA vs WGCNA)
# ==============================================================================
# Confronta i valori di NES dei pathway in comune tra Integrazione Proteomica PLS-DA e WGCNA
df_nes_comp <- inner_join(
  top_pathways_2D_from_prot %>% select(pathway, NES_PLSDA = NES, padj_PLSDA = padj),
  top_wgcna_integ_prot      %>% select(pathway, NES_WGCNA = NES, padj_WGCNA = padj),
  by = "pathway"
) %>%
  mutate(
    In_Core_175 = ifelse(pathway %in% pathways_comuni_prot, "Core Proteomico (175)", "Altri Pathway Condivisi")
  )

cor_nes <- cor(df_nes_comp$NES_PLSDA, df_nes_comp$NES_WGCNA, method = "pearson")

p_scatter_nes <- ggplot(df_nes_comp, aes(x = NES_PLSDA, y = NES_WGCNA, color = In_Core_175)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
  geom_point(alpha = 0.75, size = 2.8) +
  geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dotted", linewidth = 0.9) +
  scale_color_manual(values = c("Core Proteomico (175)" = "#E31A1C", "Altri Pathway Condivisi" = "#7F8C8D")) +
  labs(
    title = "Concordanza Direzionale dell'Arricchimento: PLS-DA vs WGCNA",
    subtitle = paste0("Integrazione Proteomica (Correlazione di Pearson r = ", round(cor_nes, 2), ")"),
    x = "NES (Ranking Loadings PLS-DA)",
    y = "NES (Ranking Network Screening WGCNA)",
    color = "Categoria"
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title    = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray30"),
    axis.text     = element_text(size = 12, color = "black"),
    axis.title    = element_text(size = 13, face = "bold"),
    legend.position = "bottom",
    legend.title  = element_text(size = 13, face = "bold"),
    legend.text   = element_text(size = 12)
  )

print(p_scatter_nes)
ggsave(file.path(dir_plot_gsea, "GSEA_Concordanza_NES_PLSDA_WGCNA.png"), 
       plot = p_scatter_nes, width = 10, height = 8, dpi = 300)


# ==============================================================================
# E) PLOT 4: ENRICHMENT SCORE CURVES (Top 2 Pathway del Super-Core)
# ==============================================================================
top2_names <- df_core_top$pathway[1:2]

p_es_1_plsda <- plotEnrichment(pathways_list_symbol[[top2_names[1]]], 
                               gsea_integ_omica_prot$rankings[[1]]) +
  labs(title = paste0("PLS-DA: ", gsub("_", " ", top2_names[1])),
       x = "Rank (Loadings Z-score)", y = "Enrichment Score") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 11, face = "bold"))

p_es_1_wgcna <- plotEnrichment(pathways_list_symbol[[top2_names[1]]], 
                               rank_wgcna_integ_prot) +
  labs(title = paste0("WGCNA: ", gsub("_", " ", top2_names[1])),
       x = "Rank (Network Screening Score)", y = "Enrichment Score") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 11, face = "bold"))

p_es_2_plsda <- plotEnrichment(pathways_list_symbol[[top2_names[2]]], 
                               gsea_integ_omica_prot$rankings[[1]]) +
  labs(title = paste0("PLS-DA: ", gsub("_", " ", top2_names[2])),
       x = "Rank (Loadings Z-score)", y = "Enrichment Score") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 11, face = "bold"))

p_es_2_wgcna <- plotEnrichment(pathways_list_symbol[[top2_names[2]]], 
                               rank_wgcna_integ_prot) +
  labs(title = paste0("WGCNA: ", gsub("_", " ", top2_names[2])),
       x = "Rank (Network Screening Score)", y = "Enrichment Score") +
  theme_bw(base_size = 12) +
  theme(plot.title = element_text(size = 11, face = "bold"))

p_es_grid <- (p_es_1_plsda | p_es_1_wgcna) / (p_es_2_plsda | p_es_2_wgcna) +
  plot_annotation(
    title = "Profili di Running Enrichment Score (Top 2 Pathway Validati)",
    theme = theme(plot.title = element_text(size = 15, face = "bold", hjust = 0.5))
  )

print(p_es_grid)
ggsave(file.path(dir_plot_gsea, "GSEA_Enrichment_Curves_Top2.png"), 
       plot = p_es_grid, width = 14, height = 9, dpi = 300)



#######################################################
#### CREAZIONE TABELLE LATEX DA INSERIRE APPENDICE ####
#######################################################

library(dplyr)
library(purrr)
library(stringr)
library(readxl)
library(kableExtra)
percorso_cartella <- "/Users/ceciliabergamini/Desktop/TESI/File_per_Eugenio"

# Cartella di output per i file .tex generati
cartella_out <- file.path(percorso_cartella, "Tabelle_LaTeX")
dir.create(cartella_out, showWarnings = FALSE)

lista_file <- list.files(path = percorso_cartella, 
                         pattern = "^[^~].*\\.xlsx$", 
                         full.names = TRUE)

# 2. Funzione di estrazione completa di tutti i fogli e livelli omici
estrai_gsea_excel <- function(file_path) {
  nome_file <- basename(file_path)
  fogli <- excel_sheets(file_path)
  
  map_dfr(fogli, function(foglio) {
    raw <- read_excel(file_path, sheet = foglio, col_names = FALSE, .name_repair = "minimal")
    names(raw) <- paste0("V", seq_len(ncol(raw)))
    
    righe_header <- which(raw$V1 == "pathway")
    
    map_dfr(seq_along(righe_header), function(i) {
      idx_start <- righe_header[i] + 1
      idx_end <- if (i < length(righe_header)) righe_header[i + 1] - 4 else nrow(raw)
      
      # Estrae il nome del livello omico dalla riga sopra 'pathway'
      livello_omico <- raw$V1[righe_header[i] - 1] %>%
        str_remove("\\s*\\(\\d+\\s*pathway\\).*$") %>%
        str_trim()
      
      raw[idx_start:idx_end, 1:5] %>%
        setNames(c("pathway", "pval", "padj", "NES", "size")) %>%
        filter(!is.na(pathway), pathway != "pathway") %>%
        mutate(
          File = nome_file,
          Foglio = foglio,
          Livello_Omico = livello_omico,
          pval = as.numeric(pval),
          padj = as.numeric(padj),
          NES = as.numeric(NES),
          size = as.integer(size)
        )
    })
  })
}

# 3. Estrazione e filtraggio per significatività (padj <= 0.05)
dati_significativi <- map_dfr(lista_file, estrai_gsea_excel) %>%
  filter(padj <= 0.05) %>%
  arrange(File, Foglio, Livello_Omico, padj, desc(abs(NES)))

# Riepilogo in console di quanti pathway significativi sono stati estratti per ciascun livello
dati_significativi %>%
  count(File, Foglio, Livello_Omico, name = "N_significativi") %>%
  print(n = Inf)

# 4. Generazione file .tex ultra-leggeri (compilazione istantanea)
gruppi <- dati_significativi %>%
  dplyr::group_split(File, Foglio)

walk(gruppi, function(df_gruppo) {
  nome_file_clean <- str_remove(unique(df_gruppo$File), "\\.xlsx$")
  nome_foglio <- unique(df_gruppo$Foglio)
  
  nome_tex <- paste0(nome_file_clean, "__", nome_foglio, ".tex")
  label_tex <- tolower(str_replace_all(paste0(nome_file_clean, "_", nome_foglio), "[^A-Za-z0-9]", "_"))
  caption_tex <- str_replace_all(paste(nome_file_clean, "-", nome_foglio), "_", " ")
  
  df_ordinato <- df_gruppo %>%
    dplyr::arrange(Livello_Omico, padj, desc(abs(NES)))
  
  indice_gruppi <- table(factor(df_ordinato$Livello_Omico, levels = unique(df_ordinato$Livello_Omico)))
  
  tabella_formattata <- df_ordinato %>%
    dplyr::mutate(
      Pathway = str_replace_all(pathway, "_", " "),
      NES_fmt = sprintf("%.2f", NES),
      padj_fmt = formatC(padj, format = "e", digits = 2)
    ) %>%
    dplyr::select(Pathway, size, NES_fmt, padj_fmt)
  
  codice_latex <- tabella_formattata %>%
    kbl(
      format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      caption = paste0("Pathway significativi ($p_{\\text{adj}} \\le 0.05$) -- ", caption_tex),
      label = label_tex,
      col.names = c("Pathway", "Size", "NES", "$p_{\\text{adj}}$"),
      align = "lccc",
      escape = FALSE
    ) %>%
    kable_styling(
      latex_options = c("repeat_header"),
      repeat_header_text = "\\textit{(continua)}",
      font_size = 5.5  # <-- Carattere compatto (~scriptsize/tiny)
    ) %>%
    # Allarga al massimo la colonna Pathway per evitare gli "a capo"
    column_spec(1, width = "16.2cm") %>%
    pack_rows(index = indice_gruppi, bold = TRUE, hline_before = TRUE, hline_after = TRUE)
  
  codice_latex <- tabella_formattata %>%
    kbl(
      format = "latex",
      booktabs = TRUE,
      longtable = TRUE,
      caption = paste0("Pathway significativi ($p_{\\text{adj}} \\le 0.05$) -- ", caption_tex),
      label = label_tex,
      col.names = c("Pathway", "Size", "NES", "$p_{\\text{adj}}$"),
      align = "lccc",
      escape = FALSE
    ) %>%
    kable_styling(
      latex_options = c("repeat_header"),
      repeat_header_text = "\\textit{(continua)}",
      font_size = 5.5
    ) %>%
    column_spec(1, width = "13.2cm") %>% # Larghezza calibrata per non uscire dal foglio
    pack_rows(index = indice_gruppi, bold = TRUE, hline_before = TRUE, hline_after = TRUE)
  
  codice_latex <- codice_latex %>%
    str_replace("\\\\begin\\{longtable\\}\\[t\\]", "\\\\begin{longtable}") %>%
    paste0(
      "\\begingroup\n",
      "\\setlength{\\LTleft}{-1.5cm}\n",         # Sposta a sinistra di 1.5 cm
      "\\setlength{\\LTright}{-1.5cm}\n",        # Bilancia a destra di 1.5 cm (centrata!)
      "\\setlength{\\tabcolsep}{2.5pt}\n",
      "\\renewcommand{\\arraystretch}{0.75}\n",
      "\\setcounter{LTchunksize}{500}\n",
      .,
      "\n\\endgroup"
    )
  
  writeLines(codice_latex, file.path(cartella_out, nome_tex))
})

message("File .tex compatti generati in: ", cartella_out)


##################################################
### ESTRAZIONE RISULTATI SPLS-DA PER APPENDICE ###
##################################################
# Funzione helper per estrarre tabella completa dei biomarcatori (Nome, Componente, Loading)
estrai_biomarcatori <- function(res_obj, nome_livello) {
  mod <- res_obj$model_ber  # Coincide con model_mcc quando list.keepX ha un solo valore
  ncomp <- mod$ncomp
  
  purrr::map_dfr(1:ncomp, function(comp) {
    sv <- mixOmics::selectVar(mod, comp = comp)
    data.frame(
      Livello_Omico = nome_livello,
      Componente = paste0("Comp ", comp),
      Biomarcatore = sv$name,
      Loading = sv$value$value.var,
      stringsAsFactors = FALSE
    )
  })
}

# Estrazione per Proteomica 2D e Proteomica 3D
bio_prot2D <- estrai_biomarcatori(res_prot2D_finali, "Proteomica 2D")
bio_prot3D <- estrai_biomarcatori(res_prot3D_finali, "Proteomica 3D")

# Controllo rapido in console
print(head(bio_prot2D))
print(head(bio_prot3D))

comuni <- intersect(bio_prot2D$Biomarcatore, bio_prot3D$Biomarcatore)
cat("Biomarcatori in comune tra Prot 2D e Prot 3D:", length(comuni), "\n")
print(comuni)


# Unisce i due dataset
tutti_biomarcatori <- dplyr::bind_rows(bio_prot2D, bio_prot3D) %>%
  dplyr::mutate(
    # Escapa eventuali underscore '_' nei nomi delle proteine/geni per evitare errori LaTeX
    Biomarcatore_tex = str_replace_all(Biomarcatore, "_", "\\\\_"),
    Loading_fmt = sprintf("%.4f", Loading)
  )

# Calcola gli indici per dividere Proteomica 2D e Proteomica 3D
indice_bio <- table(factor(tutti_biomarcatori$Livello_Omico, levels = unique(tutti_biomarcatori$Livello_Omico)))

codice_bio_tex <- tutti_biomarcatori %>%
  dplyr::select(Biomarcatore_tex, Componente, Loading_fmt) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    caption = "Biomarcatori selezionati dal modello sPLS-DA per Proteomica 2D ($\\text{keepX} = 97$) e Proteomica 3D ($\\text{keepX} = 146$).",
    label = "biomarcatori_proteomica",
    col.names = c("Biomarcatore", "Componente", "Loading"),
    align = "lcc",
    escape = FALSE
  ) %>%
  kable_styling(
    latex_options = c("repeat_header"),
    repeat_header_text = "\\textit{(continua)}",
    font_size = 7
  ) %>%
  column_spec(1, width = "7.5cm") %>%
  pack_rows(index = indice_bio, bold = TRUE, hline_before = TRUE, hline_after = TRUE) %>%
  str_replace("\\\\begin\\{longtable\\}\\[t\\]", "\\\\begin{longtable}") %>%
  paste0(
    "\\begingroup\n",
    "\\renewcommand{\\arraystretch}{0.8}\n",
    .,
    "\n\\endgroup"
  )


# 1. Definisci il percorso e crea la cartella se non esiste
percorso_cartella <- "/Users/ceciliabergamini/Desktop/Tesi_latex"
cartella_out <- file.path(percorso_cartella, "Tabelle_LaTeX")
dir.create(cartella_out, recursive = TRUE, showWarnings = FALSE)

# 2. Salva il file .tex
writeLines(codice_bio_tex, file.path(cartella_out, "Biomarcatori_Proteomica_2D_3D.tex"))



###########################################################
## Salvataggio risultati sPLS-DA integrazione proteomica ##
###########################################################
# 1. Funzione di estrazione per modelli multi-blocco (block.splsda)
estrai_biomarcatori_multiblock <- function(risultati_obj, mappa_nomi_blocchi) {
  mod <- risultati_obj$model
  ncomp <- mod$ncomp[1]
  nomi_blocchi <- names(mappa_nomi_blocchi)
  
  map_dfr(nomi_blocchi, function(blocco) {
    etichetta <- mappa_nomi_blocchi[[blocco]]
    
    map_dfr(1:ncomp, function(comp) {
      sv <- mixOmics::selectVar(mod, block = blocco, comp = comp)
      df_sub <- sv[[blocco]]
      
      data.frame(
        Blocco = etichetta,
        Componente = paste0("Comp ", comp),
        Biomarcatore = df_sub$name,
        Loading = df_sub$value$value.var,
        stringsAsFactors = FALSE
      )
    })
  })
}

# 2. Estrazione da risultati_Prot
# Nota: uso la stessa mappatura del tuo report ("transcriptomics" = "Proteomica 2D", "proteomics" = "Proteomica 3D")
mappa_prot <- c(
  "transcriptomics" = "Proteomica 2D",
  "proteomics"      = "Proteomica 3D"
)

bio_multi_prot <- estrai_biomarcatori_multiblock(risultati_Prot, mappa_prot)

# Controllo rapido in console
bio_multi_prot %>%
  dplyr::count(Blocco, Componente, name = "N_Biomarcatori") %>%
  print()

# Biomarcatori in comune tra i due blocchi nel modello integrato
comuni_multi <- intersect(
  bio_multi_prot$Biomarcatore[bio_multi_prot$Blocco == "Proteomica 2D"],
  bio_multi_prot$Biomarcatore[bio_multi_prot$Blocco == "Proteomica 3D"]
)
cat("Biomarcatori in comune tra Prot 2D e Prot 3D (integrazione):", length(comuni_multi), "\n")
print(comuni_multi)

# 3. Generazione tabella LaTeX per l'Appendice
percorso_cartella <- "/Users/ceciliabergamini/Desktop/Tesi_latexo"
cartella_out <- file.path(percorso_cartella, "Tabelle_LaTeX")
dir.create(cartella_out, recursive = TRUE, showWarnings = FALSE)

tab_bio_multi <- bio_multi_prot %>%
  dplyr::mutate(
    Biomarcatore_tex = str_replace_all(Biomarcatore, "_", "\\\\_"),
    Loading_fmt = sprintf("%.4f", Loading)
  )

indice_multi <- table(factor(tab_bio_multi$Blocco, levels = unique(tab_bio_multi$Blocco)))

codice_multi_tex <- tab_bio_multi %>%
  dplyr::select(Biomarcatore_tex, Componente, Loading_fmt) %>%
  kbl(
    format = "latex",
    booktabs = TRUE,
    longtable = TRUE,
    caption = "Biomarcatori selezionati dal modello DIABLO (block.sPLS-DA) per l'integrazione Proteomica 2D -- Proteomica 3D.",
    label = "biomarcatori_diablo_prot",
    col.names = c("Biomarcatore", "Componente", "Loading"),
    align = "lcc",
    escape = FALSE
  ) %>%
  kable_styling(
    latex_options = c("repeat_header"),
    repeat_header_text = "\\textit{(continua)}",
    font_size = 7
  ) %>%
  column_spec(1, width = "7.5cm") %>%
  pack_rows(index = indice_multi, bold = TRUE, hline_before = TRUE, hline_after = TRUE) %>%
  str_replace("\\\\begin\\{longtable\\}\\[t\\]", "\\\\begin{longtable}") %>%
  paste0(
    "\\begingroup\n",
    "\\renewcommand{\\arraystretch}{0.8}\n",
    .,
    "\n\\endgroup"
  )

writeLines(codice_multi_tex, file.path(cartella_out, "Biomarcatori_DIABLO_Prot2D_Prot3D.tex"))
message("Tabella salvata in: ", file.path(cartella_out, "Biomarcatori_DIABLO_Prot2D_Prot3D.tex"))


