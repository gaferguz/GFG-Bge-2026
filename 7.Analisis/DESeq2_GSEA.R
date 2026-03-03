library(dplyr)
library(ggplot2)
library(data.table)
library(MatrixGenerics)
library(matrixStats)
library(tximport)
library(DESeq2)
library(btools)
library(devtools)
library(tidyverse)
library(clusterProfiler)
library(ontologyIndex)

# perform ORA
term2gene <- read_tsv("./term2gene_GO.tsv")
term2name <- read_tsv("./term2name_GO.tsv")

###############
Functional_Enrichment <- function(deseq2_table, deseq2_name, outdir, outdir_imm){
  # Obtener nombre del plot reformateado y en el orden de grupo correcto (el 2 es el baseline)
  report = deseq2_table
  title = deseq2_name
  title = str_replace(title, "DESeq2_", "")
  title = str_replace(title, "_plot", "")
  
  reportf = unique(report)
  
  ############################################################################
  # GSEA
  # This time I am not filtering the changes
  gsea_list <- reportf %>%
    # dplyr::filter(log2FoldChange >= 1 & padj <= 0.05) %>%
    dplyr::arrange(desc(log2FoldChange))
  
  gsea_input <- gsea_list %>%
    dplyr::select(log2FoldChange) %>%
    unlist() %>%
    as.vector()
  names(gsea_input) <- rownames(gsea_list)
  
  # Gene enrichment analysis using the whole dataset sorted by FC
  enrichment_gsea <- GSEA(geneList = gsea_input,
                          TERM2GENE = term2gene,
                          TERM2NAME = term2name,
                          minGSSize = 10,
                          maxGSSize = 500,
                          eps = 1e-10,
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH")
  #save the enrichment result
  write.csv(file = paste0(outdir, "/", title, "_GSEA.tsv"),
            x = enrichment_gsea@result, quote = FALSE, row.names = FALSE, sep = "\t")
  
  if (any(enrichment_gsea@result$p.adjust <= 0.05)){
    p <- ridgeplot(enrichment_gsea,
                   core_enrichment= TRUE, # Cambiado, en FALSE no existe
                   fill="p.adjust",
                   orderBy = "NES",
                   showCategory=100) +
      ggtitle("Ridge plot for GSEA")
    
    ggsave(filename = paste0(outdir, "/", title, "_GSEA_ridgeplot.pdf"),
           plot =  p,  dpi = 300, width = 21, height = 70, units = "cm")
  }
  # Filtering
  enrichment_gsea_immuno <- enrichment_gsea@result
  enrichment_gsea_immuno$core_enrichment <- gsub("/","|",enrichment_gsea_immuno$core_enrichment)

  enrichment_gsea_immuno <- enrichment_gsea_immuno[grepl("defens|immun| amp |antimicrobial|infection|toll|imd",
                                                         enrichment_gsea_immuno$Description,
                                                         ignore.case = TRUE),]
  
  enrichment_gsea_immuno <- enrichment_gsea_immuno[!grepl("t cell|b cell|lymphocyte|natural killer|mast cell|recombination|neutrophil|leukocyte|myeloid",
                                                          enrichment_gsea_immuno$Description,
                                                          ignore.case = TRUE),]
  
  enrichment_gsea_immuno <- enrichment_gsea_immuno %>% filter(p.adjust < 0.05)
  
  if (nrow(enrichment_gsea_immuno) > 0) {
    write.csv(file = paste0(outdir_imm, "/", title,"_Immuno_GSEA.tsv"), x = enrichment_gsea_immuno, quote = FALSE, row.names = FALSE, sep = "\t")
  }
}

# Listamos los archivos y generamos el path
# IMPORTANTE LAS SAMPLES ESTAN PUESTAS SEGUN VAN ORDENADAS ALFABETICAMENTE EN LA VARIABLE
# REVISAR SIEMPRE BIEN SI LOS COUNTS COINCIDEN con las abundancias individuales de cada archivo

path <- "./"
x <- list.files(path)
x <- x[grepl("output", x)]
files <- file.path(path, x, "abundance.h5")

samples <- read.table("./samples.txt", header=TRUE)
rownames(samples) <- samples$sample

# Modificacion del 18.02.2026 (quitamos el grupo G2R)
samples <- samples[!grepl("G2r",samples$group),]
files <- files[!grepl("G2r", files)]
samples$path <- files # AÃ±adimos la columna a modo comprobatorio
# Vamos a renombrar los grupos
samples$group <- gsub("^GFaFeces$","Reinfectados", samples$group)
samples$group <- gsub("^GFa$","Axenicos", samples$group)
samples$group <- gsub("^Ctrl$","Untreated", samples$group)

samples$condition <- gsub("GFaFeces","Reinfected", samples$condition)
samples$condition <- gsub("GFa","Axenic", samples$condition)
samples$condition <- gsub("Ctrl","Untreated", samples$condition)

txi.kallisto <- tximport(samples$path, type = "kallisto", txOut = TRUE) # Esta da error

## AGREGACION DE LOS CONTEOS DE TRANSCRITOS A NIVEL DE GENES
tx2gene <- data.frame(names(txi.kallisto$abundance[,1]), gsub("t[0-9].*","",(names(txi.kallisto$abundance[,1]))))
colnames(tx2gene) <- c("TxID","geneID")
tx2gene$geneID <- gsub("_sm","_smt3", tx2gene$geneID) #Hacemos esta correccion para la excepcion de smt3
tx2gene$geneID <- gsub("-$","", tx2gene$geneID) #Hacemos esta correccion para la excepcion de smt3
tx2gene$geneID <- gsub("_$","", tx2gene$geneID) #Hacemos esta correccion para la excepcion de smt3

txi.kallisto <- tximport(samples$path, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

dds <- DESeqDataSetFromTximport(txi.kallisto, colData = samples, design = ~ condition) 
dds <- DESeq(dds)
vsd <- vst(dds)
PCAtissue <- plotPCA3D(vsd, intgroup = "tissue", ntop = nrow(assay(vsd)))
PCAcond <- plotPCA3D(vsd, intgroup = "condition", ntop = nrow(assay(vsd)))
PCAgroup = plotPCA3D(vsd, intgroup = "group", ntop = nrow(assay(vsd)))

saveRDS(PCAcond, paste0(path,"/PCAcond.rds"))
saveRDS(PCAgroup, paste0(path,"/PCAgroup.rds"))
saveRDS(PCAtissue,  paste0(path,"/PCAtissue.rds"))

# Modificacion de las anotaciones de cada transcrito colapsandolas por gen
Bge_annotations <- read.csv2("../Bge_Annotations.tsv", sep = "\t")
Bge_annotations <- Bge_annotations[,c("query","Description","Preferred_name")]
Bge_annotations$GeneID <- gsub("t[0-9].*", "", Bge_annotations$query)
Bge_annotations_col <- Bge_annotations
Bge_annotations_col$Notes <- paste0(Bge_annotations_col$query," : ", 
                                    Bge_annotations_col$Description, " (", 
                                    Bge_annotations_col$Preferred_name,")")

# Quitamos anotaciones vacias
Bge_annotations_col$Notes <- gsub("\\(-\\)", "", Bge_annotations_col$Notes)

Bge_annotations_col <-  Bge_annotations_col %>% 
  group_by(GeneID) %>% 
  reframe(Annotations = paste(Notes, collapse = " | ")) %>%
  ungroup()


# Vamos a iterar por tejidos y crear un nuevo objeto deseq por cada subconjunto
# argumentamos que de los contrario se inflan los valores de logFC en base a 
# una basemean y un factor de normalizacion con grupos muy distintos, sabiendo
# que las muestras se segregan bien por tejido y que hay suficientes replicas

dir.create("./Deseq2_tables")
dir.create("./GSEA")
dir.create("./GSEA_imm")

path_deseqtables <- "./Deseq2_tables/" 
path_GSEA <- "./GSEA/"
path_GSEA_imm <- "./GSEA_imm/"

deseq2_tables <- list()
deseq2_plots <- list()

# Tenemos que iterar por tejidos
lev <- levels(as.factor(samples[,2])) # get the variables

contador = 0
for (l in lev) {
  subsamples <- samples[grepl(l, samples$tissue),]
  sf <- paste(subsamples$sample, collapse="|")
  subfiles <- files[grepl(sf, files)]
  txi.kallisto <- tximport(subsamples$path, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)
  dds <- DESeqDataSetFromTximport(txi.kallisto, colData = subsamples, design = ~ condition) 
  dds <- DESeq(dds)
  
  lev2 <- levels(as.factor(subsamples[,4])) # get the variables
  L.pairs <- combn(seq_along(lev2), 2, simplify = FALSE, FUN = function(i) lev2[i])
  # sublist <- L.pairs[grepl("_Ctrl", L.pairs)]
  
  for (j in 1:length(L.pairs)) {
    contador = contador + 1
    factor1 = L.pairs[[j]][1]
    factor2 = L.pairs[[j]][2]
    nam <- paste0("DESeq2_",factor1,"_vs_",factor2)
    plot_n <- paste0(nam,"_plot")
    # dds_df <- data.frame(matrix(ncol = 7, nrow = 0))
    # x <- c("GeneID","StudiedGroup","BaseGroup","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj")
    # colnames(dds_df) <- x
    dds_all = dds
    
    factor1vs2 <- results(dds_all, contrast=c("condition", factor1, factor2), alpha = 0.05)
    factor1vs2 <- lfcShrink(dds_all, contrast=c("condition", factor1, factor2), res=factor1vs2, type = "ashr")
    factor1vs2 = data.frame(factor1vs2)
    # Guardamos el dataframe en este punto para pasarlo a bucles inferiores
    dds_df = factor1vs2
    # Formateamos las tablas ahora
    factor1vs2 = unique(factor1vs2)
    factor1vs2 = factor1vs2 %>% filter(padj < 0.05)
    factor1vs2 = factor1vs2 %>% filter(abs(log2FoldChange) > 1) # Filtro de Silva et al 2024
    factor1vs2$BaseGroup <- rep(factor2, nrow(factor1vs2))
    factor1vs2$StudiedGroup <- rep(factor1, nrow(factor1vs2))
    factor1vs2$GeneID <- rownames(factor1vs2)
    
    factor1vs2 <- factor1vs2[, c("GeneID","StudiedGroup","BaseGroup","baseMean","log2FoldChange","padj")]
    # anotamos
    factor1vs2 <- base::merge(factor1vs2, Bge_annotations_col, by = "GeneID", all.x = TRUE)
    
    assign(nam,dds_df)
    deseq2_tables[[contador]] <- dds_df
    deseq2_plots[[contador]] <- as.character(plot_n)
    write.csv2(factor1vs2, file = paste0(path_deseqtables, nam, ".tsv"), 
               quote = FALSE, col.names = TRUE, 
               sep = "\t", row.names = FALSE)
    # Functional enrichment con funcion
    Functional_Enrichment(dds_df, nam, path_GSEA, path_GSEA_imm)
    # command <- paste0("echo ",nam," >> ./DE.txt")
    # system(command)
    # system("echo \n >> ./DE.txt")
  }
}

#####################
# VULCANOS - DESeq2 #
#####################
library(tidyverse)
library(phyloseq)
library(EnhancedVolcano)
library(ggpubr)

deseq2_plots_objects <- list()

for(i in 1:length(deseq2_plots)){ 
  # Obtener nombre del plot reformateado y en el orden de grupo correcto (el 2 es el baseline)
  report = data.frame(deseq2_tables[[i]])
  title_plot = deseq2_plots[[i]]
  title_plot = str_replace(title_plot, "DESeq2_", "")
  title_plot = str_replace(title_plot, "_vs_", "(UP) vs ")
  title_plot = str_replace(title_plot, "_plot", "(DOWN)")
  
  # Filtrdo por significancia 
  reportf = unique(report)
  reportf = report %>% filter(padj < 0.05)
  reportf = reportf %>% filter(abs(log2FoldChange) > 1) # Filtro de Silva et al 2024

  # Sacamos nombres de cabecera
  ups <- reportf %>% filter(log2FoldChange > 1)
  ups <- nrow(ups)
  down <- reportf %>% filter(log2FoldChange < -1)
  down <- nrow(down)
  
  if (nrow(reportf) > 0) {
    subtitle_plot <- paste0("DEGs: ", nrow(reportf), " - Up: ", ups, " - Down: ", down)

    plot = EnhancedVolcano(report,
                           title = title_plot,
                           lab = rownames(report),
                           selectLab = selected_tx,
                           titleLabSize = 16,
                           subtitle = subtitle_plot,
                           subtitleLabSize = 18,
                           captionLabSize = 20,
                           x = "log2FoldChange",
                           y = "padj",
                           FCcutoff = 1, # LogFoldChange mayor a 2 (referencia de Silva et al 2024)
                           labSize = 0,
                           boxedLabels = FALSE,
                           pointSize = 2.0,
                           max.overlaps = 5,
                           legendPosition = "none",
                           gridlines.major = TRUE,
                           gridlines.minor = TRUE,
                           drawConnectors = FALSE,
                           widthConnectors = 0,
                           pCutoff = 0.05) 
    
    # Assing all objects to its variable name
    assign(as.character(deseq2_plots[[i]]), plot)
    deseq2_plots_objects[[i]] <-  plot
    ggsave(filename = paste0(path_deseqtables, as.character(deseq2_plots[[i]]), ".png"), 
           plot = plot, width = 6, height = 6)
  }

}
