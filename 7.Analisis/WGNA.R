# Carga de liberias
library(WGCNA)
library(DESeq2)
library(GEOquery)
library(tidyverse)
library(CorLevelPlot)
library(gridExtra)
library(tximport)


#############################################################################################
# Listamos los archivos y generamos el path, en este caso la tabla se ha hecho en excel
x <- list.files("./")
files <- file.path("./", x, "abundance.h5")
files <- files[grepl("output_", files)]
samples <- read.table(file.path("./samples.txt"), header=TRUE)
rownames(samples) <- samples$sample

# Quitamos el grupo G2R de este analisis
samples <- samples[!grepl("G2r",samples$group),]
files <- files[!grepl("G2r", files)]
samples$path <- files # Añadimos la columna a modo comprobatorio
# Vamos a renombrar los grupos
samples$group <- gsub("^GFaFeces$","Reinfected", samples$group)
samples$group <- gsub("^GFa$","Axenic", samples$group)
samples$group <- gsub("^Ctrl$","Untreated", samples$group)

samples$condition <- gsub("GFaFeces","Reinfected", samples$condition)
samples$condition <- gsub("GFa","Axenic", samples$condition)
samples$condition <- gsub("Ctrl","Untreated", samples$condition)

txi.kallisto <- tximport(samples$path, type = "kallisto", txOut = TRUE) # Esta da error

## PARTE DE AGREGACION DE COUNTS DE TRANSCRITOS EN GENES
tx2gene <- data.frame(names(txi.kallisto$abundance[,1]), gsub("t[0-9].*","",(names(txi.kallisto$abundance[,1]))))
colnames(tx2gene) <- c("TxID","geneID")
tx2gene$geneID <- gsub("_sm","_smt3", tx2gene$geneID) #Hacemos esta correccion para la excepcion de smt3
tx2gene$geneID <- gsub("-$","", tx2gene$geneID) #Hacemos esta correccion para la excepcion de smt3
tx2gene$geneID <- gsub("_$","", tx2gene$geneID) #Hacemos esta correccion para la excepcion de smt3

txi.kallisto <- tximport(samples$path, type = "kallisto", tx2gene = tx2gene, ignoreTxVersion = TRUE)

# Objeto DESeq2
dds <- DESeqDataSetFromTximport(txi.kallisto, colData = samples, design = ~ condition)
dds <- DESeq(dds)
data <- assay(dds)

#############################################################################################
# Deteccion de outliers
gsg <- WGCNA::goodSamplesGenes(t(data))
summary(gsg)
gsg$allOK

# Eliminar los outliers
data <- data[gsg$goodGenes == TRUE,]

# Representar muestras outliers con clustering jerarquico
htree <- hclust(dist(t(data)), method = "average")

# Normalización de datos
# create a deseq2 dataset

# exclude outlier samples
colData <- samples

# making the rownames and column names identical
all(rownames(colData) %in% colnames(data.subset)) # TRUE
all(rownames(colData) == colnames(data.subset)) # TRUE


# Creacion del objeto DESeq2 definitivo (despues de descartar outliers)
dds <- DESeqDataSetFromMatrix(countData = data.subset,
                              colData = colData,
                              design = ~ condition) # not spcifying model



# Transformacion para lidiar con los 0 y la composicionalidad
dds_norm <- vst(dds) 

# Extraccion de conteos normalizados
norm.counts <- assay(dds_norm) %>% t() # Transponer

#############################################################################################
# Contruccion de la red
# Se crea vector para utilizar a la hora establecer un hombral de "thresholding"
power <- c(c(1:10), seq(from = 12, to = 50, by = 2))

# Funcion de topologia de la red
sft <- pickSoftThreshold(norm.counts,
                  powerVector = power,
                  networkType = "signed",
                  verbose = 5)


sft.data <- sft$fitIndices

# Visualizacion del poder de "thresholding"
a1 <- ggplot(sft.data, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

a2 <- ggplot(sft.data, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

grid.arrange(a1, a2, nrow = 2) 

# converion de los conteos a matriz numerica
norm.counts[] <- sapply(norm.counts, as.numeric)

soft_power <- 8 # Aunque el fit es malo, por debajo de 0.8
temp_cor <- cor
cor <- WGCNA::cor

# Aplicacion de la funcion "all-in-one" para generar una matriz de topologia
bwnet <- blockwiseModules(norm.counts,
                 maxBlockSize = 14000,
                 TOMType = "signed",
                 power = soft_power,
                 mergeCutHeight = 0.25,
                 numericLabels = FALSE,
                 randomSeed = 1234,
                 verbose = 3)

cor <- temp_cor

#############################################################################################
# Eigengenes de los modulos
module_eigengenes <- bwnet$MEs

# Dendograma de agregacion de modulos
plotDendroAndColors(bwnet$dendrograms[[1]], 
                    cbind(bwnet$unmergedColors[bwnet$blockGenes[[1]]], bwnet$colors[bwnet$blockGenes[[1]]]),
                    c("unmerged", "merged"),
                    dendroLabels = FALSE,
                    addGuide = TRUE,
                    hang= 0.03,
                    guideHang = 0.05)

#############################################################################################

# Vincular modulos a grupos de metadatos

colData$path <- NULL

# Binarizacion
traits <- colData %>%
  mutate(FB_state_bin = ifelse(grepl('FB', tissue), 1, 0)) %>%
  select(5)

traits_group <- colData %>%
  mutate(Ctrl_state_bin = ifelse(grepl('Untreated', group), 1, 0)) %>%
  select(5)

traits_tisgroup <- colData %>%
  mutate(FBCtrl_state_bin = ifelse(grepl('FB_Untreated', condition), 1, 0)) %>%
  select(5)


# Binarizacion de variables categoricas

colData$tissue <- factor(colData$tissue, levels = c("FB","FG","HG","MT","MG","SG"))
colData$group <- factor(colData$group, levels = c("Untreated","Reinfected","Axenic"))
colData$condition <- factor(colData$condition, levels = unique(colData$condition))


tissue.out <- binarizeCategoricalColumns(colData$tissue,
                           includePairwise = FALSE,
                           includeLevelVsAll = TRUE,
                           minCount = 1)

group.out <- binarizeCategoricalColumns(colData$group,
                                         includePairwise = FALSE,
                                         includeLevelVsAll = TRUE,
                                         minCount = 1)

condition.out <- binarizeCategoricalColumns(colData$condition,
                                        includePairwise = FALSE,
                                        includeLevelVsAll = TRUE,
                                        minCount = 1)

traits <- cbind(traits, tissue.out)
traits_group <- cbind(traits_group, group.out)
traits_tisgroup <- cbind(traits_tisgroup, condition.out)


# Definir numero de genes y muestras
nSamples <- nrow(norm.counts)
nGenes <- ncol(norm.counts)

# Calculo de correlaciones y extraccion de p-valores
module.trait.corr <- cor(module_eigengenes, traits, use = 'p')
module.trait.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

module.group.corr <- cor(module_eigengenes, traits_group, use = 'p')
module.group.corr.pvals <- corPvalueStudent(module.trait.corr, nSamples)

module.traits_tg.corr <- cor(module_eigengenes, traits_tisgroup, use = 'p')
module.traits_tg.corr.pvals <- corPvalueStudent(module.traits_tg.corr, nSamples)


#############################################################################################
# Generacion de modulos en heatmaps

# TISSUES + GROUP
heatmap.data.tg <- merge(module_eigengenes, traits_tisgroup, by = 'row.names')

head(heatmap.data.tg)

heatmap.data.tg <- heatmap.data.tg %>%
  column_to_rownames(var = 'Row.names')


CorLevelPlot(heatmap.data.tg,
             x = names(heatmap.data.tg)[37:53],
             y = names(heatmap.data.tg)[1:36],
             col = c("blue1", "skyblue", "white", "pink", "red"))


# Extraemos manualmente 1 a 1 todos los modulos con |R| >= 0.6 (todo positivos)
FBCtrl_pos_re <- module.gene.mapping %>%
  filter(`bwnet$colors` == 'red') %>%
  rownames()
FBCtrl_pos_pt <- module.gene.mapping %>%
  filter(`bwnet$colors` == 'paleturquoise') %>%
  rownames()
FBCtrl_pos_s3 <- module.gene.mapping %>%
  filter(`bwnet$colors` == 'sienna3') %>%
  rownames()

FBAxenic_pos_pi <- module.gene.mapping %>%
  filter(`bwnet$colors` == 'pink') %>%
  rownames()
FBAxenic_pos_re <- module.gene.mapping %>%
  filter(`bwnet$colors` == 'red') %>%
  rownames()


HGCtrl_pos_tan <- module.gene.mapping %>%
  filter(`bwnet$colors` == 'tan') %>%
  rownames()
HGCtrl_pos_sb <- module.gene.mapping %>%
  filter(`bwnet$colors` == 'saddlebrown') %>%
  rownames()

HGReinfected_pos_cy <- module.gene.mapping %>%
  filter(`bwnet$colors` == 'cyan') %>%
  rownames()

MTCtrl_pos_tu <- module.gene.mapping %>%
  filter(`bwnet$colors` == 'turquoise') %>%
  rownames()

Module_list <- list(FBCtrl_pos_re, FBCtrl_pos_pt,FBCtrl_pos_s3,
                    FBAxenic_pos_pi, FBAxenic_pos_re, 
                    HGCtrl_pos_tan, HGCtrl_pos_sb, 
                    HGReinfected_pos_cy, MTCtrl_pos_tu)

Module_name <- list("FBCtrl_pos_re", "FBCtrl_pos_pt", "FBCtrl_pos_s3",
                    "FBAxenic_pos_pi", "FBAxenic_pos_re", 
                    "HGCtrl_pos_tan", "HGCtrl_pos_sb", 
                    "HGReinfected_pos_cy", "MTCtrl_pos_tu")

# Codigo de 18.12.2025
for (i in 1:length(Module_list)) {
  modulo <- Module_list[[i]]
  write.table(modulo, file = paste0(Module_name[[i]],".csv"), row.names = FALSE, sep = ",", quote = FALSE, col.names = FALSE)
  modulo <- modulo[grepl("-", modulo)]
  if (length(modulo) >= 1) {
    print(Module_name[[i]])
    print(length(Module_list[[i]]))
    print(paste(modulo, collapse = "|"))
  }
}

# [1] "FBCtrl_pos_re"
# [1] 2222
# [1] "AMP-blattellicin-g4|AMP-drosomycin-g3|AMP-termicin-g1-g2-g3|IMD-nemo|Toll-Cact|Toll-spz6"
# [1] "FBAxenic_pos_pi"
# [1] 2118
# [1] "AMP-defensin-g1|AMP-defensin-g2|IMD-FADD"
# [1] "FBAxenic_pos_re"
# [1] 2222
# [1] "AMP-blattellicin-g4|AMP-drosomycin-g3|AMP-termicin-g1-g2-g3|IMD-nemo|Toll-Cact|Toll-spz6"
# [1] "HGReinfected_pos_cy"
# [1] 687
# [1] "AMP-drosomycin-g9"
# [1] "MTCtrl_pos_tu"
# [1] 4751
# [1] "Toll-Lwr-2"


#############################################################################################
# Analisis intramodular:


# Calculo de MM (membresia del modulo) y los p-valores asociados

# The module membership/intramodular connectivity is calculated as the correlation of the eigengene and the gene expression profile.
# This quantifies the similarity of all genes on the array to every module.

module.membership.measure <- cor(module_eigengenes, norm.counts, use = 'p')
module.membership.measure.pvals <- corPvalueStudent(module.membership.measure, nSamples)


MMcors <- as.data.frame(module.membership.measure)

# Calculate the gene significance and associated p-values

gene.signf.corr <- cor(norm.counts, traits$data.MT.vs.all, use = 'p')
gene.signf.corr.pvals <- corPvalueStudent(gene.signf.corr, nSamples)

gene.signf.corr_filt <- as.data.frame(gene.signf.corr)
gene.signf.corr_filt <-  as.character(rownames(gene.signf.corr_filt %>% dplyr::filter(.[[1]] > 0.6)))


#############################################################################################

gene.signf.corr.pvals %>%
  as.data.frame() %>%
  arrange(V1) %>%
  head(100)

# Elegir gen "hub" con mas conectividad de cada modulo
colores <- bwnet$colors
Hubs <- chooseTopHubInEachModule(
  norm.counts,
  colorh = colores,
  # omitColors = "grey",
  power = 9,
  type = "signed")

TopH <- topHubs(norm.counts, colorh = colores, power = 9, omitColors = NA)

gene.signf.corr_filt <- as.data.frame(gene.signf.corr)
gene.signf.corr_filt <-  as.character(rownames(gene.signf.corr_filt %>% dplyr::filter(.[[1]] > 0.6)))

#############################################################################################
# ADAPTACION DE FUNCION PARA EXTRAER TODOS LOS VALORES DE CONECTIVIDAD DE UN MODULO

# Using the gene significance you can identify genes that have a high significance for trait of interest
# Using the module membership measures you can identify genes with high module membership in interesting modules.

# the grey module is omitted
topHubs <- function (datExpr, colorh, omitColors = "grey", power = 2, type = "signed",
                     ...)
{
  # modified from chooseTopHubInEachModule, but return the table of all genes connectivity
  isIndex = FALSE
  modules = names(table(colorh))
  if (!is.na(omitColors)[1])
    modules = modules[!is.element(modules, omitColors)]
  if (is.null(colnames(datExpr))) {
    colnames(datExpr) = 1:dim(datExpr)[2]
    isIndex = TRUE
  }

  connectivity_table <- data.frame(matrix(ncol = 3)) %>% setNames(c('gene', 'connectivity_rowSums_adj', 'module'))
  hubs = rep(NA, length(modules))
  names(hubs) = modules
  for (m in modules) {
    adj = adjacency(datExpr[, colorh == m], power = power,
                    type = type, ...)

    hub = which.max(rowSums(adj))

    hubs[m] = colnames(adj)[hub]

    sorted_genes <- rowSums(adj) %>% sort(decreasing = T) %>% as.data.frame()  %>%
      tibble::rownames_to_column() %>% setNames(c('gene', 'connectivity_rowSums_adj')) %>% mutate(module = m)
    connectivity_table <- connectivity_table %>% rbind(sorted_genes)



  }
  if (isIndex) {
    hubs = as.numeric(hubs)
    names(hubs) = modules
  }
  return(connectivity_table %>% na.omit)
}
