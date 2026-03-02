GFG-Bge-2026
Repositorio donde se implementa el código completo implementado para la memoria 
*"Reconstrucción del transcriptoma de Blattella germanica y analisis de expresión génica diferencial entre individuos control, axénicos y con microbiota reimplantada en distintos tejidos"*

------
Pipeline:

<img width="495" height="703" alt="image" src="https://github.com/user-attachments/assets/d1fa2b47-3094-4cc4-9f3b-c074b60d3803" />

------

Carpetas:

**1.Saneaminto**: control de calidad de Trimmomatic + Rcorrector + fastQC

**2.Ensamblaje**: Trinity + RNASpades (Trimming/Sanear_[X]_TFM.sh)

**3.Pseudo-alineamiento**: calculos de conteos/TPM con Kallisto (Alinear_kallisto.sh)

**4.Filtrado**: identificación de transcritos por umbral TPM  (Filtrado_TPM/TPM_Filtering.R)

**5.Transcriptoma**: (Generar_Consenso.sh + Evaluar_transcriptoma.sh)

**6.Anotacion**: Extraccion de datos Eggnog (DB_functional.R) + listas AMP/IMD/Toll 

**7.Analisis**: Expresion diferencial con DESeq2 + GSEA (Analisis_DESeq2_GSEA.R) + WGCNA_TMF.R

RESULTADOS: Excels con genes diferencialmente expresados para cada tejido y comparacion
