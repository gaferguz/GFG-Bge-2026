GFG-Bge-2026
Repositorio donde se implementa el código completo implementado para la memoria 
*"Reconstruccion del transcriptoma de blattella germanica y analisis DE EXPRESION GENICA DIFERENCIAL ENTRE INDIVIDUOS CONTROL, AXENICOS Y CON MICROBIOTA REIMPLANTADA EN DISTINTOS TEJIDOS"*

Pipeline:
<img width="495" height="703" alt="image" src="https://github.com/user-attachments/assets/d1fa2b47-3094-4cc4-9f3b-c074b60d3803" />

Carpetas:

1 - Saneaminto: control de calidad de Trimmomatic + Rcorrector + fastQC

2 - Ensamblaje: Trinity + RNASpades (Trimming/Sanear_[X]_TFM.sh)

3 - Pseudo-alineamiento: calculos de conteos/TPM usando Kallisto (Alinear_kallisto.sh)

4 - Filtrado: identificación de transcritos por TPM >= 1 en 75% de 1 grupo (Filtrado_TPM/TPM_Filtering.R)

5 - Transcriptoma: (Generar_Consenso.sh + Evaluar_transcriptoma.sh)

6 - Anotacion: Extraccion de informacion Eggnog (DB_functional.R) + anotacion listas AMP/IMD/Toll 

7 - Analisis: Expresion diferencial con DESeq2 + GSEA con clusterProfiler (Analisis_DET_GSEA.R) + WGCNA_TMF.R

RESULTADOS: Excels con genes diferencialmente expresados para cada tejido y comparacion
