GFG-Bge-2026
------
Repositorio donde se implementa el código completo implementado para la memoria *"Reconstrucción del transcriptoma de Blattella germanica y analisis de expresión génica diferencial entre individuos control, axénicos y con microbiota reimplantada en distintos tejidos"*

------
Pipeline:

<p align="center">
  <img width="495" height="703" alt="image" src="https://github.com/user-attachments/assets/d1fa2b47-3094-4cc4-9f3b-c074b60d3803" />
</p>
------

Carpetas:

**1.Saneamiento**: control de calidad de Trimmomatic + Rcorrector + fastQC

**2.Ensamblaje**: Trinity + RNASpades (Sanear_[X]_TFM.sh)

**3.Pseudo-alineamiento**: cálculos de conteos/TPM (Alinear_kallisto.sh)

**4.Filtrado**: identificación de transcritos por umbral TPM  (TPM_Filtering.R)

**5.Transcriptoma**: (Generar_Consenso.sh + Evaluar_transcriptoma.sh)

**6.Anotacion**: Extracción de datos Eggnog (DB_functional.R) + listas AMP/IMD/Toll 

**7.Analisis**: Expresión diferencial/GSEA/WGCNA (DESeq2_GSEA.R + WGCNA.R)

**RESULTADOS**: 
- Excels con genes diferencialmente expresados en cada tejido
- Transcriptoma consenso anotado (.mrna)
