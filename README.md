GFG-Bge-2026
------
Repositorio donde se recoge todo el código implementado para la memoria *"Reconstrucción del transcriptoma de Blattella germanica y análisis de expresión génica diferencial entre individuos control, axénicos y con microbiota reimplantada en distintos tejidos"*

------

Pipeline (los numeros indican los pasos de implementacion de codigo):

<p align="center">
<img width="530" height="892" alt="image" src="https://github.com/user-attachments/assets/6a0b2918-8fe2-4b79-adfd-11f97dac187b" />
</p>

------

DIRECTORIOS (numerados según el paso del pipeline):

**1.Saneamiento**: control de calidad de Trimmomatic + Rcorrector + fastQC

**2.Ensamblaje**: Trinity + RNASpades

**3.Pseudo-alineamiento**: cálculos de conteos/TPM (Alinear_kallisto.sh)

**4.Filtrado**: identificación de transcritos por umbral TPM  (TPM_Filtering.R)

**5.Transcriptoma**: (Generar_Consenso.sh + Evaluar_transcriptoma.sh)

**6.Anotacion**: Extracción de datos EggNOG (DB_functional.R) + listas AMP/IMD/Toll 

**7.Analisis**: Expresión diferencial/GSEA/WGCNA (DESeq2_GSEA.R + WGCNA.R)
