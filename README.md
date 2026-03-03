GFG-Bge-2026
------
Repositorio donde se implementa el código completo implementado para la memoria *"Reconstrucción del transcriptoma de Blattella germanica y analisis de expresión génica diferencial entre individuos control, axénicos y con microbiota reimplantada en distintos tejidos"*

------

Pipeline (los numeros indican los pasos de implementacion de codigo):

<p align="center">
<img width="520" height="886" alt="image" src="https://github.com/user-attachments/assets/150ba544-bd9e-4ff3-9130-8a1c9a2830e2" />
</p>

------

Carpetas (numerados con los pasos del pipeline):

**1.Saneamiento**: control de calidad de Trimmomatic + Rcorrector + fastQC

**2.Ensamblaje**: Trinity + RNASpades

**3.Pseudo-alineamiento**: cálculos de conteos/TPM (Alinear_kallisto.sh)

**4.Filtrado**: identificación de transcritos por umbral TPM  (TPM_Filtering.R)

**5.Transcriptoma**: (Generar_Consenso.sh + Evaluar_transcriptoma.sh)

**6.Anotacion**: Extracción de datos Eggnog (DB_functional.R) + listas AMP/IMD/Toll 

**7.Analisis**: Expresión diferencial/GSEA/WGCNA (DESeq2_GSEA.R + WGCNA.R)

**RESULTADOS**: 
- Excels con genes diferencialmente expresados en cada tejido
- Transcriptoma consenso anotado (.mrna)
