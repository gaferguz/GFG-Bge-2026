#!/bin/bash

##############
# AMPs e IMD #
##############

# Creamos la DB de blastn para la version final no anotada del genoma (Transcritos completos = .mrna)
makeblastdb -in ../Final_TPM1_okaycull.mrna -parse_seqids -out Bge_mrna -dbtype nucl

# Building a new DB, current time: 11/14/2025 08:46:19
# New DB name:   /media/bioinformatica/SEAGATE_GFG1/9-TFM-2026/4.Homology/AMP/Bge_mrna
# New DB title:  ../Final_TPM1_okaycull.mrna
# Sequence type: Nucleotide
# Keep MBits: T
# Maximum file size: 1000000000B
# Adding sequences from FASTA; added 61058 sequences in 1.46978 seconds.

# BUSQUEDA DE AMPs (facilitados por Carlos García - de Silva et. al 2020) 
# Buscamos inicialmente con una similitud y query coverage del 90% (y luego seleccionamos manualmente los de > 97%)
blastn -query ./APMs-cds.fasta \
       -db Bge_mrna \
       -out AMP_TPM1sg.out \
       -perc_identity 90 \ # umbral de identidad
       -task blastn-short -word_size 7 \ # Comando especifico para secuencias cortas
       -evalue 0.1 \ # umbral de e-valor
       -qcov_hsp_perc 90 \ # umbral de cobertura minima de la query en el alineamiento
       -outfmt  "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"

# BUSQUEDA DE CDS DE IMD (facilitados por Carlos García - de Silva et. al 2024) 
# Buscamos inicialmente con una similitud y query coverage del 80% (y luego seleccionamos manualmente los de > 97%)
blastn -query ./CDS_IMD_path.fasta \
       -db Bge_mrna \
       -out IMD_TPM1sg.out \
       -perc_identity 90 \ # umbral de identidad
       -evalue 0.1 \
       -qcov_hsp_perc 80 \
       -outfmt  "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp"


########
# TOLL #
########

# Creamos la DB de blastp para la version final no anotada del genoma (Secuencias peptidicas = .aa)
makeblastdb -parse_seqids -in ../Final_TPM1_okaycull.aa -dbtype prot -out Bge_aa

Building a new DB, current time: 11/14/2025 08:51:10
New DB name:   /media/bioinformatica/SEAGATE_GFG1/9-TFM-2026/4.Homology/TOLL/Bge_aa
New DB title:  ../Final_TPM1_okaycull.aa
Sequence type: Protein
Keep MBits: T
Maximum file size: 1000000000B
Adding sequences from FASTA; added 61058 sequences in 0.841504 seconds.


# BUSQUEDA DE PEPTIDOS DE PERIPLANETA AMERICANA (RECUPERADOS MANUALMENTE DEL NCBI USANDO EL LISTADO 4IN) 
# En este caso no se puede usar el subcomando de identidad, solo filtramos por cobertura del 50% de la query y recuperamos manualmente los que tienen un PI >= 50 o aproximado
blastp -query ./XP_list.fa \
       -db Bge_aa \
       -out TOLL_okaycull1TPMsg.out \
       -outfmt "7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs qcovhsp" \
       -num_threads 8 \
       -evalue 0.000000001 \ # Umbral de significancia (irrelevante, acaba saliendo muy bajo siempre por el tamaño de la DB)
       -qcov_hsp_perc 50 # umbral de porcentaje de la query alineada

       
