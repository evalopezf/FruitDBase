rule sortmerna_paired:
    input:
        r1 = f"{FASTP_DIR}/{{biosample}}_1.fastp.fastq.gz",
        r2 = f"{FASTP_DIR}/{{biosample}}_2.fastp.fastq.gz"
    output:
        r1 = temp(f"{REMOVE_RNA_DIR}/{{biosample}}_nonrRNA_1.fq.gz"),
        r2 = temp(f"{REMOVE_RNA_DIR}/{{biosample}}_nonrRNA_2.fq.gz")
    params:
        ref_18s = REF_18S,
        ref_28s = REF_28S,
        workdir_18s = "/lustre/home/cebas/emlopez/rnaSeqAnalysis/snakeMake/index/18S_workdir",
        workdir_28s = "/lustre/home/cebas/emlopez/rnaSeqAnalysis/snakeMake/index/28S_workdir",
        workdir_biosample = f"{REMOVE_RNA_DIR}/temp_sortmerna/{{biosample}}",
        prefix = f"{REMOVE_RNA_DIR}/{{biosample}}_nonrRNA"

    threads: 12
    resources:
        mem_mb = sortmerna_mem_mb

    shell:
        """
        set -e
        set +u
        module purge
        module load rama0.3
        module load Miniconda3/22.11.1-1
        source /dragofs/sw/foss/0.2/software/Miniconda3/4.9.2/etc/profile.d/conda.sh
        conda activate sortmerna
        set -u
        
        mkdir -p {params.workdir_biosample}
        mkdir -p {REMOVE_RNA_DIR}
        
    
        
        # Limpiar workdir previo
        if [ -d {params.workdir_biosample} ]; then
            echo "Limpiando workdir previo..."
            rm -rf {params.workdir_biosample}
        fi
        mkdir -p {params.workdir_biosample}
        
        # Ejecutar SortMeRNA con --out2 para generar archivos paired separados
        # Según documentación oficial: https://github.com/sortmerna/sortmerna/issues/295
        # --out2 genera: other_paired_fwd.fq.gz y other_paired_rev.fq.gz
        echo "Ejecutando sortmerna con --out2..."
        
        sortmerna --ref {params.ref_18s} \
                  --ref {params.ref_28s} \
                  --reads {input.r1} \
                  --reads {input.r2} \
                  --workdir {params.workdir_biosample} \
                  --fastx \
                  --other {params.prefix} \
                  --paired_out \
                  --out2 \
                  --threads {threads} \
                  -v
        
        EXIT_CODE=$?
        echo "==================================="
        echo "SortMeRNA completado (exit code: $EXIT_CODE)"
        echo "==================================="
        
        if [ $EXIT_CODE -ne 0 ]; then
            echo "ERROR: SortMeRNA falló con código $EXIT_CODE"
            exit $EXIT_CODE
        fi
        
        # Con --out2, SortMeRNA genera archivos con sufijo _paired_fwd y _paired_rev
        echo "Archivos generados por SortMeRNA:"
        ls -lh {REMOVE_RNA_DIR}/ | grep {wildcards.biosample} || true
        
        # Renombrar a formato estándar _1.fq.gz y _2.fq.gz
        # Probar primero el formato con _paired_ (más probable con --out2)

        if [ -f {params.prefix}_fwd.fq.gz ]; then
            echo "✓ Encontrados archivos _fwd/_rev"
            mv {params.prefix}_fwd.fq.gz {output.r1}
            mv {params.prefix}_rev.fq.gz {output.r2}
        else
            echo "ERROR: No se encontraron archivos paired generados por SortMeRNA"
            echo "Se esperaban:"
            echo "  - {params.prefix}_paired_fwd.fq.gz / {params.prefix}_paired_rev.fq.gz"
            echo "  - {params.prefix}_fwd.fq.gz / {params.prefix}_rev.fq.gz"
            echo ""
            echo "Archivos presentes:"
            ls -lh {REMOVE_RNA_DIR}/ | grep {wildcards.biosample} || true
            exit 1
        fi
        
        # Verificar outputs finales
        if [ ! -f {output.r1} ] || [ ! -f {output.r2} ]; then
            echo "ERROR: No se generaron los archivos de salida finales"
            exit 1
        fi
        
        echo "Archivos finales:"
        ls -lh {output.r1} {output.r2}
        
        # Verificación LIGERA usando test -s (más confiable que stat)
        if [ ! -s {output.r1} ] || [ ! -s {output.r2} ]; then
            echo "ERROR: Uno o ambos archivos están vacíos"
            exit 1
        fi
        
        # Mostrar tamaños en formato legible
        echo "✓ Verificación completada - Archivos válidos:"
        du -h {output.r1} {output.r2}
        
        # Limpiar workdir temporal
        echo "Limpiando workdir temporal..."
        rm -rf {params.workdir_biosample}
        
        echo "✓ Completado: {wildcards.biosample}"
        date
        echo "==================================="
        """
        



