rule fastp_paired:
    input:
        r1 = f"{FASTQ_BIOSAMPLE_DIR}/{{biosample}}_1.fastq",
        r2 = f"{FASTQ_BIOSAMPLE_DIR}/{{biosample}}_2.fastq"
    output:
        r1 = temp(f"{FASTP_DIR}/{{biosample}}_1.fastp.fastq.gz"),
        r2 = temp(f"{FASTP_DIR}/{{biosample}}_2.fastp.fastq.gz"),
        html = f"{QC_DIR}/{{biosample}}_report_fastp_pre.html",
        json = f"{QC_DIR}/{{biosample}}_report_fastp_pre.json"
    params:
        keep_fastq = KEEP_FASTQ
    threads: 6
    resources:
        mem_mb = 8000
    run:

        keep_flag = "true" if params.keep_fastq else "false"
        
        shell("""
        set -e
        set +u  # Deshabilitar error en variables no definidas para conda
        module purge
        module load rama0.3
        module load Miniconda3/22.11.1-1
        source $(conda info --base)/etc/profile.d/conda.sh
        

        #source /dragofs/sw/foss/0.2/software/Miniconda3/4.9.2/etc/profile.d/conda.sh
        conda activate fastp
        fastp --version
        set -u  # Reactivar después de cargar conda

        mkdir -p {FASTP_DIR}
        mkdir -p {QC_DIR}

        fastp -i {input.r1} -I {input.r2} \
            -o {output.r1} -O {output.r2} \
            --detect_adapter_for_pe \
            --qualified_quality_phred 20 \
            -w {threads} \
            --compression 1 \
            --html {output.html} \
            --json {output.json}
        
        # Los archivos temp() se eliminan automáticamente por Snakemake
        echo "✓ Fastp completado para {wildcards.biosample}"
        """)
