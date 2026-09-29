rule kallisto_quant_paired:
    input:
        r1 = f"{REMOVE_RNA_DIR}/{{biosample}}_nonrRNA_1.fq.gz",
        r2 = f"{REMOVE_RNA_DIR}/{{biosample}}_nonrRNA_2.fq.gz",
        index = KALLISTO_INDEX
    output:
        abundance = protected(f"{KALLISTO_DIR}/{{biosample}}/abundance.h5"),
        tsv = protected(f"{KALLISTO_DIR}/{{biosample}}/abundance.tsv"),
        run_info = f"{KALLISTO_DIR}/{{biosample}}/run_info.json"
    params:
        outdir = f"{KALLISTO_DIR}/{{biosample}}"

    threads: 4
    resources:
        mem_mb = 8000
    shell:
        """
        set -e
        module load rama0.4 GCC/12.3.0 OpenMPI/4.1.5 kallisto/0.51.1
        
        mkdir -p {params.outdir}
        
        kallisto quant -i {input.index} \
                       -o {params.outdir} \
                       -t {threads} \
                       {input.r1} {input.r2}
        
        # Los archivos temp() se eliminan automáticamente por Snakemake
        echo "✓ Cuantificación completada para {wildcards.biosample}"
        """

