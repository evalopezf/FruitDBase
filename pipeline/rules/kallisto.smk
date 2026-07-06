rule kallisto_quant_paired:
    input:
        r1 = f"{REMOVE_RNA_DIR}/{{run}}_nonrRNA_1.fq.gz",
        r2 = f"{REMOVE_RNA_DIR}/{{run}}_nonrRNA_2.fq.gz",
        index = KALLISTO_INDEX
    output:
        abundance = protected(f"{KALLISTO_DIR}/{{run}}/abundance.h5"),
        tsv = protected(f"{KALLISTO_DIR}/{{run}}/abundance.tsv"),
        run_info = f"{KALLISTO_DIR}/{{run}}/run_info.json"
    params:
        outdir = f"{KALLISTO_DIR}/{{run}}"
    wildcard_constraints:
        run = "[^_/]+"
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
        

        echo "✓ Quantification completed for {wildcards.run}"
        """

