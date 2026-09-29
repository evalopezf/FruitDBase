rule kallisto_report:
    input:
        expand(
            f"{KALLISTO_DIR}/{{biosample}}/run_info.json",
            biosample=ALL_BIOSAMPLES
        )
    output:
        rds = f"{KALLISTO_DIR}/kallisto_report_{TODAY}.rds"
    resources:
        mem_mb = 4096
    shell:
        """
        module load R

        Rscript scripts/kallisto_report.R \
            {input} \
            {output.rds}
        """