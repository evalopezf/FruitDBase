rule kallisto_report:
    input:
        expand(
            f"{KALLISTO_DIR}/{{run}}/run_info.json",
            run=ALL_RUNS
        )
    output:
        rds = f"{KALLISTO_DIR}/kallisto_report.rds"
    shell:
        """
        module load R

        Rscript scripts/kallisto_report.R \
            {input} \
            {output.rds}
        """