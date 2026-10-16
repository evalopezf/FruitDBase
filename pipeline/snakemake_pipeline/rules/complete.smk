rule pipeline_summary:
    input:
        fastp_reports = expand(
            f"{QC_DIR}/{{biosample}}_report_fastp_pre.json",
            biosample=ALL_BIOSAMPLES
        ),
        kallisto = expand(
            f"{KALLISTO_DIR}/{{biosample}}/run_info.json",
            biosample=ALL_BIOSAMPLES
        )

    output:
        f"{BASE_DIR}/pipeline_summary.txt"

    params:
        n_samples = len(ALL_BIOSAMPLES),
        qc_dir = QC_DIR,
        kallisto_dir = KALLISTO_DIR

    shell:
        """
        echo "PIPELINE SUMMARY - $(date)" > {output}
        echo "BioSamples: {params.n_samples}" >> {output}

        echo "Fastp reports: $(ls {params.qc_dir}/*_report_fastp_pre.json 2>/dev/null | wc -l)" >> {output}
        echo "Kallisto runs: $(ls {params.kallisto_dir}/*/run_info.json 2>/dev/null | wc -l)" >> {output}
        """