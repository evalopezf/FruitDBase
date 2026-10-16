rule merge_biosample:
    input:
        r1 = lambda wc: [
            f"{FASTQ_RUN_DIR}/{run}_1.fastq"
            for run in BIOSAMPLE_RUNS[wc.biosample]
        ],
        r2 = lambda wc: [
            f"{FASTQ_RUN_DIR}/{run}_2.fastq"
            for run in BIOSAMPLE_RUNS[wc.biosample]
        ]

    output:
        r1 = temp(f"{FASTQ_BIOSAMPLE_DIR}/{{biosample}}_1.fastq"),
        r2 = temp(f"{FASTQ_BIOSAMPLE_DIR}/{{biosample}}_2.fastq")
    resources:
        mem_mb = 8000

    shell:
        """
        set -e

        echo "Merging BioSample {wildcards.biosample}"
        mkdir -p "{FASTQ_BIOSAMPLE_DIR}" 
        cat {input.r1} > {output.r1}
        cat {input.r2} > {output.r2}

        test -s {output.r1}
        test -s {output.r2}

        echo "? BioSample {wildcards.biosample} merged"
        """