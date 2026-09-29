rule merge_biosample:
    input:
            r1 = lambda wc: [
                f"{FASTQ_DIR}/{run}_1.fastq.gz"
                for run in BIOSAMPLE_RUNS[wc.biosample]
            ],
            r2 = lambda wc: [
                f"{FASTQ_DIR}/{run}_2.fastq.gz"
                for run in BIOSAMPLE_RUNS[wc.biosample]
            ]

        output:
            r1 = temp(f"{FASTQ_DIR}/{{biosample}}_1.fastq.gz"),
            r2 = temp(f"{FASTQ_DIR}/{{biosample}}_2.fastq.gz")

        shell("""
            set -e

            echo "Merging BioSample {wildcards.biosample}"

            zcat {input.r1} | pigz > {output.r1}
            zcat {input.r2} | pigz > {output.r2}

            test -s {output.r1}
            test -s {output.r2}

            echo "✓ BioSample {wildcards.biosample} merged"
        """)