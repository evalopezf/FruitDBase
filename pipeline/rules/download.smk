rule download:
    output:
        r1 = temp(f"{FASTQ_RUN_DIR}/{{run}}_1.fastq"),
        r2 = temp(f"{FASTQ_RUN_DIR}/{{run}}_2.fastq")

    wildcard_constraints:
        run = "[^_/]+"

    threads: 6
    resources:
        mem_mb = 2048

    shell:
        """
        set -e

        mkdir -p "{SRA_DIR}" "{FASTQ_RUN_DIR}"

        echo "========================================"
        echo "Downloading {wildcards.run}"
        echo "========================================"

        "{SRA_TOOLKIT}/prefetch" \
            "{wildcards.run}" \
            --output-directory "{SRA_DIR}"

        echo "✓ SRA downloaded"

        "{SRA_TOOLKIT}/fasterq-dump" \
            "{SRA_DIR}/{wildcards.run}/{wildcards.run}.sra" \
            --split-files \
            --threads {threads} \
            --outdir "{FASTQ_RUN_DIR}"

        echo "✓ FASTQ generated"


        test -s "{FASTQ_RUN_DIR}/{wildcards.run}_1.fastq"
        test -s "{FASTQ_RUN_DIR}/{wildcards.run}_2.fastq"

        rm -rf "{SRA_DIR}/{wildcards.run}"

        echo "✓ SRA deleted"
        echo "✓ Download completed for {wildcards.run}"
        """