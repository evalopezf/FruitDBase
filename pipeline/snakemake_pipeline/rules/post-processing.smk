# Post-processing: pseudoalignment filter -> integration (BgeeCall)
# Run with:  snakemake post_processing

rule filter_by_pseudoalignment:
    input:
        metadata = METADATA_FILE,
        run_info = expand(
            f"{KALLISTO_DIR}/{{biosample}}/run_info.json",
            biosample=ALL_BIOSAMPLES
        )
    output:
        metadata_filtered = f"{POST_DIR}/metadata_filtered.tsv",
        removed = f"{POST_DIR}/removed_samples.tsv"
    params:
        kallisto_dir = KALLISTO_DIR,
        min_pct = MIN_PSEUDOALIGN_PCT
    resources:
        mem_mb = 2000
    script:
        "../scripts/filter_by_pseudoalignment.py"


rule integration:
    input:
        metadata = f"{POST_DIR}/metadata_filtered.tsv",
        fasta = TRANSCRIPTOME_FASTA,
        gtf = ANNOTATION_GTF
    output:
        tpm_gene = f"{INTEGRATION_DIR}/tpm_gene.rds",
        counts_gene = f"{INTEGRATION_DIR}/counts_gene.rds",
        calls = f"{INTEGRATION_DIR}/presence_absence_df.rds",
        presence_tissue = f"{INTEGRATION_DIR}/presence_by_tissue.rds",
        presence_dev = f"{INTEGRATION_DIR}/presence_by_dev.rds",
        score_tissue = f"{INTEGRATION_DIR}/ExpressionScoreByTissue.rds",
        score_dev = f"{INTEGRATION_DIR}/ExpressionScoreByDev.rds"
    params:
        script = f"{workflow.basedir}/scripts/integration_pipeline_v2.R",
        kallisto_dir = KALLISTO_DIR,
        species_id = SPECIES_ID,
        outdir = INTEGRATION_DIR,
        work_dir = f"{POST_DIR}/tmp_bgeecall"
    resources:
        mem_mb = 32000
    shell:
        """
        set -e
        module load R

        Rscript {params.script} \
            --metadata {input.metadata} \
            --kallisto-dir {params.kallisto_dir} \
            --fasta {input.fasta} \
            --gtf {input.gtf} \
            --species-id {params.species_id} \
            --output-dir {params.outdir} \
            --work-dir {params.work_dir}
        """


rule post_processing:
    input:
        rules.integration.output
