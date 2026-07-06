rule fastp_report:
    input:
        reports = expand(f"{QC_DIR}/{{run}}_report_fastp_pre.json", run=ALL_RUNS)
    output:
        report_fastp = f"{QC_DIR}/all_report_fastp.csv"
    shell:
        """
        set -e
        module purge
        module load rama0.3
        module load Miniconda3/22.11.1-1
        set +u
        source /dragofs/sw/foss/0.2/software/Miniconda3/4.9.2/etc/profile.d/conda.sh
        conda activate jq

        echo "Run,Total_reads_before,Total_reads_after,Reads_removed,Percent_retained,Q20_before,Q20_after,Q30_before,Q30_after,GC_content_after,DuplicationRate,Insert_size_peak" > {output.report_fastp}

        for f in {input.reports}; do
            Run=$(basename $f _report_fastp_pre.json)
            total_before=$(jq '.summary.before_filtering.total_reads' $f)
            total_after=$(jq '.summary.after_filtering.total_reads' $f)
            low_quality=$(jq '.filtering_result.low_quality_reads' $f)
            too_many_N=$(jq '.filtering_result.too_many_N_reads' $f)
            too_short=$(jq '.filtering_result.too_short_reads' $f)
            reads_removed=$((low_quality + too_many_N + too_short))
            percent_retained=$(awk -v after="$total_after" -v before="$total_before" 'BEGIN{{printf "%.2f", (after/before)*100}}')
            q20_before=$(jq '.summary.before_filtering.q20_rate' $f)
            q20_after=$(jq '.summary.after_filtering.q20_rate' $f)
            q30_before=$(jq '.summary.before_filtering.q30_rate' $f)
            q30_after=$(jq '.summary.after_filtering.q30_rate' $f)
            duplication=$(jq '.duplication.rate' $f)
            insertsize=$(jq '.insert_size.peak' $f)
            gc=$(jq '.summary.after_filtering.gc_content' $f)
            echo "$Run,$total_before,$total_after,$reads_removed,$percent_retained,$q20_before,$q20_after,$q30_before,$q30_after,$gc,$duplication,$insertsize" >> {output.report_fastp}
        done
        """