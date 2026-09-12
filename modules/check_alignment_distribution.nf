process CHECK_ALIGNMENT_DISTRIBUTION {
    tag "$sample_id"

    conda "conda-forge::python conda-forge::numpy conda-forge::matplotlib bioconda::samtools"
    publishDir "results/viruses/virus_alignment_dist/${sample_id}", mode: 'symlink'

    input:
    tuple val(sample_id), path(alignment_table), path(alignment_raw)

    output:
    tuple val(sample_id), path("${sample_id}_viruses.stats.txt"), path("${sample_id}_idxstats_raw.txt"), emit: results
    path("${sample_id}_viruses.percent_row.txt"), emit: percent_row
    path("${sample_id}_first_ref_depth_distribution.png"), emit: first_ref_depth_distribution

    script:
    """
    set -euo pipefail

    samtools view -b ${sample_id}_viruses.sam > ${sample_id}_viruses.bam

    samtools sort ${sample_id}_viruses.bam -o ${sample_id}_viruses.sorted.bam

    samtools index ${sample_id}_viruses.sorted.bam

    samtools idxstats ${sample_id}_viruses.sorted.bam > ${sample_id}_idxstats_raw.txt

    total=\$(awk '{sum+=\$3} END{print sum+0}' ${sample_id}_idxstats_raw.txt)

    awk -v total="\$total" 'BEGIN{OFS="\\t"; print "ref_name","mapped_reads","percent_mapped"}
         {percent=(total>0)?(\$3/total)*100:0; print \$1,\$3,percent}' \\
         ${sample_id}_idxstats_raw.txt > ${sample_id}_viruses.stats.txt

    header=\$(awk 'BEGIN{OFS="\\t"} {printf "%s\\t", \$1} END{print ""}' ${sample_id}_idxstats_raw.txt)
    percents=\$(awk -v total="\$total" 'BEGIN{OFS="\\t"} {percent=(total>0)?(\$3/total)*100:0; printf "%.2f\\t", percent} END{print ""}' ${sample_id}_idxstats_raw.txt)

    (echo -e "sample_id\\t\$header"; echo -e "${sample_id}\\t\$percents") > ${sample_id}_viruses.percent_row.txt

    first_ref=\$(awk 'NR==1{print \$1}' ${sample_id}_idxstats_raw.txt)

    export SAMPLE_ID="${sample_id}"
    export FIRST_REF="\$first_ref"

    samtools depth -a -r "\$first_ref" ${sample_id}_viruses.sorted.bam > ${sample_id}_first_ref.depth_by_pos.tsv

    python - <<'PY'
    import os
    import numpy as np
    import matplotlib.pyplot as plt

    sample_id = os.environ["SAMPLE_ID"]
    ref_name  = os.environ["FIRST_REF"]

    tsv_path = f"{sample_id}_first_ref.depth_by_pos.tsv"
    out_png  = f"{sample_id}_first_ref_depth_distribution.png"

    try:
        arr = np.loadtxt(tsv_path, usecols=(1,2))
    except Exception:
        arr = np.empty((0,2))

    if arr.size == 0:
        plt.figure(figsize=(10,4))
        plt.title(f"coverage across length\\nfirst ref: {ref_name}")
        plt.xlabel("Position (bp)")
        plt.ylabel("Depth")
        plt.text(0.5, 0.5, "No depth data", ha="center", va="center", transform=plt.gca().transAxes)
        plt.tight_layout()
        plt.savefig(out_png, dpi=200)
        raise SystemExit(0)

    arr = np.atleast_2d(arr)
    pos = arr[:,0].astype(int)
    dep = arr[:,1].astype(float)

    pct_cov = (dep >= 1).mean() * 100.0
    max_points = 20000
    step = max(1, len(pos) // max_points)
    pos_ds = pos[::step]
    dep_ds = dep[::step]

    plt.figure(figsize=(10,4))
    plt.plot(pos_ds, dep_ds, linewidth=0.8)
    plt.title(
        f"{sample_id}: coverage across length\\n"
        f"first ref: {ref_name} plotted every {step} bp"
    )
    plt.xlabel("Position (bp)")
    plt.ylabel("Depth")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    PY

    rm -f ${sample_id}_first_ref.depth_by_pos.tsv
    """
}
