process PLOT_SHORTSTACK_CLUSTERS {
  tag "$sample_id"
  publishDir "results/plots/cluster_positions/${sample_id}", mode: 'copy'

  conda "shortstack_plot.yaml"

  input:
  tuple val(sample_id), path(shortstack_clusters_tsv)
  path(genome_fasta)

  output:
  path("${sample_id}_shortstack_clusters.png"), emit: png
  path("${sample_id}_shortstack_clusters.pdf"), emit: pdf

  script:
  """
  python - << 'PY'
  import re
  import sys
  import pandas as pd
  import numpy as np
  import matplotlib.pyplot as plt
  from Bio import SeqIO
  from matplotlib.patches import Rectangle

  sample_id = "${sample_id}"
  tsv_path  = "${shortstack_clusters_tsv}"
  fa_path   = "${genome_fasta}"

  # -------------------------
  # Helpers
  # -------------------------
  def natural_key(s):
      return [int(t) if t.isdigit() else t.lower()
              for t in re.split(r'(\\d+)', str(s))]

  # -------------------------
  # Genome lengths
  # -------------------------
  chrom_lengths = {rec.id: len(rec.seq) for rec in SeqIO.parse(fa_path, "fasta")}
  if not chrom_lengths:
      raise SystemExit("[ERROR] No contigs read from FASTA")

  # -------------------------
  # Read ShortStack clusters
  # -------------------------
  try:
      df = pd.read_csv(tsv_path, sep='\\t', comment='#')
      if df.shape[1] == 1:
          df = pd.read_csv(tsv_path, sep=r'\\s+', engine='python', comment='#')
  except Exception:
      df = pd.read_csv(tsv_path, sep=r'\\s+', engine='python', comment='#')

  cols = {c.lower(): c for c in df.columns}
  def pick(k):
      for v in [k, *{"chrom":["chr","chromosome","seqid"],
                     "start":["begin"],
                     "end":["stop"]}.get(k, [])]:
          if v in cols:
              return cols[v]
      return None

  c_chrom = pick("chrom")
  c_start = pick("start")
  c_end   = pick("end")
  if not (c_chrom and c_start and c_end):
      raise SystemExit("[ERROR] Need Chrom/Start/End columns")

  df[c_start] = pd.to_numeric(df[c_start], errors="coerce")
  df[c_end]   = pd.to_numeric(df[c_end], errors="coerce")
  df = df.dropna(subset=[c_chrom, c_start, c_end])
  df = df[df[c_chrom].isin(chrom_lengths)]

  s = df[c_start].astype(int).to_numpy()
  e = df[c_end].astype(int).to_numpy()
  df[c_start] = np.minimum(s, e)
  df[c_end]   = np.maximum(s, e)

  # -------------------------
  # Plot
  # -------------------------
  contigs = sorted(chrom_lengths, key=natural_key)
  n = len(contigs)
  genome_max = max(chrom_lengths[c] for c in contigs)

  fig_h = max(6.0, 0.45 * n + 2.0)
  fig_w = 16.0
  fig, ax = plt.subplots(figsize=(fig_w, fig_h))

  y_pos = {c: (n - 1 - i) for i, c in enumerate(contigs)}

  # -------- VISIBILITY CONTROLS --------
  CHR_LW        = 12.0      # chromosome line thickness
  CHR_ALPHA     = 0.35

  CLUSTER_HALF  = 0.3     # vertical half-height (HEIGHT)
  VISIBLE_BP    = 5_000   # minimum visible WIDTH in bp (KEY PARAMETER)
  # ------------------------------------

  # Chromosome center lines
  for c in contigs:
      ax.hlines(
          y_pos[c],
          0,
          chrom_lengths[c],
          linewidth=CHR_LW,
          alpha=CHR_ALPHA,
          zorder=1
      )

  # Cluster bands (width + height controlled)
  for c, sub in df.groupby(c_chrom, sort=False):
      y = y_pos[c]
      for x1, x2 in zip(sub[c_start], sub[c_end]):
          length = x2 - x1
          if length < VISIBLE_BP:
              mid = (x1 + x2) / 2
              x1 = max(0, mid - VISIBLE_BP / 2)
              x2 = min(chrom_lengths[c], mid + VISIBLE_BP / 2)

          ax.add_patch(
              Rectangle(
                  (x1, y - CLUSTER_HALF),
                  x2 - x1,
                  2 * CLUSTER_HALF,
                  facecolor="black",
                  edgecolor="none",
                  zorder=3
              )
          )

  # Axes
  ax.set_yticks([y_pos[c] for c in contigs])
  ax.set_yticklabels(contigs, fontsize=10)

  ax.set_xlim(0, genome_max * 1.01)
  ax.set_xlabel("Genomic position (bp)")
  ax.set_title(f"{sample_id}: ShortStack cluster distribution")

  ax.spines["top"].set_visible(False)
  ax.spines["right"].set_visible(False)
  ax.grid(axis="x", linestyle=":", linewidth=0.6, alpha=0.4)

  ax.text(
      0.995, 0.01,
      f"Clusters plotted: {df.shape[0]}",
      transform=ax.transAxes,
      ha="right", va="bottom",
      fontsize=10
  )

  plt.tight_layout()

  out_png = f"{sample_id}_shortstack_clusters.png"
  out_pdf = f"{sample_id}_shortstack_clusters.pdf"
  fig.savefig(out_png, dpi=300)
  fig.savefig(out_pdf)
  plt.close(fig)

  print(f"[OK] Wrote {out_png} and {out_pdf}", file=sys.stderr)
  PY
  """
}
