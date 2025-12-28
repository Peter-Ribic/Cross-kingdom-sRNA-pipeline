process PLOT_SRNA_LENGTH_DISTRIBUTION {
  tag "$sample_id"
  publishDir "results/plots/length_distribution/${sample_id}", mode: 'copy'

  conda "shortstack_plot.yaml"

  input:
  tuple val(sample_id), path(length_distribution_tsv)

  output:
  path("${sample_id}_length_distribution.png"), emit: png
  path("${sample_id}_length_distribution.pdf"), emit: pdf

  script:
  """
  python - << 'PY'
  import sys
  import pandas as pd
  import matplotlib.pyplot as plt

  sample_id = "${sample_id}"
  tsv_path  = "${length_distribution_tsv}"

  df = pd.read_csv(tsv_path, sep='\\t', header=None, names=['length', 'count'])
  df['length'] = pd.to_numeric(df['length'], errors='coerce')
  df['count']  = pd.to_numeric(df['count'], errors='coerce')
  df = df.dropna().sort_values('length')

  if df.empty:
      raise SystemExit(f"[ERROR] Empty/invalid TSV: {tsv_path}")

  # -------- Display cap --------
  MAX_LEN = 50
  df = df[df['length'] <= MAX_LEN]
  # -----------------------------

  lengths = df['length'].astype(int).tolist()
  counts  = df['count'].astype(int).tolist()

  min_len = min(lengths)

  fig, ax = plt.subplots(figsize=(10, 4.5))
  ax.bar(lengths, counts)

  ax.set_xlabel("Read length (nt)")
  ax.set_ylabel("Count")
  ax.set_title(f"{sample_id}: sRNA length distribution (≤ {MAX_LEN} nt)")

  # Start plot at first visible bar
  ax.set_xlim(min_len - 0.5, MAX_LEN + 0.5)
  ax.set_xticks(lengths)

  ax.grid(axis='y', linestyle=':', linewidth=0.6, alpha=0.6)

  plt.tight_layout()

  out_png = f"{sample_id}_length_distribution.png"
  out_pdf = f"{sample_id}_length_distribution.pdf"
  fig.savefig(out_png, dpi=250)
  fig.savefig(out_pdf)
  plt.close(fig)

  print(f"[OK] Wrote {out_png} and {out_pdf}", file=sys.stderr)
  PY
  """
}
