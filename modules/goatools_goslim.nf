process GOATOOLS_GOSLIM {
  tag "${sample_id}"
  publishDir "results/go/goslim/${sample_id}", mode: 'copy'

  conda "go_slim.yaml"

  input:
  tuple val(sample_id), path(target_go), path(background_go)

  output:
  path("${sample_id}_*"), emit: results
  path("${sample_id}_goslim_summary.txt"), emit: summary

  script:
  """
  #!/bin/bash
  set -euo pipefail

  # 1. Download GO slim (swap this link to test other slims)
  wget -q https://current.geneontology.org/ontology/subsets/goslim_generic.obo -O plant_goslim.obo

  # 2. Download full GO ontology
  wget -q http://purl.obolibrary.org/obo/go/go-basic.obo -O go-basic.obo

  # 3. Run GOATOOLS goslim analysis + p-values per category
  #    (base functionality preserved; we only add stats columns and later filter P<0.05 for plots)
  python3 << 'EOF'
import math
import pandas as pd
from collections import Counter
from goatools.base import get_godag

# Load GO slim and full ontology
godag = get_godag("go-basic.obo", optional_attrs={'relationship'})
goslim = get_godag("plant_goslim.obo", optional_attrs={'relationship'})

# Define root terms to exclude
EXCLUDE_ROOTS = {'GO:0009987', 'GO:0005622', 'GO:0008152', 'GO:0005737'}

# Read GO annotations (expects: gene<TAB>GO:ID)
def read_go_annotations(filename):
    annotations = {}
    with open(filename, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\\t')
                if len(parts) == 2:
                    gene, go_term = parts
                    if go_term in EXCLUDE_ROOTS:
                        continue
                    annotations.setdefault(gene, []).append(go_term)
    return annotations

target_annotations = read_go_annotations("$target_go")
background_annotations = read_go_annotations("$background_go")

print(f"Target genes: {len(target_annotations)}")
print(f"Background genes: {len(background_annotations)}")

# Map specific GO terms to GO slim categories
def map_to_goslim(go_terms, goslim_dag):
    categories = set()
    for go_id in go_terms:
        if go_id in godag:
            go_term = godag[go_id]
            if go_id in goslim_dag and go_id not in EXCLUDE_ROOTS:
                categories.add(go_id)
            for ancestor in go_term.get_all_parents():
                if ancestor in goslim_dag and ancestor not in EXCLUDE_ROOTS:
                    categories.add(ancestor)
    return list(categories)

# Classify each gene
target_categories = {}
for gene, go_terms in target_annotations.items():
    cats = map_to_goslim(go_terms, goslim)
    if cats:
        target_categories[gene] = cats

background_categories = {}
for gene, go_terms in background_annotations.items():
    cats = map_to_goslim(go_terms, goslim)
    if cats:
        background_categories[gene] = cats

# Count categories (gene-level presence)
target_counts = Counter()
for cats in target_categories.values():
    for cat in cats:
        target_counts[cat] += 1

background_counts = Counter()
for cats in background_categories.values():
    for cat in cats:
        background_counts[cat] += 1

# --- Fisher's exact test helpers (pure python) ---
def log_comb(n, k):
    if k < 0 or k > n:
        return float('-inf')
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)

def hypergeom_p(a, r1, c1, n):
    return math.exp(log_comb(c1, a) + log_comb(n - c1, r1 - a) - log_comb(n, r1))

def fisher_exact_two_sided(a, b, c, d):
    r1 = a + b
    c1 = a + c
    c2 = b + d
    n = a + b + c + d

    amin = max(0, r1 - c2)
    amax = min(r1, c1)

    p_obs = hypergeom_p(a, r1, c1, n)
    p = 0.0
    for a2 in range(amin, amax + 1):
        p2 = hypergeom_p(a2, r1, c1, n)
        if p2 <= p_obs + 1e-15:
            p += p2
    return min(p, 1.0)

def odds_ratio(a, b, c, d):
    aa, bb, cc, dd = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    return (aa * dd) / (bb * cc)

def bh_fdr(pvals):
    m = len(pvals)
    order = sorted(range(m), key=lambda i: (pvals[i] if pvals[i] is not None else 1.0))
    q = [None] * m
    prev = 1.0
    for rank, i in enumerate(reversed(order), start=1):
        p = pvals[i] if pvals[i] is not None else 1.0
        val = min(prev, (p * m) / (m - rank + 1))
        q[i] = val
        prev = val
    return q

# Get GO term names
def get_term_info(go_id):
    if go_id in godag:
        term = godag[go_id]
        return {'id': go_id, 'name': term.name, 'namespace': term.namespace}
    return None

# Build results with p-values
results = []
T = len(target_annotations)
B = len(background_annotations)

all_ids = set(target_counts.keys()) | set(background_counts.keys())
for go_id in all_ids:
    if go_id in EXCLUDE_ROOTS:
        continue
    info = get_term_info(go_id)
    if not info:
        continue

    a = int(target_counts.get(go_id, 0))
    c = int(background_counts.get(go_id, 0))
    b = int(T - a)
    d = int(B - c)

    pval = fisher_exact_two_sided(a, b, c, d) if (T > 0 and B > 0) else None
    orat = odds_ratio(a, b, c, d) if (T > 0 and B > 0) else None

    results.append({
        'GO_Slim_ID': go_id,
        'GO_Slim_Term': info['name'],
        'Category': info['namespace'],
        'Target_Genes': a,
        'Background_Genes': c,
        'Target_Percent': (a / T * 100) if T else 0,
        'Background_Percent': (c / B * 100) if B else 0,
        'Odds_Ratio': orat,
        'P_Value': pval,
    })

df = pd.DataFrame(results)

# FDR correction
if not df.empty and 'P_Value' in df.columns:
    df['FDR_BH'] = bh_fdr(df['P_Value'].tolist())
else:
    df['FDR_BH'] = []

# Sort by target genes (descending) (base behavior preserved)
df = df.sort_values('Target_Genes', ascending=False)

# Save results (UNFILTERED: keeps base functionality & full table)
df.to_csv("${sample_id}_goslim_categories.tsv", sep='\\t', index=False)

# Save summary (still shows top 10 by Target_Genes; includes stats)
with open("${sample_id}_goslim_summary.txt", 'w') as f:
    f.write(f"Plant GO Slim Analysis - ${sample_id}\\n")
    f.write("="*60 + "\\n\\n")
    f.write(f"Target genes with GO slim annotations: {len(target_categories)}/{len(target_annotations)} ({(len(target_categories)/len(target_annotations)*100) if len(target_annotations) else 0:.1f}%)\\n")
    f.write(f"Background genes with GO slim annotations: {len(background_categories)}/{len(background_annotations)} ({(len(background_categories)/len(background_annotations)*100) if len(background_annotations) else 0:.1f}%)\\n")
    f.write(f"Unique GO slim categories found: {len(df)} (root terms excluded)\\n\\n")

    sig = df[(df['P_Value'].notna()) & (df['P_Value'] < 0.05)]
    f.write(f"Significant categories (p<0.05): {len(sig)}/{len(df)}\\n\\n")

    f.write("Top 10 GO Slim Categories in Target (with Fisher p-values):\\n")
    f.write("-"*60 + "\\n")
    for _, row in df.head(10).iterrows():
        p = row.get('P_Value', None)
        q = row.get('FDR_BH', None)
        orat = row.get('Odds_Ratio', None)

        p_str = "NA" if pd.isna(p) else f"{p:.2g}"
        q_str = "NA" if pd.isna(q) else f"{q:.2g}"
        or_str = "NA" if pd.isna(orat) else f"{orat:.2g}"

        f.write(f"{row['GO_Slim_ID']}: {row['GO_Slim_Term']}\\n")
        f.write(f"  Target: {int(row['Target_Genes'])} genes ({row['Target_Percent']:.1f}%)\\n")
        f.write(f"  Background: {int(row['Background_Genes'])} genes ({row['Background_Percent']:.1f}%)\\n")
        f.write(f"  OddsRatio: {or_str}  p: {p_str}  FDR: {q_str}\\n\\n")

print("GO Slim analysis complete (including p-values)")
EOF

  # 4. Plot A: extremes by BG-multiplier (x-axis = log2(Target%/Background%))
  #    FILTER: only categories with P_Value < 0.05 are plotted
  python3 << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

try:
    df = pd.read_csv("${sample_id}_goslim_categories.tsv", sep='\\t')
    if df.empty:
        raise ValueError("No categories found in TSV (empty dataframe).")

    df = df[(df['P_Value'].notna()) & (df['P_Value'] < 0.05)].copy()
    if df.empty:
        raise ValueError("No significant categories (p < 0.05) to plot.")

    eps = 0.1
    df['Fold_Change'] = (df['Target_Percent'] + eps) / (df['Background_Percent'] + eps)
    df['Log2_FC'] = np.log2(df['Fold_Change'])

    N = 10
    over_df = df.sort_values('Log2_FC', ascending=False).head(N)
    under_df = df.sort_values('Log2_FC', ascending=True).head(N)

    plot_df = pd.concat([under_df, over_df], axis=0).copy()
    plot_df = plot_df.sort_values('Log2_FC').drop_duplicates(subset=['GO_Slim_ID'], keep='first')
    plot_df = plot_df.sort_values('Log2_FC')

    fig, ax = plt.subplots(figsize=(12, 9))
    y_pos = range(len(plot_df))
    ax.barh(y_pos, plot_df['Log2_FC'])

    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(plot_df['GO_Slim_Term'].astype(str).str.wrap(45))

    ax.axvline(0, linewidth=1)
    ax.set_xlabel('BG multiplier (log2(Target% / Background%))')
    ax.set_title('Significant GO Slim Categories (p<0.05): BG Multiplier - ${sample_id}')

    tick_vals = np.array([-3, -2, -1, 0, 1, 2, 3], dtype=float)
    ax.set_xticks(tick_vals)
    ax.set_xticklabels([f"{2**v:g}×" for v in tick_vals])

    for i, row in enumerate(plot_df.itertuples(index=False)):
        p = float(row.P_Value)
        q = float(row.FDR_BH) if not pd.isna(row.FDR_BH) else np.nan
        p_str = f"{p:.2g}"
        q_str = "NA" if np.isnan(q) else f"{q:.2g}"
        label = f"T:{row.Target_Percent:.1f}%  B:{row.Background_Percent:.1f}%  ({row.Fold_Change:.2g}× bg)  p:{p_str}  FDR:{q_str}"
        ax.text(row.Log2_FC, i, "  " + label, va='center')

    plt.tight_layout()
    plt.savefig("${sample_id}_goslim_plot_multiplier.png", dpi=300, bbox_inches='tight')
    plt.close()

    print("Significant multiplier plot created successfully")

except Exception as e:
    print(f"Could not create multiplier plot: {e}")
EOF

  # 5. Plot B: volcano-style using ONLY significant categories? (usually you want ALL points)
  #    Here: show ALL categories, but visually it highlights significance by y-axis.
  python3 << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

try:
    df = pd.read_csv("${sample_id}_goslim_categories.tsv", sep='\\t')
    if df.empty:
        raise ValueError("No categories found in TSV (empty dataframe).")

    eps = 0.1
    df['Fold_Change'] = (df['Target_Percent'] + eps) / (df['Background_Percent'] + eps)
    df['Log2_FC'] = np.log2(df['Fold_Change'])

    p = df['P_Value'].astype(float).copy()
    p = p.fillna(1.0).clip(lower=1e-300)
    df['NegLog10P'] = -np.log10(p)

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(df['Log2_FC'], df['NegLog10P'], alpha=0.7)

    ax.axvline(0, linewidth=1)
    ax.set_xlabel('log2(Target% / Background%)')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('GO Slim Volcano: Enrichment vs Significance - ${sample_id}')

    plt.tight_layout()
    plt.savefig("${sample_id}_goslim_plot_volcano.png", dpi=300, bbox_inches='tight')
    plt.close()

    print("Volcano plot created successfully")

except Exception as e:
    print(f"Could not create volcano plot: {e}")
EOF

  # 6. Plot C: top categories by Target% (FILTER: p<0.05)
  python3 << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt

try:
    df = pd.read_csv("${sample_id}_goslim_categories.tsv", sep='\\t')
    if df.empty:
        raise ValueError("No categories found in TSV (empty dataframe).")

    df = df[(df['P_Value'].notna()) & (df['P_Value'] < 0.05)].copy()
    if df.empty:
        raise ValueError("No significant categories (p < 0.05) to plot.")

    top_n = 15
    top_df = df.sort_values('Target_Percent', ascending=False).head(top_n).copy()
    top_df = top_df.iloc[::-1]

    fig, ax = plt.subplots(figsize=(12, 8))
    y_pos = range(len(top_df))

    ax.barh(y_pos, top_df['Target_Percent'], label='Target')
    ax.barh(y_pos, top_df['Background_Percent'], label='Background', alpha=0.7)

    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(top_df['GO_Slim_Term'].astype(str).str.wrap(45))
    ax.set_xlabel('Percentage of Genes (%)')
    ax.set_title('Top Significant GO Slim Categories by Target % (p<0.05) - ${sample_id}')
    ax.legend()

    plt.tight_layout()
    plt.savefig("${sample_id}_goslim_plot_top_target.png", dpi=300, bbox_inches='tight')
    plt.close()

    print("Top Target% (significant only) plot created successfully")

except Exception as e:
    print(f"Could not create top-target plot: {e}")
EOF

  # Keep compatibility: write the default plot name expected downstream
  cp -f "${sample_id}_goslim_plot_multiplier.png" "${sample_id}_goslim_plot.png"
  """
}
