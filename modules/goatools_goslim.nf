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

  wget -q https://current.geneontology.org/ontology/subsets/goslim_generic.obo -O plant_goslim.obo
  wget -q http://purl.obolibrary.org/obo/go/go-basic.obo -O go-basic.obo

  # NOTE:
  # Without an external gene/protein annotation source (UniProt/Ensembl/TAIR GFF/GTF, etc.)
  # we cannot reliably attach true "protein name/role" beyond GO-term-derived descriptions.
  # The protein report produced includes GO-derived details.

  python3 << 'EOF'
import math
import pandas as pd
from collections import Counter, defaultdict
from goatools.base import get_godag

FDR_LIMIT = 0.1
TOP_PLOT_N = 15
EPS_FC = 0.1

EXCLUDE_ROOTS = {'GO:0009987', 'GO:0005622', 'GO:0008152', 'GO:0005737'}
PLOT_EXCLUDE_GOID = "GO:0048856"  # anatomical structure development (plots only)

godag = get_godag("go-basic.obo", optional_attrs={'relationship'})
goslim = get_godag("plant_goslim.obo", optional_attrs={'relationship'})

def read_go_annotations(filename):
    annotations = {}
    with open(filename, 'r') as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.rstrip("\\n").split('\\t')
            if len(parts) != 2:
                continue
            gene, go_term = parts
            if go_term in EXCLUDE_ROOTS:
                continue
            annotations.setdefault(gene, []).append(go_term)
    return annotations

target_annotations = read_go_annotations("$target_go")
background_annotations = read_go_annotations("$background_go")

print(f"Target genes: {len(target_annotations)}")
print(f"Background genes: {len(background_annotations)}")

def map_to_goslim(go_terms, goslim_dag):
    categories = set()
    for go_id in go_terms:
        if go_id not in godag:
            continue
        go_term = godag[go_id]
        if go_id in goslim_dag and go_id not in EXCLUDE_ROOTS:
            categories.add(go_id)
        for anc in go_term.get_all_parents():
            if anc in goslim_dag and anc not in EXCLUDE_ROOTS:
                categories.add(anc)
    return list(categories)

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

target_counts = Counter()
for cats in target_categories.values():
    for cat in cats:
        target_counts[cat] += 1

background_counts = Counter()
for cats in background_categories.values():
    for cat in cats:
        background_counts[cat] += 1

# Fisher exact helpers
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

def get_term_info(go_id):
    if go_id in godag:
        term = godag[go_id]
        return {'id': go_id, 'name': term.name, 'namespace': term.namespace}
    return None

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

    tp = (a / T * 100) if T else 0.0
    bp = (c / B * 100) if B else 0.0
    fc = (tp + EPS_FC) / (bp + EPS_FC)

    results.append({
        'GO_Slim_ID': go_id,
        'GO_Slim_Term': info['name'],
        'Category': info['namespace'],
        'Target_Genes': a,
        'Background_Genes': c,
        'Target_Percent': tp,
        'Background_Percent': bp,
        'Fold_Change': fc,
        'Odds_Ratio': orat,
        'P_Value': pval,
    })

df = pd.DataFrame(results)

if not df.empty and 'P_Value' in df.columns:
    df['FDR_BH'] = bh_fdr(df['P_Value'].tolist())
else:
    df['FDR_BH'] = []

df = df.sort_values('Target_Genes', ascending=False)
df.to_csv("${sample_id}_goslim_categories.tsv", sep='\\t', index=False)

# Top categories used for protein report (same selection logic as multiplier plot)
df_plot = df[df['GO_Slim_ID'] != PLOT_EXCLUDE_GOID].copy()
df_over = df_plot[df_plot['Fold_Change'] > 1.0].copy()

top_mode = "FDR"
df_sig_over = df_over[(df_over['FDR_BH'].notna()) & (df_over['FDR_BH'] <= FDR_LIMIT)].copy()

if df_sig_over.empty:
    top_mode = "FALLBACK_LOWEST_FDR"
    df_pool = df_over[df_over['FDR_BH'].notna()].copy()
    if df_pool.empty:
        top_categories = pd.DataFrame(columns=df_plot.columns)
    else:
        df_pool = df_pool.sort_values('FDR_BH', ascending=True).head(20)
        top_categories = df_pool.sort_values('Fold_Change', ascending=False).head(TOP_PLOT_N).copy()
else:
    top_categories = df_sig_over.sort_values('Fold_Change', ascending=False).head(TOP_PLOT_N).copy()

top_category_ids = top_categories['GO_Slim_ID'].tolist()

# Protein/gene report
cat_to_genes = defaultdict(list)
for gene, cats in target_categories.items():
    for cat in cats:
        cat_to_genes[cat].append(gene)

def go_name(go_id: str) -> str:
    if go_id in godag:
        return godag[go_id].name
    return ""

rows = []
for go_id in top_category_ids:
    cat_info = df.loc[df['GO_Slim_ID'] == go_id].head(1)
    if cat_info.empty:
        continue
    cat_term = cat_info['GO_Slim_Term'].iloc[0]
    cat_ns = cat_info['Category'].iloc[0]
    cat_fc = float(cat_info['Fold_Change'].iloc[0])
    cat_fdr = cat_info['FDR_BH'].iloc[0]
    cat_p = cat_info['P_Value'].iloc[0]
    cat_t = int(cat_info['Target_Genes'].iloc[0])
    cat_b = int(cat_info['Background_Genes'].iloc[0])

    genes = sorted(set(cat_to_genes.get(go_id, [])))
    for gene in genes:
        orig_gos = target_annotations.get(gene, [])
        orig_go_names = [go_name(g) for g in orig_gos]
        mapped_slims = sorted(set(target_categories.get(gene, [])))
        mapped_slim_names = [go_name(g) for g in mapped_slims]

        rows.append({
            'Top_Category_Mode': top_mode,
            'Top_GO_Slim_ID': go_id,
            'Top_GO_Slim_Term': cat_term,
            'Top_GO_Slim_Namespace': cat_ns,
            'Top_Category_Target_Genes': cat_t,
            'Top_Category_Background_Genes': cat_b,
            'Top_Category_Fold_Change': cat_fc,
            'Top_Category_P_Value': cat_p,
            'Top_Category_FDR_BH': cat_fdr,
            'Gene_or_Protein_ID': gene,
            'Original_GO_Terms': ";".join(orig_gos),
            'Original_GO_Term_Names': ";".join([n for n in orig_go_names if n]),
            'Mapped_GO_Slim_Terms': ";".join(mapped_slims),
            'Mapped_GO_Slim_Term_Names': ";".join([n for n in mapped_slim_names if n]),
        })

prot_df = pd.DataFrame(rows)
prot_df.to_csv("${sample_id}_top_overrepresented_proteins.tsv", sep='\\t', index=False)

# NEW: same format as *_top_overrepresented_proteins.txt but for ALL categories (not only overrepresented/top)
# Build a per-category gene/protein listing for all categories present in target (ordered by Target_Genes desc)
all_cat_ids_ordered = df_plot.sort_values('Target_Genes', ascending=False)['GO_Slim_ID'].tolist()

rows_all = []
for go_id in all_cat_ids_ordered:
    cat_info = df.loc[df['GO_Slim_ID'] == go_id].head(1)
    if cat_info.empty:
        continue

    cat_term = cat_info['GO_Slim_Term'].iloc[0]
    cat_ns = cat_info['Category'].iloc[0]
    cat_fc = float(cat_info['Fold_Change'].iloc[0])
    cat_fdr = cat_info['FDR_BH'].iloc[0]
    cat_p = cat_info['P_Value'].iloc[0]
    cat_t = int(cat_info['Target_Genes'].iloc[0])
    cat_b = int(cat_info['Background_Genes'].iloc[0])

    genes = sorted(set(cat_to_genes.get(go_id, [])))
    for gene in genes:
        orig_gos = target_annotations.get(gene, [])
        orig_go_names = [go_name(g) for g in orig_gos]
        mapped_slims = sorted(set(target_categories.get(gene, [])))
        mapped_slim_names = [go_name(g) for g in mapped_slims]

        rows_all.append({
            'GO_Slim_ID': go_id,
            'GO_Slim_Term': cat_term,
            'GO_Slim_Namespace': cat_ns,
            'Category_Target_Genes': cat_t,
            'Category_Background_Genes': cat_b,
            'Category_Fold_Change': cat_fc,
            'Category_P_Value': cat_p,
            'Category_FDR_BH': cat_fdr,
            'Gene_or_Protein_ID': gene,
            'Original_GO_Terms': ";".join(orig_gos),
            'Original_GO_Term_Names': ";".join([n for n in orig_go_names if n]),
            'Mapped_GO_Slim_Terms': ";".join(mapped_slims),
            'Mapped_GO_Slim_Term_Names': ";".join([n for n in mapped_slim_names if n]),
        })

prot_all_df = pd.DataFrame(rows_all)

with open("${sample_id}_all_categories_proteins.txt", "w") as f:
    f.write(f"All GO-slim categories (same per-category report format) - ${sample_id}\\n")
    f.write("="*90 + "\\n")
    f.write(f"FDR limit reference: {FDR_LIMIT} (not used for filtering here)\\n")
    f.write(f"Plot-only excluded GO ID: {PLOT_EXCLUDE_GOID}\\n")
    f.write(f"Fold_Change definition: (Target%+{EPS_FC})/(Background%+{EPS_FC})\\n\\n")

    f.write("Per-category gene/protein IDs and GO-derived details (ALL categories):\\n")
    f.write("-"*90 + "\\n\\n")

    if prot_all_df.empty:
        f.write("No genes mapped to any GO-slim categories.\\n")
    else:
        for go_id in all_cat_ids_ordered:
            block = prot_all_df[prot_all_df['GO_Slim_ID'] == go_id].copy()
            if block.empty:
                continue

            cat_term = block['GO_Slim_Term'].iloc[0]
            fc = float(block['Category_Fold_Change'].iloc[0])
            p = block['Category_P_Value'].iloc[0]
            q = block['Category_FDR_BH'].iloc[0]
            p_str = "NA" if pd.isna(p) else f"{p:.2g}"
            q_str = "NA" if pd.isna(q) else f"{q:.2g}"

            f.write(f"{go_id} - {cat_term}\\n")
            f.write(f"  Fold-change: {fc:.2g}×   p: {p_str}   FDR: {q_str}\\n")
            f.write(f"  Genes/proteins ({len(block)}):\\n")

            for _, rr in block.sort_values('Gene_or_Protein_ID').iterrows():
                f.write(f"    - {rr['Gene_or_Protein_ID']}\\n")
                if rr['Original_GO_Terms']:
                    f.write(f"      Original GO: {rr['Original_GO_Terms']}\\n")
                if rr['Original_GO_Term_Names']:
                    f.write(f"      Original GO names: {rr['Original_GO_Term_Names']}\\n")
                if rr['Mapped_GO_Slim_Terms']:
                    f.write(f"      Mapped slim: {rr['Mapped_GO_Slim_Terms']}\\n")
                if rr['Mapped_GO_Slim_Term_Names']:
                    f.write(f"      Mapped slim names: {rr['Mapped_GO_Slim_Term_Names']}\\n")
            f.write("\\n")

with open("${sample_id}_top_overrepresented_proteins.txt", "w") as f:
    f.write(f"Top overrepresented GO-slim categories (used for multiplier plot) - ${sample_id}\\n")
    f.write("="*90 + "\\n")
    f.write(f"Selection mode: {top_mode}\\n")
    f.write(f"FDR limit: {FDR_LIMIT}\\n")
    f.write(f"Plot-only excluded GO ID: {PLOT_EXCLUDE_GOID}\\n")
    f.write(f"Overrepresented definition: Fold_Change = (Target%+{EPS_FC})/(Background%+{EPS_FC}) > 1\\n")
    f.write(f"Max categories listed: {TOP_PLOT_N}\\n\\n")

    if top_categories.empty:
        f.write("No overrepresented categories available for reporting.\\n")
    else:
        f.write("Top categories:\\n")
        for _, r in top_categories.sort_values('Fold_Change', ascending=False).iterrows():
            p = r['P_Value']
            q = r['FDR_BH']
            p_str = "NA" if pd.isna(p) else f"{p:.2g}"
            q_str = "NA" if pd.isna(q) else f"{q:.2g}"
            f.write(f"- {r['GO_Slim_ID']} | {r['GO_Slim_Term']} | FC:{r['Fold_Change']:.2g}× | p:{p_str} | FDR:{q_str}\\n")
        f.write("\\n")

    f.write("Per-category gene/protein IDs and GO-derived details:\\n")
    f.write("-"*90 + "\\n\\n")

    if prot_df.empty:
        f.write("No genes mapped to the selected top categories.\\n")
    else:
        for go_id in top_category_ids:
            block = prot_df[prot_df['Top_GO_Slim_ID'] == go_id].copy()
            if block.empty:
                continue
            cat_term = block['Top_GO_Slim_Term'].iloc[0]
            fc = float(block['Top_Category_Fold_Change'].iloc[0])
            p = block['Top_Category_P_Value'].iloc[0]
            q = block['Top_Category_FDR_BH'].iloc[0]
            p_str = "NA" if pd.isna(p) else f"{p:.2g}"
            q_str = "NA" if pd.isna(q) else f"{q:.2g}"

            f.write(f"{go_id} - {cat_term}\\n")
            f.write(f"  Fold-change: {fc:.2g}×   p: {p_str}   FDR: {q_str}\\n")
            f.write(f"  Genes/proteins ({len(block)}):\\n")

            for _, rr in block.sort_values('Gene_or_Protein_ID').iterrows():
                f.write(f"    - {rr['Gene_or_Protein_ID']}\\n")
                if rr['Original_GO_Terms']:
                    f.write(f"      Original GO: {rr['Original_GO_Terms']}\\n")
                if rr['Original_GO_Term_Names']:
                    f.write(f"      Original GO names: {rr['Original_GO_Term_Names']}\\n")
                if rr['Mapped_GO_Slim_Terms']:
                    f.write(f"      Mapped slim: {rr['Mapped_GO_Slim_Terms']}\\n")
                if rr['Mapped_GO_Slim_Term_Names']:
                    f.write(f"      Mapped slim names: {rr['Mapped_GO_Slim_Term_Names']}\\n")
            f.write("\\n")

# Summary
with open("${sample_id}_goslim_summary.txt", 'w') as f:
    f.write(f"Plant GO Slim Analysis - ${sample_id}\\n")
    f.write("="*60 + "\\n\\n")
    f.write(f"Target genes with GO slim annotations: {len(target_categories)}/{len(target_annotations)} ({(len(target_categories)/len(target_annotations)*100) if len(target_annotations) else 0:.1f}%)\\n")
    f.write(f"Background genes with GO slim annotations: {len(background_categories)}/{len(background_annotations)} ({(len(background_categories)/len(background_annotations)*100) if len(background_annotations) else 0:.1f}%)\\n")
    f.write(f"Unique GO slim categories found: {len(df)} (root terms excluded)\\n\\n")

    sig = df[(df['FDR_BH'].notna()) & (df['FDR_BH'] <= FDR_LIMIT)].copy()
    f.write(f"Significant categories (FDR_BH<={FDR_LIMIT}): {len(sig)}/{len(df)}\\n\\n")

    f.write(f"All categories with FDR_BH <= {FDR_LIMIT} (sorted by FDR):\\n")
    f.write("-"*100 + "\\n")
    if sig.empty:
        f.write("  None\\n\\n")
    else:
        sig = sig.sort_values('FDR_BH', ascending=True)
        for _, row in sig.iterrows():
            p = row.get('P_Value', None)
            q = row.get('FDR_BH', None)
            p_str = "NA" if pd.isna(p) else f"{p:.2g}"
            q_str = "NA" if pd.isna(q) else f"{q:.2g}"
            f.write(
                f"{row['GO_Slim_ID']} | {row['GO_Slim_Term']} | "
                f"T:{int(row['Target_Genes'])} ({row['Target_Percent']:.1f}%) | "
                f"BG:{int(row['Background_Genes'])} ({row['Background_Percent']:.1f}%) | "
                f"FC:{row['Fold_Change']:.2g}× | p:{p_str} | FDR:{q_str}\\n"
            )
        f.write("\\n")

    f.write("Top 10 GO Slim Categories in Target (with Fisher p-values and BH FDR):\\n")
    f.write("-"*70 + "\\n")
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
        f.write(f"  OddsRatio: {or_str}  p: {p_str}  FDR_BH: {q_str}\\n\\n")

print("GO Slim analysis complete (including p-values, FDR, and top-category protein report)")
EOF

  # 4. Multiplier plot (RAW FC) - ONLY overrepresented - excludes GO:0048856 from plots
  python3 << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt

try:
    df_all = pd.read_csv("${sample_id}_goslim_categories.tsv", sep='\\t')
    if df_all.empty:
        raise ValueError("No categories found in TSV (empty dataframe).")

    df_all = df_all[df_all['GO_Slim_ID'] != "GO:0048856"].copy()

    df_over = df_all[df_all['Fold_Change'] > 1.0].copy()
    if df_over.empty:
        raise ValueError("No overrepresented categories (Fold_Change > 1) to plot.")

    df_sig = df_over[(df_over['FDR_BH'].notna()) & (df_over['FDR_BH'] <= 0.1)].copy()

    fallback = False
    df = df_sig
    if df.empty:
        fallback = True
        df = df_over[df_over['FDR_BH'].notna()].copy()
        if df.empty:
            raise ValueError("FDR_BH is missing for all overrepresented categories.")
        df = df.sort_values('FDR_BH', ascending=True).head(20)

    plot_df = df.sort_values('Fold_Change', ascending=False).head(15).copy()
    plot_df = plot_df.sort_values('Fold_Change', ascending=True)

    fig, ax = plt.subplots(figsize=(12, 9))
    y_pos = range(len(plot_df))
    ax.barh(y_pos, plot_df['Fold_Change'])

    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(plot_df['GO_Slim_Term'].astype(str).str.wrap(45))

    ax.axvline(1.0, linewidth=1)
    ax.set_xlabel('BG multiplier (Target% / Background%)')

    if fallback:
        ax.set_title('No overrepresented terms pass FDR≤0.1; showing lowest-FDR overrepresented - ${sample_id}')
    else:
        ax.set_title('Overrepresented GO Slim Categories (FDR≤0.1): BG Multiplier - ${sample_id}')

    ax.set_xticks([1, 1.5, 2, 3, 5])
    ax.set_xticklabels(['1×', '1.5×', '2×', '3×', '5×'])

    for i, row in enumerate(plot_df.itertuples(index=False)):
        p = row.P_Value
        q = row.FDR_BH
        p_str = "NA" if pd.isna(p) else f"{p:.2g}"
        q_str = "NA" if pd.isna(q) else f"{q:.2g}"
        label = f"T:{row.Target_Percent:.1f}%  B:{row.Background_Percent:.1f}%  ({row.Fold_Change:.2g}× bg)  p:{p_str}  FDR:{q_str}"
        ax.text(row.Fold_Change, i, "  " + label, va='center')

    plt.tight_layout()
    plt.savefig("${sample_id}_goslim_plot_multiplier.png", dpi=300, bbox_inches='tight')
    plt.close()

    print("Overrepresented multiplier plot created successfully")
except Exception as e:
    print(f"Could not create multiplier plot: {e}")
EOF

  # 5. Dot plot (volcano): annotate each dot with category name
  python3 << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

try:
    df = pd.read_csv("${sample_id}_goslim_categories.tsv", sep='\\t')
    if df.empty:
        raise ValueError("No categories found in TSV (empty dataframe).")

    df = df.copy()
    df['Log2_FC'] = np.log2(df['Fold_Change'].astype(float).clip(lower=1e-300))
    p = df['P_Value'].astype(float).fillna(1.0).clip(lower=1e-300)
    df['NegLog10P'] = -np.log10(p)

    fig, ax = plt.subplots(figsize=(10, 9))
    ax.scatter(df['Log2_FC'], df['NegLog10P'], alpha=0.7)
    ax.axvline(0, linewidth=1)

    # Annotate: for readability, annotate only the most "extreme" points
    # (otherwise this plot becomes unreadable with many categories).
    # We'll label: top 30 by -log10(p) + |log2FC|
    df['LabelScore'] = df['NegLog10P'] + df['Log2_FC'].abs()
    label_df = df.sort_values('LabelScore', ascending=False).head(30)

    for _, r in label_df.iterrows():
        ax.text(float(r['Log2_FC']), float(r['NegLog10P']), " " + str(r['GO_Slim_Term']), fontsize=7)

    ax.set_xlabel('log2(Target% / Background%)')
    ax.set_ylabel('-log10(p-value)')
    ax.set_title('GO Slim Volcano (Dot Plot): Enrichment vs Significance - ${sample_id}')

    plt.tight_layout()
    plt.savefig("${sample_id}_goslim_plot_volcano.png", dpi=300, bbox_inches='tight')
    plt.close()

    print("Volcano (dot) plot with labels created successfully")
except Exception as e:
    print(f"Could not create volcano plot: {e}")
EOF

  # 6. Top-by-Target% plot (ONLY overrepresented), excludes GO:0048856 from plots
  python3 << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt

try:
    df_all = pd.read_csv("${sample_id}_goslim_categories.tsv", sep='\\t')
    if df_all.empty:
        raise ValueError("No categories found in TSV (empty dataframe).")

    df_all = df_all[df_all['GO_Slim_ID'] != "GO:0048856"].copy()

    df_over = df_all[df_all['Fold_Change'] > 1.0].copy()
    if df_over.empty:
        raise ValueError("No overrepresented categories (Fold_Change > 1) to plot.")

    df_sig = df_over[(df_over['FDR_BH'].notna()) & (df_over['FDR_BH'] <= 0.1)].copy()

    fallback = False
    df = df_sig
    if df.empty:
        fallback = True
        df = df_over.copy()

    top_df = df.sort_values('Target_Percent', ascending=False).head(15).copy()
    top_df = top_df.iloc[::-1]

    fig, ax = plt.subplots(figsize=(12, 8))
    y_pos = range(len(top_df))

    ax.barh(y_pos, top_df['Target_Percent'], label='Target')
    ax.barh(y_pos, top_df['Background_Percent'], label='Background', alpha=0.7)

    ax.set_yticks(list(y_pos))
    ax.set_yticklabels(top_df['GO_Slim_Term'].astype(str).str.wrap(45))
    ax.set_xlabel('Percentage of Genes (%)')

    if fallback:
        ax.set_title('No overrepresented terms pass FDR≤0.1; showing top overrepresented by Target% - ${sample_id}')
    else:
        ax.set_title('Top Overrepresented GO Slim Categories by Target% (FDR≤0.1) - ${sample_id}')

    ax.legend()

    plt.tight_layout()
    plt.savefig("${sample_id}_goslim_plot_top_target.png", dpi=300, bbox_inches='tight')
    plt.close()

    print("Top Target% (overrepresented) plot created successfully")
except Exception as e:
    print(f"Could not create top-target plot: {e}")
EOF

  # Keep compatibility: default plot name expected downstream
  cp -f "${sample_id}_goslim_plot_multiplier.png" "${sample_id}_goslim_plot.png"
  """
}
