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
  
  # 1. Download plant-specific GO slim
  wget -q https://current.geneontology.org/ontology/subsets/goslim_plant.obo -O plant_goslim.obo
  
  # 2. Download full GO ontology
  wget -q http://purl.obolibrary.org/obo/go/go-basic.obo -O go-basic.obo
  
  # 3. Run GOATOOLS goslim analysis
  python3 << 'EOF'
import pandas as pd
from goatools.gosubdag.gosubdag import GoSubDag
from goatools.base import get_godag
from goatools.gosubdag.rpt.write_hierarchy import WrHierGO
  
# Load GO slim and full ontology
godag = get_godag("go-basic.obo", optional_attrs={'relationship'})
goslim = get_godag("plant_goslim.obo", optional_attrs={'relationship'})
  
# Define root terms to exclude
EXCLUDE_ROOTS = {'GO:0009987', 'GO:0005622', 'GO:0008152', 'GO:0005737'}
  
# Read GO annotations
def read_go_annotations(filename):
    annotations = {}
    with open(filename, 'r') as f:
        for line in f:
            if line.strip():
                parts = line.strip().split('\t')
                if len(parts) == 2:
                    gene, go_term = parts
                    # Skip root terms immediately
                    if go_term in EXCLUDE_ROOTS:
                        continue
                    if gene not in annotations:
                        annotations[gene] = []
                    annotations[gene].append(go_term)
    return annotations
  
target_annotations = read_go_annotations("$target_go")
background_annotations = read_go_annotations("$background_go")
  
print(f"Target genes: {len(target_annotations)}")
print(f"Background genes: {len(background_annotations)}")
  
# Map specific GO terms to plant GO slim categories
def map_to_goslim(go_terms, goslim_dag):
    categories = set()
    for go_id in go_terms:
        if go_id in godag:
            go_term = godag[go_id]
            # Check if this term is in plant GO slim (and not a root)
            if go_id in goslim_dag and go_id not in EXCLUDE_ROOTS:
                categories.add(go_id)
            # Also check ancestors (excluding roots)
            for ancestor in go_term.get_all_parents():
                if ancestor in goslim_dag and ancestor not in EXCLUDE_ROOTS:
                    categories.add(ancestor)
    return list(categories)
  
# Classify each gene
target_categories = {}
for gene, go_terms in target_annotations.items():
    categories = map_to_goslim(go_terms, goslim)
    if categories:
        target_categories[gene] = categories
  
background_categories = {}
for gene, go_terms in background_annotations.items():
    categories = map_to_goslim(go_terms, goslim)
    if categories:
        background_categories[gene] = categories
  
# Count categories
from collections import Counter
  
target_counts = Counter()
for categories in target_categories.values():
    for cat in categories:
        target_counts[cat] += 1
  
background_counts = Counter()
for categories in background_categories.values():
    for cat in categories:
        background_counts[cat] += 1
  
# Get GO term names
def get_term_info(go_id):
    if go_id in godag:
        term = godag[go_id]
        return {
            'id': go_id,
            'name': term.name,
            'namespace': term.namespace
        }
    return None
  
# Create results table (excluding root terms)
results = []
for go_id in set(list(target_counts.keys()) + list(background_counts.keys())):
    # Skip root terms in output
    if go_id in EXCLUDE_ROOTS:
        continue
    term_info = get_term_info(go_id)
    if term_info:
        results.append({
            'GO_Slim_ID': go_id,
            'GO_Slim_Term': term_info['name'],
            'Category': term_info['namespace'],
            'Target_Genes': target_counts.get(go_id, 0),
            'Background_Genes': background_counts.get(go_id, 0),
            'Target_Percent': (target_counts.get(go_id, 0) / len(target_annotations) * 100) if target_annotations else 0,
            'Background_Percent': (background_counts.get(go_id, 0) / len(background_annotations) * 100) if background_annotations else 0
        })
  
# Sort by target genes (descending)
df = pd.DataFrame(results)
df = df.sort_values('Target_Genes', ascending=False)
  
# Save results
df.to_csv("${sample_id}_goslim_categories.tsv", sep='\t', index=False)
  
# Save summary
with open("${sample_id}_goslim_summary.txt", 'w') as f:
    f.write(f"Plant GO Slim Analysis - ${sample_id}\\n")
    f.write("="*60 + "\\n\\n")
    f.write(f"Target genes with GO slim annotations: {len(target_categories)}/{len(target_annotations)} ({len(target_categories)/len(target_annotations)*100:.1f}%)\\n")
    f.write(f"Background genes with GO slim annotations: {len(background_categories)}/{len(background_annotations)} ({len(background_categories)/len(background_annotations)*100:.1f}%)\\n")
    f.write(f"Unique plant GO slim categories found: {len(df)} (root terms excluded)\\n\\n")
    
    # Top categories
    f.write("Top 10 GO Slim Categories in Target:\\n")
    f.write("-"*40 + "\\n")
    for _, row in df.head(10).iterrows():
        f.write(f"{row['GO_Slim_ID']}: {row['GO_Slim_Term']}\\n")
        f.write(f"  Target: {row['Target_Genes']} genes ({row['Target_Percent']:.1f}%)\\n")
        f.write(f"  Background: {row['Background_Genes']} genes ({row['Background_Percent']:.1f}%)\\n\\n")
  
print("GO Slim analysis complete")
EOF
  
  # 4. Create simple plot
  python3 << 'EOF'
import pandas as pd
import matplotlib.pyplot as plt
  
try:
    df = pd.read_csv("${sample_id}_goslim_categories.tsv", sep='\t')
    
    # Top 15 categories
    top_df = df.head(15).copy()
    
    # Create bar plot
    fig, ax = plt.subplots(figsize=(12, 8))
    
    y_pos = range(len(top_df))
    ax.barh(y_pos, top_df['Target_Percent'], color='steelblue', label='Target')
    ax.barh(y_pos, top_df['Background_Percent'], color='lightgray', label='Background', alpha=0.7)
    
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_df['GO_Slim_Term'].str.wrap(40))
    ax.set_xlabel('Percentage of Genes (%)')
    ax.set_title('Plant GO Slim Categories - ${sample_id}')
    ax.legend()
    
    plt.tight_layout()
    plt.savefig("${sample_id}_goslim_plot.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    print("Plot created successfully")
except Exception as e:
    print(f"Could not create plot: {e}")
EOF
  """
}