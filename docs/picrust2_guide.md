# Running PICRUSt2 for MicrobialExplorer Functional Analysis

PICRUSt2 predicts the functional potential of your microbiome from 16S
amplicon data. This guide walks you through installation and running the
analysis to produce files compatible with MicrobialExplorer's Functional
Prediction module.

---

## What you need before starting

| File | Where to get it |
|------|-----------------|
| ASV sequences (FASTA) | `Analysis/dada2/dada2_ASV.fasta` from your Nextflow output |
| OTU/ASV count table | `otu_table_clean.tsv` — prepared with `reformat_nextflow.R` |

> **Important:** The ASV IDs in the FASTA file must match the row names
> in your OTU table. If you used `reformat_nextflow.R`, this is already
> the case.

---

## Step 1 — Install PICRUSt2

PICRUSt2 is installed using conda. If you do not have conda, download
Miniconda from https://docs.conda.io/en/latest/miniconda.html first.

```bash
# Create a dedicated environment
conda create -n picrust2 -c bioconda -c conda-forge picrust2=2.5.2

# Activate the environment
conda activate picrust2

# Verify installation
picrust2_pipeline.py --version
```

You will also need the BIOM tools to convert your OTU table:

```bash
conda install -n picrust2 -c bioconda biom-format
```

---

## Step 2 — Convert OTU table to BIOM format

PICRUSt2 requires the OTU table in BIOM format, not TSV.

```bash
biom convert \
  -i otu_table_clean.tsv \
  -o otu_table.biom \
  --table-type "OTU table" \
  --to-hdf5
```

---

## Step 3 — Run PICRUSt2

```bash
picrust2_pipeline.py \
  --study_fasta  dada2_ASV.fasta \
  --input        otu_table.biom \
  --output       picrust2_output \
  --processes    4
```

- `--processes` sets the number of CPU cores. Increase for faster
  analysis on a server (e.g. `--processes 16`).
- Runtime is typically 20–60 minutes depending on your dataset size
  and machine.

---

## Step 4 — Check the output

When PICRUSt2 finishes, your output folder will contain:

```
picrust2_output/
├── KO_metagenome_out/
│   └── pred_metagenome_unstrat.tsv.gz     ← KEGG Orthology predictions
├── pathways_out/
│   └── path_abun_unstrat.tsv.gz           ← Pathway abundance predictions
└── marker_predicted_and_nsti.tsv.gz       ← NSTI quality scores
```

---

## Step 5 — Upload to MicrobialExplorer

In the Functional Prediction module, select **Upload PICRUSt2 files**
and provide the following:

| MicrobialExplorer field | PICRUSt2 file |
|-------------------------|---------------|
| Pathway abundance file | `picrust2_output/pathways_out/path_abun_unstrat.tsv.gz` |
| KO abundance file | `picrust2_output/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz` |
| NSTI file (optional) | `picrust2_output/marker_predicted_and_nsti.tsv.gz` |

---

## Understanding the outputs

**NSTI (Nearest Sequenced Taxon Index)** measures how closely your ASVs
match reference genomes. Lower values mean more reliable predictions.
Values below 0.15 are considered acceptable; values above 0.15 indicate
the prediction for that sample may be unreliable.

**Pathway abundance** (`path_abun_unstrat.tsv.gz`) contains predicted
MetaCyc pathway abundances per sample. This drives the Pathway Abundance,
Heatmap, and Differential Abundance tabs in MicrobialExplorer.

**KO abundance** (`pred_metagenome_unstrat.tsv.gz`) contains predicted
KEGG Orthology gene family abundances per sample. This drives the KEGG
Categories and Functional Diversity tabs.

---

## Troubleshooting

**Error: sequences not in reference database**
Some ASVs may not place in the reference phylogeny. PICRUSt2 will
exclude these automatically and report how many were excluded. A small
number of excluded ASVs is normal.

**Very high NSTI values**
High NSTI values suggest your samples contain taxa that are poorly
represented in reference databases. This is common for non-human or
environmental microbiomes. Results should be interpreted with caution
for samples with mean NSTI above 0.15.

**Memory errors**
Reduce the number of processes with `--processes 2` or run on a machine
with more RAM (16 GB minimum recommended).
