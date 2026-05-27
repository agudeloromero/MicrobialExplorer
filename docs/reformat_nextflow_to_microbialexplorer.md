# Reformatting Nextflow Output for MicrobialExplorer

This guide explains how to prepare files from a Nextflow 16S pipeline
for upload into MicrobialExplorer. A standalone R script is provided
to automate the reformatting steps.

---

## Files you need from your Nextflow output

| File | Location in Nextflow output | Purpose |
|------|-----------------------------|---------|
| OTU/ASV table | `Analysis/results/best_tax_merged_freq_tax.tsv` | Feature counts per sample |
| Taxonomy | `Analysis/results/tax_export/taxonomy.tsv` | Taxonomic classification |
| Phylogenetic tree | `Analysis/results/phylogeny_diversity/phylotree_mafft_rooted.nwk` | Required for Faith's PD and UniFrac |
| Metadata | Created by the user | Sample grouping variables |

---

## Metadata file

MicrobialExplorer requires a metadata file you create yourself.
Download the template below and fill it in with your sample information.

**Required format:**

- Tab-separated (`.tsv`)
- First column must be named `sample_id`
- Sample IDs must match the column names in your OTU table exactly
- Do not include the `Control_reagents` sample

**Example:**

```
sample_id	group	timepoint	sex
A17	Group1	T1	F
A18	Group1	T1	M
B16	Group2	T1	F
```

A metadata template is available to download from the MicrobialExplorer
upload page.

---

## Reformatting with the R script

The script `reformat_nextflow.R` automates all file preparation steps.
It requires R with the `dplyr` and `tidyr` packages installed.

### Install required packages (once only)

```r
install.packages(c("dplyr", "tidyr"))
```

### Run the script

```bash
Rscript reformat_nextflow.R \
  --otu   path/to/best_tax_merged_freq_tax.tsv \
  --tax   path/to/taxonomy.tsv \
  --outdir path/to/output_folder
```

### Output files

The script will create the following files in your output folder:

| File | Upload as |
|------|-----------|
| `otu_table_clean.tsv` | OTU / ASV table |
| `taxonomy_ranks.tsv` | Taxonomy table |

Upload these two files along with your metadata and the
`phylotree_mafft_rooted.nwk` tree file directly — no reformatting
needed for the tree.

---

## What the script does

1. **OTU table** — removes the `Sequence`, `Taxon`, and `Confidence`
   columns from the merged file, removes the `#q2:types` row, cleans
   sample names by removing the `.fastq` suffix, and removes the
   `Control_reagents` column.

2. **Taxonomy** — splits the single `Taxon` string into separate rank
   columns (`Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`,
   `Species`) and removes rank prefixes (`d__`, `p__`, etc.).

---

## Upload order in MicrobialExplorer

1. OTU / ASV table → `otu_table_clean.tsv`
2. Taxonomy table → `taxonomy_ranks.tsv`
3. Sample metadata → your completed metadata file
4. Phylogenetic tree → `phylotree_mafft_rooted.nwk`

Click **Import data** once all four files are uploaded.
