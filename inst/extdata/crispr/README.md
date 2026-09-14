# CRISPR reference gene sets

Reference gene sets for pooled CRISPR-screen QC (used by ExploreMageckCounts).

- `CEGv2.txt` — core-essential genes (positive controls for essentiality).
  Hart T. et al., "Evaluation and Design of Genome-Wide CRISPR/SpCas9 Knockout
  Screens", G3 2017; 7(8):2719-2727. https://doi.org/10.1534/g3.117.041277
- `NEGv1.txt` — non-essential genes (negative controls).
  Hart T. et al., "High-Resolution CRISPR Screens Reveal Fitness Genes and
  Genotype-Specific Cancer Liabilities", Cell 2015; 163(6):1515-1526.
  https://doi.org/10.1016/j.cell.2015.11.015

Source: https://github.com/hart-lab/bagel (CEGv2.txt, NEGv1.txt).
Columns: GENE, HGNC_ID, ENTREZ_ID (tab-separated). Human symbols; for mouse,
match is attempted by upper-casing gene symbols (ortholog mapping is future work).
