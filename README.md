# Oncofetal microenvironment modeling in assembloids

The vascular, immune, and stromal compartments of tumors can reactivate fetal gene programs, creating an oncofetal microenvironment. This oncofetal microenvironment regulates cancer growth and therapy resistance, yet human cancer models capture it poorly. Here, we developed tumor immune vascular assembloids (TIVAs) by differentiating human pluripotent stem cells into vascular-immune organoids and assembling them with patient-derived tumor organoids. We generated TIVAs across lung, colorectal, and hepatocellular carcinoma models, endowing tumor organoids with vascular, innate immune, and stromal compartments. Single-cell and spatial transcriptomics revealed an IL1B+ immunosuppressive macrophage niche and an oncofetal microenvironment characterized by a PLVAP+ vascular niche dependent on VEGF-KDR signaling. The PLVAP+ vascular niche conferred resistance to tyrosine kinase inhibitor therapy, as seen in corresponding lung adenocarcinoma patients. Targeting VEGF-KDR signaling restored sensitivity to therapy and disrupted the oncofetal microenvironment. TIVAs thus provide a platform for modeling the oncofetal microenvironment and reveal it as a targetable vulnerability.

## Download links for processed files

- [VIO](https://www.dropbox.com/scl/fi/z38qp5aglwuiwxqviob6r/F1_VIO.rds?rlkey=6jz62fpx38uikgtr589amfdja&dl=1)
- [PDO](https://www.dropbox.com/scl/fi/t0owxr2jovh01ugofor0v/lam.rds?rlkey=5ussjc56q39mwjblgk3hmq58r&dl=1)
- [CIVA (scRNAseq)](https://www.dropbox.com/scl/fi/l95375duse854fqfqypcm/Flex.rds?rlkey=68xx629wemlah9na412jknowl&dl=1)
- [CIVA (VisiumHD)](https://www.dropbox.com/scl/fi/pmp4k67r9owuvag45c8fk/F4_VisiumHD.h5ad?rlkey=74crcopktkabr1phpaiz1moje&dl=1)
- [Public lung cancer samples](https://www.dropbox.com/scl/fi/efxs6vumvv4ffdyuih2o2/zemin_lung.rds?rlkey=66bnxdo3rb1jsve9mlrpkgqyc&dl=1)

The files for the Fig 3 cell–cell interaction analysis can be downloaded from [cellphonefiles](https://www.dropbox.com/scl/fo/wped4hee8ew7g7lydxel3/AKwpdcIvg4DFrbtTugDOi5o?rlkey=69hav7i9gq8af7r8xtnot91n2&dl=1):

```text
cellphonefiles/
├── metadata.tsv
├── normalised_log_counts.h5ad
└── cellphonedb.zip
```

## Analysis scripts in this repository

The repository contains analysis workflows for vascular-immune organoids, tumor immune vascular assembloids, public lung cancer datasets, cell–cell interactions, gene set enrichment, and spatial transcriptomics. The scripts below are annotated with their corresponding manuscript panels.

```text
TIVA/
│
├── Fig1_VIO.ipynb                 # Figure 1E–M, together with Fig1_EC.ipynb
├── Fig1_EC.ipynb                  # Figure 1E–M, together with Fig1_VIO.ipynb
├── Fig1-4_Correlation.ipynb        # Figure 1G
│
├── Fig3_CIVA_.ipynb               # Figure 3A–B
├── Fig3_Maynard_.ipynb            # Figure 3C–D
├── Fig3_Interaction.ipynb         # Figure 3I
├── Fig3_GSEA.ipynb                # Figure 3; gene set enrichment analysis
│
├── Fig4_2umHD (1).ipynb           # Figure 4B and 4D
├── Fig4_COMMOT_v2_clean.ipynb     # Figure 4F
├── Fig4_LUAD.py                  # Figure 4E
├── Fig4_VisiumToR.ipynb
├── Fig4_violin.ipynb              # Figure 4G
│
├── Fig6_LUAD.py                  # Figure 6A
└── Fig6_PLVAP.R                  # Figure 6A
```

Updated: 18 September 2026
