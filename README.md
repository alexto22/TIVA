# Oncofetal microenvironment modeling in assembloids

The vascular, immune, and stromal compartments of tumors can reactivate fetal gene programs, creating an oncofetal microenvironment. This oncofetal microenvironment regulates cancer growth and therapy resistance, yet human cancer models capture it poorly. Here, we developed tumor immune vascular assembloids (TIVAs) by differentiating human pluripotent stem cells into vascular-immune organoids and assembling them with patient-derived tumor organoids. We generated TIVAs across lung, colorectal, and hepatocellular carcinoma models, endowing tumor organoids with vascular, innate immune, and stromal compartments. Single-cell and spatial transcriptomics revealed an IL1B+ immunosuppressive macrophage niche and an oncofetal microenvironment characterized by a PLVAP+ vascular niche dependent on VEGF-KDR signaling. The PLVAP+ vascular niche conferred resistance to tyrosine kinase inhibitor therapy, as seen in corresponding lung adenocarcinoma patients. Targeting VEGF-KDR signaling restored sensitivity to therapy and disrupted the oncofetal microenvironment. TIVAs thus provide a platform for modeling the oncofetal microenvironment and reveal it as a targetable vulnerability.

## Download links for processed files

- [VIO](https://1drv.ms/u/c/afe99e57b70de78c/IQDyyciB7x3QSb4qkRV3j8uPAd_yj6K6t6HxFEKWXQbO4JA?e=BYxcHs)
- [PDO](https://1drv.ms/u/c/afe99e57b70de78c/IQAebqMHS6anRKpiA7DcaOKPAd6HNqHVSIOyN6_pq7Fd13g?e=Vpv1V5)
- [CIVA (scRNAseq)](https://1drv.ms/u/c/afe99e57b70de78c/IQDTIhTg3gtzTqeOS04442opAZ5HXZnHa_4uoY1BGZhD5SM?e=yX50jz)
- [CIVA (VisiumHD)](https://1drv.ms/u/c/afe99e57b70de78c/IQBvPV-VjDRySYdzQnUoKzycAca5My1AQoRCklRVXKJw2Ms?e=8hLNMv)
- [Public lung cancer samples](https://1drv.ms/u/c/afe99e57b70de78c/IQD4zizjt9I3SqkSF5YAOV6HAbIYfqgPzpIgOWS0vsh0IXE?e=xoV4c0)

The files for the cell–cell interaction analysis can be downloaded from [cellphonefiles](https://1drv.ms/f/c/afe99e57b70de78c/IgAI2aysTshgS6ib-ZKALlvBAVoSyAAHw9Fkd-ckQGQ1FTU?e=Mouj0t):

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
