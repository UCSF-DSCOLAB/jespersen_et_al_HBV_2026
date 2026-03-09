This repository contains original code used for the analysis of CITE-seq data generated for Jespersen & Avanesyan et al., Clinical cure of chronic hepatitis B is dependent on priming and perpetuation of robust CD4⁺ T cell responses (2026). (bioRxiv: https://www.biorxiv.org/content/10.1101/2025.09.29.677401v1.full)

```
├── data: auxiliary data files
│   ├── citeseq_processing: contains pseudobulked gene expression matrix
│   ├── dge_result: output of differential gene expression analysis
│   ├── diff_freq_analysis: output of differential abundance analysis
│   ├── figS16_related: input and output files for Figs16 code
│   └── sample_metadata.csv: sample metadata
├── LICENSE
├── README.md
└── src: original analysis code
    ├── citeseq_processing: code for CITE-seq data processing, including DSB-based normalization of ADT data, batch-integration of RNA and ADT data, multimodal data integration, clustering, and pseudobulk expression calculation
    ├── dge_analysis: code for differential gene expression analysis
    ├── diff_freq_analysis: code for differential abundance analysis
    └── figS16_related: code for comparing gene signature scores between outcome groups at different timepoints
```
