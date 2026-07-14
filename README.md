This repository contains original code used for the analysis of CITE-seq data generated for Jespersen & Avanesyan et al., Clinical cure of chronic hepatitis B is associated with priming and perpetuation of hepatic CD4+ T cell responses (Science Translational Medicine, 2026) https://www.science.org/doi/full/10.1126/scitranslmed.adx1523.

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
