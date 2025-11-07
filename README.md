# Parkinson's Disease RNA-seq Interactive Dashboard

Interactive dashboard for exploring RNA-seq results from Parkinson's disease research.

## 🔗 Live Dashboard
- **Static Dashboard**: https://sheila-adams-sapper.github.io/parkinsons-rnaseq/
- **Interactive Features**: Volcano plots, pathway analysis, gene expression tables

## 📊 Key Results
- Directional pathway analysis showing up/down-regulated pathways
- Differential expression between Parkinson's vs controls  
- Exercise training effects (pre vs post-training)

## 📖 Study Information
- **Publication**: doi: 10.3389/fphys.2020.00653
- **Title**: Rehabilitative Impact of Exercise Training on Human Skeletal Muscle Transcriptional Programs in Parkinson's Disease
- **Data**: 12 Parkinson's patients, 12 controls, 5 patients pre/post exercise training

## ��️ Technologies
- **Analysis**: R/Bioconductor, Snakemake, DESeq2, clusterProfiler
- **Dashboard**: Shiny, plotly, DT
- **Deployment**: GitHub Pages (static HTML)

## 📁 Repository Contents
- `dashboard/` - Interactive Shiny application
- `docs/` - Static HTML dashboard for GitHub Pages
- `results/` - Analysis results (CSV files)
- `config/` - Sample metadata
