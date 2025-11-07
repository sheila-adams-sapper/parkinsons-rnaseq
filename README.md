# Parkinson's Disease RNA-seq Interactive Dashboard

Interactive dashboard for exploring RNA-seq results from Parkinson's disease research.

## 🚀 Live Dashboard
- **Production**: https://parkinsons-rnaseq.vercel.app
- **Repository**: https://github.com/sheila-adams-sapper/parkinsons-rnaseq
- **Deployment**: Vercel (auto-deploy from main branch)
- **Interactive Features**: Volcano plots, PCA plots, pathway analysis, gene expression tables

## 🛠️ Deployment
This dashboard is automatically deployed to Vercel from the `docs/` directory. 
Any push to the main branch triggers a new deployment.

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
