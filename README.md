# 🧠 MIMIC Analysis Pipeline

This repository contains an **end-to-end data analysis workflow** for extracting, preprocessing, and analyzing MIMIC-IV data using R and Python.  
It is designed for reproducible comorbidity and embedding-based analyses.

---

## 📚 Overview

**Main Steps:**
1. **Data Extraction** (`R/01_extract_data.R`): SQL-based queries from MIMIC-IV.
2. **Preprocessing** (`R/02_preprocess.R`): Cleaning, feature derivation, comorbidity flags.
3. **Embedding Analysis** (`R/03_embedding_analysis.R`): Word2Vec / UMAP / clustering.
4. **Visualization** (`R/04_visualization.R`): Network and embedding plots.

---

## ⚙️ Environment Setup

### R Dependencies
```r
install.packages(c(
  "DBI", "RPostgreSQL", "tidyverse", "lubridate",
  "data.table", "text2vec", "uwot", "ggplot2", "igraph"
))
