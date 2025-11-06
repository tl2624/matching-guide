# matchingguide <a href="https://tl2624.github.io/matching-guide"><img src="https://img.shields.io/badge/docs-pkgdown-blue" align="right" /></a>

Companion R package & website for  
**_Building a Design-Based Matching Pipeline: From Principles to Practical Implementation in R_**  
by Thomas Leavitt  
Assistant Professor, Marxe School of Public and International Affairs  
Baruch College, City University of New York (CUNY)

---

## Overview

This repository provides:

- 📘 **The full guide** rendered as an HTML vignette (and PDF in `vignettes/`).
- 🧰 **Reusable helper functions** under `R/`:
  - `fine_strat_var_est()` – conservative variance estimator for finely stratified designs  
  - `worst_case_IPW()` – sensitivity-weighted IPW estimator for matched sets
- 🗂️ **Example dataset** `peace_pre_match.rds` in `inst/extdata/`:
  ```r
  install.packages("remotes")
  remotes::install_github("tl2624/matching-guide")
  library(matchingguide)

  # Load example data
  d <- peace_pre_match()

  # Access raw path
  peace_pre_match_path()