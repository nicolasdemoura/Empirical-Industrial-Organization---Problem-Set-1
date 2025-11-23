# Empirical Industrial Organization - Problem Set 1

**Authors:** Michel Wachsmann and Nícolas de Moura  
**Date:** November 22, 2025

## Overview

Parallel implementation of demand estimation models in R and Python, including multinomial logit and random coefficients logit (BLP) models with comprehensive elasticity computation and comparison.

## Project Structure

```
├── main.R                    # Main R estimation script
├── blp.R                     # Custom BLP estimation functions
├── elasticities.R            # Elasticity computation and visualization
├── tables.R                  # LaTeX table generation functions
├── data/
│   └── dataset_ps1.xlsx      # Input data (50 products × 20 markets)
├── python/
│   ├── main.py                   # PyBLP implementation
│   ├── beta_estimates.csv        # Python parameter estimates
│   ├── sigma_estimates.csv       # Python random coefficients
│   ├── elasticities_market1.csv  # Python elasticities
│   └── blp_results.txt           # Full Python output
└── results/
    ├── logit_estimates.tex
    ├── blp_parameter_estimates.tex
    ├── blp_parameter_estimates_comparison.tex
    ├── elasticity_comparison_table.tex
    ├── elasticities_heatmap_logit.png
    ├── elasticities_heatmap_blp.png
    ├── elasticities_heatmap_pyblp.png
    └── full_estimation_environment.RData
```

## Data

- 50 products across 20 markets (1,000 observations)
- Variables: `share`, `price`, `x1-x4` (characteristics), `iv1-iv6` (instruments)
- Instruments: 11 total (1 constant + 4 exogenous characteristics + 6 external IVs)

## Models

### Multinomial Logit
- OLS and IV (2SLS) estimation
- Cluster-robust standard errors at market level
- Elasticity computation: own-price (diagonal) and cross-price (off-diagonal)

### Random Coefficients Logit (BLP)

#### R Implementation
- Nested fixed-point algorithm with contraction mapping
- One-step GMM estimation (optimized for speed)
- Random coefficient on x4 (σ₁)
- Monte Carlo integration (500 draws)
- Convergence tolerance: 1e-8

#### Python Implementation (PyBLP)
- Official PyBLP package implementation
- Monte Carlo integration (1000 draws)
- Same random seed for reproducibility
- Validates R implementation

## Usage

### Requirements

#### R Packages
```r
dplyr, tidyr, ggplot2, gmm, sandwich, AER, readxl, stargazer, progress
```

#### Python Packages
```python
pyblp, pandas, numpy
```

### Run

#### R Estimation
```r
source("main.R")
```

#### Python Estimation
```bash
python main.py
```

## Functions

### blp.R
- `blp()`: Main estimation function with one-step GMM
- `get_delta()`: Contraction mapping for share inversion
- `blp_share()`: Predicted shares with random coefficients
- `get_theta1()`: Linear parameter estimation via 2SLS
- `gmm_moments_sigma()`: GMM moment conditions for σ estimation
- `estimate_elasticities()`: BLP elasticity matrix computation

### elasticities.R
- `compute_logit_elasticities()`: Logit own-price and cross-price elasticities
- `compute_blp_elasticities()`: BLP elasticities with random coefficients
- `create_elasticity_heatmap()`: Upper-triangular heatmap visualization
- `create_elasticity_comparison_table()`: Three-way comparison (Logit, BLP R, PyBLP)

### tables.R
- `create_logit_table()`: OLS vs IV comparison with cluster-robust SEs
- `create_blp_table()`: BLP parameter estimates with optional Python comparison
- `create_combined_table()`: Side-by-side R and Python estimates

## Output

- **Tables**: LaTeX format with standard errors, significance stars (*, **, ***), and R vs PyBLP comparison
- **Heatmaps**: Upper-triangular elasticity matrices for Market 1 (300 DPI PNG)
  - Logit elasticities (homogeneous preferences)
  - BLP elasticities from R (preference heterogeneity)
  - PyBLP elasticities from Python (validation)
- **Interpretation**: Comprehensive LaTeX document explaining results and economic intuition

## Key Results

- **BLP Estimates**: σ₁ = 4.598, α = -36.909
- **Validation**: R and PyBLP produce almost identical parameter estimates
- **Elasticities**: BLP median own-price elasticity (-0.50) substantially lower than Logit (-0.26)
- **Economic Interpretation**: Preference heterogeneity reveals lower aggregate price sensitivity; Logit IIA property confounds average preferences with substitution patterns

## References

Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile prices in market equilibrium. *Econometrica*, 63(4), 841-890.

## Contact

- Nícolas de Moura: nicolasgoulartdemoura@gmail.com
- Michel Wachsmann: michel@wachsmann.com
