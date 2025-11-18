# Empirical Industrial Organization - Problem Set 1

**Authors:** Michel Wachsmann and Nícolas de Moura  
**Date:** November 18, 2025

## Overview

Implementation of demand estimation models including multinomial logit and random coefficients logit (BLP) models with elasticity computation.

## Project Structure

```
empirical/
├── main.R                    # Main estimation script
├── blp.R                     # BLP estimation functions
├── elasticities.R            # Elasticity computation and visualization
├── tables.R                  # Table generation functions
├── data/
│   └── dataset_ps1.xlsx     # Input data (50 products × 20 markets)
└── results/
    ├── logit_estimates.tex
    ├── blp_parameter_estimates.tex
    ├── elasticity_comparison_table.tex
    ├── elasticities_heatmap_logit.png
    ├── elasticities_heatmap_blp.png
    └── full_estimation_environment.RData
```

## Data

- 50 products across 20 markets (1,000 observations)
- Variables: `share`, `price`, `x1-x4` (characteristics), `iv1-iv6` (instruments)

## Models

### Multinomial Logit
- OLS and IV (2SLS) estimation
- Cluster-robust standard errors at market level
- Elasticity computation

### Random Coefficients Logit (BLP)
- Nested fixed-point algorithm with contraction mapping
- Two-step GMM estimation
- Random coefficient on x4 (σ₁)
- Monte Carlo integration (1,000 draws)

## Usage

### Requirements
```r
dplyr, tidyr, ggplot2, gmm, sandwich, AER, readxl, stargazer, progress
```

### Run
```r
source("main.R")
```

## Functions

### blp.R
- `blp()`: Main estimation function
- `get_delta()`: Contraction mapping for share inversion
- `blp_share()`: Predicted shares with random coefficients
- `get_theta1()`: Linear parameter estimation

### elasticities.R
- `compute_logit_elasticities()`: Logit elasticities
- `compute_blp_elasticities()`: BLP elasticities
- `create_elasticity_heatmap()`: Elasticity visualization
- `create_elasticity_comparison_table()`: Summary statistics

### tables.R
- `create_logit_table()`: OLS vs IV comparison
- `create_blp_table()`: BLP parameter estimates
- `create_combined_table()`: Side-by-side comparison

## Output

- **Tables**: LaTeX format with standard errors and significance levels
- **Heatmaps**: Upper-triangular elasticity matrices for Market 1

## References

Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile prices in market equilibrium. *Econometrica*, 63(4), 841-890.

## Contact

- Nícolas de Moura: nicolasgoulartdemoura@gmail.com
- Michel Wachsmann: michel@wachsmann.com
