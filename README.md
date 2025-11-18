# Empirical Industrial Organization - Problem Set 1

**Authors:** Nícolas de Moura and Michel Wachsmann  
**Date:** November 18, 2025

## Overview

This repository contains the implementation of demand estimation models for a problem set in Empirical Industrial Organization. The code estimates both multinomial logit and random coefficients logit (BLP) models, computes price elasticities, and generates publication-ready tables and visualizations.

## Project Structure

```
empirical/
├── main.R                    # Main estimation script (run this)
├── blp.R                     # BLP estimation functions
├── elasticities.R            # Elasticity computation and visualization
├── tables.R                  # LaTeX table generation functions
├── data/
│   └── dataset_ps1.xlsx     # Input data (50 products × 20 markets)
└── results/
    ├── logit_estimates.tex              # OLS vs IV parameter table
    ├── blp_parameter_estimates.tex      # BLP parameter table
    ├── elasticity_comparison_table.tex  # Logit vs BLP elasticities
    ├── elasticities_heatmap_logit.png   # Logit elasticity heatmap
    ├── elasticities_heatmap_blp.png     # BLP elasticity heatmap
    └── full_estimation_environment.RData # Saved workspace
```

## Data Description

- **Products:** 50 products across 20 geographic markets (1,000 observations)
- **Variables:**
  - `share`: Market share for each product-market
  - `price`: Product price (endogenous)
  - `x1-x4`: Product characteristics (x4 has random coefficient)
  - `iv1-iv6`: Instrumental variables for price

## Models Implemented

### 1. Multinomial Logit Model

Estimates demand using both OLS and IV (2SLS) methods:
- Accounts for price endogeneity using instrumental variables
- Computes cluster-robust standard errors at market level
- Generates own-price and cross-price elasticities

**Key Features:**
- IIA (Independence of Irrelevant Alternatives) property
- Constant elasticity patterns across products

### 2. Random Coefficients Logit (BLP) Model

Implements the Berry, Levinsohn, and Pakes (1995) estimator:
- Nested fixed-point algorithm with contraction mapping
- Two-step GMM with optimal weighting matrix
- Random coefficient on characteristic x4 (σ₁)
- Monte Carlo integration (N=1,000 draws)

**Estimation Steps:**
1. Inner loop: Contraction mapping to invert shares → δ
2. Middle loop: Linear GMM to estimate θ₁ (price and characteristics coefficients)
3. Outer loop: Nonlinear GMM to estimate σ (random coefficient std. dev.)

## Usage

### Requirements

Required R packages:
```r
dplyr, tidyr, ggplot2, gmm, sandwich, AER, readxl, stargazer, progress
```

### Running the Code

1. Set working directory to the `empirical/` folder
2. Run the main script:

```r
source("main.R")
```

The script will:
- Load and prepare data
- Estimate multinomial logit (OLS and IV)
- Estimate BLP model with random coefficients
- Compute elasticities for both models
- Generate tables and heatmaps
- Save all results to `results/` directory

### Key Functions

#### BLP Estimation (`blp.R`)
- `blp()`: Main BLP estimation function
- `get_delta()`: Contraction mapping for δ inversion
- `blp_share()`: Predicted market shares with random coefficients
- `get_theta1()`: Linear parameter estimation via GMM

#### Elasticity Analysis (`elasticities.R`)
- `compute_logit_elasticities()`: Logit model elasticities
- `compute_blp_elasticities()`: BLP model elasticities
- `create_elasticity_heatmap()`: Visualization with colorblind-safe palette
- `create_elasticity_comparison_table()`: Summary statistics table

#### Table Generation (`tables.R`)
- `create_logit_table()`: OLS vs IV comparison with cluster-robust SE
- `create_blp_table()`: BLP parameter estimates with significance stars
- `create_combined_table()`: Side-by-side model comparison

## Output Files

### Tables (LaTeX)

1. **logit_estimates.tex**: OLS and IV estimates with cluster-robust standard errors, model fit statistics
2. **blp_parameter_estimates.tex**: BLP coefficients (θ₁) and random coefficient (σ₁) with cluster-robust SE
3. **elasticity_comparison_table.tex**: Min/median/max of own-price and cross-price elasticities for both models

### Visualizations (PNG)

1. **elasticities_heatmap_logit.png**: Upper-triangular heatmap showing elasticity matrix for Market 1 (Logit model)
2. **elasticities_heatmap_blp.png**: Upper-triangular heatmap showing elasticity matrix for Market 1 (BLP model)

Both heatmaps use:
- Colorblind-safe palette (red for negative, blue for positive)
- Dual-scale normalization for own-price vs cross-price elasticities
- Numeric labels for all elasticity values

## Key Results

### Elasticity Patterns

**Multinomial Logit:**
- Constant substitution patterns (IIA property)
- Own-price elasticities depend on price and share: ε_jj = -α·p_j·(1-s_j)
- Cross-price elasticities uniform within columns: ε_jk = α·p_k·s_k

**BLP Model:**
- Flexible substitution patterns (relaxes IIA)
- Heterogeneous elasticities due to random coefficient σ₁
- Generally more elastic than logit (accommodates unobserved preference heterogeneity)

### Computational Notes

- **Convergence tolerance:** 1e-6 for contraction mapping
- **Maximum iterations:** 1,000 for δ inversion
- **Monte Carlo draws:** 1,000 individuals per market for share prediction
- **Optimization:** Brent method for σ estimation (bounded: 0.01 to 10)

## References

Berry, S., Levinsohn, J., & Pakes, A. (1995). Automobile prices in market equilibrium. *Econometrica*, 63(4), 841-890.

## License

Academic use only - Problem Set 1, Empirical Industrial Organization course, 2025.

## Contact

- Nícolas de Moura: nicolasgoulartdemoura@gmail.com
- Michel Wachsmann: michel@wachsmann.com
