###############################################################################
# Topic: Empirical Industrial Organization - Problem Set 1
# Goal: Estimate BLP model with PyBLP
# Keywords: BLP, PyBLP, Demand Estimation
# Autor: Nícolas de Moura and Michel Wachsmann
# Date: 2025-11-22
###############################################################################

import os
import numpy as np
import pandas as pd
import pyblp
import matplotlib.pyplot as plt
import seaborn as sns

# Set random seed
np.random.seed(20251115)

# Get script directory for relative paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_PATH = os.path.join(SCRIPT_DIR, "..", "data", "dataset_ps1.xlsx")
OUTPUT_DIR = SCRIPT_DIR

os.makedirs(OUTPUT_DIR, exist_ok=True)

###############################################################################
# Load Data
###############################################################################

df = pd.read_excel(DATA_PATH)

# Rename columns to PyBLP format
df = df.rename(columns={
    'market': 'market_ids',
    'product': 'product_ids', 
    'share': 'shares',
    'price': 'prices'
})

# Ensure proper dtypes
df['market_ids'] = df['market_ids'].astype(int)
df['product_ids'] = df['product_ids'].astype(int)

# Add demand instruments columns with proper naming for PyBLP
for i, iv in enumerate(['iv1', 'iv2', 'iv3', 'iv4', 'iv5', 'iv6']):
    df[f'demand_instruments{i}'] = df[iv]

###############################################################################
# BLP Estimation
###############################################################################

# Formulations according to the utility specification:
# U_ijm = -α*p_jm + β_1*x_jm1 + β_2*x_jm2 + β_3*x_jm3 + (β_4 + σ_1*v_i)*x_jm4 + ξ_jm + ε_ijm

# X1: Linear parameters (constant, price, x1, x2, x3, x4)
X1_formula = pyblp.Formulation('1 + prices + x1 + x2 + x3 + x4')

# X2: Random coefficients (only on x4 with individual heterogeneity)
X2_formula = pyblp.Formulation('0 + x4')

# Integration for Monte Carlo simulation (1000 draws like in R)
integration = pyblp.Integration('monte_carlo', size=5000)

# Create problem
# Note: PyBLP automatically uses exogenous variables (x1-x4) as instruments
# We add additional instruments (iv1-iv6) via the demand_instruments columns
problem = pyblp.Problem(
    product_formulations=(X1_formula, X2_formula),
    product_data=df,
    integration=integration
)

# Solve with initial sigma guess
initial_sigma = np.array([[1.0]])

results = problem.solve(
    sigma=initial_sigma,
    method='1s',
    optimization=pyblp.Optimization('l-bfgs-b')
)

###############################################################################
# Save Results
###############################################################################

# Save summary
with open(os.path.join(OUTPUT_DIR, "blp_results.txt"), "w") as f:
    f.write(str(results))

# Extract and save parameters
beta = results.beta
sigma = results.sigma

pd.DataFrame({
    'parameter': ['constant', 'price', 'x1', 'x2', 'x3', 'x4'],
    'estimate': beta.flatten()
}).to_csv(os.path.join(OUTPUT_DIR, "beta_estimates.csv"), index=False)

pd.DataFrame({
    'parameter': ['sigma_x4'],
    'estimate': sigma.flatten()
}).to_csv(os.path.join(OUTPUT_DIR, "sigma_estimates.csv"), index=False)

###############################################################################
# Compute Elasticities for Market 1
###############################################################################

elasticities = results.compute_elasticities()

# Get Market 1 data
market1_indices = df['market_ids'] == 1
n_products_m1 = market1_indices.sum()

# Extract elasticity submatrix for Market 1
start_idx = 0
for m in sorted(df['market_ids'].unique()):
    n_products = (df['market_ids'] == m).sum()
    if m == 1:
        elas_m1 = elasticities[start_idx:start_idx+n_products, start_idx:start_idx+n_products]
        break
    start_idx += n_products

# Save elasticity matrix
market1_products = df.loc[market1_indices, 'product_ids'].values
elas_df = pd.DataFrame(
    elas_m1,
    index=market1_products,
    columns=market1_products
)
elas_df.to_csv(os.path.join(OUTPUT_DIR, "elasticities_market1.csv"))

###############################################################################
# Create Elasticity Heatmap
###############################################################################

# Create upper triangular mask (keep diagonal and upper triangle)
mask = np.tril(np.ones_like(elas_m1, dtype=bool), k=-1)

# Create figure
plt.figure(figsize=(12, 10))

# Separate scaling for own-price (diagonal) and cross-price (off-diagonal)
own_price = np.diag(elas_m1)
cross_price = elas_m1[~np.eye(len(elas_m1), dtype=bool)]

# Create normalized values for coloring
elas_normalized = np.copy(elas_m1)
own_max = np.abs(own_price).max()
cross_max = np.abs(cross_price).max()

for i in range(len(elas_m1)):
    for j in range(len(elas_m1)):
        if i == j:
            # Normalize own-price elasticities
            elas_normalized[i, j] = elas_m1[i, j] / own_max if own_max > 0 else 0
        else:
            # Normalize cross-price elasticities
            elas_normalized[i, j] = elas_m1[i, j] / cross_max if cross_max > 0 else 0

# Apply mask to normalized values
elas_masked = np.ma.masked_where(mask, elas_normalized)

# Create heatmap with colorblind-safe colors (red for negative, blue for positive)
sns.heatmap(
    elas_masked,
    cmap='RdBu_r',
    center=0,
    vmin=-1,
    vmax=1,
    square=True,
    linewidths=0.5,
    linecolor='white',
    cbar_kws={'label': 'Elasticity (rescaled)'},
    annot=np.ma.masked_where(mask, elas_m1),
    fmt='.2f',
    annot_kws={'size': 7}
)

plt.title('Own- and Cross-Price Elasticities (BLP Model, Market 1)', 
          fontsize=12, fontweight='bold', pad=20)
plt.xlabel('Product k', fontsize=10)
plt.ylabel('Product j', fontsize=10)
plt.xticks(rotation=90)
plt.yticks(rotation=0)
plt.tight_layout()

# Save figure
plt.savefig(os.path.join(OUTPUT_DIR, "elasticities_heatmap_blp.png"), 
            dpi=300, bbox_inches='tight')
plt.close()

print("\nEstimation complete. Results saved to:", OUTPUT_DIR)
print("\nBeta estimates:")
print(beta)
print("\nSigma estimates:")
print(sigma)

