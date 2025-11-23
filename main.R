###############################################################################
# Topic: Empirical Industrial Organization - Problem Set 1
# Goal: Estimate Multinomial Logit and BLP models with elasticities
# Keywords: Logit, BLP, Empirical IO, Demand Estimation, Elasticities
# Autor: Nícolas de Moura (<nicolasgoulartdemoura@gmail.com>) and Michel Wachsmann (<michel@wachsmann.com>)
# Date: 2025-11-18
###############################################################################

###############################################################################
# Organize the working environment
###############################################################################

# Clean the working environment
rm(list = ls())
load.lib <- c("dplyr", "tidyr", "data.table", "ggplot2", "readr", "haven", "progress", "gmm",
              "MASS", "numDeriv", "optimx", "nloptr", "stargazer", "modelsummary", "readxl", "AER")
install.lib <- load.lib[!load.lib %in% installed.packages()]
for(lib in install.lib) install.packages(lib, dependencies = TRUE)
sapply(load.lib, require, character = TRUE)
rm(load.lib, install.lib, lib)
gc()

# Set the random seed for reproducibility
set.seed(20251115)

# Set working directory
setwd("C:/! PROJETOS/200 Graduate/2025.2/Empirical Industrial Organization - Problem Set 1/Empirical-Industrial-Organization---Problem-Set-1")

# Source BLP and elasticity functions
source("blp.R")
source("elasticities.R")
source("tables.R")

# Safe Colorblind Palette
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733",
                             "#332288", "#AA4499", "#44AA99", "#999933", 
                             "#882255", "#661100", "#6699CC", "#888888")

###############################################################################
# Load Data
###############################################################################

# Read from data/dataset_ps1.xlsx
data <- readxl::read_excel("data/dataset_ps1.xlsx")
data <- as.data.frame(data)

###############################################################################
# QUESTION 1: MULTINOMIAL LOGIT MODEL
###############################################################################

# Prepare data for logit estimation
data_logit <- data %>%
  group_by(market) %>%
  mutate(
    outside_share = 1 - sum(share),
    log_share_ratio = log(share) - log(outside_share)
  ) %>%
  ungroup()

# 1.4 OLS estimation
ols_model <- lm(
  log_share_ratio ~ price + x1 + x2 + x3 + x4,
  data = data_logit
)
summary(ols_model)

# IV estimation (2SLS)
iv_model <- ivreg(
  log_share_ratio ~ price + x1 + x2 + x3 + x4 |
    x1 + x2 + x3 + x4 + iv1 + iv2 + iv3 + iv4 + iv5 + iv6,
  data = data_logit
)
summary(iv_model)

create_logit_table(ols_model, iv_model, data = data_logit, market = "market", 
                   output_file = "results/logit_estimates.tex")

# 1.5 Own- and Cross-price Elasticities for Market 1
alpha_iv <- -coef(iv_model)["price"]

# Filter data for Market 1
market1_data <- data %>%
  filter(market == 1) %>%
  arrange(product)

# Compute elasticities for logit model
elasticity_matrix_logit <- compute_logit_elasticities(
  alpha = alpha_iv,
  data = market1_data,
  product = "product",
  price = "price",
  share = "share"
)

# Create and save heatmap
gg_elasticities_logit <- create_elasticity_heatmap(
  elasticity_matrix = elasticity_matrix_logit,
  title = "Own- and Cross-Price Elasticities (Logit Model, Market 1)",
  safe_colorblind_palette = safe_colorblind_palette
)

if (!dir.exists("results")) dir.create("results")
ggsave("results/elasticities_heatmap_logit.png", gg_elasticities_logit, 
       width = 12, height = 12, dpi = 300)

###############################################################################
# QUESTION 2: RANDOM COEFFICIENTS LOGIT (BLP) MODEL
###############################################################################

# 2.2 Estimate BLP parameters
blp_params <- blp(data, product = "product", market = "market", price = "price",
                  characteristics = c("x1", "x2", "x3", "x4"),
                  instruments = c("iv1", "iv2", "iv3", "iv4", "iv5", "iv6"),
                  share = "share", N = 1000)

create_blp_table(blp_params, output_file = "results/blp_parameter_estimates.tex")

# Load Python PyBLP estimates for comparison
python_beta <- read.csv("python/beta_estimates.csv")
python_sigma <- read.csv("python/sigma_estimates.csv")

# Python standard errors (from blp_results.txt)
python_se_beta <- c(1.791313, 10.18339, 1.299349, 1.184628, 0.08664314, 0.1958329)
python_se_sigma <- 2.830655

python_estimates <- list(
  beta = python_beta$estimate,
  sigma = python_sigma$estimate[1],
  se_beta = python_se_beta,
  se_sigma = python_se_sigma
)

# Create comparison table with Python
create_blp_table(blp_params, 
                 output_file = "results/blp_parameter_estimates_comparison.tex",
                 python_estimates = python_estimates)

# 2.3 Own- and Cross-price Elasticities for Market 1
elasticity_matrix_blp <- compute_blp_elasticities(
  blp_params = blp_params,
  data = market1_data,
  product = "product",
  market = "market",
  price = "price",
  characteristics = c("x1", "x2", "x3", "x4"),
  indVarying_chars = c("x4"),
  N = 500
)

# Create and save heatmap
gg_elasticities_blp <- create_elasticity_heatmap(
  elasticity_matrix = elasticity_matrix_blp,
  title = "Own- and Cross-Price Elasticities (BLP Model, Market 1)",
  safe_colorblind_palette = safe_colorblind_palette
)

ggsave("results/elasticities_heatmap_blp.png", gg_elasticities_blp, 
       width = 12, height = 12, dpi = 300)

# Load and visualize Python PyBLP elasticities
python_elasticities <- read.csv("python/elasticities_market1.csv", row.names = 1)
elasticity_matrix_pyblp <- as.matrix(python_elasticities)

# Create heatmap for PyBLP elasticities
gg_elasticities_pyblp <- create_elasticity_heatmap(
  elasticity_matrix = elasticity_matrix_pyblp,
  title = "Own- and Cross-Price Elasticities (PyBLP Model, Market 1)",
  safe_colorblind_palette = safe_colorblind_palette
)

ggsave("results/elasticities_heatmap_pyblp.png", gg_elasticities_pyblp, 
       width = 12, height = 12, dpi = 300)

###############################################################################
# COMPARISON: LOGIT VS BLP ELASTICITIES
###############################################################################

# Create comparison table with PyBLP
comparison_results <- create_elasticity_comparison_table(
  elasticity_matrix_logit = elasticity_matrix_logit,
  elasticity_matrix_blp = elasticity_matrix_blp,
  elasticity_matrix_pyblp = elasticity_matrix_pyblp
)

# Write to file
writeLines(comparison_results$table, "results/elasticity_comparison_table.tex")

###############################################################################
# Save workspace
###############################################################################

save.image("results/full_estimation_environment.RData")
