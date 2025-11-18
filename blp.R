###############################################################################
# Topic: Empirical Industrial Organization
# Goal: Create functions to do our own BLP estimator
# Keywords: BLP, Empirical IO, Demand Estimation, Structural
# Autor: Nícolas de Moura (<nicolasgoulartdemoura@gmail.com>)
# Date: 2025-11-15
###############################################################################

# Get initial_delta function
get_initial_delta <- function(alpha, beta, data, product = "product", market = "market", price = "price", characteristics = NULL) {
    # Compute delta_j as alpha * price_j + X_j * beta
    delta <- alpha * data[[price]] + as.matrix(data[, characteristics]) %*% beta
    return(as.vector(delta))
}

# Get initial_alpha and beta function via multinomial logit
get_initial_alpha_beta <- function(data, product = "product", market = "market", price = "price", characteristics = NULL, share = "share") {
    # Prepare data
    data$outside_share <- 1 - ave(data[[share]], data[[market]], FUN = sum)
    data$log_share_ratio <- log(data[[share]]) - log(data$outside_share)
    
    # Create formula for linear regression
    formula_str <- paste("log_share_ratio ~", paste(c(price, characteristics), collapse = " + "))
    formula <- as.formula(formula_str)
    
    # Run linear regression
    model <- lm(formula, data = data)
    coefs <- coef(model)
    
    alpha <- coefs[price]
    beta <- coefs[characteristics]
    
    return(list(alpha = alpha, beta = beta))
}


# Estimate Share Function
blp_share <- function(delta, sigma, data, indVarying_chars = NULL, N = 100, product = "product", market = "market", price = "price") {
    # Number of products and individuals
    M <- length(unique(data[[market]]))
    I <- N

    # Draw shocks v iid from N(0,1)
    v <- matrix(rnorm(I * length(sigma)), nrow = I, ncol = length(sigma))

    # Initialize shares vector
    shares <- numeric(nrow(data))
    
    # Loop over each market
    for (m in unique(data[[market]])) {
        market_data <- data[data[[market]] == m, ]
        J_m <- nrow(market_data)
        market_indices <- which(data[[market]] == m)
        
        # Compute mu_{ij} for market m
        mu <- matrix(0, nrow = I, ncol = J_m)
        for (k in seq_along(indVarying_chars)) {
            char_k <- indVarying_chars[k]
            sigma_k <- sigma[k]
            mu <- mu + sigma_k * v[, k] %*% t(market_data[[char_k]])
        }
        
        # Compute utility u_{ij} = delta_j + mu_{ij}
        delta_market <- delta[market_indices]
        utility <- matrix(rep(delta_market, each = I), nrow = I, ncol = J_m) + mu
        exp_utility <- exp(utility)
        sum_exp_utility <- rowSums(exp_utility) + 1  # +1 for outside good
        market_shares <- colMeans(exp_utility / sum_exp_utility) 
        shares[market_indices] <- market_shares
    }
    
    return(shares)
}

# get_delta function via contraction mapping ln(observed_shares) - ln(estimated_shares)
get_delta <- function(delta_init, sigma, data, indVarying_chars, N, observed_shares, tol = 1e-6, max_iter = 1000, silent = FALSE) {
    delta <- delta_init
    if (!silent) {
        # Using progress package to show progress bar
        pb <- progress::progress_bar$new(total = max_iter, format = "  Contraction Mapping [:bar] :percent in :elapsed")
    }
    for (iter in 1:max_iter) {
        if (!silent) pb$tick()
        estimated_shares <- blp_share(delta, sigma, data, indVarying_chars, N)
        share_diff <- log(observed_shares) - log(estimated_shares)
        delta_new <- delta + share_diff
        if (max(abs(delta_new - delta)) < tol) {
            break
        }
        delta <- delta_new
    }
    return(delta)
}

# get_theta1 function to estimate linear parameters via IV/GMM using gmm package
get_theta1 <- function(delta, data, characteristics, instruments, price = "price", market = "market") {
    library(gmm)
    library(sandwich)
    
    # Prepare moment conditions
    g <- function(theta, x) {
        X <- x$X
        Z <- x$Z
        delta <- x$delta
        residuals <- delta - X %*% theta
        moments <- Z * as.vector(residuals)
        return(moments)
    }
    
    X <- as.matrix(data[, c(price, characteristics)])
    Z <- as.matrix(data[, instruments])
    
    # Initial values from 2SLS (overidentified system)
    PZ <- Z %*% solve(t(Z) %*% Z) %*% t(Z)
    init_theta <- solve(t(X) %*% PZ %*% X) %*% t(X) %*% PZ %*% delta
    
    # Estimate using gmm package with two-step optimal weighting
    gmm_result <- gmm(g, x = list(X = X, Z = Z, delta = delta), 
                      t0 = as.vector(init_theta), 
                      type = "twoStep",
                      wmatrix = "optimal",
                      vcov = "iid",
                      centeredVcov = TRUE,
                      data = data)
    
    # Compute cluster robust vcov
    vcov_cluster <- vcovCL(gmm_result, cluster = data[[market]])
    
    return(list(theta = as.vector(coef(gmm_result)), 
                gmm_result = gmm_result,
                vcov_cluster = vcov_cluster))
}


# GMM moment conditions for nonlinear parameters (sigma)
gmm_moments_sigma <- function(sigma, x) {
    data <- x$data
    characteristics <- x$characteristics
    instruments <- x$instruments
    price <- x$price
    share <- x$share
    indVarying_chars <- x$indVarying_chars
    N <- x$N
    market <- x$market
    
    # Step 1: Get initial delta using simple logit
    init_params <- get_initial_alpha_beta(data, characteristics = characteristics, price = price, share = share)
    alpha_init <- init_params$alpha
    beta_init <- init_params$beta
    delta_init <- get_initial_delta(alpha_init, beta_init, data, characteristics = characteristics, price = price)
    
    # Step 2: Get delta via contraction mapping given sigma
    observed_shares <- data[[share]]
    delta <- get_delta(delta_init, sigma, data, indVarying_chars = indVarying_chars, N = N, observed_shares = observed_shares, silent = TRUE)
    
    # Step 3: Estimate theta1 (alpha, beta) via linear IV/GMM given delta and sigma
    theta1_result <- get_theta1(delta, data, characteristics, instruments, price = price, market = market)
    theta1 <- theta1_result$theta
    
    # Step 4: Compute residuals
    X <- as.matrix(data[, c(price, characteristics)])
    Z <- as.matrix(data[, instruments])
    predicted_delta <- X %*% theta1
    xi <- delta - predicted_delta
    
    # Step 5: Return moment conditions
    moments <- Z * as.vector(xi)
    
    return(moments)
}

# Create blp function to estimate BLP parameters using nested GMM with gmm package
blp <- function(data, product = "product", market = "market", price = "price", 
                characteristics = c("x1", "x2", "x3", "x4"), 
                instruments = c("iv1", "iv2", "iv3", "iv4", "iv5", "iv6"), 
                share = "share", indVarying_chars = c("x4"), N = 1000) {
    
    library(gmm)
    library(sandwich)
    
    x_data <- list(
        data = data,
        characteristics = characteristics,
        instruments = instruments,
        price = price,
        share = share,
        indVarying_chars = indVarying_chars,
        N = N,
        market = market
    )
    
    gmm_sigma <- gmm(
        g = gmm_moments_sigma,
        x = x_data,
        t0 = 0.5,
        type = "twoStep",
        wmatrix = "optimal",
        method = "Brent",
        lower = 0.01,
        upper = 10,
        vcov = "iid",
        data = data, 
        control = list(trace = 1)
    )
    
    # Compute cluster robust vcov for sigma
    vcov_sigma_cluster <- vcovCL(gmm_sigma, cluster = data[[market]])
    
    sigma_est <- as.numeric(coef(gmm_sigma))
    se_sigma <- sqrt(as.numeric(vcov_sigma_cluster))
    
    # Step 2: Get final delta using estimated sigma
    init_params <- get_initial_alpha_beta(data, characteristics = characteristics, price = price, share = share)
    alpha_init <- init_params$alpha
    beta_init <- init_params$beta
    delta_init <- get_initial_delta(alpha_init, beta_init, data, characteristics = characteristics, price = price)
    observed_shares <- data[[share]]
    delta_final <- get_delta(delta_init, sigma_est, data, indVarying_chars = indVarying_chars, N = N, observed_shares = observed_shares, silent = FALSE)
    
    # Step 3: Estimate theta1 (alpha, beta) via gmm package
    theta1_result <- get_theta1(delta_final, data, characteristics, instruments, price = price, market = market)
    theta1_est <- theta1_result$theta
    gmm_theta1 <- theta1_result$gmm_result
    vcov_theta1_cluster <- theta1_result$vcov_cluster
    se_theta1 <- sqrt(diag(vcov_theta1_cluster))
    X <- as.matrix(data[, c(price, characteristics)])

    # Compute residuals
    xi <- delta_final - X %*% theta1_est
    
    # Hansen J-test for theta1
    J_stat_theta1 <- specTest(gmm_theta1)$test[1]
    J_df_theta1 <- specTest(gmm_theta1)$test[2]
    J_pvalue_theta1 <- specTest(gmm_theta1)$test[3]
    
    # Hansen J-test for sigma
    J_stat_sigma <- specTest(gmm_sigma)$test[1]
    J_df_sigma <- specTest(gmm_sigma)$test[2]
    J_pvalue_sigma <- specTest(gmm_sigma)$test[3]
    
    return(list(
        sigma = sigma_est,
        se_sigma = se_sigma,
        theta1 = theta1_est,
        se_theta1 = se_theta1,
        delta = delta_final,
        xi = as.vector(xi),
        vcov_theta1 = vcov_theta1_cluster,
        vcov_sigma = as.numeric(vcov_sigma_cluster),
        gmm_theta1 = gmm_theta1,
        gmm_sigma = gmm_sigma,
        J_stat_theta1 = as.numeric(J_stat_theta1),
        J_df_theta1 = as.numeric(J_df_theta1),
        J_pvalue_theta1 = as.numeric(J_pvalue_theta1),
        J_stat_sigma = as.numeric(J_stat_sigma),
        J_df_sigma = as.numeric(J_df_sigma),
        J_pvalue_sigma = as.numeric(J_pvalue_sigma)
    ))
}

# Estimate elasticities function
estimate_elasticities <- function(blp_params, data, product = "product", market = "market", price = "price", characteristics = c("x1", "x2", "x3", "x4"), indVarying_chars = c("x4"), N = 1000) {
    sigma <- blp_params$sigma
    delta <- blp_params$delta
    alpha <- blp_params$theta1[1]
    
    # Number of products and individuals
    M <- length(unique(data[[market]]))
    I <- N
    
    # Draw shocks v iid from N(0,1)
    v <- matrix(rnorm(I * length(sigma)), nrow = I, ncol = length(sigma))
    
    # Initialize elasticity matrix
    elasticity_matrix <- matrix(0, nrow = nrow(data), ncol = nrow(data))
    rownames(elasticity_matrix) <- data[[product]]
    colnames(elasticity_matrix) <- data[[product]]
    
    # Loop over each market
    for (m in unique(data[[market]])) {
        market_data <- data[data[[market]] == m, ]
        J_m <- nrow(market_data)
        market_indices <- which(data[[market]] == m)
        
        # Compute mu_{ij} for market m
        mu <- matrix(0, nrow = I, ncol = J_m)
        for (k in seq_along(indVarying_chars)) {
            char_k <- indVarying_chars[k]
            sigma_k <- sigma[k]
            mu <- mu + sigma_k * v[, k] %*% t(market_data[[char_k]])
        }
        
        # Compute utility u_{ij} = delta_j + mu_{ij}
        delta_market <- delta[market_indices]
        utility <- matrix(rep(delta_market, each = I), nrow = I, ncol = J_m) + mu
        exp_utility <- exp(utility)
        sum_exp_utility <- rowSums(exp_utility) + 1 # +1 for outside good
        choice_probabilities <- exp_utility / sum_exp_utility
        
        # Compute shares
        shares <- colMeans(choice_probabilities)
        
        # Compute elasticities
        for (j in 1:J_m) {
            for (k in 1:J_m) {
                if (j == k) {
                    # Own-price elasticity (negative)
                    elasticity_matrix[market_data[[product]][j], market_data[[product]][k]] <- 
                        alpha * mean(choice_probabilities[, j] * (1 - choice_probabilities[, j])) * (market_data[[price]][j] / shares[j])
                } else {
                    # Cross-price elasticity (positive)
                    elasticity_matrix[market_data[[product]][j], market_data[[product]][k]] <- 
                        -alpha * mean(choice_probabilities[, j] * choice_probabilities[, k]) * (market_data[[price]][k] / shares[j])
                }
            }
        }
    }

    return(elasticity_matrix)
}
