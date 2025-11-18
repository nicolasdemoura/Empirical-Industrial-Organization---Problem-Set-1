###############################################################################
# Topic: Empirical Industrial Organization
# Goal: Functions to compute elasticities for Logit and BLP models
# Keywords: Elasticities, Logit, BLP, Demand Estimation
# Autor: Nícolas de Moura (<nicolasgoulartdemoura@gmail.com>) and Michel Wachsmann (<michel@wachsmann.com>)
# Date: 2025-11-18
###############################################################################

# Compute elasticities for multinomial logit model
compute_logit_elasticities <- function(alpha, data, product = "product", price = "price", share = "share") {
    # Number of products
    J <- nrow(data)
    
    # Extract shares and prices
    s <- data[[share]]
    p <- data[[price]]
    
    # Initialize elasticity matrix
    elasticity_matrix <- matrix(NA_real_, nrow = J, ncol = J)
    
    # Fill matrix with own- and cross-price elasticities
    for (k in 1:J) {
        for (j in 1:J) {
            if (j == k) {
                # Own-price elasticity: epsilon_jj = -alpha * p_j * (1 - s_j)
                elasticity_matrix[j, k] <- -alpha * p[j] * (1 - s[j])
            } else {
                # Cross-price elasticity: epsilon_jk = alpha * p_k * s_k
                elasticity_matrix[j, k] <- alpha * p[k] * s[k]
            }
        }
    }
    
    # Add product names as row/column labels
    prod_ids <- data[[product]]
    rownames(elasticity_matrix) <- prod_ids
    colnames(elasticity_matrix) <- prod_ids
    
    return(elasticity_matrix)
}

# Compute elasticities for BLP model
compute_blp_elasticities <- function(blp_params, data, product = "product", market = "market", 
                                     price = "price", characteristics = c("x1", "x2", "x3", "x4"), 
                                     indVarying_chars = c("x4"), N = 1000) {
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

# Create elasticity heatmap
create_elasticity_heatmap <- function(elasticity_matrix, title = "Own- and Cross-Price Elasticities", 
                                      safe_colorblind_palette = c("#88CCEE", "#CC6677", "#DDCC77", "#117733",
                                                                  "#332288", "#AA4499", "#44AA99", "#999933", 
                                                                  "#882255", "#661100", "#6699CC", "#888888")) {
    library(ggplot2)
    library(dplyr)
    
    # Number of products
    J <- nrow(elasticity_matrix)
    
    # Turn matrix into long data frame
    elasticity_df <- data.frame(
        j = rep(1:J, times = J),   # demand index (rows)
        k = rep(1:J, each = J),    # price index (columns)
        elasticity = as.vector(elasticity_matrix)
    )
    
    # Keep upper triangle (including diagonal)
    elasticity_df <- elasticity_df %>%
        mutate(
            elasticity_plot = ifelse(j <= k, elasticity, NA_real_)
        )
    
    # Max positive and most negative (for separate scaling)
    pos_max <- max(elasticity_df$elasticity_plot[elasticity_df$elasticity_plot > 0],
                   na.rm = TRUE)
    neg_min <- min(elasticity_df$elasticity_plot[elasticity_df$elasticity_plot < 0],
                   na.rm = TRUE)  # this is negative
    
    # Rescale:
    # positives in [0, 1], negatives in [-1, 0], zero -> 0
    elasticity_df <- elasticity_df %>%
        mutate(
            scaled_fill = case_when(
                is.na(elasticity_plot) ~ NA_real_,
                elasticity_plot > 0    ~  elasticity_plot / pos_max,
                elasticity_plot < 0    ~  elasticity_plot / abs(neg_min),
                TRUE                   ~  0
            ),
            label = ifelse(is.na(elasticity_plot),
                           "",
                           sprintf("%.2f", elasticity_plot))
        )
    
    # Choose distinct colors for negative and positive
    neg_col <- safe_colorblind_palette[2]  # e.g., reddish
    pos_col <- safe_colorblind_palette[11]  # e.g., bluish
    
    gg_plot <- ggplot(
        elasticity_df,
        aes(x = k, y = j, fill = scaled_fill)
    ) +
        geom_tile(color = "white") +
        geom_text(aes(label = label), size = 2) +
        scale_x_continuous(breaks = 1:J) +
        scale_y_continuous(breaks = 1:J) +
        scale_fill_gradient2(
            low       = neg_col,
            mid       = "white",
            high      = pos_col,
            limits    = c(-1, 1),
            na.value  = "white",
            name      = "Elasticity\n(rescaled)"
        ) +
        coord_fixed() +
        labs(
            title = title,
            x = "Product k",
            y = "Product j") +
        theme_classic() +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
            panel.grid  = element_blank()
        )
    
    return(gg_plot)
}

# Create elasticity comparison table
create_elasticity_comparison_table <- function(elasticity_matrix_logit, elasticity_matrix_blp) {
    # Extract own-price elasticities (diagonal)
    own_logit <- diag(elasticity_matrix_logit)
    own_blp <- diag(elasticity_matrix_blp)
    
    # Extract cross-price elasticities (off-diagonal)
    cross_logit <- elasticity_matrix_logit[row(elasticity_matrix_logit) != col(elasticity_matrix_logit)]
    cross_blp <- elasticity_matrix_blp[row(elasticity_matrix_blp) != col(elasticity_matrix_blp)]
    
    # Compute summary statistics
    logit_stats <- c(
        min(own_logit), median(own_logit), max(own_logit),
        min(cross_logit), median(cross_logit), max(cross_logit)
    )
    
    blp_stats <- c(
        min(own_blp), median(own_blp), max(own_blp),
        min(cross_blp), median(cross_blp), max(cross_blp)
    )
    
    # Create LaTeX table
    comparison_table <- paste0(
        "\\begin{table}[H] \\centering \n",
        "  \\caption{Summary of Elasticities in Market 1: Logit vs.\\ BLP} \n",
        "  \\label{tab:elasticity_comparison} \n",
        "\\begin{tabular}{@{\\extracolsep{5pt}}lcccccc} \n",
        "\\\\[-1.8ex]\\hline \n",
        "\\hline \\\\[-1.8ex] \n",
        " & \\multicolumn{3}{c}{Own-Price} & \\multicolumn{3}{c}{Cross-Price} \\\\ \n",
        "\\cline{2-4} \\cline{5-7}\n",
        "Model & Min & Median & Max & Min & Median & Max \\\\ \n",
        "\\hline \\\\[-1.8ex] \n",
        sprintf("Multinomial Logit & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\ \n",
                logit_stats[1], logit_stats[2], logit_stats[3],
                logit_stats[4], logit_stats[5], logit_stats[6]),
        sprintf("BLP & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ & $%.2f$ \\\\ \n",
                blp_stats[1], blp_stats[2], blp_stats[3],
                blp_stats[4], blp_stats[5], blp_stats[6]),
        "\\hline \n",
        "\\hline \\\\[-1.8ex] \n",
        "\\end{tabular} \n",
        "\\end{table}"
    )
    
    return(list(
        table = comparison_table,
        logit_stats = logit_stats,
        blp_stats = blp_stats
    ))
}
