###############################################################################
# Topic: Empirical Industrial Organization
# Goal: Functions to create parameter estimate tables
# Keywords: Tables, Logit, BLP, Parameter Estimates
# Autor: Nícolas de Moura (<nicolasgoulartdemoura@gmail.com>) and Michel Wachsmann (<michel@wachsmann.com>)
# Date: 2025-11-18
###############################################################################

# Create parameter estimates table for multinomial logit (OLS vs IV)
create_logit_table <- function(ols_model, iv_model, data, market = "market", output_file = NULL) {
    library(sandwich)
    library(lmtest)
    
    # Compute cluster-robust standard errors
    vcov_ols_cl <- vcovCL(ols_model, cluster = data[[market]])
    vcov_iv_cl <- vcovCL(iv_model, cluster = data[[market]])
    
    # Extract coefficients
    coef_ols <- coef(ols_model)
    coef_iv <- coef(iv_model)
    
    # Extract cluster-robust standard errors
    se_ols_cl <- sqrt(diag(vcov_ols_cl))
    se_iv_cl <- sqrt(diag(vcov_iv_cl))
    
    # Compute t-statistics and p-values
    t_ols <- coef_ols / se_ols_cl
    t_iv <- coef_iv / se_iv_cl
    p_ols <- 2 * pt(-abs(t_ols), df = nrow(data) - length(coef_ols))
    p_iv <- 2 * pt(-abs(t_iv), df = nrow(data) - length(coef_iv))
    
    # Function to add significance stars
    add_stars <- function(p) {
        ifelse(p < 0.01, "***",
               ifelse(p < 0.05, "**",
                      ifelse(p < 0.1, "*", "")))
    }
    
    # Model statistics
    n_obs <- nobs(ols_model)
    r2_ols <- summary(ols_model)$r.squared
    r2_iv <- summary(iv_model)$r.squared
    adj_r2_ols <- summary(ols_model)$adj.r.squared
    adj_r2_iv <- summary(iv_model)$adj.r.squared
    rse_ols <- summary(ols_model)$sigma
    rse_iv <- summary(iv_model)$sigma
    df_ols <- ols_model$df.residual
    df_iv <- iv_model$df.residual
    
    # Build LaTeX table
    param_names <- c("Price", "X1", "X2", "X3", "X4", "Constant")
    
    latex_table <- paste0(
        "\\begin{table}[H] \\centering \n",
        "  \\label{tab:ols_iv_logshare} \n",
        "\\begin{tabular}{@{\\extracolsep{5pt}}lcc} \n",
        "\\\\[-1.8ex]\\hline \n",
        "\\hline \\\\[-1.8ex] \n",
        " & \\textbf{OLS} & \\textbf{IV} \\\\ \n",
        "\\hline \\\\[-1.8ex] \n"
    )
    
    for (i in seq_along(param_names)) {
        latex_table <- paste0(latex_table,
            sprintf("%s & $%.3f^{%s}$ & $%.3f^{%s}$ \\\\ \n",
                    param_names[i],
                    coef_ols[i],
                    add_stars(p_ols[i]),
                    coef_iv[i],
                    add_stars(p_iv[i])),
            sprintf(" & (%.3f) & (%.3f) \\\\ \n", se_ols_cl[i], se_iv_cl[i]),
            " & & \\\\ \n"
        )
    }
    
    latex_table <- paste0(latex_table,
        "\\hline \\\\[-1.8ex] \n",
        sprintf("Observations & %d & %d \\\\ \n", n_obs, n_obs),
        sprintf("R$^{2}$ & %.3f & %.3f \\\\ \n", r2_ols, r2_iv),
        sprintf("Adjusted R$^{2}$ & %.3f & %.3f \\\\ \n", adj_r2_ols, adj_r2_iv),
        sprintf("Residual Std. Error & %.3f (df = %d) & %.3f (df = %d) \\\\ \n", 
                rse_ols, df_ols, rse_iv, df_iv),
        "\\hline \n",
        "\\hline \\\\[-1.8ex] \n",
        "\\textit{Note:} & \\multicolumn{2}{r}{$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01} \\\\ \n",
        "\\end{tabular} \n",
        "  \\caption{OLS and IV estimates for $\\log(\\text{share})$ ratio} \n",
        "\\end{table}"
    )
    
    if (!is.null(output_file)) {
        writeLines(latex_table, output_file)
    }
    
    return(latex_table)
}

# Create parameter estimates table for BLP
create_blp_table <- function(blp_params, output_file = NULL) {
    # Function to add significance stars
    add_stars <- function(coef, se) {
        t_stat <- coef / se
        p_val <- 2 * pt(-abs(t_stat), df = 1000)
        ifelse(p_val < 0.01, "***",
               ifelse(p_val < 0.05, "**",
                      ifelse(p_val < 0.1, "*", "")))
    }
    
    # Extract coefficients and standard errors
    coefs <- blp_params$theta1
    ses <- blp_params$se_theta1
    param_names <- c("Price", "X1", "X2", "X3", "X4", "Constant")
    
    # Build LaTeX table
    latex_table <- paste0(
        "\\begin{table}[H] \\centering \n",
        "  \\label{tab:blp_estimates} \n",
        "\\begin{tabular}{@{\\extracolsep{5pt}}lc} \n",
        "\\\\[-1.8ex]\\hline \n",
        "\\hline \\\\[-1.8ex] \n",
        " & \\textbf{BLP} \\\\ \n",
        "\\hline \\\\[-1.8ex] \n"
    )
    
    for (i in 1:6) {
        latex_table <- paste0(latex_table,
            sprintf("%s & $%.3f^{%s}$ \\\\ \n",
                    param_names[i],
                    coefs[i],
                    add_stars(coefs[i], ses[i])),
            sprintf(" & (%.3f) \\\\ \n", ses[i]),
            " & \\\\ \n"
        )
    }
    
    # Add sigma_1 separately
    latex_table <- paste0(latex_table,
        sprintf("$\\\\sigma_1$ & $%.3f^{%s}$ \\\\ \n",
                blp_params$sigma,
                add_stars(blp_params$sigma, blp_params$se_sigma)),
        sprintf(" & (%.3f) \\\\ \n", blp_params$se_sigma),
        " & \\\\ \n"
    )
    
    latex_table <- paste0(latex_table,
        "\\hline \n",
        "\\hline \\\\[-1.8ex] \n",
        "\\textit{Note:} & \\multicolumn{1}{r}{$^{*}$p$<$0.1; $^{**}$p$<$0.05; $^{***}$p$<$0.01} \\\\ \n",
        "\\end{tabular} \n",
        "  \\caption{BLP Parameter Estimates} \n",
        "\\end{table}"
    )
    
    if (!is.null(output_file)) {
        writeLines(latex_table, output_file)
    }
    
    return(latex_table)
}

# Create comparison table for both models side by side
create_combined_table <- function(ols_model, iv_model, blp_params, output_file = NULL) {
    library(modelsummary)
    
    ols_coefs <- coef(ols_model)
    ols_se <- sqrt(diag(vcov(ols_model)))
    
    iv_coefs <- coef(iv_model)
    iv_se <- sqrt(diag(vcov(iv_model)))
    
    blp_coefs <- c(blp_params$theta1, blp_params$sigma)
    blp_se <- c(blp_params$se_theta1, blp_params$se_sigma)
    
    param_names <- c("Price", "x1", "x2", "x3", "x4", "Constant", "sigma1")
    
    results_df <- data.frame(
        Parameter = param_names,
        OLS_Est = c(ols_coefs, NA),
        OLS_SE = c(paste0("(", sprintf("%.3f", ols_se), ")"), ""),
        IV_Est = c(iv_coefs, NA),
        IV_SE = c(paste0("(", sprintf("%.3f", iv_se), ")"), ""),
        BLP_Est = c(blp_coefs[-6], NA, blp_coefs[6]),
        BLP_SE = c(paste0("(", sprintf("%.3f", blp_se[-6]), ")"), "", 
                   paste0("(", sprintf("%.3f", blp_se[6]), ")"))
    )
    
    latex_table <- paste0(
        "\\begin{table}[H] \\centering \n",
        "  \\caption{Parameter Estimates: Logit vs BLP} \n",
        "  \\label{tab:combined_estimates} \n",
        "\\begin{tabular}{@{\\extracolsep{5pt}}lcccccc} \n",
        "\\\\[-1.8ex]\\hline \n",
        "\\hline \\\\[-1.8ex] \n",
        " & \\multicolumn{2}{c}{OLS} & \\multicolumn{2}{c}{IV (2SLS)} & \\multicolumn{2}{c}{BLP} \\\\ \n",
        "\\cline{2-3} \\cline{4-5} \\cline{6-7} \n",
        "Parameter & Est. & SE & Est. & SE & Est. & SE \\\\ \n",
        "\\hline \\\\[-1.8ex] \n"
    )
    
    for (i in 1:nrow(results_df)) {
        if (!is.na(results_df$OLS_Est[i])) {
            latex_table <- paste0(latex_table,
                sprintf("%s & %.3f & %.3f & %.3f & %.3f & ",
                        results_df$Parameter[i],
                        results_df$OLS_Est[i],
                        as.numeric(gsub("[()]", "", results_df$OLS_SE[i])),
                        results_df$IV_Est[i],
                        as.numeric(gsub("[()]", "", results_df$IV_SE[i]))))
            if (!is.na(results_df$BLP_Est[i])) {
                latex_table <- paste0(latex_table,
                    sprintf("%.3f & %.3f \\\\ \n",
                            results_df$BLP_Est[i],
                            as.numeric(gsub("[()]", "", results_df$BLP_SE[i]))))
            } else {
                latex_table <- paste0(latex_table, " & \\\\ \n")
            }
        } else if (!is.na(results_df$BLP_Est[i])) {
            latex_table <- paste0(latex_table,
                sprintf("%s & & & & & %.3f & %.3f \\\\ \n",
                        results_df$Parameter[i],
                        results_df$BLP_Est[i],
                        as.numeric(gsub("[()]", "", results_df$BLP_SE[i]))))
        }
    }
    
    latex_table <- paste0(latex_table,
        "\\hline \n",
        "\\hline \\\\[-1.8ex] \n",
        "\\end{tabular} \n",
        "\\end{table}"
    )
    
    if (!is.null(output_file)) {
        writeLines(latex_table, output_file)
    }
    
    return(latex_table)
}
