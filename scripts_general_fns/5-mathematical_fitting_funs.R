# Mathematically related functions, curve fits etc.

# Linear regressions (~plots) ----


# getting regression line and R2 values to put into the standard curve plot 
# Source: https://stackoverflow.com/a/7549819/9049673

# least square fitting - linear regression (typical use: for qPCR standard curve)
lm_std_curve <- function(df, trig = 0)
{
  x = df %>% pull(Quantity) %>% log10 ; y = df %>% pull(CT)
  m <- lm(y ~ x, df)
  lm_eqn(m, trig)
}


# Extract linear regression equation and parameters from the fit by lm
lm_eqn <- function(m, trig = 0){
  
  eq <- substitute(italic(y) == b %.% italic(x)+ a*","~~italic(r)^2~":"~r2, 
                   list(a = format(unname(coef(m)[1]), digits = 4), 
                        b = format(unname(coef(m)[2]), digits = 3), 
                        r2 = format(summary(m)$r.squared, digits = 3)))
  # if(trig == 'coeff') c(round(coef(m)[2], 2), round(summary(m)$r.squared, 2))
  if(trig == 'coeff') tibble(slope = round(coef(m)[2], 2), y_intercept = round(coef(m)[1], 2), r_square = round(summary(m)$r.squared, 3))
  else as.character(as.expression(eq)); 
}



# Exponential fitting ----


#' safe fitting : ignore errors to fit other samples
#' SSasymp provides a self-starting asymptotic function (exponential fit) with decent initialization
#' @returns Asym+(R0-Asym)*exp(-exp(lrc)*input)
#' https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/SSasymp
safe_expfit <- safely(.f = ~ nls(flipped_fraction ~ SSasymp(day, ys, y0, log_alpha), 
                                 data = .x)
)



#' Calculate the t-half from the log_alpha rate const from exponential curve fit within nested data
get_t_half_from_nls_exp_fits <- function(.df)
  
{
  .df %>% 
    
    mutate( # extract parameters from fit, attach to data
      tidied = map(.fit, ~ broom::tidy(.x$result)), # extracting fitting parameters
      augmented = map(.fit, ~ broom::augment(.x$result)), # extrapolating fitting data, for plotting
      # extrapolated = map(.fit, ~ broom::augment(.x$result, newdata = extrapol_tibble)) # extrapolate when fit didn't work
    ) %>% 
    
    # unnest the model parameters
    unnest(tidied) %>% 
    
    # arrange the parameter estimate, std. error and other stuff for each paramameter in each column
    pivot_wider(names_from = term,
                values_from = c(estimate, std.error, statistic, p.value)) %>% 
    
    # produce t1/2 estimates
    mutate(t.half = log(2)* exp(-estimate_log_alpha), 
           std.error_t.half = log(2) * exp(-estimate_log_alpha) * std.error_log_alpha,
           
           t.half.text = str_c( format(t.half, digits = 2), 
                                '+/-', 
                                format(std.error_t.half, digits = 2),
                                sep = ' ')
    ) # using error propagation - https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example
  
}



#' Access fitting parameters, calculate t half and attach augmented data for plotting
#' Fitting simpler expoonential with lm for robust fits - no singular gradients!
#' formula used : ~ lm(log(flipped_fraction) ~ day, data = .x)
#' @param .df dataframe with list column `.fit` that has the lm fit (typically nested)
get_t_half_from_lm_exp_fits <- function(.df)
  
{
  .df %>% 
    
    mutate( # extract parameters from fit, attach to data
      tidied = map(.fit, ~ broom::tidy(.x)), # extracting fitting parameters
      augmented = map(.fit, ~ broom::augment(.x)), # extrapolating fitting data, for plotting
      
    ) %>% 
    
    # unnest the model parameters
    unnest(tidied) %>% 
    
    # arrange the parameter estimate, std. error and other stuff for each paramameter in each column
    pivot_wider(names_from = term,
                values_from = c(estimate, std.error, statistic, p.value)) %>% 
    
    # produce t1/2 estimates
    mutate(t.half = -log(2)/ estimate_day, 
           std.error_t.half = -t.half / estimate_day * std.error_day,
           
           t.half.text = str_c( format(t.half, digits = 2), 
                                '+/-', 
                                format(std.error_t.half, digits = 2),
                                sep = ' ')
    ) # using error propagation - https://en.wikipedia.org/wiki/Propagation_of_uncertainty#Example
  
}

