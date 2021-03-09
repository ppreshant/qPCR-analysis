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




