calc_skewness = function(vec) {
  coeff = (3 * (mean(vec) - median(vec))) / sd(vec)
  return(coeff)
}

get_expected_value = function(x, px) {
  if (length(x) != length(px)) {
    stop("lengths are different. They must be same")
  }
  
  return(x %*% px)
}

get_expected_variance = function(x, px) {
  mu_x = get_expected_value(x, px)
  variance = sum(((x - mu_x)**2) * px)
  writeClipboard("" + variance)
  return(variance)
}


get_t_percentile = function(x, x_bar, s, n) {
  t_score = ( (x-x_bar) / (s / sqrt(n))   )
  return (pt(t_score, n-1))
}

calculateSDConfidenceInterval <- function(n, s, confidenceLevel = 0.95) {
  # Degrees of freedom
  df <- n - 1
  
  # Alpha level (1 - confidence level)
  alpha <- 1 - confidenceLevel
  
  # Chi-squared critical values for the lower and upper tails
  chiSquareLower <- qchisq(alpha / 2, df, lower.tail = TRUE)
  chiSquareUpper <- qchisq(1 - alpha / 2, df, lower.tail = TRUE)
  
  # Calculate the lower and upper bounds of the confidence interval for variance
  lowerBoundVariance <- ((n - 1) * s^2) / chiSquareUpper
  upperBoundVariance <- ((n - 1) * s^2) / chiSquareLower
  
  # Calculate the square root of variance bounds to get the confidence interval for standard deviation
  lowerBoundSD <- sqrt(lowerBoundVariance)
  upperBoundSD <- sqrt(upperBoundVariance)
  
  # Return the confidence interval
  return(c(lower = lowerBoundSD, upper = upperBoundSD))
}



easy_pie = function(data_vec, labels_vec, colors_vec) {
  new_labels_vec <- character(length = length(labels_vec))
  for (i in seq_along(data_vec)) {
    percent = round( (data_vec[i]/sum(data_vec)) * 100, digits=1)
                     
    new_labels_vec[i] = paste(
              paste(labels_vec[i], " "),
              paste(percent, "%")
    )
    
  }
  pie(data_vec, labels = new_labels_vec, col = colors_vec)
}

# requires qcc library
easy_pareto = function(data_vec, labels_vec, xlab="", ylab="") {
  library(qcc)
  
  if (length(data_vec) != length(labels_vec)) {
    stop("vectors are of unequal length. Cannot make pareto chart")
  }
  
  names(data_vec) = labels_vec
  
  pareto.chart(
                data_vec,
                cumperc = seq(0, 100, by = 20),
                ylab2 = "Cumulative Percentage",
                xlab=xlab,
                ylab=ylab
              )
}


get_binomial_distribution = function (n, p) {
  distribution_vec = c()
  for (k in 0: n) {
    mass_function = dbinom(k, n, p)
    distribution_vec = c(distribution_vec, mass_function)
  }
  return(distribution_vec)
}

# 

get_cumulative_binomial_distribution = function (n, p) {
  distribution_vec = c()
  for (k in 0: n) {
    mass_function = pbinom(k, n, p)
    distribution_vec = c(distribution_vec, mass_function)
  }
  return(distribution_vec)
}

format_double <- function(x) {
  # Function to format double values
  format_value <- function(value) {
    if (abs(value) < 0.001 || abs(value) >= 1000) {
      # Use scientific notation with 3 decimal places
      return(sprintf("%.3e", value))
    } else {
      # Use default formatting for other values
      return(as.character(round(value, digits=3)))
    }
  }
  
  # Apply formatting to each value in the vector
  formatted_values <- sapply(x, format_value)
  
  # Join formatted values into a comma-separated string
  formatted_string <- paste(formatted_values, collapse = ", ")
  
  # Print the formatted string
  cat(formatted_string, "\n")
}


pretty_print_and_copy <- function(vector) {
  if (!requireNamespace("clipr", quietly = TRUE)) {
    install.packages("clipr")
  }
  library(clipr)
  # Convert vector to character and create comma-separated string
  vector_str <- paste(format_double(vector), collapse = ", ")
  write_clip(vector_str)
  cat(vector_str, "\n")
}


Modes <- function(x) {
  ux <- unique(x)
  tab <- tabulate(match(x, ux))
  ux[tab == max(tab)]
}

plot_exponential = function(lambda, num_data_points = 20) {
  x = 0:num_data_points
  y = lambda * exp(-lambda * x)
  
  # Plot the function
  plot(x, y, type = "l", col = "blue", lwd = 2, 
       main = paste("Exponential distribution with lambda = ", 
                    lambda),
       xlab = "x", ylab = "f(x)")
  
  # Add gridlines for better visualization
  grid()
}


get_z_score = function(x, mu, sigma) {
  return((x-mu)/sigma)
}

get_confidence_interval = function(x_bar, n, sigma, confidence_level = 0.95) {
  sample_sigma = sigma / sqrt(n)
  z_star = abs(qnorm((1 - confidence_level)/2, mean = 0, sd = 1))
  interval = c(x_bar - z_star*sample_sigma, x_bar + z_star*sample_sigma)
  return(interval)
}

get_confidence_interval_difference_means = function(
    x_bar_1, 
    x_bar_2, 
    n_1,
    n_2,
    sigma_1,
    sigma_2,
    confidence_level = 0.95
) {
  standard_error = sqrt( (sigma_1**2/n_1) + (sigma_2**2/n_2)    )
  
  alpha = 1 - confidence_level
  z_star = qnorm(1 - (alpha/2), mean = 0, sd = 1)
  
  lower = abs(x_bar_1 - x_bar_2) - z_star*standard_error
  upper = abs(x_bar_1 - x_bar_2) + z_star*standard_error
  
  return(c(lower = lower, upper = upper))
  
}

get_confidence_interval_difference_means_s = function(
    x_bar_1, 
    x_bar_2, 
    n_1,
    n_2,
    s_1,
    s_2,
    confidence_level = 0.95
) {
  standard_error = sqrt( (s_1**2/n_1) + (s_2**2/n_2)    )
  
  alpha = 1 - confidence_level
  t_star = qt(1 - (alpha/2), df=(n_1 - 1) + (n_2 - 1))
  
  lower = abs(x_bar_1 - x_bar_2) - t_star*standard_error
  upper = abs(x_bar_1 - x_bar_2) + t_star*standard_error
  
  return(c(lower = lower, upper = upper))
  
}


get_difference_means_test_statistic_z = function(
    x_bar_1, 
    x_bar_2, 
    n_1,
    n_2,
    s_1,
    s_2
) {
  standard_error = sqrt( (s_1**2/n_1) + (s_2**2/n_2)    )
  
  z_ts =  ( (x_bar_1 - x_bar_2) - 0) / standard_error
  return ( z_ts )
  
}



get_confidence_interval_proportion = function(p_hat, n, confidence_level = 0.95) {
  sigma = sqrt(p_hat*(1-p_hat))
  sample_sigma = sigma / sqrt(n)
  
  z_star = abs(qnorm((1 - confidence_level)/2, mean = 0, sd = 1))
  interval = c(p_hat - z_star*sample_sigma, p_hat + z_star*sample_sigma)
  print(paste("(", round(interval[1], 4), ",", round(interval[2], 4), ")"))
  return(interval)
}



calculate_t_confidence_interval <- function(x_bar, s, n, confidence_level = 0.95) {
  alpha <- 1 - confidence_level
  df <- n - 1
  
  t_critical <- qt(1 - alpha/2, df)
  
  margin_of_error <- t_critical * (s/sqrt(n))
  
  lower_bound <- x_bar - margin_of_error
  upper_bound <- x_bar + margin_of_error
  return(c(lower_bound, upper_bound))
}

calculate_t_confidence_interval_sample_sigma <- function(x_bar, sample_sigma, n, confidence_level = 0.95) {
  alpha <- 1 - confidence_level
  df <- n - 1
  
  t_critical <- qt(1 - alpha/2, df)
  
  margin_of_error <- t_critical * sample_sigma
  
  lower_bound <- x_bar - margin_of_error
  upper_bound <- x_bar + margin_of_error
  return(c(lower_bound, upper_bound))
}


shapiro_test = function(x) {
  # test whether the passed in data follows a normal distribution or not. 
  library("ggpubr")
  return(shapiro.test(x))
}

matched_pair_interval_z_score = function(mu_d, sigma_d, n, confidence_level = 0.95) {
  alpha <- 1 - confidence_level
  z_critical <- qnorm(1 - alpha/2)
  
  ME = z_critical * (sigma_d / sqrt(n))
  lower_bound = mu_d - ME
  upper_bound = mu_d + ME
  print(paste("at confidence level of", confidence_level*100, "%, interval is:"))
  return(c(lower_bound, upper_bound))
}

get_matched_pair_difference_of_means = function (x1, x2) {
  return ( mean(x2 - x1) )
}

get_matched_pair_difference_of_standard_deviations = function (x1, x2) {
  return (
    sqrt( sum( (x2 - x1 - mean(x2 - x1) )**2 )  / (length(x1)-1 ) )
  )
}

get_matched_pair_difference_of_standard_deviations_d = function (d) {
  return (
    sqrt( sum( (d - mean(d) )**2 )  / (length(d) - 1 ) )
  )
}


matched_pair_testing_t_score_with_data = function(x1, x2, mu_d, alternative = "two.sided", confidence_level = .95) {
  
  d_bar = get_matched_pair_difference_of_means(x1, x2)
  sigma_d = get_matched_pair_difference_of_standard_deviations(x1, x2)
  n = length(x1)
  alpha <- 1 - confidence_level
  
  t_ts = (d_bar - mu_d) / (sigma_d / sqrt(n))
  print(paste("test statistic t-score:", t_ts))
  
  p_value = -1
  t_crit = -100
  
  if (alternative == "two.sided") {
    t_crit = qt(1-alpha/2, n-1)
    p_value = 2*pt( -abs(t_ts), n-1)
    print("for two sided test\n")
  }
  
  else if (alternative == "upper.tail") {
    t_crit = qt(1-alpha, n-1)
    p_value = pt( -abs(t_ts), n-1)
    print("for upper tail test\n")
  }
  
  else if (alternative == "lower.tail") {
    t_crit = qt(alpha, n-1)
    p_value = pt( -abs(t_ts), n-1)
    print("for lower tail test\n")
  }
  
  else {
    stop("idek")
  }
  
  print(paste("critical t-score:", "at", confidence_level*100, "%:", t_crit))
  print(paste("p-value:", p_value))
}

matched_pair_testing_t_score = function(differences, mu_d = 0, 
                                        alternative = "two.sided", alpha=0.05) {
  
  n = length(differences)
  d_bar = mean(differences)
  sigma_d = get_matched_pair_difference_of_standard_deviations_d(differences)
  
  t_ts = (d_bar - mu_d) / (sigma_d / sqrt(n))
  print(paste("test statistic t-score:", t_ts))
  
  p_val = pt( -abs(t_ts), n-1 )
  
  p_val = testing_helper_t(alternative, p_val, alpha, n)
  return(p_val)
}

easy_matched_pair_data = function(x1, x2) {
  t.test(x1, x2, paired = TRUE, alternative = "two.sided")
}

easy_matched_pair = function(d) {
  t.test(d, mu=0)
}

matched_pair_interval_t_score = function(x_bar_d, s_d, n, confidence_level = 0.95) {
  alpha <- 1 - confidence_level
  df <- n - 1
  
  t_critical <- qt(1 - alpha/2, df)
  
  ME = t_critical * (s_d / sqrt(n))
  lower_bound = x_bar_d - ME
  upper_bound = x_bar_d + ME
  print(paste("at confidence level of", confidence_level*100, "%, interval is:"))
  return(c(lower_bound, upper_bound))
}


get_independent_pairs_unequal_variance_score = function(n1, n2, x_bar_1, x_bar_2, s_1, s_2, mu_1=0, mu_2=0) {
  SE= sqrt( s_1**2/n1 + s_2**2/n2)
  print(paste("standard error", SE))
  return ( ((x_bar_1 - x_bar_2) - (mu_1 - mu_2))/SE )
}

get_f_test_statistic = function(s1, s2, n1, n2) {
  s_x1 = s1 / sqrt(n1)
  s_x2 = s2 / sqrt(n2)
  return ( s_x1**2 / s_x2**2 )
}

# used for hypothesis testing on the proportion
proportions_test_z = function(p_hat, p, n, alternative = "two.sided", alpha = 0.05) {
  sigma_p_hat = sqrt( (p_hat*(1-p_hat))/n )
  z_test =  (  (p_hat - p) / sigma_p_hat )
  print(paste("test statistic is", z_test))
  p_val = pnorm(  -abs(z_test) )
  
  p_val = testing_helper_z(alternative, p_val, alpha)
  return(p_val)
}

# used for hypothesis testing around the mean, with T-distribution
means_test_t = function(x_bar, mu, s, n, alternative = "two.sided", alpha = 0.05) {
  s_x = s / sqrt(n)
  t_test =   (x_bar - mu) / s_x
  print(paste("test statistic is", t_test))
  
  p_val = pnorm( -abs(t_test) )
  p_val = testing_helper_t(alternative, p_val, alpha, n)
  return(p_val)
}


testing_helper_t = function(alternative, p_val, alpha, n) {
  if (alternative == "two.sided") {
    p_val = 2*p_val
    t_crit_1 = qt(alpha/2, n - 1)
    t_crit_2 = abs(qt(alpha/2, n - 1))
    print(paste("crit scores:", t_crit_1, t_crit_2))
    print("for two sided test")
  }
  
  else if (alternative == "upper.tail") {
    t_crit = abs(qt(alpha, n - 1))
    print(paste("crit score", t_crit))
    print("for upper tail test")
  }
  
  else if (alternative == "lower.tail") {
    t_crit = qt(alpha, n - 1)
    print(paste("crit score", t_crit))
    print("for lower tail test")
  }
  
  print("p=")  
  return (p_val)
}

testing_helper_z = function(alternative, p_val, alpha) {
  if (alternative == "two.sided") {
    p_val = 2*p_val
    z_crit_1 = qnorm(alpha/2)
    z_crit_2 = abs(qnorm(alpha/2))
    print(paste("crit scores:", z_crit_1, z_crit_2))
    print("for two sided test")
  }
  
  else if (alternative == "upper.tail") {
    z_crit = abs(qnorm(alpha))
    print(paste("crit score", z_crit))
    print("for upper tail test")
  }
  
  else if (alternative == "lower.tail") {
    z_crit = qnorm(alpha)
    print(paste("crit score", z_crit))
    print("for lower tail test")
  }
  
  print("p=")  
  return (p_val)
}


difference_of_proportions_test_z = function(p_hat_1, p_hat_2, n1, n2, mu_d = 0, 
                                            alternative = "two.sided", alpha = 0.05) {
  
  # 1. calculate standard error
  p_hat = (p_hat_1*n1 + p_hat_2*n2)/(n1+n2)
  sigma_p_hat = sqrt(
    (p_hat*(1-p_hat))*(1/n1 + 1/n2) 
  )  
  
  # 2. get test statistic
  z_test =((p_hat_2 - p_hat_1) - mu_d)/sigma_p_hat
  print(paste("test statistic is", z_test))
  
  # 3. get p-value
  p_val = pnorm(  -abs(z_test) )
  
  p_val = testing_helper_z(alternative, p_val, alpha)
  return(p_val)
  
  
}


difference_of_means_test_both_large = function(mu_1, mu_2, n1, n2, s1, s2, mu_d = 0,
                                               alternative = "two.sided", alpha = 0.05){
  standard_error = sqrt( (s1**2/n1) + (s2**2/n2)    )
  z_test = (mu_1 - mu_2 - (mu_d))/standard_error
  print(paste("test statistic is", z_test))
  
  # 3. get p-value
  p_val = pnorm(  -abs(z_test) )
  
  p_val = testing_helper_z(alternative, p_val, alpha)
  return(p_val)
  
}

difference_of_means_test_one_large_one_small = function(x_bar_1, x_bar_2, n1, n2, s1, s2, mu_d = 0,
                                               alternative = "two.sided", alpha = 0.05){
  standard_error = sqrt( (1/n1)+(1/n2)   )
  
  df1 = n1 - 1
  df2 = n2 -1
  pooled_variance = sqrt(
    ((df1* (s1**2))+(df2* (s2**2))) / (df1 + df2)
  )
  t_test = (x_bar_1 - x_bar_2 - (mu_d))/ (pooled_variance*standard_error)
  print(paste("test statistic is", t_test))
  
  # 3. get p-value
  p_val = pt(  -abs(t_test), df1+df2 )
  
  p_val = testing_helper_t(alternative, p_val, alpha, df1+df2+1)
  return(p_val)
  
}


difference_of_variances_test = function(s1, s2, n1, n2) {
  F_test = s1**2/s2**2
  print(paste("f test statistic", F_test))
  
  p_val = pf(F_test, n1 - 1, n2 -1) * 2 # for 2 tailed test
  print(paste("p val:",  p_val))
}

f = function(x, y) {
  return ( 1+ x**2 + 3*y)
}