## t-value
t_val <- function(x, mu, sigma) {
  t_value = (x - mu) / sigma
  return(t_value)
}

## Population variance
p_variance <- function(c) {
  n = length(c)
  e_of_squared = (1/n) * sum(c^2)
  e_squared = ((1/n) * sum(c))^2 
  population_variance = e_of_squared - e_squared
  return(population_variance)
}

## Sample variance
s_variance <- function(c) {
  n = length(c)
  mean = (1/n) * sum(c)
  sum_of_squared_differences = 0
  for (i in 1:n) {
    sum_of_squared_differences = sum_of_squared_differences + (c[i] - mean)^2 
  }
  sample_variance = (1/(n - 1)) * sum_of_squared_differences
  return(sample_variance)
}

## Population covariance
p_covariance <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  n <- (nx + ny) / 2
  e_of_x <- (1/nx) * sum(x) 
  e_of_y <- (1/ny) * sum(y)
  differenced_x <- x - e_of_x
  differenced_y <- y - e_of_y
  
  population_covariance <- (1/n) * sum(differenced_x * differenced_y) 
  return(population_covariance)
}

## Sample covariance
s_covariance <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  n <- (nx + ny) / 2
  e_of_x <- (1/nx) * sum(x) 
  e_of_y <- (1/ny) * sum(y)
  differenced_x <- x - e_of_x
  differenced_y <- y - e_of_y
  
  sample_covariance = (1/(n - 1)) * sum(differenced_x * differenced_y)
  return(sample_covariance)
}

## Pearson's chi-squared
chi_squared <- function(table) {
  s = 0
  for (i in 1:nrow(table)) {
    for (j in 1:ncol(table)) {
      e = (sum(table[i,]) * sum(table[,j])) / sum(table[,])
      summand = ((table[i, j] - e)^2) / e
      s = s + summand
    }
  }  
  return(s)
}

## Cramer's V
cramers_v <- function(table) {
  chi <- chi_squared(table)
  n <- sum(table[,])
  k <- min(nrow(table), ncol(table))
  chi_max <- n * (k - 1)
  v <- sqrt(chi / chi_max)
  
  return(v)
}
