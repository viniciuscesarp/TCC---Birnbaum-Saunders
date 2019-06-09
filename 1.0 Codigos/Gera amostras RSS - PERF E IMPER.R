
# Códigos para obtenção de amostras via RSS sob ORDENAÇÃO PERFEITA - BS

AMOSTRA_RSS_PERF <- function(m, n, alpha, beta){
  amostra <- rbs(n = m*(n^2), alpha = alpha, beta = beta)
  amostraT<- matrix(amostra, nrow = n)
  amostra_ordenada <- apply(amostraT, 2, sort)
  pos <- cbind(rep(1:n, times = m), 1:(n*m))
  amostra_2 <- amostra_ordenada[pos]
  amostra_final <- matrix(amostra_2, nrow = m, byrow = TRUE)
  return(amostra_final)
}
 
AMOSTRA_RSS_PERF(3, 6, 1, 3) 

# ============================================================================ #
# Código para obtenção de amostras via RSS sob ORDENAÇÃO IMPERFEITA - BS

AMOSTRA_RSS_IMPERF <- function(m, n, sigma_xy, alpha, beta){
  library(MASS)
  library(VGAM)
  medias <- c(0, 0)
  sigmas <- matrix(c(1, sigma_xy, sigma_xy, 1), nrow = 2)
  amostra <- mvrnorm(n = m*(n^2), mu = medias, Sigma = sigmas)

  u <- pnorm(amostra[,2], 0, 1)
  bs <- qbisa(u, scale = beta, shape = alpha)
  
  matrix_x <- matrix(amostra[, 1], nrow = n)
  matrix_y <- matrix(bs, nrow = n)

  y_atuali <- matrix(0, ncol = n*m, nrow = n)
  x_ordenado <- apply(matrix_x, 2, order)

  for(i in 1:(m*n)){
    y_atuali[,i] <- matrix_y[x_ordenado[,i],i]
  }

  indices <- rep(1:n, times = m)
  posicao <- cbind(indices, 1:(m*n))
  y_final <- y_atuali[posicao]
  return(y_final)
}
AMOSTRA_RSS_IMPERF(6, 3, 1, 1, 1)
amem
plot(density(amem))


AM_RSS_IMPERF <- function(m, n, sigma_xy, alpha, beta){
  library(MASS)
  library(VGAM)
  medias <- c(0, 0)
  sigmas <- matrix(c(1, sigma_xy, sigma_xy, 1), nrow = 2)
  amostra <- mvrnorm(n = m*(n^2), mu = medias, Sigma = sigmas)
  
  u <- pnorm(amostra[,2], 0, 1)
  bs <- qbs(u, alpha  = alpha, beta = beta)
  
  matrix_x <- matrix(amostra[, 1], nrow = n)
  matrix_y <- matrix(bs, nrow = n)
  
  y_atuali <- matrix(0, ncol = n*m, nrow = n)
  x_ordenado <- apply(matrix_x, 2, order)
  
  for(i in 1:(m*n)){
    y_atuali[,i] <- matrix_y[x_ordenado[,i],i]
  }
  
  indices <- rep(1:n, times = m)
  posicao <- cbind(indices, 1:(m*n))
  y_final <- y_atuali[posicao]
  return(y_final)
}

