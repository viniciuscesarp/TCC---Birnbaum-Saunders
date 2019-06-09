
# ============================================================================ #
# Código para geração de amostras BS via AAS sob Ordenação Perfeita

########################### BS via AAS ###########################

##' @title Função para obtenção de uma amostra de tamanho n via AAS
##' @description Essa função permite obter uma amostra de tamanho n 
##' via AAS para a distribuição BS.
##' @param n: Tamanho da amostra.
##' @param alpha: Parâmetro da distribuição BS.
##' @param beta: Parâmetro da distribuição BS.
##' @return Retorna uma amostra de tamanho n via AAS.
##' 
AMOSTRA_AAS <- function(n, alpha=1.0, beta=1.0){
  if(!is.numeric(n)||!is.numeric(alpha)||!is.numeric(beta)) 
  {stop("non-numeric argument to mathematical function")}
  if(n<=0){stop("value of n must be greater or equal then 1")}
  if(alpha<=0){stop("alpha must be positive")}
  if(beta<=0){stop("beta must be positive")}
  method<-c("ig")
  if(method=="ig"){
    rig<-function(n,mu=1.0,lambda=1.0){
      z<-rnorm(n,0,1)
      v0<-z^2
      x1<-mu+(((mu^2)*v0)/(2*lambda))-(mu/(2*lambda))*sqrt(4*mu*lambda*v0+((mu^2)*(v0^2)))
      x2<-(mu^2)/x1
      p0<-mu/(mu*(mu+x1))
      u0<-runif(n,0,1)
      y<-seq(1,n,by=1)
      for(i in 1:n)
      {
        val1<-u0[i]
        val2<-p0[i]
        if(val1<=val2){y[i]<-x1[i]}else{y[i]<-x2[i]}
      }
      return(y)
    }
    x1<-rig(n,beta,(alpha^(-2))*beta)
    x2<-rig(n,beta^(-1),(alpha^(-2))*(beta^(-1)))
    s<-1/x2
    w<-rbinom(n,1,0.5)
    t<-(w*x1)+((1-w)*s)
  }
  return(t)
}

AMOSTRA_AAS(n = 10000, alpha = 2, beta = 3)

# ============================================================================ #
# Código para geração de amostras BS via RSS 

########################### BS via RSS ###########################

##' @title Função para obtenção de uma amostra de tamanho m*(n^2) via RSS
##' @description Essa função permite obter uma amostra de tamanho m*(n^2) 
##' via RSS pela ordenação perfeita e imperfeita.
##' @param m: Número de ciclos.
##' @param n: Tamanho da amostra.
##' @param cov_xy: covariância entre X e Y.
##' @param alpha: Parâmetro da distribuição BS.
##' @param beta: Parâmetro da distribuição BS.
##' @return Retorna uma amostra de tamanho m*(n^2) via RSS por ordenação
##' perfeita e imperfeita.

AMOSTRA_RSS <- function(m, n, cov_xy, alpha, beta){
  library(MASS)
  medias <- c(0, 0)
  sigmas <- matrix(c(1, cov_xy, cov_xy, 1), nrow = 2)
  amostra <- mvrnorm(n = m*(n^2), mu = medias, Sigma = sigmas)
  
  u <- pnorm(amostra[,2], 0, 1)
  bs <- qbs(u, alpha = alpha, beta = beta)
  
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

AMOSTRA_RSS(10, 3, 1, 2, 2)





