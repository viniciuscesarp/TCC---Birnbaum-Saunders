
# ============================================================================ #
# # Código para obtenção de estimativas via ESTIMADOR DE MÁXIMA VEROSSIMILHANÇA
# POR RSS

########################### EMV via RSS ###########################

##' @title Função para obtenção de estimativas para BS via EMV
##' @description Essa função permite obter estimativas para os parâmetros da 
##' distribuição BS.
##' @param amostra: Conjunto de dados.
##' @param alpha: Parâmetro da distribuição BS.
##' @param beta: Parâmetro da distribuição BS.
##' @return Retorna estimativas para os parâmetros da BS.
##' 
# EMV_RSS <- function(dad, alpha, beta, m, n){
#   indices <- rep(1:n, times = m)
#   lalpha <- log(alpha)
#   lbeta <- log(beta)
#   u0 <- function(alpha, beta){
#     -sum(factorial(n)/(factorial(indices-1) * factorial(n-1)) * 
#            log(dbs(dad, exp(alpha), exp(beta))*((pbs(dad, exp(alpha), exp(beta)))^(indices - 1)) *
#                  ((1 - pbs(dad, alpha = exp(alpha), beta = exp(beta)))^(n - 1))))
#   }
#   senhor <- mle2(u0, start = list(alpha = lalpha, beta = lbeta), method = "L-BFGS-B")
#   return(coef(senhor))
# }
# dados <- AMOSTRA_RSS(6, 3, 1, 5, 2)
# EMV_RSS(dados, 5, 2, 6, 3)

EMV_RSS <- function(dados, alpha, beta, m, n){
  library(VGAM)
  indices <- rep(1:n, times = m)  
  x <- data.frame(indices = indices, dados = dados)
  
  est_ord <- function(x, indice, n, alpha, beta){
    (indice - 1) * pbisa(x, beta, alpha, log.p = TRUE) + (n - indice) * 
      pbisa(x, beta, alpha, log.p = TRUE, lower.tail = FALSE) + 
      dbisa(x, beta, alpha, log = TRUE)
  }
  lik_RSS <- function(x, n, alpha, beta){
    p1 <- exp(alpha)
    p2 <- exp(beta)
    soma <- vector(mode = "numeric", length = length(unique(x[,1])))
    for (j in unique(x[,1])){
      soma[j] <- sum(est_ord(x = x[x[,1] == j,2], j , n, alpha = p1, beta = p2))}
    logs <- - sum(soma)
    return(logs)
  }
  m_lik_RSS <- mle2(lik_RSS, start=list(alpha = alpha, beta = beta),
                    data = list(x = x, n = n), optimizer="nlminb")
  final <- exp(coef(m_lik_RSS))
  sd_rss <- sqrt(diag(vcov(m_lik_RSS)))
  return(c(final, sd_rss))
}






# ============================================================================ #
# Código para obtenção de estimativas via ESTIMADOR DE MOMENTOS POR RSS

########################### EMM via RSS ###########################

##' @title Função para obtenção de estimativas para BS via EMM
##' @description Essa função permite obter estimativas para os parâmetros da 
##' distribuição BS.
##' @param x: Conjunto de dados.
##' @return Retorna estimativas para os parâmetros da BS.
##' 
agora_vai <- function(x){
  n <-  length(x)
  s <- mean(x)
  v <- var(x)
  alph <- ((-2*((s^2)-4)+2*sqrt((((s^2)-4)^2)+v*(5*(s^2)-v)))/(5*(s^2) - v))^0.5
  bet <- (2*s)/ ((alph ^2) + 2)
  cv <- (sd(x)/s) 
  return(c(alph,bet, cv))
}
dados <- AMOSTRA_RSS(10, 3, 1, 2, 2)
agora_vai(dados)


TESTE_EMM <- function(x){
  n <-  length(x)
  s <- mean(x)
  v <- var(x) * (n - 1)/n
  u <- function(alpha){
    (alpha^4)*(5*(s^2) - v) + 4*(alpha^2)*((s^2)-4) - 4*v
  }
  # senhor <- optim(par = list(alpha = 1), u, method =  "L-BFGS-B")
  senhor <- uniroot(u, c(0.01, 100))
  bet <- (2*s)/ ((senhor$root ^2) + 2)
  return(c(senhor$root, bet))
}
dados <- AMOSTRA_RSS(10, 3, 1, 2, 2)
TESTE_EMM(dados)

# ============================================================================ #
# # Código para obtenção de estimativas via ESTIMADOR DE MOMENTOS MODIFICADOS
# POR RSS

########################### EMMM via RSS ###########################

##' @title Função para obtenção de estimativas para BS via EMMM
##' @description Essa função permite obter estimativas para os parâmetros da 
##' distribuição BS.
##' @param X: Conjunto de dados.
##' @return Retorna estimativas para os parâmetros da BS.
##' 
EMMM_RSS <- function(x){
  s <- mean(x)
  r <- harmonic.mean(x)
  alph <- (2 * ((s / r)^0.5 - 1))^0.5
  bet <- (s*r)^(1/2)
  return(c(alph, bet))
}
dados <- AMOSTRA_RSS(6, 3, 1, 2, 2)
EMMM_RSS(dados)
