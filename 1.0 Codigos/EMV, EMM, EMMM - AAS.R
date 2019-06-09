
# ============================================================================ #
# # Código para obtenção de estimativas via ESTIMADOR DE MÁXIMA VEROSSIMILHANÇA
# POR AAS

########################### EMV via AAS ###########################

##' @title Função para obtenção de estimativas para BS via EMV
##' @description Essa função permite obter estimativas para os parâmetros da 
##' distribuição BS.
##' @param amostra: Conjunto de dados.
##' @param alpha: Parâmetro da distribuição BS.
##' @param beta: Parâmetro da distribuição BS.
##' @return Retorna estimativas para os parâmetros da BS.
##' 
EMV_AAS <- function(amostra, alpha, beta){
  library(bbmle)
  library(VGAM)
  u <- function(alpha, beta){- sum(dbisa(amostra, beta, alpha, log = T))}
  senhor <- mle2(u, start = list(alpha = alpha, beta = beta), 
                 method = 'L-BFGS-B')
  sd_aas <- sqrt(diag(vcov(senhor)))
  return(c(coef(senhor), sd_aas))
}

dados <- AMOSTRA_AAS(180, 3, 1)
EMV_AAS(dados, alpha = 3, beta = 1)

# ============================================================================ #
# Código para obtenção de estimativas via ESTIMADOR DE MOMENTOS POR AAS

########################### EMM via AAS ###########################

##' @title Função para obtenção de estimativas para BS via EMM
##' @description Essa função permite obter estimativas para os parâmetros da 
##' distribuição BS.
##' @param amostra: Conjunto de dados.
##' @param alpha: Parâmetro da distribuição BS.
##' @param beta: Parâmetro da distribuição BS.
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
dados <- AMOSTRA_AAS(1000, alpha = 2, beta = 3)
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
dados <- AMOSTRA_AAS(1000, alpha = 2, beta = 3)
TESTE_EMM(dados)




 
# ============================================================================ #
# # Código para obtenção de estimativas via ESTIMADOR DE MOMENTOS MODIFICADOS
# POR AAS

########################### EMMM via AAS ###########################

##' @title Função para obtenção de estimativas para BS via EMMM
##' @description Essa função permite obter estimativas para os parâmetros da 
##' distribuição BS.
##' @param X: Conjunto de dados.
##' @return Retorna estimativas para os parâmetros da BS.
##' 
EMMM_AAS <- function(x, n){
  library(psych)
  s <- mean(x)
  r <- harmonic.mean(x)
  alpha <- (2 * ((s / r)^0.5 - 1))^0.5
  beta <- (s*r)^(1/2)
  sd_alpha <- sqrt((alpha^2) / (2*n))
  
  sd_beta <- sqrt((((alpha*beta)^2)/n) * ((1 + ((3/4) *(alpha^2))) / ((1+ (0.5*(alpha^2)))^2)))
  
  return(c(alpha, beta, sd_alpha, sd_beta))
}

dados <- AMOSTRA_AAS(180, 1, 2)
EMMM_AAS(dados,18)
