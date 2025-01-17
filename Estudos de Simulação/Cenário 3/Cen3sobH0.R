library(future.apply)
library(DirichletReg)

# Funcao a ser executada para cada tamanho amostral
cen3h0 <- function(n, rep=5000, B=1000){
  # Definindo a semente principal
  set.seed(8261)
  
  # Definindo o ambiente de paralelizacao
  future::plan(multisession)
  
  # Definindo os betas arbitrarios
  betas <- matrix(data=c(-0.5,1.1,0,
                         -0.5,1.1,0,
                         2.3,0,0), 
                  nrow = 3,ncol = 3,byrow = TRUE)
  # Definindo a matriz X com as covariaveis
  X <- cbind(1, runif(n), runif(n))
  colnames(X) <- c('x11', 'x12', 'x13')
  
  # Calculando etas do modelo
  etas <- X %*% t(betas)
  colnames(etas) <- c('eta1', 'eta2', 'eta3')
  # Calculando as medias de cada componente
  denom <- 1 + exp(etas[, "eta1"]) + exp(etas[, "eta2"])
  medias <- cbind(1 / denom, 
                  exp(etas[, "eta1"]) / denom, 
                  exp(etas[, "eta2"]) / denom)
  # Calculando o parametro de precisao
  phi <- exp(etas[,"eta3"])
  # Calculando os parametros alpha
  alphas <- medias * phi
  
  #########################################################
  ## Loop externo paralelizado (replicas de Monte Carlo) ##
  #########################################################
  resultados <- future_lapply(1:rep, function(i) {
    # Guia para ver se o codigo esta rodando ou travou
    if (i %% 500 == 0) {
      print(paste("Simulacao", i))
    }
    
    # Gerando a variavel resposta
    y <- rdirichlet(n, alphas)
    # Criando um dataframe com a resposta e as covariaveis
    dados <- data.frame(y, x12 = X[, "x12"], x13 = X[, "x13"])
    # Tratando a resposta para o ajuste
    dados$y <- DR_data(dados[, 1:3])
    
    # Ajustando modelos completo e reduzido
    mod_completo <- DirichReg(y ~ x12 + x13, data = dados, model = 'alternative')
    mod_reduzido <- DirichReg(y ~ x12, data = dados, model = 'alternative')
    # Calculando valor-p do TRV
    veros_reduzido <- mod_reduzido$logLik
    veros_completo <- mod_completo$logLik
    estat_teste <- 2 * (veros_completo - veros_reduzido)
    valorp_trv <- 1 - pchisq(estat_teste, df = 2)
    
    # Obtendo betas do modelo reduzido
    b <- as.vector(mod_reduzido$coefficients)
    b <- c(b[1:2], 0, b[3:4], 0, b[5], 0, 0)
    # Definindo os betas para bootstrap
    betas_boot <- matrix(data=b, nrow = 3, ncol = 3, byrow = TRUE)
    
    # Calculando etas para bootstrap
    etas_boot <- X %*% t(betas_boot)
    colnames(etas_boot) <- c('eta1', 'eta2', 'eta3')
    # Calculando as medias para bootstrap
    denom_boot <- 1 + exp(etas_boot[, "eta1"]) + exp(etas_boot[, "eta2"])
    medias_boot <- cbind(1 / denom_boot, 
                         exp(etas_boot[, "eta1"]) / denom_boot, 
                         exp(etas_boot[, "eta2"]) / denom_boot)
    # Calculando o parametro de precisao para bootstrap
    phi_boot <- exp(etas_boot[,"eta3"])
    # Calculando os parametros alpha para bootstrap
    alphas_boot <- medias_boot * phi_boot
    
    ##########################################
    ## Loop interno (replicas de Bootstrap) ##
    ##########################################
    estatisticas_boot <- sapply(1:B, function(j) {
      # Gerando a variavel resposta
      y_boot <- rdirichlet(n, alphas_boot)
      # Criando um dataframe com a resposta e as covariaveis
      dados_boot <- data.frame(y_boot, x12 = X[, "x12"], x13 = X[, "x13"])
      # Tratando a resposta para o ajuste
      dados_boot$y_boot <- DR_data(dados_boot[, 1:3])
      
      # Ajustando modelos completo e reduzido
      mod_completo_boot <- DirichReg(y_boot ~ x12 + x13, data = dados_boot, model = 'alternative')
      mod_reduzido_boot <- DirichReg(y_boot ~ x12, data = dados_boot, model = 'alternative')
      # Calculando a estatistica de teste do TRV e armazenando
      veros_reduzido_boot <- mod_reduzido_boot$logLik
      veros_completo_boot <- mod_completo_boot$logLik
      
      2 * (veros_completo_boot - veros_reduzido_boot)
    })
    # Calculando valor-p do Boot1
    estat_corrigida <- (estat_teste * 2) / mean(estatisticas_boot)
    valorp_boot1 <- 1 - pchisq(estat_corrigida, df = 2)
    
    # Calculando valor-p do Boot2
    k <- sum(estatisticas_boot > estat_teste)
    valorp_boot2 <- (k + 1) / (B + 1)
    
    # Retornando uma lista para cada iteracao de Monte Carlo
    list(valorp_trv = valorp_trv, 
         valorp_boot1 = valorp_boot1, 
         valorp_boot2 = valorp_boot2)
  },
  # Definindo parametro que controla as sementes dentro dos nucleos paralelos
  future.seed = TRUE)
  
  # Juntando todos os resultados das iteracoes paralelizadas
  valoresp_trv   <- sapply(resultados, `[[`, "valorp_trv")
  valoresp_boot1 <- sapply(resultados, `[[`, "valorp_boot1")
  valoresp_boot2 <- sapply(resultados, `[[`, "valorp_boot2")
  
  # Calculando taxas de rejeicao
  niveis_rejeicao_trv <- c(mean(valoresp_trv   <= 0.01) * 100, 
                           mean(valoresp_trv   <= 0.05) * 100, 
                           mean(valoresp_trv   <= 0.1) * 100)
  niveis_rejeicao_boot1 <- c(mean(valoresp_boot1 <= 0.01) * 100, 
                             mean(valoresp_boot1 <= 0.05) * 100, 
                             mean(valoresp_boot1 <= 0.1) * 100)
  niveis_rejeicao_boot2 <- c(mean(valoresp_boot2 <= 0.01) * 100, 
                             mean(valoresp_boot2 <= 0.05) * 100, 
                             mean(valoresp_boot2 <= 0.1) * 100)
  
  # Organizando resultados finais em uma tabela
  resultados_finais <- data.frame(rbind(niveis_rejeicao_trv, 
                                        niveis_rejeicao_boot1, 
                                        niveis_rejeicao_boot2))
  colnames(resultados_finais) <- c("1%", "5%", "10%")
  rownames(resultados_finais) <- c("TRV", "Boot1", "Boot2")
  return(resultados_finais)
}

# Executando diferentes tamanhos amostrais
cen3h0(20,5000,1000)
cen3h0(30,5000,1000)
cen3h0(40,5000,1000)
cen3h0(50,5000,1000)
cen3h0(100,5000,1000)