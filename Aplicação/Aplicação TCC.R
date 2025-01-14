library(tidyverse)
library(DirichletReg)
library(gamlss)

# Carregando o banco de dados
dados <- read_table("C:/Users/vini-/OneDrive/Área de Trabalho/TCC/dados_aplicacao.txt")

# Removendo a coluna de identificação
dados <- dados %>% 
  select(-ID)

########################
## Análise Descritiva ##
########################
# Transformando os dados para o formato long
dados_long <- dados %>%
  pivot_longer(cols = c("N1", "N2", "N3", "REM"), names_to = "Componente", values_to = "Valor")
# Gráfico de dispersão para variável TST
ggplot(dados_long, aes(x = TST, y = Valor, color = Componente)) +
  geom_point(alpha = 0.6) +
  facet_wrap(~ Componente) +
  theme_minimal() +
  labs(x = "Tempo de Sono Total", y = "Percentual do Componente",
       title = "Percentual de Qualidade do Sono por TST") +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73"))

# Categorizando TST em 4 classes
dados_long <- dados_long %>%
  mutate(
    TST_cat = cut(
      TST, breaks = 4,
      labels = c("Baixo", "Médio-Baixo", "Médio-Alto", "Alto"),
      include.lowest = TRUE))
# Calculando as médias para cada categoria
media_categorias <- dados_long %>%
  group_by(TST_cat, Componente) %>%
  summarise(Media_Valor = mean(Valor, na.rm = TRUE)) %>%
  ungroup()
# Gráfico de perfil para variável TST categorizada
ggplot(media_categorias, aes(x = TST_cat, y = Media_Valor, color = Componente, group = Componente)) +
  geom_point(size = 3) +
  geom_line(alpha = 0.7) +
  theme_minimal() +
  labs(
    x = "Categorias de TST",
    y = "Percentual Médio do Componente",
    color = "Componente",
    title = "Percentual Médio de Qualidade do Sono por Categorias de TST"
  ) +
  scale_color_manual(values = c("#000000", "#E69F00", "#56B4E9", "#009E73"))

# Boxplots para variável Caso-Controle por componente da resposta
ggplot(dados_long, 
       aes(x = interaction(CaseControl, Componente), 
           y = Valor, fill = factor(CaseControl))) +
  geom_boxplot() +
  theme_minimal() +
  labs(x = "Caso-Controle por Componente",
       y = "Percentual da Combinação",
       fill = "Caso-Controle",
       title = "Boxplots para as Combinações de Caso-Controle e Componente") +
  scale_x_discrete(labels = function(x) gsub("\\.", " - ", x)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
  scale_fill_discrete(labels = c("0" = "Controle = 0", "1" = "Caso = 1"))


#########################
## Testes de Hipóteses ##
#########################
# Preparando os dados para ajustar o modelo
y = dados[, 3:6]
x12 = dados$TST
x13 = dados$CaseControl
dados$y = DR_data(dados[, 3:6])
colnames(dados$y) = c(c("y1","y2","y3","y4"))

### Testando covariável Caso-Controle
# Semente para reprodutibilidade
set.seed(8261)

# Aplicando o TRV para testar todos os betas (media e precisao)
mod_completo <- DirichReg(y ~ x12 + x13 | x12 + x13, data = dados, model = 'alternative')
mod_reduzido <- DirichReg(y ~ x12 | x12, data = dados, model = 'alternative')

veros_reduzido <- mod_reduzido$logLik
veros_completo <- mod_completo$logLik
estat_teste <- 2 * (veros_completo - veros_reduzido)
valorp_trv1 <- 1 - pchisq(estat_teste, df = 4)

# Aplicando o TRV para testar os betas relacionados as medias
mod_completo <- DirichReg(y ~ x12 + x13 | x12 + x13, data = dados, model = 'alternative')
mod_reduzido <- DirichReg(y ~ x12 | x12 + x13, data = dados, model = 'alternative')

veros_reduzido <- mod_reduzido$logLik
veros_completo <- mod_completo$logLik
estat_teste <- 2 * (veros_completo - veros_reduzido)
valorp_trv2 <- 1 - pchisq(estat_teste, df = 3)

# Aplicando o TRV para testar o beta do parametro de precisao
mod_completo <- DirichReg(y ~ x12 + x13 | x12 + x13, data = dados, model = 'alternative')
mod_reduzido <- DirichReg(y ~ x12 + x13 | x12, data = dados, model = 'alternative')

veros_reduzido <- mod_reduzido$logLik
veros_completo <- mod_completo$logLik
estat_teste <- 2 * (veros_completo - veros_reduzido)
valorp_trv3 <- 1 - pchisq(estat_teste, df = 1)

# Tabela com os resultados
data.frame(TRV_x13 = c(valorp_trv1, valorp_trv2, valorp_trv3))


### Testando covariável TST
## Teste para todos os betas (media e precisao)
# Semente para reprodutibilidade
set.seed(8261)

# Definindo tamanho amostral e quantidade de réplicas de bootstrap
n <- nrow(dados)
B = 1000

# Aplicando o TRV para testar todos os betas (media e precisao)
mod_completo <- DirichReg(y ~ x12 | x12, data = dados, model = 'alternative')
mod_reduzido <- DirichReg(y ~ 1, data = dados, model = 'alternative')

veros_reduzido <- mod_reduzido$logLik
veros_completo <- mod_completo$logLik
estat_teste <- 2 * (veros_completo - veros_reduzido)
valorp_trv <- 1 - pchisq(estat_teste, df = 4)

# Obtendo betas do modelo reduzido
b <- as.vector(mod_reduzido$coefficients)
b <- c(b[1], 0, b[2], 0, b[3], 0, b[4], 0)
# Definindo os betas para os testes via bootstrap
betas_boot <- matrix(data=b, 
                     nrow = 4,ncol = 2,byrow = TRUE)
# Definindo a matriz X com as covariaveis
X <- cbind(1, dados$TST)
colnames(X) <- c('x11', 'x12')

# Calculando etas para os testes via bootstrap
etas_boot <- X %*% t(betas_boot)
colnames(etas_boot) <- c('eta1', 'eta2', 'eta3', 'eta4')
# Calculando as medias para os testes via bootstrap
denom_boot <- 1 + exp(etas_boot[, "eta1"]) + exp(etas_boot[, "eta2"]) + exp(etas_boot[, "eta3"])
medias_boot <- cbind(1 / denom_boot, 
                     exp(etas_boot[, "eta1"]) / denom_boot, 
                     exp(etas_boot[, "eta2"]) / denom_boot, 
                     exp(etas_boot[, "eta3"]) / denom_boot)
# Calculando o parametro de precisao para os testes via bootstrap
phi_boot <- exp(etas_boot[,"eta4"])
# Calculando os parametros alpha para os testes via bootstrap
alphas_boot <- medias_boot * phi_boot

# Criando vetor para armazenar as réplicas de bootstrap
estatisticas_boot <- numeric(B)
for (j in 1:B){
  if (j %% 100 == 0) {
    print(paste("Simulação", j))
  }
  # Gerando a variavel resposta
  y_boot <- rdirichlet(n, alphas_boot)
  # Criando um dataframe com a resposta e as covariaveis
  dados_boot <- data.frame(y_boot, x12 = X[, "x12"])
  # Tratando a resposta para o ajuste
  dados_boot$y_boot <- DR_data(dados_boot[, 1:4])
  
  # Ajustando modelos completo e reduzido
  mod_completo_boot <- DirichReg(y_boot ~ x12 | x12, data = dados_boot, model = 'alternative')
  mod_reduzido_boot <- DirichReg(y_boot ~ 1, data = dados_boot, model = 'alternative')
  # Calculando a estatistica de teste do TRV e armazenando
  veros_reduzido_boot <- mod_reduzido_boot$logLik
  veros_completo_boot <- mod_completo_boot$logLik
  estatisticas_boot[j] <- 2 * (veros_completo_boot - veros_reduzido_boot)
}
# Calculando valor-p do Boot1
estat_corrigida <- (estat_teste*4) / mean(estatisticas_boot)
valorp_boot1 <- 1 - pchisq(estat_corrigida, df = 4)

# Calculando valor-p do Boot2
k <- sum(estatisticas_boot > estat_teste)
valorp_boot2 <- (k+1) / (B+1)

# Organizando resultados em uma tabela
todos <- data.frame(Valores_p = c(valorp_trv, valorp_boot1, valorp_boot2))


## Teste para os betas relacionados as medias
# Semente para reprodutibilidade
set.seed(8261)

# Aplicando o TRV para testar os betas relacionados as medias
mod_completo <- DirichReg(y ~ x12 | x12, data = dados, model = 'alternative')
mod_reduzido <- DirichReg(y ~ 1 | x12, data = dados, model = 'alternative')

veros_reduzido <- mod_reduzido$logLik
veros_completo <- mod_completo$logLik
estat_teste <- 2 * (veros_completo - veros_reduzido)
valorp_trv <- 1 - pchisq(estat_teste, df = 3)

# Obtendo betas do modelo reduzido
b <- as.vector(mod_reduzido$coefficients)
b <- c(b[1], 0, b[2], 0, b[3], 0, b[4:5])
# Definindo os betas para os testes via bootstrap
betas_boot <- matrix(data=b, 
                     nrow = 4,ncol = 2,byrow = TRUE)

# Calculando etas para os testes via bootstrap
etas_boot <- X %*% t(betas_boot)
colnames(etas_boot) <- c('eta1', 'eta2', 'eta3', 'eta4')
# Calculando as medias para os testes via bootstrap
denom_boot <- 1 + exp(etas_boot[, "eta1"]) + exp(etas_boot[, "eta2"]) + exp(etas_boot[, "eta3"])
medias_boot <- cbind(1 / denom_boot, 
                     exp(etas_boot[, "eta1"]) / denom_boot, 
                     exp(etas_boot[, "eta2"]) / denom_boot, 
                     exp(etas_boot[, "eta3"]) / denom_boot)
# Calculando o parametro de precisao para os testes via bootstrap
phi_boot <- exp(etas_boot[,"eta4"])
# Calculando os parametros alpha para os testes via bootstrap
alphas_boot <- medias_boot * phi_boot

# Criando vetor para armazenar as réplicas de bootstrap
estatisticas_boot <- numeric(B)
for (j in 1:B){
  if (j %% 100 == 0) {
    print(paste("Simulação", j))
  }
  # Gerando a variavel resposta
  y_boot <- rdirichlet(n, alphas_boot)
  # Criando um dataframe com a resposta e as covariaveis
  dados_boot <- data.frame(y_boot, x12 = X[, "x12"])
  # Tratando a resposta para o ajuste
  dados_boot$y_boot <- DR_data(dados_boot[, 1:4])
  
  # Ajustando modelos completo e reduzido
  mod_completo_boot <- DirichReg(y_boot ~ x12 | x12, data = dados_boot, model = 'alternative')
  mod_reduzido_boot <- DirichReg(y_boot ~ 1 | x12, data = dados_boot, model = 'alternative')
  # Calculando a estatistica de teste do TRV e armazenando
  veros_reduzido_boot <- mod_reduzido_boot$logLik
  veros_completo_boot <- mod_completo_boot$logLik
  estatisticas_boot[j] <- 2 * (veros_completo_boot - veros_reduzido_boot)
}
# Calculando valor-p do Boot1
estat_corrigida <- (estat_teste*3) / mean(estatisticas_boot)
valorp_boot1 <- 1 - pchisq(estat_corrigida, df = 3)

# Calculando valor-p do Boot2
k <- sum(estatisticas_boot > estat_teste)
valorp_boot2 <- (k+1) / (B+1)

# Organizando resultados em uma tabela
mus <- data.frame(Valores_p = c(valorp_trv, valorp_boot1, valorp_boot2))


## Teste para os betas relacionados as medias
# Semente para reprodutibilidade
set.seed(8261)

# Aplicando o TRV para testar o beta do parametro de precisao
mod_completo <- DirichReg(y ~ x12 | x12, data = dados, model = 'alternative')
mod_reduzido <- DirichReg(y ~ x12 | 1, data = dados, model = 'alternative')

veros_reduzido <- mod_reduzido$logLik
veros_completo <- mod_completo$logLik
estat_teste <- 2 * (veros_completo - veros_reduzido)
valorp_trv <- 1 - pchisq(estat_teste, df = 1)

# Obtendo betas do modelo reduzido
b <- as.vector(mod_reduzido$coefficients)
b <- c(b, 0)
# Definindo os betas para os testes via bootstrap
betas_boot <- matrix(data=b, 
                     nrow = 4,ncol = 2,byrow = TRUE)

# Calculando etas para os testes via bootstrap
etas_boot <- X %*% t(betas_boot)
colnames(etas_boot) <- c('eta1', 'eta2', 'eta3', 'eta4')
# Calculando as medias para os testes via bootstrap
denom_boot <- 1 + exp(etas_boot[, "eta1"]) + exp(etas_boot[, "eta2"]) + exp(etas_boot[, "eta3"])
medias_boot <- cbind(1 / denom_boot, 
                     exp(etas_boot[, "eta1"]) / denom_boot, 
                     exp(etas_boot[, "eta2"]) / denom_boot, 
                     exp(etas_boot[, "eta3"]) / denom_boot)
# Calculando o parametro de precisao para os testes via bootstrap
phi_boot <- exp(etas_boot[,"eta4"])
# Calculando os parametros alpha para os testes via bootstrap
alphas_boot <- medias_boot * phi_boot

# Criando vetor para armazenar as réplicas de bootstrap
estatisticas_boot <- numeric(B)
for (j in 1:B){
  if (j %% 100 == 0) {
    print(paste("Simulação", j))
  }
  # Gerando a variavel resposta
  y_boot <- rdirichlet(n, alphas_boot)
  # Criando um dataframe com a resposta e as covariaveis
  dados_boot <- data.frame(y_boot, x12 = X[, "x12"])
  # Tratando a resposta para o ajuste
  dados_boot$y_boot <- DR_data(dados_boot[, 1:4])
  
  # Ajustando modelos completo e reduzido
  mod_completo_boot <- DirichReg(y_boot ~ x12 | x12, data = dados_boot, model = 'alternative')
  mod_reduzido_boot <- DirichReg(y_boot ~ x12 | 1, data = dados_boot, model = 'alternative')
  # Calculando a estatistica de teste do TRV e armazenando
  veros_reduzido_boot <- mod_reduzido_boot$logLik
  veros_completo_boot <- mod_completo_boot$logLik
  estatisticas_boot[j] <- 2 * (veros_completo_boot - veros_reduzido_boot)
}
# Calculando valor-p do Boot1
estat_corrigida <- (estat_teste*1) / mean(estatisticas_boot)
valorp_boot1 <- 1 - pchisq(estat_corrigida, df = 1)

# Calculando valor-p do Boot2
k <- sum(estatisticas_boot > estat_teste)
valorp_boot2 <- (k+1) / (B+1)

# Organizando resultados em uma tabela
phis <- data.frame(Valores_p = c(valorp_trv, valorp_boot1, valorp_boot2))


### Análise adicional
# Semente para reprodutibilidade
set.seed(8261)

# Criando vetores para armazenar os valores-p da analise
valoresp_trv <- numeric(500)
valoresp_boot1 <- numeric(500)
valoresp_boot2 <- numeric(500)

# Analise: taxas de rejeicao se amostra tivesse tamanho 15
for(i in 1:500){
  # Tirando uma amostra sem reposicao de tamanho 15 do banco de dados original
  amostra <- dados[sample(1:nrow(dados), size = 15, replace = FALSE), ]
  
  # Preparando os dados da amostra para ajustar o modelo
  y = amostra[, 3:6]
  x12 = amostra$TST
  amostra$y = DR_data(amostra[, 3:6])
  colnames(amostra$y) = c(c("y1","y2","y3","y4"))
  
  # Definindo tamanho amostral e quantidade de réplicas de bootstrap
  n <- nrow(amostra)
  B = 1000
  
  # Aplicando o TRV para testar os betas relacionados as medias
  mod_completo <- DirichReg(y ~ x12 | x12, data = amostra, model = 'alternative')
  mod_reduzido <- DirichReg(y ~ 1 | x12, data = amostra, model = 'alternative')
  
  veros_reduzido <- mod_reduzido$logLik
  veros_completo <- mod_completo$logLik
  estat_teste <- 2 * (veros_completo - veros_reduzido)
  valoresp_trv <- 1 - pchisq(estat_teste, df = 3)
  
  # Obtendo betas do modelo reduzido
  b <- as.vector(mod_reduzido$coefficients)
  b <- c(b[1], 0, b[2], 0, b[3], 0, b[4:5])
  # Definindo os betas para os testes via bootstrap
  betas_boot <- matrix(data=b, 
                       nrow = 4,ncol = 2,byrow = TRUE)
  # Definindo a matriz X com as covariaveis da amostra
  X <- cbind(1, amostra$TST)
  colnames(X) <- c('x11', 'x12')
  
  # Calculando etas para os testes via bootstrap
  etas_boot <- X %*% t(betas_boot)
  colnames(etas_boot) <- c('eta1', 'eta2', 'eta3', 'eta4')
  # Calculando as medias para os testes via bootstrap
  denom_boot <- 1 + exp(etas_boot[, "eta1"]) + exp(etas_boot[, "eta2"]) + exp(etas_boot[, "eta3"])
  medias_boot <- cbind(1 / denom_boot, 
                       exp(etas_boot[, "eta1"]) / denom_boot, 
                       exp(etas_boot[, "eta2"]) / denom_boot, 
                       exp(etas_boot[, "eta3"]) / denom_boot)
  # Calculando o parametro de precisao para os testes via bootstrap
  phi_boot <- exp(etas_boot[,"eta4"])
  # Calculando os parametros alpha para os testes via bootstrap
  alphas_boot <- medias_boot * phi_boot
  
  # Criando vetor para armazenar as réplicas de bootstrap
  estatisticas_boot <- numeric(B)
  for (j in 1:B){
    if (j %% 100 == 0) {
      print(paste("Simulação", j))
    }
    # Gerando a variavel resposta
    y_boot <- rdirichlet(n, alphas_boot)
    # Criando um dataframe com a resposta e as covariaveis
    dados_boot <- data.frame(y_boot, x12 = X[, "x12"])
    # Tratando a resposta para o ajuste
    dados_boot$y_boot <- DR_data(dados_boot[, 1:4])
    
    # Ajustando modelos completo e reduzido
    mod_completo_boot <- DirichReg(y_boot ~ x12 | x12, data = dados_boot, model = 'alternative')
    mod_reduzido_boot <- DirichReg(y_boot ~ 1 | x12, data = dados_boot, model = 'alternative')
    # Calculando a estatistica de teste do TRV e armazenando
    veros_reduzido_boot <- mod_reduzido_boot$logLik
    veros_completo_boot <- mod_completo_boot$logLik
    estatisticas_boot[j] <- 2 * (veros_completo_boot - veros_reduzido_boot)
  }
  # Calculando valor-p do Boot1
  estat_corrigida <- (estat_teste*3) / mean(estatisticas_boot)
  valoresp_boot1 <- 1 - pchisq(estat_corrigida, df = 3)
  
  # Calculando valor-p do Boot2
  k <- sum(estatisticas_boot > estat_teste)
  valoresp_boot2 <- (k+1) / (B+1)
}
# Calculando taxas de rejeicao
niveis_rejeicao_trv <- c(mean(valoresp_trv <= 0.01) * 100, mean(valoresp_trv <= 0.05) * 100, 
                         mean(valoresp_trv <= 0.1) * 100)
niveis_rejeicao_boot1 <- c(mean(valoresp_boot1 <= 0.01) * 100, mean(valoresp_boot1 <= 0.05) * 100, 
                           mean(valoresp_boot1 <= 0.1) * 100)
niveis_rejeicao_boot2 <- c(mean(valoresp_boot2 <= 0.01) * 100, mean(valoresp_boot2 <= 0.05) * 100, 
                           mean(valoresp_boot2 <= 0.1) * 100)

# Organizando resultados em uma tabela
resultados <- data.frame(rbind(niveis_rejeicao_trv, niveis_rejeicao_boot1, niveis_rejeicao_boot2))
colnames(resultados) <- c("1%", "5%", "10%")

#################
## Diagnóstico ##
#################
# Semente para reprodutibilidade
set.seed(8261)

# Definindo tamanho amostral e quantidade de réplicas de bootstrap
n <- nrow(dados)
B = 1000

# Preparando os dados para ajustar o modelo
y = dados[, 2:5]
x12 = dados$TST
dados$y = DR_data(dados[, 2:5])
colnames(dados$y) = c(c("y1","y2","y3","y4"))
dados$x12 = dados$TST

# Ajustando o modelo final
ajuste <- DirichReg(y ~ x12 | x12, data = dados, model = 'alternative')

# Criando uma matriz para guardar os resíduos de cada componente
residuos = matrix(0,n,4)
# Calculando resíduos quantílicos aleatorizados
residuos[,1] = qnorm(pBE(unlist(y[,1]), ajuste$fitted.values$mu[,1], 1/(1+ajuste$fitted.values$phi)^0.5), 0, 1)
residuos[,2] = qnorm(pBE(unlist(y[,2]), ajuste$fitted.values$mu[,2], 1/(1+ajuste$fitted.values$phi)^0.5), 0, 1)
residuos[,3] = qnorm(pBE(unlist(y[,3]), ajuste$fitted.values$mu[,3], 1/(1+ajuste$fitted.values$phi)^0.5), 0, 1)
residuos[,4] = qnorm(pBE(unlist(y[,4]), ajuste$fitted.values$mu[,4], 1/(1+ajuste$fitted.values$phi)^0.5), 0, 1)

# Tirando o valor absoluto dos resíduos
abs_residuos <- abs(residuos)
# Identificando qual coluna tem o resíduo de maior valor absoluto e retornando o índice dessa coluna
max_indices <- max.col(abs_residuos, ties.method = "first")
# Extraindo o resíduo (com sinal) correspondente ao índice máximo
max_residuos <- residuos[cbind(1:n, max_indices)]
# Calculando a funcao sinal para os resíduos dominantes
hi <- sign(max_residuos)
# Calculando li: somando os valores absolutos de todos os resíduos de cada obs., mas com sinal definido pelo resíduo dominante
li <- hi * rowSums(abs_residuos)

# Obtendo os parametros alpha com base no modelo final ajustado
alphas_boot <- ajuste$fitted.values$alpha
# Criando uma matriz para armazenar os li para cada réplica
li_boot <- matrix(NA, nrow = n, ncol = B)
for(i in 1:B){
  # Gerando nova amostra com parametros do modelo ajustado
  y_boot <- rdirichlet(n, alphas_boot)
  # Criando um dataframe com a resposta e as covariaveis
  dados_boot <- as.data.frame(y_boot)
  colnames(dados_boot) <- c("y1", "y2", "y3", "y4")
  dados_boot$x12 <- x12
  # Tratando a resposta para o ajuste
  dados_boot$y_boot <- DR_data(dados_boot[, 1:4])
  
  # Ajustando um modelo com a nova amostra bootstrap
  ajuste_boot <- DirichReg(y_boot ~ x12 | x12, data = dados_boot, model = 'alternative')
  
  # Calculando e armazenando os resíduos para a nova amostra
  residuos_boot = matrix(0,n,4)
  residuos_boot[,1] = qnorm(pBE(unlist(y_boot[,1]), ajuste_boot$fitted.values$mu[,1], 1/(1+ajuste_boot$fitted.values$phi)^0.5), 0, 1)
  residuos_boot[,2] = qnorm(pBE(unlist(y_boot[,2]), ajuste_boot$fitted.values$mu[,2], 1/(1+ajuste_boot$fitted.values$phi)^0.5), 0, 1)
  residuos_boot[,3] = qnorm(pBE(unlist(y_boot[,3]), ajuste_boot$fitted.values$mu[,3], 1/(1+ajuste_boot$fitted.values$phi)^0.5), 0, 1)
  residuos_boot[,4] = qnorm(pBE(unlist(y_boot[,4]), ajuste_boot$fitted.values$mu[,4], 1/(1+ajuste_boot$fitted.values$phi)^0.5), 0, 1)
  
  # Mesmo procedimento descrito anteriormente até calcular li
  abs_residuos_boot <- abs(residuos_boot)
  max_indices_boot <- max.col(abs_residuos_boot, ties.method = "first")
  max_residuos_boot <- residuos_boot[cbind(1:n, max_indices_boot)]
  hi_boot <- sign(max_residuos_boot)
  li_boot[,i] <- hi_boot * rowSums(abs_residuos_boot)
}

# Criando vetores ai e si a serem preenchidos
ai <- rep(-1, n)
si <- rep(-1, n)
for(j in 1:n){
  # Calculando ai: quantas réplicas li_boot é menor que li_original
  ai[j] <- sum(li_boot[j,] < li[j])
  # Resíduo final proposto por Pereira e Cai (2024)
  si[j] = qnorm(runif(1, ai[j]/(B+1), (ai[j]+1)/(B+1)), 0, 1) 
}


### Repetição do processo acima 19 vezes para criar envelope simulado
n_sim <- 19

# Matriz para armazenar os resíduos simulados
residuos_simulados <- matrix(NA, nrow = n, ncol = n_sim)
# Obtendo os parametros alpha com base no modelo final ajustado
alphas_env <- ajuste$fitted.values$alpha
for(e in 1:n_sim){
  # Verificacao de qual simulacao esta
  print(e)
  
  y_env <- rdirichlet(n, alphas_env)
  dados_env <- as.data.frame(y_env)
  colnames(dados_env) <- c("y1", "y2", "y3", "y4")
  dados_env$x12 <- x12
  dados_env$y_env <- DR_data(dados_env[, 1:4])
  
  ajuste_env <- DirichReg(y_env ~ x12 | x12, data = dados_env, model = 'alternative')
  
  residuos_env = matrix(0,n,4)
  residuos_env[,1] = qnorm(pBE(unlist(y_boot[,1]), ajuste_env$fitted.values$mu[,1], 1/(1+ajuste_env$fitted.values$phi)^0.5), 0, 1)
  residuos_env[,2] = qnorm(pBE(unlist(y_boot[,2]), ajuste_env$fitted.values$mu[,2], 1/(1+ajuste_env$fitted.values$phi)^0.5), 0, 1)
  residuos_env[,3] = qnorm(pBE(unlist(y_boot[,3]), ajuste_env$fitted.values$mu[,3], 1/(1+ajuste_env$fitted.values$phi)^0.5), 0, 1)
  residuos_env[,4] = qnorm(pBE(unlist(y_boot[,4]), ajuste_env$fitted.values$mu[,4], 1/(1+ajuste_env$fitted.values$phi)^0.5), 0, 1)
  
  abs_residuos_env <- abs(residuos_env)
  max_indices_env <- max.col(abs_residuos_env, ties.method = "first")
  max_residuos_env <- residuos_env[cbind(1:n, max_indices_env)]
  hi_env <- sign(max_residuos_env)
  li_env <- hi_env * rowSums(abs_residuos_env)
  
  alphas_boot_env <- ajuste_env$fitted.values$alpha
  li_boot_env <- matrix(NA, nrow = n, ncol = B)
  for(ii in 1:B){
    y_boot_env <- rdirichlet(n, alphas_boot_env)
    dados_boot_env <- as.data.frame(y_boot_env)
    colnames(dados_boot_env) <- c("y1", "y2", "y3", "y4")
    dados_boot_env$x12 <- x12
    dados_boot_env$y_boot_env <- DR_data(dados_boot_env[, 1:4])
    
    ajuste_boot_env <- DirichReg(y_boot_env ~ x12 | x12, data = dados_boot_env, model = 'alternative')
    
    residuos_boot_env = matrix(0,n,4)
    residuos_boot_env[,1] = qnorm(pBE(unlist(y_boot_env[,1]), ajuste_boot_env$fitted.values$mu[,1], 1/(1+ajuste_boot_env$fitted.values$phi)^0.5), 0, 1)
    residuos_boot_env[,2] = qnorm(pBE(unlist(y_boot_env[,2]), ajuste_boot_env$fitted.values$mu[,2], 1/(1+ajuste_boot_env$fitted.values$phi)^0.5), 0, 1)
    residuos_boot_env[,3] = qnorm(pBE(unlist(y_boot_env[,3]), ajuste_boot_env$fitted.values$mu[,3], 1/(1+ajuste_boot_env$fitted.values$phi)^0.5), 0, 1)
    residuos_boot_env[,4] = qnorm(pBE(unlist(y_boot_env[,4]), ajuste_boot_env$fitted.values$mu[,4], 1/(1+ajuste_boot_env$fitted.values$phi)^0.5), 0, 1)
    
    abs_residuos_boot_env <- abs(residuos_boot_env)
    max_indices_boot_env <- max.col(abs_residuos_boot_env, ties.method = "first")
    max_residuos_boot_env <- residuos_boot_env[cbind(1:n, max_indices_boot_env)]
    hi_boot_env <- sign(max_residuos_boot_env)
    li_boot_env[,ii] <- hi_boot_env * rowSums(abs_residuos_boot_env)
  }
  
  ai_env <- rep(-1, n)
  si_env <- rep(-1, n)
  for(jj in 1:n){
    ai_env[jj] <- sum(li_boot[jj,] < li_env[jj])
    si_env[jj] = qnorm(runif(1, ai_env[jj]/(B+1), (ai_env[jj]+1)/(B+1)), 0, 1) 
  }
  residuos_simulados[,e] <- si_env
}
# Combinando os residuos originais com os simulados
residuos_completos <- cbind(si, residuos_simulados)
# Ordenando cada coluna individualmente
residuos_ordenados <- apply(residuos_completos, 2, sort)

# Calculando o mínimo de cada linha desconsiderando a primeira coluna
minimos <- apply(residuos_ordenados[, -1], 1, min)
# Calculando o máximo de cada linha desconsiderando a primeira coluna
maximos <- apply(residuos_ordenados[, -1], 1, max)

# Montando eixo x com quantis teóricos
x_env = rep(0, n)
for (a in 1:n){
  x_env[a] = qnorm((a-3/8)/(n+1/4),0,1)
}

# Definindo as bandas de confianca combinando min e max de todas as colunas
faixa <- range(residuos_ordenados[,1],minimos,maximos)
# Gráfico de envelope simulado
plot(x_env,residuos_ordenados[,1],xlab="Expected normal order statistic",
     ylab="", ylim=faixa, pch=16)
par(new=T)
# Linha inferior do envelope
lines(x_env,minimos)
# Linha superior do envelope
lines(x_env,maximos)




















