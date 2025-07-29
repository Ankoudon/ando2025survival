########## library ##########

library(tidyverse)
library(latex2exp)

########## function to calculate equation (15) ##########

# p : per-contact transmissibility 
# ve : vaccine efficacy

exposure_prob_ratio <- function(p, ve) {
  
  p_ve <- p * (1 - ve)
  
  sum_num <- 0
  sum_deno <- 0
  
  sum_num <- sum_num + 1
  sum_deno <- sum_deno + 1
  
  sum_num <- sum_num + pgeom(0, 1 / 3, lower.tail = FALSE) * (1 - p_ve) ^ (1)
  sum_deno <- sum_deno + pgeom(0, 1 / 3, lower.tail = FALSE) * (1 - p) ^ (1)
  
  sum_num <- sum_num + pgeom(1, 1 / 3, lower.tail = FALSE) * (1 - p_ve) ^ (2)
  sum_deno <- sum_deno + pgeom(1, 1 / 3, lower.tail = FALSE) * (1 - p) ^ (2)
  
  for (i in 2:10) {
    
    probability <- pgeom(i, 1 / 3, lower.tail = FALSE) * (1 - p_ve) ^ (i + 1)
    sum_num <- sum_num + probability
    
    probability <- pgeom(i, 1 / 3, lower.tail = FALSE) * (1 - p) ^ (i + 1)
    sum_deno <- sum_deno + probability
    
  }
  
  return(sum_num / sum_deno)
  
}

########## function to calculate equation (16) and (17) ##########

# p : per-contact transmissibility
# ve_per : per-contact vaccine efficacy
# ve_cox : cox based vaccine efficacy

ve_eq <- function(p, ve_per, ve_cox) {
  equ <- 1 - (1 - ve_per) * exposure_prob_ratio(p, ve_per) - ve_cox
  return(equ)
}


########## function to read the simulation data ##########

read_surv_data <- function(p, ve) {
  
  data <- read_csv(
    paste0(
      "~/desktop/research/cox-correlated-exposure-bias/result/p_", p,
      "_ve_", ve, "_cox_bias.csv"))
  
  # data for ve_cox
  data_ve_cox <- data |> 
    mutate("Per-contact transmisibility" = p,
           "ve" = ve_cox,
           "model" = "ve_cox") |>
    select("Per-contact transmisibility",
           "ve",
           "model")
  
  ve_cox_vec <- data_ve_cox$ve
  
  ve_per_vec <- c()
  
  for (i in 1:length(ve_cox_vec)) {
    
    ve_per_value <- uniroot(
      ve_eq, c(0, 1), p = p, ve_cox = ve_cox_vec[i])$root
    ve_per_vec <- c(ve_per_value, ve_per_vec)
    
  }
  
  # data for ve_per
  data_ve_per <- data |> 
    mutate("Per-contact transmisibility" = p,
           "ve" = ve_per_vec,
           "model" = "ve_per") |>
    select("Per-contact transmisibility",
           "ve",
           "model")
  
  # bind the two data
  data_ve <- bind_rows(
    data_ve_cox,
    data_ve_per) 
  
  # convert the model to factor
  data_ve <- data_ve |>
    mutate(model = factor(model, levels = c("ve_cox", "ve_per")))
  
  return(data_ve)
  
}

########## Fig. 4, 5, 6 ##########

# ve = 0.3
data_0.05_0.3 <- read_surv_data(0.05, 0.3)
data_0.1_0.3 <- read_surv_data(0.1, 0.3)
data_0.15_0.3 <- read_surv_data(0.15, 0.3)

data_0.3 <- bind_rows(
  data_0.05_0.3,
  data_0.1_0.3,
  data_0.15_0.3)

# 10 by 10
fig.4 <- ggplot() +
  geom_boxplot(
    data = data_0.3,
    aes(x = factor(`Per-contact transmisibility`),
        y = ve, fill = model), color = "black") +
  scale_fill_manual(
    values = c("ve_cox" = "#8491B4FF", "ve_per" = "#4DBBD599"),
    labels = c(
      "ve_cox" = TeX("$\\widehat{v}^*$"),
      "ve_per" = TeX("$\\widehat{v}$"))) +
  labs(x = "Per-contact transmisibility (p)",
       y = "VE estimate",
       fill = NULL) +
  geom_hline(yintercept = 0.3, color = "#DC0000FF") +
  theme_minimal() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 30),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)) 

# ve = 0.6
data_0.05_0.6 <- read_surv_data(0.05, 0.6)
data_0.1_0.6 <- read_surv_data(0.1, 0.6)
data_0.15_0.6 <- read_surv_data(0.15, 0.6)

data_0.6 <- bind_rows(
  data_0.05_0.6,
  data_0.1_0.6,
  data_0.15_0.6)

# 10 by 10
fig.5 <- ggplot() +
  geom_boxplot(
    data = data_0.6,
    aes(x = factor(`Per-contact transmisibility`),
        y = ve, fill = model), color = "black",) +
  scale_fill_manual(
    values = c("ve_cox" = "#00A087FF", "ve_per" = "#91D1C299"),
    labels = c(
      "ve_cox" = TeX("$\\widehat{v}^*$"),
      "ve_per" = TeX("$\\widehat{v}$"))) +
  labs(x = "Per-contact transmisibility (p)",
       y = "VE estimate",
       fill = NULL) +
  geom_hline(yintercept = 0.6, color = "red") +
  theme_minimal() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 30),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)) 

# ve = 0.9
data_0.05_0.9 <- read_surv_data(0.05, 0.9)
data_0.1_0.9 <- read_surv_data(0.1, 0.9)
data_0.15_0.9 <- read_surv_data(0.15, 0.9)

data_0.9 <- bind_rows(
  data_0.05_0.9,
  data_0.1_0.9,
  data_0.15_0.9)

# 10 by 10
fig.6 <- ggplot() +
  geom_boxplot(
    data = data_0.9,
    aes(x = factor(`Per-contact transmisibility`),
        y = ve, fill = model), color = "black") +
  scale_fill_manual(
    values = c("ve_cox" = "#7E6148FF", "ve_per" = "#B09C8599"),
    labels = c(
      "ve_cox" = TeX("$\\widehat{v}^*$"),
      "ve_per" = TeX("$\\widehat{v}$"))) +
  labs(x = "Per-contact transmisibility (p)",
       y = "VE estimate",
       fill = NULL) +
  geom_hline(yintercept = 0.9, color = "red") +
  theme_minimal() + 
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 30),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)) 

########## Fig. 7 ##########

data_005 <- tibble(
  ve_per = seq(0.05, 0.95, by = 0.01),
  r = exposure_prob_ratio(0.05, ve_per),
  ve_cox = 1 - (1 - ve_per) * r)

data_01 <- tibble(
  ve_per = seq(0.05, 0.95, by = 0.01),
  r = exposure_prob_ratio(0.1, ve_per),
  ve_cox = 1 - (1 - ve_per) * r)

data_015 <- tibble(
  ve_per = seq(0.05, 0.95, by = 0.01),
  r = exposure_prob_ratio(0.15, ve_per),
  ve_cox = 1 - (1 - ve_per) * r)

data <- bind_rows(
  data_005 %>% mutate(p = 0.05),
  data_01 %>% mutate(p = 0.1),
  data_015 %>% mutate(p = 0.15)
)

# 12 by 8
fig.7 <- data |>
  ggplot(aes(x = ve_per, y = ve_cox / ve_per, color = factor(p))) +
  geom_line(linewidth = 2) +
  scale_color_manual(
    name = "Per-contact transmissibility (p)",
    values = c("#AD002AFF", "#00468BFF", "#FDAF91FF"),
    labels = c("0.05", "0.1", "0.15")
  ) +
  labs(
    x = TeX("${v}$\\in$ (0.05, 0.95)"),
    y = TeX("${v}^{*}$/${v}$"),
    color = NULL
  ) +
  scale_x_continuous(
    limits = c(0.05, 0.95),
    breaks = seq(0.05, 0.95, by = 0.1)  
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.7, 0.2),
    legend.title = element_text(size = 25),
    legend.text = element_text(size = 20),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)
  )















