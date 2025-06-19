########## library ##########
library(tidyverse)
library(latex2exp)

########## function to calculate the formula (4) without m ##########
r <- function(p, ve) {
  
  # p: per-contact transmissibility after vaccination
  p_ve <- p * (1 - ve)
  
  sum <- 0
  
  sum <- sum + 1
  sum <- sum + pgeom(0, 1 / 3, lower.tail = FALSE) * (1 - p_ve) ^ (1)
  sum <- sum + pgeom(1, 1 / 3, lower.tail = FALSE) * (1 - p_ve) ^ (2)
  
  for (i in 2:50) {
    
    probability <- pgeom(i, 1 / 3, lower.tail = FALSE) * (1 - p_ve) ^ (i + 1)
    sum <- sum + probability
    
  }
  
  return(sum)
  
}

########## function to read the simlation data ##########

read_surv_data <- function(p, ve) {
  
  data <- read_csv(paste0("~/desktop/p_", p, "_ve_", ve, ".csv"))
  
  # data for ve_hr
  data_ve_hr <- data |> 
    mutate("Per-contact transmisibility" = p,
           "ve" = 1 - VE_HR,
           "model" = "ve_hr") |>
    select("Per-contact transmisibility",
           "ve",
           "model")
  
  # data for ve_hr_adj
  data_ve_adj <- data |> 
    mutate("Per-contact transmisibility" = p,
           "ve" = 1 - VE_HR * (r(p, 0)/r(p, ve)),
           "model" = "ve_adj") |>
    select("Per-contact transmisibility",
           "ve",
           "model")
  
  # bind the two data
  data_ve <- bind_rows(
    data_ve_hr,
    data_ve_adj) 
  
  # convert the model to factor
  data_ve <- data_ve |>
    mutate(model = factor(model, levels = c("ve", "ve_adj")))
  
  return(data_ve)
  
}

########## Fig. 5, 6, 7 ##########

# ve = 0.3
data_0.05_0.3 <- read_surv_data(0.05, 0.3)
data_0.1_0.3 <- read_surv_data(0.1, 0.3)
data_0.15_0.3 <- read_surv_data(0.15, 0.3)

data_0.3 <- bind_rows(
  data_0.05_0.3,
  data_0.1_0.3,
  data_0.15_0.3)

# 10 by 10
fig.5 <- ggplot() +
  geom_boxplot(
    data = data_0.3,
    aes(x = factor(`Per-contact transmisibility`),
        y = ve, fill = model), color = "black") +
  scale_fill_manual(
    values = c("ve_hr" = "#8491B4FF", "ve_adj" = "#4DBBD599"),
    labels = c(
      "ve_hr" = TeX("$\\hat{VE}_{S, HR}$"),
      "ve_adj" = TeX(
        "$1-(\\hat{r}_{pla} / \\hat{r}_{vac}) \\cdot (1-\\hat{VE}_{S, HR})$"))) +
  labs(x = "Per-contact transmisibility (p)",
       y = "VE estimate",
       fill = NULL) +
  geom_hline(yintercept = 0.3, color = "#DC0000FF") +
  theme_minimal() + 
  theme(
    legend.position = c(0.8, 0.9),
    legend.text = element_text(size = 20),
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
fig.6 <- ggplot() +
  geom_boxplot(
    data = data_0.6,
    aes(x = factor(`Per-contact transmisibility`),
        y = ve, fill = model), color = "black",) +
  scale_fill_manual(
    values = c("ve_hr" = "#00A087FF", "ve_adj" = "#91D1C299"),
    labels = c(
      "ve_hr" = TeX("$\\hat{VE}_{S, HR}$"),
      "ve_adj" = TeX(
        "$1-(\\hat{r}_{pla} / \\hat{r}_{vac}) \\cdot (1-\\hat{VE}_{S, HR})$"))) +
  labs(x = "Per-contact transmisibility (p)",
       y = "VE estimate",
       fill = NULL) +
  geom_hline(yintercept = 0.6, color = "red") +
  theme_minimal() + 
  theme(
    legend.position = c(0.8, 0.9),
    legend.text = element_text(size = 20),
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
    values = c("ve_hr" = "#7E6148FF", "ve_adj" = "#B09C8599"),
    labels = c(
      "ve_hr" = TeX("$\\hat{VE}_{S, HR}$"),
      "ve_adj" = TeX(
        "$1-(\\hat{r}_{pla} / \\hat{r}_{vac}) \\cdot (1-\\hat{VE}_{S, HR})$"))) +
  labs(x = "Per-contact transmisibility (p)",
       y = "VE estimate",
       fill = NULL) +
  geom_hline(yintercept = 0.9, color = "red") +
  theme_minimal() + 
  theme(
    legend.position = c(0.8, 0.9),
    legend.text = element_text(size = 20),
    axis.line = element_line(colour = "black"),
    panel.background = element_blank(),
    axis.title.x = element_text(size = 25),
    axis.title.y = element_text(size = 25),
    axis.text.x = element_text(size = 20),
    axis.text.y = element_text(size = 20)) 

########## Fig. 4 ##########

data_005 <- tibble(
  ve = seq(0.01, 0.99, by = 0.01),
  r = r(0.05, ve)/r(0.05, 0),
  ve_hr = 1 - (1 - ve) * r)

data_01 <- tibble(
  ve = seq(0.01, 0.99, by = 0.01),
  r = r(0.1, ve)/r(0.1, 0),
  ve_hr = 1 - (1 - ve) * r)

data_015 <- tibble(
  ve = seq(0.01, 0.99, by = 0.01),
  r = r(0.15, ve)/r(0.15, 0),
  ve_hr = 1 - (1 - ve) * r)

data <- bind_rows(
  data_005 %>% mutate(p = 0.05),
  data_01 %>% mutate(p = 0.1),
  data_015 %>% mutate(p = 0.15)
)

fig.4 <- data |>
  ggplot(aes(x = ve, y = ve_hr / ve, color = factor(p))) +
  geom_line(linewidth = 2) +
  scale_color_manual(
    name = "Per-contact transmissibility (p)",
    values = c("#AD002AFF", "#00468BFF", "#FDAF91FF"),
    labels = c("0.05", "0.1", "0.15")
  ) +
  labs(
    x = TeX("${VE}$\\in$ (0, 1)"),
    y = TeX("${VE}^{*}_{S,HR}$/${VE}$"),
    color = NULL
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
