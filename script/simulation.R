########## library ##########

library(tidyverse)

########## parameters ##########

# 0.01 (fixed)
m <- 0.01
# 10000 (fixed)
n <- 10000
# {0.05, 0.1, 0.15}
p <- 0.15
# {0.3, 0.6, 0.9}
ve <- 0.9
# 180 (fixed)
trial_length <- 180
# 500 (fixed)
session <- 500

########## tibble for the result ##########

result <- tibble(
  ve_cox = as.numeric(),
)

########## simulation ##########

# 500 simulations

for (i in 1:session) {
  
  print(i)
  
  set.seed(i)
  
  # participants data
  participants_data <- tibble(
    state = 0,
    num_exposures = 0,
    vaccine = sample(c(rep(0, n/2), rep(1, n/2)), n, replace = FALSE),
    day = 1,
    participant_id = 1:n,
    inf_window = 0,
    inf_window_periods = vector("list", n)
  )
  
  # 180 days o
  for (t in 1:trial_length) {
    
    # update the survival day 
    participants_data <- participants_data |>
      mutate(
        day = if_else(state == 0, t, day))
    
    # update the number of exposures
    participants_data <- participants_data |> mutate(
      inf_window_periods = map2(inf_window_periods,
                                participant_id, ~{.x[.x >= t]}),
      inf_window = map_int(inf_window_periods, length)  
    )
    
    # extract the IDs of participants who are susceptible
    suseptible_id <- participants_data |>
      filter(state == 0) |>
      pull(participant_id)
    
    # randomly select participants who get into the infection window
    inf_window_flag <- rbinom(length(suseptible_id), size=1, prob=m)
    # extract the IDs of participants  who get into the infection window
    inf_window_ids <- suseptible_id[inf_window_flag == 1]
    
    # update the infection window for the selected participants
    participants_data$inf_window_periods[inf_window_ids] <- 
      map2(
        participants_data$inf_window_periods[inf_window_ids],
        inf_window_ids,
        ~ append(.x, t + rgeom(1, 1/3))
      )
    
    participants_data <- participants_data |> 
      mutate(
        infection_prob = case_when(
          vaccine == 0 ~ 1 - (1 - p) ^ inf_window,
          vaccine == 1 ~ 1 - (1 - p * (1 - ve)) ^ inf_window,
          TRUE ~ 0),
        newly_infected = rbinom(n(), 1, infection_prob),
        state = if_else(state == 0 & newly_infected == 1, 1, state)
      )
    
  }
  
  # calculate the hazard ratio using Cox hazard models
  cox_model <- survival::coxph(
    survival::Surv(day, state) ~ vaccine,
    data = participants_data)
  
  # extract the hazard ratio for the vaccine group
  hr <- exp(coef(cox_model))[[1]]
  
  # add the hazard ratio to the result tibble
  result <- result |> 
    add_row(ve_cox = 1 - hr)
  
}

# save the result
write_csv(result, paste0("~/desktop/research/cox-correlated-exposure-bias/result/p_",
                         p, "_ve_", ve, "_cox_bias.csv"))
