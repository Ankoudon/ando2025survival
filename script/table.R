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

# per-contact VE estimator for table 1.

uniroot(ve_eq, c(0, 1), p = 0.05, ve_cox = 0.575)$root 
uniroot(ve_eq, c(0, 1), p = 0.1, ve_cox = 0.575)$root 
uniroot(ve_eq, c(0, 1), p = 0.15, ve_cox = 0.575)$root 

# CIs for per-contact VE for table 1.

u_0.05 <- uniroot(ve_eq, c(0, 1), p = 0.05, ve_cox = 0.748)$root
l_0.05 <- uniroot(ve_eq, c(0, 1), p = 0.05, ve_cox = 0.282)$root
print(c(l_0.05, u_0.05))

u_0.1 <- uniroot(ve_eq, c(0, 1), p = 0.1, ve_cox = 0.748)$root
l_0.1 <- uniroot(ve_eq, c(0, 1), p = 0.1, ve_cox = 0.282)$root
print(c(l_0.1, u_0.1))

u_0.15 <- uniroot(ve_eq, c(0, 1), p = 0.15, ve_cox = 0.748)$root
l_0.15 <- uniroot(ve_eq, c(0, 1), p = 0.15, ve_cox = 0.282)$root
print(c(l_0.15, u_0.15))

