# install.packages("tmvtnorm",dependencies=TRUE, repos='http://cran.rstudio.com/')

library(tmvtnorm)
 
source('case_2/case_2_preamble.R')
source("case_2/abcSMC.R")

msl_data <- read.csv(file.path("data", "Measles_data_time.csv"))

# Length of observed outbreak (weeks)
n_obs<-nrow(msl_data)

# Age classes: 0.5-1, 1-2,..., 19-20, 20-100
age<-data.frame(lower=c(0.5, seq(1,20)), upper=c(seq(1,20),100)) 
n_age<-nrow(age)  

#  Observed percentage in three age groups
msl_age_perc<-c(42, 30, 28) 
cl.to.test<-list(1:5, 6:14, 15:21)

#  Latent period
nu<- 7 

# Infectious period
mu<- 7 

# Max duration of an outbreak
t_end<-1000 

# Number of particles
N <- 1000

# Epsilon values for temporal data 
epsilon_T <- c(50000, 27500, 25000, 22500, 20000, 17500, 15000, 12500, 11000, 10000)  

# Epsilon values for age data 
epsilon_A <- c(15, 10, 9, 8, 7, 6, 5, 4, 3.5, 3)

#  Lower and upper boundaries for priors
lm.low <- c(140*10^3, 0.8, 0.1, 0, 0)
lm.upp <- c(300*10^3, 1.4, 0.4, 5*10^(-6), 0.001)


abcSMC(data = msl_data, N = N, epsilon_A = epsilon_A, epsilon_T = epsilon_T, fileName = "results_case_2_ABC_SMC_gen_")
