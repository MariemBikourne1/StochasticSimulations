library(readxl)
library(pracma) 
library(tseries)
library(dplyr)
library(stats)
library(zoo)
library(nloptr) 
library(glmnet)
library(plyr)
library(readr)
library(dplyr)
library(Rcpp)
library(caret)
library(ggplot2)
library(repr)

file <- read_excel("Data.xlsx", sheet = "Sheet1")

# Convert 'date_column' to Date class if it's not already
file$Quarter <- as.Date(file$Quarter)

# Filter out rows where the year is greater than or equal to 1990
file <- file[file$Quarter > as.Date("1990-01-01"), ]

# Define data variables
GDP <- ts(file$"GDP",start=c(1990,1),frequency=4)
Consumption <-  ts(file$"Consumption",start=c(1990,1),frequency=4)
Quarter <-  ts(file$"Quarter",start=c(1990,1),frequency=4)

Time <- 136
dt <- 0.25 
time <- seq(1, Time, by = dt)

# Define the likelihood empty vectors
f <- numeric(length(GDP))
L <- numeric(length(GDP))
M <- numeric(length(GDP))
expectation <- numeric(length(GDP))
stdev <- numeric(length(GDP))
expectation_J <- numeric(length(GDP))
stdev_J <- numeric(length(GDP))
ll <- numeric(length(GDP))

# Define NLL functions
G_Y <- function(gamma)
{
  return(lambda *( ((1+gamma)^(2-alpha) - 1)/(2-alpha) -1 - gamma/2))
}

G_K <- function(gamma)
{
  return(lambda *( ((1+gamma)^(2-alpha) - 1)/(2-alpha) -1 - gamma*(1-alpha)/2))
}

# Define the log-likelihood function for a normal distribution
log_likelihood <- function(parameters){
  data <- GDP
  t1 <- parameters[1]
  t2 <- parameters[2]
  t3 <- parameters[3]
  t4 <- parameters[4]
  mu_J <- parameters[5]
  sigma_J <- parameters[6]
  lambda <- parameters[7]
  
  ll[1] <- 0
  
  for (i in seq(1, Time, by = dt)) {
    f[i] <- t1 * data[i] + t2 * data[i] ^ t3 
    L[i] <- t1 + t2*t3*(data[i]^(t3-1))
    M[i] <- ((t4*(2-t3))^2/2)*t2*t3*(t3-1)*data[i]^(t3-2)
    
    expectation[i] <- data[i] + f[i]* exp(L[i]*dt-1)/(L[i]) + M[i]* (exp(L[i]*dt)-1-L[i]*dt)/(L[i]^2) 
    expectation_J[i] <- data[i] + f[i]* exp(L[i]*dt-1)/(L[i]) + M[i]* (exp(L[i]*dt)-1-L[i]*dt)/(L[i]^2)  + mu_J
    
    
    stdev[i] <- (t4*(2-t3)) * sqrt(exp(2*L[i]*dt-1)/(2*L[i]))
    stdev_J[i] <- sqrt((t4*(2-t3))^2 * (exp(2*L[i]*dt-1)/(2*L[i])) + sigma_J^2)
    
    # Calculate the log-likelihood ll
    ll[i+dt] <- ll[i] + (1-lambda*dt)*dnorm(data[i+dt], mean = expectation[i], sd = stdev[i], log = TRUE)
    + lambda*dt * dnorm(data[i+dt], mean = expectation_J[i], sd = stdev_J[i], log = TRUE)
  }
  print(sum(ll))
  return(-sum(ll))
}

# Constraint for positivity for the mean of GDP equation
constraint <- function(x)
{
  data <- GDP
  return( c(x[1] + x[2]*x[3]*(data^(x[3]-1)))
  )
}  

# Define an initial guess for the parameters
initial_parameters <- c(6,0.09,-0.2,0.006, 0.1, 0.01, 0.2)

lower_bounds <-  c(-Inf, -Inf,-Inf, -Inf,-Inf, -Inf)
upper_bounds <- c(Inf, Inf, Inf, Inf, Inf, Inf)

opt <- list(
  algorithm = "NLOPT_LN_COBYLA",  # Choose the optimization algorithm (COBYLA in this case)
  x0 = initial_parameters,        # Initial parameter values
  lb = lower_bounds,              # Lower bounds
  ub = upper_bounds,              # Upper bounds
  eval_f = log_likelihood,        # Objective function
  xtol_rel = 1e-5, 
  maxeval = 20000
)

# Perform the optimisation method
result <- nloptr(x0 = initial_parameters, eval_f = log_likelihood, opts = opt, eval_g_ineq = constraint)

# Extract the estimated parameters
print(result$message)
print(result$solution)

# store to the parameter vector Theta = (t1, t2, t3, t4, mu_J, sigma_J, lambda)
t1 <- result$solution[1]
t2 <- result$solution[2]
t3 <- result$solution[3]
t4 <- result$solution[4]
mu_J <- result$solution[5]
sigma_J <- result$solution[6]
lambda <- result$solution[7]

# calculate GDP SDE parameters
alpha <- 1/(2-t3) 
A <- (t2/alpha)^alpha
sigma <- t4*(2-t3)
gamma <- 1 #choosed

#find the consumption coefficients using OLS method 
dGDP <- diff(GDP)
dConsumption <- diff(Consumption)
dta <- data.frame(dConsumption,dGDP)
lm_model <- lm(dConsumption ~ dGDP-1 , data = dta)
summary(lm_model)

#find the consumption coefficients using MLE method 

# Define the negative log-likelihood function
neg_log_likelihood <- function(a, x, y) {
  n <- length(y)
  sigma_squared <- var(y - a * x)
  -(-n/2 * log(2 * pi * sigma_squared) - 1/(2 * sigma_squared) * sum((y - a * x)^2))
}

# Initial guess for the slope parameter 
initial_guess <- 0.1
# Perform MLE using optim function
mle_result <- optim(initial_guess, neg_log_likelihood, x = dGDP, y = dConsumption, method = "BFGS")
print(mle_result)

# Perform execution time (microseconds) for OLS and MLE
microbenchmark(
  OLS = {
    lm_model <- lm(dConsumption ~ dGDP-1 , data = dta)
    
  },
  MLE = {
    mle_result <- optim(initial_guess, neg_log_likelihood, x = dGDP, y = dConsumption, method = "BFGS")
  },
  times = 1000
)

# Store to the parameter of equation c(t)=Y(t)*omega2
omega2 <- lm_model$coefficients[1]/A

# Calculate the remaining variables
delta <- 0.5*(alpha-1)*sigma^2 -omega2 - t1/alpha - G_K(gamma)
rho <- alpha*(alpha-1)*sigma^2 - t1 - delta - G_Y(gamma)


# Graphical representation
Y <- numeric(Time)
Z <- numeric(Time)
Y[1] <- GDP[1]
Z[1] <- 1

N <- 1000
trajectories <- matrix(0, ncol = length(1:Time), nrow = N)

# Simulate the jump sizes H_i
H <- rnorm(Time, mean = mu_J, sd = sigma_J)
plot(1:Time, H, type = "h")
mean(H)

# Perform Rejection Sampling method
J <- numeric(Time)
for (t in 1:Time){
  U1 <- runif(1, min = 0, max = 1)
  if(U1 < lambda*dt)
  {
    J[t] <-1
  }
  else
  {
    J[t] <- 0
    
  }
}
plot(1:Time, J , type="s")

# Simulate Monte Carlo trajectories
for (i in 1:N) {
  dW <- numeric(Time)
  WN <- randn(Time,1) #white noise
  dW <- sqrt(dt) * WN  
  
  for (t in 1:(length(1:Time)-1) ){
    Y[t+1] <- Y[t] +( (- delta - rho - alpha*(1-alpha)*(sigma^2) - G_Y(gamma) )* Y[t] + alpha * A^(1/alpha)*Y[t]^(2-1/alpha) )* dt + alpha* sigma * Y[t]* dW[t+1] + Y[t]*J[t+1]*((1+gamma*H[t+1])^alpha-1)
    #Y[t+1] <- Y[t] + (- (delta + rho) - alpha*(1-alpha)*(sigma^2) +1/2 )*dt+alpha * A^(1/alpha)*(Y[t])^(1-1/alpha)* dt + alpha* sigma * dW[t]+alpha*log(X_t[t])
  }
  trajectories[i, ] <- Y[]
}

g <- t(trajectories)
date <- seq.Date(as.Date("1990-01-01"), by = "quarter", length.out = length(GDP)) 

time <- as.yearmon(date)
g <-  zoo(g, order.by = time)
h <- rowMeans(g)
label <- " Quarter"
paste(label)

matplot(time, g/1000, type = "l", col = "lightgray", main="With Jumps", lty = 1 , xaxt = "n", alpha=0.5, xlab = "Quarter", ylab = "GDP ")
lines(time,h/1000, lty=2)
lines(time, GDP/1000, lwd=2)
legend("topleft", legend = c("GDP", "Random paths", "Mean path"), col = c("black", "lightgray", "black"), lty = c(1,1,2), lwd=c(2,2,1), cex=0.7)
axis(1, at = seq(1990, 2025, by = 5))


c(alpha, delta , rho, sigma, A, lambda, mu_J, sigma_J)

# Calculate the value function and its derivative
v <- function(x, alpha)
{
  return(A/(rho*omega2^alpha)+omega2^(-alpha)*(x/A)^(1-alpha)/(1-alpha))
}

vp<- function(x, A)
{
  return(omega2^(-alpha)*(x/A)^(-1))
}

y <- seq(0,3000,by=0.01)
plot(y,v(y,alpha ), type="l", lwd=2)
lines(y,v(y,0.4))
lines(y,v(y,0.5))
lines(y,v(y,0.6))

plot(y,vp(A*y, A), type="l", lwd=2, ylab = "Derivative of the value function", xlab="Capital value")
lines(y,vp(A*y,1) , type="l", lwd=2, lty=2, col="blue")
lines(y,vp(A*y,2) , type="l", lwd=2 , lty=2, col="red")
#lines(y,vp(y,3) , type="l", lwd=2, , lty=2, col="green")
legend("topright", legend = c("A = 0.5","A = 1", "A = 2"), col = c("black", "blue", "red"), lty = c(1,2,2), lwd=c(2,2,2), cex=0.7)


