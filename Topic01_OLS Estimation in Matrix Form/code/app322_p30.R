# Clear workspace
rm(list = ls())

# Clear plots
check_dev <- dev.list()
if(!is.null(check_dev)){
  dev.off(dev.list()["RStudioGD"])  
}

# Get location of current script
fileloc <- dirname(rstudioapi::getSourceEditorContext()$path)

# Set working directory to script location
setwd(fileloc)

# Remove fileloc object
rm(fileloc)

# Load packages
library(tidyverse)

# Load data on Application 3.2.2
table_f3.1 <- read.csv("TableF3-1.csv")
mod_tab3.1 <- lm(RealInv ~ Trend + RealGNP + Interest + Infl, 
                 data = table_f3.1)
summary(mod_tab3.1)

# Verify "manually" 
# (Greene gets the values wrong on p. 33, 
# even the online "correction" is wrong)
X <- as.matrix(select(table_f3.1, c("Trend", "RealGNP", "Interest", "Infl")))
X
X <- cbind(rep(1,15), X)
colnames(X)[1] <- "Intercept"
X

XprimeX <- t(X) %*% X
XprimeX

Y <- as.matrix(select(table_f3.1, "RealInv"))
Y
XprimeY <- t(X) %*% Y
colnames(XprimeY) <- "XprimeY"
XprimeY

b <- solve(XprimeX) %*% XprimeY

b

summary(mod_tab3.1)

# Construct the M matrix ("residual maker")
M <- diag(15) - X %*% solve(XprimeX) %*% t(X)
M

M - t(M) # Symmetry
M %*% M - M # Idempotency

as.matrix(mod_tab3.1$residuals)

e <- M %*% Y
colnames(e) <- "resid"
e

# Check whether MX = 0
M %*% X

# Yhat
Yhat = X %*% b
colnames(Yhat) <- "Yhat"
Yhat

# Check Y = Yhat + e
Y - Yhat - e

# Check Yhat'*e = 0
t(Yhat) %*% e

# Define the projection matrix P
P <- diag(15) - M
P

P - t(P) # Symmetry
P %*% P - P # Idempotency

Yhat2 <- P %*% Y
Yhat - Yhat2

# Check PM = MP = 0
P %*% M
M %*% P

# Check PX = X
P %*% X - X

# Check y = Py + My
Y - P %*% Y - M %*% Y

# Store summary of estimated model
smr <- summary(mod_tab3.1)
smr

# Get K (the number of estimated regression coefficients, incl. intercept)
K <- length(coef(mod_tab3.1))

# Get n (number of obs.)
n <- length(resid(mod_tab3.1))

# Sum of squared residuals (SSE)
SSE <- sum(resid(mod_tab3.1)^2)

# Residual (regression) standard error
s <- sqrt(SSE/(n - K))

# Get diagonal elements of (X'X)^{-1}
S_kk <- diag(solve(XprimeX))

# Set level of significance
alpha <- 0.05

# Calculate critical t-value
t_crit <- qt(1 - alpha/2, n - K)

# Calculate t-stats
t_stats <- smr$coefficients[,1]/(s * sqrt(S_kk))

confint_left <- smr$coefficients[,1] - t_crit * s * sqrt(S_kk)
confint_right <- smr$coefficients[,1] + t_crit * s * sqrt(S_kk)

confint(mod_tab3.1)

# Calculate p-values of estimated coefficients
p_vals <- 2*pt(-abs(t_stats), n - K)
p_vals

summary(mod_tab3.1)
