# Library
library(dplyr)

# Data with forest cover and population from 1990 to 2017
# see: http://dx.doi.org/10.18167/DVN1/AUBRRC
df_par <- read.table("data/scenarios/barnes.txt", header=TRUE, sep=",")

# Time-interval
int <- df_par$Year[-c(1)]-df_par$Year[-length(df_par$Year)]

## Deforestation
area_defor <- df_par$For[-length(df_par$For)]-df_par$For[-c(1)]
ann_defor <- area_defor/int

# New columns
df_par$D <- c(ann_defor,NA)
# Coefficients from Barnes model
df_par$lX <- log(df_par$D) - 0.607 * log(df_par$For) - 0.493 * log(df_par$Pop)

# Linear model to estimate beta0
mod <- lm(lX~1,data=df_par)
beta0 <- mod$coefficients[1] # -5.941
sigma2 <- var(mod$residuals) # 0.219

# Function to predict annual deforestation
D.func <- function (Forest, Pop, par) {
	beta0 <- par[1]
	beta1 <- par[2]
	beta2 <- par[3]
	V <- par[4]
	D <- exp(beta0+(V/2)+beta1*log(Forest)+beta2*log(Pop))
	return (D)
}
