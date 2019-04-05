# Library
library(dplyr)
library(readr)

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

# Projecting deforestation and forest cover
df_un <- read_csv("data/scenarios/un_pop.txt")
Year <- c(2017, seq(2020, 2100, by=5))
niter <- length(Year)-1
For <- For_np <- c(df_par$For[6], rep(NA,17))
Pop <- c(df_par$Pop[6], df_un$Pop[15:31])
D <- D_np <- rep(NA, length(Year))
par <- c(-5.941,0.607,0.493,0.219)
for (i in 1:niter) {
	D[i] <- D.func(For[i], Pop[i], par)
	D_np[i] <- D.func(For[i], Pop[1], par) # constant pop
	interval <- ifelse(Year[i]==2017, 3, 5)
	D_int <- D[i]*interval
	D_int_np <- D_np[i]*interval
	For[i+1] <- For[i]-D_int
	For_np[i+1] <- For_np[i]-D_int_np
}
D[niter+1] <- D.func(For[niter+1], Pop[niter+1], par)
D_np[niter+1] <- D.func(For[niter+1], Pop[1], par)
df_forest <- tibble(Year, For, Pop, For_np, D, D_np)
write_csv(df_forest, "output/for_proj.csv")



