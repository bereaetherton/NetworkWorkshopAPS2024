README.md
================
Jacobo Robledo
2024-07-23

## SIR Modeling:

- *sir_model function:* Defines the differential equations for the SIR
  model.
- *Initial conditions and parameters:* Sets the initial state of the
  population and the model parameters (transmission and recovery rates).
- *Time frame:* Defines the period over which the simulation runs.
- *Solving the SIR model:* Uses the ode function from the deSolve
  package to solve the differential equations.
- *Plotting:* Visualizes the proportions of the susceptible, infected,
  and recovered populations over time.
- *Basic Reproduction Number (R0):* Calculates the basic reproduction
  number to understand the spread potential of the disease.

``` r
#------------SIR Modeling----------#

# First, download deSolve from https://desolve.r-forge.r-project.org/
# which is a differential equations package
#install.packages("deSolve") 
library("deSolve")
```

    ## Warning: package 'deSolve' was built under R version 4.3.3

``` r
# Define the SIR model function
sir_model <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    dS <- -beta * S * I # Rate of change of the Susceptible class
    dI <- beta * S * I - gamma * I # Rate of change of the Infected class
    dR <- gamma * I # Rate of change of the Recovered class
    
    # Return the rate of change for each class
    list(c(dS, dI, dR))
  })
}

# Initial conditions
initial_state <- c(S = 0.99, # 99% susceptible
                   I = 0.01, # 1% infected
                   R = 0.00) # 0% recovered

# Parameters
parameters <- c(beta = 0.3,  # Transmission rate
                gamma = 0.1) # Recovery rate

# Time frame
times <- seq(0, 160, by = 1) # Simulate for 160 days

# Solve the SIR problem
out <- ode(y = initial_state, times = times, func = sir_model, parms = parameters)

# After 130 days, for example...
out[130,]
```

    ##         time            S            I            R 
    ## 1.290000e+02 5.882805e-02 1.425075e-04 9.410294e-01

``` r
# This shows the proportion of the population in each class

# Convert the output to a data frame for plotting
out_df <- as.data.frame(out)
plot(out_df$time, out_df$S, type = "l", col = "blue", ylim = c(0, 1), xlab = "Time", ylab = "Proportion", lwd = 2,
     main = "SIR Model")
lines(out_df$time, out_df$I, col = "red", lwd = 2)
lines(out_df$time, out_df$R, col = "green", lwd = 2)
legend("right", legend = c("Susceptible", "Infected", "Recovered"), col = c("blue", "red", "green"), lwd = 2)
```

![](README_files/figure-gfm/cars-1.png)<!-- -->

``` r
# Basic Reproduction number
R_0 <- parameters[1] / parameters[2] # R0 = beta / gamma
R_0
```

    ## beta 
    ##    3

``` r
# This shows how many individuals (or plants) an individual will spread the disease to
```

## Epidemic Networks:

*Random Network:* Creates and plots a random network using the
Erdős–Rényi model. *Small World Network:* Creates and plots a
small-world network and calculates its average local transitivity (a
measure of the clustering coefficient). *Scale-Free Network:* Creates
and plots a scale-free network using the Barabási–Albert model and
calculates its average degree. *Degree Distribution:* Plots the degree
distribution histograms for the small-world and scale-free networks to
analyze their structure.

    ## 
    ## Attaching package: 'igraph'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     decompose, spectrum

    ## The following object is masked from 'package:base':
    ## 
    ##     union

![](README_files/figure-gfm/pressure-1.png)<!-- -->![](README_files/figure-gfm/pressure-2.png)<!-- -->

    ## [1] 0.5904762

![](README_files/figure-gfm/pressure-3.png)<!-- -->

    ## [1] 3.7

![](README_files/figure-gfm/pressure-4.png)<!-- -->![](README_files/figure-gfm/pressure-5.png)<!-- -->
