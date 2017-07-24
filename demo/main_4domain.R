library(DNAr)

#
# 4-domain Examples
#

run_ApBeC_4d <- function() {
    behavior <- react_4domain(
        species   = c('A', 'B', 'C'),
        ci        = c(1e3, 1e3, 0),
        reactions = c('A + B -> C'),
        ki        = c(1e-7),
        qmax      = 1e-3,
        cmax      = 1e5,
        alpha     = 1,
        beta      = 1,
        t         = seq(0, 72000, 10)
    )
}

run_AeB_4d <- function() {
    behavior <- react_4domain(
        species   = c('A', 'B'),
        ci        = c(1e4, 0),
        reactions = c('A -> B'),
        ki        = c(5e-5 / 1e5),
        qmax      = 1e-5,
        cmax      = 1e5,
        alpha     = 1,
        beta      = 1,
        t         = seq(0, 72000, 10)
    )
}

run_Lotka_4d <- function() {
    behavior <- react_4domain(
        species   = c('X1', 'X2'),
        ci        = c(20e-9, 10e-9),
        reactions = c('X1 + X2 -> 2X2',
                      'X1 -> 2X1',
                      'X2 -> 0'),
        ki        = c(5e5,
                      1/300,
                      1/300),
        qmax      = 1e6,
        cmax      = 10e-6,
        alpha     = 1,
        beta      = 1,
        t         = seq(0, 12600, 1)
    )
}

#behavior <- run_ApBeC_4d()
#behavior <- run_AeB_4d()
behavior <- run_Lotka_4d()

plot_behavior(behavior)
