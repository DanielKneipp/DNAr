library(DNAr)

#
# 4-domain Examples
#

run_ApBeC_4d <- function() {
    behaviour <- react_4domain(
        species   = c('A', 'B', 'C'),
        ci        = c(1e3, 1e3, 0),
        reactions = c('A + B -> C'),
        ki        = c(1e-7),
        qmax      = 1e-3,
        cmax      = 1e5,
        t         = seq(0, 72000, 10)
    )
}

run_AeB_4d <- function() {
    behaviour <- react_4domain(
        species   = c('A', 'B'),
        ci        = c(1e4, 0),
        reactions = c('A -> B'),
        ki        = c(5e-5 / 1e5),
        qmax      = 1e-5,
        cmax      = 1e5,
        t         = seq(0, 72000, 10)
    )
}

run_Lotka_4d <- function() {
    behaviour <- react_4domain(
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
        t         = seq(0, 12600, 1)
    )
}

run_neuron_4d <- function() {
    behaviour <- react_4domain(
        species   = c('X1', 'X2', 'X3', 'X4', 'E1', 'E2', 'I1', 'I2'),
        ci        = c(0, 0, 0, 0, 0.5 * 5e-7, 0.5 * 5e-7, 10 * 5e-7, 5 * 5e-7),
        reactions = c('I2 -> 2I2',
                      'I2 -> X3',
                      'I1 -> 2I1',
                      'I1 -> X1',
                      'X1 + E1 -> X2 + E2',
                      'X3 + E2 -> X4 + E1',
                      'X2 -> 0',
                      'X4 -> 0'),
        ki        = c(0.001, 0.001, 0.001, 0.001, 1e6, 1e6, 0.01, 0.01),
        qmax      = 10e6,
        cmax      = 10e-5,
        t         = seq(0, 1000, 10)
    )
}

#behaviour <- run_ApBeC_4d()
#behaviour <- run_AeB_4d()
#behaviour <- run_Lotka_4d()
behaviour <- run_neuron_4d()

plot_behaviour(behaviour)
