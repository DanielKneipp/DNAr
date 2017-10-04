library(DNAr)

#
# Examples
#

run_ApBeC <- function() {
    behavior <- react(
        species   = c('A', 'B', 'C'),
        ci        = c(1e3, 1e3, 0),
        reactions = c('A + B -> C'),
        ki        = c(1e-7),
        t         = seq(0, 72000, 10)
    )
    return(behavior)
}

run_lotka <- function() {
    behavior <- react(
        species   = c('A', 'B'),
        ci        = c(2, 1),
        reactions = c('A + B -> 2B',
                      'A -> 2A',
                      'B -> 0'),
        ki        = c(1.5,
                      1,
                      1),
        t         = seq(0, 45, 0.1)
    )
    return(behavior)
}

run_lotka_scaled <- function() {
    behavior <- react(
        species   = c('A', 'B'),
        ci        = c(20e-9, 10e-9),
        reactions = c('A + B -> 2B',
                      'A -> 2A',
                      'B -> 0'),
        ki        = c(5e5,
                      1/300,
                      1/300),
        t         = seq(0, 12600, 1)
    )
    return(behavior)
}

run_ApBeC_4domain <- function() {
    behavior <- react(
        species   = c('A', 'B', 'C', 'L', 'H', 'W', 'O', 'T'),
        ci        = c(1e3, 1e3, 0.0, 1e5, 0.0, 1e5, 0.0, 1e5),
        reactions = c('A + L -> H + W',
                      'H + W -> A + L',
                      'B + H -> O',
                      'O + T -> C'),
        ki        = c(1e-7,
                      1e-3,
                      1e-3,
                      1e-3),
        t         = seq(0, 72000, 1)
    )
    return(behavior)
}

run_origonator <- function() {
    # It is not working properly yet
    behavior <- react(
        species   = c('A', 'B', 'C'),
        ci        = c(0.8e-9, 4.5e-9, 0.8e-9),
        reactions = c('B -> A',
                      'B + A -> 0',
                      'A -> 2A + C',
                      '2A -> 0',
                      'C -> B',
                      'C -> 0'),
        ki        = c(3.889e-7,
                      5e5,
                      0.00232,
                      12500,
                      0.00198,
                      1.195e-5),
        t         = seq(0, 360000, 1)
    )
    return(behavior)
}

run_rossler <- function() {
    # It is not working properly yet
    behavior <- react(
        species   = c('A', 'B', 'C'),
        ci        = c(1.8e-9, 1.8e-9, 1.8e-9),
        reactions = c('A -> 2A',
                      '2A -> A',
                      'B + A -> 2B',
                      'B -> 0',
                      'A + C -> 0',
                      'C -> 2C',
                      '2C -> C'),
        ki        = c(0.002,
                      200000,
                      400000,
                      0.000667,
                      400000,
                      0.001,
                      200000),
        t         = seq(0, 144000, 1)
    )
    return(behavior)
}

run_consensus <- function() {
    behavior <- react(
        species   = c('X', 'Y', 'B'),
        ci        = c(0.7 * 80e-9, 0.3 * 80e-9, 0.0),
        reactions = c('X + Y -> 2B',
                      'B + X -> 2X',
                      'B + Y -> 2Y'),
        ki        = c(2e3, 2e3, 2e3),
        t         = seq(0, 54000, 5)
    )
}

#behavior <- run_ApBeC()
#behavior <- run_lotka()
#behavior <- run_lotka_scaled()
#behavior <- run_ApBeC_4domain()
#behavior <- run_origonator()
#behavior <- run_rossler()
behavior <- run_consensus()

plot_behavior(behavior)
