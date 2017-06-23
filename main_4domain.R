source('4domain_reactor.R')

#
# 4-domain Examples
#

run_ApBeC_4d <- function() {
    behaviour <- react_4domain(
        species   = c('A', 'B', 'C'),
        ci        = c(1e3, 1e3, 0),
        reactions = c('A + B -> C'),
        qi        = c(1e-7),
        qmax      = c(1e-3),
        cmax      = c(1e5),
        t         = seq(0, 72000, 10)
    )
}

run_AeB_4d <- function() {
    behaviour <- react_4domain(
        species   = c('A', 'B'),
        ci        = c(1e4, 0),
        reactions = c('A -> B'),
        qi        = c(5e-5 / 1e5),
        qmax      = c(1e-5),
        cmax      = c(1e5),
        t         = seq(0, 72000, 10)
    )
}

#behaviour <- run_ApBeC_4d()
behaviour <- run_AeB_4d()

plot_behaviour(behaviour)
