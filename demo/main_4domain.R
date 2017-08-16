library(DNAr)

#
# 4-domain Examples
#

run_ApBeC_4d <- function() {
    result <- react_4domain(
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
    behavior <- result$behavior[,1:(3 + 1)]
}

run_AeB_4d <- function() {
    result <- react_4domain(
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
    behavior <- result$behavior[,1:(2 + 1)]
}

run_Lotka_4d <- function() {
    result <- react_4domain(
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
    behavior <- result$behavior[,1:(2 + 1)]
}

run_consensus_4d <- function() {
    result <- react_4domain(
        species   = c('X', 'Y', 'B'),
        ci        = c(0.7 * 80e-9, 0.3 * 80e-9, 0.0),
        reactions = c('X + Y -> 2B',
                      'B + X -> 2X',
                      'B + Y -> 2Y'),
        ki        = c(2e3, 2e3, 2e3),
        qmax      = 1e6,
        cmax      = 10e-6,
        alpha     = 1,
        beta      = 1,
        t         = seq(0, 54000, 5)
    )

    # save_behavior_csv(result$behavior, '../consensus')
    #
    # save_reactions_txt(
    #     species   = result$species,
    #     cis       = result$ci,
    #     reactions = result$reactions,
    #     kis       = result$ki,
    #     filename  = '../consensus'
    # )
    #
    # save_dsd_script(
    #     species   = c('X', 'Y', 'B'),
    #     ci        = c(0.7 * 80e-9, 0.3 * 80e-9, 0.0),
    #     reactions = c('X + Y -> 2B',
    #                   'B + X -> 2X',
    #                   'B + Y -> 2Y'),
    #     ki        = c(2e3, 2e3, 2e3),
    #     qmax      = 1e6,
    #     cmax      = 10e-6,
    #     alpha     = 1,
    #     beta      = 1,
    #     t         = seq(0, 54000, 5),
    #     filename  = '../consensus.dsd'
    # )

    behavior <- result$behavior[,1:(3 + 1)]
    return(behavior)
}

#behavior <- run_ApBeC_4d()
#behavior <- run_AeB_4d()
#behavior <- run_Lotka_4d()
behavior <- run_consensus_4d()

plot_behavior(behavior)
