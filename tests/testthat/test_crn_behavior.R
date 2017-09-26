library(DNAr)
context('CRN behavior')

run_reaction <- function(parms, filename, save = FALSE) {
    behavior <- do.call(react, parms)
    if(save) {
        save_behavior_csv(behavior, filename)
    } else {
        behavior_ok <- load_behavior_csv(filename)
        return(list(behavior, behavior_ok))
    }
}

test_that(
    'react reproduces correctly the behavior of the reaction A + B -> C',
    {
        parms <- list(
            species   = c('A', 'B', 'C'),
            ci        = c(1e3, 1e3, 0),
            reactions = c('A + B -> C'),
            ki        = c(1e-7),
            t         = seq(0, 72000, 10)
        )
        behaviors <- run_reaction(parms, 'behavior_ApBeC')
        expect_equal(behaviors[[1]], behaviors[[2]])
    }
)

test_that(
    'react reproduces correctly the behavior of the Lotka CRN',
    {
        parms <- list(
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
        behaviors <- run_reaction(parms, 'behavior_lotka')
        expect_equal(behaviors[[1]], behaviors[[2]])
    }
)

test_that(
    'react reproduces correctly the behavior of the CRN that represents
    an A + B -> C in the 4-domain approach',
    {
        parms <- list(
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
            t         = seq(0, 72000, 10)
        )
        behaviors <- run_reaction(parms, 'behavior_ApBeC_4d')
        expect_equal(behaviors[[1]], behaviors[[2]])
    }
)

test_that(
    'react reproduces correctly the behavior of the Consensus CRN ',
    {
        parms <- list(
            species   = c('X', 'Y', 'B'),
            ci        = c(0.7 * 80e-9, 0.3 * 80e-9, 0.0),
            reactions = c('X + Y -> 2B',
                          'B + X -> 2X',
                          'B + Y -> 2Y'),
            ki        = c(2e3, 2e3, 2e3),
            t         = seq(0, 54000, 5)
        )
        behaviors <- run_reaction(parms, 'behavior_consensus')
        expect_equal(behaviors[[1]], behaviors[[2]])
    }
)
