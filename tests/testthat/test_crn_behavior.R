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
    'react reproduces correctly the behavior of the reaction 2A -> 2B',
    {
        parms <- list(
            species   = c('A', 'B'),
            ci        = c(1e3, 0),
            reactions = c('2A -> 2B'),
            ki        = c(1e-7),
            t         = seq(0, 72000, 10)
        )
        behaviors <- run_reaction(parms, 'behavior_2Ae2B')
        expect_equal(behaviors[[1]], behaviors[[2]])
    }
)

test_that(
    'react reproduces correctly the behavior of the reaction 2A -> 2B + A',
    {
        parms <- list(
            species   = c('A', 'B'),
            ci        = c(1e3, 0),
            reactions = c('2A -> 2B + A'),
            ki        = c(1e-7),
            t         = seq(0, 72000, 10)
        )
        behaviors <- run_reaction(parms, 'behavior_2Ae2BpA')
        expect_equal(behaviors[[1]], behaviors[[2]])
    }
)

test_that(
    'react reproduces correctly the behavior of the reaction A -> 2B',
    {
        parms <- list(
            species   = c('A', 'B'),
            ci        = c(1e3, 0),
            reactions = c('A -> 2B'),
            ki        = c(1e-7),
            t         = seq(0, 72000, 10)
        )
        behaviors <- run_reaction(parms, 'behavior_Ae2B')
        expect_equal(behaviors[[1]], behaviors[[2]])
    }
)

test_that(
    'react reproduces correctly the behavior of the reaction 0 -> A',
    {
        parms <- list(
            species   = c('A'),
            ci        = c(0),
            reactions = c('0 -> A'),
            ki        = c(2),
            t         = seq(0, 100, 1)
        )
        behaviors <- run_reaction(parms, 'behavior_0eA')
        expect_equal(behaviors[[1]], behaviors[[2]])
    }
)

test_that(
    'react reproduces correctly the behavior of the reaction A -> 0',
    {
        parms <- list(
            species   = c('A'),
            ci        = c(1e3),
            reactions = c('A -> 0'),
            ki        = c(2),
            t         = seq(0, 100, 1)
        )
        behaviors <- run_reaction(parms, 'behavior_Ae0')
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
    'react reproduces correctly the behavior of the Consensus CRN',
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

test_that(
    'analyze_behavior correctly returns the derivatives of A + B -> C',
    {
        d <- analyze_behavior(
            species   = c('A', 'B', 'C'),
            ci        = c(1e3, 1e3, 0),
            reactions = c('A + B -> C'),
            ki        = c(1e-7)
        )

        expect_equal(d[['A']], 'd[A]/dt = (-1e-07 * [A] * [B])')
        expect_equal(d[['B']], 'd[B]/dt = (-1e-07 * [A] * [B])')
        expect_equal(d[['C']], 'd[C]/dt = (1e-07 * [A] * [B])')
    }
)

test_that(
    'analyze_behavior correctly returns the derivatives of 0 -> A',
    {
        d <- analyze_behavior(
            species   = c('A'),
            ci        = c(0),
            reactions = c('0 -> A'),
            ki        = c(2)
        )

        expect_equal(d[['A']], 'd[A]/dt = (2)')
    }
)

test_that(
    'analyze_behavior correctly returns the derivatives of A -> 0',
    {
        d <- analyze_behavior(
            species   = c('A'),
            ci        = c(1e3),
            reactions = c('A -> 0'),
            ki        = c(2)
        )

        expect_equal(d[['A']], 'd[A]/dt = (-2 * [A])')
    }
)

test_that(
    'analyze_behavior correctly returns the derivatives of 2A -> 2B',
    {
        d <- analyze_behavior(
            species   = c('A', 'B'),
            ci        = c(1e3, 0),
            reactions = c('2A -> 2B'),
            ki        = c(1e-7)
        )

        expect_equal(d[['A']], 'd[A]/dt = (-2e-07 * [A]^2)')
        expect_equal(d[['B']], 'd[B]/dt = (2e-07 * [A]^2)')
    }
)

test_that(
    'analyze_behavior correctly returns the derivatives of 2A -> B',
    {
        d <- analyze_behavior(
            species   = c('A', 'B'),
            ci        = c(1e3, 0),
            reactions = c('2A -> B'),
            ki        = c(1e-7)
        )

        expect_equal(d[['A']], 'd[A]/dt = (-2e-07 * [A]^2)')
        expect_equal(d[['B']], 'd[B]/dt = (1e-07 * [A]^2)')
    }
)

test_that(
    'analyze_behavior correctly returns the derivatives of 2A -> 2B + A',
    {
        d <- analyze_behavior(
            species   = c('A', 'B'),
            ci        = c(1e3, 0),
            reactions = c('2A -> 2B + A'),
            ki        = c(1e-7)
        )

        expect_equal(d[['A']], 'd[A]/dt = (-1e-07 * [A]^2)')
        expect_equal(d[['B']], 'd[B]/dt = (2e-07 * [A]^2)')
    }
)

test_that(
    'analyze_behavior correctly returns the derivatives of 2A -> 2B + 2A',
    {
        d <- analyze_behavior(
            species   = c('A', 'B'),
            ci        = c(1e3, 0),
            reactions = c('2A -> 2B + 2A'),
            ki        = c(1e-7)
        )

        expect_equal(d[['A']], 'd[A]/dt = 0')
        expect_equal(d[['B']], 'd[B]/dt = (2e-07 * [A]^2)')
    }
)

test_that(
    'analyze_behavior correctly returns the derivatives of a CRN
    which the first reaction don\'t have all the species',
    {
        d <- analyze_behavior(
            species   = c('A', 'B', 'C'),
            ci        = c(1e3, 0, 0),
            reactions = c('A -> 2B',
                          'B -> 2C'),
            ki        = c(1e-7, 1e-7)
        )

        expect_equal(d[['A']], 'd[A]/dt = (-1e-07 * [A])')
        expect_equal(d[['B']], 'd[B]/dt = (2e-07 * [A]) + (-1e-07 * [B])')
        expect_equal(d[['C']], 'd[C]/dt = (2e-07 * [B])')
    }
)

test_that(
    'analyze_behavior correctly returns the derivatives of Consensus',
    {
        d <- analyze_behavior(
            species   = c('X', 'Y', 'B'),
            ci        = c(0.7 * 80e-9, 0.3 * 80e-9, 0.0),
            reactions = c('X + Y -> 2B',
                          'B + X -> 2X',
                          'B + Y -> 2Y'),
            ki        = c(2e3, 2e3, 2e3)
        )

        expect_equal(d[['X']], 'd[X]/dt = (-2000 * [X] * [Y]) + (2000 * [B] * [X])')
        expect_equal(d[['Y']], 'd[Y]/dt = (-2000 * [X] * [Y]) + (2000 * [B] * [Y])')
        expect_equal(d[['B']], 'd[B]/dt = (4000 * [X] * [Y]) + (-2000 * [B] * [X]) + (-2000 * [B] * [Y])')
    }
)

test_that(
    'analyze_behavior correctly returns the derivatives of A + B -> C
    with the concentration values at the time point 0',
    {
        behavior <- react(
            species   = c('A', 'B', 'C'),
            ci        = c(1e3, 1e3, 0),
            reactions = c('A + B -> C'),
            ki        = c(1e-7),
            t         = seq(0, 72000, 10)
        )
        d <- analyze_behavior(
            species   = c('A', 'B', 'C'),
            ci        = c(1e3, 1e3, 0),
            reactions = c('A + B -> C'),
            ki        = c(1e-7),
            time_point = 1,
            behavior = behavior
        )

        expect_equal(d[['A']], 'd[A]/dt = (-1e-07 * 1000[A] * 1000[B])')
        expect_equal(d[['B']], 'd[B]/dt = (-1e-07 * 1000[A] * 1000[B])')
        expect_equal(d[['C']], 'd[C]/dt = (1e-07 * 1000[A] * 1000[B])')
    }
)

test_that(
    'analyze_behavior throws an exception when the time_point parameter is
    passed but behavior isn\'t',
    {
        expect_error(
            analyze_behavior(
                species   = c('A', 'B', 'C'),
                ci        = c(1e3, 1e3, 0),
                reactions = c('A + B -> C'),
                ki        = c(1e-7),
                time_point = 1
            )
        )
    }
)
