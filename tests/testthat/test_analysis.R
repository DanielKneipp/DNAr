library(DNAr)
context('CRN analysis')

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
            time_points = 1,
            behavior = behavior
        )

        expect_equal(d[['A']], 'd[A]/dt = (-1e-07 * 1000[A] * 1000[B])')
        expect_equal(d[['B']], 'd[B]/dt = (-1e-07 * 1000[A] * 1000[B])')
        expect_equal(d[['C']], 'd[C]/dt = (1e-07 * 1000[A] * 1000[B])')
    }
)

test_that(
    'analyze_behavior throws an exception when the time_points parameter is
    passed but behavior isn\'t',
    {
        expect_error(
            analyze_behavior(
                species   = c('A', 'B', 'C'),
                ci        = c(1e3, 1e3, 0),
                reactions = c('A + B -> C'),
                ki        = c(1e-7),
                time_points = 1
            )
        )
    }
)

test_that(
    'eval_derivative correct calculates the derivatives of A + B -> C',
    {
        parms <- list(
            species   = c('A', 'B', 'C'),
            ci        = c(1e3, 1e3, 0),
            reactions = c('A + B -> C'),
            ki        = c(1e-7),
            t         = seq(0, 72000, 10)
        )
        behavior <- do.call(react, parms)

        parms$t <- NULL
        parms$time_points <- 1
        parms$behavior <- behavior
        d <- do.call(analyze_behavior, parms)

        expect_equal(eval_derivative(d[['A']]), -0.1)
        expect_equal(eval_derivative(d[['B']]), -0.1)
        expect_equal(eval_derivative(d[['C']]), 0.1)
})

test_that(
    'eval_derivative correct calculates the derivatives of 2A + B -> C',
    {
        parms <- list(
            species   = c('A', 'B', 'C'),
            ci        = c(1e3, 1e3, 0),
            reactions = c('2A + B -> C'),
            ki        = c(1e-7),
            t         = seq(0, 72000, 10)
        )
        behavior <- do.call(react, parms)

        parms$t <- NULL
        parms$time_points <- 1
        parms$behavior <- behavior
        d <- do.call(analyze_behavior, parms)

        expect_equal(eval_derivative(d[['A']]), -200)
        expect_equal(eval_derivative(d[['B']]), -100)
        expect_equal(eval_derivative(d[['C']]), 100)
    }
)

test_that(
    'eval_derivative_part',
    {
        parms <- list(
            species   = c('A', 'B', 'C'),
            ci        = c(1e5, 1e3, 0),
            reactions = c('2A + B -> C',
                          '2A -> A'),
            ki        = c(1e-7,
                          1e-7),
            t         = seq(0, 72000, 10)
        )
        behavior <- do.call(react, parms)

        parms$t <- NULL
        parms$time_points <- 1
        parms$behavior <- behavior
        d <- do.call(analyze_behavior, parms)

        r <- eval_derivative_part(d[['A']])
        expect_equal(r[[1]], -2e+06)
        expect_equal(r[[2]], -1e+03)
        expect_equal(names(r)[[1]], '(-2e-07 * 1e+10[A]^2 * 1000[B])')
        expect_equal(names(r)[[2]], '(-1e-07 * 1e+10[A]^2)')

        r <- eval_derivative_part(d[['B']])
        expect_equal(r[[1]], -1e+06)
        expect_equal(names(r)[[1]], '(-1e-07 * 1e+10[A]^2 * 1000[B])')

        r <- eval_derivative_part(d[['C']])
        expect_equal(r[[1]], 1e+06)
        expect_equal(names(r)[[1]], '(1e-07 * 1e+10[A]^2 * 1000[B])')
    }
)

test_that(
    'compare_behaviors_nrmse correctly returns 0 if there is no difference',
    {
        b <- react(
            species   = c('A', 'B', 'C'),
            ci        = c(1e3, 1e3, 0),
            reactions = c('A + B -> C'),
            ki        = c(1e-7),
            t         = seq(0, 72000, length.out = 10)
        )

        diff <- compare_behaviors_nrmse(b, b)

        expect_equal(diff[['A']], 0)
        expect_equal(diff[['B']], 0)
        expect_equal(diff[['C']], 0)
    }
)

test_that(
    'compare_behaviors_nrmse correctly returns greater than 0
    when there is a difference',
    {
        b1 <- react(
            species   = c('A', 'B', 'C'),
            ci        = c(1e3, 1e3, 0),
            reactions = c('A + B -> C'),
            ki        = c(1e-7),
            t         = seq(0, 72000, length.out = 10)
        )

        b2 <- react(
            species   = c('A', 'B', 'C'),
            ci        = c(1e4, 1e3, 0),
            reactions = c('A + B -> C'),
            ki        = c(1e-7),
            t         = seq(0, 72000, length.out = 10)
        )

        diff <- round(compare_behaviors_nrmse(b1, b2), 2)

        expect_equal(diff[['A']], 8.78)
        expect_equal(diff[['B']], 0.27)
        expect_equal(diff[['C']], 0.27)
    }
)
