library(DNAr)
context('DSD script generation')

temp_dsd_filename <- 'temp_dsd.dsd'

EXPORT_MODE <- list(FIRST_TIME = 'first_time', NORMAL = 'normal')

export_dsd <- function(parms, ok_filename, mode = EXPORT_MODE$NORMAL) {
    switch (mode,
        'normal' = {
            do.call(save_dsd_script, parms)
            out <- readChar(parms$filename, file.info(parms$filename)$size)
            out_ok <- readChar(ok_filename, file.info(ok_filename)$size)
            file.remove(parms$filename)
            return(list(out, out_ok))
        },
        'first_time' = {
            parms$filename <- ok_filename
            do.call(save_dsd_script, parms)
            return(list(TRUE, TRUE))
        },
        stop('Mode does not exists')
    )
}

test_that(
    'save_dsd_script exports a dsd script of the reaction A + B -> C',
    {
        parms <- list(
            species   = c('A', 'B', 'C'),
            ci        = c(1e3, 1e3, 0),
            reactions = c('A + B -> C'),
            ki        = c(1e-7),
            qmax      = 1e-3,
            cmax      = 1e5,
            alpha     = 1,
            beta      = 1,
            t         = seq(0, 72000, 10),
            filename  = temp_dsd_filename
        )
        outs <- export_dsd(parms, 'data/ApBeC.dsd')
        expect_equal(outs[[1]], outs[[2]])
    }
)

test_that(
    'save_dsd_script exports a dsd script of the reaction A -> B + C',
    {
        parms <- list(
            species   = c('A', 'B', 'C'),
            ci        = c(1e3, 0, 0),
            reactions = c('A -> B + C'),
            ki        = c(1e-4),
            qmax      = 1e1,
            cmax      = 1e5,
            alpha     = 1,
            beta      = 1,
            t         = seq(0, 72000, 10),
            filename  = temp_dsd_filename
        )
        outs <- export_dsd(parms, 'data/AeBpC.dsd')
        expect_equal(outs[[1]], outs[[2]])
    }
)

test_that(
    'save_dsd_script exports a dsd script of the Consensus CRN',
    {
        parms <- list(
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
            t         = seq(0, 54000, 5),
            filename  = temp_dsd_filename
        )
        outs <- export_dsd(parms, 'data/consensus.dsd')
        expect_equal(outs[[1]], outs[[2]])
    }
)

test_that(
    'save_dsd_script correctly exports a dsd script of a CRN called dummy1:
    A + B -> 0; A -> 0; A + B -> C; C -> B',
    {
        parms <- list(
            species   = c('A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'),
            ci        = c(1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3, 1e3),
            reactions = c('A + B -> 0',
                          'C -> 0',
                          'D + E -> F',
                          'G -> H'),
            ki        = c(1e-7, 1e-4, 1e-7, 1e-4),
            qmax      = 1e1,
            cmax      = 1e5,
            alpha     = 1,
            beta      = 1,
            t         = seq(0, 35000, length.out = 100),
            filename  = temp_dsd_filename
        )
        outs <- export_dsd(parms, 'data/dummy1.dsd')
        expect_equal(outs[[1]], outs[[2]])
    }
)

test_that(
    'save_dsd_script correctly exports a dsd scripts of a CRN of
    formation and degradation reactions (called dummy2)',
    {
        parms <- list(
            species   = c('A', 'B', 'C', 'D', 'E', 'F'),
            ci        = c(1e3,  0,  1e3, 1e3, 1e3,  0),
            reactions = c('A -> 0',
                          '0 -> B',
                          '2C -> 0',
                          'D + E -> 0',
                          '0 -> 2F'),
            ki =        c(1, 2, 1e-4, 1e-4, 2),
            qmax      = 1e4,
            cmax      = 1e6,
            alpha     = 1,
            beta      = 1,
            t         = seq(0, 100, 1),
            filename  = temp_dsd_filename
        )
        outs <- export_dsd(parms, 'data/dummy2.dsd')
        expect_equal(outs[[1]], outs[[2]])
    }
)

test_that(
    'save_dsd_script correctly exports a dsd scripts of a CRN (called
     3prod) with reactions with 3 products.',
    {
        parms <- list(
            species   = c('sA', 'sB', 'sC', 'sD', 'sE', 'sF', 'sG', 'sH'),
            ci        = c( 0,    0,    0,    1e3,  0,    0,    1e3,  0),
            reactions = c('0 -> sA + sB + sC',
                          'sD -> 2sE + sF',
                          '2sG -> 3sH'),
            ki        = c(1e1, 1e-1, 1e-4),
            qmax      = 1e3,
            cmax      = 1e6,
            alpha     = 1,
            beta      = 1,
            t         = seq(0, 100, 1),
            filename  = temp_dsd_filename
        )
        outs <- export_dsd(parms, 'data/3prod.dsd')
        expect_equal(outs[[1]], outs[[2]])
    }
)
