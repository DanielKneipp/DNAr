library(DNAr)
context('Import and export tools')

temp_txt_filename <- 'data/temp_txt'
temp_csv_filename <- 'data/temp_csv'
dummy_behavior_csv_filename <- 'data/dummy_behavior'
dummy_behavior <- data.frame(
    'A' = seq(1, 10, length.out = 10),
    'B' = seq(0, 1, length.out = 10)
)

export_txt <- function(parms, ok_filename, first_time = FALSE) {
    do.call(save_reactions_txt, parms)
    if(!first_time) {
        filename_with_ext <- paste(parms$filename, '.txt', sep = '')
        out <- readChar(filename_with_ext, file.info(filename_with_ext)$size)
        out_ok <- readChar(ok_filename, file.info(ok_filename)$size)
        file.remove(filename_with_ext)
        return(list(out, out_ok))
    }
}

test_that('save_reactions_txt correctly exports the Consensus CRN', {
    parms <- list(
        species   = c('X', 'Y', 'B'),
        ci        = c(0.7 * 80e-9, 0.3 * 80e-9, 0.0),
        reactions = c('X + Y -> 2B',
                      'B + X -> 2X',
                      'B + Y -> 2Y'),
        kis        = c(2e3, 2e3, 2e3),
        filename   = temp_txt_filename
    )
    outs <- export_txt(parms, 'data/consensus.txt')
    expect_equal(outs[[1]], outs[[2]])
})

test_that('load_behavior_csv correctly imports a dummy behavior', {
    out_loaded <- load_behavior_csv(dummy_behavior_csv_filename)
    expect_equal(out_loaded, dummy_behavior)
})

test_that('save_behavior_csv correctly exports a dummy behavior', {
    save_behavior_csv(dummy_behavior, temp_csv_filename)
    out_loaded <- load_behavior_csv(temp_csv_filename)
    file.remove(paste(temp_csv_filename, '.csv', sep = ''))
    expect_equal(out_loaded, dummy_behavior)
})
